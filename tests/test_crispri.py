#!/usr/bin/env python

from collections import Counter

import pandas as pd
import pytest
from Bio.SeqFeature import FeatureLocation, SeqFeature
from pydna.dseqrecord import Dseqrecord

from teemi.design.crispr_cas import SgRNAargs, revcomp
from teemi.design.crispri import (
    extract_sgRNAs_for_crispri,
    filter_crispri_guides,
    find_sgrna_hits_cas9_crispri,
)


def _cds(location, locus_tag):
    return SeqFeature(location, type="CDS", qualifiers={"locus_tag": [locus_tag]})


def _hitframe(rows):
    """Build a hitframe with the columns filter_crispri_guides expects."""
    return pd.DataFrame(
        rows,
        columns=[
            "sgrna_loc",
            "gene_strand",
            "sgrna_strand",
            "pam",
            "sgrna",
            "downstream",
            "gc",
            "off_target_count",
        ],
    )


# ---------------------------------------------------------------------------
# find_sgrna_hits_cas9_crispri
# ---------------------------------------------------------------------------


# A gene with exactly one Cas9 guide in its CDS and one in its promoter.
# Each region holds a single "CC", and neither region contains "GG", so the
# reverse-complement scan adds no further hits.
PRE = "T" * 10
UP = "CC" + "ATATATATATATATATATATA" + "TTTTTT"  # 29 bp: a 23 nt window + 6 bp
BODY = "CC" + "AAATTTAAATTTAAATTTAAA"  # 23 bp: exactly one window
POST = "A" * 10
GENOME = PRE + UP + BODY + POST
GENE_START, GENE_END = len(PRE + UP), len(PRE + UP + BODY)

UP_GUIDE = "TATATATATATATATATATA"  # revcomp of the upstream window
CDS_GUIDE = "TTTAAATTTAAATTTAAATT"  # revcomp of the CDS window


def _genome(sequence, start, end, strand, name="strain"):
    record = Dseqrecord(sequence)
    record.id = name
    record.features = [_cds(FeatureLocation(start, end, strand=strand), "G")]
    return record


def _plus_strand_genome():
    return _genome(GENOME, GENE_START, GENE_END, 1, "crispri_plus")


def _minus_strand_genome():
    """The same gene, in the reverse-complemented genome, on the minus strand."""
    return _genome(
        revcomp(GENOME),
        len(GENOME) - GENE_END,
        len(GENOME) - GENE_START,
        -1,
        "crispri_minus",
    )


def _find(record):
    return find_sgrna_hits_cas9_crispri(
        record, record.id, ["G"], Counter(), 13, revcomp,
        extension_to_promoter_region=len(UP),
    )


def test_find_sgrna_hits_crispri_reports_cds_and_promoter_guides():
    df = _find(_plus_strand_genome())

    assert list(df["region"]) == ["CDS", "upstream"]
    cds_row, upstream_row = df.iloc[0], df.iloc[1]

    # The CDS guide is measured from the start of the gene, counting into it:
    # match.end() + protospacer_len + pam_len == 0 + 20 + 3.
    assert cds_row["sgrna"] == CDS_GUIDE
    assert cds_row["pam"] == "TGG"
    assert cds_row["sgrna_loc"] == 23
    assert cds_row["locus_tag"] == "G"

    # The promoter guide's window ends 6 bp before the gene, so it is negative.
    assert upstream_row["sgrna"] == UP_GUIDE
    assert upstream_row["pam"] == "TGG"
    assert upstream_row["sgrna_loc"] == -6
    # the synthetic "_upstream" feature name is not leaked
    assert upstream_row["locus_tag"] == "G"

    # Cas9 guides pair with a CCN on the strand scanned, so they lie opposite
    # the gene's own strand.
    assert list(df["gene_strand"]) == [1, 1]
    assert list(df["sgrna_strand"]) == [-1, -1]


def test_find_sgrna_hits_crispri_sgrna_loc_is_strand_symmetric():
    # The same gene placed on either strand must give the same guides at the
    # same positions relative to the gene, with the strands mirrored.
    plus, minus = _find(_plus_strand_genome()), _find(_minus_strand_genome())

    for column in ["region", "sgrna", "pam", "sgrna_loc"]:
        assert list(plus[column]) == list(minus[column])
    assert list(plus["sgrna_strand"]) == [-1, -1]
    assert list(minus["sgrna_strand"]) == [1, 1]


def test_find_sgrna_hits_crispri_without_promoter_extension():
    df = find_sgrna_hits_cas9_crispri(
        _plus_strand_genome(), "s", ["G"], Counter(), 13, revcomp,
        extension_to_promoter_region=0,
    )

    assert list(df["region"]) == ["CDS"]
    assert list(df["sgrna"]) == [CDS_GUIDE]


def test_find_sgrna_hits_crispri_clamps_promoter_region_at_the_sequence_ends():
    # The gene starts 10 bp into the record, so a 1000 bp promoter region
    # cannot be taken in full; the search is clamped instead of failing.
    record = _genome(PRE + BODY, len(PRE), len(PRE + BODY), 1)

    df = find_sgrna_hits_cas9_crispri(
        record, "s", ["G"], Counter(), 13, revcomp,
        extension_to_promoter_region=1000,
    )

    # 10 bp of promoter is too short to hold a 23 nt guide window
    assert list(df["region"]) == ["CDS"]

    # same at the other end of the record, for a gene on the minus strand
    minus = _genome(BODY + PRE, 0, len(BODY), -1)
    df_minus = find_sgrna_hits_cas9_crispri(
        minus, "s", ["G"], Counter(), 13, revcomp,
        extension_to_promoter_region=1000,
    )
    assert list(df_minus["region"]) == ["CDS"]


def test_find_sgrna_hits_crispri_ignores_other_locus_tags_and_features():
    record = _plus_strand_genome()
    record.features.append(
        _cds(FeatureLocation(0, len(PRE), strand=1), "OTHER_LOCUS")
    )
    record.features.append(SeqFeature(FeatureLocation(0, 20, strand=1), type="gene"))

    df = _find(record)

    assert set(df["locus_tag"]) == {"G"}


def test_find_sgrna_hits_crispri_matches_all_locus_tags():
    df = find_sgrna_hits_cas9_crispri(
        _plus_strand_genome(), "s", ["all"], Counter(), 13, revcomp,
        extension_to_promoter_region=len(UP),
    )

    assert list(df["region"]) == ["CDS", "upstream"]


def test_find_sgrna_hits_crispri_counts_off_targets():
    # The CDS guide's seed occurs twice in the genome, so one hit is an
    # off-target; the upstream guide's seed occurs once.
    counter = Counter({"TTAAATTTAAATT": 2, "ATATATATATATA": 1})

    df = find_sgrna_hits_cas9_crispri(
        _plus_strand_genome(), "s", ["G"], counter, 13, revcomp,
        extension_to_promoter_region=len(UP),
    )

    assert list(df["off_target_count"]) == [1, 0]


def test_filter_crispri_guides_keeps_only_tss_window_and_non_template_strand():
    hitframe = _hitframe(
        [
            (-150, 1, -1, "AGG", "AAA", "TTT", 0.5, 0),  # upstream of the window
            (150, 1, -1, "AGG", "AAA", "TTT", 0.5, 0),  # downstream of the window
            (50, 1, 1, "AGG", "AAA", "TTT", 0.5, 0),  # template strand
            (-50, 1, -1, "AGG", "AAA", "TTT", 0.5, 0),  # keeper
        ]
    )
    args = SgRNAargs(
        Dseqrecord("ATGC"),
        ["LOCUS"],
        upstream_tss=100,
        dwstream_tss=100,
        target_non_template_strand=True,
        gc_lower=0.0,
        gc_upper=1.0,
        off_target_upper=10,
    )

    filtered = filter_crispri_guides(args, hitframe)

    assert list(filtered.index) == [3]
    assert filtered.iloc[0]["sgrna_loc"] == -50


def test_filter_crispri_guides_applies_sequence_gc_and_off_target_filters():
    hitframe = _hitframe(
        [
            (0, 1, 1, "AGG", "AAA", "TTT", 0.5, 0),  # keeper
            (0, 1, 1, "TGG", "AAA", "TTT", 0.5, 0),  # removed by pam_remove
            (0, 1, 1, "AGG", "GGGG", "TTT", 0.5, 0),  # removed by sgrna_remove
            (0, 1, 1, "AGG", "AAA", "CCCC", 0.5, 0),  # removed by downstream_remove
            (0, 1, 1, "AGG", "AAA", "TTT", 0.95, 0),  # gc above gc_upper
            (0, 1, 1, "AGG", "AAA", "TTT", 0.1, 0),  # gc below gc_lower
            (0, 1, 1, "AGG", "AAA", "TTT", 0.5, 5),  # too many off-targets
        ]
    )
    args = SgRNAargs(
        Dseqrecord("ATGC"),
        ["LOCUS"],
        pam_remove=["TGG"],
        sgrna_remove=["GGGG"],
        downstream_remove=["CCCC"],
        gc_lower=0.2,
        gc_upper=0.8,
        off_target_upper=2,
        upstream_tss=100,
        dwstream_tss=100,
        target_non_template_strand=False,  # gene/sgrna strand filter switched off
    )

    filtered = filter_crispri_guides(args, hitframe)

    assert list(filtered.index) == [0]


# ---------------------------------------------------------------------------
# extract_sgRNAs_for_crispri
# ---------------------------------------------------------------------------


def test_extract_sgrnas_for_crispri_find_step_only():
    df = extract_sgRNAs_for_crispri(
        SgRNAargs(
            _plus_strand_genome(),
            ["G"],
            cas_type="cas9",
            step=["find"],
            extension_to_promoter_region=len(UP),
        )
    )

    # sorted by position, so the promoter guide comes first
    assert list(df["sgrna_loc"]) == [-6, 23]
    assert list(df["sgrna"]) == [UP_GUIDE, CDS_GUIDE]
    assert list(df["region"]) == ["upstream", "CDS"]


def test_extract_sgrnas_for_crispri_filter_step_applies_the_tss_window():
    args = dict(
        cas_type="cas9",
        step=["find", "filter"],
        extension_to_promoter_region=len(UP),
    )

    # Both guides sit within +-100 bp of the gene start
    wide = extract_sgRNAs_for_crispri(
        SgRNAargs(_plus_strand_genome(), ["G"], upstream_tss=100, dwstream_tss=100, **args)
    )
    assert list(wide["sgrna_loc"]) == [-6, 23]

    # Narrowing the window downstream of the start drops the CDS guide
    narrow = extract_sgRNAs_for_crispri(
        SgRNAargs(_plus_strand_genome(), ["G"], upstream_tss=100, dwstream_tss=10, **args)
    )
    assert list(narrow["sgrna_loc"]) == [-6]


def test_extract_sgrnas_for_crispri_filter_step_keeps_non_template_strand():
    # Cas9 guides lie opposite the gene strand, which is the non-template
    # strand, so this filter keeps them both.
    filtered = extract_sgRNAs_for_crispri(
        SgRNAargs(
            _plus_strand_genome(),
            ["G"],
            cas_type="cas9",
            step=["find", "filter"],
            extension_to_promoter_region=len(UP),
            target_non_template_strand=True,
        )
    )

    assert list(filtered["sgrna"]) == [UP_GUIDE, CDS_GUIDE]


def test_extract_sgrnas_for_crispri_rejects_unsupported_cas_type():
    with pytest.raises(ValueError, match="not supported for CRISPRi"):
        extract_sgRNAs_for_crispri(
            SgRNAargs(_plus_strand_genome(), ["G"], cas_type="cas3", step=["find"])
        )


def test_extract_sgrnas_for_crispri_requires_the_find_step():
    with pytest.raises(ValueError, match="step must include 'find'"):
        extract_sgRNAs_for_crispri(
            SgRNAargs(_plus_strand_genome(), ["G"], cas_type="cas9", step=["filter"])
        )


def test_extract_sgrnas_for_crispri_requires_dseqrecord():
    with pytest.raises(ValueError):
        SgRNAargs("ATGC", ["LOCUS"])


def test_filter_crispri_guides_rejects_a_filter_without_its_column():
    # downstream_remove filters a "downstream" column, which none of the sgRNA
    # hit finders produce; say so instead of raising a bare KeyError.
    hitframe = _hitframe([(0, 1, 1, "AGG", "AAA", "TTT", 0.5, 0)]).drop(
        columns=["downstream"]
    )
    args = SgRNAargs(
        Dseqrecord("ATGC"), ["LOCUS"], downstream_remove=["CCCC"],
        upstream_tss=100, dwstream_tss=100,
    )

    with pytest.raises(ValueError, match="do not produce a 'downstream' column"):
        filter_crispri_guides(args, hitframe)


def test_find_sgrna_hits_crispri_gene_with_no_room_upstream():
    # A gene starting at the first base has no promoter region to search.
    record = _genome(BODY + POST, 0, len(BODY), 1)

    df = find_sgrna_hits_cas9_crispri(
        record, "s", ["G"], Counter(), 13, revcomp,
        extension_to_promoter_region=len(UP),
    )

    assert list(df["region"]) == ["CDS"]
