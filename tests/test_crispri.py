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


def test_find_sgrna_hits_crispri_positive_strand_gene():
    # Gene on the + strand: "CC" + 21 A's + 10 T's, preceded by 30 T's.
    # The only CC-match sits at offset 0 of the coding sequence, so the
    # 23 nt protospacer+PAM window is "CC" + "A" * 21 and its reverse
    # complement is "T" * 21 + "GG" -> sgrna "T" * 20, pam "TGG".
    genome = Dseqrecord("T" * 30 + "CC" + "A" * 21 + "T" * 10)
    genome.id = "crispri_plus_strand"
    genome.features = [_cds(FeatureLocation(30, 63, strand=1), "LOCUS_PLUS")]

    df = find_sgrna_hits_cas9_crispri(
        genome,
        genome.id,
        ["LOCUS_PLUS"],
        Counter({"T" * 13: 1}),
        13,
        revcomp,
        extension_to_promoter_region=0,
    )

    assert len(df) == 1
    row = df.iloc[0]
    assert row["strain_name"] == "crispri_plus_strand"
    assert row["locus_tag"] == "LOCUS_PLUS"
    assert row["gene_loc"] == 31  # 0-based feature start 30 -> 1-based
    assert row["gene_strand"] == 1
    assert row["sgrna_strand"] == 1  # watson guide for a + strand gene
    assert row["sgrna"] == "T" * 20
    assert row["pam"] == "TGG"
    assert row["gc"] == 0.0
    # sgrna_loc for + strand genes is covered by the mirror-image test below.
    assert row["sgrna_seed_sequence"] == "T" * 13
    assert row["off_target_count"] == 0  # counter value 1, minus the guide itself


@pytest.mark.xfail(
    strict=True,
    reason="sgrna_loc formula is picked by strand label, not by which sequence "
    "is scanned, so + strand genes are measured from the 3' end (crispri.py:165-169)",
)
def test_find_sgrna_hits_crispri_sgrna_loc_is_strand_symmetric():
    # The same gene placed on the + strand and, in the reverse-complemented
    # genome, on the - strand has an identical promoter-extended coding
    # sequence, so the guide must get the same position relative to the TSS.
    plus_genome = Dseqrecord("T" * 30 + "CC" + "A" * 21 + "T" * 10)
    plus_genome.features = [_cds(FeatureLocation(30, 63, strand=1), "GENE")]
    minus_genome = Dseqrecord(revcomp(str(plus_genome.seq)))
    minus_genome.features = [_cds(FeatureLocation(0, 33, strand=-1), "GENE")]

    def hits(genome):
        return find_sgrna_hits_cas9_crispri(
            genome, "strain", ["GENE"], Counter(), 13, revcomp,
            extension_to_promoter_region=5,
        )

    plus_hits, minus_hits = hits(plus_genome), hits(minus_genome)

    assert list(plus_hits["sgrna"]) == list(minus_hits["sgrna"]) == ["T" * 20]
    assert list(plus_hits["sgrna_loc"]) == list(minus_hits["sgrna_loc"])


def test_find_sgrna_hits_crispri_negative_strand_gene_uses_promoter_extension():
    # Gene on the - strand at [10, 30). With a 5 nt promoter extension the
    # genomic slice is [10, 35) = "AA" + "T" * 21 + "GG"; the coding sequence
    # is its reverse complement, "CC" + "A" * 21 + "TT".
    genome = Dseqrecord("C" * 10 + "AA" + "T" * 21 + "GG" + "C" * 5)
    genome.id = "crispri_minus_strand"
    genome.features = [_cds(FeatureLocation(10, 30, strand=-1), "LOCUS_MINUS")]

    df = find_sgrna_hits_cas9_crispri(
        genome,
        genome.id,
        ["LOCUS_MINUS"],
        Counter({"T" * 13: 1}),
        13,
        revcomp,
        extension_to_promoter_region=5,
    )

    assert len(df) == 1
    row = df.iloc[0]
    assert row["locus_tag"] == "LOCUS_MINUS"
    assert row["gene_loc"] == 11
    assert row["gene_strand"] == -1
    assert row["sgrna_strand"] == -1  # watson guide for a - strand gene
    assert row["sgrna"] == "T" * 20
    assert row["pam"] == "TGG"
    # match.end() + protospacer_len + pam_len - extension == 0 + 20 + 3 - 5:
    # the 3' end of the guide window sits 18 nt downstream of the TSS.
    assert row["sgrna_loc"] == 18
    assert row["off_target_count"] == 0


def test_find_sgrna_hits_crispri_ignores_other_locus_tags():
    genome = Dseqrecord("T" * 30 + "CC" + "A" * 21 + "T" * 10)
    genome.id = "crispri_other_locus"
    genome.features = [
        _cds(FeatureLocation(30, 63, strand=1), "LOCUS_PLUS"),
        SeqFeature(FeatureLocation(0, 30, strand=1), type="gene"),
    ]

    df = find_sgrna_hits_cas9_crispri(
        genome,
        genome.id,
        ["SOME_OTHER_LOCUS"],
        Counter(),
        13,
        revcomp,
        extension_to_promoter_region=0,
    )

    assert df.empty
    assert list(df.columns)[:2] == ["strain_name", "locus_tag"]


def test_find_sgrna_hits_crispri_skips_guide_with_truncated_pam(capsys):
    # The coding sequence is only 22 nt, so the window starting at the CC match
    # yields a full 20 nt protospacer but a 2 nt PAM -> the hit is discarded.
    genome = Dseqrecord("CC" + "A" * 20)
    genome.id = "crispri_truncated"
    genome.features = [_cds(FeatureLocation(0, 22, strand=1), "LOCUS_SHORT")]

    df = find_sgrna_hits_cas9_crispri(
        genome,
        genome.id,
        ["LOCUS_SHORT"],
        Counter(),
        13,
        revcomp,
        extension_to_promoter_region=0,
    )

    assert df.empty
    assert "outside the designated border in LOCUS_SHORT" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# filter_crispri_guides
# ---------------------------------------------------------------------------


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
    genome = Dseqrecord("T" * 30 + "CC" + "A" * 21 + "T" * 10)
    genome.id = "crispri_pipeline_find"
    genome.features = [_cds(FeatureLocation(30, 63, strand=1), "LOCUS_PLUS")]

    args = SgRNAargs(
        genome,
        ["LOCUS_PLUS"],
        cas_type="cas9",
        step=["find"],
        off_target_seed=13,
    )

    df = extract_sgRNAs_for_crispri(args)

    assert len(df) == 1
    row = df.iloc[0]
    assert row["sgrna"] == "T" * 20
    assert row["pam"] == "TGG"
    assert row["sgrna_loc"] == 30
    assert row["off_target_count"] == 0


def test_extract_sgrnas_for_crispri_filter_step_drops_template_strand_guides():
    genome = Dseqrecord("T" * 30 + "CC" + "A" * 21 + "T" * 10)
    genome.id = "crispri_pipeline_filter"
    genome.features = [_cds(FeatureLocation(30, 63, strand=1), "LOCUS_PLUS")]

    find_only = extract_sgRNAs_for_crispri(
        SgRNAargs(genome, ["LOCUS_PLUS"], cas_type="cas9", step=["find"])
    )
    assert len(find_only) == 1
    # The single guide lies on the same strand as the gene ...
    assert find_only.iloc[0]["gene_strand"] == find_only.iloc[0]["sgrna_strand"]

    filtered = extract_sgRNAs_for_crispri(
        SgRNAargs(
            genome,
            ["LOCUS_PLUS"],
            cas_type="cas9",
            step=["find", "filter"],
            target_non_template_strand=True,
        )
    )

    # ... so requiring the non-template strand removes it.
    assert filtered.empty


def test_extract_sgrnas_for_crispri_sorts_by_sgrna_location():
    # Two CC matches in the coding sequence give two guides at different
    # distances from the TSS; the pipeline returns them sorted ascending.
    genome = Dseqrecord("CC" + "A" * 21 + "CC" + "A" * 21 + "T" * 6)
    genome.id = "crispri_sorted"
    genome.features = [_cds(FeatureLocation(0, len(genome.seq), strand=1), "LOCUS_TWO")]

    df = extract_sgRNAs_for_crispri(
        SgRNAargs(genome, ["LOCUS_TWO"], cas_type="cas9", step=["find"])
    )

    assert len(df) == 2
    assert list(df["sgrna_loc"]) == sorted(df["sgrna_loc"])
    # Both guides are the reverse complement of "CC" + 21 A's.
    assert set(df["sgrna"]) == {"T" * 20}


def test_extract_sgrnas_for_crispri_promoter_extension_adds_upstream_guides():
    # The CC match lies 10 nt upstream of the CDS start, so it is only seen
    # once the promoter region is included.
    genome = Dseqrecord("T" * 20 + "CC" + "A" * 21 + "G" * 30)
    genome.id = "crispri_extension"
    genome.features = [_cds(FeatureLocation(33, len(genome.seq), strand=1), "LOCUS_EXT")]

    without_extension = extract_sgRNAs_for_crispri(
        SgRNAargs(
            genome,
            ["LOCUS_EXT"],
            cas_type="cas9",
            step=["find"],
            extension_to_promoter_region=0,
        )
    )
    with_extension = extract_sgRNAs_for_crispri(
        SgRNAargs(
            genome,
            ["LOCUS_EXT"],
            cas_type="cas9",
            step=["find"],
            extension_to_promoter_region=20,
        )
    )

    # Without the extension only the CDS body is scanned, so the guide sitting
    # in the promoter is never seen; every guide found there ends in G's.
    assert not without_extension.empty
    assert all(sgrna.endswith("G") for sgrna in without_extension["sgrna"])
    assert "T" * 20 not in set(without_extension["sgrna"])
    # Including the promoter region surfaces the upstream guide as well.
    assert "T" * 20 in set(with_extension["sgrna"])
    upstream = with_extension[with_extension["sgrna"] == "T" * 20].iloc[0]
    assert upstream["pam"] == "TGG"
    assert upstream["gene_loc"] == 34


def test_extract_sgrnas_for_crispri_requires_dseqrecord():
    with pytest.raises(ValueError):
        SgRNAargs("ATGC", ["LOCUS"])
