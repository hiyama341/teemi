#!/usr/bin/env python

from collections import Counter

import pandas as pd
import pytest
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from pydna.dseqrecord import Dseqrecord

from teemi.design.crispr_cas import (
    SgRNAargs,
    extract_sgRNAs,
    filter_guides,
    find_off_target_hits,
    find_sgrna_hits_cas3,
    find_sgrna_hits_cas9,
    find_sgrna_hits_cas12a,
    revcomp,
)


def test_extract_sgrnas_handles_compound_feature_locations():
    genome = Dseqrecord("ATGCGG" * 80)
    genome.id = "test_compound_genome"

    compound_feature = SeqFeature(
        CompoundLocation(
            [
                FeatureLocation(30, 90, strand=-1),
                FeatureLocation(120, 180, strand=-1),
            ]
        ),
        type="CDS",
        qualifiers={"locus_tag": ["TEST_LOCUS_001"]},
    )
    genome.features = [compound_feature]

    args = SgRNAargs(
        genome,
        ["TEST_LOCUS_001"],
        cas_type="cas9",
        step=["find", "filter"],
        gc_upper=1.0,
        gc_lower=0.0,
        off_target_seed=13,
        off_target_upper=10,
    )

    sgrna_df = extract_sgRNAs(args)

    assert "locus_tag" in sgrna_df.columns


def test_extract_sgrnas_cas9_off_target_count_is_non_negative():
    genome = Dseqrecord("ATG" + "A" * 40 + "CC" + "G" * 40 + "ATGCGG" * 40)
    genome.id = "test_cas9_genome"

    feature = SeqFeature(
        FeatureLocation(0, len(genome.seq), strand=1),
        type="CDS",
        qualifiers={"locus_tag": ["TEST_LOCUS_002"]},
    )
    genome.features = [feature]

    args = SgRNAargs(
        genome,
        ["TEST_LOCUS_002"],
        cas_type="cas9",
        step=["find"],
        gc_upper=1.0,
        gc_lower=0.0,
        off_target_seed=13,
        off_target_upper=10,
    )

    sgrna_df = extract_sgRNAs(args)

    assert not sgrna_df.empty
    assert (sgrna_df["off_target_count"] >= 0).all()


def test_extract_sgrnas_cas9_reports_genomic_guide_strand():
    genome = Dseqrecord("A" * 40 + "CC" + "G" * 30 + "A" * 40)
    genome.id = "test_cas9_strand_genome"

    feature = SeqFeature(
        FeatureLocation(0, len(genome.seq), strand=-1),
        type="CDS",
        qualifiers={"locus_tag": ["TEST_LOCUS_005"]},
    )
    genome.features = [feature]

    args = SgRNAargs(
        genome,
        ["TEST_LOCUS_005"],
        cas_type="cas9",
        step=["find"],
        gc_upper=1.0,
        gc_lower=0.0,
        off_target_seed=13,
        off_target_upper=10,
    )

    sgrna_df = extract_sgRNAs(args)

    assert not sgrna_df.empty
    assert 1 in set(sgrna_df["sgrna_strand"])


def test_extract_sgrnas_cas12a_uses_requested_protospacer_length():
    genome = Dseqrecord("A" * 60 + "TTTA" + "C" * 40 + "G" * 40)
    genome.id = "test_cas12a_genome"

    feature = SeqFeature(
        FeatureLocation(0, len(genome.seq), strand=1),
        type="CDS",
        qualifiers={"locus_tag": ["TEST_LOCUS_003"]},
    )
    genome.features = [feature]

    args = SgRNAargs(
        genome,
        ["TEST_LOCUS_003"],
        cas_type="cas12a",
        step=["find"],
        gc_upper=1.0,
        gc_lower=0.0,
        off_target_seed=13,
        off_target_upper=10,
        protospacer_len=23,
    )

    sgrna_df = extract_sgRNAs(args)

    assert not sgrna_df.empty
    assert (sgrna_df["sgrna"].str.len() == 23).all()


def test_extract_sgrnas_cas12a_reports_genomic_guide_strand():
    genome = Dseqrecord("A" * 40 + "TTTA" + "C" * 30 + "A" * 40)
    genome.id = "test_cas12a_strand_genome"

    feature = SeqFeature(
        FeatureLocation(0, len(genome.seq), strand=-1),
        type="CDS",
        qualifiers={"locus_tag": ["TEST_LOCUS_004"]},
    )
    genome.features = [feature]

    args = SgRNAargs(
        genome,
        ["TEST_LOCUS_004"],
        cas_type="cas12a",
        step=["find"],
        gc_upper=1.0,
        gc_lower=0.0,
        off_target_seed=13,
        off_target_upper=10,
        protospacer_len=23,
    )

    sgrna_df = extract_sgRNAs(args)

    assert not sgrna_df.empty
    assert -1 in set(sgrna_df["sgrna_strand"])


def _cds(location, locus_tag):
    return SeqFeature(location, type="CDS", qualifiers={"locus_tag": [locus_tag]})


def _single_cds_genome(sequence, locus_tag, strand, record_id):
    genome = Dseqrecord(sequence)
    genome.id = record_id
    genome.features = [_cds(FeatureLocation(0, len(sequence), strand=strand), locus_tag)]
    return genome


# ---------------------------------------------------------------------------
# SgRNAargs
# ---------------------------------------------------------------------------


def test_sgrnaargs_rejects_input_that_is_not_a_dseqrecord():
    with pytest.raises(ValueError, match="Dseqrecord"):
        SgRNAargs("ATGCGG", ["TEST_LOCUS"])


def test_sgrnaargs_wraps_a_bare_locus_tag_and_reads_the_strain_name():
    genome = Dseqrecord("ATGCGG")
    genome.id = "strain_x"

    args = SgRNAargs(genome, "TEST_LOCUS")

    assert args.locus_tag == ["TEST_LOCUS"]
    assert args.strain_name == "strain_x"


# ---------------------------------------------------------------------------
# find_off_target_hits
# ---------------------------------------------------------------------------


def test_find_off_target_hits_cas3_takes_the_seed_downstream_of_the_pam():
    # Two TTC PAMs, each followed by the same 5 nt seed.
    contig = "AAA" + "TTC" + "GGGGG" + "TTC" + "GGGGG"

    counter = find_off_target_hits([contig], off_target_seed=5, cas_type="cas3")

    assert counter == Counter({"GGGGG": 2})


def test_find_off_target_hits_skips_seeds_truncated_by_the_contig_end():
    # The PAM is 2 nt from the end, so no full 5 nt seed can be read.
    counter = find_off_target_hits(["AAATTCGG"], off_target_seed=5, cas_type="cas3")

    assert counter == Counter()


def test_find_off_target_hits_cas9_takes_the_seed_upstream_of_the_pam():
    counter = find_off_target_hits(["AAAAAGG" + "TTTTTGG"], off_target_seed=5, cas_type="cas9")

    # GG at index 5 is preceded by "AAAAA"; GG at index 12 by "TTTTT".
    assert counter == Counter({"AAAAA": 1, "TTTTT": 1})


# ---------------------------------------------------------------------------
# find_sgrna_hits_cas9 / cas12a border cases
# ---------------------------------------------------------------------------


def test_find_sgrna_hits_cas9_discards_guide_whose_pam_runs_past_the_cds(capsys):
    # The CDS is only 22 nt, so the window starting at the CC match yields a
    # full 20 nt protospacer but only a 2 nt PAM.
    genome = _single_cds_genome("CC" + "A" * 20, "TEST_LOCUS_006", 1, "cas9_short_pam")

    df = find_sgrna_hits_cas9(genome, genome.id, ["TEST_LOCUS_006"], Counter(), 13, revcomp)

    assert df.empty
    assert "Pam was found outside designated locus_tag" in capsys.readouterr().out


def test_find_sgrna_hits_cas12a_discards_pam_at_the_very_end_of_the_cds():
    # The TTTA PAM ends exactly at the CDS end, leaving no protospacer at all.
    genome = _single_cds_genome("G" * 10 + "TTTA", "TEST_LOCUS_007", 1, "cas12a_no_spacer")

    df = find_sgrna_hits_cas12a(
        genome, genome.id, ["TEST_LOCUS_007"], Counter(), 13, revcomp, protospacer_len=23
    )

    assert df.empty


def test_find_sgrna_hits_cas12a_discards_protospacer_truncated_by_the_cds_end():
    # Only 5 nt follow the PAM, short of the requested 23 nt protospacer.
    genome = _single_cds_genome(
        "G" * 10 + "TTTA" + "C" * 5, "TEST_LOCUS_008", 1, "cas12a_short_spacer"
    )

    df = find_sgrna_hits_cas12a(
        genome, genome.id, ["TEST_LOCUS_008"], Counter(), 13, revcomp, protospacer_len=23
    )

    assert df.empty


# ---------------------------------------------------------------------------
# cas3
# ---------------------------------------------------------------------------


def test_extract_sgrnas_cas3_forward_strand_guide():
    # "TTC" PAM at the CDS start followed by a 34 nt protospacer.
    genome = _single_cds_genome("TTC" + "A" * 34 + "GGG", "TEST_LOCUS_009", 1, "cas3_plus")

    args = SgRNAargs(
        genome, ["TEST_LOCUS_009"], cas_type="cas3", step=["find"], off_target_seed=13
    )
    df = extract_sgRNAs(args)

    assert len(df) == 1
    row = df.iloc[0]
    assert row["strain_name"] == "cas3_plus"
    assert row["locus_tag"] == "TEST_LOCUS_009"
    assert row["pam"] == "TTC"
    assert row["sgrna"] == "A" * 34
    assert row["gc"] == 0.0
    assert row["gene_loc"] == 1
    assert row["gene_strand"] == 1
    assert row["sgrna_strand"] == 1
    # match.end() + protospacer_len + pam_len == 3 + 34 + 3; the same as the
    # mirror-image minus-strand genome in the next test
    assert row["sgrna_loc"] == 40
    assert row["sgrna_seed_sequence"] == "A" * 13
    # The seed occurs once in the genome, which is this guide itself.
    assert row["off_target_count"] == 0


def test_extract_sgrnas_cas3_reverse_strand_guide():
    # Reverse complement of the forward-strand genome above, with the CDS
    # annotated on the minus strand, so the same guide is found.
    genome = _single_cds_genome("CCC" + "T" * 34 + "GAA", "TEST_LOCUS_010", -1, "cas3_minus")

    args = SgRNAargs(
        genome, ["TEST_LOCUS_010"], cas_type="cas3", step=["find"], off_target_seed=13
    )
    df = extract_sgRNAs(args)

    assert len(df) == 1
    row = df.iloc[0]
    assert row["pam"] == "TTC"
    assert row["sgrna"] == "A" * 34
    assert row["gene_strand"] == -1
    assert row["sgrna_strand"] == -1
    # match.end() + protospacer_len + pam_len == 3 + 34 + 3
    assert row["sgrna_loc"] == 40
    assert row["off_target_count"] == 0


def test_extract_sgrnas_cas3_matches_all_locus_tags():
    genome = _single_cds_genome("TTC" + "A" * 34 + "GGG", "TEST_LOCUS_011", 1, "cas3_all")

    df = extract_sgRNAs(
        SgRNAargs(genome, ["all"], cas_type="cas3", step=["find"], off_target_seed=13)
    )

    assert list(df["locus_tag"]) == ["TEST_LOCUS_011"]


def test_find_sgrna_hits_cas3_discards_pam_at_the_very_end_of_the_cds(capsys):
    genome = _single_cds_genome("A" * 10 + "TTC", "TEST_LOCUS_012", 1, "cas3_no_spacer")

    df = find_sgrna_hits_cas3(genome, genome.id, ["TEST_LOCUS_012"], Counter(), 13, revcomp)

    assert df.empty
    assert "No sgRNA found for locus tag TEST_LOCUS_012" in capsys.readouterr().out


def test_find_sgrna_hits_cas3_discards_protospacer_truncated_by_the_cds_end():
    genome = _single_cds_genome(
        "A" * 10 + "TTC" + "G" * 5, "TEST_LOCUS_013", 1, "cas3_short_spacer"
    )

    df = find_sgrna_hits_cas3(genome, genome.id, ["TEST_LOCUS_013"], Counter(), 13, revcomp)

    assert df.empty


# ---------------------------------------------------------------------------
# filter_guides / pipeline steps
# ---------------------------------------------------------------------------


def test_filter_guides_excludes_rows_matching_the_removal_patterns():
    hitframe = pd.DataFrame(
        {
            "pam": ["AGG", "TGG", "AGG", "AGG", "AGG", "AGG", "AGG"],
            "sgrna": ["AAA", "AAA", "GGGG", "AAA", "AAA", "AAA", "AAA"],
            "downstream": ["TTT", "TTT", "TTT", "CCCC", "TTT", "TTT", "TTT"],
            "gc": [0.5, 0.5, 0.5, 0.5, 0.95, 0.1, 0.5],
            "off_target_count": [0, 0, 0, 0, 0, 0, 5],
        }
    )
    args = SgRNAargs(
        Dseqrecord("ATGCGG"),
        ["TEST_LOCUS"],
        pam_remove=["TGG"],
        sgrna_remove=["GGGG"],
        downstream_remove=["CCCC"],
        gc_lower=0.2,
        gc_upper=0.8,
        off_target_upper=2,
    )

    filtered = filter_guides(args, hitframe)

    # Only the first row survives every filter.
    assert list(filtered.index) == [0]


def test_filter_guides_keeps_everything_when_no_patterns_are_given():
    hitframe = pd.DataFrame(
        {
            "pam": ["AGG", "TGG"],
            "sgrna": ["AAA", "GGGG"],
            "downstream": ["TTT", "CCCC"],
            "gc": [0.5, 0.5],
            "off_target_count": [0, 1],
        }
    )
    args = SgRNAargs(Dseqrecord("ATGCGG"), ["TEST_LOCUS"], gc_lower=0.0, gc_upper=1.0)

    filtered = filter_guides(args, hitframe)

    assert list(filtered.index) == [0, 1]


def test_extract_sgrnas_predict_step_is_a_placeholder(capsys):
    genome = _single_cds_genome(
        "ATG" + "A" * 40 + "CC" + "G" * 40, "TEST_LOCUS_014", 1, "cas9_predict"
    )

    df = extract_sgRNAs(
        SgRNAargs(genome, ["TEST_LOCUS_014"], cas_type="cas9", step=["find", "predict"])
    )

    assert not df.empty
    assert "We will try to implement this part later" in capsys.readouterr().out


def test_find_sgrna_hits_cas3_guide_on_the_reverse_complement_scan():
    # The CDS itself holds no TTC, but its reverse complement starts with one,
    # so the guide is found while scanning the opposite strand.
    genome = _single_cds_genome("T" * 34 + "GAA", "TEST_LOCUS_014", 1, "cas3_revcomp")

    df = find_sgrna_hits_cas3(genome, genome.id, ["TEST_LOCUS_014"], Counter(), 13, revcomp)

    assert len(df) == 1
    row = df.iloc[0]
    assert row["pam"] == "TTC"
    assert row["sgrna"] == "A" * 34
    assert row["gene_strand"] == 1
    # found on the reverse complement of a plus-strand gene
    assert row["sgrna_strand"] == -1
    # len(coding seq) - match.start() - pam_len == 37 - 0 - 3
    assert row["sgrna_loc"] == 34
