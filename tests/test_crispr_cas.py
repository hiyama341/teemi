#!/usr/bin/env python

from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from pydna.dseqrecord import Dseqrecord

from teemi.design.crispr_cas import SgRNAargs, extract_sgRNAs


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
