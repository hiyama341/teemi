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
