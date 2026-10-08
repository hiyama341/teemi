#!/usr/bin/env python

# Tests for the gibson_cloning module
import random

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer
from pydna.tm import tm_default

import teemi.design.gibson_cloning as gibson_cloning
from teemi.design.gibson_cloning import (
    assemble_multiple_plasmids_with_repair_templates_for_deletion,
    assemble_single_plasmid_with_repair_templates,
    extract_locus_tag_homology_arms,
    find_up_dw_repair_templates,
    update_primer_names,
)


def random_dna(length, seed):
    """Deterministic pseudo random DNA used to build small synthetic genomes."""
    rng = random.Random(seed)
    return "".join(rng.choice("ATGC") for _ in range(length))


GENOME_SEQ = random_dna(4000, 7)
VECTOR_SEQ = random_dna(1200, 99)


def make_genome():
    """A 4 kb genome with one CDS (SCO_0001) spanning 1500..1800."""
    cds = SeqFeature(FeatureLocation(1500, 1800, strand=1), type="CDS")
    cds.qualifiers["locus_tag"] = ["SCO_0001"]
    other_cds = SeqFeature(FeatureLocation(2500, 2700, strand=-1), type="CDS")
    other_cds.qualifiers["locus_tag"] = ["SCO_9999"]
    gene = SeqFeature(FeatureLocation(1500, 1800, strand=1), type="gene")
    gene.qualifiers["locus_tag"] = ["SCO_0001"]
    return SeqRecord(
        Seq(GENOME_SEQ),
        id="synthetic_chromosome",
        name="synthetic_chromosome",
        features=[gene, cds, other_cds],
    )


def test_find_up_dw_repair_templates():
    genome = make_genome()

    records = find_up_dw_repair_templates(
        genome, ["SCO_0001"], target_tm=55, primer_calc_function=tm_default
    )

    # only the CDS feature with a matching locus tag gives a record
    assert len(records) == 1
    record = records[0]
    assert record["name"] == "SCO_0001"

    # the templates are the 1000 bp directly up- and downstream of the CDS
    assert record["location_up_start"] == 500
    assert record["location_up_end"] == 1500
    assert record["location_dw_start"] == 1800
    assert record["location_dw_end"] == 2800
    assert str(record["up_repair"].seq) == GENOME_SEQ[500:1500]
    assert str(record["dw_repair"].seq) == GENOME_SEQ[1800:2800]
    assert record["up_repair"].template.name == "Repair_Template_UPSTREAMSCO_0001"
    assert record["dw_repair"].template.name == "Repair_Template_DownstreamSCO_0001"

    # primers anneal at the very ends of their template
    up_fw = str(record["up_forwar_p"].seq)
    up_rv = str(record["up_reverse_p"].seq)
    assert GENOME_SEQ[500:1500].startswith(up_fw)
    assert str(Seq(up_rv).reverse_complement()) == GENOME_SEQ[1500 - len(up_rv) : 1500]

    dw_fw = str(record["dw_forwar_p"].seq)
    dw_rv = str(record["dw_reverse_p"].seq)
    assert GENOME_SEQ[1800:2800].startswith(dw_fw)
    assert str(Seq(dw_rv).reverse_complement()) == GENOME_SEQ[2800 - len(dw_rv) : 2800]

    # the reported Tms come from the supplied Tm function
    assert record["tm_up_forwar_p"] == pytest.approx(tm_default(up_fw))
    assert record["tm_up_reverse_p"] == pytest.approx(tm_default(up_rv))
    assert record["tm_dw_forwar_p"] == pytest.approx(tm_default(dw_fw))
    assert record["tm_dw_reverse_p"] == pytest.approx(tm_default(dw_rv))

    # the primers are close to the requested target Tm
    for tm_key in (
        "tm_up_forwar_p",
        "tm_up_reverse_p",
        "tm_dw_forwar_p",
        "tm_dw_reverse_p",
    ):
        assert 52 <= record[tm_key] <= 58


def test_find_up_dw_repair_templates_without_matching_locus_tag():
    genome = make_genome()

    records = find_up_dw_repair_templates(
        genome, ["NOT_IN_GENOME"], target_tm=55, primer_calc_function=tm_default
    )

    assert records == []


def test_extract_locus_tag_homology_arms_skips_other_features():
    cds = SeqFeature(FeatureLocation(100, 160, strand=-1), type="CDS")
    cds.qualifiers["locus_tag"] = ["WANTED"]
    other_cds = SeqFeature(FeatureLocation(200, 260, strand=1), type="CDS")
    other_cds.qualifiers["locus_tag"] = ["NOT_WANTED"]
    # a non CDS feature carrying the wanted locus tag must be ignored
    gene = SeqFeature(FeatureLocation(100, 160, strand=-1), type="gene")
    gene.qualifiers["locus_tag"] = ["WANTED"]
    # a CDS without a locus_tag qualifier must not raise
    unnamed_cds = SeqFeature(FeatureLocation(10, 20, strand=1), type="CDS")

    genome = SeqRecord(
        Seq(random_dna(300, 3)),
        id="chrom1",
        name="chrom1",
        features=[gene, cds, other_cds, unnamed_cds],
    )

    arms = extract_locus_tag_homology_arms(genome, ["WANTED"], arm_length=20)

    assert len(arms) == 1
    assert arms.loc[0, "locus_tag"] == "WANTED"
    assert arms.loc[0, "gene_strand"] == -1
    assert arms.loc[0, "upstream_arm"] == str(genome.seq[80:100])
    assert arms.loc[0, "downstream_arm"] == str(genome.seq[160:180])
    assert arms.loc[0, "repair_oligo"] == str(genome.seq[80:100]) + str(
        genome.seq[160:180]
    )
    assert arms.loc[0, "repair_oligo_length"] == 40


def test_extract_locus_tag_homology_arms_clamps_to_genome_ends():
    first_cds = SeqFeature(FeatureLocation(5, 20, strand=1), type="CDS")
    first_cds.qualifiers["locus_tag"] = ["B_EDGE"]
    last_cds = SeqFeature(FeatureLocation(30, 45, strand=1), type="CDS")
    last_cds.qualifiers["locus_tag"] = ["A_EDGE"]
    genome = SeqRecord(
        Seq(random_dna(50, 4)),
        id="small",
        name="small",
        features=[first_cds, last_cds],
    )

    arms = extract_locus_tag_homology_arms(genome, ["A_EDGE", "B_EDGE"], arm_length=45)

    # rows are sorted by locus tag, not by genome position
    assert list(arms["locus_tag"]) == ["A_EDGE", "B_EDGE"]
    assert list(arms.index) == [0, 1]

    # the arms are clamped at the start and the end of the genome
    b_edge = arms[arms["locus_tag"] == "B_EDGE"].iloc[0]
    assert b_edge["upstream_arm_start"] == 0
    assert b_edge["upstream_arm_length"] == 5
    a_edge = arms[arms["locus_tag"] == "A_EDGE"].iloc[0]
    assert a_edge["downstream_arm_end"] == 50
    assert a_edge["downstream_arm_length"] == 5


def test_extract_locus_tag_homology_arms_returns_empty_dataframe():
    genome = SeqRecord(Seq(random_dna(100, 5)), id="chrom1", name="chrom1")

    arms = extract_locus_tag_homology_arms(genome, ["MISSING"])

    assert isinstance(arms, pd.DataFrame)
    assert arms.empty


def make_primer_record(gene_name, offset):
    """Four distinct primers wrapped in the dict layout update_primer_names expects."""
    return {
        "gene_name": gene_name,
        "up_forwar_p": Primer(random_dna(20, 100 + offset), id="f", name="f"),
        "up_reverse_p": Primer(random_dna(20, 200 + offset), id="r", name="r"),
        "dw_forwar_p": Primer(random_dna(20, 300 + offset), id="f", name="f"),
        "dw_reverse_p": Primer(random_dna(20, 400 + offset), id="r", name="r"),
    }


def test_update_primer_names_renames_unique_primers():
    records = [make_primer_record("geneA", 0), make_primer_record("geneB", 1)]

    update_primer_names(records)

    assert records[0]["up_forwar_p_name"] == "geneA_up_F0"
    assert records[0]["up_reverse_p_name"] == "geneA_up_R0"
    assert records[0]["dw_forwar_p_name"] == "geneA_dw_F0"
    assert records[0]["dw_reverse_p_name"] == "geneA_dw_R0"
    assert records[1]["up_forwar_p_name"] == "geneB_up_F1"
    assert records[1]["up_reverse_p_name"] == "geneB_up_R1"
    assert records[1]["dw_forwar_p_name"] == "geneB_dw_F1"
    assert records[1]["dw_reverse_p_name"] == "geneB_dw_R1"

    # the primer objects themselves are renamed in place
    assert records[0]["up_forwar_p"].name == "geneA_up_F0"
    assert records[0]["up_forwar_p"].id == "geneA_up_F0"
    assert records[1]["dw_reverse_p"].name == "geneB_dw_R1"
    assert records[1]["dw_reverse_p"].id == "geneB_dw_R1"


def test_update_primer_names_reuses_name_of_identical_primer():
    first = make_primer_record("geneA", 0)
    # the second record shares all four primer sequences with the first one
    second = {
        "gene_name": "geneB",
        "up_forwar_p": Primer(str(first["up_forwar_p"].seq), id="f", name="f"),
        "up_reverse_p": Primer(str(first["up_reverse_p"].seq), id="r", name="r"),
        "dw_forwar_p": Primer(str(first["dw_forwar_p"].seq), id="f", name="f"),
        "dw_reverse_p": Primer(str(first["dw_reverse_p"].seq), id="r", name="r"),
    }
    records = [first, second]

    update_primer_names(records)

    # the duplicates are not renamed, they point at the already ordered primer
    assert second["up_forwar_p_name"] == "geneA_up_F0"
    assert second["up_reverse_p_name"] == "geneA_up_R0"
    assert second["dw_forwar_p_name"] == "geneA_dw_F0"
    assert second["dw_reverse_p_name"] == "geneA_dw_R0"
    assert second["up_forwar_p"].id == "f"
    assert second["dw_reverse_p"].id == "r"


def overlap_length(left, right, max_overlap=80):
    """Length of the longest suffix of ``left`` that is a prefix of ``right``."""
    for length in range(min(max_overlap, len(left), len(right)), 0, -1):
        if left[-length:] == right[:length]:
            return length
    return 0


def test_assemble_single_plasmid_with_repair_templates():
    genome = make_genome()
    records = find_up_dw_repair_templates(
        genome, ["SCO_0001"], target_tm=55, primer_calc_function=tm_default
    )
    repair_templates = [records[0]["up_repair"], records[0]["dw_repair"]]
    vector = Dseqrecord(VECTOR_SEQ, name="pDEL_SCO_0001")

    parts = assemble_single_plasmid_with_repair_templates(
        repair_templates, vector, overlap=35
    )

    # the vector is added both in front of and behind the repair templates
    assert len(parts) == 4
    assert str(parts[0].seq) == VECTOR_SEQ
    assert str(parts[3].seq) == VECTOR_SEQ

    # the two amplicons grew beyond their 1000 bp template because of the tails
    up_amplicon, dw_amplicon = parts[1], parts[2]
    assert len(up_amplicon) > 1000
    assert len(dw_amplicon) > 1000
    assert str(up_amplicon.seq).startswith(VECTOR_SEQ[-35:])
    assert str(dw_amplicon.seq).endswith(VECTOR_SEQ[:35])

    # every junction shares at least the requested overlap
    assert overlap_length(VECTOR_SEQ, str(up_amplicon.seq)) >= 35
    assert overlap_length(str(up_amplicon.seq), str(dw_amplicon.seq)) >= 35
    assert overlap_length(str(dw_amplicon.seq), VECTOR_SEQ) >= 35

    assert str(up_amplicon.forward_primer.footprint) in GENOME_SEQ[500:1500]
    assert str(dw_amplicon.forward_primer.footprint) in GENOME_SEQ[1800:2800]


def test_assemble_multiple_plasmids_with_repair_templates_for_deletion(monkeypatch):
    # keep the Tm calculation offline and deterministic
    monkeypatch.setattr(
        gibson_cloning, "primer_tm_neb", lambda seq: round(tm_default(str(seq)), 2)
    )

    genome = make_genome()
    repair_dna_templates = find_up_dw_repair_templates(
        genome, ["SCO_0001"], target_tm=55, primer_calc_function=tm_default
    )
    vector = Dseqrecord(VECTOR_SEQ, name="pDEL_SCO_0001")

    records = assemble_multiple_plasmids_with_repair_templates_for_deletion(
        ["SCO_0001"], [vector], repair_dna_templates, overlap=35
    )

    assert len(records) == 1
    record = records[0]
    assert record["gene_name"] == "SCO_0001"
    assert record["name"] == "pDEL_SCO_0001"

    # vector + both 1000 bp repair templates, circularised
    contig = record["contig"]
    assert len(contig) == 1200 + 1000 + 1000
    assert contig.circular is True

    # the two CDS features mark the repair templates behind the vector
    cds_features = [f for f in contig.features if f.type == "CDS"]
    assert len(cds_features) == 2
    assert cds_features[0].qualifiers["label"] == "UP_repairSCO_0001"
    assert int(cds_features[0].location.start) == 1200
    assert int(cds_features[0].location.end) == 2200
    assert cds_features[1].qualifiers["label"] == "DW_repairSCO_0001"
    assert int(cds_features[1].location.start) == 2200
    assert int(cds_features[1].location.end) == 3200

    # primer footprints anneal on the genome, the full primers carry Gibson tails
    assert str(record["up_forwar_p_anneal"]) in GENOME_SEQ[500:1500]
    assert str(record["dw_forwar_p_anneal"]) in GENOME_SEQ[1800:2800]
    assert str(record["up_forwar_p"].seq).endswith(str(record["up_forwar_p_anneal"]))
    assert len(record["up_forwar_p"]) > len(record["up_forwar_p_anneal"])
    assert record["up_forwar_p_name"] == record["up_forwar_p"].id
    assert record["dw_reverse_p_name"] == record["dw_reverse_p"].id

    # Tms are computed on the annealing part only
    assert record["tm_up_forwar_p"] == pytest.approx(
        round(tm_default(str(record["up_forwar_p_anneal"])), 2)
    )
    assert record["tm_dw_reverse_p"] == pytest.approx(
        round(tm_default(str(record["dw_reverse_p_anneal"])), 2)
    )


def test_assemble_multiple_plasmids_skips_unmatched_input(monkeypatch):
    monkeypatch.setattr(
        gibson_cloning, "primer_tm_neb", lambda seq: round(tm_default(str(seq)), 2)
    )

    genome = make_genome()
    repair_dna_templates = find_up_dw_repair_templates(
        genome, ["SCO_0001"], target_tm=55, primer_calc_function=tm_default
    )

    # the plasmid name does not contain the gene name
    unrelated_vector = Dseqrecord(VECTOR_SEQ, name="pUNRELATED")
    assert (
        assemble_multiple_plasmids_with_repair_templates_for_deletion(
            ["SCO_0001"], [unrelated_vector], repair_dna_templates
        )
        == []
    )

    # the plasmid matches but there is no repair template for that gene
    matching_vector = Dseqrecord(VECTOR_SEQ, name="pDEL_SCO_0002")
    assert (
        assemble_multiple_plasmids_with_repair_templates_for_deletion(
            ["SCO_0002"], [matching_vector], repair_dna_templates
        )
        == []
    )
