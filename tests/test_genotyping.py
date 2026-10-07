#!/usr/bin/env python

# Test genotyping module

import pandas as pd
from Bio import SeqIO

# Importing the module we are  testing
from teemi.test.genotyping import *


def test_pairwise_alignment_of_templates():
    # templates
    small_prom = []
    for seq_record in SeqIO.parse('../teemi/tests/files_for_testing/templates_for_pairwise_alignment.fasta', format= 'fasta'):
        small_prom.append(seq_record)

    # sequencing reads
    reads = []
    for seq_record in SeqIO.parse('../teemi/tests/files_for_testing/sequencing_reads.fasta', format= 'fasta'):
        reads.append(seq_record)    

    # primers
    pad_pG8H_fw = SeqIO.read('../teemi/tests/files_for_testing/pad_pG8H_fw.fasta', format = 'fasta')
    pad_pCPR_fw = SeqIO.read('../teemi/tests/files_for_testing/pad_pCPR_fw.fasta', format = 'fasta')
    primers_for_seq = [pad_pG8H_fw, pad_pCPR_fw]

    df_alignment = pairwise_alignment_of_templates(reads,small_prom, primers_for_seq)

    assert df_alignment.iloc[0]['inf_part_name'] == 'pTPI1'
    assert df_alignment.iloc[1]['inf_part_name'] == 'pTPI1'
    assert df_alignment.iloc[2]['inf_part_name'] == 'pCYC1'
    assert df_alignment.iloc[3]['inf_part_name'] == 'pCYC1'
    assert df_alignment.iloc[4]['inf_part_name'] == 'pCCW12'




    





     
def test_pairwise_alignment_of_templates_trims_ns_and_finds_primer():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    primer_seq = "GTTCCAGAGGCAAGCTTGAC"
    part1_body = "ATGACCGTTAAGGCTTCAGGATCCTTAGCCGGTATTCAAGTCGACGGAATTCCTGCAGTA"
    part2_seq = "CATCGGTAGCTTACGGAATCGATCCGTTAGAGCTACGGTTACCAGTGCATGAATCCGGTAAGCTTCGAATTCGGACT"

    # part2 is listed first, so part1 has to score strictly higher to be inferred
    templates = [
        SeqRecord(Seq(part2_seq), id="part2", name="part2", description="2"),
        SeqRecord(Seq(primer_seq + part1_body), id="part1", name="part1", description="1"),
    ]
    primers = [
        SeqRecord(Seq(primer_seq), id="seq_fw", name="seq_fw"),
        SeqRecord(Seq("AAAAAAAAAAAAAAAAAAAA"), id="unused", name="unused"),
    ]
    # read 1 starts with the sequencing primer and is an exact copy of part1 (after N removal)
    read1_core = primer_seq + part1_body[:40]
    read1 = SeqRecord(Seq("NN" + read1_core + "N"), id="read1", name="read1")
    # read 2 contains no primer and is an exact copy of a stretch of part2
    read2_core = part2_seq[10:60]
    read2 = SeqRecord(Seq(read2_core[:20] + "N" + read2_core[20:]), id="read2", name="read2")

    df = pairwise_alignment_of_templates([read1, read2], templates, primers)

    assert list(df.columns) == ["Sample-Name", "inf_part_name", "align_score", "inf_part_number"]
    assert df["Sample-Name"].tolist() == ["read1", "read2"]
    assert df["inf_part_name"].tolist() == ["part1", "part2"]
    assert df["inf_part_number"].tolist() == ["1", "2"]
    # localxx scores 1 per identical base, so an exact match scores the read length
    assert df["align_score"].tolist() == [float(len(read1_core)), float(len(read2_core))]


def _genotyping_fixture():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    primer_seq = "GTTCCAGAGGCAAGCTTGAC"
    body = "ATGACCGTTAAGGCTTCAGGATCCTTAGCCGGTATTCAAGTCGACGGAATTCCTGCAGTA"
    template = SeqRecord(Seq(primer_seq + body), id="part1", name="part1", description="1")
    primer = SeqRecord(Seq(primer_seq), id="seq_fw", name="seq_fw")
    return template, primer, primer_seq + body[:40]


def test_pairwise_alignment_of_templates_can_trim_at_primer():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    template, primer, core = _genotyping_fixture()
    upstream = "CAGGATCCTTAGCC"  # vector sequence in front of the primer site
    template.seq = upstream + template.seq
    read = SeqRecord(Seq(upstream + core), id="read", name="read")

    untrimmed = pairwise_alignment_of_templates([read], [template], [primer])
    trimmed = pairwise_alignment_of_templates([read], [template], [primer], trim_at_primer=True)

    # default: the whole read is aligned (published notebook behaviour)
    assert untrimmed["align_score"][0] == len(upstream + core)
    # trimmed: only the read from the primer onwards
    assert trimmed["align_score"][0] == len(core)


def test_pairwise_alignment_of_templates_short_reads_are_not_called():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    template, primer, core = _genotyping_fixture()
    short1 = SeqRecord(Seq("NNACGTACGTACGTACGTACGTANN"), id="short1", name="short1")
    long = SeqRecord(Seq(core), id="long", name="long")
    short2 = SeqRecord(Seq("ACGTACGT"), id="short2", name="short2")

    df = pairwise_alignment_of_templates([short1, long, short2], [template], [primer])

    assert df["Sample-Name"].tolist() == ["short1", "long", "short2"]
    assert df["align_score"].tolist() == [0.0, float(len(core)), 0.0]
    assert df["inf_part_name"].isna().tolist() == [True, False, True]
    assert df["inf_part_name"][1] == "part1"
    assert df["inf_part_number"].isna().tolist() == [True, False, True]
