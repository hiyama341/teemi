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




    





     
