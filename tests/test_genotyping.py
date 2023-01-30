#!/usr/bin/env python

# Test PCR module

import pandas as pd
from Bio import SeqIO

# Importing the module we are  testing
from teemi.test.genotyping import *

sequencing_plates = pd.read_csv('../teemi/tests/files_for_testing/Plate2Seq.csv',index_col=False)
df = [sequencing_plates]

#Plate2SeqFunctions
sliced = slicing_and_naming_seq_plates(df)
list_of_dfs = plat_seq_data_wrangler(sliced)
filtered_plates = plate_AvgQual(list_of_dfs)
split_df = split_df_names(filtered_plates)
all_data_frames = concatenating_list_of_dfs(split_df)

def test_slicing_and_naming_seq_plates():
    
    assert len(sliced[0]) == 81


def test_plat_seq_data_wrangler():

    assert len(list_of_dfs[0]) == 81
    assert type(list_of_dfs[0].iloc[3]['AvgQual']) == type(np.float64(0))


def test_plate_AvgQual():

    # do we actually filter on our parameters?
    true_false = (filtered_plates[0].iloc[:]['AvgQual']>=50).any()
    true_false1 = (filtered_plates[0].iloc[:]['used']>=25).any()

    assert len(filtered_plates[0]) == 72
    assert true_false == True
    assert true_false1 == True


def test_split_df_names():
    assert  split_df[0].columns[7] == 'plate'
    assert  split_df[0].columns[8] == 'well'


def test_concatenating_list_of_dfs():
    assert len(all_data_frames) == 72
    assert len(all_data_frames.columns) == 9

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




    





     
