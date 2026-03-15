#!/usr/bin/env python

# Test fetch_sequences module
from teemi.design.fetch_sequences import *
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

def test_retrieve_sequences_from_ncbi():
    # random acc_number
    acc_numbers = ['Q05001']
    # call the function: 
    assert retrieve_sequences_from_ncbi(acc_numbers, '../tests/files_for_testing/test_fetch.fasta') == None

    # check if it is the correct sequences
    Q05001_seq = SeqIO.read('../teemi/tests/files_for_testing/test_fetch.fasta', 'fasta')
    assert str(Q05001_seq.seq[:10]) == 'MDSSSEKLSP'
    assert Q05001_seq.id == 'sp|Q05001.1|NCPR_CATRO'
 
def test_read_fasta_files(): 
    sequences = read_fasta_files('../teemi/tests/files_for_testing/test_fetch.fasta')
    assert sequences[0].seq[:10] == 'MDSSSEKLSP'
    assert sequences[0].id == 'sp|Q05001.1|NCPR_CATRO'

def test_retrieve_sequences_from_PDB(): 
    acc_numbers = ['Q1PQK4']
    sequences = retrieve_sequences_from_PDB(acc_numbers)
    assert str(sequences[0][0].seq[:10]) == 'MQSTTSVKLS'
    assert str(sequences[0][0].id) == 'sp|A0A2U1LIM9|NCPR1_ARTAN'


def test_read_genbank_files():
    test_gb = read_genbank_files('../teemi/tests/files_for_testing/MIA-HA-1.gb')
    assert str(test_gb[0].seq[20:30]) == 'AGTTATATAG'
    assert test_gb[0].name == 'MIA-HA-1'


def test_concatenate_genbank_records():
    feature_a = SeqFeature(FeatureLocation(5, 10, strand=1), type='CDS')
    feature_a.qualifiers['locus_tag'] = ['GENE_A']
    feature_b = SeqFeature(FeatureLocation(2, 6, strand=-1), type='CDS')
    feature_b.qualifiers['locus_tag'] = ['GENE_B']

    record_a = SeqRecord(Seq('A' * 20), id='chrom1', name='chrom1', features=[feature_a])
    record_b = SeqRecord(Seq('T' * 15), id='chrom2', name='chrom2', features=[feature_b])

    combined_record = concatenate_genbank_records(
        [record_a, record_b],
        record_id='aspergillus_combined',
        record_name='aspergillus_combined',
    )

    assert len(combined_record.seq) == 35
    assert combined_record.id == 'aspergillus_combined'
    assert combined_record.name == 'aspergillus_combined'
    assert len(combined_record.features) == 2
    assert int(combined_record.features[0].location.start) == 5
    assert int(combined_record.features[0].location.end) == 10
    assert int(combined_record.features[1].location.start) == 22
    assert int(combined_record.features[1].location.end) == 26
    assert combined_record.features[1].qualifiers['locus_tag'] == ['GENE_B']


## TODO: I cannot make these work with github actions yet. It lacks some packages(intermine). Have to figure this out later. 
# def test_fetch_promoter(): 
#     cyc1 = fetch_promoter('CYC1')
#     assert cyc1[:20] == 'GAGGCACCAGCGTCAGCATT'

# def test_fetch_multiple_promoters(): 
#     list_of_promoters = ['YAR035C-A', 'YGR067C']
#     seqs = fetch_multiple_promoters(list_of_promoters)

#     assert str(seqs[0].seq[:10]) == 'CCCTGGTGGC'
#     assert str(seqs[1].seq[:10])== 'AGACAACCTA'
#     assert seqs[0].id == 'YAR035C-A'
#     assert seqs[1].id == 'YGR067C'
