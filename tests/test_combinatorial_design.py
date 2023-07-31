#!/usr/bin/env python

# Test Combinatorial_design
import pandas as pd
from pandas.testing import assert_frame_equal

# For getting the sequences
from Bio import SeqIO
from pydna.dseqrecord import Dseqrecord
from Bio.Seq import Seq
from pydna.amplicon import Amplicon
import pytest

# Import the DesignAssembly modules
from teemi.design.combinatorial_design import (
    DesignAssembly,
    count_unique_parts,
    get_combinatorial_list,
    get_systematic_names,
    simple_amplicon_maker,
    get_primers,
    assembly_maker,
    unique_primers,
    get_assembly_figure

)
# First some data is imported
promoter = []
cds = []
terminator = []

# GEtting the objects from the fasta files
for seq_record in SeqIO.parse('../teemi/tests/files_for_testing/promoter_test_data.fasta', format= 'fasta'):
    promoter.append(seq_record)
for seq_record in SeqIO.parse('../teemi/tests/files_for_testing/CDS_test_data.fasta', format= 'fasta'):
    cds.append(seq_record)
for seq_record in SeqIO.parse('../teemi/tests/files_for_testing/terminator_test_data.fasta', format= 'fasta'):
    terminator.append(seq_record)

# Slicing for smaller dataset
test_prom = promoter[0:2]
test_cds = cds[0:2]
test_term = terminator[0:2]
pad = [Dseqrecord('')]


# Getting their sequences
test_prom_seqrecord = [Dseqrecord(seq) for seq in test_prom]
test_cds_seqrecord = [Dseqrecord(seq) for seq in test_cds]
test_term_seqrecord = [Dseqrecord(seq) for seq in test_term]

#Getting the objects names
test_prom_names = [names.name for names in test_prom]
test_cds_names = [names.name for names in test_cds]
test_term_names = [names.name for names in test_term]

# making listsoflists
list_of_seqs = [test_prom_seqrecord, test_cds_seqrecord, test_term_seqrecord]
list_of_names =[test_prom_names, test_cds_names,test_term_names]

# Initializing the design
test_assembly = DesignAssembly(list_of_seqs, pad, [2])
index0 = test_assembly.list_of_amplicons

import io
import sys

def test_show_contigs():
    # Redirect stdout to a buffer
    stdout = sys.stdout
    sys.stdout = io.StringIO()

    # Call the method
    test_assembly.show_contigs()

    # Get the output
    output = sys.stdout.getvalue()

    # Restore stdout
    sys.stdout = stdout

    # Remove empty lines from the output
    output_lines = [line for line in output.splitlines() if line.strip()]

    # Compare the first four lines of the output to the first four lines of the expected output
    expected_output = """\
Contig(1, 1, 1)
Template:  MLS1
Template:  AKF02530
Template:  ADH1
"""
    assert output_lines[:4] == expected_output.splitlines()[:4]


def test_DesignAssembly_lengths():
    # test lenght of assembly, primers, PCR
    assert len(test_assembly.show_variants_lib_df()) == int(8)
    assert len(test_assembly.primer_list_to_dataframe()) == int(16)
    assert len(test_assembly.pcr_list_to_dataframe()) == int(14)


def test_DesignAssembly__primer_print():
    ''' Test the print of primers '''

    from Bio.Seq import Seq
    test1 = ['P001', 'MLS1', Seq('TTTAATCTTTAGGGAGGG'), 54.72, 18, 32.4, 'Anneals to MLS1', Seq('TTTAATCTTTAGGGAGGG'),18 ]
    test2 = ['P002','MLS1', Seq('TTCCATTTCATTATCCATTTTCTTAATTCTTTTATGTGCTTTT'), 54.65, 43, 77.4, 'Anneals to MLS1, overlaps to 1611bp_PCR_prod', Seq('TTTCTTAATTCTTTTATGTGCTTTT'),25]
    primers = test_assembly.primer_list()
    index0 = primers[0]
    index1 = primers[1]
    assert test1== index0
    assert test2== index1

    assert len(primers) == 16


def test_primer_list_to_dataframe():
    # Call the method
    df = test_assembly.primer_list_to_dataframe()

    # Only keep the first four rows
    df = df.head(4)

    # Convert the "sequence" column to string type
    df['sequence'] = df['sequence'].apply(''.join)
    
    # Create a DataFrame for the expected output
    expected_df = pd.DataFrame({
        'id': ['P001', 'P002', 'P003', 'P004'],
        'anneals to': ['MLS1', 'MLS1', 'AKF02530', 'AKF02530'],
        'sequence': [
            ''.join(('T', 'T', 'T', 'A', 'A', 'T', 'C', 'T', 'T', 'T', 'A', 'G', 'G', 'G', 'A', 'G', 'G', 'G')),
            ''.join(('T', 'T', 'C', 'C', 'A', 'T', 'T', 'T', 'C', 'A', 'T', 'T', 'A', 'T', 'C', 'C', 'A', 'T', 'T', 'T', 'T', 'C', 'T', 'T', 'A', 'A', 'T', 'T', 'C', 'T', 'T', 'T', 'T', 'A', 'T', 'G', 'T', 'G', 'C', 'T', 'T', 'T', 'T')),
            ''.join(('C', 'A', 'T', 'A', 'A', 'A', 'A', 'G', 'A', 'A', 'T', 'T', 'A', 'A', 'G', 'A', 'A', 'A', 'A', 'T', 'G', 'G', 'A', 'T', 'A', 'A', 'T', 'G', 'A', 'A', 'A', 'T', 'G', 'G', 'A', 'A', 'A', 'C', 'T', 'A')),
            ''.join(('A', 'A', 'A', 'T', 'C', 'A', 'T', 'A', 'A', 'G', 'A', 'A', 'A', 'T', 'T', 'C', 'G', 'C', 'A', 'T', 'A', 'A', 'G', 'A', 'A', 'T', 'C', 'T', 'G', 'G', 'A', 'T', 'T', 'A', 'T', 'T', 'T', 'T', 'A', 'C', 'A', 'T', 'A', 'A'))
        ],
        'annealing temperature': [54.72, 54.65, 55.07, 54.50],
        'length': [18, 43, 40, 44],
        'price(DKK)': [32.4, 77.4, 72.0, 79.2],
        'description': [
            'Anneals to MLS1',
            'Anneals to MLS1, overlaps to 1611bp_PCR_prod',
            'Anneals to AKF02530, overlaps to MLS1',
            'Anneals to AKF02530, overlaps to 518bp_PCR_prod'
        ],
        'footprint': [
            ''.join(('T', 'T', 'T', 'A', 'A', 'T', 'C', 'T', 'T', 'T', 'A', 'G', 'G', 'G', 'A', 'G', 'G', 'G')),
            ''.join(('T', 'T', 'T', 'C', 'T', 'T', 'A', 'A', 'T', 'T', 'C', 'T', 'T', 'T', 'T', 'A', 'T', 'G', 'T', 'G', 'C', 'T', 'T', 'T', 'T')),
            ''.join(('A', 'T', 'G', 'G', 'A', 'T', 'A', 'A', 'T', 'G', 'A', 'A', 'A', 'T', 'G', 'G', 'A', 'A', 'A', 'C', 'T', 'A')),
            ''.join(('A', 'T', 'A', 'A', 'G', 'A', 'A', 'T', 'C', 'T', 'G', 'G', 'A', 'T', 'T', 'A', 'T', 'T', 'T', 'T', 'A', 'C', 'A', 'T', 'A', 'A'))
        ],
        'len_footprint': [18, 25, 22, 26]
    })

    # Convert the "sequence" column to string type in the expected DataFrame
    expected_df['sequence'] = expected_df['sequence'].apply(''.join)

    # Check that the DataFrame is as expected
    assert_frame_equal(df, expected_df)



def test_pcr_list_to_dataframe():
    # Call the method
    df = test_assembly.pcr_list_to_dataframe()

    # Create a DataFrame for the expected output
    expected_df = pd.DataFrame({
        'pcr_number': ['PCR1', 'PCR2', 'PCR3', 'PCR4', 'PCR5', 'PCR6', 'PCR7', 'PCR8', 'PCR9', 'PCR10', 'PCR11', 'PCR12', 'PCR13', 'PCR14'],
        'template': ['MLS1', 'AKF02530', 'ADH1', 'AKF02530', 'CYC1', 'MLS1', 'AAA17732', 'AAA17732', 'URA1', 'AKF02530', 'AKF02530', 'URA1', 'AAA17732', 'AAA17732'],
        'forward_primer': ['P001', 'P003', 'P005', 'P003', 'P008', 'P001', 'P011', 'P011', 'P012', 'P014', 'P014', 'P012', 'P016', 'P016'],
        'reverse_primer': ['P002', 'P004', 'P006', 'P007', 'P009', 'P010', 'P004', 'P007', 'P013', 'P004', 'P007', 'P015', 'P004', 'P007'],
        'f_tm': [54.72, 55.07, 55.15, 55.07, 56.02, 54.72, 55.15, 55.15, 54.63, 55.07, 55.07, 54.63, 55.15, 55.15],
        'r_tm': [54.65, 54.50, 53.60, 54.50, 56.37, 54.65, 54.50, 54.50, 54.85, 54.50, 54.50, 54.85, 54.50, 54.50]
    })

    # Check that the DataFrame is as expected
    assert_frame_equal(df, expected_df)

def test_show_variants_lib_df():
    # Call the method
    df = test_assembly.show_variants_lib_df()

    # Create a DataFrame for the expected output
    expected_df = pd.DataFrame({
        0: ['MLS1', 'MLS1', 'MLS1', 'MLS1', 'URA1', 'URA1', 'URA1', 'URA1'],
        1: ['AKF02530', 'AKF02530', 'AAA17732', 'AAA17732', 'AKF02530', 'AKF02530', 'AAA17732', 'AAA17732'],
        2: ['ADH1', 'CYC1', 'ADH1', 'CYC1', 'ADH1', 'CYC1', 'ADH1', 'CYC1'],
        'Systematic_name': [(1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2), (2, 1, 1), (2, 1, 2), (2, 2, 1), (2, 2, 2)],
        'Variant': [0, 1, 2, 3, 4, 5, 6, 7]
    })

    # Check that the DataFrame is as expected
    assert_frame_equal(df, expected_df)

def test_DesignAssembly_Combinatorial_correct_names():
    # test correct names
    assert test_assembly.combinatorial_list_of_names[0] == ('MLS1', 'AKF02530', 'ADH1')
    assert test_assembly.combinatorial_list_of_names[7] == ('URA1','AAA17732' ,'CYC1' )

    # test the lenght is correct
    assert len(test_assembly.combinatorial_list_of_names[0]) == 3
    assert len(test_assembly.combinatorial_list_of_names[1]) == 3
    assert len(test_assembly.combinatorial_list_of_names[2]) == 3
    assert len(test_assembly.combinatorial_list_of_names[3]) == 3
    assert len(test_assembly.combinatorial_list_of_names[4]) == 3
    assert len(test_assembly.combinatorial_list_of_names[5]) == 3
    assert len(test_assembly.combinatorial_list_of_names[6]) == 3
    assert len(test_assembly.combinatorial_list_of_names[7]) == 3
    # for all of them
    assert len(test_assembly.combinatorial_list_of_names) == 8


def test_DesignAssembly_correct_amplicons():
    assert len(test_assembly.list_of_amplicons) == 3
    assert len(test_assembly.list_of_amplicon_primers) == 3
    assert len(test_assembly.list_of_amplicon_primer_temps) == 3

def test_DesignAssembly_combinatorial_lenght():
    ''' To test that the contig lenght remain the same when the the amplicons get tail on'''

    # Native amplicons before their get tails on
    first_prom = test_assembly.list_of_amplicons[0][0]
    first_cds = test_assembly.list_of_amplicons[1][0]
    first_term = test_assembly.list_of_amplicons[2][0]

    # amplicons - These have ovellaps to the next part
    Amp1 = test_assembly.list_of_assemblies[0][0]
    Amp2 = test_assembly.list_of_assemblies[0][1]
    Amp3 = test_assembly.list_of_assemblies[0][2]

    # Assmebling the amplicons to 1 contig
    from pydna.assembly import Assembly
    assemblyobj = Assembly([Amp1,Amp2,Amp3])
    assembled = assemblyobj.assemble_linear()
    
    assert len(assembled[0].seq) == len(first_prom) + len(first_cds)+ len(first_term)


def test_get_combinatorial_list():
    # Define the input list of lists
    input_list = [[1, 2], ['a', 'b']]

    # Call the function
    result = get_combinatorial_list(input_list)

    # Define the expected output
    expected_output = [(1, 'a'), (1, 'b'), (2, 'a'), (2, 'b')]

    # Assert that the function output is as expected
    assert result == expected_output



def test_get_systematic_names():
    # Define the input list of lists
    input_list = [[1, 2], ['a', 'b'], [True, False, None]]

    # Call the function
    result = get_systematic_names(input_list)

    # Define the expected output
    expected_output = [(1, 1, 1), (1, 1, 2), (1, 1, 3), (1, 2, 1), (1, 2, 2), (1, 2, 3), (2, 1, 1), (2, 1, 2), (2, 1, 3), (2, 2, 1), (2, 2, 2), (2, 2, 3)]

    # Assert that the function output is as expected
    assert result == expected_output


def test_simple_amplicon_maker():

    amplicons , amplicon_primers, amplicon_primer_temps = simple_amplicon_maker(test_assembly.list_of_seqs, test_assembly.list_of_names, target_tm=55, limit=10)

    # Define the expected output
    expected_amplicons = [
        [len(amplicons[0][0]), len(amplicons[0][1])],
        [len(amplicons[1][0]),len(amplicons[1][1])],
        [len(amplicons[2][0]),len(amplicons[2][1])]
    ]
    expected_amplicon_primers = [
        [(Seq('TTTAATCTTTAGGGAGGG'), Seq('TTTCTTAATTCTTTTATGTGCTTTT')), (Seq('GTTGTATTAATTTTCTCGAAGG'), Seq('GTTTGGTACGGAAGTTC'))],
        [(Seq('ATGGATAATGAAATGGAAACTA'), Seq('ATAAGAATCTGGATTATTTTACATAA')), (Seq('ATGGATATGGAAATGTAAACTA'), Seq('ATAAGAATCTGGATTATTTTACATAA'))],
        [(Seq('GCGAATTTCTTATGATTTATGATTT'), Seq('CGTAAAAAAAGCATGCAC')), (Seq('ACAGGCCCCTTTTC'), Seq('GTCGACAACTAAACTGGAA'))]
    ]
    expected_amplicon_primer_temps = [
        [(54.723946785693045, 54.645978169312286), (54.63044764231597, 54.850554260067895)],
        [(55.07036270899346, 54.496272626779955), (55.14528208347963, 54.496272626779955)],
        [(55.148984180455045, 53.59627752094883), (56.01845012485285, 56.37058397752662)]
    ]

    # Assert that the output is as expected
    for i in range(len(amplicons)):
        for j in range(len(amplicons[i])):
            assert len(amplicons[i][j]) == expected_amplicons[i][j]
            assert str(amplicon_primers[i][j][0]) == str(expected_amplicon_primers[i][j][0])
            assert str(amplicon_primers[i][j][1]) == str(expected_amplicon_primers[i][j][1])
            assert amplicon_primer_temps[i][j][0] == pytest.approx(expected_amplicon_primer_temps[i][j][0], 0.01)
            assert amplicon_primer_temps[i][j][1] == pytest.approx(expected_amplicon_primer_temps[i][j][1], 0.01)

def test_get_primers():

    primers = get_primers(test_assembly.list_of_assemblies, test_assembly.combinatorial_list_of_names, test_assembly.combinatorial_list_of_primer_tm)

    # there are 8 constructs
    assert len(primers) == 8
    # three parts per construct and thus 3 primer pairs
    assert len(primers[0]) == 3
    # each is a pair
    assert len(primers[0][0]) == 2

    # sequence of those pairs
    assert primers[0][0][0].seq == Seq('TTTAATCTTTAGGGAGGG')
    assert primers[0][0][1].seq == Seq('TTCCATTTCATTATCCATTTTCTTAATTCTTTTATGTGCTTTT')

    # len 
    assert len(primers[0][0][0]) ==18
    assert len(primers[0][0][1]) == 43


def test_assembly_maker():
    list_of_assemblies = assembly_maker(test_assembly.combinatorial_list_of_amplicons, overlap=35)

    # Define the expected output
    expected_list_of_assemblies = [
        [len(list_of_assemblies[0][0]), len(list_of_assemblies[0][1]), len(list_of_assemblies[0][2])],
        [len(list_of_assemblies[1][0]),len(list_of_assemblies[1][1]), len(list_of_assemblies[1][2])],
        [len(list_of_assemblies[2][0]),len(list_of_assemblies[2][1]), len(list_of_assemblies[2][2])], 
        [len(list_of_assemblies[3][0]), len(list_of_assemblies[3][1]), len(list_of_assemblies[3][2])],
        [len(list_of_assemblies[4][0]),len(list_of_assemblies[4][1]), len(list_of_assemblies[4][2])],
        [len(list_of_assemblies[5][0]),len(list_of_assemblies[5][1]), len(list_of_assemblies[5][2])], 
        [len(list_of_assemblies[6][0]), len(list_of_assemblies[6][1]), len(list_of_assemblies[6][2])],
        [len(list_of_assemblies[7][0]),len(list_of_assemblies[7][1]), len(list_of_assemblies[7][2])],
    ]

    # Assert that the output is as expected
    for i in range(len(list_of_assemblies)):
        for j in range(len(list_of_assemblies[i])):
            assert len(list_of_assemblies[i][j]) == expected_list_of_assemblies[i][j]



def test_unique_primers():

    # Call the function
    primer_info = unique_primers(test_assembly.primers, test_assembly.list_of_assemblies)
    # lets just check the 4 first
    primer_info = primer_info[:4]

    # Define the expected output
    expected_primer_info = [['P001',
  'MLS1',
  Seq('TTTAATCTTTAGGGAGGG'),
  54.72,
  18,
  32.4,
  'Anneals to MLS1',
  Seq('TTTAATCTTTAGGGAGGG'),
  18],
 ['P002',
  'MLS1',
  Seq('TTCCATTTCATTATCCATTTTCTTAATTCTTTTATGTGCTTTT'),
  54.65,
  43,
  77.4,
  'Anneals to MLS1, overlaps to AKF02530',
  Seq('TTTCTTAATTCTTTTATGTGCTTTT'),
  25],
 ['P003',
  'AKF02530',
  Seq('CATAAAAGAATTAAGAAAATGGATAATGAAATGGAAACTA'),
  55.07,
  40,
  72.0,
  'Anneals to AKF02530, overlaps to MLS1',
  Seq('ATGGATAATGAAATGGAAACTA'),
  22],
 ['P004',
  'AAA17732',
  Seq('AAATCATAAGAAATTCGCATAAGAATCTGGATTATTTTACATAA'),
  54.5,
  44,
  79.2,
  'Anneals to AAA17732, overlaps to ADH1',
  Seq('ATAAGAATCTGGATTATTTTACATAA'),
  26]]

    # Assert that the output is as expected
    assert len(primer_info) == len(expected_primer_info)
    for i in range(len(primer_info)):
        assert primer_info[i][0] == expected_primer_info[i][0]
        assert primer_info[i][1] == expected_primer_info[i][1]
        assert str(primer_info[i][2]) == expected_primer_info[i][2]
        assert pytest.approx(primer_info[i][3], 0.01) == expected_primer_info[i][3]
        assert primer_info[i][4] == expected_primer_info[i][4]
        assert pytest.approx(primer_info[i][5], 0.01) == expected_primer_info[i][5]
        assert primer_info[i][6] == expected_primer_info[i][6]
        assert str(primer_info[i][7]) == expected_primer_info[i][7]
        assert primer_info[i][8] == expected_primer_info[i][8]



def test_count_unique_parts():
# Define a DataFrame for testing
    data = {
        'G8H': ['Rsep', 'Smus', 'Rsep', 'Vmin', 'Smus'], # 3 
        'pG8H': ['ENO2', 'ENO2', 'ENO2', 'ENO2', 'ENO2'], # 1
        'CPR': ['Clo', 'Cac', 'Ara', 'Ara', 'Aan'], # 3
        'pCPR': ['TPI1', 'TPI1', 'TPI1', 'TPI1', 'CCW12'] # 2
    }
    df = pd.DataFrame(data)

    # Define the expected output
    expected_output = {'G8H': ['Rsep', 'Smus', 'Vmin'],
 'pG8H': ['ENO2'],
 'pCPR': ['TPI1'],
 'CPR': ['Clo', 'Cac', 'Ara'],
 'sum_of_parts': '6',
 'prediction_number': '3'}

    # Call the function with the test DataFrame and a max_combinations of 625
    result = count_unique_parts(df,  6)


    # Assert that the function output is as expected
    assert result == expected_output

def test_get_assembly_figure():
    # Define the input
    first_10 = test_assembly.list_of_assemblies[:3]
    input_assembly = first_10[2]

    # Define the expected output
    expected_output = """MLS1|36
     \\/
     /\\
     36|AAA17732|36
                 \\/
                 /\\
                 36|ADH1"""

    # Call the function with the test input
    result = get_assembly_figure(input_assembly)

    # Assert that the function output is as expected
    assert result == expected_output