#!/usr/bin/env python

# Test Combinatorial_design
import pandas as pd

# For getting the sequences
from Bio import SeqIO
from pydna.dseqrecord import Dseqrecord

# Import the DesignAssembly modules
from teemi.design.combinatorial_design import DesignAssembly, count_unique_parts

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
pad = Dseqrecord('')


# Getting their sequences
test_prom_seqrecord = [Dseqrecord(seq.seq) for seq in test_prom]
test_cds_seqrecord = [Dseqrecord(seq.seq) for seq in test_cds]
test_term_seqrecord = [Dseqrecord(seq.seq) for seq in test_term]

#Getting the objects names
test_prom_names = [names.name for names in test_prom]
test_cds_names = [names.name for names in test_cds]
test_term_names = [names.name for names in test_term]

# making listsoflists
list_of_seqs = [test_prom_seqrecord, test_cds_seqrecord, test_term_seqrecord]
list_of_names =[test_prom_names, test_cds_names,test_term_names]

# Initializing the design
test_assembly = DesignAssembly(list_of_seqs, list_of_names, pad, 2)

indec0 = test_assembly.list_of_amplicons

def test_DesignAssembly_lengths():
    # test lenght of assembly, primers, PCR
    assert len(test_assembly.ShowVariantsLibDF()) == int(8)
    assert len(test_assembly.primer_list_to_dataframe()) == int(20)
    assert len(test_assembly.PCR_list_to_dataframe()) == int(14)

# def test_DesignAssembly_df():
#     ''' To test if it makes a correct dataframe
#         Note: something goes wrong here
#     '''
#     dataframe = {0: ['MLS1','MLS1','MLS1','MLS1','URA1','URA1','URA1','URA1'],
#                 1: ['AKF02530', 'AKF02530', 'AAA17732','AAA17732','AKF02530','AKF02530', 'AAA17732', 'AAA17732' ],
#                 2: ['ADH1', 'CYC1','ADH1', 'CYC1','ADH1', 'CYC1','ADH1', 'CYC1'],
#                 'Systematic_name': [(1, 1, 1), (1, 1, 2), (1, 2, 1),(1, 2, 2),(2, 1, 1), (2, 1, 2), (2, 2, 1),(2, 2, 2)],
#                  'Variant': [np.int64(0),np.int64(1),np.int64(2),np.int64(3),np.int64(4),np.int64(5),np.int64(6),np.int64(7)]}
#
#     df = pd.DataFrame(dataframe)
#     df1 = pd.DataFrame(test_assembly.ShowVariantsLibDF())
#
#     assert df == df1

def test_DesignAssembly__primer_print():
    ''' Test the print of primers '''

    from Bio.Seq import Seq
    test = ['F001', 'MLS1', Seq('TTTAATCTTTAGGGAGGGT'), 56.91, 19, 34.2]
    primers = test_assembly.primer_list()
    index0 = primers[0]
    assert test== index0
    assert len(primers) == 20

def test_DesignAssembly__PCR_print():
    # Test pcr_print
    test = ['PCR1', 'MLS1', 'F001', 'R011', 56.91, 56.35]
    pcr_list = test_assembly.PCR_list()
    index0_pcr = pcr_list[0]

    assert test == index0_pcr
    assert len(pcr_list) == 14

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



# def test_count_unique_parts():
#     # Create a test DataFrame
#     df = pd.DataFrame({'G8H': [1, 2, 3, 4, 5],
#                       'pG8H': [6, 7, 8, 9, 10],
#                       'CPR': [11, 12, 13, 14, 15],
#                       'pCPR': [16, 17, 18, 19, 20]})
    
#     # Convert values to strings
#     df = df.astype(str)
    
#     # Test that the function returns the correct output
#     assert count_unique_parts(df, max_combinations=10) == {'G8H': ['1', '2'],
#                                                            'pG8H': ['6', '7'],
#                                                            'pCPR': ['16', '17'],
#                                                            'CPR': ['11', '12'],
#                                                            'Sum of parts': '16',
#                                                            'Predictions': '2'}


#if __name__ == "__main__":
#     pytest.main([__file__, "-vv", "-s"])
