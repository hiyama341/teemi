#!/usr/bin/env python
# Test Transformation module


# Importing the module we are  testing
from teemi.build.transformation import *
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pydna.amplify import pcr
import io
import sys
import pytest

def test_ng_to_nmol():
    ''' Testing that it calculates correctly'''
    
    #increasing ng
    assert ng_to_nmol(100, 1000) == 0.00015384615384615385
    assert ng_to_nmol(1000, 1000) == 0.0015384615384615385
    
    #increaseing basepair
    assert ng_to_nmol(100, 5000) == 3.076923076923077e-05
    assert ng_to_nmol(100, 10000) == 1.5384615384615384e-05
    
    # Non-valid input
    assert ng_to_nmol(-500, 8000) == 'non-valid_input'
    assert ng_to_nmol(500, -9000) == 'non-valid_input'
    
    # division by zero
    assert ng_to_nmol(0, 8000) == 'non-valid_input'
    assert ng_to_nmol(500, 0) == 'non-valid_input'



def test_ODtime(): 
    
    assert ODtime(1, 1) == 1.414

    assert ODtime(0, 1) == 0.0

    assert ODtime(1, 0) == 1.0
    assert ODtime(1, 10) == 32.0
    
    #non-valid input
    assert ODtime(-1, 1) == 'non-valid_input'
    assert ODtime(1, -1) == 'non-valid_input'
    
    
def test_wanted_mass(): 
    
    assert wanted_mass(10, 20) == 130000
    assert wanted_mass(10, 150) == 975000
    assert wanted_mass(10, 250) == 1625000

    
def test_wanted_volume():
    assert wanted_volume(100, 54) == 1.9
    assert wanted_volume(100, 153) == 0.7
    assert wanted_volume(143, 89) == 1.6


def test_transformation_mix(): 
    # load some files
    vector = SeqIO.read('../teemi/tests/files_for_testing/MIA-HA-1.gb', format = 'gb')
    gRNA1_pcr_prod = SeqIO.read('../teemi/tests/files_for_testing/gRNA_test1.fasta', format = 'fasta')
    gRNA2_pcr_prod = SeqIO.read('../teemi/tests/files_for_testing/gRNA_test2.fasta', format = 'fasta')
    LEU_plasmid = SeqIO.read('../teemi/tests/files_for_testing/MIA-HA-1.gb', format = 'gb')[1000:5000]


    #names 
    vector.name = 'p0056\\(pESC-LEU-ccdB-USER)'
    gRNA1_pcr_prod.name = 'ATF1' 
    gRNA2_pcr_prod.name  = 'CroCPR' 
    LEU_plasmid.name = 'LEU_plasmid'

    #annott
    vector.annotations['batches'] = [{'location':'Freezer1', 'concentration':123}]
    gRNA1_pcr_prod.annotations['batches']  = [{'location':'Freezer2','concentration':274}]
    gRNA2_pcr_prod.annotations['batches']   = [{'location':'Freezer3','concentration':124}]
    LEU_plasmid.annotations['batches']  = [{'location':'Freezer4','concentration':187}]

    # 1. Mention which reacion names you have
    reaction_names = ["insert", "n.ctr", "n.ctr", "n.ctr", "p. ctr"]

    # 2. Add reaction reaction_participants
    reaction_participants = [[vector, gRNA1_pcr_prod,gRNA2_pcr_prod], #the insert we want
                    [vector],                                        #negative control
                    [gRNA1_pcr_prod],                                #negative control
                    [gRNA2_pcr_prod],                                #negative control
                    [LEU_plasmid]]                                   #positive control

    # 2. Calculate nmol:
    nmol_vector = ng_to_nmol(ng = 15, bp = len(vector))
    nmol_gRNA = ng_to_nmol(ng = 30, bp = len(gRNA1_pcr_prod))
    nmol_pctr = ng_to_nmol(ng = 10, bp = len(LEU_plasmid))

    # 3. Add the concentrations
    wanted_concentrations = {'p0056\\(pESC-LEU-ccdB-USER)' : nmol_vector,
                    'ATF1'                        : nmol_gRNA,
                    'CroCPR'                      : nmol_gRNA,
                    'LEU_plasmid'                 : nmol_pctr}

    # 4. what media the transformants are plated on (5 transformations here)
    media = ['LB_AMP'] * 5

    # 5. initate the function
    trans_df = transformation_mix(reaction_names, reaction_participants, wanted_amounts = wanted_concentrations, water_dna_p_reac = 7, media = media)
    
    # insert 
    assert trans_df.iloc[0]['name'] == 'insert'
    assert trans_df.iloc[0]['p0056\(pESC-LEU-ccdB-USER)'] == 0.1
    assert trans_df.iloc[0]['ATF1'] == 0.1
    assert trans_df.iloc[0]['CroCPR'] == 0.1
    assert trans_df.iloc[0]['water'] == 6.7
    assert trans_df.iloc[0]['plate on'] == 'LB_AMP'

    # n.ctr 
    assert trans_df.iloc[1]['name'] == 'n.ctr'
    assert trans_df.iloc[1]['p0056\(pESC-LEU-ccdB-USER)'] == 0.1
    assert trans_df.iloc[1]['water'] == 6.9
    assert trans_df.iloc[1]['plate on'] == 'LB_AMP'

    # p.ctr 
    assert trans_df.iloc[4]['name'] == 'p. ctr'
    assert trans_df.iloc[4]['LEU_plasmid'] == 0.1
    assert trans_df.iloc[4]['water'] == 6.9
    assert trans_df.iloc[4]['plate on'] == 'LB_AMP'


def test_transformation_participants():
    # Create some sample reaction participants
    part1 = SeqRecord(Seq("ATCG"), name="part1")
    part2 = SeqRecord(Seq("GCTA"), name="part2")
    part3 = SeqRecord(Seq("CTGA"), name="part3")
    part4 = SeqRecord(Seq("TGAC"), name ="part4")
    reaction_participants = [[part1, part2], [part3, part4]]

    # Test with default values
    expected_output = {"part1": 0.0005, "part2": 0.0005, "part3": 0.0005, "part4": 0.0005}
    assert transformation_partitipants(reaction_participants) == expected_output

    # Test with sgRNA_plasmid_name provided
    expected_output = {"part1": 0.0005, "part2": 0.01, "part3": 0.0005, "part4": 0.0005}
    assert transformation_partitipants(reaction_participants, sgRNA_plasmid_name="part2", sgRNA_plasmid_conc=0.01) == expected_output

    # Test with different amnt value
    expected_output = {"part1": 0.001, "part2": 0.001, "part3": 0.001, "part4": 0.001}
    assert transformation_partitipants(reaction_participants, amnt=0.001) == expected_output
     
 



def test_pool_parts(example_input):
    amplicons, part_names, part_amounts, pool_names, pool_lengths = example_input
    
    expected_output = {
        'Part1': {
            'Amplicon1': {
                'volume_to_mix': 50.8,
                'location': 'Tube1',
                'concentration': 100
            },
            'Amplicon2': {
                'volume_to_mix': 25.4,
                'location': 'Tube2',
                'concentration': 200
            }
        },
        'Part2': {
            'Amplicon1': {
                'volume_to_mix': 101.6,
                'location': 'Tube1',
                'concentration': 100
            },
            'Amplicon2': {
                'volume_to_mix': 50.8,
                'location': 'Tube3',
                'concentration': 200
            }
        }
    }
    
    assert pool_parts(amplicons, part_names, part_amounts, pool_names, pool_lengths) == expected_output


def test_pool_parts():
    # the ones we want to pool
    pool_names = ['OeuG8H_tADH1','RsepG8H_tADH1','VminG8H_tADH1']
    part_names = ['OeuG8H_tADH1','RsepG8H_tADH1','VminG8H_tADH1']
    part_amounts = [0.002, 0.0014, 0.0021]
    pool_lengths = [1845,1851,1845]

    template1 = Dseqrecord('ATGATATATGGCTCGACTGCAGGGGGATTTTTCCTCTGCTGCTATATTGCCGGATCGCGGTCGATGACTGATACTACTACGACTACTAG')
    primer_fw1 = Dseqrecord('ATGATATATGGCTCGAC')
    primer_rv1 = Dseqrecord('TACTACGACTACTAG').reverse_complement()
    amplicon1 = pcr(primer_fw1, primer_rv1, template1)
    amplicon1.template.name = 'OeuG8H_tADH1'
    amplicon1.annotations['batches']= [{'location': 'l5_B02', 'concentration': 23}]


    template2 = Dseqrecord('ATGATATATGGCTCGACTGCAGGGGGATGGGGGGCGCGGCGCCTTTTCCGGATCGCGGTCGATGACTGATACTACTACGACTACTAG')
    primer_fw2 =Dseqrecord('ATGATATATGGCTCGAC')
    primer_rv2 = Dseqrecord('TACTACGACTACTAG').reverse_complement()
    amplicon2 = pcr(primer_fw2, primer_rv2,template2)
    amplicon2.template.name = 'RsepG8H_tADH1'
    amplicon2.annotations['batches']= [{'location': 'l5_B02', 'concentration': 55.0}]

    template3 =Dseqrecord('ATGATATATGGCTCGACTGCAGGGGGATTTTTGACTCAGCATGTCCGGATCGCGGTCGATGACTGATACTACTACGACTACTAG')
    primer_fw3 = Dseqrecord('ATGATATATGGCTCGAC')
    primer_rv3 = Dseqrecord('TACTACGACTACTAG').reverse_complement()
    amplicon3 = pcr( primer_fw3, primer_rv3,template3)
    amplicon3.template.name = 'VminG8H_tADH1'
    amplicon3.annotations['batches']= [{'location': 'l5_B02', 'concentration': 194.0}]


    amplicons = [amplicon1, amplicon2, amplicon3]

    # test the function
    dict_with_pooled_volumes = pool_parts(amplicons,part_names,part_amounts,  pool_names, pool_lengths)

    expected_output = {'OeuG8H_tADH1': {'89bp_PCR_prod': {'volume_to_mix': 104.3,
   'location': 'l5_B02',
   'concentration': 23}},
 'RsepG8H_tADH1': {'87bp_PCR_prod': {'volume_to_mix': 30.6,
   'location': 'l5_B02',
   'concentration': 55.0}},
 'VminG8H_tADH1': {'84bp_PCR_prod': {'volume_to_mix': 13.0,
   'location': 'l5_B02',
   'concentration': 194.0}}}

    assert dict_with_pooled_volumes == expected_output



def test_print_pooled_parts():
    # Define input parameters
    pooled_volumes = {
        'part1': {
            'sample1': {'volume_to_mix': 10, 'concentration': 100},
            'sample2': {'volume_to_mix': 20, 'concentration': 200},
        },
        'part2': {
            'sample3': {'volume_to_mix': 30, 'concentration': 300},
            'sample4': {'volume_to_mix': 40, 'concentration': 400},
        }
    }

    # Redirect stdout to a buffer
    stdout = io.StringIO()
    sys.stdout = stdout

    # Call function with input parameters
    print_pooled_parts(pooled_volumes)

    # Capture output
    output = stdout.getvalue().strip()

    # Test output
    expected_output = '''To be pooled together
part1
sample1 {'volume_to_mix': 10, 'concentration': 100}
sample2 {'volume_to_mix': 20, 'concentration': 200}
vol 30
calculated con 166.66666666666666 

part2
sample3 {'volume_to_mix': 30, 'concentration': 300}
sample4 {'volume_to_mix': 40, 'concentration': 400}
vol 70
calculated con 357.14285714285717'''
    assert output == expected_output


def test_time_to_inculate(capsys):
    # Test 1
    time_to_inoculate(initialOD=0.0025, td=0.4, verbose=True, transformation_time=12,target_OD =1, plot=False )
    captured = capsys.readouterr()
    assert "GOAL: to get enough cells in exponential phase for transformation" in captured.out
    assert "Assumptions:" in captured.out
    assert "Hours to target OD: 	" in captured.out
    assert "Time of inoculation: 	" in captured.out
    assert "NB: If you inoculated now, the cells will have reached the target OD by:   " in captured.out


def test_calculate_volume_and_total_concentration():
    template1 = Dseqrecord('ATGATATATGGCTCGACTGCAGGGGGATTTTTCCTCTGCTGCTATATTGCCGGATCGCGGTCGATGACTGATACTACTACGACTACTAG')
    primer_fw1 = Dseqrecord('ATGATATATGGCTCGAC')
    primer_rv1 = Dseqrecord('TACTACGACTACTAG').reverse_complement()
    amplicon1 = pcr(primer_fw1, primer_rv1, template1)
    amplicon1.name = 'OeuG8H_tADH1'
    amplicon1.annotations['batches']= [{'location': 'l5_B02', 'concentration': 23}]


    template2 = Dseqrecord('ATGATATATGGCTCGACTGCAGGGGGATGGGGGGCGCGGCGCCTTTTCCGGATCGCGGTCGATGACTGATACTACTACGACTACTAG')
    primer_fw2 =Dseqrecord('ATGATATATGGCTCGAC')
    primer_rv2 = Dseqrecord('TACTACGACTACTAG').reverse_complement()
    amplicon2 = pcr(primer_fw2, primer_rv2,template2)
    amplicon2.name = 'RsepG8H_tADH1'
    amplicon2.annotations['batches']= [{'location': 'l5_B02', 'concentration': 55.0}]

    template3 =Dseqrecord('ATGATATATGGCTCGACTGCAGGGGGATTTTTGACTCAGCATGTCCGGATCGCGGTCGATGACTGATACTACTACGACTACTAG')
    primer_fw3 = Dseqrecord('ATGATATATGGCTCGAC')
    primer_rv3 = Dseqrecord('TACTACGACTACTAG').reverse_complement()
    amplicon3 = pcr( primer_fw3, primer_rv3,template3)
    amplicon3.name = 'VminG8H_tADH1'
    amplicon3.annotations['batches']= [{'location': 'l5_B02', 'concentration': 194.0}]

    amplicons = [amplicon1, amplicon2, amplicon3]

    amplicon_parts_amounts_total = {'OeuG8H_tADH1': 100, 'RsepG8H_tADH1': 200, 'VminG8H_tADH1': 300}

    n = 2

    # Define expected outputs
    expected_volumes = [503043.4, 411272.8, 168866.0]
    expected_ngs = [11569998.200000001, 22620004.0, 32760004.0]
    expected_total_conc = 61.80862850220397

    # Call function and check outputs
    volumes, ngs, total_conc = calculate_volume_and_total_concentration(amplicons, amplicon_parts_amounts_total, n=n)

    assert volumes == expected_volumes
    assert ngs == expected_ngs
    assert total_conc == expected_total_conc

