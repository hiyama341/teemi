#!/usr/bin/env python

# Test PCR module

from pydna.dseqrecord import Dseqrecord
from pydna.amplify import pcr
from pydna.primer import Primer
from Bio.SeqRecord import SeqRecord


# Importing the module we are  testing
from teemi.build.PCR import *

def test_amplicon_by_name(): 
    
    # initialize
    middle = 'a'*2000
    template = Dseqrecord("tacactcaccgtctatcattatcagcgacgaagcgagcgcgaccgcgagcgcgagcgca"+middle+"caggagcgagacacggcgacgcgagcgagcgagcgatactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatcattatc")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon = pcr(p1, p2, template)
    amplicon.name = 'AMPICON_FOR_TESTING_amplicon_byname_function'

    amplicon_list = [amplicon]
    my_amplicon = amplicon_by_name('AMPICON_FOR_TESTING_amplicon_byname_function', amplicon_list)

    assert my_amplicon.name == 'AMPICON_FOR_TESTING_amplicon_byname_function'
    assert len(my_amplicon) == 2123
    
def test_calculate_processing_speed(): 
    
    # initialize
    template = Dseqrecord("tacactcaccgtctatcattatctactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatcattatc")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon = pcr(p1, p2, template)
    
    # the tests
    amplicon.annotations['polymerase'] = "OneTaq Hot Start"
    calculate_processing_speed(amplicon)
    assert amplicon.annotations['proc_speed'] == 60
    
    amplicon.annotations['polymerase'] = "Q5 Hot Start"
    calculate_processing_speed(amplicon)
    assert amplicon.annotations['proc_speed'] == 30
    

    amplicon.annotations['polymerase'] = "Phusion"
    calculate_processing_speed(amplicon)
    assert amplicon.annotations['proc_speed'] == 30
    

def test_calculate_elongation_time():
    
    # initialize
    template = Dseqrecord("tacactcaccgtctatcattatctactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatcattatc")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon = pcr(p1, p2, template)
    amplicon.annotations['polymerase'] = "OneTaq Hot Start"
    calculate_processing_speed(amplicon)
    #the test
    calculate_elongation_time(amplicon)
    
    assert amplicon.annotations['elongation_time'] == 4
    
    
    # initialize 2
    middle = 'a'*2000
    template = Dseqrecord("tacactcaccgtctatcattatcagcgacgaagcgagcgcgaccgcgagcgcgagcgca"+middle+"caggagcgagacacggcgacgcgagcgagcgagcgatactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatcattatc")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon = pcr(p1, p2, template)
    amplicon.annotations['polymerase'] = "OneTaq Hot Start"
    calculate_processing_speed(amplicon)
    calculate_elongation_time(amplicon)

    #tests
    assert amplicon.annotations['elongation_time'] == 128
    
    
    # initialize 2
    middle = 'a'*3000
    template = Dseqrecord("tacactcaccgtctatcattatcagcgacgaagcgagcgcgaccgcgagcgcgagcgca"+middle+"caggagcgagacacggcgacgcgagcgagcgagcgatactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatcattatc")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon = pcr(p1, p2, template)
    amplicon.annotations['polymerase'] = "OneTaq Hot Start"
    calculate_processing_speed(amplicon)
    calculate_elongation_time(amplicon)

    #test
    assert amplicon.annotations['elongation_time'] == 188

    

def test_PCR_program(): 
    # intialize an amplicon 
    middle = 'a'*2000
    template = Dseqrecord("tacactcaccgtctatcattatcagcgacgaagcgagcgcgaccgcgagcgcgagcgca"+middle+"caggagcgagacacggcgacgcgagcgagcgagcgatactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatcattatc")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon = pcr(p1, p2, template)
    amplicon.name = 'AMPICON_FOR_TESTING_PCR_program_function'
    amplicon.annotations['polymerase'] = "OneTaq Hot Start"
    calculate_processing_speed(amplicon)

    # initialize a string object
    program = Q5_NEB_PCR_program(amplicon)

    assert program[1:5] == '98°C'
    assert program[35:39] == '61.0'
    assert program[61:65] == '72°C'


def test_grouper():
    elong_times = [60,60, 46, 60, 45, 30, 200, 100]
    elong_times.sort()
    elong_time_max_diff = 10
    groups = dict(enumerate(grouper(elong_times,elong_time_max_diff), 1))

    assert groups == {1: [30], 2: [45, 46], 3: [60, 60, 60], 4: [100], 5: [200]}



def test_calculate_required_thermal_cyclers():
    # amplicon1
    middle = 'c'*100
    template = Dseqrecord("tacactcaccgtctatcattatcagcgacgaagcgagcgcgaccgcgagcgcgagcgca"+middle+"caggagcgagacacggcgacgcgagcgagcgagcgatactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatca")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon1 = pcr(p1, p2, template)
    amplicon1.name = 'AMPLICON1'
    amplicon1.annotations['polymerase'] = "OneTaq Hot Start"
    calculate_processing_speed(amplicon1)
    calculate_elongation_time(amplicon1)

    # amplicon2
    middle = 'a'*200
    template = Dseqrecord("tacactcaccgtctatcattatcagcgacgaagcgagcgcgaccgcgagcgcgagcgca"+middle+"caggagcgagacacggcgacgcgagcgagcgagcgatactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatca")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon2 = pcr(p1, p2, template)
    amplicon2.name = 'AMPLICON2'
    amplicon2.annotations['polymerase'] = "OneTaq Hot Start"
    calculate_processing_speed(amplicon2)
    calculate_elongation_time(amplicon2)

    # adding extra annotations
    amplicon1_program = Q5_NEB_PCR_program(amplicon1)
    amplicon2_program = Q5_NEB_PCR_program(amplicon2)

    amplicons = [amplicon1, amplicon2]
    # running the function
    thermal_cyclers = calculate_required_thermal_cyclers(amplicons, polymerase='Q5 Hot Start') # pol

    assert thermal_cyclers.iloc[0]['tas'] == 59
    assert thermal_cyclers.iloc[0]['elong_times'] == 20
    assert thermal_cyclers.iloc[0]['amplicons'] == 'AMPLICON1, AMPLICON2'


def test_pcr_locations():
    from pydna.dseqrecord import Dseqrecord
    from pydna.amplify import pcr


    dna = Dseqrecord('ATGATATATGGCTCGACTGCAGGGGGATTTTTCCGGATCGCGGTCGATGACTGATACTACTACGACTACTAG')
    primer1 = Dseqrecord('ATGATATATGGCTCGAC')
    primer2 = Dseqrecord('TACTACGACTACTAG').reverse_complement()

    #names 
    dna.name = 'dna'
    primer1.name = 'primer1' 
    primer2.name  = 'primer2' 

    #annott
    dna.annotations['batches'] = [{'location':'Freezer1', 'concentration':123}]
    primer1.annotations['batches']  = [{'location':'Freezer2','concentration':274}]
    primer2.annotations['batches']   = [{'location':'Freezer3','concentration':124}]

    # make a pcr_prod
    gRNA1_pcr_prod = pcr(primer1,primer2, dna)

    # run the function
    pcr_locations_df = pcr_locations([gRNA1_pcr_prod])

    assert pcr_locations_df.iloc[0]['location'] == 'Freezer1'
    assert pcr_locations_df.iloc[0]['name'] == '72bp_PCR_prod'
    assert pcr_locations_df.iloc[0]['template'] == 'Freezer1'
    assert pcr_locations_df.iloc[0]['fw'] == 'Freezer2'
    assert pcr_locations_df.iloc[0]['rv'] == 'Freezer3'


def test_nanophotometer_concentrations(): 
    list_of_conc = nanophotometer_concentrations(path = '../teemi/tests/files_for_testing/2021-03-29_G8H_CPR_library_part_concentrations.tsv')

    assert list_of_conc[0] ==142.8
    assert list_of_conc[1] ==134.5

    assert list_of_conc[-2] ==17.5
    assert list_of_conc[-1] ==39.9


def test_calculate_volumes(): 

    calculate_volumes_df = calculate_volumes(vol_p_reac = 20, 
            no_of_reactions = 3,
            standard_reagents = ["Template", "Primer 1", "Primer 2", "H20", "Pol"],
            standard_volumes = [1, 2.5, 2.5, 19, 25])
    
    # template
    assert calculate_volumes_df.iloc[0]['vol_p_reac'] == 0.4
    assert round(calculate_volumes_df.iloc[0]['vol_p_3_reac'],2) == 1.2

    # h2o
    assert calculate_volumes_df.iloc[3]['vol_p_reac'] == 7.6
    assert round(calculate_volumes_df.iloc[3]['vol_p_3_reac'],2) == 22.8

    # total
    assert calculate_volumes_df.iloc[5]['vol_p_reac'] == 20.0
    assert round(calculate_volumes_df.iloc[5]['vol_p_3_reac'],2) == 60.0