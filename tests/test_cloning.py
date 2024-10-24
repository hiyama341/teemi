#!/usr/bin/env python

# Importing the module we are  testing
from pydna.amplify import pcr
from teemi.design.cloning import *
from pydna.dseqrecord import Dseqrecord
from pydna.amplify import pcr
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from teemi.design.fetch_sequences import read_fasta_files, read_genbank_files


def test_USER_enzyme(): 
    # inititalize
    template = 'TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT'
    primerF = 'CGTGCGAUTCTTTGAAAAGATAATGTATGA'
    primerR = 'ACCTGCACUTAACTAATTACATGACTCGA'

    gRNA1_pcr_prod = pcr(primerF,primerR, template)

    Digested_USER = USER_enzyme(gRNA1_pcr_prod)
    assert len(Digested_USER.seq.watson) == 417
    assert Digested_USER.seq.watson == 'TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT'


def test_remove_features_with_negative_loc():
    template = 'TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT'
    primerF = 'CGTGCGAUTCTTTGAAAAGATAATGTATGA'
    primerR = 'ACCTGCACUTAACTAATTACATGACTCGA'
    rec_vector = pcr(primerF,primerR, template)
    rec_vector.add_feature(-40,60, label =['gRNA'] )
    assert len(rec_vector.features) == 3

    # remove negative featureS
    remove_features_with_negative_loc(rec_vector)
    assert len(rec_vector.features) == 2


def test_CAS9_cutting():

    template = Dseqrecord('TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT')
    gRNA = Dseqrecord('TCTAGATTTTGTAGTGCCCT')
    up, dw = CAS9_cutting(gRNA, template)

    assert len(up)== 125
    assert up.seq.watson == 'TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGC'

    assert len(dw) == 292
    assert dw.seq.watson == 'CCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT'


# def test_CRIPSR_knockout():
#     # initialize
#     insertion_site = Dseqrecord('TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT')
#     gRNA = Dseqrecord('TCTAGATTTTGTAGTGCCCT')
#     repair_template = Dseqrecord('CGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGT')

#     # call the function
#     Knock_out = CRIPSR_knockout(gRNA,insertion_site, repair_template)
#     assert Knock_out.seq.watson == 'TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT'

def test_extract_gRNAs():
    template = 'TCTTTGAAAAGATAATGTATGATTATGCTTTCACTCATATTTATACAGAAACTTGATGTTTTCTTTCGAGTATATACAAGGTGATTACATGTACGTTTGAAGTACAACTCTAGATTTTGTAGTGCCCTCTTGGGCTAGCGGTAAAGGTGCGCATTTTTTCACACCCTACAATGTTCTGTTCAAAAGATTTTGGTCAAACGCTGTAGAAGTGAAAGTTGGTGCGCATGTTTCGGCGTTCGAAACTTCTCCGCAGTGAAAGATAAATGATCGCCGTAGTAACGTCGCTGTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGGTGCTTTTTTTGTTTTTTATGTCTTCGAGTCATGTAATTAGTTAAGTGCAGGT'
    primerF = 'CGTGCGAUTCTTTGAAAAGATAATGTATGA'
    primerR = 'ACCTGCACUTAACTAATTACATGACTCGA'

    rec_vector = pcr(primerF,primerR, template)

    # Adding a feature
    rec_vector.add_feature(40,60, name ='gRNA' , label = ['gRNA'])
    gRNA = extract_gRNAs(rec_vector, 'gRNA')

    assert gRNA[0].seq.watson == 'ACTCATATTTATACAGAAAC'


def test_extract_template_amplification_sites():
    # Initializing
    features1 = [SeqFeature(FeatureLocation(1, 100, strand=1), type='CDS')]
    features2 = [SeqFeature(FeatureLocation(101, 300, strand=1), type='terminator')]
    dict1 = {'name': 'ATF'}
    dict2 = {'name': 'ATF_tTEF'}
    features1[0].qualifiers = dict1
    features2[0].qualifiers = dict2

    TEMPLATES = [SeqRecord(seq = Seq('ATGATGA'*1000), id='seq_DPrPIMvy', name='ATF', description='<unknown description>', dbxrefs=[], features = features1+ features2)]

    # the test
    names_of_sites_we_want_to_amplify = ['ATF']
    name_of_terminator_site_to_be_incorporated= 'tTEF'

    extractions_sites = extract_template_amplification_sites(TEMPLATES, names_of_sites_we_want_to_amplify, name_of_terminator_site_to_be_incorporated)
    
    assert len(extractions_sites[0].seq) == 299
    assert len(extractions_sites[0].features) == 2
    assert str(extractions_sites[0].seq) == 'TGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATGAATGATG'
    assert extractions_sites[0].name == 'ATF'
    

def test_extract_sites(): 
    # Test 1
    features1 = [SeqFeature(FeatureLocation(2, 100, strand=1), type='promoter1')]
    features2 = [SeqFeature(FeatureLocation(200, 300, strand=1), type='promoter2')]
    features1[0].qualifiers['name'] =  ['pCYC1']
    features2[0].qualifiers['name'] =  ['pPCK']

    prom_seq = Seq('CAGCATTTTCAAAGGTGTGTTCTTCGTCAGACATGTTTTAGTGTGTGAATGAAAATATAGATAGATGATGATGATA'*50)

    prom_template = [SeqRecord(seq = prom_seq, id='seq_DPrPIMvy', name='Promoter_test', description='<unknown description>', dbxrefs=[], features = features1+ features2)]
    p_annotations = ['pCYC1']
    p_names = ['pCYC1']

    promoter_sites =  extract_sites(p_annotations, prom_template, p_names)
    assert promoter_sites[0].seq == 'GCATTTTCAAAGGTGTGTTCTTCGTCAGACATGTTTTAGTGTGTGAATGAAAATATAGATAGATGATGATGATACAGCATTTTCAAAGGTGTGTTCTT'
    assert promoter_sites[0].name == 'pCYC1'

    ##  Test 2
    test_plasmid = [SeqIO.read('../teemi/tests/files_for_testing/MIA-HA-1.gb', format = 'genbank')]
    test_plasmid

    p_annotations = ['XI-2']
    p_names = ['XI-2']

    site_on_plasmid =  extract_sites(p_annotations, test_plasmid, p_names)
    assert site_on_plasmid[0].name == 'XI-2'
    assert len(site_on_plasmid[0].seq) == 1573


def test_seq_to_annotation(): 
    test_plasmid = SeqIO.read('../teemi/tests/files_for_testing/MIA-HA-1.gb', 'gb')
    test_sequence = test_plasmid[0:100]
    
    # The actual test
    seq_to_annotation(test_sequence, test_plasmid, 'AMPLICON2')

    ## Assertions to verify the annotation
    assert test_plasmid.features[-1].type == 'AMPLICON2', "Feature type should be 'AMPLICON2'"
    
    # Access the start and end positions correctly
    assert int(test_plasmid.features[-1].location.start) == 0, "Feature start position should be 0"
    assert int(test_plasmid.features[-1].location.end) == 100, "Feature end position should be 100"


def test_casembler():
    # name 
    assembly_name = 'test'
    
    # strain and sgRNA
    HA1 = SeqIO.read('../teemi/tests/files_for_testing/MIA-HA-1.gb', 'gb')
    XI2_2_gRNA = pydna.dseqrecord.Dseqrecord("ACCCCCCTCAACTGATCAAC", name = "XI2-2_gRNA")

    # primers
    forward_primers = read_fasta_files('../teemi/tests/files_for_testing/casembler_test/MIA-HA-1_forward_primers.fasta')
    reverse_primers = read_fasta_files('../teemi/tests/files_for_testing/casembler_test/MIA-HA-1_reverse_primers.fasta')

    # amplicon sequences
    amplicons = read_genbank_files('../teemi/tests/files_for_testing/casembler_test/MIA-HA-1_casembler.gb')

    # Making amplicon objects
    pcr_amplicons = []
    for i in range( len(amplicons)): 
        pcr_amplicons.append(pydna.amplify.pcr(forward_primers[i],  reverse_primers[i],amplicons[i]))

    # adding parameters    
    parameters = {
        'bg_strain': HA1,
        'site_names': ["XI-2"],
        'gRNAs': [XI2_2_gRNA],
        'assembly_limits':[30],
        'verbose': False,
        'to_benchling': False  
        }       
    # Assembling the new strain
    assembly = [casembler(**{**parameters,
                            'parts'          : [pcr_amplicons],
                            'assembly_names' : [assembly_name]
                            })]

    assert len(assembly[0]) == 8993
    assert str(assembly[0].seq[:20]) =='GTTTGTAGTTGGCGGTGGAG'


def test_find_sequence_location(): 
    
    test_plasmid = SeqIO.read('../teemi/tests/files_for_testing/MIA-HA-1.gb', format = 'genbank')
    sgRNA = SeqRecord(Seq('TGACGAATCGTTAGGCACAG'), name = 'random_sgRNA', id = '1483', description = 'This is a test sgRNA')
    start_end_location = find_sequence_location(sgRNA, test_plasmid)
    
    assert start_end_location[0] == 2
    assert start_end_location[1] == 22
    assert start_end_location[2] == 1

    # reverse complement
    sgRNA = sgRNA.reverse_complement()
    start_end_location = find_sequence_location(sgRNA, test_plasmid)

    assert start_end_location[0] == 22
    assert start_end_location[1] == 2
    assert start_end_location[2] == -1


def test_crispr_db_break_location(): 
    
    assert crispr_db_break_location(220,200, -1) == 217
    assert crispr_db_break_location(200,220, 1) == 217


def test_add_feature_annotation_to_seqrecord(): 
    test_plasmid = SeqIO.read('../teemi/tests/files_for_testing/MIA-HA-1.gb', 'gb')
    add_feature_annotation_to_seqrecord(test_plasmid,label=f'This a test')
    
    assert test_plasmid.features[0].qualifiers['label'] == 'This a test'
    # Directly use the ExactPosition as an integer
    assert int(test_plasmid.features[-1].location.start) == 0
    assert int(test_plasmid.features[-1].location.end) == len(test_plasmid)



def test_find_all_occurences_of_a_sequence(): 
    test_plasmid = SeqIO.read('../teemi/tests/files_for_testing/MIA-HA-1.gb', 'gb')
    sgRNA = SeqRecord(Seq('attcattaccatagtattact'), name = 'random_sgRNA', id = '1483', description = 'This is a test sgRNA')
    do_we_have_multiple_hits_on_the_genome = find_all_occurrences_of_a_sequence(sgRNA, test_plasmid)

    assert do_we_have_multiple_hits_on_the_genome == 1
    
    sgRNA = SeqRecord(Seq('TGACGAATCGTTAGGCACAG'), name = 'random_sgRNA', id = '1483', description = 'This is a test sgRNA')
    do_we_have_multiple_hits_on_the_genome = find_all_occurrences_of_a_sequence(sgRNA, test_plasmid)

    assert do_we_have_multiple_hits_on_the_genome == 1

    sgRNA = SeqRecord(Seq('CTATTTTTTTCTGCTTACGCGAGAGAGAGATAGATAGA'), name = 'random_sgRNA', id = '1483', description = 'This is a test sgRNA')
    do_we_have_multiple_hits_on_the_genome = find_all_occurrences_of_a_sequence(sgRNA, test_plasmid)

    assert do_we_have_multiple_hits_on_the_genome == 0

    sgRNA = SeqRecord(Seq('CTAT'), name = 'random_sgRNA', id = '1483', description = 'This is a test sgRNA')
    do_we_have_multiple_hits_on_the_genome = find_all_occurrences_of_a_sequence(sgRNA, test_plasmid)

    assert do_we_have_multiple_hits_on_the_genome == 27
