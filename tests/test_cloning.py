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
from teemi.design.gibson_cloning import extract_locus_tag_homology_arms


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


def test_extract_locus_tag_homology_arms():
    feature = SeqFeature(FeatureLocation(100, 160, strand=1), type="CDS")
    feature.qualifiers["locus_tag"] = ["AO_TEST_001"]
    genome = SeqRecord(Seq("A" * 300), id="chrom1", name="chrom1", features=[feature])

    homology_arms = extract_locus_tag_homology_arms(
        genome, ["AO_TEST_001"], arm_length=45
    )

    assert homology_arms.shape[0] == 1
    assert homology_arms.loc[0, "locus_tag"] == "AO_TEST_001"
    assert homology_arms.loc[0, "upstream_arm_name"] == "AO_TEST_001_UP45"
    assert homology_arms.loc[0, "downstream_arm_name"] == "AO_TEST_001_DW45"
    assert homology_arms.loc[0, "repair_oligo_name"] == "AO_TEST_001_DEL_45bp_arms"
    assert homology_arms.loc[0, "upstream_arm_start"] == 55
    assert homology_arms.loc[0, "upstream_arm_end"] == 100
    assert homology_arms.loc[0, "downstream_arm_start"] == 160
    assert homology_arms.loc[0, "downstream_arm_end"] == 205
    assert homology_arms.loc[0, "upstream_arm_length"] == 45
    assert homology_arms.loc[0, "downstream_arm_length"] == 45
    assert homology_arms.loc[0, "repair_oligo_length"] == 90
    assert homology_arms.loc[0, "upstream_arm"] == "A" * 45
    assert homology_arms.loc[0, "downstream_arm"] == "A" * 45
    assert homology_arms.loc[0, "repair_oligo"] == "A" * 90



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


# ---------------------------------------------------------------------------
# Additional tests
# ---------------------------------------------------------------------------
import random
import pandas as pd
import pytest
from pydna.dseq import Dseq


def _random_dna(length, seed):
    """Deterministic pseudo random DNA."""
    rng = random.Random(seed)
    return ''.join(rng.choice('ATGC') for _ in range(length))


def test_CAS9_cutting_warns_when_gRNA_is_not_found(capsys):
    background = Dseqrecord('ACGT' * 20, name='background')
    gRNA = Dseqrecord('TTTTTTTTTTTTTTTTTTTT', name='missing_gRNA')

    CAS9_cutting(gRNA, background)

    assert "CAN'T FIND THE CUT SITE IN YOUR SEQUENCE" in capsys.readouterr().out


def test_CAS9_cutting_warns_when_gRNA_cuts_twice(capsys):
    gRNA_seq = 'TTTTGGGGCCCCAAAATTTT'
    background = Dseqrecord('AAAA' + gRNA_seq + 'CGCGCGTATA' + gRNA_seq + 'GGGG',
                            name='background')
    gRNA = Dseqrecord(gRNA_seq, name='double_gRNA')

    up, dw = CAS9_cutting(gRNA, background)

    # cut 17 bp into the first occurrence of the gRNA
    assert len(up) == 4 + 17
    assert len(dw) == len(background) - len(up)
    assert 'cuts more than one location' in capsys.readouterr().out


def test_extract_template_amplification_sites_without_strand():
    cds = SeqFeature(FeatureLocation(1, 100), type='CDS')
    cds.qualifiers = {'name': 'ATF'}
    terminator = SeqFeature(FeatureLocation(101, 300), type='terminator')
    terminator.qualifiers = {'name': 'ATF_tTEF'}
    template = Dseqrecord(Seq('ATGATGA' * 100), name='ATF')
    template.features = [cds, terminator]

    # features without an explicit strand are treated as being on the plus strand
    assert cds.location.strand is None
    sites = extract_template_amplification_sites([template], ['ATF'], 'tTEF')

    assert len(sites) == 1
    assert str(sites[0].seq) == str(template.seq)[1:300]


def test_extract_template_amplification_sites_on_minus_strand():
    cds = SeqFeature(FeatureLocation(200, 300, strand=-1), type='CDS')
    cds.qualifiers = {'name': 'BTF'}
    terminator = SeqFeature(FeatureLocation(1, 100, strand=-1), type='terminator')
    terminator.qualifiers = {'name': 'BTF_tTEF'}
    template = Dseqrecord(Seq('ATGCATG' * 100), name='BTF')
    template.features = [cds, terminator]

    sites = extract_template_amplification_sites([template], ['BTF'], 'tTEF')

    # on the minus strand the CDS gives the end and the terminator the start
    assert len(sites) == 1
    assert len(sites[0].seq) == 299
    assert str(sites[0].seq) == str(template.seq)[1:300]


def test_extract_template_amplification_sites_appends_to_existing_batches():
    cds = SeqFeature(FeatureLocation(1, 100, strand=1), type='CDS')
    cds.qualifiers = {'name': 'ATF'}
    terminator = SeqFeature(FeatureLocation(101, 300, strand=1), type='terminator')
    terminator.qualifiers = {'name': 'ATF_tTEF'}
    template = Dseqrecord(Seq('ATGATGA' * 100), name='ATF')
    template.features = [cds, terminator]
    template.annotations['batches'] = ['batch_one']

    sites = extract_template_amplification_sites([template], ['ATF'], 'tTEF')

    # the slice inherits 'batches' and the template batch is appended to it
    assert sites[0].annotations['batches'] == ['batch_one', 'batch_one']


def test_nicking_enzyme():
    middle = _random_dna(30, 21)
    watson = 'CGCGTG' + middle + 'CGTGCG'
    crick = str(Seq(watson).reverse_complement()) + 'AT'
    vector = Dseqrecord(Dseq(watson=watson, crick=crick, ovhg=2))

    nicked = nicking_enzyme(vector)

    assert isinstance(nicked, Dseq)
    # the 6 bp recognition site is removed from both strands
    assert nicked.watson == watson[6:]
    assert nicked.crick == crick[6:]
    assert nicked.ovhg == 8


def test_nicking_enzyme_without_recognition_site(capsys):
    vector = Dseqrecord(_random_dna(40, 22))

    assert nicking_enzyme(vector) is None
    assert 'No nicking sequnce' in capsys.readouterr().out


def test_Nt_Bbc_CI():
    watson = 'GCGATCGCT' + _random_dna(26, 23) + 'TTTAAT'
    crick = 'TAAACC' + str(Seq(watson).reverse_complement()) + 'GC'
    linearized_vector = Dseqrecord(Dseq(watson=watson, crick=crick, ovhg=2))

    nicked = Nt_Bbc_CI(linearized_vector)

    assert isinstance(nicked, Dseq)
    assert nicked.watson == watson[6:]
    assert nicked.crick == crick[6:]
    assert nicked.ovhg == 8


def test_Nt_Bbc_CI_raises_without_nicking_site():
    linearized_vector = Dseqrecord(_random_dna(40, 24))

    with pytest.raises(ValueError, match='No nicking sequence found'):
        Nt_Bbc_CI(linearized_vector)


def test_casembler_verbose_writes_genbank_files(tmp_path, monkeypatch):
    site_seq = _random_dna(300, 11)
    genome_seq = _random_dna(50, 12) + site_seq + _random_dna(50, 13)
    site_feature = SeqFeature(FeatureLocation(50, 350, strand=1), type='misc_feature')
    site_feature.qualifiers['name'] = ['X-1']
    bg_strain = SeqRecord(Seq(genome_seq), id='bg', name='bg', features=[site_feature])

    # cut 17 bp into the gRNA, i.e. at position 117 of the site
    gRNA = Dseqrecord(site_seq[100:120], name='gRNA1')
    up_seq, dw_seq = site_seq[:117], site_seq[117:]
    insert = _random_dna(60, 14)
    repair_template = Dseqrecord(up_seq[-30:] + insert + dw_seq[:30], name='part1')

    monkeypatch.chdir(tmp_path)
    assembly = casembler(
        bg_strain=bg_strain,
        site_names=['X-1'],
        gRNAs=[gRNA],
        parts=[[repair_template]],
        assembly_limits=[30],
        assembly_names=['testasm'],
        verbose=True,
    )

    assert assembly.name == 'testasm'
    assert str(assembly.seq) == up_seq + insert + dw_seq
    assert len(assembly) == len(site_seq) + len(insert)
    assert assembly.features[-1].qualifiers['name'] == 'X-1'
    assert assembly.features[-1].qualifiers['label'] == 'X-1'
    assert int(assembly.features[-1].location.end) == len(assembly)

    # verbose=True writes the up, assembly and down sequences as genbank files
    written = sorted(p.name for p in tmp_path.glob('*.gb'))
    assert written == ['DW_gRNA1_X-1.gb', 'UP_gRNA1_X-1.gb', 'testasm.gb']
    reread = SeqIO.read(tmp_path / 'testasm.gb', 'gb')
    assert str(reread.seq) == str(assembly.seq)


def test_plate_plot():
    amplicon_df = pd.DataFrame(
        {
            'name': ['PCR_01', 'PCR_02', 'PCR_03', 'PCR_04'],
            'prow': ['A', 'A', 'B', 'B'],
            'pcol': [1, 2, 1, 2],
            'ignored': [0, 0, 0, 0],
        }
    )

    plate = plate_plot(amplicon_df, 'name')

    assert list(plate.index) == ['A', 'B']
    assert plate.index.name == 'prow'
    assert list(plate.columns.get_level_values(0).unique()) == ['name']
    assert list(plate.columns.get_level_values(1)) == [1, 2]
    assert plate.loc['A', ('name', 1)] == 'PCR_01'
    assert plate.loc['B', ('name', 2)] == 'PCR_04'


def test_seq_to_annotation_on_reverse_strand():
    chromosome = SeqRecord(Seq('ATGCTAGCCAGTCGTAAACCGGTT'), id='chrom1', name='chrom1')
    # the reverse complement of chromosome[4:11]
    reverse_hit = SeqRecord(Seq('CTGGCTA'), id='reverse_hit', name='reverse_hit')

    seq_to_annotation(reverse_hit, chromosome, 'misc_feature')

    feature = chromosome.features[-1]
    assert feature.type == 'misc_feature'
    # start/end are swapped so that the feature location stays ascending
    assert int(feature.location.start) == 4
    assert int(feature.location.end) == 11
    assert feature.qualifiers['label'] == 'reverse_hit'
