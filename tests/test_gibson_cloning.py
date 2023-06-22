import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from teemi.build.PCR import primer_tm_neb
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO 

from teemi.design.gibson_cloning import find_up_dw_repair_templates, update_primer_names, assemble_single_plasmid_with_repair_templates, assemble_multiple_plasmids_with_repair_templates_for_deletion
from Bio.SeqUtils import MeltingTemp as mt
from pydna.dseqrecord import Dseqrecord
from Bio.Restriction import StuI


# read in a geneome to use as repair temolate for all tests
genome_path = '../teemi/tests/files_for_testing/gibson/Streptomyces_coelicolor_A3_chromosome.gb'
StrepA3 = SeqIO.read(genome_path, "gb")
repair_templates = ['SCO5892']

# Call the function to test
repair_DNA_templates = find_up_dw_repair_templates(StrepA3, repair_templates, target_tm = 55, primer_calc_function= mt.Tm_NN)


def test_find_up_dw_repair_templates():

    # Assert the expected output
    assert len(repair_DNA_templates) == 1

    # Assert the content of the first record
    record1 = repair_DNA_templates[0]
    assert record1["name"] == "SCO5892"
    assert len(record1["up_repair"]) == 1000

    assert int(record1["tm_up_forwar_p"]) == 54
    assert int(record1["tm_up_reverse_p"]) == 55
    assert record1["location_up_start"] == 6448548
    assert record1["location_up_end"] == 6449548

    assert len(record1["dw_repair"]) == 1000
    assert int(record1["tm_dw_forwar_p"]) == 55
    assert int(record1["tm_dw_reverse_p"]) == 53
    assert record1["location_dw_start"] == 6456442
    assert record1["location_dw_end"] == 6457442


def test_update_primer_names():
    # Create mock data for testing
    list_of_records = [
        {
            "gene_name": "gene1",
            "up_forwar_p": SeqRecord(Seq("ATCG"), id="up_F0"),
            "up_reverse_p": SeqRecord(Seq("GCAT"), id="up_R0"),
            "dw_forwar_p": SeqRecord(Seq("GCTA"), id="dw_F0"),
            "dw_reverse_p": SeqRecord(Seq("TACG"), id="dw_R0"),
        },
        {
            "gene_name": "gene2",
            "up_forwar_p": SeqRecord(Seq("AATT"), id="up_F1"),
            "up_reverse_p": SeqRecord(Seq("TTAA"), id="up_R1"),
            "dw_forwar_p": SeqRecord(Seq("CGCG"), id="dw_F1"),
            "dw_reverse_p": SeqRecord(Seq("GCGC"), id="dw_R1"),
        },
        # Add more mock data as needed
    ]

    # Call the function to test
    update_primer_names(list_of_records)

    # Assert the content of the first updated record
    record1 = list_of_records[0]
    assert record1["gene_name"] == "gene1"
    assert record1["up_forwar_p_name"] == "gene1_up_F0"
    assert record1["up_reverse_p_name"] == "gene1_up_R0"
    assert record1["dw_forwar_p_name"] == "gene1_dw_F0"
    assert record1["dw_reverse_p_name"] == "gene1_dw_R0"
    # Add more assertions for the remaining fields

    # Assert the content of the second updated record
    record2 = list_of_records[1]
    assert record2["gene_name"] == "gene2"
    assert record2["up_forwar_p_name"] == "gene2_up_F1"
    assert record2["up_reverse_p_name"] == "gene2_up_R1"
    assert record2["dw_forwar_p_name"] == "gene2_dw_F1"
    assert record2["dw_reverse_p_name"] == "gene2_dw_R1"
 


def test_assemble_single_plasmid_with_repair_templates():
    #import a plasmid
    CRISPR_plasmids_wo_repair = SeqIO.read('../teemi/tests/files_for_testing/gibson/SCO5087_oligo1.gb', "gb") 

    vector = Dseqrecord(CRISPR_plasmids_wo_repair, circular = True)
    vector.name = 'pCRISPR–1_SCO5892'

    # restriction digest
    vector_StulI = vector.cut(StuI)[0]
    vector_StulI.name = 'pCRISPR–1_SCO5892'

    #lets_assemble our constructs
    repair_names = ['up_SCO5892', 'dw_SCO5892']
    repair_templates_sc000x = [repair_DNA_templates[0]['up_repair'], repair_DNA_templates[0]['dw_repair']]
    repair_templates_sc000x[0].name, repair_templates_sc000x[1].name = repair_names[0], repair_names[1]

    # ASSEMBLE it
    first_vector = assemble_single_plasmid_with_repair_templates(repair_templates_sc000x,vector_StulI , overlap = 35)

    assert len(first_vector) == 4

    assert len(first_vector[0]) == 11279
    assert first_vector[0].name == 'pCRISPR–1_SCO5892'
    assert first_vector[0].seq[:10].watson =='CCTATATGCG'

    assert len(first_vector[1]) == 1053
    assert first_vector[1].name == '1053bp_PCR_prod'
    assert first_vector[1].seq[:10].watson == 'CGCTTTTGCG'

    assert len(first_vector[2]) == 1053
    assert first_vector[2].name == '1053bp_PCR_prod'
    assert first_vector[2].seq[:10].watson == 'CTCGGAACGG'

    assert len(first_vector[3]) == 11279
    assert first_vector[0].name == 'pCRISPR–1_SCO5892'



def test_assemble_multiple_plasmids_with_repair_templates_for_deletion():
    # import plasmids
    CRISPR_plasmids_wo_repair1 = [SeqIO.read('../teemi/tests/files_for_testing/gibson/SCO5892_oligo4.gb', "gb") ]
    CRISPR_plasmids_wo_repair2 = [SeqIO.read('../teemi/tests/files_for_testing/gibson/SCO5892_oligo5.gb', "gb") ]
    CRISPR_plasmids_wo_repair = CRISPR_plasmids_wo_repair1 + CRISPR_plasmids_wo_repair2

    # digest the plasmids
    digested_plasmids = [Dseqrecord(digest, circular = True).cut(StuI)[0] for digest in CRISPR_plasmids_wo_repair]
    # rename them appropriatly
    for i in range(len( digested_plasmids)):
        digested_plasmids[i].name = CRISPR_plasmids_wo_repair[i].name
    
    # run function
    list_of_gene_names = ["SCO5892"]
    list_of_digested_plasmids = digested_plasmids

    list_of_records = assemble_multiple_plasmids_with_repair_templates_for_deletion(list_of_gene_names,list_of_digested_plasmids, repair_DNA_templates, overlap = 40 )

    # first plasmid
    record1 = list_of_records[0]
    assert record1["name"] == "SCO5892_oligo4_pCRISPR"
    assert len(record1["contig"]) == 13279

    assert int(record1["tm_up_forwar_p"]) == 69
    assert int(record1["tm_up_reverse_p"]) == 71
    assert record1["up_forwar_p_anneal"] ==  Seq('CGACGAGCTGGACGTCG')
    assert record1["up_reverse_p_anneal"] == Seq('CTACCGGGCCGTTCCGA')

    assert int(record1["tm_up_forwar_p"]) == 69
    assert int(record1["tm_up_reverse_p"]) == 71
    assert record1["up_forwar_p_anneal"] ==  Seq('CGACGAGCTGGACGTCG')
    assert record1["up_reverse_p_anneal"] == Seq('CTACCGGGCCGTTCCGA')

    # second plasmid
    record2 = list_of_records[1]
    assert record2["name"] == "SCO5892_oligo5_pCRISPR"
    assert len(record2["contig"]) == 13279
    assert int(record2["tm_up_forwar_p"]) == 69
    assert int(record2["tm_up_reverse_p"]) == 71
    assert record2["up_forwar_p_anneal"] ==  Seq('CGACGAGCTGGACGTCG')
    assert record2["up_reverse_p_anneal"] == Seq('CTACCGGGCCGTTCCGA')
    assert int(record2["tm_up_forwar_p"]) == 69
    assert int(record2["tm_up_reverse_p"]) == 71
    assert record2["up_forwar_p_anneal"] ==  Seq('CGACGAGCTGGACGTCG')
    assert record2["up_reverse_p_anneal"] == Seq('CTACCGGGCCGTTCCGA')