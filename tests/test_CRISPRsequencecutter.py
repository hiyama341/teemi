#!/usr/bin/env python
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dataclasses import asdict

from teemi.design.CRISPRsequencecutter import CRISPRSequenceCutter

# def test_crispr_cutting_instance():
#     chromosome = SeqRecord(Seq("CTAACTATTCCCTATGTCTTATAGGGGCCTACGTTATCTGCCTGTCGAAC"))
#     sequence = SeqRecord(Seq("TATAGGGGCCTACGTTATCT"))

#     crispr_cutting_instance = CRISPRSequenceCutter(sequence=sequence, chromosome=chromosome)
#     dict_with_cut = asdict(crispr_cutting_instance)


#     assert str(dict_with_cut['sequence'].seq) =='TATAGGGGCCTACGTTATCT'    
#     assert str(dict_with_cut['chromosome'].seq) =='CTAACTATTCCCTATGTCTTATAGGGGCCTACGTTATCTGCCTGTCGAAC'


# def test_get_upstream_and_downstream_sequences(crispr_cutting_instance):
#     upstream, downstream = crispr_cutting_instance.get_upstream_and_downstream_sequences()
#     assert upstream.seq == Seq("ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC")
#     assert downstream.seq == Seq("TAGCTAGCTAGCTAGCTAGCT")


# def test_add_annotation_to_upstream_and_downstream_sequences(crispr_cutting_instance):
#     crispr_cutting_instance.add_annotation_to_upstream_and_downstream_sequences()
#     assert crispr_cutting_instance.upstream_crispr_sequence.features[0].qualifiers["label"] == [
#         "upstream_sequence_of_CRISPR_cut"
#     ]
#     assert crispr_cutting_instance.downstream_crispr_sequence.features[0].qualifiers["label"] == [
#         "downstream_sequence_of_CRISPR_cut"
#     ]


# def test_display_cut(crispr_cutting_instance):
#     assert crispr_cutting_instance.display_cut() == "CTAGCTAGCTAGCTAGCTAG---x----TAGCTAGCTAGCTAGCTAGC"

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

@pytest.fixture
def mock_chromosome():
    # Create a mock chromosome to use in testing
    seq = Seq("GCGCTTGGTAAGGCCATCGCGAATACCAGGTATCGTGTAAGTAGCGTAGGCCCGCACGCAAGATAAACTGCTAGGGAACCGCGTTTCCACGACCGGTGCACGATTTAATTTCGCCGACGTGATGACATTCCAGGCAGTGCCTCTGCCGCCGGACCCCTCTCGTGATTGGGTAGCTGGACATGCCCTTGTAAGATATAACA")
    return SeqRecord(seq, name="chromosome1", id = "chromosome1")

def test_valid_sequence_length(mock_chromosome):
    # Test that CRISPRSequenceCutter raises a ValueError if the sgRNA is not 20 bp long
    with pytest.raises(ValueError):
        CRISPRSequenceCutter(Seq("GCACGATTTAATTTCGCCGACGTGATGAC"), mock_chromosome)

def test_valid_sgRNA_site(mock_chromosome):
    # Test that CRISPRSequenceCutter raises a ValueError if there are not exactly 1 sgRNA site in the chromosome
    with pytest.raises(ValueError):
        CRISPRSequenceCutter("AGCT" * 10, mock_chromosome)

def test_valid_crispr_db_break_location(mock_chromosome):
    # Test that the CRISPR break location is calculated correctly
    cutter = CRISPRSequenceCutter(SeqRecord(Seq("AATTTCGCCGACGTGATGAC")), mock_chromosome)
    assert cutter.get_crispr_db_break_location() == 123

def test_valid_upstream_and_downstream_sequences(mock_chromosome):
    # Test that the upstream and downstream sequences are extracted correctly
    cutter = CRISPRSequenceCutter(SeqRecord(Seq('AATTTCGCCGACGTGATGAC')), mock_chromosome)
    upstream, downstream = cutter.upstream_crispr_sequence, cutter.downstream_crispr_sequence
    assert str(upstream.seq) == 'GCGCTTGGTAAGGCCATCGCGAATACCAGGTATCGTGTAAGTAGCGTAGGCCCGCACGCAAGATAAACTGCTAGGGAACCGCGTTTCCACGACCGGTGCACGATTTAATTTCGCCGACGTGAT'
    assert str(downstream.seq) == 'GACATTCCAGGCAGTGCCTCTGCCGCCGGACCCCTCTCGTGATTGGGTAGCTGGACATGCCCTTGTAAGATATAACA'


def test_valid_display_cut(mock_chromosome):
    # Test that the display_cut function returns the correct string
    cutter = CRISPRSequenceCutter(SeqRecord(Seq('AATTTCGCCGACGTGATGAC')), mock_chromosome)
    assert cutter.display_cut() == 'TTTAATTTCGCCGACGTGAT---x----GACATTCCAGGCAGTGCCTC'
