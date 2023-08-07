# import pytest
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from teemi.build.PCR import primer_tm_neb
# from Bio.SeqFeature import SeqFeature, FeatureLocation
# from Bio import SeqIO 

from teemi.design.golden_gate_cloning import (retrieve_primers_from_golden_gate,
                         make_GG_overhangs_forward_primers, make_GG_overhangs_reverse_primers,
                         make_amplicons, GoldenGateCloning, digest_amplicons_w_BsaI)
import pytest
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydna.primer import Primer
import pytest
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer
from pydna.amplicon import Amplicon
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt



def test_make_GG_overhangs_reverse_primers():
    # Simple example
    amplicon1 = Dseqrecord(Dseq("AAAAGTCTAGAGGATCC"))
    amplicon1.reverse_primer = Primer(Dseq("GGATCC"))
    amplicon2 = Dseqrecord(Dseq("AAAAGTCTAGAGGATCC"))
    amplicon2.reverse_primer = Primer(Dseq("GGATCC"))
    amplicons = [amplicon1, amplicon2]

    sgRNA1 = Dseqrecord(Dseq("GGTCTAGA"))
    sgRNA2 = Dseqrecord(Dseq("GGTCTAGA"))
    sgRNAs = [sgRNA1, sgRNA2]

    modified_amplicons = make_GG_overhangs_reverse_primers(amplicons, sgRNAs)

    assert len(modified_amplicons) == 2
    assert str(modified_amplicons[0].reverse_primer.seq) == "gATCAGGTCTCGGACCGGATCC"
    assert str(modified_amplicons[1].reverse_primer.seq) == "gATCAGGTCTCGCTAgGGATCC"


    # more complex
    # Generate amplicons
    sgRNA_handle_cys4_sites = [Dseqrecord('GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT', name = 'sgRNA_handle_cys4')]*5
    sgRNA_handle_cys4_sites_as_amplicons =make_amplicons(sgRNA_handle_cys4_sites,target_tm=55 , tm_function=mt.Tm_Wallace )
    # Add sgRNAs
    sgRNA = [Dseqrecord('GAGCAGTTCCCAGAACTGCC'), Dseqrecord('GTCGACCTCCCAGTCACGGC'), Dseqrecord('GCCGACCGCCCAGGCGACCT'),Dseqrecord('CCGTTCACAGGTCGCGGCGG'), Dseqrecord('ACCGCCCAGGCGACCTCGGC')]


    f_primers = make_GG_overhangs_forward_primers(sgRNA_handle_cys4_sites_as_amplicons, sgRNA, restriction_overhang_f = "GATCGggtctcc",     backbone_overhang_f = "cATG")
    f_r_primers = make_GG_overhangs_reverse_primers(f_primers, sgRNA,     backbone_overhang_r = "cTAG", restriction_overhang_r = "GATCAGGTCTCg")

    assert str(f_r_primers[1].reverse_primer.seq) == 'GATCAGGTCTCgCGGCAAAAAAGCACCGACTCG'
    assert len(str(f_r_primers[1].reverse_primer.seq)) == 33
    assert str(f_r_primers[4].reverse_primer.seq) == 'GATCAGGTCTCgCTAgAAAAAAGCACCGACTCG'
    assert len(str(f_r_primers[4].reverse_primer.seq)) == 33

def test_make_GG_overhangs_forward_primers():

    # Generate amplicons
    sgRNA_handle_cys4_sites = [Dseqrecord('GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT', name = 'sgRNA_handle_cys4')]*5
    sgRNA_handle_cys4_sites_as_amplicons =make_amplicons(sgRNA_handle_cys4_sites,target_tm=55 , tm_function=mt.Tm_Wallace )
    # Add sgRNAs
    sgRNA = [Dseqrecord('GAGCAGTTCCCAGAACTGCC'), Dseqrecord('GTCGACCTCCCAGTCACGGC'), Dseqrecord('GCCGACCGCCCAGGCGACCT'),Dseqrecord('CCGTTCACAGGTCGCGGCGG'), Dseqrecord('ACCGCCCAGGCGACCTCGGC')]

    # run funciton
    f_primers = make_GG_overhangs_forward_primers(sgRNA_handle_cys4_sites_as_amplicons, sgRNA, restriction_overhang_f = "GATCGggtctcc",     backbone_overhang_f = "cATG")

    assert len(f_primers[0].forward_primer) == 85
    assert str(f_primers[0].forward_primer.seq) == 'GATCGggtctcccATGgTTCACTGCCGTATAGGCAGCTAAGAAAGAGCAGTTCCCAGAACTGCCGTTTTAGAGCTAGAAATAGCA'
    assert len(f_primers[1].forward_primer) == 53
    assert str(f_primers[1].forward_primer.seq) == 'GATCGggtctccGTCGACCTCCCAGTCACGGCGTTTTAGAGCTAGAAATAGCA'
    
    assert len(f_primers) == 5


def test_retrieve_primers_from_golden_gate():
    primers = ['ATCG', 'GCTA', 'CGAT']
    primers_footprint = ['A', 'GC', 'TAC']
    result = retrieve_primers_from_golden_gate(primers, primers_footprint)
    assert len(result) == 3
    assert result[0].seq == 'ATCG'
    assert len(result[0].footprint) == len(str(primers_footprint[0]))
    assert result[1].seq == 'GCTA'
    assert len(result[1].footprint) == len(primers_footprint[1])
    assert result[2].seq == 'CGAT'
    assert len(result[2].footprint) == len(primers_footprint[2])

def test_make_amplicons(): 
    # genrate some random amplicons
    sgRNA_handle_cys4_sites = [Dseqrecord('GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT', name = 'sgRNA_handle_cys4')]*5
    sgRNA_handle_cys4_sites_as_amplicons =make_amplicons(sgRNA_handle_cys4_sites,target_tm=55 , tm_function=mt.Tm_Wallace )


    assert sgRNA_handle_cys4_sites_as_amplicons[0].seq.watson == 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT'
    assert len(sgRNA_handle_cys4_sites_as_amplicons[0].seq) == 82
    assert len(sgRNA_handle_cys4_sites_as_amplicons) == 5


def test_GoldenGate(): 
    
    sgRNA_handle_cys4_sites = [Dseqrecord('GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT', name = 'sgRNA_handle_cys4')]*5
    sgRNA_list = [Dseqrecord('GAGCAGTTCCCAGAACTGCC'), Dseqrecord('GTCGACCTCCCAGTCACGGC'), Dseqrecord('GCCGACCGCCCAGGCGACCT'),Dseqrecord('CCGTTCACAGGTCGCGGCGG'), Dseqrecord('ACCGCCCAGGCGACCTCGGC')]

    golden_gate = GoldenGateCloning(sgRNA_list,
                                sgRNA_handle_cys4_sites, 
                                target_tm= 60, tm_function=mt.Tm_Wallace  ) #, backbone_overhang_r = 'ctag' )

# def test_make_amplicons():
#     amplicons = [Dseqrecord(Seq("ATCG"))]
#     result = make_amplicons(amplicons)
#     assert len(result) == 1
#     assert isinstance(result[0], Amplicon)

# # def test_GoldenGateCloning():
# #     sgRNAs = [Dseqrecord(Seq("GCTA"))]
# #     amplicons = [Dseqrecord(Seq("ATCG"))]
# #     gg_cloning = GoldenGateCloning(sgRNAs, amplicons)
# #     assert len(gg_cloning.amplicons) == 1
# #     assert isinstance(gg_cloning.amplicons[0], Amplicon)
# #     assert len(gg_cloning.amplicons_w_f_primer_overhang) == 1
# #     assert isinstance(gg_cloning.amplicons_w_f_primer_overhang[0], Amplicon)
# #     assert len(gg_cloning.amplicons_w_f_r_primer_overhang) == 1
# #     assert isinstance(gg_cloning.amplicons_w_f_r_primer_overhang[0], Amplicon)

# def test_digest_amplicons_w_BsaI():
#     amplicons = [Dseqrecord(Seq("ATCG"))]
#     result = digest_amplicons_w_BsaI(amplicons)
#     assert len(result) == 1
#     assert isinstance(result[0], Dseqrecord)


