#!/usr/bin/env python
# MIT License
# Copyright (c) 2023, Technical University of Denmark (DTU)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer
from teemi.build.PCR import primer_ta_neb, primer_tm_neb
from pydna.design import primer_design
from Bio.Seq import Seq
import pandas as pd


def retrieve_primers_from_golden_gate(primers, primers_footprint):
    """Makes primers with footprints from two list of strings"""
    primers_list = []
    for i in range(len(primers)):
        primers_list.append(Primer(primers[i], footprint=len(primers_footprint[i])))
    return primers_list


def make_GG_overhangs_forward_primers(
    list_of_amplicons: list,
    sgRNA_list: list,
    restriction_overhang_f="GATCGggtctcc",
    backbone_overhang_f="cATG",
    cys4="gTTCACTGCCGTATAGGCAGCTAAGAAA",
):
    """Generates forward primer overhangs to accomodate GoldenGate cloning.

    Parameters
    ----------
    list_of_amplicons: list
        list of pydna.Dseqrecords
    sgRNA_list: list
        list of pydna.Dseqrecords
    - restriction_overhang_f: str
        representing the restriction overhang sequence for the forward primer (default='GATCGggtctcc')
    - backbone_overhang_f: str
        representing the backbone overhang sequence for the forward primer (default='cATG')
    - cys4: str
        representing the Cys4 sequence (default='gTTCACTGCCGTATAGGCAGCTAAGAAA')

    Returns
    -------
    list_of_amplicons: list
        lsit of pydna.amplicon objects with modified forward primer sequences

    """
    # First overhang
    list_of_amplicons[0].forward_primer.seq = (
        restriction_overhang_f
        + backbone_overhang_f
        + cys4
        + sgRNA_list[0].seq.watson
        + list_of_amplicons[0].forward_primer.seq
    )
    list_of_amplicons[0].forward_primer.name, list_of_amplicons[0].forward_primer.id = (
        f"F{0}",
        f"F{0}",
    )

    # the rest - we start at index 1.
    for i in range(1, len(list_of_amplicons)):
        list_of_amplicons[i].forward_primer.seq = (
            restriction_overhang_f
            + sgRNA_list[i].seq.watson
            + list_of_amplicons[i].forward_primer.seq
        )
        # change_name
        (
            list_of_amplicons[i].forward_primer.name,
            list_of_amplicons[i].forward_primer.id,
        ) = (f"F{i}", f"F{i}")

    return list_of_amplicons


def make_GG_overhangs_reverse_primers(
    list_of_amplicons: list,
    sgRNA_list: list,
    backbone_overhang_r="cTAG",
    restriction_overhang_r="gATCAGGTCTCG",
):
    """Generates reverse primer overhangs to accomodate GoldenGate cloning.

    Parameters
    ----------
    list_of_amplicons: list
        list of pydna.Dseqrecords
    sgRNA_list: list
        list of pydna.Dseqrecords

    backbone_overhang_r: str
        representing the backbone overhang sequence for the reverse primer (default='cTAG')
    restriction_overhang_r: str
        representing the restriction overhang sequence for the reverse primer (default='gATCAGGTCTCG')

    Returns
    -------
     -------
    list_of_amplicons: list
        list of pydna.amplicon objects with modified reverse primer sequences

    """
    golden_gate_overhangs = []
    # Find overhangs  - last is fixed
    for i in range(1, len(sgRNA_list)):
        golden_gate_overhangs.append(sgRNA_list[i].seq.watson[0:4])
    # last amplicon has fixed overhnag
    golden_gate_overhangs.append(backbone_overhang_r)

    # We add restriction site, overhangs, primer seq
    for i in range(0, len(list_of_amplicons)):
        list_of_amplicons[i].reverse_primer.seq = (
            restriction_overhang_r
            + Seq(golden_gate_overhangs[i]).reverse_complement()
            + list_of_amplicons[i].reverse_primer.seq
        )
        # change_name/id
        (
            list_of_amplicons[i].reverse_primer.name,
            list_of_amplicons[i].reverse_primer.id,
        ) = (f"R{i}", f"R{i}")

        list_of_amplicons[i].name = sgRNA_list[i].name

    return list_of_amplicons


def make_amplicons(
    list_of_amplicons: list, target_tm=55, limit=10, tm_function=primer_tm_neb
):
    """Generates pydna.amplicons which contains primers with a target temperature.

    Parameters
    ----------
    list_of_amplicons : list
        list of pydna.Dseqrecords
    target_tm : int
        representing the target melting temperature for the primers (default=55)
    limit: int
        representing the maximum primer size (default=5)
    tm_function : function
        for calculating primer melting temperature (default=primer_tm_neb)

    Returns:
    amplicons: list
        list of amplicon objects with designed primer sequences
    """
    amplicons = []
    for i in range(len(list_of_amplicons)):
        amplicon = primer_design(
            list_of_amplicons[i],
            target_tm=target_tm,
            limit=limit,
            tm_function=tm_function,
        )

        amplicons.append(amplicon)

    return amplicons


from dataclasses import dataclass, field
from pydna.amplify import pcr
from Bio.Restriction import BsaI
from typing import Callable


@dataclass
class GoldenGateCloning:
    """Data class for making golden gate cloning a brease.

    Parameters:
        sgRNAs: list
            list of sgRNA sequences
        list_of_amplicons: list
            list of amplicon sequences
        target_tm: int, default = 55
            target melting temperature for the primers
        tm_function : function, default = primer_tm_neb
            function to calculate melting temperature
        restriction_overhang_f: str, default = 'GATCGggtctcc'
            overhang to be added to the forward primer
        restriction_overhang_r: str, default = 'GATCAGGTCTCg'
            overhang to be added to the reverse primer
        backbone_overhang_f:str, default = 'cATG'
            overhang to be added to the forward primer
        backbone_overhang_r:str, default = 'cTAG'
            overhang to be added to the reverse primer
        cys4: str, default = 'gTTCACTGCCGTATAGGCAGCTAAGAAA'
            to be apended to the forward primers
    Attributes:
        amplicons : list
            list of amplicon objects
        amplicons_w_f_primer_overhang : list
            list of amplicon objects with forward primer overhangs
        amplicons_w_f_r_primer_overhang : list
            list of amplicon objects with forward and reverse primer overhangs

    Methods:
        simulate_pcrs()
            simulates PCR reactions and returns list of PCR products
        calculate_melting_temperature()
            calculates melting temperature of the primers and returns list of tuples
        make_pcr_df()
            makes a dataframe of PCR amplicon details
        make_primer_df()
            makes a dataframe of primer details
    """

    sgRNAs: list
    list_of_amplicons: list
    target_tm: int = 55
    tm_function: Callable = primer_tm_neb  # u can use your own

    restriction_overhang_f: str = "GATCGggtctcc"
    restriction_overhang_r: str = "GATCAGGTCTCg"

    backbone_overhang_f: str = "cATG"
    backbone_overhang_r: str = "cTAG"

    # to be appended to the forward primers
    cys4: str = "gTTCACTGCCGTATAGGCAGCTAAGAAA"

    # for post_init
    amplicons: list = field(init=False)
    amplicons_w_f_primer_overhang: list = field(init=False)
    amplicons_w_f_r_primer_overhang: list = field(init=False)

    def __post_init__(self) -> None:
        # generate the amplicons
        self.amplicons = make_amplicons(
            self.list_of_amplicons,
            target_tm=self.target_tm,
            tm_function=self.tm_function,
        )

        # make f primer overhangs
        self.amplicons_w_f_primer_overhang = make_GG_overhangs_forward_primers(
            self.amplicons,
            self.sgRNAs,
            restriction_overhang_f=self.restriction_overhang_f,
            backbone_overhang_f=self.backbone_overhang_f,
            cys4=self.cys4,
        )
        # make f and r primer overhangs
        self.amplicons_w_f_r_primer_overhang = make_GG_overhangs_reverse_primers(
            self.amplicons_w_f_primer_overhang,
            self.sgRNAs,
            backbone_overhang_r=self.backbone_overhang_r,
            restriction_overhang_r=self.restriction_overhang_r,
        )

    def simulate_pcrs(self):
        pcr_products = []
        for amplicon in self.amplicons_w_f_r_primer_overhang:
            pcr_product = pcr(
                amplicon.forward_primer.seq,
                amplicon.reverse_primer.seq,
                amplicon.template,
            )
            pcr_products.append(pcr_product)

        return pcr_products

    def calculate_melting_temperature(self):
        tm_list = []
        for amplicon in self.amplicons_w_f_r_primer_overhang:
            forward_tm = self.tm_function(str(amplicon.forward_primer.footprint))
            reverse_tm = self.tm_function(str(amplicon.reverse_primer.footprint))

            tm_list.append((forward_tm, reverse_tm))

        return tm_list

    def make_pcr_df(self):
        list_of_dicts = []
        for amplicon in self.amplicons_w_f_r_primer_overhang:
            # MAKE A DICT
            record = {
                "name": amplicon.name,
                # up repair
                "f_primer": str(amplicon.forward_primer.seq),
                "r_primer": str(amplicon.reverse_primer.seq),
                "f_tm": self.tm_function(str(amplicon.forward_primer.footprint)),
                "r_tm": self.tm_function(str(amplicon.reverse_primer.footprint)),
                "ta": primer_ta_neb(
                    str(amplicon.forward_primer.footprint),
                    str(amplicon.reverse_primer.footprint),
                ),
                "template": amplicon.template.name,
            }
            list_of_dicts.append(record)
        df = pd.DataFrame(list_of_dicts)

        return df

    def make_primer_df(self):
        # protospacer seq
        # Golden gate overhang seq

        list_of_dicts = []
        for amplicon in self.amplicons_w_f_r_primer_overhang:
            # MAKE A DICT
            record = {
                "name": amplicon.name,
                # up repair
                "f_primer": str(amplicon.forward_primer.seq),
                "r_primer": str(amplicon.reverse_primer.seq),
                "f_tm": self.tm_function(str(amplicon.forward_primer.footprint)),
                "r_tm": self.tm_function(str(amplicon.reverse_primer.footprint)),
                "f_primer_annealing": str(amplicon.forward_primer.footprint),
                "r_primer_annealing": str(amplicon.reverse_primer.footprint),
                "ta": primer_ta_neb(
                    str(amplicon.forward_primer.footprint),
                    str(amplicon.reverse_primer.footprint),
                ),
                "template": amplicon.template.name,
            }
            list_of_dicts.append(record)

        df = pd.DataFrame(list_of_dicts)

        return df


def digest_amplicons_w_BsaI(list_of_amplicons: list) -> None:
    """Generates pydna.amplicons which contains primers with a target temperature.

    Parameters
    ----------
    list_of_amplicons : list
        list of pydna.Dseqrecords


    Returns
    -------
    None
        The list of amplicons have been digested with BsaI.
    """

    list_of_digested_amplicons = []
    for amplicon in list_of_amplicons:
        vector_BsaI = Dseqrecord(amplicon.cut(BsaI)[1])
        list_of_digested_amplicons.append(vector_BsaI)

    return list_of_digested_amplicons

    # def retrieve_primers_from_golden_gate(primers,primers_footprint ):


#    '''Makes primers with footprints from two list of strings'''
#    primers_list = []
#    for i in range(len(primers)):
#        primers_list.append( Primer(primers[i], footprint=len(primers_footprint[i]) ))
#    return primers_list
