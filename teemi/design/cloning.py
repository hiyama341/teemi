#!/usr/bin/env python
# MIT License
# Copyright (c) 2022, Technical University of Denmark (DTU)
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


""" Module used for cloning of microbial strains."""
import functools
import pandas as pd
import pydna
import Bio
import Bio.SeqFeature
import Bio.SeqRecord
import Bio.SeqIO
from Bio.Seq import Seq
from pydna.assembly import Assembly
from pydna.dseq import Dseq


def CAS9_cutting(gRNA_record, background_record):

    """Simulates double-stranded-break by CAS9 given a gRNA.

    Parameters
    ----------
    gRNA_record: pydna.dseqrecord.
        A 20 bp DNA sequence
    background_record: pydna.dseqrecord.
        The sequence of interest for CRISPR mediated DSB

    Returns
    -------
    1.pydna.dseqrecord.
        Sequence upstream of the DSB: pydna.dseqrecord.
    2. pydna.dseqrecord.
        Sequence downstream of the DSB: pydna.dseqrecord.

    """
    gRNA_sequence = gRNA_record.seq.watson.upper()

    background_sequence = background_record.seq.upper()

    gRNA_strand = 1

    if background_sequence.find(gRNA_sequence) == -1:

        gRNA_strand = -1
        gRNA_sequence = Seq(gRNA_sequence)

        gRNA_sequence = (gRNA_sequence).reverse_complement()
        if background_sequence.find(gRNA_sequence) == -1:
            print("not on -1, CAN'T FIND THE CUT SITE IN YOUR SEQUENCE")

    if gRNA_strand == 1:
        cut_pos_rel_start = 17
    else:
        cut_pos_rel_start = 3

    gRNA_start = background_sequence.find(gRNA_sequence)

    cut_site = gRNA_start + cut_pos_rel_start

    up = background_record[0:cut_site]
    dw = background_record[cut_site : len(background_sequence)]

    up.name = "UP" + "_" + gRNA_record.name + "_" + background_record.name
    dw.name = "DW" + "_" + gRNA_record.name + "_" + background_record.name

    up_feature = Bio.SeqFeature.SeqFeature(
        Bio.SeqFeature.FeatureLocation(0, len(up)), type="misc_feature", strand=+1
    )
    up_feature.qualifiers["label"] = up.name
    up.features.append(up_feature)

    dw_feature = Bio.SeqFeature.SeqFeature(
        Bio.SeqFeature.FeatureLocation(0, len(dw)), type="misc_feature", strand=+1
    )
    dw_feature.qualifiers["label"] = dw.name
    dw.features.append(dw_feature)

    up = pydna.dseqrecord.Dseqrecord(
        up,
    )
    dw = pydna.dseqrecord.Dseqrecord(dw)

    # UPS more than one cut site?
    if dw.seq.find(gRNA_sequence) != -1:
        print("OBS", gRNA_sequence, "cuts more than one location!")

    return (up, dw)


def extract_gRNAs(template, name):
    """Extracts gRNAs from a template.

    Parameters
    ----------
    template: pydna.dseqrecord or pydna.amplicon.Amplicon
        a plasmid or piece of DNA

    name: str
        a string that would include the feature name for example: gRNA

    Returns
    -------
    list of pydna.dseqrecord or pydna.amplicon.Amplicon
        list of with the found features and their sequences

    """
    gRNAs = []
    for feature in template.features:
        if name in feature.qualifiers.get("name", ""):

            gRNA = template[feature.location.start : feature.location.end]
            gRNA.name = feature.qualifiers.get("name", "")
            gRNAs.append(gRNA)

    return gRNAs


def remove_features_with_negative_loc(record):
    """Removes a SeqFeatures if negative.

    Parameters
    ----------
    record: pydna.amplicon.Amplicon.
        A amplicon with SeqFeature and locations

    Returns
    -------
    record: pydna.amplicon.Amplicon.
        With the negative features deleted

    """

    for i in range(len(record.features)):
        if record.features[i].location.start < 0 or record.features[i].location.end < 0:
            del record.features[i]


def extract_template_amplification_sites(templates, names, terminator):
    """Extracts amplifications sites from a templates features


    Parameters
    ----------
    templates: list of Bio.SeqRecord.SeqRecord
        list of Bio.SeqRecord.SeqRecord objects with SeqFeatures

    names: list of strings
        list of strings to be extracted

    terminator: str
        a string with the name of upstream terminator


    Returns
    -------
    record: list of Bio.SeqRecord.SeqRecord
        list of extracted elements

    """
    template_amplification_sites = []
    for name, template in zip(names, templates):
        for feature in template.features:
            if feature.qualifiers["name"] == name:

                CDS_strand = feature.strand
                if CDS_strand == 1:
                    start = feature.location.start
                    print(start)
                else:
                    end = feature.location.end
                    print(end)

            if feature.qualifiers["name"].endswith(terminator):
                terminator_strand = feature.strand
                if terminator_strand == 1:
                    end = feature.location.end
                else:
                    start = feature.location.start
                    print(start)

        template_amplification_site = template[start:end]

        if "batches" in template_amplification_site.annotations.keys():
            template_amplification_site.annotations["batches"].append(
                template.annotations["batches"][0]
            )  # {'location': template.annotations['batches'][0]['location']}
        else:
            template_amplification_site.annotations["batches"] = []
            template_amplification_site.annotations["batches"].append(
                template.annotations
            )  # {'location': template.annotations['batches'][0]['location']}

        template_amplification_sites.append(template_amplification_site)
    return template_amplification_sites


def extract_sites(annotations, templates, names):
    """This function extracts the sequences from annotated sequences based
    on their names

    Parameters
    ----------
    annotations: list
        list of annotations for sequences that will be extracted

    templates: list of Bio.SeqRecord.SeqRecord
        A list of Bio.SeqRecord.SeqRecord with SeqFeatures

    names: str
        name of the sequence that will be extracted

    Returns
    -------
    record: list of Bio.SeqRecord.SeqRecord
        list of extracted sites
    """

    sites = []
    for anno, template, name in zip(annotations, templates, names):
        for feature in template.features:
            # site = template[feature.location.start : feature.location.end]
            if str(feature.qualifiers["name"][0]) == anno:
                site = template[feature.location.start : feature.location.end]
                site.name = name

                # # If there is an batch anotation we can save it
                # if "batches" in site.annotations.keys():
                #     site.annotations["batches"].append(
                #         template.annotations["batches"][0]
                #     )
                # else:
                #     site.annotations["batches"] = []
                #     site.annotations["batches"].append(template.annotations)
                site.annotations = template.annotations

                sites.append(site)
    return sites


def seq_to_annotation(seqrec_from, seqrec_onto, aType):
    """Anotate an amplicon object from another amplicon object.

    Parameters
    ----------
    seqrec_from: str
        annotation sequence that will be extracted

    seqrec_onto: list of Bio.SeqRecord.SeqRecord
        A list of Bio.SeqRecord.SeqRecord with SeqFeatures

    aType: str
        name of the sequence that will be extracted

    Returns
    -------
    record: list of Bio.SeqRecord.SeqRecord
        list of extracted sites
    """

    seq_from = seqrec_from.seq.watson.upper()
    seq_onto = seqrec_onto.seq.watson.upper()

    strand = 1
    match_index = seq_onto.find(seq_from)

    # if there is match
    if match_index != -1:
        start = match_index
        end = start + len(seq_from)
    else:
        # If no match we look at the reverse complement
        seq_onto = Seq(seq_onto)
        rev_match_index = seq_onto.reverse_complement().find(seq_from)

        # if we get a match here
        if rev_match_index != -1:
            strand = -1
            reclength = len(str(seqrec_onto.seq))
            end = reclength - rev_match_index

            start = end - len(seq_from)

        else:
            print(
                "no match! seq:"
                + str(seqrec_from.name)
                + "\nnot annealing to:"
                + str(seqrec_onto.name)
            )

    # add the feature to the amplicon
    feature = Bio.SeqFeature.SeqFeature(
        Bio.SeqFeature.FeatureLocation(start, end), type=aType, strand=strand
    )
    feature.qualifiers["label"] = seqrec_from.id
    seqrec_onto.features.append(feature)


def USER_enzyme(amplicon):
    """Simulates digestion with USER enzyme.

    Parameters
    ----------
    amplicon: pydna.amplicon.Amplicon
        An pydna.amplicon.Amplicon to with Uracil integrated

    Returns
    -------
    Dseqrecord
        USER digested Dseqrecord with USER tails
    """
    fw_U_idx = amplicon.forward_primer.seq.find("U")
    # fw_U_idx

    rv_U_idx = amplicon.reverse_primer.seq.find("U")
    # rv_U_idx

    digested_watson = amplicon.seq.watson[fw_U_idx + 1 :]
    # digested_watson

    digested_crick = amplicon.seq.crick[rv_U_idx + 1 :]
    # digested_crick

    digested_pcr = pydna.dseqrecord.Dseqrecord(
        pydna.dseq.Dseq(watson=digested_watson, crick=digested_crick),
        features=amplicon.features,
        annotations=amplicon.annotations,
    )
    # digested_pcr.seq
    return digested_pcr


def nicking_enzyme(vector):
    """Nt.Bbc.CI (nicking enzyme, Nicks) a vector with the
    sequence 'CGCGTG' on watson and 'CGCACG' on crick strand.
    Parameters
    ----------
    vector: Dseq
        digested Dseqrecord - usually with AsiSI or similar overhang

    Returns
    -------
    Dseq with nick - ready for USER cloning
    """
    if vector.seq[0:8].watson == "CGCGTG" and vector.seq.crick[:6] == "CGCACG":
        return Dseq(watson=vector.seq.watson[6:], crick=vector.seq.crick[6:], ovhg=8)
    else:
        print("No nicking sequnce")


def casembler(
    bg_strain,
    site_names=None,
    gRNAs=None,
    parts=None,
    assembly_limits=None,
    assembly_names=None,
    verbose=False,
    to_benchling=False,
):

    """Simulate in vivo assembly and integration with the possibility
    of printing to gb files or send it directly to benchling.

    Parameters
    ----------
    bg_strain : GenBank
        strain of choice eg. genbank file

    site_names: list
          list of names e.g. [X-3, XI-3]

    gRNAs: Seqrecords
        list of 20 bp seqrecords e.g. [ATF1_gRNA, CroCPR_gRNA]

    parts: list
        list of list of parts e.g. [[ATF1_repair_template],[CPR_repair_template]]

    assembly_limits: list
        list of numbers of bp assembly limits e.g. [200,400]

    assembly_names: list
        list of names of DNA post assembly e.g.["X_3_tADH1_P2_pPGK1", "XI_3_UP_DW"]

    verbose: bool
        write DNA e.g. False

    to_benchling: bool
        upload DNA to benchling e.g. False

    Returns
    -------
    One dseqrecord
        of assembled contig

    """
    assemblies = []
    for int_no in range(0, len(gRNAs)):

        gRNA = gRNAs[int_no]

        site_name = site_names[int_no]

        site = extract_sites([site_name], [bg_strain], [site_name])[0]

        UP, DW = CAS9_cutting(gRNA, site)

        fragments = [UP] + parts[int_no] + [DW]

        assembly = pydna.assembly.Assembly(
            fragments, limit=assembly_limits[int_no]
        ).assemble_linear()[0]

        # sometimes pydna.assembly.Assembly distorts the start, end location of features to become negative which produces and error when printing.
        # This function is created as a workaround.
        # CPR assembly gives DW_XI_3 annotation called "DW_XI_3" a negative start location.
        # A quick workaround is to remove featuere
        remove_features_with_negative_loc(assembly)

        assembly.name = assembly_names[int_no]

        assembly_feat = Bio.SeqFeature.SeqFeature(
            Bio.SeqFeature.FeatureLocation(0, len(assembly), strand=1),
            type="misc_feature",
        )
        assembly_feat.qualifiers["name"] = site_names[int_no]
        assembly_feat.qualifiers["label"] = site_names[int_no]
        assembly.features.append(assembly_feat)

        if verbose:
            DNAs = [UP] + [assembly] + [DW]
            for DNA in DNAs:
                DNA.write("./" + DNA.name + ".gb")  # "../data/processed/"

        if to_benchling:
            to_benchling(assembly, "to_benchling")

        assemblies.append(assembly)

    return functools.reduce(lambda x, y: x + y, assemblies)


def UPandDW(strain, isite_name, path_to_gRNA_table="../data/raw/gRNAtable.csv"):
    """Finds upstream and downstream sequences based on genome and site name.

    Parameters
    ----------
    strain : str
        name of the strain eg. CENPK113-7d
        (you should specify path to the chromosome)

    isite_name : str
        a string of the site chomosomal site you want to retrieve

    Returns
    -------
    UP_sites : list
        list of pydna.dseqrecord or pydna.amplicon.Amplicon

    DW_sites : list
        list of pydna.dseqrecord or pydna.amplicon.Amplicon

    """

    # load lookup table
    gRNAtable = pd.read_csv(path_to_gRNA_table, index_col="name")

    chromosome_no = gRNAtable.loc[isite_name, "chromosome"]

    # load chromosome
    PathToChromosomeSeq = (
        "../data/raw/" + strain + "/" + str(chromosome_no).zfill(2) + ".fa"
    )
    ChromosomeSeq = Bio.SeqIO.read(PathToChromosomeSeq, "fasta").seq

    # define homology region location with respect to gRNA sequence
    # f_hom is the length of the UP homology
    # e_hom is the length of the DW homology
    # f_dist is distance from end of UP homology to the first base in isite_sequence
    # e_dist is distance from end of the isite_sequence to DW homology
    f_dist = gRNAtable.loc[isite_name, "f_dist"]
    e_dist = gRNAtable.loc[isite_name, "e_dist"]
    f_hom = gRNAtable.loc[isite_name, "f_hom"]
    e_hom = gRNAtable.loc[isite_name, "e_hom"]

    isite_sequence = gRNAtable.loc[isite_name, "sequence"]
    isite_sequence = Bio.Seq.Seq(isite_sequence)

    # Determine gRNA sequence strand
    gRNA_strand = 1
    if ChromosomeSeq.find(isite_sequence) == -1:
        print("not on +1")
        gRNA_strand = -1
        isite_sequence = isite_sequence.reverse_complement()
        if ChromosomeSeq.find(isite_sequence) == -1:
            print("not on -1")
            print("CAN'T FIND THE CUT SITE IN YOUR SEQUENCE")

    # Locate UP and DW
    StartIndex = ChromosomeSeq.find(isite_sequence)
    EndIndex = StartIndex

    UPseq = ChromosomeSeq[StartIndex + f_dist - f_hom : StartIndex + f_dist]
    DWseq = ChromosomeSeq[EndIndex + e_dist : EndIndex + e_dist + e_hom]

    UPrec = Bio.SeqRecord.SeqRecord(UPseq, name=isite_name + "UP")
    DWrec = Bio.SeqRecord.SeqRecord(DWseq, name=isite_name + "DW")

    # Annotate
    UP_feature = Bio.SeqFeature.SeqFeature(
        Bio.SeqFeature.FeatureLocation(0, len(UPseq)), type="misc_feature", strand=+1
    )
    UP_feature.qualifiers["label"] = UPrec.name
    UPrec.features.append(UP_feature)

    DW_feature = Bio.SeqFeature.SeqFeature(
        Bio.SeqFeature.FeatureLocation(0, len(DWseq)), type="misc_feature", strand=+1
    )
    DW_feature.qualifiers["label"] = DWrec.name
    DWrec.features.append(DW_feature)

    return ([UPrec], [DWrec])


def multiply_list(myList):
    """Multiplies elements one by one.

    Parameters
    ----------
    myList: list
        list of integers to be multiplied

    Returns
    -------
    result : int

    """
    result = 1
    for x in myList:
        result = result * x
    return result


def remove_tuple_duplicates(lst: list) -> list:
    """Removes tuple duplicates

    Parameters
    ----------
    lst: list
        list with duplicated elements

    Returns
    -------
    list
        without duplicates
    """
    return [t for t in (set(tuple(i) for i in lst))]


def recs_no_duplicates(recs_with_duplicates: list) -> list:
    """Removes duplicate sequences from a list.

    Parameters
    ----------
    recs_with_duplicates: list
        list with duplicated elements

    Returns
    -------
    list
        without duplicates
    """
    seen_sequences = set()
    recs_no_dup = []
    for rec in recs_with_duplicates:
        if rec.seq.watson not in seen_sequences:
            recs_no_dup.append(rec)
            seen_sequences.add(rec.seq.watson)
    return recs_no_dup


def recs_no_duplicates_names(recs_with_duplicates):
    """Removes duplicate names from a list

    Parameters
    ----------
    recs_with_duplicates: list
        list with duplicated elements

    Returns
    -------
    list
        without duplicates


    """
    seen_names = set()
    recs_no_dup = []
    for rec in recs_with_duplicates:
        if rec.name not in seen_names:
            recs_no_dup.append(rec)
            seen_names.add(rec.name)
    return recs_no_dup


def plate_plot(df, value):
    """Plots a 96 well plate as a pandas df.

    Parameters
    ----------
    df: pd.Dataframe
        A pandas dataframe with

    value: pandas dataframe column name
        The name of the pandas dataframe coloumn that you want to display
    Returns
    -------
    pd.Dataframe
        in a 96 well plate format of the chosen column

    Example
    -------

    1. Initialize:
    Amplicon_df = {
        name	location	template_name	fw_name	fw_location	rv_name	rv_location	prow	pcol
    29	PCR_G8H_01	l5_A03	VminG8H_tADH1	PR_G8H_01	op4_A10	PR_G8H_02	op4_A01	A	1
    25	PCR_G8H_05	l5_A07	SmusG8H_tADH1	PR_G8H_01	op4_A10	PR_G8H_06	op4_A02	A	2
    21	PCR_G8H_09	l5_B02	RsepG8H_tADH1	PR_G8H_01	op4_A10	PR_G8H_10	op4_A03	A	3
    17	PCR_G8H_13	l5_B06	CacuG8H_tADH1	PR_G8H_01	op4_A10	PR_G8H_14	op4_A04	A	4
    13	PCR_G8H_17	l5_C01	OpumG8H_tADH1	PR_G8H_01	op4_A10	PR_G8H_18	op4_A05	A	5
    }

    2. Call the function:  # here we call the name coloumn
    plate_plot(amplicon_df, 'name')

    <<Result:
    name
    pcol	1	2	3	4	5	6	7	8	9	10	11	12
    prow
    A	PCR_G8H_01	PCR_G8H_05	PCR_G8H_09	PCR_G8H_13	PCR_G8H_17	PCR_G8H_21	PCR_G8H_25	PCR_G8H_29	PCR_G8H_33	PCR_UP_tADH1_01	PCR_PRO_01	NaN1
    B	PCR_G8H_02	PCR_G8H_06	PCR_G8H_10	PCR_G8H_14	PCR_G8H_18	PCR_G8H_22	PCR_G8H_26	PCR_G8H_30	PCR_G8H_34	PCR_TRP1-DW_02	PCR_PRO_02	NaN2
    C	PCR_G8H_03	PCR_G8H_07	PCR_G8H_11	PCR_G8H_15	PCR_G8H_19	PCR_G8H_23	PCR_G8H_27	PCR_G8H_31	PCR_G8H_35	PCR_TRP1-DW_01	PCR_PRO_03	NaN3
    D	PCR_G8H_04	PCR_G8H_08	PCR_G8H_12	PCR_G8H_16	PCR_G8H_20	PCR_G8H_24	PCR_G8H_28	PCR_G8H_32	PCR_G8H_36	NaN4	PCR_PRO_04	NaN5
    E	PCR_CPR_01	PCR_CPR_10	PCR_CPR_03	PCR_CPR_09	PCR_CPR_02	PCR_CPR_06	PCR_CPR_07	PCR_CPR_08	PCR_CPR_04	PCR_CPR_05	PCR_PRO_05	NaN6
    F	PCR_CPR_11	PCR_CPR_20	PCR_CPR_13	PCR_CPR_19	PCR_CPR_12	PCR_CPR_16	PCR_CPR_17	PCR_CPR_18	PCR_CPR_14	PCR_CPR_15	PCR_PRO_06	NaN7
    G	PCR_CPR_21	PCR_CPR_30	PCR_CPR_23	PCR_CPR_29	PCR_CPR_22	PCR_CPR_26	PCR_CPR_27	PCR_CPR_28	PCR_CPR_24	PCR_CPR_25	PCR_PRO_07	NaN8
    H	PCR_CPR_31	PCR_CPR_40	PCR_CPR_33	PCR_CPR_39	PCR_CPR_32	PCR_CPR_36	PCR_CPR_37	PCR_CPR_38	PCR_CPR_34	PCR_CPR_35	PCR_PRO_08	NaN9

    """

    cols = [value, "prow", "pcol"]
    return df[cols].set_index(["prow", "pcol"]).unstack(level=-1)


def CRIPSR_knockout(gRNA_record, insertion_site, repair_DNA):

    """Simple version of casembler - Cuts the insertion site with
     CAS9_cutting and assembles knockout with a repair template.

    Parameters
    ----------
    gRNA_record: pydna.dseqrecord.
        A 20 bp DNA sequence

    insertion_site: pydna.dseqrecord.
        The site to knock out

    repair_DNA: pydna.dseqrecord.
        Repair template. Typucally 90 bp or longer

    Returns
    -------
    pydna.dseqrecord.
        Of assembled contig after CRISPR-mediated KO
    """
    # Create fragments after CAS9 cut
    IS_UP, IS_DW = CAS9_cutting(gRNA_record, insertion_site)

    # create list of parts and assemble to knockout sequence
    assmeble_parts = IS_UP, repair_DNA, IS_DW
    assembled_knockout = Assembly(assmeble_parts).assemble_linear()[0]

    return assembled_knockout
