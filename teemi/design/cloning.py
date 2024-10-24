#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)
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
import pydna
import Bio
import Bio.SeqFeature
import Bio.SeqRecord
import Bio.SeqIO
from Bio.Seq import Seq
from pydna.assembly import Assembly
from pydna.dseq import Dseq
from math import fabs
import re


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
    1. up : pydna.dseqrecord.
        Sequence upstream of the DSB: pydna.dseqrecord.
    2. dw : pydna.dseqrecord.
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
        Bio.SeqFeature.FeatureLocation(0, len(up)), type="misc_feature")
    up_feature.qualifiers["label"] = up.name
    up.features.append(up_feature)

    dw_feature = Bio.SeqFeature.SeqFeature(
        Bio.SeqFeature.FeatureLocation(0, len(dw)), type="misc_feature")
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
            )
        else:
            template_amplification_site.annotations["batches"] = []
            template_amplification_site.annotations["batches"].append(
                template.annotations
            )

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
            if str(feature.qualifiers["name"][0]) == anno:
                site = template[feature.location.start : feature.location.end]
                site.name = name
                site.annotations = template.annotations

                sites.append(site)
    return sites


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
    rv_U_idx = amplicon.reverse_primer.seq.find("U")

    digested_watson = amplicon.seq.watson[fw_U_idx + 1 :]
    digested_crick = amplicon.seq.crick[rv_U_idx + 1 :]

    digested_pcr = pydna.dseqrecord.Dseqrecord(
        pydna.dseq.Dseq(watson=digested_watson, crick=digested_crick),
        features=amplicon.features,
        annotations=amplicon.annotations,
    )
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
            Bio.SeqFeature.FeatureLocation(0, len(assembly)),
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


def plate_plot(df, value):
    """
    Plots a 96 well plate as a pandas df.

    Parameters
    ----------
    df : pd.Dataframe
        A pandas dataframe.
    value : str
        The name of the pandas dataframe column that you want to display.

    Returns
    -------
    pd.Dataframe
        in a 96 well plate format of the chosen column.

    Example
    -------
    .. code-block:: python

        # Initialize:
        Amplicon_df = {
            'name': ['PCR_G8H_01', 'PCR_G8H_05', ...],
            'location': ['l5_A03', 'l5_A07', ...],
            ...
        }

        # Call the function:
        plate_plot(amplicon_df, 'name')

    .. code-block:: none

        name
        pcol	1	2	3	4	5	6	7	8	9	10	11	12
        prow
        A	PCR_G8H_01	PCR_G8H_05	PCR_G8H_09	PCR_G8H_13	PCR_G8H_17	...
        ...

    """
    cols = [value, "prow", "pcol"]
    return df[cols].set_index(["prow", "pcol"]).unstack(level=-1)


def seq_to_annotation(
    seq_record_from: Bio.SeqRecord, seq_record_onto: Bio.SeqRecord, type_name: str
):
    """Anotate an Bio.SeqRecord object from another
    bio.seqrecord object.

    Parameters
    ----------
    seqrec_from: Bio.SeqRecord
        annotation sequence that will be extracted
    seqrec_onto: Bio.SeqRecord
    type_name: str
        name of the sequence that will be extracted

    Returns
    -------
    None
    """
    match_index = find_sequence_location(seq_record_from, seq_record_onto)

    if match_index[2] > 0:
        start_location, end_location, strand = (
            match_index[0],
            match_index[1],
            match_index[2],
        )
    else:
        start_location, end_location, strand = (
            match_index[1],
            match_index[0],
            match_index[2],
        )

    feature = Bio.SeqFeature.SeqFeature(
        Bio.SeqFeature.FeatureLocation(start_location, end_location),
        type=type_name)

    feature.qualifiers["label"] = seq_record_from.id
    seq_record_onto.features.append(feature)


def find_sequence_location(
    sequence: Bio.SeqRecord, sequence_to_search_in: Bio.SeqRecord
) -> tuple:
    """Finds start and end location of a mathced sequence.

    Parameters
    ----------
    sequence : str
    sequence_to_search_in : Bio.SeqRecord

    Returns
    -------
    (start_index,end_index) : tuple
    """
    strand = +1
    start_index = sequence_to_search_in.seq.find(sequence.seq)
    end_index = start_index + len(sequence)

    if start_index == -1:
        # search reverse_comp
        start_index = len(
            sequence_to_search_in
        ) - sequence_to_search_in.reverse_complement().seq.find(sequence.seq)
        end_index = start_index - len(sequence)
        strand = -1

        if start_index == -1:
            raise ValueError("ValueERROR - couldnt find a match")

    return (start_index, end_index, strand)


def crispr_db_break_location(start_location, end_location, strand):
    """
    Determine the CRISPR cut location in the genome.

    Parameters
    ----------
    start_location : int
        Start position of the sgRNA sequence in the chromosome.
    end_location : int
        End position of the sgRNA sequence in the chromosome.
    strand : int
        Strand of the sgRNA sequence in the chromosome, +1 for positive strand, -1 for negative strand.

    Returns
    -------
    crispr_db_break : int
        CRISPR cut location in the genome.
    """
    if strand == +1:
        crispr_db_break = start_location + 17
    if strand == -1:
        crispr_db_break = start_location - 3

    return crispr_db_break


def add_feature_annotation_to_seqrecord(
    sequence: Bio.SeqRecord, label="", type_name="misc_feature", strand=0
) -> None:
    """Adds feature, label and name to a Bio.Seqrecord sequence.
    Parameters
    ----------
    sequence : Bio.SeqRecord
    label : str (optional)
    type_name : str (default: "misc_feature")
    strand : int (default 0)

    Returns
    -------
    None
    """
    bio_feature = Bio.SeqFeature.SeqFeature(
        Bio.SeqFeature.FeatureLocation(0, len(sequence)), type=type_name)

    # label
    sequence.features.append(bio_feature)
    sequence.features[0].qualifiers["label"] = label
    sequence.features[0].qualifiers["name"] = sequence.name


def find_all_occurrences_of_a_sequence(
    sequence: Bio.SeqRecord, sequence_to_search_in: Bio.SeqRecord
) -> tuple:
    """
    Searches for all occurrences of a given sequence in a given string.

    Parameters
    ----------
    sequence : Bio.SeqRecord
        Sequence to search for.
    sequence_to_search_in : Bio.SeqRecord
        Sequence to search in.

    Returns
    -------
    tuple
        Number of occurrences of `sequence` in `sequence_to_search_in`.

    """
    finder = re.finditer(
        str(sequence.seq.upper()), str(sequence_to_search_in.seq.upper())
    )
    matches_watson = [(match.start(), match.end()) for match in finder]

    if len(matches_watson) < 2:
        finder = re.finditer(
            str(sequence.seq.upper()),
            str(sequence_to_search_in.seq.reverse_complement().upper()),
        )
        matches_crick = [(match.start(), match.end()) for match in finder]

        return len(matches_watson) + len(matches_crick)
    else:
        return len(matches_watson)


