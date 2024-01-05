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

from pydna.dseqrecord import Dseqrecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from pydna.assembly import Assembly
from Bio.SeqRecord import SeqRecord
from pydna.design import primer_design
from pydna.design import assembly_fragments
from typing import List
from teemi.build.PCR import primer_tm_neb


def find_up_dw_repair_templates(
    genome: SeqRecord,
    repair_templates: List[str],
    target_tm: int = 65,
    primer_calc_function: callable = primer_tm_neb,
) -> List[dict]:
    """
    Find repair templates upstream and downstream of coding sequences in a genome.

    Parameters
    ----------
    genome : Bio.SeqRecord.SeqRecord
        A SeqRecord object containing the genome sequence.
    repair_templates : List[str]
        A list of locus tags of the coding sequences for which repair templates need to be found.
    target_tm : int, optional
        The target melting temperature for the PCR primers. Default is 65°C.
    primer_calc_function : callable, optional
        A function to calculate the melting temperature of primers. Default is the NEB formula.

    Returns
    -------
    List[dict]
        A list of dictionaries containing information about the repair templates. Each dictionary has the following keys:
            - name: The locus tag of the coding sequence for which the repair template was found.
            - up_repair: A Dseqrecord object containing the upstream repair template.
            - up_forwar_p: A Dseqrecord object containing the forward primer for the upstream repair template.
            - up_reverse_p: A Dseqrecord object containing the reverse primer for the upstream repair template.
            - tm_up_forwar_p: The melting temperature of the forward primer for the upstream repair template.
            - tm_up_reverse_p: The melting temperature of the reverse primer for the upstream repair template.
            - location_up_start: The starting location of the upstream repair template on the genome.
            - location_up_end: The ending location of the upstream repair template on the genome.
            - dw_repair: A Dseqrecord object containing the downstream repair template.
            - dw_forwar_p: A Dseqrecord object containing the forward primer for the downstream repair template.
            - dw_reverse_p: A Dseqrecord object containing the reverse primer for the downstream repair template.
            - tm_dw_forwar_p: The melting temperature of the forward primer for the downstream repair template.
            - tm_dw_reverse_p: The melting temperature of the reverse primer for the downstream repair template.
            - location_dw_start: The starting location of the downstream repair template on the genome.
            - location_dw_end: The ending location of the downstream repair template on the genome.
    """
    repair_DNA_templates = []
    for feature in genome.features:
        if feature.type == "CDS":
            # getting the locus tag if we see cds
            if feature.qualifiers["locus_tag"][0] in repair_templates:

                # get start and end locations
                start_location = int(feature.location.start)
                end_location = int(feature.location.end)

                # Fetch 500 upstream to start + end to 500 downstream
                repair_up = primer_design(
                    Dseqrecord(
                        str(genome.seq[start_location - 1000 : start_location]),
                        name=f"Repair_Template_UPSTREAM{feature.qualifiers['locus_tag'][0]}",
                    ),
                    target_tm=target_tm,
                    tm_func=primer_calc_function,
                )

                # Fetch 500 downstream to start + end to 500 downstream
                repair_dw = primer_design(
                    Dseqrecord(
                        str(genome.seq[end_location : end_location + 1000]),
                        name=f"Repair_Template_Downstream{feature.qualifiers['locus_tag'][0]}",
                    ),
                    target_tm=target_tm,
                    tm_func=primer_calc_function,
                )
                # MAKE A DICT
                record = {
                    "name": feature.qualifiers["locus_tag"][0],
                    # up repair
                    "up_repair": repair_up,
                    "up_forwar_p": repair_up.forward_primer,
                    "up_reverse_p": repair_up.reverse_primer,
                    "tm_up_forwar_p": primer_calc_function(
                        str(repair_up.forward_primer.seq)
                    ),
                    "tm_up_reverse_p": primer_calc_function(
                        str(repair_up.reverse_primer.seq)
                    ),
                    "location_up_start": start_location - 1000,
                    "location_up_end": start_location,
                    ""
                    # dw repair
                    "dw_repair": repair_dw,
                    "dw_forwar_p": repair_dw.forward_primer,
                    "dw_reverse_p": repair_dw.reverse_primer,
                    "tm_dw_forwar_p": primer_calc_function(
                        str(repair_dw.forward_primer.seq)
                    ),
                    "tm_dw_reverse_p": primer_calc_function(
                        str(repair_dw.reverse_primer.seq)
                    ),
                    "location_dw_start": end_location,
                    "location_dw_end": end_location + 1000,
                }

                repair_DNA_templates.append(record)

    return repair_DNA_templates


def update_primer_names(list_of_records: list) -> list:
    """Update primer names from dict.
    Parameters
    ----------
    list_of_records : list of dictionaries
        List of dictionaries where each dictionary contains information about a gene and its primers.

    Returns
    -------
    list of dictionaries
        The updated list of records with updated primer names.
    """
    # Making sure we only get unique sequences by checking if we saw a sequence before
    unique_primer_seqs = []
    unique_primers = []

    # changing names
    for i in range(len(list_of_records)):
        # UP forward primers
        if list_of_records[i]["up_forwar_p"].seq not in unique_primer_seqs:
            name = list_of_records[i]["gene_name"]
            list_of_records[i]["up_forwar_p"].name = f"{name}_up_F{i}"
            list_of_records[i]["up_forwar_p"].id = f"{name}_up_F{i}"
            list_of_records[i]["up_forwar_p_name"] = f"{name}_up_F{i}"

            unique_primers.append(list_of_records[i]["up_forwar_p"])
            unique_primer_seqs.append(list_of_records[i]["up_forwar_p"].seq)

        else:
            for seq in unique_primers:
                if seq.seq == list_of_records[i]["up_forwar_p"].seq:
                    list_of_records[i]["up_forwar_p_name"] = seq.id

        # UP reverse primers
        if list_of_records[i]["up_reverse_p"].seq not in unique_primer_seqs:
            name = list_of_records[i]["gene_name"]
            list_of_records[i]["up_reverse_p"].name = f"{name}_up_R{i}"
            list_of_records[i]["up_reverse_p"].id = f"{name}_up_R{i}"
            list_of_records[i]["up_reverse_p_name"] = f"{name}_up_R{i}"

            unique_primers.append(list_of_records[i]["up_reverse_p"])
            unique_primer_seqs.append(list_of_records[i]["up_reverse_p"].seq)
        else:
            for seq in unique_primers:
                if seq.seq == list_of_records[i]["up_reverse_p"].seq:
                    list_of_records[i]["up_reverse_p_name"] = seq.id

        # DW forward primers
        if list_of_records[i]["dw_forwar_p"].seq not in unique_primer_seqs:
            name = list_of_records[i]["gene_name"]
            list_of_records[i]["dw_forwar_p"].name = f"{name}_dw_F{i}"
            list_of_records[i]["dw_forwar_p"].id = f"{name}_dw_F{i}"
            list_of_records[i]["dw_forwar_p_name"] = f"{name}_dw_F{i}"

            unique_primers.append(list_of_records[i]["dw_forwar_p"])
            unique_primer_seqs.append(list_of_records[i]["dw_forwar_p"].seq)
        else:
            for seq in unique_primers:
                if seq.seq == list_of_records[i]["dw_forwar_p"].seq:
                    list_of_records[i]["dw_forwar_p_name"] = seq.id

        # DW reverse primers
        if list_of_records[i]["dw_reverse_p"].seq not in unique_primer_seqs:
            name = list_of_records[i]["gene_name"]
            list_of_records[i]["dw_reverse_p"].name = f"{name}_dw_R{i}"
            list_of_records[i]["dw_reverse_p"].id = f"{name}_dw_R{i}"
            list_of_records[i]["dw_reverse_p_name"] = f"{name}_dw_R{i}"

            unique_primers.append(list_of_records[i]["dw_reverse_p"])
            unique_primer_seqs.append(list_of_records[i]["dw_reverse_p"].seq)
        else:
            for seq in unique_primers:
                if seq.seq == list_of_records[i]["dw_reverse_p"].seq:
                    list_of_records[i]["dw_reverse_p_name"] = seq.id


def assemble_single_plasmid_with_repair_templates(
    repair_templates: list, digested_vector: Dseqrecord, overlap: int = 35
) -> list:
    """Assembles plasmids based on homology.

    Parameters
    ----------
    repair_templates : list of Dseqrecord
        List of Dseqrecord objects representing the repair templates to be assembled.
    digested_vector : Dseqrecord
        Dseqrecord object representing the digested vector.
    overlap : int, optional
        Length of the overlap between the repair templates. Default is 35.

    Returns
    -------
    list of Dseqrecord
        List of assembled Dseqrecord objects representing the repaired plasmids.
    """

    # initialize the circular vector
    all_parts = [digested_vector] + repair_templates + [digested_vector]

    # assemble them in a circular fashion - we actually do that by adding the digested vector first and last aboved
    new_vector = assembly_fragments(all_parts, overlap=overlap)

    return new_vector


def assemble_multiple_plasmids_with_repair_templates_for_deletion(
    list_of_gene_names: list,
    list_of_digested_plasmids: list,
    repair_DNA_templates: list,
    overlap=30,
) -> dict:
    """Assemble plasmids with repair templates via Gibson cloning.

    Parameters
    ----------
    list_of_gene_names : list
        A list of gene names.
    list_of_digested_plasmids : list
        A list of digested plasmids.
    repair_DNA_templates : list
        A list of templates to repair the plasmids
    overlap : int, optional
        The overlap between the repair templates and the digested plasmids. Default is 30.

    Returns
    -------
    dict
        A dictionary containing the assembled plasmids, primers and other information.
    """
    # initialize
    list_of_records = []

    # iterate through gene_names
    for gene_name in list_of_gene_names:
        # iterate through plasmids and find out if the gene name is in the name of the plasmid so we can assemble correctly
        for plasmid in list_of_digested_plasmids:
            plasmid_name = str(plasmid.name)
            x = plasmid_name.find(gene_name)
            if x != -1:
                # we got a match - get repair templates
                for repair_template_name in repair_DNA_templates:
                    # double check that the names are working
                    if repair_template_name["name"] == gene_name:
                        # making a list of the repair templates
                        repair_templates = [
                            repair_template_name["up_repair"],
                            repair_template_name["dw_repair"],
                        ]

                        # Assemble vector with the
                        assembled_vector = (
                            assemble_single_plasmid_with_repair_templates(
                                repair_templates, plasmid, overlap=overlap
                            )
                        )

                        # Make this to a contig and make genbank files of them
                        assemblyobj = Assembly(assembled_vector)
                        contig = Dseqrecord(assemblyobj.assemble_circular()[0])

                        ## change features
                        contig.features.append(
                            SeqFeature(
                                FeatureLocation(len(plasmid), len(plasmid) + 1000),
                                type="CDS",
                                qualifiers={"label": f"UP_repair{gene_name}"},
                            )
                        )
                        contig.features.append(
                            SeqFeature(
                                FeatureLocation(
                                    len(plasmid) + 1000, len(plasmid) + 2000
                                ),
                                type="CDS",
                                qualifiers={"label": f"DW_repair{gene_name}"},
                            )
                        )

                        # Retrieve information
                        # MAKE A DICT
                        record = {
                            "gene_name": gene_name,
                            "name": plasmid_name,
                            "contig": contig,
                            # up repair
                            "up_forwar_p_name": assembled_vector[1].forward_primer.id,
                            "up_reverse_p_name": assembled_vector[1].reverse_primer.id,
                            "up_forwar_p": assembled_vector[1].forward_primer,
                            "up_reverse_p": assembled_vector[1].reverse_primer,
                            "up_forwar_p_anneal": assembled_vector[
                                1
                            ].forward_primer.footprint,
                            "up_reverse_p_anneal": assembled_vector[
                                1
                            ].reverse_primer.footprint,
                            "tm_up_forwar_p": primer_tm_neb(
                                str(assembled_vector[1].forward_primer.footprint)
                            ),
                            "tm_up_reverse_p": primer_tm_neb(
                                str(assembled_vector[1].reverse_primer.footprint)
                            ),
                            # dw repair
                            "dw_forwar_p_name": assembled_vector[2].forward_primer.id,
                            "dw_reverse_p_name": assembled_vector[2].reverse_primer.id,
                            "dw_forwar_p": assembled_vector[2].forward_primer,
                            "dw_reverse_p": assembled_vector[2].reverse_primer,
                            "dw_forwar_p_anneal": assembled_vector[
                                2
                            ].forward_primer.footprint,
                            "dw_reverse_p_anneal": assembled_vector[
                                2
                            ].reverse_primer.footprint,
                            "tm_dw_forwar_p": primer_tm_neb(
                                str(assembled_vector[2].forward_primer.footprint)
                            ),
                            "tm_dw_reverse_p": primer_tm_neb(
                                str(assembled_vector[2].reverse_primer.footprint)
                            ),
                        }
                        list_of_records.append(record)

    return list_of_records
