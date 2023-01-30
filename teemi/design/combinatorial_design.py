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

""" This part of the design module is used for making combinatorial libraries from DNA fragments."""

# standard libraries
import itertools
import numpy as np
import pandas as pd

# Pydna for the molecular bio
from pydna.design import primer_design
from pydna.design import assembly_fragments
from pydna.assembly import Assembly
from pydna.tm import tm_default as _tm_default


def combinatorial_list_maker(
    listOflist_that_is_being_made_into_all_combinations: list,
) -> list:
    """Makes all possible combinations from a list of list.

    Parameters
    ----------
    listOflist_that_is_being_made_into_all_combinations: list[list[any_type]]
        can be of any type inside the list of lists

    Returns
    -------
    combinations: list[tuple(any_type)]
        all possible combinations of the list of lists

    """
    combinations = list(
        itertools.product(*listOflist_that_is_being_made_into_all_combinations)
    )

    return combinations


def systematic_names_function(List_of_list_parts: list) -> list:
    """Returns a list of list with systematic names i.e [1,1,1], [1,2,1]... etc

    Parameters
    ----------
    List_of_list_parts: list of list
        can have anny type withing the list[list[any_type]]

    Returns
    -------
    combinatorial_list_of_indexes
        list of tuples with the systematic names eg. [(1,1,1),(1,2,1)]

    """
    # The number of parts of each fragment
    no_parts = [int(len(l)) for l in List_of_list_parts]

    ### For naming the strains systematically ### basicly making a list from the number of parts with indexes
    list_of_systematic = []
    midlertidiglist = []
    for parts in no_parts:
        for j in range(0, parts):
            midlertidiglist.append(j + 1)
        list_of_systematic.append(midlertidiglist)
        midlertidiglist = []

    # Then we use itertools to make the right combinations
    combinatorial_list_of_indexes = list(itertools.product(*list_of_systematic))

    return combinatorial_list_of_indexes


def empty_list_maker(list_of_sequences: list):
    """returns empty list in the length of seqs

    Parameters
    ----------
    list_of_sequences: list
        could be any list with any types

    Returns
    -------
    EmptyList:list
        an empty list with the same dimensions

    """
    EmptyList = [[] for i in range(len(list_of_sequences))]

    return EmptyList


def simple_amplicon_maker(
    list_of_seqs: list, list_of_names: list, target_tm=56.0, limit=13
):
    """Creates amplicons, updates their names

    Parameters
    ----------
    list_of_seqs : list[list[pydna.dseqrecord.Dseqrecord]]
        List of the pydna.dseqrecord import Dseqrecord elements u want to made into amplicons

    list_of_names : list[list[str]]
        provide names for the sequences since pydna changes their names to amplicon

    Returns
    -------
    list_of_amplicons : list[pydna.amplicon.Amplicon]
        list with the pydna.amplicon.Amplicon objects that have been made

    list_of_amplicon_primers : list[list[(pydna.seq.Seq, pydna.seq.Seq)]]
        a list of all the generated primers in tuples where index0 = forward primer
        and index1=reverse primer. Both are pydna.seq.Seq objects

    list_of_amplicon_primer_temps : list[list[(float, float)]]
        a list of melting temperatures in tuples where index0 = forward primer melting temp
        and index1=reverse primer melting temp.

    """
    # Start by making an empty list
    list_of_amplicons = [[] for i in range(len(list_of_seqs))]
    list_of_amplicon_primers = [[] for i in range(len(list_of_seqs))]
    list_of_amplicon_primer_temps = [[] for i in range(len(list_of_seqs))]

    ### HERE WE CALCULATE Amplicons, primers, and their temperatures
    # Then we calculate the primers with the NEB calculator
    for i in range(0, len(list_of_seqs)):
        for j in range(0, len(list_of_seqs[i])):
            # Append Amplicons
            amplicons = primer_design(
                list_of_seqs[i][j],
                tm_func=_tm_default,
                target_tm=target_tm,
                limit=limit,
            )  ############## Can add NEB Calculator here: primer_TM ################# _tm_default i.e tm_func = _tm_default,

            # Updating names
            amplicons.name = list_of_names[i][j]
            list_of_amplicons[i].append(amplicons)

            # Save the primers
            primers = (amplicons.forward_primer.seq, amplicons.reverse_primer.seq)
            list_of_amplicon_primers[i].append(primers)

            # Save melting temps
            ############## Can add NEB Calculator here: primer_TM #############################
            melting_temps = (
                _tm_default(amplicons.forward_primer.seq),
                _tm_default(amplicons.reverse_primer.seq),
            )
            list_of_amplicon_primer_temps[i].append(melting_temps)

    return list_of_amplicons, list_of_amplicon_primers, list_of_amplicon_primer_temps


def get_primers(
    List_of_assemblies: list,
    combinatorial_list_of_names: list,
    combinatorial_list_of_primer_tm: list,
):
    """Returns a list of ALL primers from the combinatorial library,
    updates names and what they anneal to.

    Parameters
    ----------
    List_of_assemblies : list[list[pydna.amplicon.Amplicon]]
    combinatorial_list_of_names : list[(str)]
    combinatorial_list_of_primer_tm : list[(float, float),..)...]

    Returns
    -------
    primers : list[list[[pydna.primer.Primer, pydna.primer.Primer]]
        All primers that have been made for all assemblies
    """

    primers_temporary = []
    primers = []

    counter = 0
    for i in range(0, len(List_of_assemblies)):
        for j in range(0, len(List_of_assemblies[i])):
            counter += 1
            # Names
            List_of_assemblies[i][j].name = combinatorial_list_of_names[i][j]
            # Primers
            # description ------ DESCRIBES what other part it overlaps-------------
            if j == 0:  # START OF THE ASSEMBLY
                List_of_assemblies[i][
                    j
                ].forward_primer.description = "Anneals to " + str(
                    List_of_assemblies[i][j].name
                )
                List_of_assemblies[i][j].reverse_primer.description = (
                    "Anneals to "
                    + str(List_of_assemblies[i][j].name)
                    + ", overlaps to "
                    + str(List_of_assemblies[i][j + 1].name)
                )
            if j > 0 and j < len(List_of_assemblies[i]) - 1:  #      # THE rest:
                List_of_assemblies[i][
                    j
                ].forward_primer.description = "Anneals to " + str(
                    List_of_assemblies[i][j].name
                    + ", overlaps to "
                    + str(List_of_assemblies[i][j - 1].name)
                )
                List_of_assemblies[i][
                    j
                ].reverse_primer.description = "Anneals to " + str(
                    List_of_assemblies[i][j].name
                    + ", overlaps to "
                    + str(List_of_assemblies[i][j + 1].name)
                )
            if j == len(List_of_assemblies[i]) - 1:  # THE END OF THE ASSEMBLY
                List_of_assemblies[i][j].forward_primer.description = (
                    "Anneals to "
                    + str(List_of_assemblies[i][j].name)
                    + ", overlaps to "
                    + str(List_of_assemblies[i][j - 1].name)
                )
                List_of_assemblies[i][
                    j
                ].reverse_primer.description = "Anneals to " + str(
                    List_of_assemblies[i][j].name
                )

            # template it aneals to
            List_of_assemblies[i][j].forward_primer.name = str(
                List_of_assemblies[i][j].name
            )
            List_of_assemblies[i][j].reverse_primer.name = str(
                List_of_assemblies[i][j].name
            )

            # Primer tm
            List_of_assemblies[i][j].forward_primer.features = round(
                float(combinatorial_list_of_primer_tm[i][j][0]), 2
            )
            List_of_assemblies[i][j].reverse_primer.features = round(
                float(combinatorial_list_of_primer_tm[i][j][1]), 2
            )

            fwd_rev_primers = [
                List_of_assemblies[i][j].forward_primer,
                List_of_assemblies[i][j].reverse_primer,
            ]
            primers_temporary.append(fwd_rev_primers)

        primers.append(primers_temporary)
        primers_temporary = []

    return primers


def assembly_maker(combinatorial_list_of_amplicons: list, overlap=35):
    """Assembles Amplicons with pad and makes new overlapping primers.

    Parameters
    ----------
    combinatorial_list_of_amplicons : list[[pydna.amplicon.Amplicon]]
        the list of pydna.amplicon.Amplicon that you want generate
        overlapping primers for.
    overlap : int (default set to 35)
        How many basepair overlaps

    Returns
    -------
    List_of_assemblies : list[[pydna.amplicon.Amplicon]]
        amplicons that overlaps eachother with the specified overlap value.

    """

    List_of_assemblies = []
    for i in range(0, len(combinatorial_list_of_amplicons)):
        List_of_assemblies.append(
            assembly_fragments(combinatorial_list_of_amplicons[i], overlap, maxlink=40)
        )

    return List_of_assemblies


def unique_primers(primers: list, list_of_assemblies):
    """Finds unique primers from a list of assemblies
    Parameters
    ----------
    primers : list[list[[pydna.primer.Primer, pydna.primer.Primer]]
        a list of all the primers made for the combinatorial library

    list_of_assemblies: list[[pydna.amplicon.Amplicon]]
        used here to update the names of the primers

    Returns
    -------
    unique_primers : list[list(ID,Anneals_to,Sequence,Annealing_temp,Length,Price(DKK))]
        Relevant metrics for the unique primers of the combinatorial library.

    """

    unikke_F_primers = []
    unikke_R_primers = []
    length_of_unique_primers = 0
    counter = 0
    primer_list = []

    for i in range(0, len(primers)):
        for j in range(0, len(primers[i])):
            counter += len(primers[i][j])
            if primers[i][j][0] not in unikke_F_primers:
                unikke_F_primers.append(primers[i][j][0])
            if primers[i][j][1] not in unikke_R_primers:
                unikke_R_primers.append(primers[i][j][1])

    counter = 0
    unique_forward_primers = []
    unique_reverse_primers = []

    ### CHANGING THE NAMES OF THE PRIMERS
    # Forward primers
    for i in range(len(unikke_F_primers)):
        counter += 1
        unikke_F_primers[i].id = "F{number:03}".format(number=counter)
        length_of_unique_primers += len(unikke_F_primers[i].seq)
        U_f_primers = [
            unikke_F_primers[i].id,
            unikke_F_primers[i].name,
            unikke_F_primers[i].seq,
            unikke_F_primers[i].features,  # anealing temp
            len(unikke_F_primers[i].seq),  # lenght
            len(unikke_F_primers[i].seq) * 1.8,  # price
        ]
        unique_forward_primers.append(U_f_primers)
    # Reverse primers
    for i in range(len(unikke_R_primers)):
        counter += 1
        unikke_R_primers[i].id = "R{number:03}".format(number=counter)
        length_of_unique_primers += len(unikke_R_primers[i].seq)
        U_r_primers = [
            unikke_R_primers[i].id,
            unikke_R_primers[i].name,
            unikke_R_primers[i].seq,
            unikke_R_primers[i].features,
            len(unikke_R_primers[i].seq),
            len(unikke_R_primers[i].seq) * 1.8,  # cost
        ]
        unique_reverse_primers.append(U_r_primers)

    primer_list = (
        unique_forward_primers + unique_reverse_primers
    )  # COULD CONCATONATE THEM INTO: unique_forward_primers + unique_reverse_primers

    ### Updating primer names and removing duplicates
    for i in range(0, len(list_of_assemblies)):
        for j in range(0, len(list_of_assemblies[i])):
            for l in range(0, len(unikke_F_primers)):
                if (
                    list_of_assemblies[i][j].forward_primer.seq
                    == unikke_F_primers[l].seq
                ):
                    list_of_assemblies[i][j].forward_primer = unikke_F_primers[l]
            for m in range(0, len(unique_reverse_primers)):
                if (
                    list_of_assemblies[i][j].reverse_primer.seq
                    == unikke_R_primers[m].seq
                ):
                    list_of_assemblies[i][j].reverse_primer = unikke_R_primers[m]

    return primer_list


def unique_amplicons(list_of_assemblies: list):

    """Finds Unique amplicons from a list of assemblies
    Parameters
    ----------
    list_of_assemblies: list[[pydna.amplicon.Amplicon]]
        list of the combinatorial libarary with overlapping ends

    Returns
    -------
        unique_amplicons: list[pydna.amplicon.Amplicon]
            returns a list of unique amplicons where relavant metrics
            are added to the objects.
    """
    ### Unique amplicons
    unique_amplicons = []
    for i in range(0, len(list_of_assemblies)):
        for j in range(0, len(list_of_assemblies[i])):
            if list_of_assemblies[i][j] not in unique_amplicons:
                unique_amplicons.append(list_of_assemblies[i][j])

    return unique_amplicons


def making_assembly_objects(list_of_assemblies: list):
    """Assembling amplicons into assembling class that shows
    fragments, limit,nodes and which algorithm that was used
    for assembling.

    Parameters
    ----------
    list_of_assemblies: list[[pydna.amplicon.Amplicon]]
        list of the combinatorial libarary with overlapping ends

    Returns
    -------
        list_of_assembly_objects: list[pydna.assembly.Assembly]
            shows which algorithm that was used, nodes, limit and fragments

    """
    list_of_assembly_objects = []
    for i in range(0, len(list_of_assemblies)):
        list_of_assembly_objects.append(Assembly((list_of_assemblies[i]), limit=35))

    return list_of_assembly_objects


def making_assembled_contigs(list_of_assembly_objects: list):
    """Assembles a list of assembly object into
    linear contigs.

    Parameters
    ----------
    list_of_assembly_objects : list[pydna.assembly.Assembly]
        these objects can be assembled into contigs

    Returns
    -------
    list_of_assembly_objects : list[]
        list_of_assembly_objects have been assembled into contigs
    """
    contigs_assembled = []
    for j in range(0, len(list_of_assembly_objects)):
        contigs_assembled.append(list_of_assembly_objects[j].assemble_linear())

    return list_of_assembly_objects


class DesignAssembly:
    """Class able to make a combinatorial library from DNA fragments.

    Parameters
    ----------
    list_of_seqs : list
        A list of list of a constructs of choice.
    list_of_names : list
        A list of list of the names wanted for the construct of choice.
    pad : pydna.Dseqrecord
        A nucleotide sequence to be incorporated into the primers (Max is 40 bp)
    position_of_pad : int
        the position in the list of seqs where the pad is incorporated (zero indexed)

    Returns
    -------
    teemi.design.combinatorial_design.DesignAssembly object
        A powerfull class and a lot of information can be retrieved.
        Such as: showing all the amplicons needed to construct a combinatorial library
        with the simple method --> PCR_list_to_dataframe or Primer_list_to_dataframe.

    """

    def __init__(
        self,
        list_of_seqs: list,
        list_of_names: list,
        pad: str,
        position_of_pad: int,
        target_tm=56.0,
        limit=13,
        overlap=35,
    ):

        ###  1.INITIALIZING ##
        self.list_of_seqs = list_of_seqs
        self.list_of_names = list_of_names
        self.pad = pad
        self.position_of_pad = position_of_pad

        ### 2. Amplicons, primers, and their temperatures
        (
            self.list_of_amplicons,
            self.list_of_amplicon_primers,
            self.list_of_amplicon_primer_temps,
        ) = simple_amplicon_maker(
            self.list_of_seqs, self.list_of_names, target_tm=target_tm, limit=limit
        )

        # Systematic names
        self.systematic_names = systematic_names_function(self.list_of_seqs)

        ### 3. COMBINATORIAL LISTS
        self.combinatorial_list_of_amplicons = combinatorial_list_maker(
            self.list_of_amplicons
        )
        self.combinatorial_list_of_names = combinatorial_list_maker(self.list_of_names)
        self.combinatorial_list_of_primer_tm = combinatorial_list_maker(
            self.list_of_amplicon_primer_temps
        )

        # Making the combinations into a list so we can insert PADS later (They are tuples at this stage, and insert doesnt work for tuples)
        for i in range(0, len(self.combinatorial_list_of_amplicons)):
            self.combinatorial_list_of_amplicons[i] = list(
                self.combinatorial_list_of_amplicons[i]
            )

        #### 4. Adding PAD ###
        for i in range(0, len(self.combinatorial_list_of_amplicons)):
            self.combinatorial_list_of_amplicons[i].insert(
                self.position_of_pad, self.pad
            )

        ### 5. Assembling and making overlapping primers
        self.list_of_assemblies = assembly_maker(
            self.combinatorial_list_of_amplicons, overlap=overlap
        )

        ### 6. GETTING all primers, annotating, adding features
        self.primers = get_primers(
            self.list_of_assemblies,
            self.combinatorial_list_of_names,
            self.combinatorial_list_of_primer_tm,
        )

        ### 7. Getting Unique primers and re-annotating list_assemblies to get right names
        self.unique_primers = unique_primers(self.primers, self.list_of_assemblies)

        ### 8. Unique amplicons
        self.unique_amplicons = unique_amplicons(self.list_of_assemblies)

    def ShowContigs(self):
        """Returns a string of the contigs generated by the assembly"""
        print("Template, Primer, tm")
        for i in range(0, len(self.list_of_assemblies)):
            print("\nContig" + str(self.systematic_names[i]))
            for j in range(0, len(self.list_of_assemblies[i])):
                print(
                    "Template: ", self.list_of_assemblies[i][j].name[0:15]
                )  # , '\t', self.primers[i][j][0].name,'\t',self.primers[i][j][0].features)
        return

    def ShowVariantsLibDF(self):
        """Returns a dataframe of all the variants"""
        combinatorial_lib_variants_df = pd.DataFrame(self.combinatorial_list_of_names)
        systematic_names = self.systematic_names
        combinatorial_lib_variants_df["Systematic_name"] = systematic_names
        combinatorial_lib_variants_df["Variant"] = np.arange(
            len(combinatorial_lib_variants_df)
        )

        return combinatorial_lib_variants_df

    def print_primer_list(self):
        """Return the list of transfers in human-readable format."""
        for primers in self.unique_primers:
            print(primers)

    def primer_list(self):
        """Return the list of transfers in human-readable format."""
        primer_list = []
        for primers in self.unique_primers:
            primer_list.append(primers)

        return primer_list

    def primer_list_to_dataframe(self):
        """Return a pandas dataframe with list of primers."""
        df = pd.DataFrame(self.unique_primers)
        df.columns = [
            "ID",
            "Anneals to",
            "Sequence",
            "Annealing temperature",
            "Length",
            "Price(DKK)",
        ]

        return df

    def print_PCR_list(self):
        """Prints PCR_list"""
        print("PCR#, Template,forward_primer, reverse primer, F_tm, R_tm")
        for i in range(0, len(self.unique_amplicons)):
            print(
                "PCR{number}".format(number=i + 1),
                ",",
                self.unique_amplicons[i].name,
                ",",
                self.unique_amplicons[i].forward_primer.id,
                ",",
                self.unique_amplicons[i].reverse_primer.id,
                ",",
                self.unique_amplicons[i].forward_primer.features,
                ",",
                self.unique_amplicons[i].reverse_primer.features,
            )

    def PCR_list(self):
        """Returns a PCR_list"""
        pcr_list = []
        for i in range(0, len(self.unique_amplicons)):
            PCR = [
                "PCR{number}".format(number=i + 1),
                self.unique_amplicons[i].name,
                self.unique_amplicons[i].forward_primer.id,
                self.unique_amplicons[i].reverse_primer.id,
                self.unique_amplicons[i].forward_primer.features,
                self.unique_amplicons[i].reverse_primer.features,
            ]
            pcr_list.append(PCR)

        return pcr_list

    def PCR_list_to_dataframe(self):
        """Prints PCR_list into a pandas dataframe"""
        dataframe_list = []
        for i in range(0, len(self.unique_amplicons)):
            lst = [
                "PCR{number}".format(number=i + 1),
                self.unique_amplicons[i].name,
                self.unique_amplicons[i].forward_primer.id,
                self.unique_amplicons[i].reverse_primer.id,
                self.unique_amplicons[i].forward_primer.features,
                self.unique_amplicons[i].reverse_primer.features,
            ]
            dataframe_list.append(lst)

        df = pd.DataFrame(dataframe_list)
        df.columns = [
            "PCR#",
            "Template",
            "forward_primer",
            "reverse_primer",
            "F_tm",
            "R_tm",
        ]

        return df

    def graphical_representation_of_assemblies(self):
        """
        Takes in the assembly object and returns graphical report of the
        fragments assembled
        """
        graphical_representation = [
            self.assembly_object[x].assemble_linear()[0].figure()
            for x in range(0, len(self.assembly_object))
        ]

        return graphical_representation


def count_unique_parts(df, max_combinations: int):
    """Iterate through the list of predictions and save new encountered parts.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing predictions

    Returns:
    --------
    parts_encounteres : dict
        A dictionary containing the unique parts encountered in 'G8H','pG8H', 'pCPR', 'CPR' columns,
        total number of unique combinations encountered in 'Sum of parts' and total predictions
        encountered in 'Predictions'

    """
    # Iterate through the list of predictions and save new encountered parts. Stop after 180 combiantions.
    # Initialisation
    parts_encounteres = {
        "G8H": [],
        "pG8H": [],
        "pCPR": [],
        "CPR": [],
        "Sum of parts": "",
        "Predictions": "",
    }
    sum_of_parts = 0
    i = 0
    g8h_count = 0
    cpr_count = 0
    pg8h_count = 0
    pcpr_count = 0

    # Loop through the predctions and save new parts.
    while sum_of_parts < max_combinations:
        sum_of_parts = g8h_count * cpr_count * pg8h_count * pcpr_count

        parts_encounteres["Sum of parts"] = str(sum_of_parts)
        parts_encounteres["Predictions"] = str(i)

        g8h = df.G8H[i]
        pg8h = df.pG8H[i]
        cpr = df.CPR[i]
        pcpr = df.pCPR[i]
        if g8h not in parts_encounteres["G8H"]:
            parts_encounteres["G8H"].append(g8h)
            g8h_count += 1
        if pg8h not in parts_encounteres["pG8H"]:
            parts_encounteres["pG8H"].append(pg8h)
            pg8h_count += 1
        if cpr not in parts_encounteres["CPR"]:
            parts_encounteres["CPR"].append(cpr)
            cpr_count += 1
        if pcpr not in parts_encounteres["pCPR"]:
            parts_encounteres["pCPR"].append(pcpr)
            pcpr_count += 1
        i += 1

    return parts_encounteres
