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

""" This part of the design module is used fetching gene homologs"""


import re
from Bio import pairwise2
import pandas as pd
from dnachisel import DnaOptimizationProblem, EnforceGCContent, CodonOptimize
from typing import List
from Bio.SeqRecord import SeqRecord


def alignment_identity(query: list, reference: str) -> list:
    """Calculates percent identity between a reference and query(s).
    Parameters
    ----------
    query : list
        list of Biopython Seqrecord objects
    reference : str

    Returns
    -------
    list of percent identeties as floats
    """

    alignment_score = []
    for alignment in query:
        alignment_score.append(
            pairwise2.align.globalmx(
                reference.seq, alignment.seq, 1, 0, score_only=True
            )
            / len(reference.seq)
        )

    return alignment_score


def filter_blast_results(
    blast_record,
    E_VALUE_THRESH=0.4,
    LOWER_PROTEIN_IDENTITY_THRESH=0.1,
    UPPER__PROTEIN_IDENTITY_THRESH=1,
    show_alignment=False,
):
    Alignments_that_follow_our_criteria = []  # These are the ACC numbers

    # saving some of the metrics
    Names = []
    Identity = []
    E_value = []
    Length = []
    ACC_number = []

    counter = 0
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            LENGTH = alignment.length  #

            IDENTITY = hsp.identities / hsp.align_length  # not super accurate
            # We add the added reqirement that the sequence should be within the AA identeties
            if (
                hsp.expect < E_VALUE_THRESH
                and IDENTITY > LOWER_PROTEIN_IDENTITY_THRESH
                and IDENTITY < UPPER__PROTEIN_IDENTITY_THRESH
            ):
                counter += 1
                if show_alignment:
                    print("\nAlignment#", counter)
                    print("Name:", alignment.hit_def)
                    print("Title:", alignment.title)
                    print("Length:", alignment.length)
                    print("E value:", hsp.expect)
                    print("Query:", hsp.query[0:75] + "...")
                    print("Match:", hsp.match[0:75] + "...")
                    print("Subjt:", hsp.sbjct[0:75] + "...")
                    print("Identitiy", "{:.2f}".format(IDENTITY))
                    print(alignment)

                # Saving the metrics we want into different lists
                Alignments_that_follow_our_criteria.append(alignment)
                Names.append(alignment.hit_def)
                Identity.append(IDENTITY)
                E_value.append(hsp.expect)
                Length.append(alignment.length)
                ACC_number.append(alignment.accession)

    if show_alignment:
        print("\nTOTAL HOMOLOGS", counter)

    data = list(zip(Names, Identity, E_value, Length, ACC_number))
    df_alignments_that_follow_our_criteria = pd.DataFrame(
        data, columns=["Name", "Identity", "E_value", "Length", "ACC_number"]
    )

    return df_alignments_that_follow_our_criteria


def codon_optimize_with_dnachisel(
    sequences: List[SeqRecord],
    lower_GC: float = 0.3,
    upper_GC: float = 0.7,
    species: str = None,
    codon_usage_table=None,
    window: int = 50,
) -> List[SeqRecord]:
    """Codon-optimize sequences with_dnachisel.

    Parameters
    ----------
    sequences : list
        list of Bio.SeqRecord objects
    lower_GC : float
        the lowest GC content in the region of 50 bp
    upper_GC : float
        the highest GC content in the region of 50 bp
    species : str
        name of the species for which to optimize the sequence.
        examples: 'e_coli, s_cerevisiae, h_sapiens, c_elegans, b_subtilis, d_melanogaster
        check python_codon_tables for more info.
    codon_usage_table:
        a codon table following the structure of:
        {'*': {'TAA': 0.0, 'TAG': 0.0, 'TGA': 1.0},...


    Returns
    -------
    list of codon optimized sequences for yeast
    """

    if not species and not codon_usage_table:
        raise ValueError(
            "At least one of `species` and `codon_usage_table` must be specified."
        )

    codon_optimized_seqs = []

    # DEFINE THE OPTIMIZATION PROBLEM
    for seq in sequences:
        if species:
            problem = DnaOptimizationProblem(
                sequence=seq.seq,
                constraints=[
                    EnforceGCContent(mini=lower_GC, maxi=upper_GC, window=window)
                ],
                objectives=[CodonOptimize(species=species)],
            )

            # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE
            problem.resolve_constraints()
            problem.optimize()

            print(problem.constraints_text_summary())
            print(problem.objectives_text_summary())

        elif codon_usage_table:
            problem = DnaOptimizationProblem(
                sequence=seq.seq,
                constraints=[
                    EnforceGCContent(mini=lower_GC, maxi=upper_GC, window=window)
                ],
                objectives=[CodonOptimize(codon_usage_table=codon_usage_table)],
            )

        # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE
        problem.resolve_constraints()
        problem.optimize()

        print(problem.constraints_text_summary())
        print(problem.objectives_text_summary())

        # GET THE FINAL SEQUENCE AS ANNOTATED BIOPYTHON RECORDS)
        final_record = problem.to_record(with_sequence_edits=True)
        final_record.id = seq.id
        final_record.name = seq.name
        final_record.description = seq.description

        codon_optimized_seqs.append(final_record)

    return codon_optimized_seqs


def find_all_starts(seq):
    """Find the starting index of all start codons in a lowercase seq
    This function was made by Justin Bois : http://justinbois.github.io/.

    """
    # Initialize array of indices of start codons
    starts = []

    # Find index of first start codon (remember, find() returns -1 if not found)
    i = seq.find("atg")

    # Keep looking for subsequence incrementing starting point of search
    while i >= 0:
        starts.append(i)
        i = seq.find("atg", i + 1)

    return tuple(starts)


def find_first_in_register_stop(seq):
    """
    Find first stop codon on lowercase seq that starts at an index
    that is divisible by three.
    This function was made by Justin Bois : http://justinbois.github.io/.

    """
    # Compile regexes for stop codons
    regex_stop = re.compile("(taa|tag|tga)")

    # Stop codon iterator
    stop_iterator = regex_stop.finditer(seq)

    # Find next stop codon that is in register
    for stop in stop_iterator:
        if stop.end() % 3 == 0:
            return stop.end()

    # Return -1 if we failed to find a stop codon
    return -1


def all_orfs(seq):
    """Return all ORFs of a sequence.
    This function was made by Justin Bois : http://justinbois.github.io/."""
    # Make sure sequence is all lower case
    seq = seq.lower()

    # Find the indices of all start codons
    start_inds = find_all_starts(seq)

    # Keep track of stops
    stop_inds = []

    # Initialze ORFs.  Each entry in list is [ORF length, ORF start, ORF stop]
    orfs = []

    # For each start codon, find the next stop codon in register
    for start in start_inds:
        relative_stop = find_first_in_register_stop(seq[start:])

        if relative_stop != -1:
            # Index of stop codon
            stop = start + relative_stop

            # If already had stop, a longer ORF contains this one
            if stop not in stop_inds:
                orfs.append((relative_stop, start, stop))
                stop_inds.append(stop)

    # Get sorted list of ORF length
    orfs = sorted(orfs, reverse=True)

    # Remove lengths
    for i, orf in enumerate(orfs):
        orfs[i] = (orf[1], orf[2])

    return tuple(orfs)


def longest_orf(seq, n=1):
    """Longest ORF of a sequence.
    This function was made by Justin Bois : http://justinbois.github.io/."""
    orfs = all_orfs(seq)

    if len(orfs) == 0:
        return ""
    elif n == 1 or len(orfs) == 1:
        return seq[orfs[0][0] : orfs[0][1]]
    else:
        return_list = []
        for i in range(min(n, len(orfs))):
            return_list.append(seq[orfs[i][0] : orfs[i][1]])

        return tuple(return_list)
