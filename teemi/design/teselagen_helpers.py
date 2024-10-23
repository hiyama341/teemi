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

import numpy as np
import pydna
from teemi.design.cloning import seq_to_annotation


def primer_matrix_teselagen(
    teselagen_primer_matrix,
    teselagen_parts_matrix,
    combinations_matrix,
    no_combs: int,
    no_frags,
):
    """Makes an primer matrix from teselagen output.

    Parameters
    ----------
    teselagen_primer_matrix : pandas DataFrame
        DataFrame with teselagen primer information, including the columns 'Name', 'Sequence', and 'Location'.
    teselagen_parts_matrix : pandas DataFrame
        DataFrame with teselagen part information, including the columns 'Sequence' and 'Overlaps_to'.
    combinations_matrix : list of list of Dseqrecord
        List of lists of Dseqrecord objects representing the combinations of parts.
    no_combs : int
        Number of combinations.

    Returns
    -------
    amplicon_matrix : list of list of tuple
        List of lists of tuples containing primer pairs for each combination.
    """
    # Initialize matrix
    primer_matrix = [[] for i in range(no_combs)]
    for comb_no in range(0, no_combs):
        for frag_no in range(0, no_frags):
            comb_seq = str(combinations_matrix[comb_no][frag_no].seq.upper())

            # index used for getting the names from overlapping parts
            if frag_no == 1:
                idx = frag_no + 1
            elif frag_no == 4:
                idx = frag_no - 1
            else:
                idx = frag_no
            #
            if idx != frag_no:
                name = combinations_matrix[comb_no][idx].name
                true_false = np.logical_and(
                    teselagen_parts_matrix["Sequence"].str.contains(comb_seq),
                    teselagen_parts_matrix["Overlaps_to"] == name,
                )
            else:
                true_false = teselagen_parts_matrix["Sequence"].str.contains(comb_seq)

            comb_seq_idx = teselagen_parts_matrix[true_false].index

            forward_primer_name = teselagen_parts_matrix.loc[
                comb_seq_idx, ["Forward Oligo Name"]
            ].values[0][0]
            reverse_primer_name = teselagen_parts_matrix.loc[
                comb_seq_idx, ["Reverse Oligo Name"]
            ].values[0][0]

            forward_primer_idx = teselagen_primer_matrix[
                teselagen_primer_matrix["Name"] == forward_primer_name
            ].index
            forward_seq = teselagen_primer_matrix.loc[
                forward_primer_idx, ["Sequence"]
            ].values[0][0]
            forward_loc = teselagen_primer_matrix.loc[
                forward_primer_idx, ["Location"]
            ].values[0][0]
            forward_primer = pydna.primer.Primer(
                pydna.dseqrecord.Dseqrecord(forward_seq, id=forward_primer_name)
            )

            forward_primer.annotations["batches"] = []
            forward_primer.annotations["batches"].append(
                {"location": forward_loc, "volume": 100, "concentration": 10}
            )
            reverse_primer_idx = teselagen_primer_matrix[
                teselagen_primer_matrix["Name"] == reverse_primer_name
            ].index
            reverse_seq = teselagen_primer_matrix.loc[
                reverse_primer_idx, ["Sequence"]
            ].values[0][0]
            reverse_loc = teselagen_primer_matrix.loc[
                reverse_primer_idx, ["Location"]
            ].values[0][0]
            reverse_primer = pydna.primer.Primer(
                pydna.dseqrecord.Dseqrecord(reverse_seq, id=reverse_primer_name)
            )
            reverse_primer.annotations["batches"] = []
            reverse_primer.annotations["batches"].append(
                {"location": reverse_loc, "volume": 100, "concentration": 10}
            )

            primer_matrix[comb_no].append((forward_primer, reverse_primer))

    return primer_matrix


def amplicon_matrix_teselagen(
    teselagen_parts_matrix,
    primer_matrix,
    combinations_matrix,
    no_combs: int,
    no_frags: int,
):
    """Makes an amplicon_matrix from tesselagen output.

    Parameters
    ----------

    teselagen_parts_matrix : pd.DataFrame
    primer_matrix : list of list
    combinations_matrix : list of list
    no_combs : int

    Returns
    -------
    amplicon_matrix : list of list"""

    # Initizalise
    amplicon_matrix = [[] for i in range(no_combs)]

    for comb_no in range(0, no_combs):
        for frag_no in range(0, no_frags):
            template = pydna.dseqrecord.Dseqrecord(
                combinations_matrix[comb_no][frag_no]
            )

            forward_primer = primer_matrix[comb_no][frag_no][0]
            reverse_primer = primer_matrix[comb_no][frag_no][1]

            amplicon = pydna.amplify.pcr(forward_primer, reverse_primer, template)

            amplicon.annotations["template_name"] = template.name

            comb_seq = amplicon.seq.watson.upper()
            true_false = teselagen_parts_matrix["Sequence"].str.contains(comb_seq)
            comb_seq_idx = teselagen_parts_matrix[true_false].index

            amplicon_name = teselagen_parts_matrix.loc[comb_seq_idx, ["Name"]].values[
                0
            ][0]
            amplicon.name = amplicon_name

            amplicon_loc = teselagen_parts_matrix.loc[
                comb_seq_idx, ["Location"]
            ].values[0][0]
            amplicon_vol = teselagen_parts_matrix.loc[comb_seq_idx, ["volume"]].values[
                0
            ][0]
            amplicon_con = teselagen_parts_matrix.loc[
                comb_seq_idx, ["concentration"]
            ].values[0][0]
            amplicon.annotations["batches"] = []
            amplicon.annotations["batches"].append(
                {
                    "location": amplicon_loc,
                    "volume": amplicon_vol,
                    "concentration": amplicon_con,
                }
            )

            amplicon.forward_primer.annotations["tm Q5 Hot Start"] = (
                teselagen_parts_matrix.loc[comb_seq_idx, ["Q5_fw_tm"]].values[0][0]
            )
            amplicon.reverse_primer.annotations["tm Q5 Hot Start"] = (
                teselagen_parts_matrix.loc[comb_seq_idx, ["Q5_rv_tm"]].values[0][0]
            )
            amplicon.annotations["ta Q5 Hot Start"] = teselagen_parts_matrix.loc[
                comb_seq_idx, ["Q5_ta"]
            ].values[0][0]

            ## Remove additional wrong primer annotation added by pydna.amplify.pcr
            amplicon.features = [
                feat for feat in amplicon.features if feat.type != "primer_bind"
            ]

            seq_to_annotation(amplicon, amplicon, "PCR_product")
            seq_to_annotation(forward_primer, amplicon, "primer_bind")
            seq_to_annotation(reverse_primer, amplicon, "primer_bind")

            amplicon_matrix[comb_no].append(amplicon)

    return amplicon_matrix
