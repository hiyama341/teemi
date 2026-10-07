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

import pandas as pd
from Bio import pairwise2


def pairwise_alignment_of_templates(
    reads: list, templates: list, primers: list, trim_at_primer: bool = False
) -> pd.DataFrame:
    """Infers relationship of templates to reads based on highest
    score from a pairwise alignment.

    Parameters
    ----------
    reads: list of Bio.SeqRecord.SeqRecord
        these are .ab1 files made into Bio.SeqRecord.SeqRecord objects
    templates: list of Bio.SeqRecord.SeqRecord
        Templates for inferring relationship with - could be plasmid fx
    primers: list of Bio.SeqRecord.SeqRecord
        list of primers to be for finding were the read should start
    trim_at_primer: bool
        if True, reads are trimmed to start at the primer before aligning.
        Defaults to False: whole reads are aligned, as in the published notebooks.

    Returns
    -------
    pd.Dataframe in the following way:

    Example
    -------
    <<<df_alignment = pairwise_alignment_of_templates(reads,templates, primers_for_seq)

    <<< df_alignment

    Sample-Name	inf_promoter_name	align_score	inf_promoter
    132	yp53re_cpr_A10_A10-pad_cpr_fw	pCCW12	634.0	5
    188	yp53re_cpr_A11_A11-pad_cpr_fw	pTPI1	904.0	6
    247	yp53re_cpr_A12_A12-pad_cpr_fw	pTPI1	851.0	6
    93	yp53re_cpr_A1_A01-pad_cpr_fw	pCCW12	543.0	5
    41	yp53re_cpr_A2_A02-pad_cpr_fw	pCCW12	636.0	5

    Notes
    -----
    Reads of 25 bp or less (after removing Ns) are too short to call; they get
    an align_score of 0 and no inferred part.

    If you want inf_part_number column then change your the description
    of the Bio.SeqRecord.SeqRecord as follows:

    pCCW12.description = '1'
    """

    best_scores = []
    read_list = []
    template_list = []
    template_number_list = []

    for i in range(len(reads)):

        sample = reads[i].seq.replace("N", "")

        # If we see the primers in the sample the alignment can start from there
        if trim_at_primer:
            for primer in primers:
                start = sample.find(primer.seq)
                if start != -1:
                    sample = sample[start:]

        score = 0.0
        temp_name = None
        temp_number = None
        read_name = reads[i].name

        # Aling with templates
        if len(sample) > 25:
            for j in range(len(templates)):
                template = templates[j].seq

                # Align # identical = 1, non-identical = -2 , gap = -2 , extending gap = -2
                # alignments = pairwise2.align.globalms(template, sample,2, -2, -3, -3)
                alignment_score = pairwise2.align.localxx(
                    template, sample, score_only=True
                )

                # Get the best alignment of them all
                if alignment_score > score:
                    score = alignment_score
                    temp_name = templates[j].name
                    temp_number = templates[j].description

        # Saving the alignmets and their names
        best_scores.append(score)
        read_list.append(read_name)
        template_list.append(temp_name)
        template_number_list.append(temp_number)

    # Making a pandas. dataframe
    df = pd.DataFrame()
    df["Sample-Name"] = read_list
    df["inf_part_name"] = template_list
    df["align_score"] = best_scores
    df["inf_part_number"] = template_number_list

    return df
