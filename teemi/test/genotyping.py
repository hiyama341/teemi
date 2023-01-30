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

import pandas as pd
import numpy as np
from Bio import pairwise2


def slicing_and_naming_seq_plates(sequencing_plates, where_to_slice=7) -> list:
    """Slices rows of a list of dataframes and changes the names.
    Is used to ease pre-processing of Plate2seq excel files

    Parameters
    ----------
    sequencing_plates : list of pd.DataFrames
        Plate2seq pd.dataframes
    where_to_slice : int
        indicate where to slice the dataframe

    Returns
    -------
    list of plates sliced pd.DataFrames

    """
    # changing column names and slicing
    for i in range(len(sequencing_plates)):
        sequencing_plates[i].columns = (
            "Number",
            "Sample-Name",
            "AvgQual",
            "Length",
            "GoodQualFrom",
            "GoodQualTo",
            "used",
        )  # sequencing_plates[0].iloc[6]
        sequencing_plates[i] = sequencing_plates[i][where_to_slice:]

    return sequencing_plates


def plat_seq_data_wrangler(sequencing_plates: list) -> list:
    """Makes list of Plate2Seq pd.DataFrames into numeric
    values and removes nan values.

    Parameters
    ----------
    sequencing_plates : list of pd.DataFrames
        Sliced Plate2seq pd.dataframes

    Returns
    -------
    Plate2Seq pd.DataFrames with numeric values

    """

    list_with_dfs = []

    for i in range(len(sequencing_plates)):
        # taking only a subset of the dataframe:
        numeric_values = sequencing_plates[i][
            ["AvgQual", "Length", "GoodQualFrom", "GoodQualTo", "used"]
        ]

        # if values are non nummeric make them NaN
        numeric_values = numeric_values.replace("n.a.", np.NaN)

        # Making them numeric
        numeric_values = numeric_values.apply(pd.to_numeric, errors="coerce")

        # Adding names column
        name_column = sequencing_plates[i]["Sample-Name"]
        number_column = sequencing_plates[i]["Number"]

        # Adding them to the dataframe
        data1 = pd.concat([number_column, name_column, numeric_values], axis=1)
        list_with_dfs.append(data1)

    return list_with_dfs


def plate_AvgQual(list_of_dfs_numeric: list, Avg_qual=50, used_bases=25) -> list:
    """Filters out rows that doesnt follow the criteria.

    Parameters
    ----------
    list_of_dfs_numeric : list of pd.DataFrames
        Sliced and Plate2seq pd.dataframes
    Avg_qual : int
    used_bases : int

    Returns
    -------
    Plate2Seq pd.DataFrames with that follows Avg_qual and used_bases criteria

    """
    # Initialize
    filtered_plates = []

    for i in range(len(list_of_dfs_numeric)):
        # Filter
        filter_Avg_qual = list_of_dfs_numeric[i][
            list_of_dfs_numeric[i]["AvgQual"] > Avg_qual
        ]
        filer_used_bases = filter_Avg_qual[filter_Avg_qual["used"] > used_bases]
        # Save the filtered plates
        filtered_plates.append(filer_used_bases)

    return filtered_plates


def split_df_names(
    df_names_column, which_column_to_split1=0, which_column_to_split2=2
) -> list:
    """Splits sample names from plate2seq pd.dataframes into plate and well columns"""
    df_with_names_split = []

    for i in range(len(df_names_column)):
        # splitting
        df_filter_plates = df_names_column[i]["Sample-Name"].str.split("_", expand=True)

        # selecting
        column1 = df_filter_plates[which_column_to_split1]
        column2 = df_filter_plates[which_column_to_split2]

        # concating
        concatenated = pd.concat(
            [df_names_column[i], column1, column2], axis=1, ignore_index=False
        )

        # changing names
        concatenated.columns = (
            "Number",
            "Sample-Name",
            "AvgQual",
            "Length",
            "GoodQualFrom",
            "GoodQualTo",
            "used",
            "plate",
            "well",
        )

        # save
        df_with_names_split.append(concatenated)

    return df_with_names_split


def concatenating_list_of_dfs(list_of_dfs: list):
    """Concatenating a list of daframes into one pd.dataframe by rows"""
    assembled_dfs = pd.concat(list_of_dfs, axis=0, ignore_index=False)

    return assembled_dfs


def pairwise_alignment_of_templates(
    reads: list, templates: list, primers: list
) -> dict:
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

        # If we see the primers in the sample we the alignment will start from there
        for k in range(len(primers)):
            start = sample.find(primers[k].seq)

            if start != -1:
                sample[start:]
            else:
                continue

        # Aling with templates
        if len(sample) > 25:
            score = 0.0
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
                    read_name = reads[i].name

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
