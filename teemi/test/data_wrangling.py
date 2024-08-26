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
import numpy as np


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
        )
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

        # if values are non numeric make them NaN
        numeric_values = numeric_values.replace("n.a.", np.nan)  # Use np.nan here

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
