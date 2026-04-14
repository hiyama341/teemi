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

"""Easy to use functions to fetch sequences and objects from local csv database."""

import os

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


_TEXT_COLUMNS = {
    "name",
    "date",
    "row",
    "seq",
    "features",
    "location",
    "reference",
    "comments",
    "description",
}


def _normalize_database_dtypes(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Make CSV-backed text columns writable across pandas versions."""
    text_columns = {
        column: "object" for column in _TEXT_COLUMNS if column in dataframe.columns
    }
    if text_columns:
        dataframe = dataframe.astype(text_columns)
    return dataframe


def get_unique_id(path="../data/csv_database") -> int:
    """Make a single unique ID from the csv_database files."""
    files = os.listdir(path)
    csv_files = [csv for csv in files if csv.endswith("csv")]

    list_of_ids = []
    for df in csv_files:
        dataframe = pd.read_csv(path + "/" + df)
        list_of_ids += dataframe["ID"].to_list()

    arr = np.array(list_of_ids)
    new_data = arr[~np.isnan(arr)]

    if new_data.size == 0:
        return 10000
    return int(np.amax(new_data)) + 1


def add_sequences_to_dataframe(list_of_DNA: list, csv_database_as_df, index=0) -> None:
    """Add sequences to a local CSV database dataframe."""
    counter = 0
    for ds_dna in list_of_DNA:
        if index != 0:
            blank_row_index = index + counter
            counter += 1
        else:
            blank_row_bool = csv_database_as_df.loc[:, "ID"].isna()
            blank_row_index = [i + index for i, x in enumerate(blank_row_bool) if x][0]

        csv_database_as_df.loc[blank_row_index, "ID"] = int(ds_dna.id)
        csv_database_as_df.loc[blank_row_index, "description"] = ds_dna.description
        csv_database_as_df.loc[blank_row_index, "size"] = len(ds_dna.seq)
        csv_database_as_df.loc[blank_row_index, "seq"] = str(ds_dna.seq)
        csv_database_as_df.loc[blank_row_index, "date"] = str(
            pd.to_datetime("today").strftime("%m-%d-%Y")
        )
        csv_database_as_df.loc[blank_row_index, "name"] = ds_dna.name
        csv_database_as_df.loc[blank_row_index, "features"] = str(ds_dna.features)
        csv_database_as_df.loc[blank_row_index, "concentration"] = ds_dna.annotations[
            "batches"
        ][0]["concentration"]
        csv_database_as_df.loc[blank_row_index, "volume"] = ds_dna.annotations[
            "batches"
        ][0]["volume"]
        csv_database_as_df.loc[blank_row_index, "location"] = ds_dna.annotations[
            "batches"
        ][0]["location"]
        csv_database_as_df.loc[blank_row_index, "comments"] = ds_dna.annotations[
            "comments"
        ]
        csv_database_as_df.loc[blank_row_index, "reference"] = ds_dna.annotations[
            "reference"
        ]


def get_plate(plate_number: int, csv_database_as_df):
    """Return the plate from a specified CSV database dataframe."""
    return csv_database_as_df.loc[csv_database_as_df["plate"] == plate_number]


def get_box(box_number: int, csv_database_as_df):
    """Return the box from a specified CSV database dataframe."""
    return csv_database_as_df.loc[csv_database_as_df["box"] == box_number]


def add_unique_ids(list_of_parts: list, path="../data/csv_database") -> None:
    """Add unique ids to a list of SeqRecords."""
    for i in range(len(list_of_parts)):
        unique_id = get_unique_id(path) + i
        list_of_parts[i].id = str(unique_id)


def add_annotations(
    list_of_parts: list,
    concentration: float = 0.0,
    reference: str = "",
    volume: float = 0.0,
    comments: str = "",
    location: str = "",
) -> list:
    """Add the necessary annotations to a list of SeqRecord objects."""
    for annotations in list_of_parts:
        annotations.annotations = {
            "reference": reference,
            "comments": comments,
            "batches": [
                {
                    "location": location,
                    "volume": float(volume),
                    "concentration": float(concentration),
                }
            ],
        }
    return list_of_parts


def update_database(dataframe, which_database: str, path="../data/csv_database/") -> None:
    """Update the chosen CSV database."""
    dataframe.to_csv(path + which_database + ".csv", index=False)


def get_dna_from_plate_name(
    name: str,
    database_name: str,
    database_path="../data/csv_database/",
    genbank_files_path="../data/genbank_files/",
    genbank=False,
) -> SeqRecord:
    """Fetch DNA based on the name from the plate database of choice."""
    path = database_path + database_name + ".csv"
    dataframe = pd.read_csv(path)
    found_the_record_df = dataframe.loc[dataframe["name"] == name].reset_index()

    if genbank:
        ID_from_df = str(int(found_the_record_df.loc[0, "ID"]))
        Record = SeqIO.read(genbank_files_path + ID_from_df + ".gb", format="gb")
    else:
        Record = SeqRecord(Seq(str(found_the_record_df.loc[0, "seq"])))
        Record.id = str(found_the_record_df.loc[0, "ID"])
        Record.name = found_the_record_df.loc[0, "name"]
        Record.description = found_the_record_df.loc[0, "description"]

    Record.annotations = {
        "plate": found_the_record_df.loc[0, "plate"],
        "row": found_the_record_df.loc[0, "row"],
        "col": found_the_record_df.loc[0, "col"],
        "batches": [
            {
                "location": str(found_the_record_df.loc[0, "location"])
                + "_"
                + str(found_the_record_df.loc[0, "plate"])
                + "_"
                + str(found_the_record_df.loc[0, "row"])
                + str(found_the_record_df.loc[0, "col"]),
                "volume": float(found_the_record_df.loc[0, "volume"]),
                "concentration": float(found_the_record_df.loc[0, "concentration"]),
            }
        ],
    }
    return Record


def get_dna_from_box_name(
    name: str,
    database_name: str,
    database_path="../data/csv_database/",
    genbank_files_path="../data/genbank_files/",
    genbank=False,
) -> SeqRecord:
    """Fetch DNA based on the name from the box database of choice."""
    path = database_path + database_name + ".csv"
    dataframe = pd.read_csv(path)
    found_the_record_df = dataframe.loc[dataframe["name"] == name].reset_index()

    if genbank:
        ID_from_df = str(int(found_the_record_df.loc[0, "ID"]))
        Record = SeqIO.read(genbank_files_path + ID_from_df + ".gb", format="gb")
    else:
        Record = SeqRecord(Seq(str(found_the_record_df.loc[0, "seq"])))
        Record.id = str(found_the_record_df.loc[0, "ID"])
        Record.name = found_the_record_df.loc[0, "name"]
        Record.description = found_the_record_df.loc[0, "description"]

    Record.annotations = {
        "box": found_the_record_df.loc[0, "box"],
        "row": found_the_record_df.loc[0, "row"],
        "col": found_the_record_df.loc[0, "col"],
        "batches": [
            {
                "location": str(found_the_record_df.loc[0, "location"])
                + "_"
                + str(found_the_record_df.loc[0, "box"])
                + "_"
                + str(found_the_record_df.loc[0, "row"])
                + str(found_the_record_df.loc[0, "col"]),
                "volume": float(found_the_record_df.loc[0, "volume"]),
                "concentration": float(found_the_record_df.loc[0, "concentration"]),
            }
        ],
    }
    return Record


def get_database(name: str, path="../data/csv_database/"):
    """Fetch the csv database as a pandas dataframe."""
    try:
        dataframe = pd.read_csv(path + name + ".csv", index_col=False)
        return _normalize_database_dtypes(dataframe)
    except Exception:
        print("Couldnt find that databse. Hack: Dont add csv extention.")


def change_row(row_index: int, csv_database_as_df, biopython_object):
    """Insert a biopython object into the database at a specific index."""
    csv_database_as_df.loc[row_index, "ID"] = int(biopython_object.id)
    csv_database_as_df.loc[row_index, "description"] = biopython_object.description
    csv_database_as_df.loc[row_index, "size"] = len(biopython_object.seq)
    csv_database_as_df.loc[row_index, "seq"] = str(biopython_object.seq)
    csv_database_as_df.loc[row_index, "date"] = str(
        pd.to_datetime("today").strftime("%m-%d-%Y")
    )
    csv_database_as_df.loc[row_index, "name"] = biopython_object.name
    csv_database_as_df.loc[row_index, "features"] = str(biopython_object.features)
    csv_database_as_df.loc[row_index, "concentration"] = biopython_object.annotations[
        "batches"
    ][0]["concentration"]
    csv_database_as_df.loc[row_index, "volume"] = biopython_object.annotations[
        "batches"
    ][0]["volume"]
    csv_database_as_df.loc[row_index, "comments"] = biopython_object.annotations[
        "comments"
    ]
    csv_database_as_df.loc[row_index, "reference"] = biopython_object.annotations[
        "reference"
    ]
    return csv_database_as_df


def delete_row_df(row_index, which_df):
    """Delete a row in the database without changing the names."""
    coloumns = [
        "ID",
        "name",
        "size",
        "seq",
        "concentration",
        "features",
        "location",
        "reference",
        "volume",
        "comments",
        "description",
        "date",
    ]
    which_df.loc[row_index, coloumns] = np.nan
