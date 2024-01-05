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

""" Easy to use benchling functions to fetch sequences and objects"""
import os
from benchlingapi import Session
import pandas as pd
import datetime
import Bio
from Bio.SeqFeature import SeqFeature
import pydna
from teemi.utils import (
    rename_dict_keys,
    split_based_on_keys,
    nest_dict,
    start_end_to_location,
)

########################################################
########## Initializing Benchling environmen t##########
########################################################

# Fetching API key from env file  - the function can find the env file in the root folder
import os
from dotenv import load_dotenv, find_dotenv

# find dotenv file
load_dotenv(find_dotenv())

# Initializing session
api_key = os.environ.get("API_KEY")

# Add your benchling home environment here
home_url = os.environ.get("HOME_url")
session = Session(api_key=api_key, home=home_url)


""" This part of the LIMS module is used for importing sequences and exporting sequences to benchling"""


def sequence_to_benchling(folder_name, oligo_name, oligo_bases, schema):
    """This function uploads sequences to Benchling.

    Parameters
    ----------
    folder_name : str
        The name of the folder in which to upload the sequence.
    oligo_name : str
        The name of the oligo.
    oligo_bases : str
        The base sequence of the oligo.
    schema : str
        The schema to use for the oligo. Should be one of the following:
        "Primer", "DNA Fragment", "Plasmid", "Gene", "gRNA", "Marker","Promoter",
        "Terminator", "Tag", "Origin of Replication".

    Returns
    -------
    None
    """

    folder = session.Folder.find_by_name(folder_name)
    dna = ""
    dna = session.DNASequence(
        name=oligo_name, bases=oligo_bases, folder_id=folder.id, is_circular=False
    )

    # save the dna to Benchling account
    dna.save()

    schemas = (
        "Primer",
        "DNA Fragment",
        "Plasmid",
        "Gene",
        "gRNA",
        "Marker",
        "Promoter",
        "Terminator",
        "Tag",
        "Origin of Replication",
    )
    if schema in schemas:
        dna.set_schema(schema)

    dna.register()


def from_benchling(bname: str, schema: str = ""):
    """Extract information of object on benchling.
    Parameters
    ----------
    bname : str
        The name of the object on Benchling.
    schema : str, optional
        The schema of the object, by default ""

    Returns
    -------
    object :
        The extracted object from Benchling.
    """
    # retrieve benchling sequence dict
    bench_dict = session.DNASequence.find_by_name(bname).dump()

    # Convert benchling sequence dict into Biopython/Genbank object
    ## rename keys from Benchling to Biopython/Genbank
    trans_dict = {
        "bases": "seq",
        "id": "id",
        "annotations": "features",
        "name": "name",
        "fields": "annotations",
    }
    trans_bench_dict = rename_dict_keys(bench_dict, trans_dict)

    ## Select most import information
    translated_bench_dict_sel, translated_bench_dict_other = split_based_on_keys(
        trans_bench_dict, trans_dict.values()
    )

    ## Add Benchling custom fields and isCircular to fields dict in Biopython annotations.
    translated_bench_dict_sel["annotations"].update(
        translated_bench_dict_other["customFields"]
    )
    translated_bench_dict_sel["annotations"].update(
        {"topology": translated_bench_dict_other["isCircular"]}
    )

    ##
    date = datetime.datetime.today().strftime("%d-%b-%Y").upper()
    # comment is renamed to commentary as Bio.IO.write provides an error when this key is in the annotation dict.
    comment = translated_bench_dict_sel["annotations"].pop("comment", None)

    translated_bench_dict_sel["annotations"].update(
        {
            "data_file_division": "PLN",
            "date": date,
            "molecule_type": "DNA",
            "location": "unknown",
            "commentary": comment,
        }
    )

    ## Convert values
    ### Seq
    translated_bench_dict_sel["seq"] = Bio.Seq.Seq(translated_bench_dict_sel["seq"])

    ### Features
    seq_length = len(translated_bench_dict_sel["seq"])

    #### Translate start end to compound locations object
    translated_bench_dict_sel["features"] = [
        start_end_to_location(feature_dict, seq_length)
        for feature_dict in translated_bench_dict_sel["features"]
    ]

    #### Add "label" feature.qualifier dict for CLC workbench visualization
    for feature in translated_bench_dict_sel["features"]:
        feature.update({"label": feature.get("name", "")})

    #### Move non Biopython features to qualifier subdict
    translated_bench_dict_sel["features"] = [
        nest_dict(
            feature_dict,
            first_order_keys=["location", "type", "strand"],
            key_for_nested_dict="qualifiers",
        )
        for feature_dict in translated_bench_dict_sel["features"]
    ]

    #### Create SeqFeatures
    translated_bench_dict_sel["features"] = [
        Bio.SeqFeature.SeqFeature(**feature_dict)
        for feature_dict in translated_bench_dict_sel["features"]
    ]

    seqRecord = Bio.SeqRecord.SeqRecord(**translated_bench_dict_sel)
    seqRecord.name = bname
    seqRecord = update_loc_vol_conc(seqRecord)

    if schema == "Primer":
        seqRecord = pydna.primer.Primer(seqRecord)

    # return object
    return seqRecord


def update_loc_vol_conc(seqRecord, DBpath: str = ""):
    """Update with location volume and concentration
    information downloaded from benchling if possible.

    Parameters
    ----------
    seqRecord : Bio.SeqRecord.SeqRecord
        The SeqRecord object to update.
    DBpath : str, optional
        The path to the csv file containing the location, volume and concentration information.

    Returns:
    --------
    seqRecord : Bio.SeqRecord.SeqRecord
        The updated SeqRecord object.
    """
    DB = pd.read_csv(DBpath)

    seqRecord.annotations["batches"] = []
    for i, row in DB.iterrows():
        if row["batchEntId"] == seqRecord.id:
            location = row["parentBoxPlateName"] + "_" + row["parentBoxPlatePos"]
            seqRecord.annotations["batches"].append(
                {
                    "box": row["parentBoxPlateName"],
                    "position": row["parentBoxPlatePos"],
                    "volume": int(row["volume"]),
                    "concentration": int(row["Concentration (ng/ul)"]),
                    "location": location,
                }
            )
        else:
            pass
    return seqRecord
