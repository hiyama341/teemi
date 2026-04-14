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

"""Easy to use Benchling functions to fetch sequences and objects."""

import datetime
import os

import Bio
import pandas as pd
import pydna
from benchlingapi import Session
from dotenv import find_dotenv, load_dotenv

from teemi.utils import (
    nest_dict,
    rename_dict_keys,
    split_based_on_keys,
    start_end_to_location,
)


load_dotenv(find_dotenv())

api_key = os.environ.get("API_KEY")
home_url = os.environ.get("HOME_url")
session = Session(api_key=api_key, home=home_url)


def sequence_to_benchling(folder_name, oligo_name, oligo_bases, schema):
    """Upload sequences to Benchling."""
    folder = session.Folder.find_by_name(folder_name)
    dna = session.DNASequence(
        name=oligo_name, bases=oligo_bases, folder_id=folder.id, is_circular=False
    )
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
    """Extract information of an object on Benchling."""
    bench_dict = session.DNASequence.find_by_name(bname).dump()

    trans_dict = {
        "bases": "seq",
        "id": "id",
        "annotations": "features",
        "name": "name",
        "fields": "annotations",
    }
    trans_bench_dict = rename_dict_keys(bench_dict, trans_dict)

    translated_bench_dict_sel, translated_bench_dict_other = split_based_on_keys(
        trans_bench_dict, trans_dict.values()
    )

    translated_bench_dict_sel["annotations"].update(
        translated_bench_dict_other["customFields"]
    )
    translated_bench_dict_sel["annotations"].update(
        {"topology": translated_bench_dict_other["isCircular"]}
    )

    date = datetime.datetime.today().strftime("%d-%b-%Y").upper()
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

    translated_bench_dict_sel["seq"] = Bio.Seq.Seq(translated_bench_dict_sel["seq"])
    seq_length = len(translated_bench_dict_sel["seq"])

    translated_bench_dict_sel["features"] = [
        start_end_to_location(feature_dict, seq_length)
        for feature_dict in translated_bench_dict_sel["features"]
    ]
    for feature in translated_bench_dict_sel["features"]:
        feature.update({"label": feature.get("name", "")})
    translated_bench_dict_sel["features"] = [
        nest_dict(
            feature_dict,
            first_order_keys=["location", "type", "strand"],
            key_for_nested_dict="qualifiers",
        )
        for feature_dict in translated_bench_dict_sel["features"]
    ]
    translated_bench_dict_sel["features"] = [
        Bio.SeqFeature.SeqFeature(**feature_dict)
        for feature_dict in translated_bench_dict_sel["features"]
    ]

    seqRecord = Bio.SeqRecord.SeqRecord(**translated_bench_dict_sel)
    seqRecord.name = bname
    seqRecord = update_loc_vol_conc(seqRecord)

    if schema == "Primer":
        seqRecord = pydna.primer.Primer(seqRecord)

    return seqRecord


def update_loc_vol_conc(seqRecord, DBpath: str = ""):
    """Update location, volume, and concentration information from a CSV export."""
    DB = pd.read_csv(DBpath)

    seqRecord.annotations["batches"] = []
    for _, row in DB.iterrows():
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
    return seqRecord
