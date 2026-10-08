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

""" This part of the design module is used fetching sequences"""

import time
from Bio import SeqIO
from Bio import Entrez
import requests as r
from io import StringIO

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature


def concatenate_genbank_records(
    records: list, record_id: str = "combined_record", record_name: str = "combined_record"
):
    """Concatenate multiple GenBank records into one sequence record.

    Parameters
    ----------
    records : list of Bio.SeqRecord.SeqRecord
        Ordered list of chromosome or contig records.
    record_id : str, optional
        Identifier for the concatenated record.
    record_name : str, optional
        Name for the concatenated record.

    Returns
    -------
    Bio.SeqRecord.SeqRecord
        A new record with concatenated sequence and shifted feature locations.
    """
    concatenated_sequence = Seq("")
    concatenated_features = []
    sequence_offset = 0

    for record in records:
        concatenated_sequence += record.seq

        for feature in record.features:
            shifted_feature = SeqFeature(
                location=feature.location._shift(sequence_offset),
                type=feature.type,
                id=feature.id,
                qualifiers=feature.qualifiers.copy(),
            )
            concatenated_features.append(shifted_feature)

        sequence_offset += len(record.seq)

    concatenated_record = SeqRecord(
        concatenated_sequence,
        id=record_id,
        name=record_name,
        description=f"Concatenated record built from {len(records)} GenBank records",
        features=concatenated_features,
    )

    return concatenated_record


def retrieve_sequences_from_ncbi(
    list_of_acc_numbers: list, out_file: str, db="protein"
):
    """Retrieves sequences from ncbi.
    Parameters
    ----------
    list_of_acc_numbers: list
        list_of_acc_numbers such as: ['Q05001', 'Q1PQK4','Q9SB48' ,'AFX82679']

    Returns
    -------
    A fasta file with your sequences
    """
    try:
        email = "youremail@gmail.com"

        with open(out_file, "w") as out_handle:
            for i in range(0, len(list_of_acc_numbers)):
                Entrez.email = email
                handle = Entrez.efetch(
                    db=db, id=list_of_acc_numbers[i], rettype="fasta", retmode="text"
                )
                out_handle.write(handle.read())

    except Exception:
        print(
            "An exception occurred, please double-check your accession numbers or connection"
        )


def read_fasta_files(path):
    """Reads FASTA files.
    Parameters
    ----------
    path: str
        path to the fasta file you want to read.

    Returns
    -------
    list of Bio.SeqRecord.SeqRecord
    """

    ncbi_hits = []
    for seq_record in SeqIO.parse(path, format="fasta"):
        ncbi_hits.append(seq_record)

    return ncbi_hits


def read_genbank_files(path):
    """Reads single Genbank files.
    Parameters
    ----------
    path: str
        path to the genbank file you want to read.

    Returns
    -------
    list of Bio.SeqRecord.SeqRecord
    """

    ncbi_hits = []
    for seq_record in SeqIO.parse(path, format="gb"):
        ncbi_hits.append(seq_record)

    return ncbi_hits


def retrieve_sequences_from_PDB(query: list):
    """Retrieves sequences from PDB.
    Parameters
    ----------
    query: list
        list of accession numbers in the form of strings

    Returns
    -------
    list of Bio.SeqRecord.SeqRecord
    """
    list_of_protein_seqs = []

    for q in query:
        cID = q

        baseUrl = "http://www.uniprot.org/uniprot/"
        currentUrl = baseUrl + cID + ".fasta"
        response = r.post(currentUrl)
        cData = "".join(response.text)

        Seq = StringIO(cData)
        Protein_sequence = list(SeqIO.parse(Seq, "fasta"))
        list_of_protein_seqs.append(Protein_sequence)

    return list_of_protein_seqs


ENSEMBL_REST_URL = "https://rest.ensembl.org"


def _ensembl_get(path: str, retries: int = 5):
    """GET a JSON resource from the Ensembl REST API, backing off when rate limited."""
    for attempt in range(retries):
        response = r.get(
            ENSEMBL_REST_URL + path,
            headers={"Content-Type": "application/json"},
            timeout=30,
        )
        if response.status_code != 429 and response.status_code < 500:
            break
        time.sleep(float(response.headers.get("Retry-After", 2**attempt)))
    return response


def fetch_promoter(promoter_name: str):
    """Retrieves a yeast promoter sequence, defined as the 1 kb upstream of the gene.

    The sequence used to come from YeastMine, which SGD retired in July 2024.
    It is now fetched from the Ensembl REST API (S. cerevisiae S288C), which
    returns the same 1 kb upstream flanking region.

    Parameters
    ----------
    promoter_name: str
        standard (e.g. ``"CYC1"``) or systematic (e.g. ``"YJR048W"``) gene name

    Returns
    -------
    promoter sequence : str
        empty if the gene is not found
    """
    gene = _ensembl_get(f"/lookup/symbol/saccharomyces_cerevisiae/{promoter_name}")
    if not gene.ok:
        # systematic names are Ensembl stable IDs
        gene = _ensembl_get(f"/lookup/id/{promoter_name}")
    if not gene.ok:
        return ""
    gene = gene.json()

    if gene["strand"] == 1:
        start, end = max(gene["start"] - 1000, 1), gene["start"] - 1
    else:
        start, end = gene["end"] + 1, gene["end"] + 1000
    region = f"{gene['seq_region_name']}:{start}..{end}:{gene['strand']}"

    sequence = _ensembl_get(f"/sequence/region/saccharomyces_cerevisiae/{region}")
    sequence.raise_for_status()
    return sequence.json()["seq"]


def fetch_multiple_promoters(List_of_promoter_names: list):
    """Retrieves yeast promoter sequences (1 kb upstream), see fetch_promoter.
    Parameters
    ----------
    List_of_promoter_names: list
        list of strings of promoter names fx : ['YAR035C-A', 'YGR067C', 'JEN1', 'YNR034W-A', 'ACH1']

    Returns
    -------
    list of Bio.SeqRecord.SeqRecord

    """
    # #initializing
    LIST_OF_BIOrecord_objects = []

    for i in range(0, len(List_of_promoter_names)):
        # fetching the seqs
        promoters_seq = SeqRecord(Seq(fetch_promoter(List_of_promoter_names[i])))
        promoters_seq.name = str(List_of_promoter_names[i]) + " Promoter"
        promoters_seq.id = str(List_of_promoter_names[i])
        promoters_seq.description = "Defined as being 1kb upstream of the TSS and fetched through the Ensembl REST API"

        # Append to list
        LIST_OF_BIOrecord_objects.append(promoters_seq)

    return LIST_OF_BIOrecord_objects
