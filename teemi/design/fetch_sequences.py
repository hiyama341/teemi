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

from Bio import SeqIO
from Bio import Entrez
import requests as r
from io import StringIO

# intermine
# from __future__ import print_function
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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

        out_handle = open(out_file, "w")

        for i in range(0, len(list_of_acc_numbers)):
            Entrez.email = email
            handle = Entrez.efetch(
                db=db, id=list_of_acc_numbers[i], rettype="fasta", retmode="text"
            )
            out_handle.write(handle.read())
        out_handle.close()

    except:
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


def fetch_promoter(promoter_name: str):
    from intermine.webservice import Service

    """Retrieves a yeast promoter sequence from intermine.
    Parameters
    ----------
    promoter_name: str

    Returns
    -------
    promoter sequence : str
    """
    seq = ""
    service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")
    query = service.new_query("Gene")
    query.add_view(
        "secondaryIdentifier",
        "symbol",
        "length",
        "flankingRegions.direction",
        "flankingRegions.sequence.length",
        "flankingRegions.sequence.residues",
    )

    query.add_constraint("Gene", "LOOKUP", promoter_name, "S. cerevisiae", code="B")
    query.add_constraint("flankingRegions.direction", "=", "upstream", code="C")
    query.add_constraint("flankingRegions.distance", "=", "1.0kb", code="A")
    query.add_constraint("flankingRegions.includeGene", "=", "false", code="D")

    for row in query.rows():
        seq = row["flankingRegions.sequence.residues"]

    return seq


def fetch_multiple_promoters(List_of_promoter_names: list):
    """Retrieves a yeast promoter sequence from intermine.
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
        promoters_seq.description = "Defined as being 1kb upstream of the TSS and fetched through Intermines API"

        # Append to list
        LIST_OF_BIOrecord_objects.append(promoters_seq)

    return LIST_OF_BIOrecord_objects
