
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

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pydna.dseqrecord import Dseqrecord
import pandas as pd
from typing import Counter
import re
from collections import Counter
from Bio.Seq import Seq
from typing import Callable, Dict, Optional, List, Tuple
from Bio.SeqFeature import SeqFeature, FeatureLocation
from teemi.design.crispr_cas import (
    find_off_target_hits,
    find_sgrna_hits_cas9,
    revcomp,
    parse_genbank_record,
    SgRNAargs,
)


def filter_crispri_guides(args: SgRNAargs, hitframe: pd.DataFrame) -> pd.DataFrame:
    """
    Filter the sgRNAs based on the specified criteria.

    Parameters
    ----------
    args : SgRNAargs
        An object of class SgRNAargs containing the parameters for sgRNA filtering.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the filtered sgRNAs.
    """

    def exclude_rows_based_on_patterns(frame: pd.DataFrame, column: str, patterns: List[str]) -> pd.DataFrame:
        """Exclude rows from frame where frame[column] contains any of the provided patterns"""
        if not patterns:
            return frame
        if column not in frame.columns:
            raise ValueError(
                f"cannot filter on {column!r}: the sgRNA hit finders do not "
                f"produce a {column!r} column"
            )
        mask = frame[column].apply(lambda x: not any(pattern in x for pattern in patterns))
        return frame[mask]

    # TODO remove this since it filters too hard. 
    # Filter based on the range relative to TSS 
    filtered_frame = hitframe.loc[
        (hitframe['sgrna_loc'] >= -args.upstream_tss) & 
        (hitframe['sgrna_loc'] <= args.dwstream_tss)
    ]

    # Filter rows where gene_strand and sgrna_strand are different
    if args.target_non_template_strand:
        filtered_frame = filtered_frame[filtered_frame['gene_strand'] != filtered_frame['sgrna_strand']]

    # Apply other filters; each is a no-op when nothing was asked for
    filtered_frame = exclude_rows_based_on_patterns(
        filtered_frame, "pam", args.pam_remove
    )
    filtered_frame = exclude_rows_based_on_patterns(
        filtered_frame, "sgrna", args.sgrna_remove
    )
    filtered_frame = exclude_rows_based_on_patterns(
        filtered_frame, "downstream", args.downstream_remove
    )

    filtered_frame = filtered_frame.loc[filtered_frame.gc >= args.gc_lower]
    filtered_frame = filtered_frame.loc[filtered_frame.gc <= args.gc_upper]
    filtered_frame = filtered_frame.loc[filtered_frame.off_target_count <= args.off_target_upper]

    return filtered_frame



def find_sgrna_hits_cas9_crispri(
    record: Dseqrecord,
    strain_name: str,
    locus_tags: List[str],
    off_target_counter: Counter,
    off_target_seed: int,
    revcomp: callable,
    extension_to_promoter_region: int = 100,
) -> pd.DataFrame:
    """
    Find Cas9 sgRNA hits in each CDS and in the promoter region upstream of it.

    The guides are found with :func:`teemi.design.crispr_cas.find_sgrna_hits_cas9`,
    once for the annotated CDS and once for a synthetic feature covering the
    upstream region, so both get their positions from the same code. Upstream
    guides are reported relative to the start of the gene, i.e. with a negative
    ``sgrna_loc``, and a ``region`` column says which of the two a guide is in.

    Parameters
    ----------
    record : Dseqrecord
        The record to parse.
    strain_name : str
        Name reported in the ``strain_name`` column.
    locus_tags : List[str]
        Locus tags to look for, or ``["all"]`` for every CDS.
    off_target_counter : Counter
        Counter object containing the frequency of each off-target hit.
    off_target_seed : int
        The length of the off-target seed sequence to match.
    revcomp : callable
        Function to get the reverse complement of a sequence.
    extension_to_promoter_region : int
        How many base pairs upstream of the gene to search (default 100).

    Returns
    -------
    sgrna_df : pd.DataFrame
        A DataFrame of sgRNA hits information, with an added ``region`` column
        holding either ``"CDS"`` or ``"upstream"``.
    """
    upstream_len = extension_to_promoter_region

    # Guides inside the coding sequence
    sgrna_cds = find_sgrna_hits_cas9(
        record, strain_name, locus_tags, off_target_counter, off_target_seed, revcomp
    )
    sgrna_cds["region"] = "CDS"

    if upstream_len <= 0:
        return sgrna_cds

    # Guides in the promoter region, found through a synthetic upstream feature
    sgrna_upstream_frames = []
    sequence_length = len(record.seq)

    for feature in record.features:
        if feature.type != "CDS":
            continue
        locus_tag = feature.qualifiers.get("locus_tag", ["NA"])[0]
        if locus_tag not in locus_tags and "all" not in locus_tags:
            continue

        gene_strand = feature.location.strand
        if gene_strand == 1:
            upstream_start = max(0, feature.location.start - upstream_len)
            upstream_end = feature.location.start
        else:
            upstream_start = feature.location.end
            upstream_end = min(sequence_length, feature.location.end + upstream_len)
        if upstream_end <= upstream_start:
            continue

        upstream_tag = locus_tag + "_upstream"
        upstream_record = Dseqrecord(record.seq)
        upstream_record.features = [
            SeqFeature(
                FeatureLocation(upstream_start, upstream_end, strand=gene_strand),
                type="CDS",
                qualifiers={"locus_tag": [upstream_tag]},
            )
        ]

        upstream_df = find_sgrna_hits_cas9(
            upstream_record,
            strain_name,
            [upstream_tag],
            off_target_counter,
            off_target_seed,
            revcomp,
        )
        if upstream_df.empty:
            continue

        # Report upstream guides relative to the start of the gene
        upstream_df["sgrna_loc"] = upstream_df["sgrna_loc"] - (
            upstream_end - upstream_start
        )
        upstream_df["region"] = "upstream"
        upstream_df["locus_tag"] = locus_tag
        sgrna_upstream_frames.append(upstream_df)

    if not sgrna_upstream_frames:
        return sgrna_cds

    return pd.concat([sgrna_cds] + sgrna_upstream_frames, ignore_index=True)


def extract_sgRNAs_for_crispri(args: SgRNAargs) -> Tuple[pd.DataFrame, Counter, pd.DataFrame]:
    """
    Execute all three functions together to extract gene information, off-target hits, 
    and sgRNA hits from a given genbank file.

    Parameters
    ----------
    args : SgRNAargs
        An instance of the SgRNAargs class.

    Returns
    -------
    gene_df : pd.DataFrame
        A DataFrame of gene information containing locus tag, gene name, strand, start, end.
    off_target_counter : Counter
        Counter object containing the frequency of each off-target hit.
    sgrna_df : pd.DataFrame
        A DataFrame of sgRNA hits information containing genbank file path, locus tag, 
        gene name, strand, offset, position offset, GC content, sgrna, PAM, 
        downstream sequence, sgrna_pam_downstream, seed, off-target count.
    """
    # Extract gene information
    sequences = parse_genbank_record(args.dseqrecord)

    if "cas9" not in args.cas_type:
        raise ValueError(
            f"cas_type {args.cas_type!r} is not supported for CRISPRi; "
            "only 'cas9' guides can be designed here"
        )
    if "find" not in args.step:
        raise ValueError(
            "step must include 'find': the guides have to be found before "
            f"they can be filtered (got {args.step!r})"
        )

    # Find all potential off-target hits
    off_target_counter = find_off_target_hits(
        sequences, args.off_target_seed, cas_type=args.cas_type
    )

    # Find all potential sgRNA hits
    sgrna_df = find_sgrna_hits_cas9_crispri(
        args.dseqrecord,
        args.strain_name,
        args.locus_tag,
        off_target_counter,
        args.off_target_seed,
        revcomp,
        extension_to_promoter_region=args.extension_to_promoter_region,
    )

    # Sort sgrna_df by position relative to the start of the gene
    sgrna_df.sort_values(by="sgrna_loc", ascending=True, inplace=True)

    # Filter guides if 'filter' is in the steps
    if "filter" in args.step:
        sgrna_df = filter_crispri_guides(args, sgrna_df)

    return sgrna_df
