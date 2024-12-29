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
import pandas as pd
from typing import Counter
import re
from collections import Counter
from Bio.Seq import Seq
from typing import Callable, Dict, Optional, List, Tuple
from pydna.dseqrecord import Dseqrecord
from Bio.SeqRecord import SeqRecord
import pandas as pd
from typing import Tuple, List
from Bio import SeqIO


class SgRNAargs:
    """
    A class used to store the parameters for sgRNA generation and filtering.

    Attributes
    ----------
    dseqrecord : Dseqrecord
        Dseqrecord containing the genomic data.
    locus_tag : list of str
        List of locus tags (i.e., genes) to be analyzed.
    cas_type : str
        The type of the CRISPR-Cas system used ('cas9', 'cas12a', or 'cas13').
    step : list of str, optional
        Pipeline components to execute (default is ['find', 'filter']).
    downstream : int, optional
        Number of downstream nucleotides to consider when filtering sgRNAs (default is 15).
    off_target_seed : int, optional
        How many PAM-adjacent bases to use as seed-sequence for off-target hits (default is 9).
    pam_remove : list of str, optional
        List of PAM sequences to filter (default is None).
    downstream_remove : list of str, optional
        List of downstream sequences to filter (default is None).
    sgrna_remove : list of str, optional
        List of specific sgRNA targets to filter (default is None).
    gc_upper : float, optional
        GC cutoff filter (maximum-value) (default is 1).
    gc_lower : float, optional
        GC cutoff filter (minimum-value) (default is 0).
    off_target_upper : int, optional
        Maximum number of seed nucleotides allowed for off-target (default is 10).
    extension_to_promoter_region : int, optional
        Additional parameter to extend the promoter region (default is 0).
    upstream_tss : int, optional
        Number of upstream nucleotides to consider from the TSS (default is 100).
    dwstream_tss : int, optional
        Number of downstream nucleotides to consider from the TSS (default is 100).
    target_non_template_strand : bool, optional
        Whether to target the non-template strand (default is True).
    strain_name : str
        Strain name extracted from the Dseqrecord.
    """

    def __init__(self, dseqrecord: Dseqrecord, locus_tag: List[str], cas_type: str = 'cas9',
                 step: List[str] = ['find', 'filter'], 
                 downstream: int = 15, 
                 off_target_seed: int = 9, 
                 pam_remove: Optional[List[str]] = None, 
                 downstream_remove: Optional[List[str]] = None, 
                 sgrna_remove: Optional[List[str]] = None, 
                 gc_upper: float = 1, gc_lower: float = 0, 
                 off_target_upper: int = 10, 
                 extension_to_promoter_region: int = 0, 
                 upstream_tss: int = 100, 
                 dwstream_tss: int = 100, 
                 target_non_template_strand = True,

                 ):
        
        if not isinstance(dseqrecord, Dseqrecord):
            raise ValueError("Input must be an instance of Dseqrecord.")
        
        self.dseqrecord = dseqrecord
        self.locus_tag = locus_tag if isinstance(locus_tag, list) else [locus_tag]
        self.step = step
        self.downstream = downstream
        self.off_target_seed = off_target_seed
        self.pam_remove = pam_remove
        self.downstream_remove = downstream_remove
        self.sgrna_remove = sgrna_remove
        self.gc_upper = gc_upper
        self.gc_lower = gc_lower
        self.off_target_upper = off_target_upper
        self.cas_type = cas_type
        self.extension_to_promoter_region = extension_to_promoter_region
        self.upstream_tss = upstream_tss
        self.dwstream_tss = dwstream_tss
        self.target_non_template_strand = target_non_template_strand

        # Extract strain name from the Dseqrecord
        self.strain_name = self.dseqrecord.id

# TODO incorporate this into the code? 
def revcomp(x: str) -> str:
    """
    Compute the reverse complement of a DNA sequence.

    Parameters
    ----------
    x : str
        A string representing a DNA sequence.

    Returns
    -------
    str
        The reverse complement of the input DNA sequence.
    """
    return str(Seq(x).reverse_complement())


def parse_genbank_record(dseqrecord: SeqRecord) -> Tuple[List[str], pd.DataFrame]:
    """
    Parse a Dseqrecord and extract gene information and sequences.

    Parameters
    ----------
    dseqrecord : Dseqrecord
        The Dseqrecord to parse.

    Returns
    -------
    sequences : List[str]
        A list of sequences of positive and negative strands (reverse complement).
    gene_df : pd.DataFrame
        A DataFrame of gene information containing locus tags, gene names, strands, starts, and ends.
    """
    sequences = list()

    # Extract gene information and sequences from the Dseqrecord
    sequences.append(str(dseqrecord.seq))  # Store sequence of positive strand
    sequences.append(str(dseqrecord.reverse_complement().seq))  # Store sequence of negative strand (reverse complement)

    return sequences


def find_off_target_hits(sequences: List[str], off_target_seed: int, cas_type:str) -> Counter:
    """
    Find all potential off-target hits in both strands and count their frequencies.

    Parameters
    ----------
    sequences : List[str]
        A list of sequences in which to find off-target hits.
    off_target_seed : int
        The length of the off-target seed sequence to match.
    cas_type : str
        The type of the CRISPR-Cas system used ('cas9', 'cas12a', or 'cas3').

    Returns
    -------
    off_target_counter : Counter
        A Counter object containing the frequency of each off-target hit.

    """
    # Initialize list to store off-target hits
    off_target_hits = list()
    # Define PAM pattern based on cas_type
    pam_patterns = {
        'cas3': r"TTC",
        'cas9': r"(?=GG)",
        'cas12a': r"TTT[ACG]" }
    
    # Find all potential off-target hits in both strands
    for contig in sequences:
        for match in re.finditer(pam_patterns[cas_type], contig):
            if cas_type == 'cas12a':
                # For Cas12a, the seed sequence is immediately downstream of the PAM
                seed_sequence = contig[
                    match.end():(match.end() + off_target_seed)
                ]
            elif cas_type == 'cas3':
                # For Cas3, the seed sequence is immediately downstream of the PAM
                seed_sequence = contig[
                    match.end():(match.end() + off_target_seed)
                ]
            else: # CAS9
                seed_sequence = contig[
                    (match.start() - off_target_seed-1):(match.start() -1)
    ]    
            # Count the frequency of each off-target hit
            off_target_hits.append(seed_sequence)


    # Count the frequency of each off-target hit
    off_target_counter = Counter(off_target_hits)

    return off_target_counter


def find_sgrna_hits_cas9(record: Dseqrecord, strain_name: str, locus_tags: List[str], off_target_counter: Counter, off_target_seed: int, revcomp: callable) -> pd.DataFrame:
    """
    Parse a Dseqrecord object to find sgRNA hits.

    Parameters
    ----------
    record : Dseqrecord
        Dseqrecord object from pydna.
    locus_tags : List[str]
        List of locus tags to find in the record.
    off_target_counter : Counter
        Counter object containing the frequency of each off-target hit.
    off_target_seed : int
        The length of the off-target seed sequence to match.
    downstream : int
        The number of downstream base pairs to include in the downstream sequence.
    revcomp : callable
        Function to get the reverse complement of a sequence.

    Returns
    -------
    sgrna_df : pd.DataFrame
        A DataFrame of sgRNA hits information.

    """
    # Initialize list to store sgRNA hits
    sgrna_hits = list()

    # Counter for genes processed
    gene_counter = 1

    # Cas9-specific parameters
    pam_pattern = r"(?=CC)"
    protospacer_len = 20
    pam_len = 3

    # Parse the Dseqrecord again to find sgRNAs
    for feature in record.features:
        if feature.type == "CDS":  # Check if feature is a coding sequence
            locus_tag = feature.qualifiers.get("locus_tag", ["NA"])[0]
            if locus_tag in locus_tags or "all" in locus_tags:  # Check if gene tag is in the list of locus tags
                # Extract the sequence of the coding sequence
                coding_sequence = str(feature.extract(record.seq))
                coding_sequence_revcomp = revcomp(coding_sequence)

                # Extract location of the feature
                location = feature.location
                gene_strand = feature.location.strand

                # Find potential sgRNAs in both the coding sequence and its reverse complement
                for sequence in [(-1, coding_sequence), (1, coding_sequence_revcomp)]:

                    # Increment gene counter
                    gene_counter = gene_counter + 1
                    for match in re.finditer(pam_pattern, sequence[1]):

                        sgrna_pam = str(sequence[1][match.start():(match.start() + pam_len + protospacer_len)])
            
                        # Get reverse complement of the sgRNA and PAM sequence
                        sgrna = revcomp(sgrna_pam)[0:protospacer_len]
                        pam = revcomp(sgrna_pam)[protospacer_len:protospacer_len+pam_len]

                        if not sgrna:
                            print(f"No sgRNA found for locus tag {locus_tag}. Skipping to next locus tag.")
                            continue  # This skips the rest of the current iteration and moves to the next feature

                        if len(sgrna) != protospacer_len:  # Check if sgRNA is exactly 23 nt long
                            print(f"sgRNA generated were outside the designated border in {locus_tag}. To incorporate this extent borders. Skipping to next locus tag.")
                            continue  # This skips the rest of the current iteration and moves to the next feature
                        
                        if len(pam) != pam_len:  # Check if sgRNA is exactly 23 nt long
                            print(f"Pam was found outside designated locus_tag: {locus_tag}. To incorporate this extent borders. Skipping to next locus tag.")
                            continue  # This skips the rest of the current iteration and moves to the next feature

                        # Calculate GC content of the sgRNA
                        gc_content = len([base for base in sgrna if base in ["C", "G"]]) / len(sgrna) if len(sgrna) > 0 else 0

                        # Calculate genomic location of the sgRNA depending on the strand
                        if sequence[0] == 1:
                            genome_location = (int(location.start)) +1
                            position_sgrna = len(sequence[1]) -match.start()-3
                            strand_sgrna = 1

                        if sequence[0] == -1:
                            genome_location = int(location.start)+1
                            position_sgrna = match.end() + protospacer_len+3   
                            strand_sgrna = -1

                        sgrna_seed = sgrna[(protospacer_len - off_target_seed):protospacer_len]

                        # Get number of off-target hits for the seed sequence
                        off_target_count = (
                            off_target_counter[sgrna_seed] - 1 # if it turns out to be minus 1 it has not found the seed and there is a mistake.
                        )
                        # Store sgRNA hit
                        sgrna_hits.append(
                            (strain_name, 
                                locus_tag, 
                                genome_location, 
                                gene_strand, 
                                strand_sgrna,
                                position_sgrna,
                                gc_content,
                                pam,
                                sgrna,
                                sgrna_seed,
                                off_target_count,
                            )
                        )

    # Convert sgRNA hits into a DataFrame
    sgrna_df = pd.DataFrame(
        sgrna_hits,
        columns=['strain_name',
            "locus_tag",
            "gene_loc",
            "gene_strand",
            "sgrna_strand",
            "sgrna_loc",
            "gc",
            "pam",
            "sgrna",
            "sgrna_seed_sequence",
            "off_target_count",
        ]
    )
    return sgrna_df

#TODO input here is a filepath. WOuld be nice with consistency i.e. using a Dseqrecord
def find_sgrna_hits_cas12a(filepath: str, strain_name: str, locus_tags: List[str], off_target_counter: Counter, off_target_seed: int, revcomp: callable) -> pd.DataFrame:
    """
    Parse a genbank file to find sgRNA hits.

    Parameters
    ----------
    filepath : str
        The path to the genbank file to parse.
    locus_tags : List[str]
        List of locus tags to find in the genbank file.
    off_target_counter : Counter
        Counter object containing the frequency of each off-target hit.
    off_target_seed : int
        The length of the off-target seed sequence to match.
    downstream : int
        The number of downstream base pairs to include in the downstream sequence.
    revcomp : callable
        Function to get the reverse complement of a sequence.

    Returns
    -------
    sgrna_df : pd.DataFrame
        A DataFrame of sgRNA hits information.

    """
    # Initialize list to store sgRNA hits
    sgrna_hits = list()

    # Counter for genes processed
    gene_counter = 1

    # Cas12a-specific parameters
    pam_pattern = r"TTT[ACG]"
    protospacer_len = 23
    pam_len = 4

    # Parse the genbank file again to find sgRNAs
    for record in SeqIO.parse(filepath, "gb"):
        if record.features:
            for feature in record.features:
                if feature.type == "CDS":  # Check if feature is a coding sequence
                    locus_tag = feature.qualifiers.get("locus_tag", ["NA"])[0]
                    if locus_tag in locus_tags or "all" in locus_tags:  # Check if gene tag is in the list of locus tags
                        # Extract the sequence of the coding sequence
                        coding_sequence = str(feature.extract(record.seq))
                        coding_sequence_revcomp = revcomp(coding_sequence)

                        # Extract location of the feature
                        location = feature.location
                        gene_strand = feature.location.strand

                        # Find potential sgRNAs in both the coding sequence and its reverse complement
                        for sequence in [(-1, coding_sequence), (1, coding_sequence_revcomp)]:
                            
                            # Counter for sgRNAs found in the current gene
                            sgrna_counter = 0

                            # Increment gene counter
                            gene_counter = gene_counter + 1

                            for match in re.finditer(pam_pattern, sequence[1]):
                                
                                # Increment sgRNA counter
                                sgrna_counter = sgrna_counter + 1

                                # Extract the sgRNA sequence following the PAM
                                sgrna_start = match.start() - protospacer_len  
                                pam_start = sgrna_start + protospacer_len  
                                sgrna_with_pam = str(sequence[1][sgrna_start+protospacer_len:pam_start + pam_len+ protospacer_len])

                                # We're dealing with the reverse complement
                                sgrna = sgrna_with_pam[pam_len:]  # Extract sgRNA
                                pam = sgrna_with_pam[:pam_len]  # Extract pam

                                if not sgrna:
                                    print(f"No sgRNA found for locus tag {locus_tag}. Skipping to next locus tag.")
                                    continue  # This skips the rest of the current iteration and moves to the next feature

                                if len(sgrna) != protospacer_len:  # Check if sgRNA is exactly 23 nt long
                                    print(f"sgRNA generated were too small{locus_tag}. To incorporate this extent borders. Skipping to next locus tag.")
                                    continue  # This skips the rest of the current iteration and moves to the next feature

                                # Calculate GC content of the sgRNA
                                gc_content = len([base for base in sgrna if base in ["C", "G"]]) / len(sgrna) if len(sgrna) > 0 else 0
                                
                                # Calculate genomic location of the sgRNA depending on the strand
                                if sequence[0] == 1:
                                    genome_location = (int(location.start)) +1
                                    position_sgrna = len(sequence[1]) -match.start()-3
                                    strand_sgrna = sequence[0]
                                if sequence[0] == -1:
                                    genome_location = int(location.start)+1
                                    position_sgrna = match.end() + protospacer_len+3   
                                    strand_sgrna = sequence[0]

                                # For Cas12a, extract the seed sequence from the beginning of the sgRNA
                                sgrna_seed = sgrna[:off_target_seed]

                                # Get number of off-target hits for the seed sequence
                                off_target_count = (
                                    off_target_counter[sgrna_seed] - 1
                                )

                                # Store sgRNA hit
                                sgrna_hits.append(
                                    (strain_name, 
                                        locus_tag, 
                                        genome_location, 
                                        gene_strand, 
                                        strand_sgrna,
                                        position_sgrna,
                                        gc_content,
                                        pam,
                                        sgrna,
                                        sgrna_seed,
                                        off_target_count,
                                    )
                                )

    # Convert sgRNA hits into a DataFrame
    sgrna_df = pd.DataFrame(
        sgrna_hits,
        columns=['strain_name',
            "locus_tag",
            "gene_loc",
            "gene_strand",
            "sgrna_strand",
            "sgrna_loc",
            "gc",
            "pam",
            "sgrna",
            "sgrna_seed_sequence",
            "off_target_count",
        ]
    )
    return sgrna_df




def find_sgrna_hits_cas3(record: Dseqrecord, strain_name: str, locus_tags: List[str], off_target_counter: Counter, off_target_seed: int, revcomp: callable) -> pd.DataFrame:
    """
    Parse a genbank file to find sgRNA hits.

    Parameters
    ----------
    record: Dseqrecord
        The path to the genbank file to parse.
    locus_tags : List[str]
        List of locus tags to find in the genbank file.
    off_target_counter : Counter
        Counter object containing the frequency of each off-target hit.
    off_target_seed : int
        The length of the off-target seed sequence to match.
    downstream : int
        The number of downstream base pairs to include in the downstream sequence.
    revcomp : callable
        Function to get the reverse complement of a sequence.

    Returns
    -------
    sgrna_df : pd.DataFrame
        A DataFrame of sgRNA hits information.

    """
    # Initialize list to store sgRNA hits
    sgrna_hits = list()

    # Counter for genes processed
    gene_counter = 1

    # Cas12a-specific parameters
    pam_pattern = r"TTC"
    protospacer_len = 34
    pam_len = 3

    # Parse the genbank file again to find sgRNAs

    for feature in record.features:
        if feature.type == "CDS":  # Check if feature is a coding sequence
            locus_tag = feature.qualifiers.get("locus_tag", ["NA"])[0]
            if locus_tag in locus_tags or "all" in locus_tags:  # Check if gene tag is in the list of locus tags
                # Extract the sequence of the coding sequence
                coding_sequence = str(feature.extract(record.seq))
                coding_sequence_revcomp = revcomp(coding_sequence)

                # Extract location of the feature
                location = feature.location
                gene_strand = feature.location.strand

                # Find potential sgRNAs in both the coding sequence and its reverse complement
                for sequence in [(-1, coding_sequence), (1, coding_sequence_revcomp)]:
                    
                    # Counter for sgRNAs found in the current gene
                    sgrna_counter = 0

                    # Increment gene counter
                    gene_counter = gene_counter + 1

                    for match in re.finditer(pam_pattern, sequence[1]):
                        
                        # Increment sgRNA counter
                        sgrna_counter = sgrna_counter + 1

                        # Extract the sgRNA sequence following the PAM
                        sgrna_start = match.start() - protospacer_len  
                        pam_start = sgrna_start + protospacer_len  
                        sgrna_with_pam = str(sequence[1][sgrna_start+protospacer_len:pam_start + pam_len+ protospacer_len])

                        # We're dealing with the reverse complement
                        sgrna = sgrna_with_pam[pam_len:]  # Extract sgRNA
                        pam = sgrna_with_pam[:pam_len]  # Extract pam

                        if not sgrna:
                            print(f"No sgRNA found for locus tag {locus_tag}. Skipping to next locus tag.")
                            continue  # This skips the rest of the current iteration and moves to the next feature

                        if len(sgrna) != protospacer_len:  # Check if sgRNA is exactly 23 nt long
                            print(f"sgRNA generated were too small{locus_tag}. To incorporate this extent borders. Skipping to next locus tag.")
                            continue  # This skips the rest of the current iteration and moves to the next feature

                        if len(pam) != pam_len:  # Check if sgRNA is exactly 23 nt long
                            print(f"Pam was found outside designated locus_tag: {locus_tag}. To incorporate this extent borders. Skipping to next locus tag.")
                            continue  # This skips the rest of the current iteration and moves to the next feature

                        # Calculate GC content of the sgRNA
                        gc_content = len([base for base in sgrna if base in ["C", "G"]]) / len(sgrna) if len(sgrna) > 0 else 0

                        # Calculate genomic location of the sgRNA depending on the strand
                        if sequence[0] == 1:
                            genome_location = (int(location.start)) +1
                            position_sgrna = len(sequence[1]) -match.start()-3
                            strand_sgrna = sequence[0]
                        if sequence[0] == -1:
                            genome_location = int(location.start)+1
                            position_sgrna = match.end() + protospacer_len+3   
                            strand_sgrna = sequence[0]

                        # For Cas12a, extract the seed sequence from the beginning of the sgRNA
                        sgrna_seed = sgrna[:off_target_seed]

                        # Get number of off-target hits for the seed sequence
                        off_target_count = (
                            off_target_counter[sgrna_seed] - 1
                        )

                        # Store sgRNA hit
                        sgrna_hits.append(
                            (strain_name, 
                                locus_tag, 
                                genome_location, 
                                gene_strand, 
                                strand_sgrna,
                                position_sgrna,
                                gc_content,
                                pam,
                                sgrna,
                                sgrna_seed,
                                off_target_count,
                            )
                        )

    # Convert sgRNA hits into a DataFrame
    sgrna_df = pd.DataFrame(
        sgrna_hits,
        columns=['strain_name',
            "locus_tag",
            "gene_loc",
            "gene_strand",
            "sgrna_strand",
            "sgrna_loc",
            "gc",
            "pam",
            "sgrna",
            "sgrna_seed_sequence",
            "off_target_count",
        ]
    )
    return sgrna_df



def filter_guides(args: SgRNAargs, hitframe: pd.DataFrame) -> pd.DataFrame:
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
        mask = frame[column].apply(lambda x: not any(pattern in x for pattern in patterns))
        return frame[mask]

    return (hitframe
            .pipe(exclude_rows_based_on_patterns, "pam", args.pam_remove)
            .pipe(exclude_rows_based_on_patterns, "sgrna", args.sgrna_remove)
            .pipe(exclude_rows_based_on_patterns, "downstream", args.downstream_remove)
            .loc[hitframe.gc >= args.gc_lower]
            .loc[hitframe.gc <= args.gc_upper]
            .loc[hitframe.off_target_count <= args.off_target_upper])


def extract_sgRNAs(args: SgRNAargs) -> Tuple[pd.DataFrame, Counter, pd.DataFrame]:
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


    # If 'find' is in the steps, execute off-target and sgRNA finding
    if "find" in args.step:
        # Find all potential off-target hits
        off_target_counter = find_off_target_hits(sequences, args.off_target_seed, cas_type=args.cas_type)
        
        if 'cas3' in args.cas_type:
            # Find all potential sgRNA hits
            sgrna_df = find_sgrna_hits_cas3(args.dseqrecord, args.strain_name, args.locus_tag, off_target_counter, 
                                    args.off_target_seed, revcomp=revcomp)

        if 'cas9' in args.cas_type:
            # Find all potential sgRNA hits
            sgrna_df = find_sgrna_hits_cas9(args.dseqrecord, args.strain_name, args.locus_tag, off_target_counter, 
                                    args.off_target_seed, revcomp=revcomp)

        #TODO fix input to DSEQrecord
        if 'cas12a' in args.cas_type:
            # Find all potential sgRNA hits
            sgrna_df = find_sgrna_hits_cas12a(args.genbank,args.strain_name, args.locus_tag, off_target_counter, 
                                    args.off_target_seed, revcomp)

        # Sort sgrna_df by 'off-targets' in ascending order
        sgrna_df.sort_values(by='off_target_count', ascending=True, inplace=True)
    
    # Filter guides if 'filter' is in the steps
    if "filter" in args.step:
        sgrna_df = filter_guides(args, sgrna_df)
    
    # use ML models to predict efficient sgRNAs
    if "predict" in args.step and args.cas_type == 'cas9':
        print('We will try to implement this part later')

    return sgrna_df
