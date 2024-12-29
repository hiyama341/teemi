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
from Bio.Seq import Seq


def identify_base_editing_sites(sgrna_df: pd.DataFrame, editing_window_start: int = 3, editing_window_end: int = 10) -> pd.DataFrame:
    """
    Identify potential base editing sites for C-to-T substitutions in sgRNAs.
    
    Parameters
    ----------
    sgrna_df : pd.DataFrame
        DataFrame containing sgRNA information.
    editing_window_start : int, optional
        Start position of the editing window (1-based index), by default 3.
    editing_window_end : int, optional
        End position of the editing window (1-based index), by default 10.
        
    Returns
    -------
    pd.DataFrame
        Updated DataFrame with additional columns for base editing annotations.

    Examples
    --------
    >>> data = {'sgrna': ['GACCGT', 'CCGTGA']}
    >>> df = pd.DataFrame(data)
    >>> result = identify_base_editing_sites(df)
    >>> print(result)
    """
    def find_editable_cytosines(sgrna: str) -> str:
        """Find positions of editable cytosines within the editing window."""
        editable_positions = []
        for i in range(editing_window_start - 1, editing_window_end):
            if i < len(sgrna) and sgrna[i] == 'C':
                editable_positions.append(i + 1)  # Convert to 1-based index
                
        return ",".join(map(str, editable_positions))
    
    def find_context_dependent_seqs(sgrna: str) -> str:
        """Find positions of editable cytosines within the editing window."""
        sequence_context_bases = []
        for i in range(editing_window_start - 1, editing_window_end):
            if i < len(sgrna) and sgrna[i] == 'C':
                # editing context
                if i < len(sgrna) and sgrna[i-1] == 'G':
                    sequence_context_bases.append(1)
                else: 
                    sequence_context_bases.append(0)
                
        return ",".join(map(str, sequence_context_bases))
    
    sgrna_df = sgrna_df.copy() 
    sgrna_df['editing_context'] = sgrna_df['sgrna'].apply(find_context_dependent_seqs)
    sgrna_df['editable_cytosines'] = sgrna_df['sgrna'].apply(find_editable_cytosines)
    return sgrna_df


def filter_sgrnas_for_base_editing(sgrna_df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter sgRNAs to include only those with editable cytosines within the editing window.
    
    Parameters
    ----------
    sgrna_df : pd.DataFrame
        DataFrame containing sgRNA information.
        
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with sgRNAs suitable for base editing.

    Examples
    --------
    >>> data = {'sgrna': ['GACCGT', 'CCGTGA'], 'editable_cytosines': ['3', '']}
    >>> df = pd.DataFrame(data)
    >>> result = filter_sgrnas_for_base_editing(df)
    >>> print(result)
    """
    return sgrna_df[sgrna_df['editable_cytosines'] != '']


def process_base_editing(df: pd.DataFrame, gene_sequences: dict,
                         only_stop_codons: bool = False, 
                         editing_context: bool = True) -> pd.DataFrame:
    """
    Process the DataFrame to apply C-to-T (or C-to-A if sgrna_strand is -1) mutations at specified positions and identify amino acid changes.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing sgRNA information with editable cytosines.
    gene_sequences : dict
        Dictionary mapping locus tags to gene sequences.
    only_stop_codons : bool
        If True, filter out rows without STOP codon mutations.
        
        
    Returns
    -------
    pd.DataFrame
        Updated DataFrame with new columns for mutated sequences and amino acid changes.
    """
    
    def mutate_sequence(row):
        gene_seq = gene_sequences[row['locus_tag']]
        sgrna_start = row['sgrna_loc'] - 20  # Calculate the start of the sgRNA
        
        # Create a mutable list of the gene sequence
        mutated_seq = list(str(gene_seq))
        gene_length = len(mutated_seq)
        
        # Mutate the editable cytosines in the sgRNA
        for pos in map(int, row['editable_cytosines'].split(',')):
            if row['sgrna_strand'] == -1:
                genome_pos = sgrna_start - pos + 20
                if 0 <= genome_pos < gene_length:  # Boundary check
                    if mutated_seq[genome_pos] == 'G':
                        mutated_seq[genome_pos] = 'A'
                else:
                    print(f"Warning: genome_pos {genome_pos} out of range for gene length {gene_length}")
            else:
                genome_pos = sgrna_start + pos - 1
                if 0 <= genome_pos < gene_length:  # Boundary check
                    if mutated_seq[genome_pos] == 'C':
                        mutated_seq[genome_pos] = 'T'
                else:
                    print(f"Warning: genome_pos {genome_pos} out of range for gene length {gene_length}")
    
        # Return the mutated sequence as a string
        return ''.join(mutated_seq)
    
    def translate_sequence(sequence):
        """
        Translate a nucleotide sequence into an amino acid sequence.
        """
        return str(Seq(sequence).translate())

    def find_amino_acid_changes(row):
        original_seq = gene_sequences[row['locus_tag']]
        mutated_seq = row['mutated_sequence']
        
        original_aa_seq = translate_sequence(original_seq)
        mutated_aa_seq = translate_sequence(mutated_seq)
        
        mutations = []
        for i, (orig_aa, mut_aa) in enumerate(zip(original_aa_seq, mutated_aa_seq), start=1):
            if orig_aa != mut_aa:
                mutations.append(f"{orig_aa}{i}{mut_aa}")
        
        return ', '.join(mutations)
    
    def extract_first_mutation_position(mutations):
        if mutations:
            # Extract the position of the first mutation
            first_mutation = mutations.split(',')[0]
            position = int(''.join(filter(str.isdigit, first_mutation)))
            return position
        return float('inf')
    
    # Use .loc to avoid SettingWithCopyWarning
    df = df.copy()
    df.loc[:, 'mutated_sequence'] = df.apply(mutate_sequence, axis=1)
    df.loc[:, 'mutations'] = df.apply(find_amino_acid_changes, axis=1)
    
    # Remove rows where mutations column is empty
    df = df[df['mutations'] != '']
    # Drop the mutated_sequence column if you don't want it in the final output
    df = df.drop(columns=['mutated_sequence'])
    
    if only_stop_codons:
        # Filter out rows where there is no stop codon in the mutations column
        df = df[df['mutations'].str.contains(r'\*')]
        
        # Add a column for the position of the first mutation
        df['first_mutation_position'] = df['mutations'].apply(extract_first_mutation_position)
        
        # Sort the DataFrame by the position of the first mutation in ascending order
        df = df.sort_values(by='first_mutation_position', ascending=True)
        
        # Drop the first_mutation_position column if you don't want it in the final output
        df = df.drop(columns=['first_mutation_position'])

    if editing_context: 
        # Filter out rows where 'editing_context' contains '1'
        df = df[~df['editing_context'].str.contains(r'1')]

    
    return df