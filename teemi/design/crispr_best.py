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


# CRISPR-BEST is an SpCas9 system and nothing else. Its two editors are
# rAPOBEC1 (cBEST) and ecTadA (aBEST) fused to a Streptomyces codon-optimised
# SpCas9 D10A nickase, and the default 4-10 editing window below is the window
# cBEST was measured on, in that 20 nt geometry (Tong et al. 2019,
# doi:10.1073/pnas.1913493116, which gives aBEST's window only as a span of
# about 6 nt without stating positions, so 4-10 is a cBEST number serving
# both). Passing a Cas12a or Cas3 guide table through this module does not
# fail: it invents amino-acid changes, which is worse than failing. On SCO5087
# a Cas3 table yields five of them at these defaults.
#
# A guide table carries no nuclease label -- the finders' frame is eleven
# columns and ``cas_type`` is not one of them -- so the only discriminator in
# the data is the PAM. Every SpCas9 PAM ends in GG, and neither other nuclease
# in ``crispr_cas._NUCLEASES`` can produce one that does: Cas12a's pattern is
# ``TTT[ACG]``, so its third base is always T, and Cas3's is ``TTC``. Spacer
# length would not serve in its place, because ``SgRNAargs.protospacer_len``
# is overridable and a truncated sgRNA is real practice. ``cas9`` and
# ``dcas9`` are both accepted: they share the PAM and the 20 nt geometry, and
# the window is a property of the guide rather than of whether the nuclease
# still cuts.
#
# Two limits of deciding on the PAM, both deliberate. A ``pam`` cell is an
# ordinary string, so a caller who overwrites one is believed; the claim above
# is that no finder can emit a GG-ending PAM for another nuclease, not that
# the column is tamper-proof. And the PAM fixes the nuclease, not the
# geometry: a table found with ``cas_type="cas9"`` and ``protospacer_len``
# overridden carries real NGG PAMs, passes this check, and still has the 4-10
# window read off a spacer that window was never measured on. Refusing that
# would mean deciding spacer length here, which this does not do.
#
# The check is strict in ``process_base_editing``, which is where an
# amino-acid change gets asserted, and advisory in
# ``identify_base_editing_sites``, which only reports where the cytosines in a
# spacer are -- true of that spacer whichever nuclease produced it.
_CAS9_PAM_SUFFIX = "GG"


def _reject_non_cas9_guides(sgrna_df: pd.DataFrame) -> None:
    """Raise unless every guide in ``sgrna_df`` came from Cas9 or dCas9.

    Decided on the ``pam`` column, for the reasons given above this function.

    A table with no ``pam`` column is not judged here, because whether that is
    an error depends on the caller. ``process_base_editing`` lists ``pam`` among
    the columns it requires, since that is where a biological claim gets made;
    ``identify_base_editing_sites`` only annotates sequence, which is true of a
    spacer whatever cut it, so it checks the PAMs when they are present and does
    not demand them.

    An empty table is allowed through. A gene whose guides carry no editable
    cytosine is a legitimate result, and a table the finders returned empty
    carries no PAM to judge, so there is nothing to decide on. That is the one
    case this cannot catch.
    """
    if "pam" not in sgrna_df.columns:
        return
    # ``fillna`` so a missing PAM is judged as ``''``: left as NaN it reaches
    # ``sorted`` below as a float, which cannot be ordered against the real PAM
    # strings, so a table carrying one foreign PAM and one missing raised
    # ``TypeError`` instead of the message below. ``strip`` because a CSV round
    # trip can pad the column, and only the suffix is tested, so a trailing
    # space was refused with a message that renders identically to an accepted
    # PAM while a leading one passed. ``upper`` so a record read in lower case
    # is not refused for that alone.
    pams = sgrna_df["pam"].fillna("").astype(str).str.strip().str.upper()
    foreign = sorted(set(pams[~pams.str.endswith(_CAS9_PAM_SUFFIX)]))
    if foreign:
        raise ValueError(
            f"PAM(s) {foreign} are not supported for CRISPR-BEST; only 'cas9' "
            "guides can be base edited here, and every Cas9 PAM ends in "
            f"{_CAS9_PAM_SUFFIX!r}"
        )


def identify_base_editing_sites(sgrna_df: pd.DataFrame, editing_window_start: int = 4, editing_window_end: int = 10) -> pd.DataFrame:
    """
    Identify potential base editing sites for C-to-T substitutions in sgRNAs.

    CRISPR-BEST is an SpCas9 base editor, so a guide table carrying a ``pam``
    column from any other nuclease is refused -- see
    ``_reject_non_cas9_guides``. The column is not required here, only checked
    when it is there; ``process_base_editing`` requires it.

    Parameters
    ----------
    sgrna_df : pd.DataFrame
        DataFrame containing sgRNA information. Needs the ``sgrna`` column.
    editing_window_start : int, optional
        Start position of the editing window (1-based index), by default 4.
        Positions are counted from the PAM-distal end of the protospacer, so
        the default window 4-10 is the one reported for CRISPR-cBEST; position
        3 is outside it and has no measured conversion.
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
                # editing context: a G immediately 5' of the C. Guard i > 0 so
                # a window starting at position 1 cannot read sgrna[-1], the
                # base at the far end of the guide.
                if i > 0 and sgrna[i-1] == 'G':
                    sequence_context_bases.append(1)
                else: 
                    sequence_context_bases.append(0)
                
        return ",".join(map(str, sequence_context_bases))
    
    _reject_non_cas9_guides(sgrna_df)

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
    Process the DataFrame to apply C-to-T mutations at specified positions and identify amino acid changes.

    The deaminated cytosine is a C in the gene sequence when the guide reads
    along the CDS, and a G there when it reads against it, so the edit shows up
    as C-to-T or G-to-A respectively.

    CRISPR-BEST is an SpCas9 base editor, so this refuses a guide table from any
    other nuclease -- see ``_reject_non_cas9_guides``.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing sgRNA information with editable cytosines. Needs
        the ``locus_tag``, ``sgrna``, ``sgrna_loc``, ``sgrna_strand``,
        ``gene_strand``, ``editable_cytosines`` and ``pam`` columns, plus
        ``editing_context`` when ``editing_context`` is left True.
    gene_sequences : dict
        Dictionary mapping locus tags to gene sequences.
    only_stop_codons : bool
        If True, filter out rows without STOP codon mutations.
        
        
    Returns
    -------
    pd.DataFrame
        The input columns plus ``mutations``, with the rows that produce no
        change dropped. ``mutated_sequence`` is an intermediate and is not
        returned. An empty input comes back empty, with an empty ``mutations``.
    """
    
    def mutate_sequence(row):
        gene_seq = str(gene_sequences[row['locus_tag']])
        sgrna = row['sgrna']

        # ``sgrna_strand`` and ``gene_strand`` are both genome-relative, so the
        # guide reads along the CDS exactly when the two agree. Comparing them
        # rather than testing sgrna_strand alone is what keeps minus-strand
        # genes right.
        sense_to_cds = row['sgrna_strand'] == row['gene_strand']
        protospacer = sgrna if sense_to_cds else str(Seq(sgrna).reverse_complement())

        # ``sgrna_loc`` minus the protospacer length is the protospacer's offset
        # in the CDS, for both of the finder's scan branches and both gene
        # strands. Searching for the protospacer instead would read the wrong
        # copy inside a tandem repeat, where the same 20-mer occurs at several
        # overlapping offsets and only this coordinate says which one the guide
        # came from. Checking the sequence at the offset keeps that honest: a
        # guide table that does not belong to this gene fails loudly rather
        # than editing an arbitrary position.
        start_idx = row['sgrna_loc'] - len(sgrna)
        if start_idx < 0 or gene_seq[start_idx : start_idx + len(sgrna)] != protospacer:
            raise RuntimeError(
                f"protospacer {protospacer} is not at offset {start_idx} of "
                f"{row['locus_tag']}: the guide table and gene_sequences do not "
                "describe the same gene"
            )

        mutated_seq = list(gene_seq)
        for pos in map(int, row['editable_cytosines'].split(',')):
            if sense_to_cds:
                genome_pos = start_idx + pos - 1
            else:
                # The protospacer was reverse-complemented to match gene_seq, so
                # position pos from the guide's 5' end lands len(sgrna) - pos
                # along it.
                genome_pos = start_idx + len(sgrna) - pos
            # The protospacer was checked against gene_seq above, so a genuine
            # editable cytosine is a C here when the guide reads along the CDS
            # and a G when it reads against it. Anything else means the position
            # did not come from identify_base_editing_sites.
            base = mutated_seq[genome_pos]
            if base not in ('C', 'G'):
                raise RuntimeError(
                    f"position {pos} of {sgrna} is {base!r} in "
                    f"{row['locus_tag']}, not an editable cytosine"
                )
            mutated_seq[genome_pos] = 'T' if base == 'C' else 'A'

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
        # Rows without mutations have already been dropped below
        first_mutation = mutations.split(',')[0]
        return int(''.join(filter(str.isdigit, first_mutation)))
    
    required = {
        "locus_tag", "sgrna", "sgrna_loc", "sgrna_strand", "gene_strand",
        "editable_cytosines",
        # Required because this is where a biological claim is made, and ``pam``
        # is what says the claim is about a nuclease CRISPR-BEST exists for.
        "pam",
    }
    if editing_context:
        # Read at the end of this function, so required only when asked for.
        # Without this a frame meeting the documented contract raised KeyError
        # from pandas rather than the message below.
        required.add("editing_context")
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"process_base_editing needs the column(s) {sorted(missing)}; pass a "
            "guide table from the sgRNA finders through "
            "identify_base_editing_sites and filter_sgrnas_for_base_editing"
        )
    _reject_non_cas9_guides(df)

    # ``identify_base_editing_sites`` leaves this empty for a guide with no
    # cytosine in the window -- 2.9% of the guides over 50 CDSs of the
    # validation chromosome, and some on every one of them -- and
    # ``filter_sgrnas_for_base_editing`` is what drops them. Skipping it
    # reached ``int('')`` in the apply below and died with a bare pandas
    # ``ValueError: invalid literal for int() with base 10: ''``, which names
    # neither the column nor the step that was left out.
    unfiltered = df["editable_cytosines"].isna() | (
        df["editable_cytosines"].astype(str).str.strip() == ""
    )
    if unfiltered.any():
        raise ValueError(
            f"{int(unfiltered.sum())} of {len(df)} guides have no editable "
            "cytosine; pass the table through filter_sgrnas_for_base_editing "
            "before process_base_editing"
        )

    # Use .loc to avoid SettingWithCopyWarning
    df = df.copy()
    if df.empty:
        # A row-wise apply over a 0-row frame returns a DataFrame rather than a
        # Series, and assigning that to a single column raises ValueError from
        # pandas. A gene whose guides carry no editable cytosine is a legitimate
        # result, so return what the non-empty path returns: the input columns
        # plus ``mutations``, without the intermediate ``mutated_sequence``,
        # which the caller may have passed in.
        return df.drop(columns=["mutated_sequence"], errors="ignore").assign(
            mutations=pd.Series(dtype=object)
        )
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