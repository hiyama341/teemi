#!/usr/bin/env python

import pandas as pd

from teemi.design.crispr_best import (
    filter_sgrnas_for_base_editing,
    identify_base_editing_sites,
    process_base_editing,
)

# Genes used by the process_base_editing tests. Index 0 is the first base.
#                  0  1  2   3  4  5   6  7  8   9 10 11  12..
GENE_STOP_PLUS = "ATG" "AAA" "AAA" "CAA" "GGG" "GGG" "GGG" "GGG" "GGG"  # MKKQGGGGG
GENE_STOP_MINUS = "ATG" "AAA" "AAA" "TGG" "CCC" "CCC" "CCC"  # MKKWPPP
GENE_MISSENSE_AND_STOP = "ATG" "CCA" "AAA" "CAG" "GGG" "GGG"  # MPKQGG
GENE_EARLY_STOP = "ATG" "CAA" "GGG" "GGG" "GGG" "GGG"  # MQGGGG
GENE_DOUBLE_STOP = "ATG" "CAA" "CAA" "GGG" "GGG" "GGG"  # MQQGGG


def _row(locus_tag, sgrna_loc, sgrna_strand, editable_cytosines, editing_context="0"):
    return {
        "locus_tag": locus_tag,
        "sgrna_loc": sgrna_loc,
        "sgrna_strand": sgrna_strand,
        "editable_cytosines": editable_cytosines,
        "editing_context": editing_context,
    }


# ---------------------------------------------------------------------------
# identify_base_editing_sites
# ---------------------------------------------------------------------------


def test_identify_base_editing_sites_default_window():
    # Window is positions 3..10 (1-based), i.e. indices 2..9.
    #              idx 0123456789
    df = pd.DataFrame({"sgrna": ["AAGCAAAACA"]})

    result = identify_base_editing_sites(df)

    # C at index 3 -> position 4, C at index 8 -> position 9.
    assert result.loc[0, "editable_cytosines"] == "4,9"
    # The C at index 3 is preceded by G (index 2) -> context 1;
    # the C at index 8 is preceded by A (index 7) -> context 0.
    assert result.loc[0, "editing_context"] == "1,0"


def test_identify_base_editing_sites_ignores_cytosines_outside_window():
    # C at index 0 (position 1) and index 11 (position 12) are both outside
    # the default 3..10 window.
    df = pd.DataFrame({"sgrna": ["CAAAAAAAAAAC", "AAAAAAAAAA"]})

    result = identify_base_editing_sites(df)

    assert list(result["editable_cytosines"]) == ["", ""]
    assert list(result["editing_context"]) == ["", ""]


def test_identify_base_editing_sites_handles_sgrna_shorter_than_window():
    # Only 4 bases, so indices 4..9 of the default window do not exist.
    df = pd.DataFrame({"sgrna": ["AACA"]})

    result = identify_base_editing_sites(df)

    assert result.loc[0, "editable_cytosines"] == "3"
    assert result.loc[0, "editing_context"] == "0"


def test_identify_base_editing_sites_custom_window():
    #              idx 0123456
    df = pd.DataFrame({"sgrna": ["ACGCACG"]})

    narrow = identify_base_editing_sites(df, editing_window_start=2, editing_window_end=3)
    wide = identify_base_editing_sites(df, editing_window_start=2, editing_window_end=5)

    # Narrow window covers indices 1..2 -> only the C at index 1 (position 2).
    assert narrow.loc[0, "editable_cytosines"] == "2"
    # Wide window covers indices 1..4 -> C at index 1 and C at index 3.
    assert wide.loc[0, "editable_cytosines"] == "2,4"
    # Index 1 is preceded by A -> 0; index 3 is preceded by G -> 1.
    assert wide.loc[0, "editing_context"] == "0,1"


def test_identify_base_editing_sites_does_not_mutate_input():
    df = pd.DataFrame({"sgrna": ["AAGCAAAACA"]})

    identify_base_editing_sites(df)

    assert list(df.columns) == ["sgrna"]


# ---------------------------------------------------------------------------
# filter_sgrnas_for_base_editing
# ---------------------------------------------------------------------------


def test_filter_sgrnas_for_base_editing_drops_rows_without_editable_cytosines():
    df = pd.DataFrame(
        {"sgrna": ["GACCGT", "CCGTGA", "AAAAAA"], "editable_cytosines": ["3", "", "4,5"]}
    )

    result = filter_sgrnas_for_base_editing(df)

    assert list(result.index) == [0, 2]
    assert list(result["editable_cytosines"]) == ["3", "4,5"]


# ---------------------------------------------------------------------------
# process_base_editing
# ---------------------------------------------------------------------------


def test_process_base_editing_plus_strand_c_to_t():
    # sgrna_loc 26 -> sgrna_start 6; editable cytosine 4 -> gene index 9,
    # the C of codon 4 ("CAA" = Q), which becomes "TAA" = stop.
    df = pd.DataFrame([_row("g1", 26, 1, "4")])

    result = process_base_editing(df, {"g1": GENE_STOP_PLUS})

    assert len(result) == 1
    assert result.iloc[0]["mutations"] == "Q4*"
    # The helper column is not leaked into the output.
    assert "mutated_sequence" not in result.columns


def test_process_base_editing_minus_strand_g_to_a():
    # sgrna_strand -1 -> genome_pos = sgrna_loc - pos = 15 - 4 = 11, the last G
    # of codon 4 ("TGG" = W), which becomes "TGA" = stop.
    df = pd.DataFrame([_row("g2", 15, -1, "4")])

    result = process_base_editing(df, {"g2": GENE_STOP_MINUS})

    assert len(result) == 1
    assert result.iloc[0]["mutations"] == "W4*"


def test_process_base_editing_drops_rows_without_amino_acid_change():
    # Gene index 9 is a G on the plus strand, so no C-to-T edit happens and the
    # protein is unchanged.
    df = pd.DataFrame([_row("g1", 26, 1, "4"), _row("early", 20, 1, "7")])

    result = process_base_editing(df, {"g1": GENE_STOP_PLUS, "early": GENE_EARLY_STOP})

    assert list(result["mutations"]) == ["Q4*"]


def test_process_base_editing_reports_multiple_mutations():
    # sgrna_start 0; cytosines at positions 4 and 7 -> gene indices 3 and 6,
    # the C of codon 2 and the C of codon 3, both "CAA" = Q -> "TAA" = stop.
    df = pd.DataFrame([_row("g5", 20, 1, "4,7")])

    result = process_base_editing(df, {"g5": GENE_DOUBLE_STOP})

    assert result.iloc[0]["mutations"] == "Q2*, Q3*"


def test_process_base_editing_warns_when_plus_strand_position_out_of_range(capsys):
    # sgrna_loc 1 -> sgrna_start -19 -> genome_pos -16, before the gene start.
    df = pd.DataFrame([_row("g1", 1, 1, "4")])

    result = process_base_editing(df, {"g1": GENE_STOP_PLUS})

    assert result.empty  # nothing mutated, so no amino-acid change
    out = capsys.readouterr().out
    assert "genome_pos -16 out of range for gene length 27" in out


def test_process_base_editing_warns_when_minus_strand_position_out_of_range(capsys):
    # sgrna_strand -1 -> genome_pos = 100 - 4 = 96, past the 21 nt gene.
    df = pd.DataFrame([_row("g2", 100, -1, "4")])

    result = process_base_editing(df, {"g2": GENE_STOP_MINUS})

    assert result.empty
    out = capsys.readouterr().out
    assert "genome_pos 96 out of range for gene length 21" in out


def test_process_base_editing_only_stop_codons_filters_and_sorts():
    df = pd.DataFrame(
        [
            # "CCA" (P) codon 2 -> "TCA" (S): a missense change, no stop.
            _row("missense", 20, 1, "4"),
            # "CAG" (Q) codon 4 -> "TAG": stop at amino acid 4.
            _row("missense", 26, 1, "4"),
            # "CAA" (Q) codon 2 -> "TAA": stop at amino acid 2.
            _row("early", 20, 1, "4"),
        ]
    )
    gene_sequences = {
        "missense": GENE_MISSENSE_AND_STOP,
        "early": GENE_EARLY_STOP,
    }

    everything = process_base_editing(df, gene_sequences)
    assert list(everything["mutations"]) == ["P2S", "Q4*", "Q2*"]

    stops_only = process_base_editing(df, gene_sequences, only_stop_codons=True)

    # The missense row is gone and the remaining rows are sorted by the
    # position of their first mutation.
    assert list(stops_only["mutations"]) == ["Q2*", "Q4*"]
    assert "first_mutation_position" not in stops_only.columns


def test_process_base_editing_editing_context_filter():
    rows = pd.DataFrame(
        [
            _row("g1", 26, 1, "4", editing_context="0"),
            _row("g5", 20, 1, "4,7", editing_context="0,1"),
        ]
    )
    gene_sequences = {"g1": GENE_STOP_PLUS, "g5": GENE_DOUBLE_STOP}

    filtered = process_base_editing(rows, gene_sequences, editing_context=True)
    unfiltered = process_base_editing(rows, gene_sequences, editing_context=False)

    # A "1" in editing_context marks a G-preceded cytosine, which is discarded.
    assert list(filtered["mutations"]) == ["Q4*"]
    assert list(unfiltered["mutations"]) == ["Q4*", "Q2*, Q3*"]


def test_process_base_editing_does_not_mutate_input():
    df = pd.DataFrame([_row("g1", 26, 1, "4")])

    process_base_editing(df, {"g1": GENE_STOP_PLUS})

    assert "mutations" not in df.columns
