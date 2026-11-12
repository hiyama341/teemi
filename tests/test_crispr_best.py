#!/usr/bin/env python

import pandas as pd
import pytest
from Bio.Seq import Seq

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
GENE_SYNONYMOUS = "ATG" "GGG" "GGG" "GGG"  # MGGG; GGG -> GGA stays Gly

GUIDE_LEN = 12


def _row(
    locus_tag,
    gene,
    target_idx,
    *,
    sense,
    gene_strand=1,
    editing_context="0",
):
    """Build a guide row that places an editable cytosine on ``gene[target_idx]``.

    The guide is cut out of the gene rather than written by hand, so the row is
    always self-consistent: a guide reading along the CDS carries a C where the
    gene has one, and a guide reading against it carries a C where the gene has
    the complementary G. ``sense`` says which, and ``gene_strand`` then fixes
    ``sgrna_strand``, since the code decides orientation by comparing the two.
    """
    start = max(0, min(target_idx - GUIDE_LEN // 2, len(gene) - GUIDE_LEN))
    region = gene[start : start + GUIDE_LEN]
    if sense:
        sgrna = region
        position = target_idx - start + 1
    else:
        sgrna = str(Seq(region).reverse_complement())
        position = GUIDE_LEN - (target_idx - start)
    assert sgrna[position - 1] == "C", (
        f"gene[{target_idx}] is {gene[target_idx]!r}; a "
        f"{'sense' if sense else 'antisense'} guide needs a "
        f"{'C' if sense else 'G'} there"
    )
    return {
        "locus_tag": locus_tag,
        "sgrna": sgrna,
        "gene_strand": gene_strand,
        "sgrna_strand": gene_strand if sense else -gene_strand,
        # What the finder reports: the offset one past the protospacer's end.
        "sgrna_loc": start + GUIDE_LEN,
        "editable_cytosines": str(position),
        "editing_context": editing_context,
    }


# ---------------------------------------------------------------------------
# identify_base_editing_sites
# ---------------------------------------------------------------------------


def test_identify_base_editing_sites_default_window():
    # Window is positions 4..10 (1-based), i.e. indices 3..9.
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
    # the default 4..10 window.
    df = pd.DataFrame({"sgrna": ["CAAAAAAAAAAC", "AAAAAAAAAA"]})

    result = identify_base_editing_sites(df)

    assert list(result["editable_cytosines"]) == ["", ""]
    assert list(result["editing_context"]) == ["", ""]


def test_identify_base_editing_sites_handles_sgrna_shorter_than_window():
    # Only 4 bases, so indices 4..9 of the default window do not exist. The C at
    # index 2 is position 3, which the CRISPR-cBEST window 4..10 excludes.
    df = pd.DataFrame({"sgrna": ["AACA"]})

    result = identify_base_editing_sites(df)

    assert result.loc[0, "editable_cytosines"] == ""
    assert result.loc[0, "editing_context"] == ""


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


def test_process_base_editing_guide_along_the_cds_edits_c_to_t():
    # Gene index 9 is the C of codon 4 ("CAA" = Q), which becomes "TAA" = stop.
    df = pd.DataFrame([_row("g1", GENE_STOP_PLUS, 9, sense=True)])

    result = process_base_editing(df, {"g1": GENE_STOP_PLUS})

    assert len(result) == 1
    assert result.iloc[0]["mutations"] == "Q4*"
    # The helper column is not leaked into the output.
    assert "mutated_sequence" not in result.columns


def test_process_base_editing_guide_against_the_cds_edits_g_to_a():
    # The guide reads against the CDS, so its C is the G at gene index 11, the
    # last base of codon 4 ("TGG" = W), which becomes "TGA" = stop.
    df = pd.DataFrame([_row("g2", GENE_STOP_MINUS, 11, sense=False)])

    result = process_base_editing(df, {"g2": GENE_STOP_MINUS})

    assert len(result) == 1
    assert result.iloc[0]["mutations"] == "W4*"


@pytest.mark.parametrize("gene_strand", [1, -1])
@pytest.mark.parametrize("sense", [True, False])
def test_process_base_editing_depends_on_strand_only_through_the_comparison(
    gene_strand, sense
):
    """A gene on the minus strand must give the same answer as one on the plus.

    ``gene_sequences`` is always the CDS read 5'->3', so what decides the edit
    is whether the guide runs along that sequence -- ``sgrna_strand ==
    gene_strand`` -- and not the sign of ``sgrna_strand`` on its own. Reading
    the sign alone silently mirrors the window on every minus-strand gene.
    """
    gene = GENE_STOP_PLUS if sense else GENE_STOP_MINUS
    target = 9 if sense else 11
    expected = "Q4*" if sense else "W4*"

    df = pd.DataFrame([_row("g", gene, target, sense=sense, gene_strand=gene_strand)])

    result = process_base_editing(df, {"g": gene})

    assert list(result["mutations"]) == [expected]


def test_process_base_editing_drops_rows_without_amino_acid_change():
    # Gene index 5 is the third base of codon 2 ("GGG" = G); G-to-A gives "GGA",
    # still Gly, so the row carries no amino-acid change and is dropped.
    df = pd.DataFrame(
        [
            _row("g1", GENE_STOP_PLUS, 9, sense=True),
            _row("syn", GENE_SYNONYMOUS, 5, sense=False),
        ]
    )

    result = process_base_editing(
        df, {"g1": GENE_STOP_PLUS, "syn": GENE_SYNONYMOUS}
    )

    assert list(result["mutations"]) == ["Q4*"]


def test_process_base_editing_reports_multiple_mutations():
    # One guide over gene[0:12] carries the C of codon 2 (index 3) at position 4
    # and the C of codon 3 (index 6) at position 7. Both "CAA" (Q) -> "TAA".
    row = _row("g5", GENE_DOUBLE_STOP, 3, sense=True)
    row["editable_cytosines"] = "4,7"
    df = pd.DataFrame([row])

    result = process_base_editing(df, {"g5": GENE_DOUBLE_STOP})

    assert result.iloc[0]["mutations"] == "Q2*, Q3*"


def test_process_base_editing_raises_when_the_guide_is_not_where_it_says():
    # The sequence at sgrna_loc is checked, so a guide table that does not
    # belong to this gene fails loudly instead of editing an arbitrary position.
    row = _row("g1", GENE_STOP_PLUS, 9, sense=True)
    row["sgrna"] = "CCCCCCCCCCCC"
    df = pd.DataFrame([row])

    with pytest.raises(RuntimeError, match="is not at offset"):
        process_base_editing(df, {"g1": GENE_STOP_PLUS})


def test_process_base_editing_raises_when_a_position_is_not_a_cytosine():
    # Gene index 10 is an A, so position 8 of this guide is not editable. The old
    # code silently edited nothing and reported the guide as having no effect.
    row = _row("g1", GENE_STOP_PLUS, 9, sense=True)
    row["editable_cytosines"] = "8"
    df = pd.DataFrame([row])

    with pytest.raises(RuntimeError, match="not an editable cytosine"):
        process_base_editing(df, {"g1": GENE_STOP_PLUS})


def test_process_base_editing_requires_the_columns_it_reads():
    df = pd.DataFrame([{"locus_tag": "g1", "editable_cytosines": "4"}])

    with pytest.raises(ValueError, match="gene_strand"):
        process_base_editing(df, {"g1": GENE_STOP_PLUS})


def test_process_base_editing_only_stop_codons_filters_and_sorts():
    df = pd.DataFrame(
        [
            # "CCA" (P) codon 2 -> "TCA" (S): a missense change, no stop.
            _row("missense", GENE_MISSENSE_AND_STOP, 3, sense=True),
            # "CAG" (Q) codon 4 -> "TAG": stop at amino acid 4.
            _row("missense", GENE_MISSENSE_AND_STOP, 9, sense=True),
            # "CAA" (Q) codon 2 -> "TAA": stop at amino acid 2.
            _row("early", GENE_EARLY_STOP, 3, sense=True),
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
    multi = _row("g5", GENE_DOUBLE_STOP, 3, sense=True, editing_context="0,1")
    multi["editable_cytosines"] = "4,7"
    rows = pd.DataFrame(
        [_row("g1", GENE_STOP_PLUS, 9, sense=True, editing_context="0"), multi]
    )
    gene_sequences = {"g1": GENE_STOP_PLUS, "g5": GENE_DOUBLE_STOP}

    filtered = process_base_editing(rows, gene_sequences, editing_context=True)
    unfiltered = process_base_editing(rows, gene_sequences, editing_context=False)

    # A "1" in editing_context marks a G-preceded cytosine, which is discarded.
    assert list(filtered["mutations"]) == ["Q4*"]
    assert list(unfiltered["mutations"]) == ["Q4*", "Q2*, Q3*"]


def test_process_base_editing_does_not_mutate_input():
    df = pd.DataFrame([_row("g1", GENE_STOP_PLUS, 9, sense=True)])

    process_base_editing(df, {"g1": GENE_STOP_PLUS})

    assert "mutations" not in df.columns


# A tandem repeat long enough to contain the same protospacer twice. Editing
# position 4 of that guide hits codon 3 at the first copy and codon 7 at the
# second, so the two guides have different answers.
GENE_TANDEM_REPEAT = "ATG" "GGGCAAGGGGGG" "GGGCAAGGGGGG" "TAA"  # MGQGGGQGG*


def test_process_base_editing_distinguishes_repeated_protospacers():
    """Placement must use sgrna_loc, not a search for the protospacer.

    Inside a tandem repeat the same protospacer sits at several offsets, and
    ``gene_seq.find()`` would return the first one for both guides. On
    *S. coelicolor* that silently moves the edit for 15 of 133,168 guides,
    always onto a position the guide cannot even bind, because the repeat
    breaks inside the PAM.
    """
    guide = GENE_TANDEM_REPEAT[3:15]
    assert GENE_TANDEM_REPEAT.count(guide) == 2

    # The two rows differ only in sgrna_loc, which is the one field that says
    # which copy the guide came from. The protospacers start at 3 and 15, and
    # sgrna_loc is one past their ends.
    rows = pd.DataFrame(
        [
            {
                "locus_tag": "tandem",
                "sgrna": guide,
                "gene_strand": 1,
                "sgrna_strand": 1,
                "sgrna_loc": start + len(guide),
                "editable_cytosines": "4",
                "editing_context": "0",
            }
            for start in (3, 15)
        ]
    )

    result = process_base_editing(rows, {"tandem": GENE_TANDEM_REPEAT})

    assert list(result["mutations"]) == ["Q3*", "Q7*"]
