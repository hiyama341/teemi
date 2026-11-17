#!/usr/bin/env python

import itertools
import re

import pandas as pd
import pytest
from Bio.Seq import Seq

from teemi.design.crispr_best import (
    filter_sgrnas_for_base_editing,
    identify_base_editing_sites,
    process_base_editing,
)
from teemi.design.crispr_cas import _NUCLEASES

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
        # process_base_editing refuses a table whose PAM is not Cas9's, so the
        # rows it is given carry one. Which NGG it is makes no difference to
        # any assertion here.
        "pam": "AGG",
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
                "pam": "AGG",
            }
            for start in (3, 15)
        ]
    )

    result = process_base_editing(rows, {"tandem": GENE_TANDEM_REPEAT})

    assert list(result["mutations"]) == ["Q3*", "Q7*"]


# ---------------------------------------------------------------------------
# CRISPR-BEST is dCas9 only
#
# The editors are rAPOBEC1 and ecTadA fused to an SpCas9 D10A nickase, and the
# 4-10 editing window is a measurement on that 20 nt Cas9 geometry, so a guide
# table from another nuclease has to be refused. It used not to be: with the
# sgrna_loc positions correct, a Cas3 table passed every check and came back
# with invented amino-acid changes, five of them on SCO5087.
# ---------------------------------------------------------------------------


def _cas9_row(pam="AGG"):
    """One well-formed Cas9 row whose only interesting field is the PAM.

    Index 9 of ``GENE_STOP_PLUS`` is the C the other tests edit, so the row
    produces a real mutation and an acceptance can be told from a no-op.
    """
    row = _row("g", GENE_STOP_PLUS, 9, sense=True)
    row["pam"] = pam
    return row


def _pams_matching(pattern, length):
    """Every k-mer of ``length`` that ``pattern`` matches end to end."""
    return [
        "".join(kmer)
        for kmer in itertools.product("ACGT", repeat=length)
        if re.fullmatch(pattern, "".join(kmer))
    ]


def test_no_pam_another_supported_nuclease_can_emit_is_accepted():
    """Tied to the nuclease table, not to hand-written strings.

    This is the assertion that has to keep holding if a nuclease is ever added:
    the guard tests ``endswith("GG")``, and that is only a sound proxy for "came
    from Cas9" as long as no other entry in ``_NUCLEASES`` can produce such a
    PAM. Expanding each pattern over ACGT checks exactly that, so adding an
    entry whose PAM could end in GG fails here rather than silently widening
    what CRISPR-BEST accepts.
    """
    foreign = {
        name: _pams_matching(spec.pam_pattern, spec.pam_len)
        for name, spec in _NUCLEASES.items()
        if name not in ("cas9", "dcas9")
    }
    # Guard the guard: if a pattern expands to nothing the loop below is vacuous.
    assert foreign, "no non-Cas9 nucleases found in _NUCLEASES"
    for name, pams in foreign.items():
        assert pams, f"{name}'s pam_pattern expanded to no k-mers"
        for pam in pams:
            with pytest.raises(ValueError, match="not supported for CRISPR-BEST"):
                process_base_editing(
                    pd.DataFrame([_cas9_row(pam)]), {"g": GENE_STOP_PLUS}
                )


@pytest.mark.parametrize("pam", ["AGG", "CGG", "GGG", "TGG", "agg", "tgg"])
def test_every_cas9_pam_is_accepted(pam):
    """All four NGG, and lower case, since a record can be read either way.

    ``_NUCLEASES["cas9"].pam_pattern`` is ``(?=CC)``, the form matched on the
    scanned strand; what the finder reports in the ``pam`` column is the NGG
    form, which is what this covers.
    """
    result = process_base_editing(
        pd.DataFrame([_cas9_row(pam)]), {"g": GENE_STOP_PLUS}
    )

    assert len(result) == 1


@pytest.mark.parametrize("pam", ["GGA", "GGT", "GGTTC", "AGGAG"])
def test_a_pam_containing_gg_without_ending_in_it_is_refused(pam):
    """Pins ``endswith`` rather than ``contains``.

    Every other PAM in these tests either ends in GG or holds no GG at all, so
    relaxing the guard to ``str.contains`` passed the whole suite. No entry in
    ``_NUCLEASES`` emits a PAM like these, so this closes no reachable bypass;
    it pins the predicate that ``_CAS9_PAM_SUFFIX``, the module comment and the
    error message all state.
    """
    with pytest.raises(ValueError, match="not supported for CRISPR-BEST"):
        process_base_editing(
            pd.DataFrame([_cas9_row(pam)]), {"g": GENE_STOP_PLUS}
        )


@pytest.mark.parametrize("pam", ["AGG ", " AGG", "\tagg\n"])
def test_a_cas9_pam_with_surrounding_whitespace_is_accepted(pam):
    """A CSV round trip can pad the column; the guide is still Cas9.

    Without the ``strip`` this was asymmetric, because only the suffix is
    tested: a leading space passed and a trailing one was refused, reported as
    a PAM that renders identically to an accepted one.
    """
    result = process_base_editing(
        pd.DataFrame([_cas9_row(pam)]), {"g": GENE_STOP_PLUS}
    )

    assert len(result) == 1


def test_process_base_editing_refuses_a_cas3_guide_table():
    result = pytest.raises(
        ValueError,
        match=r"PAM\(s\) \['TTC'\] are not supported for CRISPR-BEST",
    )
    with result:
        process_base_editing(
            pd.DataFrame([_cas9_row("TTC")]), {"g": GENE_STOP_PLUS}
        )


def test_one_foreign_row_refuses_the_whole_table_and_the_message_names_it():
    """A mixed table is a scope error too, and only the odd PAM is reported."""
    rows = pd.DataFrame([_cas9_row("AGG"), _cas9_row("TTC"), _cas9_row("CGG")])

    with pytest.raises(ValueError, match=r"\['TTC'\]") as excinfo:
        process_base_editing(rows, {"g": GENE_STOP_PLUS})

    assert "AGG" not in str(excinfo.value)
    assert "CGG" not in str(excinfo.value)


def test_process_base_editing_requires_the_pam_column():
    row = _cas9_row()
    del row["pam"]

    with pytest.raises(ValueError, match=r"needs the column\(s\) \['pam'\]"):
        process_base_editing(pd.DataFrame([row]), {"g": GENE_STOP_PLUS})


@pytest.mark.parametrize("missing", [None, float("nan"), ""])
def test_a_missing_pam_is_refused_rather_than_passed_over(missing):
    """A guide with no recorded PAM carries no evidence it came from Cas9.

    Refused deliberately, not by accident. What the ``fillna`` in the guard is
    for is pinned by
    ``test_a_missing_pam_beside_a_foreign_one_reports_both`` instead.
    """
    with pytest.raises(ValueError, match="not supported for CRISPR-BEST"):
        process_base_editing(
            pd.DataFrame([_cas9_row(missing)]), {"g": GENE_STOP_PLUS}
        )


def test_a_missing_pam_among_real_ones_still_refuses():
    """A NaN mixed with valid PAMs must not slip through."""
    rows = pd.DataFrame([_cas9_row("AGG"), _cas9_row(None), _cas9_row("CGG")])

    with pytest.raises(ValueError, match="not supported for CRISPR-BEST"):
        process_base_editing(rows, {"g": GENE_STOP_PLUS})


def test_a_missing_pam_beside_a_foreign_one_reports_both():
    """What the ``fillna`` in the guard is actually for.

    Left as NaN, a missing PAM reaches ``sorted`` as a float, which cannot be
    ordered against the real PAM strings, so this table raised ``TypeError: '<'
    not supported between instances of 'str' and 'float'`` out of the guard
    rather than the refusal below.
    """
    rows = pd.DataFrame([_cas9_row("TTC"), _cas9_row(None)])

    with pytest.raises(
        ValueError, match="not supported for CRISPR-BEST"
    ) as excinfo:
        process_base_editing(rows, {"g": GENE_STOP_PLUS})

    assert "TTC" in str(excinfo.value)
    assert "''" in str(excinfo.value)


def test_identify_base_editing_sites_refuses_a_non_cas9_table_when_it_can_tell():
    df = pd.DataFrame({"sgrna": ["AAGCAAAACA"], "pam": ["TTC"]})

    with pytest.raises(ValueError, match="not supported for CRISPR-BEST"):
        identify_base_editing_sites(df)


def test_identify_base_editing_sites_does_not_require_a_pam_column():
    """It only reports where the cytosines in a spacer are.

    That is true of the spacer whichever nuclease produced it, so the column is
    checked when present and not demanded. ``process_base_editing``, which
    asserts an amino-acid change, requires it.
    """
    result = identify_base_editing_sites(pd.DataFrame({"sgrna": ["AAGCAAAACA"]}))

    assert result.loc[0, "editable_cytosines"] == "4,9"


def test_process_base_editing_returns_an_empty_frame_rather_than_raising():
    """A gene whose guides carry no editable cytosine is a legitimate result.

    A row-wise apply over a 0-row frame returns a DataFrame rather than a
    Series, so assigning it to one column used to raise a bare pandas
    ValueError naming nothing a caller could act on.

    The input carries ``mutated_sequence`` because a caller re-running the
    chain can already have it, and the non-empty path drops that column after
    using it; without it in the input the assertion below is vacuous.
    """
    empty = pd.DataFrame(columns=list(_cas9_row()) + ["mutated_sequence"])

    result = process_base_editing(empty, {"g": GENE_STOP_PLUS})

    assert result.empty
    # The same column set the non-empty path returns.
    assert "mutations" in result.columns
    assert "mutated_sequence" not in result.columns
    assert set(_cas9_row()) <= set(result.columns)


def test_process_base_editing_requires_editing_context_only_when_it_reads_it():
    """It is read at the end, and only when the flag is left True.

    Listing it unconditionally would refuse a legitimate
    ``editing_context=False`` call; not listing it at all raised ``KeyError``
    from pandas for a frame that met the documented contract.
    """
    row = _cas9_row()
    del row["editing_context"]
    df = pd.DataFrame([row])

    with pytest.raises(ValueError, match=r"needs the column\(s\) \['editing_context'\]"):
        process_base_editing(df, {"g": GENE_STOP_PLUS})

    assert len(process_base_editing(df, {"g": GENE_STOP_PLUS}, editing_context=False)) == 1


@pytest.mark.parametrize("empty", ["", "   ", "\t\n", None, float("nan")])
def test_an_unfiltered_guide_table_names_the_filter_it_is_missing(empty):
    """``identify_base_editing_sites`` leaves this empty for a guide with no
    cytosine in the window, and ``filter_sgrnas_for_base_editing`` is what
    drops those rows.

    Following the chain the column guard above prescribes, without that filter,
    used to reach ``int('')`` and die with a bare pandas ``ValueError: invalid
    literal for int() with base 10: ''`` -- naming neither the column nor the
    step that had been left out. On a real Cas9 table it is routinely half the
    rows: 15 of 30 on the 306 nt CDS this was found with.

    A whitespace-only value counts as empty, because a CSV round trip can pad
    the column and ``int('   ')`` fails exactly as ``int('')`` does. Without
    the ``strip`` in the guard those cases reach the same bare pandas error.
    """
    row = _cas9_row()
    row["editable_cytosines"] = empty

    with pytest.raises(ValueError, match="filter_sgrnas_for_base_editing") as excinfo:
        process_base_editing(pd.DataFrame([row]), {"g": GENE_STOP_PLUS})

    assert "1 of 1 guides have no editable cytosine" in str(excinfo.value)


def test_one_unfiltered_row_refuses_the_whole_table():
    """Refusing is deliberate, rather than dropping the row and carrying on.

    Dropping it would make ``filter_sgrnas_for_base_editing`` optional and
    silently return fewer rows than the caller passed, which is the shape of
    mistake this module raises about everywhere else.
    """
    rows = pd.DataFrame([_cas9_row(), {**_cas9_row(), "editable_cytosines": ""}])

    with pytest.raises(ValueError, match="1 of 2 guides have no editable cytosine"):
        process_base_editing(rows, {"g": GENE_STOP_PLUS})


def test_the_column_guard_names_every_step_of_the_chain():
    """Its message is the only instruction a caller gets, so it has to be whole.

    It named the finders and ``identify_base_editing_sites`` but not
    ``filter_sgrnas_for_base_editing``, so following it literally landed on the
    crash in the test above.
    """
    row = _cas9_row()
    del row["pam"]

    with pytest.raises(ValueError) as excinfo:
        process_base_editing(pd.DataFrame([row]), {"g": GENE_STOP_PLUS})

    message = str(excinfo.value)
    assert "identify_base_editing_sites" in message
    assert "filter_sgrnas_for_base_editing" in message
