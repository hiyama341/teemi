#!/usr/bin/env python
"""Validate the CRISPR modules against references produced outside teemi.

The unit tests elsewhere check that the code does what it is written to do.
These check that what it is written to do is right, against two references that
teemi had no hand in making:

* **CRISPy-web** output for *Streptomyces coelicolor* A3(2) SCO5087, exported
  from the web tool (``tests/files_for_testing/crispy_web/``). CRISPy-web is a
  sound oracle for Cas9 guide lists and, on this plus-strand gene, for
  base-editing mutation calls.
* **A hand curation of SCO7048**, which is on the minus strand. CRISPy-web is
  not reliable for base editing there, so those rows were worked out by hand.
  This is the reference that matters most: placing an edit on a minus-strand
  gene is what the ``sgrna_strand == gene_strand`` comparison in
  ``process_base_editing`` exists to get right, and reading ``sgrna_strand``
  on its own puts every one of these mutations on the wrong residue.

The genome file is the same one StreptoCAD validates against, byte for byte.
"""

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from pydna.dseqrecord import Dseqrecord

from teemi.design.crispr_best import (
    filter_sgrnas_for_base_editing,
    identify_base_editing_sites,
    process_base_editing,
)
from teemi.design.crispr_cas import SgRNAargs, extract_sgRNAs

GENOME = "tests/files_for_testing/gibson/Streptomyces_coelicolor_A3_chromosome.gb"
CRISPY_DIR = "tests/files_for_testing/crispy_web"

# Hand-curated base-editing results for SCO7048, a minus-strand CDS: the guide
# and the amino-acid changes its C-to-T edits produce.
CURATED_SCO7048 = [
    ("TCCCGGCGAATAGCCGGCGC", "R26W, R27*"),
    ("CTCCCGGCGAATAGCCGGCG", "R26W, R27*"),
    ("CCCGGCGAATAGCCGGCGCG", "R27*"),
    ("CTGGCTCCAGCTGTCCCCCA", "W66*, S67N"),
    ("GCTGGCTCCAGCTGTCCCCC", "W66*, S67N"),
    ("CTGGAGCCAGCTGTTCGTCG", "Q68*"),
    ("ACCGGCAGGTCTTCACGCGC", "Q93*"),
    ("GCGACCACACCAGGATGTCG", "V146I, W147*"),
    ("CGACCACACCAGGATGTCGG", "V146I, W147*"),
    ("CGCGACCACACCAGGATGTC", "W147*"),
    ("CCGCGACCACACCAGGATGT", "W147*"),
    ("GGGCAGGCGAGCCGGGCGAC", "Q151*, A152V"),
    ("CCGTCCACCAGGGACTGTTC", "H231Y, Q232*"),
    ("GCCGTCCACCAGGGACTGTT", "H231Y, Q232*"),
    ("CCGGGCCCAGGCCATCACCG", "A305V, Q306*"),
    ("CTACCACCAGACGGGGCTGC", "H321Y, Q322*"),
]


def _split(mutations):
    return [m.strip() for m in str(mutations).split(",") if m.strip()]


@pytest.fixture(scope="module")
def coelicolor():
    return Dseqrecord(SeqIO.read(GENOME, "genbank"))


def _gene_sequences(record, locus_tags):
    """Map each locus tag to its CDS read 5'->3'."""
    return {
        feature.qualifiers["locus_tag"][0]: str(feature.extract(Seq(str(record.seq))))
        for feature in record.features
        if feature.type == "CDS"
        and feature.qualifiers.get("locus_tag", [""])[0] in locus_tags
    }


def _guides(record, locus_tag):
    return extract_sgRNAs(
        SgRNAargs(
            dseqrecord=record,
            locus_tag=[locus_tag],
            cas_type="cas9",
            gc_upper=0.9999,
            gc_lower=0.0001,
            off_target_seed=13,
            off_target_upper=10,
            step=["find", "filter"],
        )
    )


def _stop_codon_guides(record, locus_tag):
    guides = identify_base_editing_sites(_guides(record, locus_tag))
    editable = filter_sgrnas_for_base_editing(guides)
    return process_base_editing(
        editable,
        gene_sequences=_gene_sequences(record, [locus_tag]),
        only_stop_codons=True,
        editing_context=False,
    ).reset_index(drop=True)


@pytest.fixture(scope="module")
def sco5087_stops(coelicolor):
    return _stop_codon_guides(coelicolor, "SCO5087")


@pytest.fixture(scope="module")
def sco7048_stops(coelicolor):
    return _stop_codon_guides(coelicolor, "SCO7048")


# ---------------------------------------------------------------------------
# Base editing against the hand curation of a MINUS-strand gene
# ---------------------------------------------------------------------------


def test_base_editing_matches_hand_curation_on_a_minus_strand_gene(sco7048_stops):
    curated = pd.DataFrame(CURATED_SCO7048, columns=["sgrna", "mutations"])

    merged = sco7048_stops[["sgrna", "mutations"]].merge(
        curated, on="sgrna", suffixes=("_teemi", "_curated")
    )

    # The curation is exhaustive -- it is every stop-codon guide for SCO7048 --
    # so assert set equality, not containment. An inner join on its own passes
    # silently both when teemi stops producing these guides and when it starts
    # producing extra ones.
    assert set(sco7048_stops["sgrna"]) == set(curated["sgrna"]), (
        f"missing: {sorted(set(curated['sgrna']) - set(sco7048_stops['sgrna']))}; "
        f"unexpected: {sorted(set(sco7048_stops['sgrna']) - set(curated['sgrna']))}"
    )
    assert len(merged) == len(curated)

    wrong = {
        row.sgrna: (_split(row.mutations_teemi), _split(row.mutations_curated))
        for row in merged.itertuples()
        if _split(row.mutations_teemi) != _split(row.mutations_curated)
    }
    assert not wrong, f"mutation calls disagree with the curation: {wrong}"


# ---------------------------------------------------------------------------
# Base editing against CRISPy-web on a PLUS-strand gene
# ---------------------------------------------------------------------------


def test_base_editing_matches_crispy_web_on_a_plus_strand_gene(sco5087_stops):
    expected = (
        pd.read_csv(f"{CRISPY_DIR}/BEST_SCO5087_only_stop_codons.csv", index_col=0)[
            ["Sequence", "C to T mutations"]
        ]
        .rename(columns={"Sequence": "sgrna", "C to T mutations": "mutations"})
        .sort_values("sgrna")
        .reset_index(drop=True)
    )
    actual = sco5087_stops[["sgrna", "mutations"]].sort_values("sgrna").reset_index(
        drop=True
    )

    assert list(actual["sgrna"]) == list(expected["sgrna"])
    for sgrna, got, want in zip(
        actual["sgrna"], actual["mutations"], expected["mutations"]
    ):
        assert _split(got) == _split(want), f"guide {sgrna}"


# ---------------------------------------------------------------------------
# Cas9 guide finding against CRISPy-web
# ---------------------------------------------------------------------------


def test_cas9_finds_every_guide_crispy_web_reports(coelicolor):
    """Every CRISPy-web guide must be found, with the same PAM.

    teemi reports a superset: 289 guides against CRISPy-web's 206 for this
    gene, and the 83 extras are not explained by GC, by position in the CDS or
    by duplicates, so they are most likely filters CRISPy-web applies by
    default and teemi does not. That is why this asserts containment rather
    than equality -- it still fails if the finder loses guides or moves a PAM.
    """
    expected = pd.read_csv(f"{CRISPY_DIR}/cas9_SCO5087_0_mismatches.csv", index_col=0)
    actual = _guides(coelicolor, "SCO5087")

    assert actual["sgrna"].is_unique
    assert not set(expected["Sequence"]) - set(actual["sgrna"])

    merged = actual[["sgrna", "pam"]].merge(
        expected[["Sequence", "PAM"]].rename(columns={"Sequence": "sgrna"}), on="sgrna"
    )
    assert len(merged) == len(expected)
    mismatched = merged[merged["pam"] != merged["PAM"]]
    assert mismatched.empty, f"PAM differs for {list(mismatched['sgrna'])}"
