#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)

"""Tests for :mod:`teemi.design.combinatorial.simulate`.

The four verdicts are exercised with real sequences rather than stubs: a design
that assembles, one whose parts are ordered differently from the amplicons, one
whose primers have no template, and one that reuses the same sequence twice so
the homology is genuinely ambiguous.
"""

import json
import random
from pathlib import Path

import pytest
from Bio.Seq import Seq

from teemi.design.combinatorial import design as D
from teemi.design.combinatorial.simulate import (
    PYDNA_ASSEMBLY_LIMIT,
    canonical_sequence,
    verify_variant,
)
from teemi.design.combinatorial.sources import SourceResolver
from teemi.design.combinatorial.spec import load_spec
from teemi.design.combinatorial.tm import FootprintTm

from tests.combinatorial_netguard import no_network  # noqa: F401

REPO = Path(__file__).resolve().parents[1]


def random_dna(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGTACGTGC") for _ in range(n))


def inline_spec(parts, topology="linear", homology_length=35) -> dict:
    """A spec whose structure constraints are wide open.

    These tests are about the assembly, not about primer3's floors, so the
    floors are set where no candidate breaches them and the split search has
    nothing to work around.
    """
    return {
        "design": {"name": "simulate fixture", "topology": topology},
        "positions": [
            {
                "id": f"slot{index}",
                "parts": [
                    {"id": part_id, "source": {"type": "inline", "sequence": seq}}
                ],
            }
            for index, (part_id, seq) in enumerate(parts)
        ],
        "primer_design": {
            "target_tm": 55.0,
            "footprint_len_min": 18,
            "junction": {"homology_length": homology_length},
        },
        "constraints": {
            "footprint_tm_window": {"min": 45.0, "max": 75.0},
            "gc_fraction_min": 0.05,
            "gc_fraction_max": 0.95,
            "hairpin_dg_min": -30.0,
            "homodimer_dg_min": -30.0,
            "heterodimer_dg_min": -30.0,
            "three_prime_dimer_dg_min": -30.0,
            "hairpin_tm_max": 99.0,
            "homodimer_tm_max": 99.0,
            "heterodimer_tm_max": 99.0,
        },
    }


def designed(raw: dict, base_dir=None):
    """Everything :func:`verify_variant` needs, for one spec."""
    spec = load_spec(raw)
    records = SourceResolver(base_dir=base_dir).resolve_spec(spec)
    footprints, _ = D.design_footprints(records, spec, FootprintTm(offline=True))
    junctions = D.junctions_for_spec(spec, records, footprints)
    variants = D.variants(spec)
    pairs = {v.id: D.variant_oligos(spec, v, footprints, junctions) for v in variants}
    return spec, records, variants, pairs


# --------------------------------------------------------------------------- #
# The invariant comparison
# --------------------------------------------------------------------------- #


def test_canonical_sequence_is_strand_invariant():
    seq = random_dna(200, 1)
    reverse = str(Seq(seq).reverse_complement())
    assert canonical_sequence(seq, circular=False) == canonical_sequence(
        reverse, circular=False
    )


def test_canonical_sequence_is_rotation_invariant_only_when_circular():
    seq = random_dna(200, 2)
    rotated = seq[57:] + seq[:57]
    assert canonical_sequence(seq, circular=True) == canonical_sequence(
        rotated, circular=True
    )
    assert canonical_sequence(seq, circular=False) != canonical_sequence(
        rotated, circular=False
    )


def test_canonical_sequence_separates_different_molecules():
    left = random_dna(120, 3)
    right = random_dna(120, 4)
    assert canonical_sequence(left, circular=True) != canonical_sequence(
        right, circular=True
    )


def test_assembly_limit_is_pydnas_default_capped_at_the_homology():
    """A limit equal to the intended homology hides every join that came out
    short; pydna's own default is the shortest homology a reaction is taken to
    recombine on."""
    assert PYDNA_ASSEMBLY_LIMIT == 25


# --------------------------------------------------------------------------- #
# The four verdicts
# --------------------------------------------------------------------------- #


def test_linear_design_verifies():
    raw = inline_spec(
        [
            ("a", random_dna(400, 10)),
            ("b", random_dna(400, 11)),
            ("c", random_dna(400, 12)),
        ]
    )
    spec, records, variants, pairs = designed(raw)
    variant = variants[0]
    verification, construct = verify_variant(variant, pairs[variant.id], records, spec)
    assert verification.verdict == "verified"
    assert verification.ok
    assert verification.observed_bp == verification.expected_bp == 1200
    assert construct is not None
    assert len(construct) == 1200
    assert construct.id == variant.id
    assert construct.name == variant.id[:16]


def test_circular_design_verifies_and_is_marked_circular():
    raw = inline_spec(
        [
            ("a", random_dna(400, 20)),
            ("b", random_dna(400, 21)),
            ("c", random_dna(400, 22)),
        ],
        topology="circular",
    )
    spec, records, variants, pairs = designed(raw)
    variant = variants[0]
    verification, construct = verify_variant(variant, pairs[variant.id], records, spec)
    assert verification.verdict == "verified"
    assert construct is not None and construct.circular
    assert len(construct) == 1200


def test_a_different_intended_order_is_a_mismatch():
    """The amplicons are the design's; the expectation is the parts in another
    order, so the product is a real molecule that is not the intended one."""
    raw = inline_spec(
        [
            ("a", random_dna(400, 30)),
            ("b", random_dna(400, 31)),
            ("c", random_dna(400, 32)),
        ]
    )
    spec, records, variants, pairs = designed(raw)
    variant = variants[0]
    reordered = D.Variant(id=variant.id, parts=variant.parts[::-1])
    verification, construct = verify_variant(
        reordered, pairs[variant.id], records, spec
    )
    assert verification.verdict == "mismatch"
    assert construct is None
    assert not verification.ok


def test_primers_with_no_template_give_no_product():
    raw = inline_spec(
        [
            ("a", random_dna(400, 40)),
            ("b", random_dna(400, 41)),
            ("c", random_dna(400, 42)),
        ]
    )
    spec, records, variants, pairs = designed(raw)
    variant = variants[0]
    from pydna.dseqrecord import Dseqrecord

    unrelated = dict(records)
    unrelated["b"] = Dseqrecord(random_dna(400, 999), circular=False)
    verification, construct = verify_variant(
        variant, pairs[variant.id], unrelated, spec
    )
    assert verification.verdict == "no_product"
    assert verification.n_products == 0
    assert construct is None
    assert "no single product" in verification.message


def test_a_repeated_part_sequence_gives_multiple_products():
    """Two parts with the same sequence make the homology ambiguous, so the real
    reaction would give a mixture.  The first product must not be reported."""
    repeated = random_dna(400, 51)
    raw = inline_spec([("a", random_dna(400, 50)), ("b", repeated), ("c", repeated)])
    spec, records, variants, pairs = designed(raw)
    variant = variants[0]
    verification, construct = verify_variant(variant, pairs[variant.id], records, spec)
    assert verification.verdict == "multiple_products"
    assert verification.n_products > 1
    assert construct is None
    assert "mixture" in verification.message


def test_pydna_refusing_a_path_is_reported_as_no_product_with_its_reason():
    """pydna raises when a join leaves the molecule split; that is no product,
    and the reason travels with the verdict rather than out of the call."""
    repeated = random_dna(400, 61)
    raw = inline_spec(
        [("a", random_dna(400, 60)), ("b", repeated), ("c", repeated)],
        topology="circular",
    )
    spec, records, variants, pairs = designed(raw)
    variant = variants[0]
    verification, construct = verify_variant(variant, pairs[variant.id], records, spec)
    assert verification.verdict == "no_product"
    assert "pydna could not assemble" in verification.message
    assert construct is None


# --------------------------------------------------------------------------- #
# The repository's own specs
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "spec_path",
    sorted((REPO / "tests" / "files_for_testing").glob("combinatorial_spec_*.json")),
    ids=lambda p: p.stem,
)
def test_every_variant_of_every_repository_spec_verifies(spec_path):
    raw = json.loads(spec_path.read_text())
    spec, records, variants, pairs = designed(raw, base_dir=REPO)
    for variant in variants:
        verification, construct = verify_variant(
            variant, pairs[variant.id], records, spec
        )
        assert verification.verdict == "verified", verification.message
        assert len(construct) == verification.expected_bp


def test_verified_construct_carries_the_part_and_primer_features():
    """pydna propagates features from the fragments into the product, which is
    what makes the exported GenBank annotated."""
    fixture = (
        REPO / "tests" / "files_for_testing" / "combinatorial_spec_fixed_overhang.json"
    )
    raw = json.loads(fixture.read_text())
    spec, records, variants, pairs = designed(raw, base_dir=REPO)
    variant = variants[0]
    _, construct = verify_variant(variant, pairs[variant.id], records, spec)
    labels = set()
    for feature in construct.features:
        label = feature.qualifiers["label"]
        labels.add(label[0] if isinstance(label, list) else label)
    assert {part.id for part in variant.parts} <= labels
    assert sum(1 for f in construct.features if f.type == "primer_bind") == 6
