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

"""Does the design assemble into the construct it claims to?

The primers are run against their own templates with ``pydna.amplify.pcr``, the
amplicons are joined with ``pydna.assembly2.Assembly``, and the product is
compared with the plain concatenation of the variant's parts.  Nothing here
designs anything; it only asks whether what was designed works.

Two comparisons are deliberately invariant, because a correct construct can come
back written either way round and, when circular, starting anywhere:

* **strand** -- a sequence and its reverse complement are the same molecule, so
  both are reduced to the lexicographically smaller of the two;
* **rotation** -- a circular construct has no first base, so it is reduced to
  its smallest rotation first (``pydna.utils.smallest_rotation``, Booth's
  algorithm).

The verdict distinguishes four outcomes, and **multiple products is not a
success**: more than one distinct product means the homology is ambiguous, so
the real library would carry a mixture.  Taking the first product would report a
plausible-looking construct for a design that scrambles.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Literal, Mapping, Optional, Sequence, Tuple

from Bio.Seq import Seq
from pydna.amplicon import Amplicon
from pydna.amplify import pcr
from pydna.assembly2 import Assembly
from pydna.dseqrecord import Dseqrecord
from pydna.utils import smallest_rotation

from teemi.design.combinatorial.design import OligoPair, Variant
from teemi.design.combinatorial.spec import DesignSpec

__all__ = [
    "PYDNA_ASSEMBLY_LIMIT",
    "Verdict",
    "Verification",
    "canonical_sequence",
    "verify_variant",
]

#: pydna's own default minimum homology, which is the shortest stretch a Gibson
#: reaction is taken to recombine on.  The limit actually used is this capped at
#: the design's ``homology_length``, so the intended joins are always visible;
#: the design's own L is deliberately *not* used as the limit, because a limit
#: equal to the intended homology hides every join that came out shorter than
#: intended and reports it as no product.
PYDNA_ASSEMBLY_LIMIT = 25

Verdict = Literal["verified", "mismatch", "no_product", "multiple_products"]


@dataclass(frozen=True, slots=True)
class Verification:
    """What the simulated assembly of one variant produced."""

    verdict: Verdict
    message: str
    n_products: int
    expected_bp: int
    observed_bp: Optional[int]

    @property
    def ok(self) -> bool:
        return self.verdict == "verified"

    def to_dict(self) -> dict:
        return {
            "verdict": self.verdict,
            "message": self.message,
            "n_products": self.n_products,
            "expected_bp": self.expected_bp,
            "observed_bp": self.observed_bp,
        }


def canonical_sequence(seq: str, *, circular: bool) -> str:
    """A strand- and (when circular) rotation-invariant key for ``seq``.

    Two sequences describe the same molecule exactly when their keys are equal.
    """
    forward = seq
    reverse = str(Seq(seq).reverse_complement())
    if circular:
        forward, reverse = smallest_rotation(forward), smallest_rotation(reverse)
    return min(forward, reverse)


def verify_variant(
    variant: Variant,
    pairs: Sequence[OligoPair],
    records: Mapping[str, Dseqrecord],
    spec: DesignSpec,
) -> Tuple[Verification, Optional[Dseqrecord]]:
    """Amplify, assemble and compare one variant.

    The second element is the assembled record -- carrying the features pydna
    propagated out of the fragments -- and is ``None`` whenever there is not
    exactly one product, so a caller cannot export a construct the verification
    did not pin down.

    The fragments are handed to pydna in the design's order, which is what the
    spec's ordered positions mean; the question answered is whether that order
    closes into exactly the intended construct and does so uniquely.
    """
    circular = spec.design.topology == "circular"
    expected = "".join(str(records[part.id].seq) for part in variant.parts)

    def outcome(
        verdict: Verdict,
        message: str,
        *,
        n_products: int = 0,
        observed_bp: Optional[int] = None,
    ) -> Verification:
        return Verification(
            verdict=verdict,
            message=f"variant {variant.id}: {message}",
            n_products=n_products,
            expected_bp=len(expected),
            observed_bp=observed_bp,
        )

    try:
        amplicons = _amplify(pairs, records)
    except ValueError as exc:
        return outcome("no_product", str(exc)), None

    limit = min(PYDNA_ASSEMBLY_LIMIT, spec.primer_design.junction.homology_length)
    assembly = Assembly(list(amplicons), limit=limit)
    try:
        products = (
            assembly.assemble_circular() if circular else assembly.assemble_linear()
        )
    except ValueError as exc:
        # pydna refuses a path whose joins leave the molecule split, which happens
        # when the homology is ambiguous enough that two joins land on the same
        # stretch.  That is no product, with a reason.
        return (
            outcome(
                "no_product",
                f"pydna could not assemble the {len(amplicons)} amplicons at "
                f"{limit} bp of homology ({exc})",
            ),
            None,
        )

    distinct: Dict[str, Dseqrecord] = {}
    for product in products:
        distinct.setdefault(
            canonical_sequence(str(product.seq), circular=circular), product
        )

    if not distinct:
        return (
            outcome(
                "no_product",
                f"the {len(amplicons)} amplicons do not assemble at {limit} bp "
                "of homology",
            ),
            None,
        )
    if len(distinct) > 1:
        return (
            outcome(
                "multiple_products",
                f"{len(distinct)} distinct products "
                f"({', '.join(str(len(p)) for p in distinct.values())} bp) -- the "
                "homology is not unique, so the real reaction would give a mixture",
                n_products=len(distinct),
            ),
            None,
        )

    key, product = next(iter(distinct.items()))
    if key != canonical_sequence(expected, circular=circular):
        return (
            outcome(
                "mismatch",
                f"the assembled product is {len(product)} bp but the intended "
                f"construct is {len(expected)} bp, or differs in sequence",
                n_products=1,
                observed_bp=len(product),
            ),
            None,
        )

    product.id = variant.id
    product.name = variant.id[:16]  # GenBank LOCUS field
    product.description = spec.design.name
    product.annotations["molecule_type"] = "DNA"
    return (
        outcome(
            "verified",
            f"{len(product)} bp, assembles to the intended construct",
            n_products=1,
            observed_bp=len(product),
        ),
        product,
    )


def _amplify(
    pairs: Sequence[OligoPair], records: Mapping[str, Dseqrecord]
) -> Tuple[Amplicon, ...]:
    """One amplicon per PCR, named after its part.

    ``pcr`` raises :class:`ValueError` when a primer pair does not give exactly
    one product on its template -- no annealing site, or several -- which is a
    real failure of the design and is reported as such rather than worked
    around.  The name is set because pydna labels the assembly figure with it.
    """
    amplicons: List[Amplicon] = []
    for pair in pairs:
        try:
            amplicon = pcr(pair.forward, pair.reverse, records[pair.part_id])
        except ValueError as exc:
            raise ValueError(
                f"primers {pair.forward.id} / {pair.reverse.id} give no single "
                f"product on part {pair.part_id!r} ({exc})"
            ) from exc
        amplicon.name = pair.part_id[:16]
        amplicons.append(amplicon)
    return tuple(amplicons)
