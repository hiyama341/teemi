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

The verdict distinguishes five outcomes, and **multiple products is not a
success**: more than one distinct product means the homology is ambiguous, so
the real library would carry a mixture.  Taking the first product would report a
plausible-looking construct for a design that scrambles.  The fifth,
``assembly_failed``, is pydna declining the assembly at all; it is a verdict
rather than an exception so that :func:`~teemi.design.combinatorial.plan`
always returns a result.
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
from teemi.design.combinatorial.spec import DesignSpec, Part

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

Verdict = Literal[
    "verified",
    "mismatch",
    "no_product",
    "multiple_products",
    #: pydna refused the assembly outright -- incompatible ends, which it
    #: reports as a bare ``TypeError``.  A fact about the design, so it is a
    #: verdict like the others rather than an exception out of ``plan``.
    "assembly_failed",
]


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

    An ``ordered`` part has no PCR, so it enters the assembly as the plain
    :class:`~pydna.dseqrecord.Dseqrecord` it will be bought as, in its place
    between the amplicons.  That is the whole point of the mode and it is what
    is checked here: the synthetic fragment is joined by its neighbours' tails
    without ever being a template, and the product still has to equal the
    intended construct base for base.
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

    pieces = _assembly_pieces(variant, amplicons, records, circular=circular)
    limit = min(PYDNA_ASSEMBLY_LIMIT, spec.primer_design.junction.homology_length)
    assembly = Assembly(list(pieces), limit=limit)
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
                f"pydna could not assemble the {len(pieces)} fragments at "
                f"{limit} bp of homology ({exc})",
            ),
            None,
        )

    # pydna's linear assembly also returns the degenerate path that is just the
    # shared region between two amplicons -- a 30 bp "product" beside the 800 bp
    # one, which made every two-fragment linear design report an ambiguity it
    # does not have.  A real join of n amplicons is at least as long as its
    # longest member, so anything shorter is that artefact and not a product.
    floor = max((len(piece) for piece in pieces), default=0)
    distinct: Dict[str, Dseqrecord] = {}
    for product in products:
        if len(product) < floor:
            continue
        distinct.setdefault(
            canonical_sequence(str(product.seq), circular=circular), product
        )

    if not distinct:
        return (
            outcome(
                "no_product",
                f"the {len(pieces)} fragments do not assemble at {limit} bp "
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


def _assembly_pieces(
    variant: Variant,
    amplicons: Sequence[Amplicon],
    records: Mapping[str, Dseqrecord],
    *,
    circular: bool,
) -> Tuple[Dseqrecord, ...]:
    """The molecules the reaction actually contains, in construct order.

    Two kinds, interleaved: an amplicon for every amplified fragment, and one
    record for every run of ``ordered`` parts, which are bought rather than
    made and so enter the tube as themselves.  An ``in_tail`` part is neither
    -- it is inside its neighbours' primers -- so it contributes no molecule
    here while still being in :attr:`Variant.parts`, which is what the product
    is compared against.

    **A run of ordered parts is ONE molecule.**  Two synthetic fragments side
    by side share no sequence, so handing them over separately asks pydna to
    join them on nothing and it correctly refuses; they are adjacent in the
    construct, so they are bought and assembled as a single piece.  That is the
    same grouping the primers were designed against, where the run's
    concatenation is what the two flanking tails reach across.

    Walked over :attr:`Variant.parts`, so the order is the construct's and not
    the order two separate lists happened to be built in.  Part ids are unique
    across a whole spec, which is what makes the membership test here
    unambiguous: one id cannot be amplified at one slot and ordered at another.
    """
    amplified = iter(amplicons)
    ordered = {part.id for part in variant.ordered}
    fragments = {part.id for part in (variant.fragments or variant.parts)}
    pieces: List[Dseqrecord] = []
    # (first part id, joined sequence) per run of ordered parts, in order.
    runs: List[Tuple[str, str]] = []
    slots: List[Optional[int]] = []
    run: List[Part] = []

    def flush() -> None:
        if not run:
            return
        slots.append(len(runs))
        runs.append(
            (run[0].id, "".join(str(records[part.id].seq) for part in run))
        )
        run.clear()

    for part in variant.parts:
        if part.id in ordered:
            run.append(part)
            continue
        flush()
        if part.id in fragments:
            slots.append(None)
    flush()

    # A circular construct has no first base, so a run at each end is one run
    # spanning the join: the trailing one comes first, as the circle reads it.
    wraps = bool(slots) and slots[0] == 0 and slots[-1] == len(runs) - 1
    if circular and len(runs) > 1 and wraps:
        last_id, last_seq = runs.pop()
        first_id, first_seq = runs[0]
        runs[0] = (last_id, last_seq + first_seq)
        slots = slots[:-1]

    for slot in slots:
        if slot is None:
            pieces.append(next(amplified))
            continue
        part_id, sequence = runs[slot]
        record = Dseqrecord(sequence)
        record.id = part_id
        record.name = part_id[:16]
        pieces.append(record)
    return tuple(pieces)


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
