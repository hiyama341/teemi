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

"""The one entry point, the result, and the three tables.

:func:`plan` runs the whole pipeline -- spec, sources, footprints, junction
splits, primers, PCRs, NEB temperatures, assembly verification -- and returns a
:class:`DesignResult`.  It writes no files: that is
:mod:`teemi.design.combinatorial.export`, which the caller invokes with
``result.designs``.

It is deterministic.  The same spec gives the same result, so ``hash(spec)``
keys a cache, and :meth:`DesignResult.to_dict` serialises under a bare
``json.dumps`` -- no pydna objects, no ``Seq``, no sets, no tuple keys.
``include_designs`` does not change it: every variant is verified either way,
and the flag only decides whether the assembled records are kept for export.

Every failure that is about the *input* -- an unparseable spec, an unreadable
part -- is caught here, at this one boundary, and becomes a critical diagnostic
rather than an exception out of :func:`plan`.  Every failure that is about the
*design* is a constraint breach, and the design still comes back carrying it.

The table columns keep the two Tm quantities apart by name, because conflating
them has caused real bugs.  A row carries machine-readable keys and the tables
render them through :data:`COLUMN_LABELS`, which is where the unit lives: a
primer row has ``footprint_Tm (C)`` and ``oligo_Tm (C)``, a PCR row has those
split ``_fwd``/``_rev``, and the annealing temperature is
``Ta_from_footprints (C)`` -- named for the pair it is allowed to come from,
which is never the tailed oligos.  Both Tm of a row come from one engine, named
once as ``tm_engine``; Ta keeps its own, because its model is a different one.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import pandas as pd
from pydantic import ValidationError
from pydna.dseqrecord import Dseqrecord

from teemi.design.combinatorial.design import (
    BLOCKING_BOUNDS,
    STRUCTURE_BOUNDS,
    Diagnostic,
    apply_slot_flags,
    plan_layout,
    FootprintError,
    Footprints,
    JunctionDesign,
    JunctionError,
    OligoPair,
    OligoStructures,
    PairStructures,
    PcrCheck,
    Variant,
    Violation,
    check_pcrs,
    design_footprints,
    junctions_for_spec,
    oligo_structures,
    oligo_violations,
    variant_oligos,
    variants,
)
from teemi.design.combinatorial.simulate import Verification, verify_variant
from teemi.design.combinatorial.sources import SourceError, SourceResolver
from teemi.design.combinatorial.spec import ID_SEPARATOR, DesignSpec, load_spec
from teemi.design.combinatorial.tm import (
    CONDITIONS,
    FootprintTm,
    PcrOligos,
    PcrTemperatures,
    Tm,
    neb_pair_temperatures,
)

__all__ = [
    "InfeasibleDesignError",
    "PrimerRecord",
    "PcrRecord",
    "VariantPlan",
    "DesignResult",
    "plan",
    "variant_overview_dataframe",
    "primer_dataframe",
    "pcr_dataframe",
]

#: The readable projection of a primer row; everything else is behind ``full``.
PRIMER_COLUMNS = (
    "primer",
    "primer_id",
    "part_id",
    "direction",
    "n_uses",
    "length_nt",
    "sequence",
    "footprint_len",
    "footprint_tm_c",
    "hairpin_folded_at_ta",
    "hairpin_3prime_paired",
    "hairpin_dg",
    "homodimer_dg",
    "primer3_clamped",
    "flagged_by",
    "price",
)

#: The readable projection of a PCR row.
PCR_COLUMNS = (
    "part_id",
    "forward_primer",
    "reverse_primer",
    "amplicon_bp",
    "ta_from_footprints",
    "footprint_tm_forward",
    "footprint_tm_reverse",
    "hairpin_folded_at_ta",
    "hairpin_3prime_paired",
    "heterodimer_dg",
    "three_prime_dimer_dg",
    "homodimer_dg_worst",
    "primer3_clamped",
    "flagged_by",
)

#: Headings for the tables: the key a row carries -> what a reader should see.
#:
#: The unit belongs in the heading.  A dG read as cal/mol, or a structure Tm
#: read as a margin against Ta, is the kind of mistake these names are here to
#: stop, and a heading saying ``heterodimer_fwd_rev`` cannot be mistaken for a
#: single primer's homodimer.  Both Tm columns of a row come from one engine,
#: so one ``tm_engine`` heading covers them and the second is dropped.
COLUMN_LABELS = {
    "footprint_tm_c": "footprint_Tm (C)",
    "oligo_tm_c": "oligo_Tm (C)",
    "footprint_tm_engine": "tm_engine",
    "footprint_tm_forward": "footprint_Tm_fwd (C)",
    "footprint_tm_reverse": "footprint_Tm_rev (C)",
    "oligo_tm_forward": "oligo_Tm_fwd (C)",
    "oligo_tm_reverse": "oligo_Tm_rev (C)",
    "ta_from_footprints": "Ta_from_footprints (C)",
    "hairpin_dg": "hairpin_deltaG (kcal/mol)",
    "hairpin_tm": "hairpin_Tm (C)",
    "homodimer_dg": "homodimer_deltaG (kcal/mol)",
    "homodimer_tm": "homodimer_Tm (C)",
    "heterodimer_dg": "heterodimer_fwd_rev_deltaG (kcal/mol)",
    "heterodimer_tm": "heterodimer_fwd_rev_Tm (C)",
    "three_prime_dimer_dg": "three_prime_dimer_deltaG (kcal/mol)",
    "homodimer_dg_worst": "worst_primer_homodimer_deltaG (kcal/mol)",
    "primer3_clamped": "not_analysable_over_60nt",
    "hairpin_folded_at_ta": "hairpin_folded_at_Ta",
    "hairpin_3prime_paired": "hairpin_blocks_3prime",
}


class InfeasibleDesignError(ValueError):
    """Raised by :meth:`DesignResult.raise_for_status`."""

    def __init__(self, diagnostics: Sequence[Diagnostic]) -> None:
        super().__init__(
            f"{len(diagnostics)} blocking diagnostic(s): "
            + "; ".join(d.message for d in diagnostics)
        )
        self.diagnostics = tuple(diagnostics)


def _shared_id(uses: Sequence[str]) -> str:
    """The name of an oligo, given every junction that asked for it.

    A primer id is ``<left>__<right>_<direction>``, and an oligo two junctions
    share differs from each of them in one token -- two parts with the same end
    give the same tail on the same footprint.  Joining the distinct tokens
    position by position names every context the molecule serves and privileges
    none of them, so a reaction on ``AAA17732`` is not told to use a tube
    labelled ``AKF02530``.  Ids whose shape does not line up keep the first
    name rather than being spliced into a wrong one.
    """
    if len(uses) == 1:
        return uses[0]
    tokens = [use.split(ID_SEPARATOR) for use in uses]
    if len({len(parts) for parts in tokens}) != 1:
        return uses[0]
    columns = list(zip(*tokens))
    varying = [index for index, column in enumerate(columns) if len(set(column)) > 1]
    # Exactly one column may be merged, and never the last one.  Two differing
    # columns would name the cross product -- junctions the design need not
    # contain -- and the last column carries the ``_fw``/``_rv`` suffix, so
    # merging it buries the direction mid-name (``A__X_rv+Y_rv``).  Either way
    # the id falls back to one real name and ``used_as`` carries the rest.
    if len(varying) != 1 or varying[0] == len(columns) - 1:
        return uses[0]
    return ID_SEPARATOR.join(
        "+".join(sorted(set(column))) if index == varying[0] else column[0]
        for index, column in enumerate(columns)
    )


def _one_engine(engines: Iterable[str]) -> str:
    """The engine behind a column, named once instead of in every row.

    A per-row engine column said the same thing 142 times, because a degraded
    run gives up NEB for everything rather than half way.  If that ever stops
    being true this says so rather than quietly picking one.
    """
    seen = sorted(set(engines))
    return "/".join(seen) if seen else ""


def _severity(violation: Violation) -> str:
    """``error`` only where the design cannot be built; otherwise a flag.

    One rule, read off :data:`BLOCKING_BOUNDS`.  A thermodynamic prediction is
    a flag however hopeless it looks, and an unorderable primer is an error
    however good its thermodynamics are.  Whether anything could have avoided
    the breach is still recorded, in the message -- it says whether re-running
    with a search would help -- but it no longer decides whether the design is
    refused.
    """
    return "error" if violation.constraint in BLOCKING_BOUNDS else "warning"


def _flag_severity(violations: Sequence[Violation]) -> str:
    """The worst severity among a row's breaches, or ``""`` for none."""
    if not violations:
        return ""
    return _worst(_severity(violation) for violation in violations)


def _worst(severities: Iterable[str]) -> str:
    """``error`` beats ``warning`` beats nothing."""
    seen = set(severities)
    return "error" if "error" in seen else ("warning" if "warning" in seen else "")


def _dg_flagged(violations: Sequence[Violation]) -> bool:
    return any(v.constraint.endswith("_dg_min") for v in violations)


def _violated(violations: Sequence[Violation]) -> str:
    return ",".join(sorted({v.constraint for v in violations}))


@dataclass(frozen=True, slots=True)
class PrimerRecord:
    """One oligo to order, with the numbers that decide whether to order it.

    An oligo is identified by its **sequence**, not by its name: two junctions
    that happen to need the same molecule share one record, one tube and one
    line on the invoice.  :attr:`uses` lists every name that asked for it, so a
    reusable primer is visible as such instead of being ordered twice.

    :attr:`footprint_tm` and :attr:`oligo_tm` are the authoritative batched
    numbers for the annealing region and for the whole tailed molecule.  They
    are different quantities from the same engine and are never combined.
    """

    id: str
    part_id: str
    direction: str
    sequence: str
    footprint: str
    footprint_tm: Tm
    oligo_tm: Tm
    structures: OligoStructures
    violations: Tuple[Violation, ...]
    price: float
    #: Position in the order sheet, from 1.  ``label`` renders it zero-padded
    #: to three digits -- ``007`` -- because that is what gets said at a bench
    #: and ``PRO_ANID_TEF1__gene_KP2612_000062_rv`` is not.
    #:
    #: Padded and a string, so it sorts and lines up in a column.  A
    #: spreadsheet will read ``007`` as the number 7 and drop the padding on
    #: paste, which costs nothing here but is worth knowing before a vendor
    #: form is filled from one.
    #:
    #: A LABEL, NOT AN IDENTITY.  It is the molecule's position in *this* plan:
    #: deterministic, so the same spec always numbers the same way, but
    #: positional, so inserting a gene renumbers everything after it.  The
    #: sequence is the identity, and it is in the row beside the number for
    #: exactly that reason -- check it before ordering against a number you
    #: wrote down from an earlier run.
    number: int = 0
    uses: Tuple[str, ...] = ()
    #: Every part this molecule amplifies, unjoined.  :attr:`part_id` is the
    #: readable join of these; anything that has to name a real part uses this.
    part_ids: Tuple[str, ...] = ()
    #: Breached bounds no candidate in the sweep that produced this molecule
    #: avoided.  For an untailed terminal primer there was no sweep and no
    #: tail, so every structure bound it breaches is in here.
    unavoidable: Tuple[str, ...] = ()
    #: The coolest annealing temperature this molecule is used at, which is the
    #: one its hairpin is most folded at.  A shared oligo can sit in reactions
    #: with different Ta, so the worst of them is the one to report.
    anneal_ta: Optional[float] = None
    #: Whether a junction split produced this molecule.  A tailed primer's
    #: breach is reported as a diagnostic against its junction, where the split
    #: that caused it lives; the row still shows it either way.
    tailed: bool = False

    @property
    def label(self) -> str:
        """``001``: short enough to say out loud."""
        return f"{self.number:03d}"

    def to_dict(self) -> dict:
        return {
            "primer": self.label,
            "primer_id": self.id,
            "part_id": self.part_id,
            "direction": self.direction,
            "n_uses": len(self.uses),
            "used_as": ",".join(self.uses),
            "sequence": self.sequence,
            "length_nt": len(self.sequence),
            "footprint": self.footprint,
            "footprint_len": len(self.footprint),
            "footprint_tm_c": self.footprint_tm.celsius,
            "footprint_tm_engine": self.footprint_tm.engine,
            "oligo_tm_c": self.oligo_tm.celsius,
            "oligo_tm_engine": self.oligo_tm.engine,
            **self.structures.to_dict(),
            # What a Tm cannot say: how much of the primer is tied up when it
            # has to anneal.  Blank rather than 0 where no Ta is known, so an
            # absent number is never read as "none folded".
            "hairpin_folded_at_ta": (
                None
                if self.anneal_ta is None
                else self.structures.folded_at(self.anneal_ta)
            ),
            "dg_flagged": _dg_flagged(self.violations),
            "flagged_by": _violated(self.violations),
            "flag_severity": _flag_severity(self.violations),
            "violations": [v.to_dict() for v in self.violations],
            "price": self.price,
        }


@dataclass(frozen=True, slots=True)
class PcrRecord:
    """One reaction to run: its two primers, its temperatures, its structures.

    :attr:`structures` holds the *pair* measurements -- the two primers against
    each other.  A single primer's own hairpin and homodimer belong to the
    primer, not to the reaction, so they arrive here as
    :attr:`primer_violations` and as :attr:`homodimer_dg_worst`: without them a
    reaction carrying a flagged oligo would read as unflagged.
    """

    part_id: str
    forward_id: str
    reverse_id: str
    amplicon_bp: int
    temperatures: PcrTemperatures
    structures: PairStructures
    violations: Tuple[Violation, ...]
    #: The two primers' order-sheet numbers, so a reaction can be read without
    #: a sixty-character name in each cell.
    forward_label: str = ""
    reverse_label: str = ""
    primer_violations: Tuple[Violation, ...] = ()
    homodimer_dg_worst: Optional[float] = None
    #: The worse of the two primers' own verdicts, each scored against its own
    #: unavoidable set.  Pooling the two sets and scoring the pooled breaches
    #: against the union let a bound that is unavoidable for one primer
    #: downgrade the other primer's avoidable breach of the same bound.
    primer_severity: str = ""
    #: The worse of the two primers' hairpins at this reaction's own Ta, and
    #: whether either of them has its 3' end inside a stem.
    hairpin_folded_at_ta: Optional[float] = None
    hairpin_3prime_paired: bool = False

    def to_dict(self) -> dict:
        return {
            "part_id": self.part_id,
            "forward_primer": self.forward_label,
            "reverse_primer": self.reverse_label,
            "forward_primer_id": self.forward_id,
            "reverse_primer_id": self.reverse_id,
            "amplicon_bp": self.amplicon_bp,
            **self.temperatures.to_dict(),
            **self.structures.to_dict(),
            "homodimer_dg_worst": self.homodimer_dg_worst,
            "hairpin_folded_at_ta": self.hairpin_folded_at_ta,
            "hairpin_3prime_paired": self.hairpin_3prime_paired,
            "dg_flagged": _dg_flagged(self.violations + self.primer_violations),
            "flagged_by": _violated(self.violations + self.primer_violations),
            # The pair bounds are errors: the sweep optimises each junction on
            # its own and never searched this forward-against-reverse pairing,
            # so "nothing avoids it" has not been established for them.
            "flag_severity": "error" if self.violations else self.primer_severity,
            # Both sets, because ``flagged_by`` names bounds from both and a
            # consumer reading it then needs the observed value and the limit.
            # ``subject`` says which it was: a primer, or the pair.
            "violations": [
                v.to_dict() for v in self.violations + self.primer_violations
            ],
        }


@dataclass(frozen=True, slots=True)
class VariantPlan:
    """One member of the library: what it is made of and whether it assembles.

    :attr:`construct` is the assembled, annotated record, kept only when
    :func:`plan` was asked for ``include_designs`` and the verification found
    exactly one product.  It is the object :mod:`.export` writes.
    """

    id: str
    choices: Tuple[Tuple[str, str], ...]
    fragment_ids: Tuple[str, ...]
    length_bp: int
    verification: Verification
    construct: Optional[Dseqrecord] = None

    def to_dict(self) -> dict:
        return {
            "id": self.id,
            "choices": [list(choice) for choice in self.choices],
            "fragment_ids": list(self.fragment_ids),
            "length_bp": self.length_bp,
            "verification": self.verification.to_dict(),
        }


@dataclass(frozen=True, slots=True)
class DesignResult:
    """Everything :func:`plan` worked out, in one JSON-serialisable object."""

    spec: Optional[DesignSpec]
    diagnostics: Tuple[Diagnostic, ...] = ()
    footprints: Tuple[Footprints, ...] = ()
    junctions: Tuple[JunctionDesign, ...] = ()
    primers: Tuple[PrimerRecord, ...] = ()
    pcrs: Tuple[PcrRecord, ...] = ()
    variants: Tuple[VariantPlan, ...] = ()

    @property
    def errors(self) -> Tuple[Diagnostic, ...]:
        return tuple(d for d in self.diagnostics if d.severity in ("error", "critical"))

    @property
    def feasible(self) -> bool:
        """No blocking diagnostic: every bound met and every variant verified."""
        return not self.errors

    @property
    def designs(self) -> Tuple[VariantPlan, ...]:
        """The variants whose assembled construct is available to export."""
        return tuple(v for v in self.variants if v.construct is not None)

    @property
    def price(self) -> float:
        """What the oligos cost at the spec's price per base."""
        return sum(p.price for p in self.primers)

    def raise_for_status(self) -> None:
        if not self.feasible:
            raise InfeasibleDesignError(self.errors)

    def summary(self) -> str:
        if self.spec is None:
            return "no design: " + "; ".join(d.message for d in self.diagnostics)
        design = self.spec.design
        verdicts = Counter(v.verification.verdict for v in self.variants)
        severities = Counter(d.severity for d in self.diagnostics)
        longest = max((len(p.sequence) for p in self.primers), default=0)
        return "\n".join(
            (
                f"{design.name} [{design.topology}] "
                f"{'feasible' if self.feasible else 'INFEASIBLE'}",
                f"  {len(self.variants)} variants from "
                f"{len(self.spec.positions)} positions, "
                f"{len(self.footprints)} parts",
                f"  {len(self.junctions)} junctions, "
                f"{sum(1 for j in self.junctions if j.feasible)} feasible, "
                f"homology {self.spec.primer_design.junction.homology_length} bp"
                + (
                    " ("
                    + ", ".join(
                        f"{n} by {who}"
                        for who, n in sorted(
                            Counter(j.designed_by for j in self.junctions).items()
                        )
                    )
                    + ")"
                    if self.junctions
                    else ""
                ),
                "  temperatures: "
                + ", ".join(
                    f"{what} {engine}"
                    for what, engine in (
                        (
                            "Tm",
                            _one_engine(
                                p.footprint_tm.engine for p in self.primers
                            ),
                        ),
                        (
                            "Ta",
                            _one_engine(
                                p.temperatures.ta_from_footprints.engine
                                for p in self.pcrs
                            ),
                        ),
                        ("structures", "primer3"),
                    )
                    if engine
                ),
                f"  {len(self.primers)} primers (longest {longest} nt), "
                f"{len(self.pcrs)} PCRs, "
                f"{self.price:.0f} {self.spec.objective.currency}",
                "  assembly: "
                + ", ".join(
                    f"{n} {verdict}" for verdict, n in sorted(verdicts.items())
                ),
                "  diagnostics: "
                + (
                    ", ".join(f"{n} {sev}" for sev, n in sorted(severities.items()))
                    or "none"
                ),
            )
        )

    def to_dict(self) -> dict:
        design = None
        cost = None
        if self.spec is not None:
            design = {
                "name": self.spec.design.name,
                "topology": self.spec.design.topology,
                "description": self.spec.design.description,
                "n_positions": len(self.spec.positions),
                "n_variants": self.spec.n_variants,
                "homology_length": self.spec.primer_design.junction.homology_length,
            }
            cost = {
                "bases": sum(len(p.sequence) for p in self.primers),
                "price": self.price,
                "price_per_base": self.spec.objective.price_per_base,
                "currency": self.spec.objective.currency,
            }
        return {
            "design": design,
            "feasible": self.feasible,
            "cost": cost,
            "variants": [v.to_dict() for v in self.variants],
            "footprints": [f.to_dict() for f in self.footprints],
            "junctions": [j.to_dict() for j in self.junctions],
            "primers": [p.to_dict() for p in self.primers],
            "pcrs": [p.to_dict() for p in self.pcrs],
            "diagnostics": [d.to_dict() for d in self.diagnostics],
        }


def plan(
    spec: Union[Mapping, DesignSpec],
    *,
    sources: Optional[SourceResolver] = None,
    offline: bool = False,
    include_designs: bool = False,
) -> DesignResult:
    """Design the whole library and verify it.

    ``sources`` is the filesystem boundary; without one, only ``inline`` part
    sources resolve.  ``offline=True`` makes zero network requests and every
    temperature comes from the local pydna models.  ``include_designs=True``
    keeps the assembled records on the result for :mod:`.export`; it changes
    nothing in :meth:`DesignResult.to_dict`.
    """
    try:
        validated = spec if isinstance(spec, DesignSpec) else load_spec(spec)
    except ValidationError as exc:
        return _stopped(None, "invalid_spec", f"the design spec is not valid: {exc}")

    resolver = SourceResolver() if sources is None else sources
    try:
        records = resolver.resolve_spec(validated)
    except SourceError as exc:
        return _stopped(validated, "unresolved_source", str(exc), part_id=exc.part_id)

    records, flag_notes = apply_slot_flags(records, validated)

    layout = plan_layout(validated, records)

    engine = FootprintTm(
        offline=offline, prodcode=validated.primer_design.polymerase
    )
    try:
        footprints, diagnostics = design_footprints(records, validated, engine, layout)
        junctions = junctions_for_spec(validated, records, footprints, layout)
    except FootprintError as exc:
        return _stopped(validated, "no_footprint", str(exc), part_id=exc.part_id)
    except JunctionError as exc:
        return _stopped(
            validated,
            "impossible_junction",
            str(exc),
            junction_id=exc.junction_id,
        )

    library = variants(validated, layout)
    oligos = {
        variant.id: variant_oligos(validated, variant, footprints, junctions)
        for variant in library
    }
    pairs = [pair for variant_pairs in oligos.values() for pair in variant_pairs]
    checks = check_pcrs(pairs, validated.constraints)
    temperatures, batch_failures = neb_pair_temperatures(
        [pair.oligos for pair in pairs],
        offline=offline,
        prodcode=validated.primer_design.polymerase,
    )

    primers, oligo_of = _primer_records(pairs, junctions, temperatures, validated)
    by_id = {primer.id: primer for primer in primers}
    pcrs = _pcr_records(pairs, checks, temperatures, records, oligo_of, by_id)
    plans, verification_diagnostics = _variant_plans(
        library,
        oligos,
        records,
        validated,
        oligo_of,
        include_designs=include_designs,
    )

    breaches = (
        [
            Diagnostic.from_violation(
                violation,
                junction_id=junction.id,
                severity=_severity(violation),
                exhausted=(
                    f"no split avoids it ({len(junction.profile)} tried)"
                    if (violation.constraint, violation.subject)
                    in junction.unavoidable
                    else None
                ),
            )
            for junction in junctions.values()
            for violation in junction.chosen.violations
        ]
        + [
            Diagnostic.from_violation(
                violation, part_id=check.part_id, severity=_severity(violation)
            )
            for check in checks
            for violation in check.violations
        ]
        + [
            # Only the primers no junction produced: a tailed primer's breach is
            # already a diagnostic against its junction above, and reporting it
            # twice would double every count.  An untailed terminal primer is
            # the footprint itself -- no split, no tail, one possible molecule --
            # so a structure bound it breaches was never avoidable and the
            # search had exactly one option.
            Diagnostic.from_violation(
                violation,
                # Only when it is one real part: the composite ``part_id`` of a
                # shared molecule matches no Part.id anywhere.
                part_id=primer.part_ids[0] if len(primer.part_ids) == 1 else None,
                primer_id=primer.id,
                severity=_severity(violation),
                exhausted=(
                    "the footprint is fixed and this primer has no tail"
                    if violation.constraint in STRUCTURE_BOUNDS
                    else None
                ),
            )
            for primer in primers
            if not primer.tailed
            for violation in primer.violations
        ]
    )
    unreachable = [
        Diagnostic(severity="warning", code="neb_unreachable", message=message)
        for message in engine.failures + batch_failures
    ]

    return DesignResult(
        spec=validated,
        diagnostics=(
            tuple(flag_notes)
            + layout.diagnostics
            + tuple(diagnostics)
            + tuple(breaches)
            + verification_diagnostics
            + tuple(unreachable)
        ),
        footprints=tuple(footprints.values()),
        junctions=tuple(junctions.values()),
        primers=primers,
        pcrs=pcrs,
        variants=plans,
    )


def _stopped(
    spec: Optional[DesignSpec], code: str, message: str, **location: str
) -> DesignResult:
    """A result that says why the pipeline could not start.

    An unparseable spec or an unreadable part is bad *input*, not an infeasible
    design, so nothing downstream ran and there is nothing to report but the
    reason.
    """
    return DesignResult(
        spec=spec,
        diagnostics=(
            Diagnostic(severity="critical", code=code, message=message, **location),
        ),
    )


def _primer_records(
    pairs: Iterable[OligoPair],
    junctions: Mapping[Tuple[str, str], JunctionDesign],
    temperatures: Mapping[PcrOligos, PcrTemperatures],
    spec: DesignSpec,
) -> Tuple[Tuple[PrimerRecord, ...], Dict[str, str]]:
    """One record per distinct oligo to order, plus the name-to-oligo map.

    Two junctions can need the same molecule -- parts that share an end give
    identical tails on identical footprints -- so records are keyed by sequence
    and every name that asked for one is collected in its ``uses``.  The second
    return value maps each of those names onto the id of the oligo that serves
    it, which is what the PCR table references.

    Every record carries the bounds its own molecule breaches, tailed or not: a
    row that showed a homodimer past the floor while reporting itself unflagged
    is the contradiction between these tables that this field exists to remove.
    What differs is where it is *reported* -- a tailed primer's breach becomes a
    diagnostic against its junction, where the split that caused it lives, and
    :attr:`PrimerRecord.tailed` is how :func:`plan` tells the two apart without
    recomputing which sequences came out of a sweep.
    """
    # Sequence -> the bounds that no sweep needing this molecule could avoid.
    # Several junctions can land on one molecule, and their sweeps are separate,
    # so the sets are intersected: if any one of them had a candidate that
    # avoided a bound, there is a choice to make and the breach stays an error.
    # A dict comprehension here would instead keep whichever junction happened
    # to be enumerated last.
    tailed: Dict[str, frozenset] = {}
    for junction in junctions.values():
        for subject, sequence in (
            ("forward_primer", junction.chosen.forward_primer),
            ("reverse_primer", junction.chosen.reverse_primer),
        ):
            # Only the pairs recorded against the role this sequence played.  A
            # junction's set covers both its primers, and scoring one molecule
            # against the other's evidence is the mistake this filter removes.
            avoided = frozenset(
                name for name, who in junction.unavoidable if who == subject
            )
            previous = tailed.get(sequence)
            tailed[sequence] = avoided if previous is None else previous & avoided
    records: Dict[str, PrimerRecord] = {}
    uses: Dict[str, List[str]] = {}
    parts: Dict[str, List[str]] = {}
    directions: Dict[str, List[str]] = {}
    coolest_ta: Dict[str, float] = {}
    for pair in pairs:
        temps = temperatures[pair.oligos]
        for direction, primer, footprint, footprint_tm, oligo_tm in (
            (
                "forward",
                pair.forward,
                pair.oligos.forward_footprint,
                temps.footprint_tm_forward,
                temps.oligo_tm_forward,
            ),
            (
                "reverse",
                pair.reverse,
                pair.oligos.reverse_footprint,
                temps.footprint_tm_reverse,
                temps.oligo_tm_reverse,
            ),
        ):
            sequence = str(primer.seq)
            this_ta = temps.ta_from_footprints.celsius
            coolest_ta[sequence] = min(
                coolest_ta.get(sequence, this_ta), this_ta
            )
            if sequence in records:
                if primer.id not in uses[sequence]:
                    uses[sequence].append(primer.id)
                if pair.part_id not in parts[sequence]:
                    parts[sequence].append(pair.part_id)
                if direction not in directions[sequence]:
                    directions[sequence].append(direction)
                continue
            uses[sequence] = [primer.id]
            parts[sequence] = [pair.part_id]
            directions[sequence] = [direction]
            structures = oligo_structures(sequence, CONDITIONS)
            records[sequence] = PrimerRecord(
                id=primer.id,
                part_id=pair.part_id,
                direction=direction,
                sequence=sequence,
                footprint=footprint,
                footprint_tm=footprint_tm,
                oligo_tm=oligo_tm,
                structures=structures,
                violations=oligo_violations(
                    f"{direction}_primer",
                    sequence,
                    structures,
                    spec.constraints,
                    footprint,
                ),
                price=len(sequence) * spec.objective.price_per_base,
                tailed=sequence in tailed,
            )
    # Numbered in the order the walk first met each molecule, which follows the
    # construct, so 001 is the first fragment's forward primer.  dict keeps
    # insertion order, and the walk is deterministic, so the same spec numbers
    # the same way every run -- which is the only reason a number is safe to
    # write down at all.
    oligos = tuple(
        replace(
            record,
            number=position,
            id=_shared_id(uses[sequence]),
            uses=tuple(uses[sequence]),
            # A shared molecule belongs to every part it amplifies, so the
            # column agrees with the row's own id and ``used_as`` and can be
            # joined against the PCR table.
            anneal_ta=coolest_ta[sequence],
            part_id="+".join(sorted(parts[sequence])),
            part_ids=tuple(sorted(parts[sequence])),
            # One molecule can be a part's forward primer and another part's
            # reverse primer -- one part's end being the reverse complement of
            # another's start is all it takes.  Keeping the first encounter's
            # direction had the order sheet call it forward while the PCR table
            # listed it under ``reverse_primer_id``.
            direction="+".join(sorted(directions[sequence])),
            # A tailed molecule inherits its sweep's verdict; an untailed one
            # had no sweep and no tail, so nothing it breaches was avoidable.
            unavoidable=(
                tuple(sorted(tailed[sequence]))
                if sequence in tailed
                else tuple(sorted({v.constraint for v in record.violations}))
            ),
        )
        for position, (sequence, record) in enumerate(records.items(), start=1)
    )
    return oligos, {
        name: record.id for record in oligos for name in record.uses
    }


def _pcr_records(
    pairs: Iterable[OligoPair],
    checks: Iterable[PcrCheck],
    temperatures: Mapping[PcrOligos, PcrTemperatures],
    records: Mapping[str, Dseqrecord],
    oligo_of: Mapping[str, str],
    primers: Mapping[str, PrimerRecord],
) -> Tuple[PcrRecord, ...]:
    """One record per distinct reaction: one template, two oligos.

    Primers are referenced by the oligo that serves them, so two junctions that
    need the same two molecules on the same template are the same reaction and
    appear once.  The template is part of the key, because the same pair of
    oligos on a different template is a different reaction.

    The amplicon length is the part plus the two 5' tails: the footprints sit at
    the part's own ends, so everything a primer adds beyond its footprint is
    added to the product.
    """
    by_key = {(c.forward_id, c.reverse_id): c for c in checks}
    out: Dict[Tuple[str, str, str], PcrRecord] = {}
    for pair in pairs:
        forward_id = oligo_of[pair.forward.id]
        reverse_id = oligo_of[pair.reverse.id]
        key = (pair.part_id, forward_id, reverse_id)
        if key in out:
            continue
        oligos = pair.oligos
        check = by_key[(pair.forward.id, pair.reverse.id)]
        breaches = tuple(
            violation
            for primer_id in (forward_id, reverse_id)
            for violation in primers[primer_id].violations
        )
        out[key] = PcrRecord(
            part_id=pair.part_id,
            forward_id=forward_id,
            reverse_id=reverse_id,
            forward_label=primers[forward_id].label,
            reverse_label=primers[reverse_id].label,
            amplicon_bp=(
                len(records[pair.part_id])
                + len(oligos.forward_oligo)
                - len(oligos.forward_footprint)
                + len(oligos.reverse_oligo)
                - len(oligos.reverse_footprint)
            ),
            temperatures=temperatures[oligos],
            structures=check.structures,
            violations=check.violations,
            primer_violations=breaches,
            primer_severity=_worst(
                _flag_severity(primers[primer_id].violations)
                for primer_id in (forward_id, reverse_id)
            ),
            homodimer_dg_worst=min(
                primers[forward_id].structures.homodimer_dg,
                primers[reverse_id].structures.homodimer_dg,
            ),
            hairpin_folded_at_ta=max(
                primers[primer_id].structures.folded_at(
                    temperatures[oligos].ta_from_footprints.celsius
                )
                for primer_id in (forward_id, reverse_id)
            ),
            hairpin_3prime_paired=any(
                primers[primer_id].structures.hairpin_3prime_paired
                for primer_id in (forward_id, reverse_id)
            ),
        )
    return tuple(out.values())


def _relabel_primer_binds(
    construct: Optional[Dseqrecord], oligo_of: Mapping[str, str]
) -> Optional[Dseqrecord]:
    """Give the ``primer_bind`` features the name the order sheet uses.

    ``pydna.amplify.pcr`` labels each feature with the ``Primer.id`` of the
    junction it came from, which for a shared molecule is only one of the names
    that molecule serves.  Someone reading the GenBank has to be able to find
    that label on the order sheet, so the feature carries the oligo's name and
    the junction context is left to be read off the map.

    Only features this design's own PCRs created are touched.  A part sourced
    from a GenBank file brings its own annotations through the assembly, and one
    of its ``primer_bind`` features could carry a label that happens to equal a
    constructed primer name; renaming that would claim an oligo binds a site it
    does not and lose the user's annotation.  ``pcr`` marks its own features
    with a ``PCR_conditions`` qualifier, which is the discriminator, and a
    feature with no label of its own is left without one rather than given an
    empty list.
    """
    if construct is None:
        return None
    for feature in construct.features:
        if feature.type != "primer_bind" or "PCR_conditions" not in feature.qualifiers:
            continue
        labels = feature.qualifiers.get("label")
        if not labels:
            continue
        feature.qualifiers["label"] = [oligo_of.get(label, label) for label in labels]
    return construct


def _variant_plans(
    library: Sequence[Variant],
    oligos: Mapping[str, Sequence[OligoPair]],
    records: Mapping[str, Dseqrecord],
    spec: DesignSpec,
    oligo_of: Mapping[str, str],
    *,
    include_designs: bool,
) -> Tuple[Tuple[VariantPlan, ...], Tuple[Diagnostic, ...]]:
    """Verify every variant; keep the assembled record only if asked to.

    Verification is unconditional, so a result never reports a design it did not
    simulate; ``include_designs`` only decides whether the records survive the
    loop.
    """
    position_ids = tuple(position.id for position in spec.positions)
    slot_ids = tuple(slot_id for slot_id, _ in spec.ordered_slots())
    plans: List[VariantPlan] = []
    diagnostics: List[Diagnostic] = []
    for variant in library:
        by_slot = dict(zip(slot_ids, variant.parts))
        try:
            verification, construct = verify_variant(
                variant, oligos[variant.id], records, spec
            )
        except Exception as exc:  # noqa: BLE001 -- pydna raises freely here
            # Every other failure in plan() becomes a diagnostic, and an
            # assembly that pydna cannot carry out is a fact about the design,
            # not a bug in the caller.  pydna signals incompatible ends with a
            # bare TypeError out of ``Dseq.__add__``, so the type is no guide
            # and the class has to be broad to keep the contract that plan()
            # returns a result.
            verification, construct = (
                Verification(
                    verdict="assembly_failed",
                    message=f"pydna could not assemble this variant: {exc}",
                    n_products=0,
                    expected_bp=0,
                    observed_bp=None,
                ),
                None,
            )
        if not verification.ok:
            diagnostics.append(
                Diagnostic(
                    severity="error",
                    code=f"assembly_{verification.verdict}",
                    message=verification.message,
                    variant_id=variant.id,
                )
            )
        plans.append(
            VariantPlan(
                id=variant.id,
                choices=tuple(
                    (position_id, by_slot[position_id].id)
                    for position_id in position_ids
                ),
                # What is amplified, which is not every part of the construct: an
                # in_tail slot reaches the product through the primer tails.
                fragment_ids=tuple(
                    part.id for part in (variant.fragments or variant.parts)
                ),
                length_bp=verification.expected_bp,
                verification=verification,
                construct=(
                    _relabel_primer_binds(construct, oligo_of)
                    if include_designs
                    else None
                ),
            )
        )
    return tuple(plans), tuple(diagnostics)


# --------------------------------------------------------------------------- #
# The three tables
# --------------------------------------------------------------------------- #


def _project(rows: List[dict], columns: Sequence[str], full: bool) -> pd.DataFrame:
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    if not full:
        frame = frame[list(columns)]
    return frame.rename(columns=COLUMN_LABELS)


#: Columns the overview always has, which a position heading must not take.
_OVERVIEW_RESERVED = ("variant_id", "length_bp", "assembly")


def _overview_headings(spec: DesignSpec) -> Dict[str, str]:
    """A distinct column heading per position.

    ``Position.name`` is free text: nothing makes it unique or keeps it clear of
    the fixed columns.  Two positions sharing a name would collapse into one
    column and lose a position's choices, and a position named ``assembly``
    would be overwritten by the fixed key that follows it in the row literal.
    A name is used only when it is safe; otherwise the position's id is, which
    validation has already made unique.
    """
    names = Counter(
        position.name for position in spec.positions if position.name is not None
    )
    return {
        position.id: (
            position.name
            if position.name
            and names[position.name] == 1
            and position.name not in _OVERVIEW_RESERVED
            and position.name not in {p.id for p in spec.positions}
            else position.id
        )
        for position in spec.positions
    }


def variant_overview_dataframe(result: DesignResult) -> pd.DataFrame:
    """One row per variant, one column per position, headed by its name.

    The column headings come from the spec, so a renamed position renames the
    column and nothing has to be kept in step by hand.
    """
    if result.spec is None:
        return pd.DataFrame()
    headings = _overview_headings(result.spec)
    return pd.DataFrame(
        [
            {
                "variant_id": variant.id,
                **{
                    headings[position_id]: part_id
                    for position_id, part_id in variant.choices
                },
                "length_bp": variant.length_bp,
                "assembly": variant.verification.verdict,
            }
            for variant in result.variants
        ]
    )


def primer_dataframe(result: DesignResult, *, full: bool = False) -> pd.DataFrame:
    """The order sheet: one row per distinct **molecule**.

    An oligo several junctions share is one row, with ``n_uses`` above 1 and
    every name that asked for it in ``used_as``, so it is ordered and paid for
    once.  ``footprint_Tm (C)`` is the annealing region and ``oligo_Tm (C)`` the
    whole tailed molecule, both from ``tm_engine``; ΔG is kcal/mol and more
    positive is better.  ``not_analysable_over_60nt`` marks an oligo longer than
    primer3 will align: it is ordered exactly as designed, but its ΔG and
    structure Tm describe only its 3'-most 60 bases, so they bound the structure
    rather than read the whole molecule.  ``flagged_by`` names the bounds it
    breaches
    and ``flag_severity`` says whether to act: ``error`` when some other split
    avoids it, ``warning`` when none does.  ``full=True`` adds the structure Tm,
    the footprint sequence, the engine of each Tm and the breached bounds in
    full.
    """
    return _project(
        [primer.to_dict() for primer in result.primers], PRIMER_COLUMNS, full
    )


def pcr_dataframe(result: DesignResult, *, full: bool = False) -> pd.DataFrame:
    """One row per distinct reaction: one template, two molecules.

    ``Ta_from_footprints (C)`` is the annealing temperature of the **footprint**
    pair, which is what anneals in cycle 1; the tailed oligos' own Tm are in
    their own columns and are never what Ta is taken from.

    The ΔG columns are the *pair* --
    ``heterodimer_fwd_rev_deltaG (kcal/mol)`` and the 3' dimer -- plus
    ``worst_primer_homodimer_deltaG (kcal/mol)``, which is the worse of the two
    primers' own homodimers.  That last one is why a reaction holding a flagged
    oligo cannot read as unflagged here while the order sheet flags it.
    """
    return _project([pcr.to_dict() for pcr in result.pcrs], PCR_COLUMNS, full)
