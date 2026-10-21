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
    STRUCTURE_BOUNDS,
    Diagnostic,
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
    "primer_id",
    "part_id",
    "direction",
    "n_uses",
    "used_as",
    "length_nt",
    "sequence",
    "footprint_len",
    "footprint_tm_c",
    "oligo_tm_c",
    "footprint_tm_engine",
    "hairpin_dg",
    "homodimer_dg",
    "flagged_by",
    "flag_severity",
    "price",
)

#: The readable projection of a PCR row.
PCR_COLUMNS = (
    "part_id",
    "forward_primer_id",
    "reverse_primer_id",
    "amplicon_bp",
    "ta_from_footprints",
    "ta_engine",
    "footprint_tm_forward",
    "footprint_tm_reverse",
    "oligo_tm_forward",
    "oligo_tm_reverse",
    "footprint_tm_engine",
    "heterodimer_dg",
    "three_prime_dimer_dg",
    "homodimer_dg_worst",
    "flagged_by",
    "flag_severity",
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


def _flag_severity(
    violations: Sequence[Violation], unavoidable: Sequence[str]
) -> str:
    """Whether a flagged row needs acting on: ``error``, ``warning`` or ``""``.

    A row is only a warning when every bound it breaches is one the split sweep
    owns *and* one no candidate avoided.  A length or GC breach stays an error
    because no split was ever going to fix it, and a structure breach some
    other split did avoid stays an error because there is a choice to make.
    """
    if not violations:
        return ""
    return (
        "warning"
        if all(
            violation.constraint in STRUCTURE_BOUNDS
            and violation.constraint in unavoidable
            for violation in violations
        )
        else "error"
    )


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
    uses: Tuple[str, ...] = ()
    #: Breached bounds no candidate in the sweep that produced this molecule
    #: avoided.  For an untailed terminal primer there was no sweep and no
    #: tail, so every structure bound it breaches is in here.
    unavoidable: Tuple[str, ...] = ()
    #: Whether a junction split produced this molecule.  A tailed primer's
    #: breach is reported as a diagnostic against its junction, where the split
    #: that caused it lives; the row still shows it either way.
    tailed: bool = False

    def to_dict(self) -> dict:
        return {
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
            "dg_flagged": _dg_flagged(self.violations),
            "flagged_by": _violated(self.violations),
            "flag_severity": _flag_severity(self.violations, self.unavoidable),
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
    primer_violations: Tuple[Violation, ...] = ()
    homodimer_dg_worst: Optional[float] = None
    #: The two primers' unavoidable bounds, pooled, for scoring
    #: :attr:`primer_violations`.
    unavoidable: Tuple[str, ...] = ()

    def to_dict(self) -> dict:
        return {
            "part_id": self.part_id,
            "forward_primer_id": self.forward_id,
            "reverse_primer_id": self.reverse_id,
            "amplicon_bp": self.amplicon_bp,
            **self.temperatures.to_dict(),
            **self.structures.to_dict(),
            "homodimer_dg_worst": self.homodimer_dg_worst,
            "dg_flagged": _dg_flagged(self.violations + self.primer_violations),
            "flagged_by": _violated(self.violations + self.primer_violations),
            # The pair bounds are errors: the sweep optimises each junction on
            # its own and never searched this forward-against-reverse pairing,
            # so "nothing avoids it" has not been established for them.
            "flag_severity": (
                "error"
                if self.violations
                else _flag_severity(self.primer_violations, self.unavoidable)
            ),
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
                f"homology {self.spec.primer_design.junction.homology_length} bp",
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

    engine = FootprintTm(offline=offline)
    try:
        footprints, diagnostics = design_footprints(records, validated, engine)
        junctions = junctions_for_spec(validated, records, footprints)
    except FootprintError as exc:
        return _stopped(validated, "no_footprint", str(exc), part_id=exc.part_id)
    except JunctionError as exc:
        return _stopped(
            validated,
            "impossible_junction",
            str(exc),
            junction_id=exc.junction_id,
        )

    library = variants(validated)
    oligos = {
        variant.id: variant_oligos(validated, variant, footprints, junctions)
        for variant in library
    }
    pairs = [pair for variant_pairs in oligos.values() for pair in variant_pairs]
    checks = check_pcrs(pairs, validated.constraints)
    temperatures, batch_failures = neb_pair_temperatures(
        [pair.oligos for pair in pairs], offline=offline
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
                severity=(
                    "warning"
                    if violation.constraint in junction.unavoidable
                    else "error"
                ),
                tried=(
                    len(junction.profile)
                    if violation.constraint in junction.unavoidable
                    else None
                ),
            )
            for junction in junctions.values()
            for violation in junction.chosen.violations
        ]
        + [
            Diagnostic.from_violation(violation, part_id=check.part_id)
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
                part_id=primer.part_id,
                severity=(
                    "warning"
                    if violation.constraint in STRUCTURE_BOUNDS
                    else "error"
                ),
                tried=1 if violation.constraint in STRUCTURE_BOUNDS else None,
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
            tuple(diagnostics)
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
        for sequence in (
            junction.chosen.forward_primer,
            junction.chosen.reverse_primer,
        ):
            avoided = frozenset(junction.unavoidable)
            previous = tailed.get(sequence)
            tailed[sequence] = avoided if previous is None else previous & avoided
    records: Dict[str, PrimerRecord] = {}
    uses: Dict[str, List[str]] = {}
    parts: Dict[str, List[str]] = {}
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
            if sequence in records:
                if primer.id not in uses[sequence]:
                    uses[sequence].append(primer.id)
                if pair.part_id not in parts[sequence]:
                    parts[sequence].append(pair.part_id)
                continue
            uses[sequence] = [primer.id]
            parts[sequence] = [pair.part_id]
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
                ),
                price=len(sequence) * spec.objective.price_per_base,
                tailed=sequence in tailed,
            )
    oligos = tuple(
        replace(
            record,
            id=_shared_id(uses[sequence]),
            uses=tuple(uses[sequence]),
            # A shared molecule belongs to every part it amplifies, so the
            # column agrees with the row's own id and ``used_as`` and can be
            # joined against the PCR table.
            part_id="+".join(sorted(parts[sequence])),
            # A tailed molecule inherits its sweep's verdict; an untailed one
            # had no sweep and no tail, so nothing it breaches was avoidable.
            unavoidable=(
                tuple(sorted(tailed[sequence]))
                if sequence in tailed
                else tuple(sorted({v.constraint for v in record.violations}))
            ),
        )
        for sequence, record in records.items()
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
            unavoidable=tuple(
                sorted(
                    {
                        name
                        for primer_id in (forward_id, reverse_id)
                        for name in primers[primer_id].unavoidable
                    }
                )
            ),
            homodimer_dg_worst=min(
                primers[forward_id].structures.homodimer_dg,
                primers[reverse_id].structures.homodimer_dg,
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
    """
    if construct is None:
        return None
    for feature in construct.features:
        if feature.type != "primer_bind":
            continue
        labels = feature.qualifiers.get("label") or []
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
        verification, construct = verify_variant(
            variant, oligos[variant.id], records, spec
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
                fragment_ids=tuple(part.id for part in variant.parts),
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


def variant_overview_dataframe(result: DesignResult) -> pd.DataFrame:
    """One row per variant, one column per position, headed by its name.

    The column headings come from the spec, so a renamed position renames the
    column and nothing has to be kept in step by hand.
    """
    if result.spec is None:
        return pd.DataFrame()
    headings = {
        position.id: position.name or position.id for position in result.spec.positions
    }
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
    positive is better.  ``flagged_by`` names the bounds the molecule breaches
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
