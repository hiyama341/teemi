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

"""Primer design: the footprints, and where the overhang is split.

Two steps, in this order, and only the second one searches.

**Footprints.**  :func:`design_footprints` makes the one
``pydna.design.primer_design`` call in this package, with
:class:`~teemi.design.combinatorial.tm.FootprintTm` as ``tm_func`` and the local
model as ``estimate_function`` -- pydna's own two-stage walk, which halves the
NEB queries (6 per footprint against 12 on this repo's fixtures).  The
footprints are fixed from then on; nothing downstream recomputes, trims or
re-searches one.

**The split.**  At a junction between upstream A and downstream B with
``homology_length = L``, the forward primer of B carries A's last ``t_f`` bases
and the reverse primer of A carries the reverse complement of B's first ``t_r``,
with ``t_f + t_r = L``.  The two amplicons then overlap by exactly L and the
construct is A followed by B with nothing duplicated.  L means L: pydna's
``assembly_fragments`` applies ``ceil(overlap/2)`` to *both* sides, so its
``overlap=35`` yields 36 shared bases, which is not what the spec asked for.
Only the split moves -- :func:`design_junction` sweeps ``t_f`` and scores the
complete tailed primers.

Structures come from primer3 and nothing else, and they are scored at the scope
each one belongs to:

* **per primer** -- hairpin and self-dimer, on the complete tailed oligo, during
  the sweep (:class:`OligoStructures`);
* **per PCR** -- heterodimer and 3'-end stability, between the two primers that
  actually share a tube, once the junctions are chosen
  (:class:`PairStructures`, :func:`check_pcrs`).

A junction's own two primers are **not** such a pair.  They are perfectly
complementary over all L bases by construction -- that complementarity is the
overlap -- so scoring them against each other reads about -30 kcal/mol and
would reject every design ever written.  The forward primer's tube-mate is the
downstream part's *other* primer, which the next junction decides.

primer3 is never used for a Tm here.  Its thermodynamic Tm is a third model, and
mixing models is the defect class this package is built to avoid.
"""

from __future__ import annotations

import functools
import itertools
import math
from dataclasses import dataclass
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Literal,
    Mapping,
    Optional,
    Sequence,
    Tuple,
)

import primer3
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from pydna.amplify import pcr
from pydna.design import assembly_fragments, primer_design
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer

from teemi.design.combinatorial.spec import (
    ID_SEPARATOR,
    Constraints,
    DesignSpec,
    Pad,
    Part,
    Position,
    Weights,
)
from teemi.design.combinatorial.tm import (
    CONDITIONS,
    Conditions,
    FootprintTm,
    PcrOligos,
    Tm,
    require_same_engine,
)

__all__ = [
    "CLAMP_NT",
    "CAL_PER_KCAL",
    "PRIMER3_TEMP_C",
    "Severity",
    "Violation",
    "Diagnostic",
    "FootprintError",
    "JunctionError",
    "Footprints",
    "design_footprints",
    "OligoStructures",
    "PairStructures",
    "oligo_structures",
    "oligo_violations",
    "Candidate",
    "Junction",
    "JunctionDesign",
    "design_junction",
    "junctions_for_spec",
    "Variant",
    "variants",
    "adjacent_part_pairs",
    "OligoPair",
    "variant_oligos",
    "PcrCheck",
    "check_pcrs",
]

#: primer3's thermodynamic alignment raises ``RuntimeError`` above 60 nt, so
#: every sequence handed to it is cut to its 3'-most :data:`CLAMP_NT` bases --
#: the end that must stay free for the polymerase to extend from.  A primer long
#: enough to be cut already breaches ``primer_len_max``, whose own ceiling is 60,
#: so the approximation only ever applies to a candidate being rejected anyway.
CLAMP_NT = 60

#: primer3 reports dG and dH in cal/mol; every constraint is in kcal/mol.  An
#: off-by-1000 here leaves the -6.0 floor either inert or impossible and looks
#: plausible either way, so the conversion happens exactly once per primer3
#: call, at that boundary.
CAL_PER_KCAL = 1000.0

#: The temperature primer3 evaluates dG at.  primer3's own default, kept because
#: the -6.0 kcal/mol floors were set against numbers measured at it.
PRIMER3_TEMP_C = 37.0

Severity = Literal["info", "warning", "error", "critical"]

#: The bounds the split sweep is there to satisfy: secondary structure on the
#: complete oligo.
#:
#: Exhausting the sweep says something about these and nothing about a bound no
#: split could ever have changed -- an unorderable length, a GC fraction, a
#: footprint Tm outside its window.  Those stay errors however hard the search
#: tried, because the remedy is a different spec, not a different split.
STRUCTURE_BOUNDS = frozenset(
    {
        "hairpin_dg_min",
        "hairpin_tm_max",
        "homodimer_dg_min",
        "homodimer_tm_max",
        "heterodimer_dg_min",
        "heterodimer_tm_max",
        "three_prime_dimer_dg_min",
    }
)


@dataclass(frozen=True, slots=True)
class Violation:
    """A breached bound: the observation, the limit, and who measured both.

    ``measured_by`` is not optional.  An observation from one model checked
    against a limit calibrated for another is the defect class this package
    exists to avoid, so a breach cannot be recorded without naming the single
    source behind both sides of the comparison.
    """

    constraint: str
    subject: str
    bound: Literal["min", "max"]
    observed: float
    limit: float
    measured_by: str
    units: str

    @property
    def overshoot(self) -> float:
        """How far past the limit; positive for a real breach."""
        if self.bound == "min":
            return self.limit - self.observed
        return self.observed - self.limit

    def to_dict(self) -> dict:
        return {
            "constraint": self.constraint,
            "subject": self.subject,
            "bound": self.bound,
            "observed": self.observed,
            "limit": self.limit,
            "measured_by": self.measured_by,
            "units": self.units,
        }


@dataclass(frozen=True, slots=True)
class Diagnostic:
    """One flat record.  No taxonomy: severity, code, message, breach."""

    severity: Severity
    code: str
    message: str
    constraint: Optional[Violation] = None
    part_id: Optional[str] = None
    #: The oligo a breach belongs to, when it belongs to one.  A molecule
    #: several parts share has no single ``part_id``, and naming a part that
    #: does not exist would break grouping worse than leaving it unset.
    primer_id: Optional[str] = None
    junction_id: Optional[str] = None
    variant_id: Optional[str] = None

    @classmethod
    def from_violation(
        cls,
        violation: Violation,
        *,
        part_id: Optional[str] = None,
        junction_id: Optional[str] = None,
        primer_id: Optional[str] = None,
        severity: Severity = "error",
        exhausted: Optional[str] = None,
    ) -> "Diagnostic":
        """The one place a breach is turned into prose, so every message names
        the quantity, the bound and the model behind both.

        ``exhausted`` is the reason there was nothing better to pick, in the
        caller's own words, and is only ever passed with ``severity="warning"``.
        It is worth saying out loud because the two reasons are different: a
        sweep that walked every split and found no escape, or a terminal primer
        that is its own footprint and has no tail to move.
        """
        because = "" if exhausted is None else f"; {exhausted}"
        return cls(
            severity=severity,
            code=violation.constraint,
            message=(
                f"{part_id or junction_id or primer_id or 'design'}: "
                f"{violation.subject} "
                f"{violation.observed:.2f} {violation.units} breaches "
                f"{violation.constraint} {violation.limit:.2f} "
                f"(measured by {violation.measured_by}){because}"
            ),
            constraint=violation,
            part_id=part_id,
            junction_id=junction_id,
            primer_id=primer_id,
        )

    def to_dict(self) -> dict:
        return {
            "severity": self.severity,
            "code": self.code,
            "message": self.message,
            "constraint": (
                None if self.constraint is None else self.constraint.to_dict()
            ),
            "part_id": self.part_id,
            "primer_id": self.primer_id,
            "junction_id": self.junction_id,
            "variant_id": self.variant_id,
        }


class FootprintError(ValueError):
    """No primer can anneal to this part at all."""

    def __init__(self, message: str, *, part_id: str) -> None:
        super().__init__(f"part {part_id!r}: {message}")
        self.part_id = part_id


class JunctionError(ValueError):
    """The two parts cannot between them supply the requested homology."""

    def __init__(self, message: str, *, junction_id: str) -> None:
        super().__init__(f"junction {junction_id!r}: {message}")
        self.junction_id = junction_id


def _check_min(
    constraint: str,
    subject: str,
    observed: float,
    limit: float,
    measured_by: str,
    units: str,
) -> Optional[Violation]:
    if observed >= limit:
        return None
    return Violation(constraint, subject, "min", observed, limit, measured_by, units)


def _check_max(
    constraint: str,
    subject: str,
    observed: float,
    limit: float,
    measured_by: str,
    units: str,
) -> Optional[Violation]:
    if observed <= limit:
        return None
    return Violation(constraint, subject, "max", observed, limit, measured_by, units)


# --------------------------------------------------------------------------- #
# Footprints
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class Footprints:
    """The two annealing regions of one part, with the Tm of each.

    Immutable by contract, not merely by ``frozen``: these are what anneal in
    cycle 1 and what the annealing temperature is set from, so the junction
    search is allowed to move the split and nothing else.
    """

    part_id: str
    forward: str
    reverse: str
    forward_tm: Tm
    reverse_tm: Tm

    @property
    def pair_delta(self) -> float:
        """``|Tm_f - Tm_r|`` of this **footprint** pair.

        :meth:`~teemi.design.combinatorial.tm.Tm.delta` refuses two engines, so
        this cannot become the cross-model difference that used to be checked
        against ``pair_tm_delta_max``.
        """
        return self.forward_tm.delta(self.reverse_tm)

    def to_dict(self) -> dict:
        return {
            "part_id": self.part_id,
            "forward_footprint": self.forward,
            "reverse_footprint": self.reverse,
            "forward_footprint_len": len(self.forward),
            "reverse_footprint_len": len(self.reverse),
            "footprint_tm_forward": self.forward_tm.celsius,
            "footprint_tm_reverse": self.reverse_tm.celsius,
            "footprint_tm_engine": require_same_engine(
                self.forward_tm, self.reverse_tm
            ),
            "footprint_tm_delta": self.pair_delta,
        }


START_CODON = "ATG"

#: The three stop codons, so ``drop_stop`` knows what it is removing.
STOP_CODONS = ("TAA", "TAG", "TGA")


def apply_slot_flags(
    records: Mapping[str, Dseqrecord], spec: DesignSpec
) -> Tuple[Dict[str, Dseqrecord], Tuple[Diagnostic, ...]]:
    """Trim the codons the fusion does not keep, then orient each part.

    The codons first, in the part's own reading: a start codon is a fact about
    the part's frame, not about the construct.  Flipping first would look for
    the ATG at the 5' end of the oriented sequence, where a reversed CDS
    instead presents the reverse complement of its stop codon -- so the flag
    would miss its target and trim three bases of the wrong end if they
    happened to match.

    A ``drop_start`` or ``drop_stop`` whose codon is not there is reported and
    nothing is removed.  Silently trimming three bases that were not a codon
    would shift the reading frame of the fusion, which is worse than saying so.
    """
    by_part = {
        part.id: position
        for position in spec.positions
        for part in position.parts
    }
    out = dict(records)
    notes: List[Diagnostic] = []
    for part_id, position in by_part.items():
        record = out.get(part_id)
        if record is None:
            continue
        for flag, codons, where in (
            (position.drop_start, (START_CODON,), "start"),
            (position.drop_stop, STOP_CODONS, "end"),
        ):
            if not flag:
                continue
            sequence = str(record.seq).upper()
            present = (
                sequence[:3] in codons if where == "start" else sequence[-3:] in codons
            )
            if not present:
                notes.append(
                    Diagnostic(
                        severity="warning",
                        code=f"drop_{where}_codon_absent",
                        message=(
                            f"{part_id}: drop_"
                            f"{'start' if where == 'start' else 'stop'} was "
                            f"asked for but the part "
                            f"{'begins' if where == 'start' else 'ends'} "
                            f"{sequence[:3] if where == 'start' else sequence[-3:]!r}, "
                            f"not a {where} codon; nothing removed"
                        ),
                        part_id=part_id,
                    )
                )
                continue
            record = record[3:] if where == "start" else record[:-3]
            record.id = part_id
            record.name = part_id[:16]
        if position.reverse:
            flipped = record.reverse_complement()
            flipped.id = part_id
            flipped.name = part_id[:16]
            record = flipped
        out[part_id] = record
    return out, tuple(notes)


def design_footprints(
    records: Mapping[str, Dseqrecord],
    spec: DesignSpec,
    engine: FootprintTm,
    layout: Optional[Layout] = None,
) -> Tuple[Dict[str, Footprints], Tuple[Diagnostic, ...]]:
    """Footprints for every amplified part and pad, plus the bounds they breach.

    Only the amplified ones.  A part riding in its neighbours' primer tails is
    never a template, so no primer anneals to it and asking for a footprint
    would refuse a 9 bp Kozak sequence for being too short to do something it
    was never going to do.

    Raises :class:`FootprintError` for a part no primer can anneal to -- shorter
    than ``footprint_len_min``, or one pydna cannot design against.  That is
    invalid input rather than an infeasible design, so it is an exception and the
    caller turns it into a diagnostic at its own boundary.
    """
    amplified = {
        part.id
        for _, parts in (layout or Layout.from_spec(spec)).fragments
        for part in parts
    }
    footprints: Dict[str, Footprints] = {}
    diagnostics: List[Diagnostic] = []
    for part_id, record in records.items():
        if part_id not in amplified:
            continue
        footprint = _design_one(part_id, record, spec, engine)
        footprints[part_id] = footprint
        diagnostics.extend(_footprint_diagnostics(footprint, spec.constraints))
    return footprints, tuple(diagnostics)


def _design_one(
    part_id: str,
    record: Dseqrecord,
    spec: DesignSpec,
    engine: FootprintTm,
) -> Footprints:
    design = spec.primer_design
    if len(record) < design.footprint_len_min:
        raise FootprintError(
            f"{len(record)} bp is shorter than footprint_len_min "
            f"{design.footprint_len_min}, so there is nothing for a primer to "
            "anneal to",
            part_id=part_id,
        )
    try:
        amplicon = primer_design(
            record,
            limit=design.footprint_len_min,
            target_tm=design.target_tm,
            tm_func=engine.tm_func,
            estimate_function=engine.estimate_function,
        )
        seed_forward = str(amplicon.forward_primer.seq)
        seed_reverse = str(amplicon.reverse_primer.seq)
    except (AttributeError, ValueError, IndexError) as exc:
        raise FootprintError(
            f"pydna could not design a footprint at target_tm="
            f"{design.target_tm} C from {len(record)} bp ({exc})",
            part_id=part_id,
        ) from exc

    # pydna designs the forward footprint to ``target_tm`` and then re-targets
    # the reverse one at the forward's *measured* Tm, which is what keeps the
    # footprint pair balanced and ``pair_tm_delta_max`` satisfiable.  Both are
    # held to the same window edge here, and neither is touched again.
    forward, forward_tm = _shortest_footprint(
        str(record.seq), seed_forward, design.target_tm, engine
    )
    reverse, reverse_tm = _shortest_footprint(
        str(record.reverse_complement().seq), seed_reverse, design.target_tm, engine
    )
    # Both Tm are read back here, after pydna's walk, so they always come from
    # the same model: degradation is one-way and ``FootprintTm.tm`` stops
    # consulting NEB entirely once it has failed.  That is what makes
    # ``pair_delta`` always well defined.
    return Footprints(
        part_id=part_id,
        forward=forward,
        reverse=reverse,
        forward_tm=forward_tm,
        reverse_tm=reverse_tm,
    )


def _shortest_footprint(
    template: str, seed: str, target_tm: float, engine: FootprintTm
) -> Tuple[str, Tm]:
    """The shortest prefix of ``template``, from ``seed`` up, that reaches
    ``target_tm``.

    ``target_tm`` is the footprint window's lower edge -- the spec guarantees
    ``footprint_tm_window.min <= target_tm <= max`` -- so a footprint that
    reaches it is inside the window, and the window's own bounds stay as the
    hard check for the case this walk cannot satisfy: a part that runs out of
    template first.

    The walk exists because pydna keeps whichever of the last two lengths is
    *nearer* the target, which can leave the footprint one base short of it.
    Lengthening from there is the whole correction, and the lengths pydna
    already visited are cache hits.

    Reading the target as an edge rather than as a point is what makes this
    robust to the 3-4 C offset between the local model and NEB: clearing the
    same target by NEB needs more bases, and a point target plus a hard
    footprint-length cap is non-monotonic in ``target_tm`` where this is not.
    """
    length = len(seed)
    tm = engine.tm(template[:length])
    while tm.celsius < target_tm and length < len(template):
        length += 1
        tm = engine.tm(template[:length])
    return template[:length], tm


def _footprint_diagnostics(
    footprint: Footprints, constraints: Constraints
) -> Tuple[Diagnostic, ...]:
    """The footprint-pair bounds: the window on each, the delta on the pair.

    Both are footprint-engine bounds, and both sides of either comparison come
    from that one engine.  ``pair_tm_delta_max`` used to be checked on
    ``|local - NEB|``, which breached at 4.21 C on a pair whose footprints are
    0.88 C apart locally and 0.36 C apart at NEB.
    """
    window = constraints.footprint_tm_window
    engine = require_same_engine(footprint.forward_tm, footprint.reverse_tm)
    breaches = [
        _check_min(
            "footprint_tm_window.min", subject, tm.celsius, window.min, tm.engine, "C"
        )
        or _check_max(
            "footprint_tm_window.max", subject, tm.celsius, window.max, tm.engine, "C"
        )
        for subject, tm in (
            ("forward_footprint", footprint.forward_tm),
            ("reverse_footprint", footprint.reverse_tm),
        )
    ]
    breaches.append(
        _check_max(
            "pair_tm_delta_max",
            "footprint_pair",
            footprint.pair_delta,
            constraints.pair_tm_delta_max,
            engine,
            "C",
        )
    )
    return tuple(
        Diagnostic.from_violation(breach, part_id=footprint.part_id)
        for breach in breaches
        if breach is not None
    )


# --------------------------------------------------------------------------- #
# primer3 structures
# --------------------------------------------------------------------------- #


@functools.lru_cache(maxsize=1 << 16)
def _structure(
    calc: Callable[..., Any], seqs: Tuple[str, ...], conditions: Conditions
) -> Tuple[float, float]:
    """One primer3 structure: its dG in kcal/mol and its Tm in C.

    primer3 answers in cal/mol; the bounds are written in kcal/mol.  The
    conversion lives here and nowhere else, so the two scales cannot meet.
    More positive dG is better.
    """
    result = calc(*seqs, temp_c=PRIMER3_TEMP_C, **conditions.primer3_kwargs)
    return result.dg / CAL_PER_KCAL, result.tm


@dataclass(frozen=True, slots=True)
class OligoStructures:
    """primer3's hairpin and self-dimer numbers for one **complete** primer.

    Complete, tail included: a tail that pairs back with its own footprint is a
    hairpin the bench sees and the footprint alone does not.  dG is kcal/mol and
    more positive is better; the Tm are absolute primer3 structure
    temperatures, never a margin against Ta.
    """

    hairpin_dg: float
    hairpin_tm: float
    homodimer_dg: float
    homodimer_tm: float
    clamped: bool

    def to_dict(self, prefix: str = "") -> dict:
        return {
            f"{prefix}hairpin_dg": self.hairpin_dg,
            f"{prefix}hairpin_tm": self.hairpin_tm,
            f"{prefix}homodimer_dg": self.homodimer_dg,
            f"{prefix}homodimer_tm": self.homodimer_tm,
            f"{prefix}primer3_clamped": self.clamped,
        }


@dataclass(frozen=True, slots=True)
class PairStructures:
    """primer3's numbers for two primers that share a tube.

    Only ever a within-PCR pair.  ``three_prime_dimer_dg`` is the worse of the
    two directions, because a duplex occluding either primer's 3' end is the one
    a polymerase extends into a primer dimer, and it is scored under different
    alignment rules from the heterodimer.
    """

    heterodimer_dg: float
    heterodimer_tm: float
    three_prime_dimer_dg: float
    clamped: bool

    def to_dict(self) -> dict:
        return {
            "heterodimer_dg": self.heterodimer_dg,
            "heterodimer_tm": self.heterodimer_tm,
            "three_prime_dimer_dg": self.three_prime_dimer_dg,
            "primer3_clamped": self.clamped,
        }


def oligo_structures(seq: str, conditions: Conditions) -> OligoStructures:
    clamped = seq[-CLAMP_NT:]
    hairpin_dg, hairpin_tm = _structure(primer3.calc_hairpin, (clamped,), conditions)
    homodimer_dg, homodimer_tm = _structure(
        primer3.calc_homodimer, (clamped,), conditions
    )
    return OligoStructures(
        hairpin_dg=hairpin_dg,
        hairpin_tm=hairpin_tm,
        homodimer_dg=homodimer_dg,
        homodimer_tm=homodimer_tm,
        clamped=len(seq) > CLAMP_NT,
    )


def _pair_structures(
    forward: str, reverse: str, conditions: Conditions
) -> PairStructures:
    f, r = forward[-CLAMP_NT:], reverse[-CLAMP_NT:]
    heterodimer_dg, heterodimer_tm = _structure(
        primer3.calc_heterodimer, (f, r), conditions
    )
    return PairStructures(
        heterodimer_dg=heterodimer_dg,
        heterodimer_tm=heterodimer_tm,
        # Either primer's 3' end being occluded is what a polymerase extends.
        three_prime_dimer_dg=min(
            _structure(primer3.calc_end_stability, (f, r), conditions)[0],
            _structure(primer3.calc_end_stability, (r, f), conditions)[0],
        ),
        clamped=max(len(forward), len(reverse)) > CLAMP_NT,
    )


def oligo_violations(
    subject: str, seq: str, structures: OligoStructures, constraints: Constraints
) -> Tuple[Violation, ...]:
    """The hard bounds on one complete primer, whatever tube it ends up in.

    ``primer_len_max`` and ``gc_fraction`` are properties of the molecule that
    gets ordered, so they are read on the whole tailed primer.  The footprint-Tm
    bounds are deliberately absent -- footprints are fixed before any split is
    tried, so they are checked once per part instead of once per candidate -- and
    there is no whole-oligo pair-Tm bound at all: the tail does not anneal in
    cycle 1, and from cycle 2 a tail only raises Tm, so with Ta set from the
    footprints an imbalance there has no failure mode.
    """
    gc = gc_fraction(seq)
    checks = (
        _check_max(
            "primer_len_max",
            subject,
            float(len(seq)),
            float(constraints.primer_len_max),
            "sequence",
            "nt",
        ),
        _check_min(
            "gc_fraction_min",
            subject,
            gc,
            constraints.gc_fraction_min,
            "sequence",
            "fraction",
        ),
        _check_max(
            "gc_fraction_max",
            subject,
            gc,
            constraints.gc_fraction_max,
            "sequence",
            "fraction",
        ),
        _check_min(
            "hairpin_dg_min",
            subject,
            structures.hairpin_dg,
            constraints.hairpin_dg_min,
            "primer3",
            "kcal/mol",
        ),
        _check_max(
            "hairpin_tm_max",
            subject,
            structures.hairpin_tm,
            constraints.hairpin_tm_max,
            "primer3",
            "C",
        ),
        _check_min(
            "homodimer_dg_min",
            subject,
            structures.homodimer_dg,
            constraints.homodimer_dg_min,
            "primer3",
            "kcal/mol",
        ),
        _check_max(
            "homodimer_tm_max",
            subject,
            structures.homodimer_tm,
            constraints.homodimer_tm_max,
            "primer3",
            "C",
        ),
    )
    return tuple(check for check in checks if check is not None)


def _pair_violations(
    structures: PairStructures, constraints: Constraints
) -> Tuple[Violation, ...]:
    """The hard bounds on two primers that share a tube."""
    checks = (
        _check_min(
            "heterodimer_dg_min",
            "primer_pair",
            structures.heterodimer_dg,
            constraints.heterodimer_dg_min,
            "primer3",
            "kcal/mol",
        ),
        _check_max(
            "heterodimer_tm_max",
            "primer_pair",
            structures.heterodimer_tm,
            constraints.heterodimer_tm_max,
            "primer3",
            "C",
        ),
        _check_min(
            "three_prime_dimer_dg_min",
            "primer_pair",
            structures.three_prime_dimer_dg,
            constraints.three_prime_dimer_dg_min,
            "primer3",
            "kcal/mol",
        ),
    )
    return tuple(check for check in checks if check is not None)


def _cost(
    forward: str,
    reverse: str,
    forward_structures: OligoStructures,
    reverse_structures: OligoStructures,
    weights: Weights,
) -> float:
    """What ``accept="argmin"`` ranks candidates by; 0 is unimprovable.

    Only quantities the split actually moves and the split alone can see.  The
    footprints are fixed before the sweep starts and ``t_f + t_r = L`` fixes the
    two primers' total length, so what is left is each tailed primer's own
    structure energy and how unevenly the tails are shared.  Each dG term is the
    shortfall below 0 kcal/mol, so a candidate with no structure pays nothing.

    ``Weights.heterodimer_dg`` has no term here on purpose: a heterodimer needs
    two primers in one tube, and a junction's two primers are in two different
    tubes (see the module docstring).  It is checked, per PCR, in
    :func:`check_pcrs`.
    """

    def penalty(dg: float) -> float:
        return max(0.0, -dg)

    return (
        weights.hairpin_dg
        * (
            penalty(forward_structures.hairpin_dg)
            + penalty(reverse_structures.hairpin_dg)
        )
        + weights.homodimer_dg
        * (
            penalty(forward_structures.homodimer_dg)
            + penalty(reverse_structures.homodimer_dg)
        )
        + weights.length_imbalance * abs(len(forward) - len(reverse))
    )


# --------------------------------------------------------------------------- #
# The junction and its split sweep
# --------------------------------------------------------------------------- #


def _suffix(seq: str, n: int) -> str:
    """The last ``n`` bases, and ``""`` for ``n == 0``.

    ``seq[-0:]`` is the whole string in Python, which at the ``t_f = 0``
    boundary would paste an entire part into a primer.
    """
    return seq[len(seq) - n :] if n else ""


@dataclass(frozen=True, slots=True)
class Candidate:
    """One scored split of a junction's homology."""

    t_f: int
    t_r: int
    forward_primer: str
    reverse_primer: str
    forward_structures: OligoStructures
    reverse_structures: OligoStructures
    violations: Tuple[Violation, ...]
    cost: float

    @property
    def feasible(self) -> bool:
        return not self.violations

    @property
    def total_overshoot(self) -> float:
        """Sum of the breach magnitudes.

        Mixed units (nt, C, kcal/mol, fraction) on purpose: it exists only to
        order infeasible candidates so the one reported is the closest to
        working, and is never compared against a bound.
        """
        return sum(violation.overshoot for violation in self.violations)

    def to_dict(self) -> dict:
        return {
            "t_f": self.t_f,
            "t_r": self.t_r,
            "forward_primer": self.forward_primer,
            "reverse_primer": self.reverse_primer,
            "forward_primer_len": len(self.forward_primer),
            "reverse_primer_len": len(self.reverse_primer),
            "cost": self.cost,
            "feasible": self.feasible,
            "violations": [v.to_dict() for v in self.violations],
            **self.forward_structures.to_dict("forward_"),
            **self.reverse_structures.to_dict("reverse_"),
        }


@dataclass(frozen=True, slots=True)
class Junction:
    """One adjacency of the library: upstream part A, downstream part B.

    The homology is taken from the two parts' own sequences -- A's last ``t_f``
    bases followed by B's first ``t_r`` -- so the construct gains no base that
    is not in the design.
    """

    upstream: Footprints
    downstream: Footprints
    upstream_seq: str
    downstream_seq: str
    homology_length: int
    interposed: str = ""
    """Sequence that sits between the two parts without being amplified.

    An ``in_tail`` slot's part.  It is *inside* the homology, not added to it:
    ``homology_length`` stays the number of bases the two amplicons share, and
    what rides in the tails spends part of that budget, exactly as pydna's
    ``assembly_fragments`` accounts for a fragment shorter than its ``maxlink``.
    Asking for 40 bp of homology gives 40, whether or not a Kozak sequence sits
    in the middle of it."""

    forced_t_f: Optional[int] = None
    """The only split allowed here, when a ``no_tail`` flag pins it.

    ``no_tail_rev`` upstream leaves the downstream forward primer carrying all
    the homology (``t_f == homology_length``); ``no_tail_fwd`` downstream leaves
    it all on the upstream reverse primer (``t_f == 0``).  The sweep then has
    one candidate, which is why :func:`_unavoidable` declines to conclude
    anything from it."""

    @property
    def id(self) -> str:
        return f"{self.upstream.part_id}{ID_SEPARATOR}{self.downstream.part_id}"

    @property
    def budget(self) -> int:
        """How much of the homology the two parts themselves must supply.

        ``homology_length`` less whatever rides in the tails, so the split is
        searched over what is left rather than over the whole.
        """
        return self.homology_length - len(self.interposed)

    def primers(self, t_f: int) -> Tuple[str, str]:
        """``(forward primer of B, reverse primer of A)`` for this split.

        The tail goes on the 5' side of each primer and the footprint stays at
        the 3' end, which is the end that has to anneal and extend.
        """
        t_r = self.budget - t_f
        forward = (
            _suffix(self.upstream_seq, t_f)
            + self.interposed
            + self.downstream.forward
        )
        reverse = (
            str(
                Seq(self.interposed + self.downstream_seq[:t_r]).reverse_complement()
            )
            + self.upstream.reverse
        )
        return forward, reverse

    def homology(self, t_f: int) -> str:
        """The bases the two amplicons will share.

        ``homology_length`` of it comes from the two parts; anything riding in
        the tails sits in the middle and lengthens it.
        """
        return (
            _suffix(self.upstream_seq, t_f)
            + self.interposed
            + self.downstream_seq[: self.budget - t_f]
        )

    def split_range(self) -> Tuple[int, int]:
        """The ``t_f`` values the two parts are long enough to supply."""
        low = max(0, self.budget - len(self.downstream_seq))
        high = min(self.budget, len(self.upstream_seq))
        return low, high


@dataclass(frozen=True, slots=True)
class JunctionDesign:
    """The split chosen for one adjacency, and the profile it came from."""

    id: str
    upstream_id: str
    downstream_id: str
    homology_length: int
    homology: str
    chosen: Candidate
    #: ``(t_f, cost)`` for every candidate evaluated, in the order tried.  Under
    #: ``first_feasible`` that is only as far as the sweep got.
    profile: Tuple[Tuple[int, float], ...]
    #: What rides in the tails between the two parts, if anything.  It sits
    #: inside :attr:`homology`, which is why that can be longer than
    #: :attr:`homology_length`.
    interposed: str = ""
    #: Who built these primers.  ``pydna`` is
    #: ``pydna.design.assembly_fragments`` and the default; ``sweep`` is this
    #: module's own search, which ``tune_overhangs`` turns on; ``forced`` is one
    #: split this module had to build because a ``no_tail`` flag pins it and
    #: pydna has no way to express that.  Worth recording because the three
    #: mean different things by the split, and because a design that asked for
    #: pydna should be able to see where it did not get it.
    designed_by: Literal["pydna", "sweep", "forced"] = "sweep"
    #: ``(bound, subject)`` pairs that no candidate in the sweep avoided.
    #:
    #: Two things make this a statement about every legal split rather than a
    #: truncated walk: ``first_feasible`` only stops early on a *feasible*
    #: candidate, so whenever the chosen split still breaches, the whole range
    #: was evaluated; and :func:`_unavoidable` refuses to conclude anything from
    #: a sweep of fewer than two candidates, which is what ``fixed_balanced``
    #: produces.
    unavoidable: Tuple[Tuple[str, str], ...] = ()

    @property
    def feasible(self) -> bool:
        return self.chosen.feasible

    def to_dict(self) -> dict:
        return {
            "junction_id": self.id,
            "upstream_id": self.upstream_id,
            "downstream_id": self.downstream_id,
            "homology_length": self.homology_length,
            "interposed": self.interposed,
            "designed_by": self.designed_by,
            "unavoidable": [list(pair) for pair in self.unavoidable],
            "homology": self.homology,
            "profile": [list(point) for point in self.profile],
            "n_evaluated": len(self.profile),
            **self.chosen.to_dict(),
        }


def _split_order(order: str, homology_length: int) -> Tuple[int, ...]:
    """The ``t_f`` values to try, in the order to try them.

    ``outward_from_balanced`` starts at ``ceil(L/2)`` and steps +/-1 outward,
    which is the order that pays off with ``first_feasible``: the balanced split
    is the one that keeps both primers shortest.  ``exhaustive`` walks 0..L, so
    it is the one to pair with ``argmin``, which scores every candidate the
    order yields.  ``fixed_balanced`` tries the balanced split and nothing else.
    """
    balanced = math.ceil(homology_length / 2)
    if order == "fixed_balanced":
        return (balanced,)
    if order == "exhaustive":
        return tuple(range(homology_length + 1))
    return tuple(
        sorted(range(homology_length + 1), key=lambda t_f: (abs(t_f - balanced), t_f))
    )


def _shared_run(left: str, right: str, most: int) -> str:
    """The longest suffix of ``left`` that is also a prefix of ``right``.

    Measured rather than assumed, because on the pydna path the homology is
    whatever ``assembly_fragments`` built -- which is ``2 * ceil(overlap/2)``
    for an odd overlap, and includes its share of anything riding in the tails.
    """
    for length in range(min(most, len(left), len(right)), 0, -1):
        if left[-length:] == right[:length]:
            return left[-length:]
    return ""


def _pydna_candidate(
    junction: Junction, spec: DesignSpec, conditions: Conditions
) -> Tuple[Candidate, str]:
    """The one candidate ``pydna.design.assembly_fragments`` produces.

    pydna is handed this adjacency on its own -- two amplicons, with whatever
    rides in the tails between them -- which gives the same primers as handing
    it the whole construct, because the tails it builds depend only on the
    immediate neighbours.  So the module keeps designing each adjacency once
    and reusing it, and pydna decides the tails.

    The split is read back off the primers pydna returned rather than computed:
    it distributes a linker as half to each side, so the tail lengths are not
    something to predict.
    """
    upstream = pcr(
        Primer(junction.upstream.forward),
        Primer(junction.upstream.reverse),
        Dseqrecord(junction.upstream_seq),
    )
    downstream = pcr(
        Primer(junction.downstream.forward),
        Primer(junction.downstream.reverse),
        Dseqrecord(junction.downstream_seq),
    )
    fragments: List[object] = [upstream]
    if junction.interposed:
        fragments.append(Dseqrecord(junction.interposed))
    fragments.append(downstream)
    tailed = assembly_fragments(
        fragments,
        overlap=junction.homology_length,
        maxlink=spec.constraints.in_tail_max_bp,
    )
    forward = str(tailed[-1].forward_primer.seq)
    reverse = str(tailed[0].reverse_primer.seq)
    homology = _shared_run(
        str(tailed[0].seq),
        str(tailed[-1].seq),
        junction.homology_length + len(junction.interposed) + 8,
    )
    forward_structures = oligo_structures(forward, conditions)
    reverse_structures = oligo_structures(reverse, conditions)
    return Candidate(
        t_f=len(forward) - len(junction.downstream.forward),
        t_r=len(reverse) - len(junction.upstream.reverse),
        forward_primer=forward,
        reverse_primer=reverse,
        forward_structures=forward_structures,
        reverse_structures=reverse_structures,
        violations=oligo_violations(
            "forward_primer", forward, forward_structures, spec.constraints
        )
        + oligo_violations(
            "reverse_primer", reverse, reverse_structures, spec.constraints
        ),
        cost=_cost(
            forward,
            reverse,
            forward_structures,
            reverse_structures,
            spec.objective.weights,
        ),
    ), homology


def _unavoidable(
    candidates: Sequence[Candidate], chosen: Candidate
) -> Tuple[Tuple[str, str], ...]:
    """The ``(bound, subject)`` pairs the chosen split breaches that no split avoided.

    A breach every candidate shares is a property of the two parts, not of the
    choice the optimizer made, so re-running the search cannot help.  Only
    :data:`STRUCTURE_BOUNDS` are considered: a length or GC breach is not
    something a split was ever going to fix.

    Pairs rather than bare names, because a candidate's ``violations``
    concatenate the forward *and* the reverse primer's breaches.  Matching on
    the bound alone lets one that merely migrates from one primer to the other
    across the sweep count as unavoidable for both, and each would then be
    downgraded on the other's evidence.

    A sweep of one candidate proves nothing, so it claims nothing.
    ``fixed_balanced`` pins the split, and ``all()`` over its single candidate
    would otherwise "prove" every breach unavoidable from evidence that was
    never gathered -- when the remedy is precisely to let the optimizer search,
    which is what a blocking error should prompt.
    """
    if len(candidates) < 2:
        return ()
    breached = {
        (violation.constraint, violation.subject)
        for violation in chosen.violations
        if violation.constraint in STRUCTURE_BOUNDS
    }
    return tuple(
        sorted(
            pair
            for pair in breached
            if all(
                any(
                    (violation.constraint, violation.subject) == pair
                    for violation in candidate.violations
                )
                for candidate in candidates
            )
        )
    )


def design_junction(
    junction: Junction, spec: DesignSpec, conditions: Conditions = CONDITIONS
) -> JunctionDesign:
    """Sweep the split and keep one candidate.

    Never raises for an infeasible design and never returns ``None``: when no
    split passes every bound, the least-violating candidate comes back with the
    bounds it breached attached, because an unbuildable junction that names its
    problem is worth more than a missing one.

    :class:`JunctionError` is for a junction that cannot exist at all -- two
    parts too short between them to supply ``homology_length``.
    """
    low, high = junction.split_range()
    if low > high:
        raise JunctionError(
            f"{len(junction.upstream_seq)} bp and {len(junction.downstream_seq)} bp "
            f"cannot supply the {junction.budget} bp of homology left after "
            f"{len(junction.interposed)} bp riding in the tails"
            if junction.interposed
            else f"cannot supply homology_length {junction.homology_length}",
            junction_id=junction.id,
        )
    settings = spec.primer_design.junction
    search = settings.split_search
    if not settings.tune_overhangs and junction.forced_t_f is None:
        # pydna's own tails.  Nothing is searched, so there is one candidate and
        # the homology is whatever it built; both are read back off the result.
        chosen, homology = _pydna_candidate(junction, spec, conditions)
        return JunctionDesign(
            id=junction.id,
            upstream_id=junction.upstream.part_id,
            downstream_id=junction.downstream.part_id,
            homology_length=junction.homology_length,
            homology=homology,
            chosen=chosen,
            profile=((chosen.t_f, chosen.cost),),
            interposed=junction.interposed,
            designed_by="pydna",
        )
    if junction.forced_t_f is not None:
        if not low <= junction.forced_t_f <= high:
            raise JunctionError(
                f"a no_tail flag pins the split at t_f={junction.forced_t_f}, "
                f"which these two parts cannot supply (only {low}..{high} of "
                f"{junction.budget} bp is available)",
                junction_id=junction.id,
            )
        order: Tuple[int, ...] = (junction.forced_t_f,)
        built_by: Literal["pydna", "sweep", "forced"] = "forced"
    else:
        order = _split_order(search.order, junction.budget)
        built_by = "sweep"
    candidates: List[Candidate] = []
    for t_f in order:
        if not low <= t_f <= high:
            continue
        forward, reverse = junction.primers(t_f)
        forward_structures = oligo_structures(forward, conditions)
        reverse_structures = oligo_structures(reverse, conditions)
        candidates.append(
            Candidate(
                t_f=t_f,
                t_r=junction.budget - t_f,
                forward_primer=forward,
                reverse_primer=reverse,
                forward_structures=forward_structures,
                reverse_structures=reverse_structures,
                violations=oligo_violations(
                    "forward_primer", forward, forward_structures, spec.constraints
                )
                + oligo_violations(
                    "reverse_primer", reverse, reverse_structures, spec.constraints
                ),
                cost=_cost(
                    forward,
                    reverse,
                    forward_structures,
                    reverse_structures,
                    spec.objective.weights,
                ),
            )
        )
        if search.accept == "first_feasible" and candidates[-1].feasible:
            break

    feasible = [candidate for candidate in candidates if candidate.feasible]
    if feasible:
        chosen = (
            min(feasible, key=lambda candidate: candidate.cost)
            if search.accept == "argmin"
            else feasible[0]
        )
    else:
        chosen = min(
            candidates,
            key=lambda candidate: (
                len(candidate.violations),
                candidate.total_overshoot,
                candidate.cost,
            ),
        )
    return JunctionDesign(
        id=junction.id,
        upstream_id=junction.upstream.part_id,
        downstream_id=junction.downstream.part_id,
        homology_length=junction.homology_length,
        homology=junction.homology(chosen.t_f),
        interposed=junction.interposed,
        chosen=chosen,
        profile=tuple((candidate.t_f, candidate.cost) for candidate in candidates),
        designed_by=built_by,
        unavoidable=_unavoidable(candidates, chosen),
    )


# --------------------------------------------------------------------------- #
# The library: variants, and the junctions they share
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class Variant:
    """One member of the library, in assembly order.

    The id names only the position choices.  Pads are in :attr:`parts` because
    they are assembled, but they are the same in every variant and so
    distinguish nothing.

    :attr:`parts` is the whole construct, :attr:`fragments` only what is
    amplified.  They differ by the ``in_tail`` slots, which end up in the
    product via the primer tails and therefore belong to the first and not the
    second.
    """

    id: str
    parts: Tuple[Part, ...]
    fragments: Tuple[Part, ...] = ()


@dataclass(frozen=True, slots=True)
class Layout:
    """Which slots are amplified, and what rides in the tails between them.

    An ``in_tail`` slot is not a fragment: its sequence is built into the two
    flanking primers instead, so the construct still contains it while one
    fewer PCR is run.  That decision needs the sequences -- a part past
    ``in_tail_max_bp`` is amplified after all -- so it is made once, here, and
    everything downstream reads it rather than re-deriving it.
    """

    #: ``(slot_id, candidate parts, rides_in_tail)`` in assembly order.
    slots: Tuple[Tuple[str, Tuple[Part, ...], bool], ...]
    diagnostics: Tuple[Diagnostic, ...] = ()

    @classmethod
    def from_spec(cls, spec: DesignSpec) -> "Layout":
        """The layout the spec asks for, before the sequences are known.

        Honours every ``in_tail`` flag as written.  :func:`plan_layout` is the
        one that can also decline a request, because only it has the lengths.
        """
        asked = {position.id for position in spec.positions if position.in_tail}
        return cls(
            slots=tuple(
                (slot_id, parts, slot_id in asked)
                for slot_id, parts in spec.ordered_slots()
            )
        )

    @property
    def fragments(self) -> Tuple[Tuple[str, Tuple[Part, ...]], ...]:
        """The slots that become their own PCR."""
        return tuple(
            (slot_id, parts) for slot_id, parts, tail in self.slots if not tail
        )

    def interposed(self, records: Mapping[str, Dseqrecord]) -> Dict[str, str]:
        """Slot id -> the sequence it contributes to its neighbours' tails."""
        return {
            slot_id: str(records[parts[0].id].seq)
            for slot_id, parts, tail in self.slots
            if tail
        }


def plan_layout(
    spec: DesignSpec, records: Mapping[str, Dseqrecord]
) -> Layout:
    """Decide which ``in_tail`` requests the primer budget can actually honour."""
    asked = {position.id for position in spec.positions if position.in_tail}
    limit = spec.constraints.in_tail_max_bp
    slots: List[Tuple[str, Tuple[Part, ...], bool]] = []
    notes: List[Diagnostic] = []
    for slot_id, parts in spec.ordered_slots():
        if slot_id not in asked:
            slots.append((slot_id, parts, False))
            continue
        length = len(records[parts[0].id])
        if length > limit:
            notes.append(
                Diagnostic(
                    severity="warning",
                    code="in_tail_too_long",
                    message=(
                        f"{slot_id}: in_tail was asked for but the part is "
                        f"{length} bp, past in_tail_max_bp {limit}; it would be "
                        "added to both flanking primers, so it is amplified as "
                        "its own fragment instead"
                    ),
                    part_id=parts[0].id,
                )
            )
            slots.append((slot_id, parts, False))
            continue
        slots.append((slot_id, parts, True))
    return Layout(slots=tuple(slots), diagnostics=tuple(notes))


def variants(spec: DesignSpec, layout: Optional[Layout] = None) -> Tuple[Variant, ...]:
    """The library: the product of the positions, pads in place.

    :attr:`Variant.parts` is the whole construct in order, including whatever
    rides in the tails, because that is what the product has to come out equal
    to.  :attr:`Variant.fragments` is only what gets amplified.
    """
    layout = layout or Layout.from_spec(spec)
    fragment_slots = layout.fragments
    return tuple(
        _variant(spec, layout, combination)
        for combination in itertools.product(*(parts for _, parts in fragment_slots))
    )


def _variant(
    spec: DesignSpec, layout: Layout, chosen: Tuple[Part, ...]
) -> Variant:
    """Put one combination of fragments back into the full construct order."""
    picks = iter(chosen)
    order: List[Part] = []
    for _, parts, tail in layout.slots:
        order.append(parts[0] if tail else next(picks))
    return Variant(
        id=ID_SEPARATOR.join(
            part.id for part in order if not isinstance(part, Pad)
        ),
        parts=tuple(order),
        fragments=chosen,
    )


def adjacent_part_pairs(
    spec: DesignSpec, layout: Optional[Layout] = None
) -> Tuple[Tuple[str, str], ...]:
    """Every ``(upstream_id, downstream_id)`` that occurs in any variant, once.

    ``sum(n_i * n_i+1)`` pairs rather than ``prod(n_i)`` variants: a junction
    between two parts is the same junction whatever else the variant contains,
    so it is designed once and reused by every variant that contains it.

    Adjacency is between *fragments*: an ``in_tail`` slot does not break one, it
    sits inside it.
    """
    slots = (layout or Layout.from_spec(spec)).fragments
    adjacencies = list(zip(slots, slots[1:]))
    if spec.design.topology == "circular":
        adjacencies.append((slots[-1], slots[0]))
    pairs: Dict[Tuple[str, str], None] = {}
    for (_, upstream), (_, downstream) in adjacencies:
        for part_a in upstream:
            for part_b in downstream:
                pairs.setdefault((part_a.id, part_b.id), None)
    return tuple(pairs)


def junctions_for_spec(
    spec: DesignSpec,
    records: Mapping[str, Dseqrecord],
    footprints: Mapping[str, Footprints],
    layout: Optional[Layout] = None,
    conditions: Conditions = CONDITIONS,
) -> Dict[Tuple[str, str], JunctionDesign]:
    """Design each unique adjacency once, keyed by its part-id pair."""
    homology_length = spec.primer_design.junction.homology_length
    layout = layout or Layout.from_spec(spec)
    between = _interposed_between(spec, layout, records)
    slot_of = {
        part.id: position
        for position in spec.positions
        for part in position.parts
    }
    return {
        (upstream_id, downstream_id): design_junction(
            Junction(
                upstream=footprints[upstream_id],
                downstream=footprints[downstream_id],
                upstream_seq=str(records[upstream_id].seq),
                downstream_seq=str(records[downstream_id].seq),
                homology_length=homology_length,
                interposed=between.get((upstream_id, downstream_id), ""),
                forced_t_f=_forced_split(
                    slot_of.get(upstream_id),
                    slot_of.get(downstream_id),
                    homology_length
                    - len(between.get((upstream_id, downstream_id), "")),
                ),
            ),
            spec,
            conditions,
        )
        for upstream_id, downstream_id in adjacent_part_pairs(spec, layout)
    }


def _interposed_between(
    spec: DesignSpec, layout: Layout, records: Mapping[str, Dseqrecord]
) -> Dict[Tuple[str, str], str]:
    """What rides in the tails of each fragment adjacency.

    Consecutive ``in_tail`` slots concatenate in order, so two short elements
    between the same pair of fragments both end up in the product once.
    """
    riding = layout.interposed(records)
    order = list(layout.slots)
    if spec.design.topology == "circular":
        order = order + order[: len(order)]
    between: Dict[Tuple[str, str], str] = {}
    pending: List[str] = []
    previous: Optional[Tuple[str, Tuple[Part, ...]]] = None
    seen = 0
    for slot_id, parts, tail in order:
        if tail:
            pending.append(riding[slot_id])
            continue
        if previous is not None:
            for part_a in previous[1]:
                for part_b in parts:
                    between[(part_a.id, part_b.id)] = "".join(pending)
        previous = (slot_id, parts)
        pending = []
        seen += 1
        if seen > len(layout.fragments):
            break
    return between


def _forced_split(
    upstream: Optional[Position], downstream: Optional[Position], budget: int
) -> Optional[int]:
    """The split a ``no_tail`` flag pins this junction to, if either does.

    A pad carries no flags, hence the ``Optional`` on both sides.  Flags facing
    each other across one junction are refused in the spec, so at most one of
    these can fire.
    """
    if upstream is not None and upstream.no_tail_rev:
        return budget
    if downstream is not None and downstream.no_tail_fwd:
        return 0
    return None


@dataclass(frozen=True, slots=True)
class OligoPair:
    """One PCR of one variant: two named primers plus the sequences the NEB
    batch needs.

    The primers are :class:`pydna.primer.Primer` objects with real ids, so the
    ``primer_bind`` features ``pydna.amplify.pcr`` adds carry those ids instead
    of a placeholder label.
    """

    forward: Primer
    reverse: Primer
    oligos: PcrOligos

    @property
    def part_id(self) -> str:
        return self.oligos.part_id


def variant_oligos(
    spec: DesignSpec,
    variant: Variant,
    footprints: Mapping[str, Footprints],
    designs: Mapping[Tuple[str, str], JunctionDesign],
) -> Tuple[OligoPair, ...]:
    """The two primers each fragment of this variant is amplified with.

    A fragment's forward primer comes from the junction with its upstream
    neighbour and its reverse primer from the junction with its downstream
    neighbour, so a primer is shared by every variant with that same pair of
    neighbours.  In a linear design the first fragment's forward primer and the
    last fragment's reverse primer are the bare footprints, with no tail; in a
    circular one the last fragment joins the first and every primer is tailed.

    Only the fragments: an ``in_tail`` part has no PCR of its own, and its
    sequence reaches the product through the tails of the two primers either
    side of it.
    """
    parts = variant.fragments or variant.parts
    count = len(parts)
    circular = spec.design.topology == "circular"
    pairs = []
    for index, part in enumerate(parts):
        footprint = footprints[part.id]
        if circular or index > 0:
            upstream_id = parts[index - 1].id
            forward_id = f"{upstream_id}{ID_SEPARATOR}{part.id}_fw"
            forward_seq = designs[(upstream_id, part.id)].chosen.forward_primer
        else:
            forward_id = f"{part.id}_fw"
            forward_seq = footprint.forward
        if circular or index + 1 < count:
            downstream_id = parts[(index + 1) % count].id
            reverse_id = f"{part.id}{ID_SEPARATOR}{downstream_id}_rv"
            reverse_seq = designs[(part.id, downstream_id)].chosen.reverse_primer
        else:
            reverse_id = f"{part.id}_rv"
            reverse_seq = footprint.reverse
        pairs.append(
            OligoPair(
                forward=Primer(forward_seq, id=forward_id, name=forward_id),
                reverse=Primer(reverse_seq, id=reverse_id, name=reverse_id),
                oligos=PcrOligos(
                    part_id=part.id,
                    forward_oligo=forward_seq,
                    reverse_oligo=reverse_seq,
                    forward_footprint=footprint.forward,
                    reverse_footprint=footprint.reverse,
                ),
            )
        )
    return tuple(pairs)


@dataclass(frozen=True, slots=True)
class PcrCheck:
    """One PCR's within-tube numbers and the bounds they breach."""

    part_id: str
    forward_id: str
    reverse_id: str
    structures: PairStructures
    violations: Tuple[Violation, ...]

    def to_dict(self) -> dict:
        return {
            "part_id": self.part_id,
            "forward_primer_id": self.forward_id,
            "reverse_primer_id": self.reverse_id,
            "violations": [v.to_dict() for v in self.violations],
            **self.structures.to_dict(),
        }


def check_pcrs(
    pairs: Iterable[OligoPair],
    constraints: Constraints,
    conditions: Conditions = CONDITIONS,
) -> Tuple[PcrCheck, ...]:
    """Score the primer pairs that actually share a tube, once each.

    This is the only cross-reactivity scope there is: the two primers of one
    PCR.  A junction's two primers are complementary by construction and belong
    to two different reactions, so they are not checked against each other --
    and primers from different reactions never meet.

    Deduplicated by primer id, because one PCR is shared by every variant that
    contains the same part with the same two neighbours.
    """
    checks: Dict[Tuple[str, str], PcrCheck] = {}
    for pair in pairs:
        key = (pair.forward.id, pair.reverse.id)
        if key in checks:
            continue
        structures = _pair_structures(
            pair.oligos.forward_oligo, pair.oligos.reverse_oligo, conditions
        )
        checks[key] = PcrCheck(
            part_id=pair.part_id,
            forward_id=pair.forward.id,
            reverse_id=pair.reverse.id,
            structures=structures,
            violations=_pair_violations(structures, constraints),
        )
    return tuple(checks.values())
