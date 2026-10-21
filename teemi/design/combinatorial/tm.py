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

"""Melting and annealing temperatures, each tagged with the engine behind it.

Two engines, and never a number from one checked against a number from the
other:

* **local** -- :func:`pydna.tm.tm_default` (Biopython ``Tm_NN``) with this
  package's reaction conditions on its salt arguments.  Fast, offline, and the
  value every degraded path lands on.
* **NEB** -- the vendor's own model, which reads about 3-4 C below the local one
  for the same oligo.  Per sequence it is :func:`pydna.tm.tm_neb`, which is what
  pydna's two-stage ``primer_design`` walk calls; for the authoritative final
  numbers it is one batched POST to ``/tm/batch`` whose rows carry ``tm1``,
  ``tm2`` and ``ta`` together.

The 3-4 C offset is why :class:`Tm` carries its engine and has no arithmetic:
``abs(a - b)`` on two :class:`Tm` is a :class:`TypeError`, and :meth:`Tm.delta`
refuses two engines.  A cross-engine subtraction checked against a same-engine
bound shipped three times in the code this replaces; here it cannot be written.

:data:`PcrTemperatures` keeps the same separation on the output side: the
footprint and whole-oligo Tm live in differently named fields that nothing
mixes, and the annealing temperature is called ``ta_from_footprints`` because
that is the only pair it may come from.  NEB's recommendation for a 55-60 nt
tailed oligo runs 13-17 C above the footprints, and a PCR run at it fails.
"""

from __future__ import annotations

import functools
import math
from dataclasses import dataclass
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import requests
from pydna.tm import ta_dbd, tm_default, tm_neb

__all__ = [
    "NEB_BATCH_ENDPOINT",
    "NEB_PRODCODE",
    "NEB_TIMEOUT_S",
    "ENGINE_LOCAL_TM",
    "ENGINE_LOCAL_TA",
    "ENGINE_NEB",
    "Conditions",
    "CONDITIONS",
    "TmEngineMismatch",
    "Tm",
    "require_same_engine",
    "local_tm",
    "local_ta",
    "FootprintTm",
    "PcrOligos",
    "PcrTemperatures",
    "neb_pair_temperatures",
    "clear_tm_caches",
]

#: NEB's batch endpoint.  ``seqpairs`` takes a list of ``[seq]`` or
#: ``[seq1, seq2]``, so N pairs are one request rather than N.
NEB_BATCH_ENDPOINT = "https://tmapi.neb.com/tm/batch"

#: Q5, the polymerase this package's conditions describe.
NEB_PRODCODE = "q5-0"

NEB_TIMEOUT_S = 15.0

ENGINE_LOCAL_TM = "pydna_tm_default"
ENGINE_LOCAL_TA = "pydna_ta_dbd"
ENGINE_NEB = "neb"


@dataclass(frozen=True, slots=True)
class Conditions:
    """The reaction's physical conditions, in the units each engine wants.

    The defaults are :func:`pydna.tm.tm_default`'s own, so an offline run
    reproduces pydna's numbers exactly.  One object rather than three literal
    dicts because the monovalent total (Biopython's ``Na + K + Tris/2``) has to
    mean the same thing in ``tm_default``, in primer3 and in NEB's ``conc``;
    pydna's defaults put both strand concentrations at 250 nM.
    """

    na_mm: float = 40.0
    k_mm: float = 0.0
    tris_mm: float = 75.0
    mg_mm: float = 1.5
    dntps_mm: float = 0.8
    primer_nm: float = 250.0

    @property
    def monovalent_mm(self) -> float:
        return self.na_mm + self.k_mm + self.tris_mm / 2.0

    @property
    def neb_conc_um(self) -> float:
        """NEB's ``conc`` is in micromolar; ``primer_nm`` is in nanomolar."""
        return self.primer_nm / 1000.0

    @property
    def pydna_tm_kwargs(self) -> Dict[str, float]:
        return {
            "Na": self.na_mm,
            "K": self.k_mm,
            "Tris": self.tris_mm,
            "Mg": self.mg_mm,
            "dNTPs": self.dntps_mm,
            "dnac1": self.primer_nm,
            "dnac2": self.primer_nm,
        }

    @property
    def primer3_kwargs(self) -> Dict[str, float]:
        return {
            "mv_conc": self.monovalent_mm,
            "dv_conc": self.mg_mm,
            "dntp_conc": self.dntps_mm,
            "dna_conc": self.primer_nm,
        }


CONDITIONS = Conditions()


class TmEngineMismatch(TypeError):
    """Two temperatures from different models were about to be compared."""


@dataclass(frozen=True, slots=True)
class Tm:
    """A temperature in C, tagged with the model that produced it.

    Deliberately not a number: it has no ``__sub__`` and no ordering, so the
    only way to relate two of these is :meth:`delta`, which refuses two
    engines.  Read :attr:`celsius` to report or to check against a bound whose
    engine is known.
    """

    celsius: float
    engine: str

    def delta(self, other: "Tm") -> float:
        """``|self - other|``, or :class:`TmEngineMismatch` across engines."""
        require_same_engine(self, other)
        return abs(self.celsius - other.celsius)


def require_same_engine(*values: Tm) -> str:
    """The one engine behind every value, or :class:`TmEngineMismatch`.

    Called wherever a bound is checked, so a constraint breach can only be
    recorded against an engine that produced both sides of the comparison.
    """
    engines = {value.engine for value in values}
    if len(engines) != 1:
        raise TmEngineMismatch(
            "cannot compare temperatures from different models: "
            + ", ".join(sorted(engines))
        )
    return engines.pop()


# --------------------------------------------------------------------------- #
# Local models
# --------------------------------------------------------------------------- #


@functools.lru_cache(maxsize=1 << 16)
def _local_celsius(seq: str, conditions: Conditions) -> float:
    return float(tm_default(seq, **conditions.pydna_tm_kwargs))


def local_tm(seq: str, conditions: Conditions = CONDITIONS) -> Tm:
    """:func:`pydna.tm.tm_default` at the reaction's salts, memoized.

    Memoized because footprints recur across a combinatorial library: the same
    part's footprint is scored once per junction it takes part in.
    """
    return Tm(_local_celsius(seq, conditions), ENGINE_LOCAL_TM)


def local_ta(forward_footprint: str, reverse_footprint: str) -> Tm:
    """The local annealing temperature of a **footprint** pair.

    :func:`pydna.tm.ta_dbd` is ``min(tm_dbd(fp), tm_dbd(rp)) + 3`` and ignores
    its ``seq`` argument.  It is the DNA-binding-domain model, which is the one
    that matches :data:`NEB_PRODCODE`: Q5 carries a DBD, which raises the usable
    annealing temperature by about 4 C over the Taq-style ``ta_default``.
    """
    return Tm(float(ta_dbd(forward_footprint, reverse_footprint, "")), ENGINE_LOCAL_TA)


# --------------------------------------------------------------------------- #
# NEB, per sequence: the engine pydna's primer_design walk calls
# --------------------------------------------------------------------------- #


@functools.lru_cache(maxsize=1 << 14)
def _neb_celsius(seq: str, conc: float, prodcode: str) -> float:
    """One GET per distinct sequence, memoized; raises on anything else.

    ``lru_cache`` does not store exceptions, so a failure is not cached -- which
    is why :class:`FootprintTm` stops asking after the first one.
    """
    return float(tm_neb(seq, conc, prodcode))


def clear_tm_caches() -> None:
    """Drop both memos.

    The NEB memo holds answers from a remote service for the lifetime of the
    process, so a long-running caller and a test both need a way to forget them.
    """
    _local_celsius.cache_clear()
    _neb_celsius.cache_clear()


class FootprintTm:
    """The Tm engine pydna's ``primer_design`` walk is driven by.

    ``offline=True`` issues no request at all.  Online, each distinct sequence
    costs one :func:`pydna.tm.tm_neb` GET, and the first unreachable call
    records a message and leaves the engine local for the rest of the run:
    ``tm_neb`` passes no ``timeout`` to ``requests``, so a dead endpoint would
    otherwise cost one hang per sequence rather than one in total.

    Every value comes out tagged with the model that actually produced it, so a
    run that degraded half way through says so value by value.
    """

    def __init__(
        self,
        *,
        offline: bool,
        conditions: Conditions = CONDITIONS,
        prodcode: str = NEB_PRODCODE,
    ) -> None:
        self.offline = bool(offline)
        self.conditions = conditions
        self.prodcode = prodcode
        self._local_only = self.offline
        self._failures: list[str] = []
        #: pydna's cheap first guess, which halves the NEB queries: measured on
        #: this repo's fixtures at 6 per footprint against 12 without it.
        #: ``None`` offline, where ``tm_func`` is already the cheap model and a
        #: second pass would only double the work.
        self.estimate_function = None if self.offline else self._local_tm_func

    @property
    def engine(self) -> str:
        """The model the next value will carry."""
        return ENGINE_LOCAL_TM if self._local_only else ENGINE_NEB

    @property
    def failures(self) -> Tuple[str, ...]:
        """Messages for the caller to raise as diagnostics."""
        return tuple(self._failures)

    def tm(self, seq: str) -> Tm:
        """Once degraded, local for everything -- including sequences NEB had
        already answered.  That keeps every footprint Tm a degraded run reports
        on one model, so a pair delta is always well defined."""
        if self._local_only:
            return local_tm(seq, self.conditions)
        try:
            celsius = _neb_celsius(seq, self.conditions.neb_conc_um, self.prodcode)
        except (requests.RequestException, ValueError, KeyError) as exc:
            self._local_only = True
            self._failures.append(
                f"pydna.tm.tm_neb could not be reached ({exc}); footprint Tm "
                f"came from the local {ENGINE_LOCAL_TM} model instead, which "
                "reads about 3-4 C above NEB for the same oligo"
            )
            return local_tm(seq, self.conditions)
        return Tm(celsius, ENGINE_NEB)

    def tm_func(self, seq: str) -> float:
        """``(str) -> float``, for ``primer_design(tm_func=...)``.

        The bare float is pydna's contract.  It is consumed by pydna's walk,
        which returns a *sequence*, not a number; the chosen footprint is then
        re-read through :meth:`tm`, so an untagged float never reaches a
        comparison.
        """
        return self.tm(seq).celsius

    def _local_tm_func(self, seq: str) -> float:
        return local_tm(seq, self.conditions).celsius


# --------------------------------------------------------------------------- #
# NEB, batched: the authoritative numbers
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class PcrOligos:
    """The sequences of one PCR, as the batch request needs them.

    Four named fields rather than a tuple: ``(oligo, footprint)`` the wrong way
    round is a silent 13-17 C error in the Ta column, not a type error.
    """

    part_id: str
    forward_oligo: str
    reverse_oligo: str
    forward_footprint: str
    reverse_footprint: str


@dataclass(frozen=True, slots=True)
class PcrTemperatures:
    """What one PCR's rows said, kept apart by what they describe.

    The oligo and footprint Tm are separate fields that nothing averages or
    subtracts, and the annealing temperature is named for the pair it is
    allowed to come from.
    """

    oligo_tm_forward: Tm
    oligo_tm_reverse: Tm
    footprint_tm_forward: Tm
    footprint_tm_reverse: Tm
    ta_from_footprints: Tm

    def to_dict(self) -> Dict[str, Any]:
        """The report's column names, fixed here rather than at each caller.

        The two Tm quantities never share a column, each carries the model that
        produced it, and the Ta column says which pair it came from -- printing
        a footprint Tm and an oligo Tm as one quantity, and labelling the
        oligos' Ta as footprint-derived, are two of the three bugs this module
        was rewritten to make unwriteable.
        """
        return {
            "oligo_tm_forward": self.oligo_tm_forward.celsius,
            "oligo_tm_reverse": self.oligo_tm_reverse.celsius,
            "oligo_tm_engine": require_same_engine(
                self.oligo_tm_forward, self.oligo_tm_reverse
            ),
            "footprint_tm_forward": self.footprint_tm_forward.celsius,
            "footprint_tm_reverse": self.footprint_tm_reverse.celsius,
            "footprint_tm_engine": require_same_engine(
                self.footprint_tm_forward, self.footprint_tm_reverse
            ),
            "ta_from_footprints": self.ta_from_footprints.celsius,
            "ta_engine": self.ta_from_footprints.engine,
        }


def neb_pair_temperatures(
    pcrs: Sequence[PcrOligos],
    *,
    offline: bool,
    conditions: Conditions = CONDITIONS,
    prodcode: str = NEB_PRODCODE,
    timeout: float = NEB_TIMEOUT_S,
) -> Tuple[Dict[PcrOligos, PcrTemperatures], Tuple[str, ...]]:
    """The authoritative Tm and Ta for every PCR, in **one** POST.

    Both row sets travel in the same request -- the tailed-oligo pairs first,
    then the footprint pairs -- so ``tm1``/``tm2`` are read off the oligo rows
    and ``ta`` off the **footprint** rows at no extra cost.  Ta is the
    temperature the footprints anneal at; NEB's recommendation for the tailed
    oligos is 13-17 C higher and would make those reactions fail.

    ``offline=True`` issues no request.  A failed request, or a reply missing a
    number, degrades the affected PCRs to the local models and returns a message
    per failure for the caller's diagnostics; it never raises.
    """
    unique = tuple(dict.fromkeys(pcrs))
    if not unique:
        return {}, ()
    if offline:
        return {p: _local_temperatures(p, conditions) for p in unique}, ()

    payload = {
        "seqpairs": [[p.forward_oligo, p.reverse_oligo] for p in unique]
        + [[p.forward_footprint, p.reverse_footprint] for p in unique],
        "conc": conditions.neb_conc_um,
        "prodcode": prodcode,
    }
    rows, failure = _post_batch(payload, timeout)
    if failure is not None:
        return {p: _local_temperatures(p, conditions) for p in unique}, (failure,)

    results: Dict[PcrOligos, PcrTemperatures] = {}
    failures: list[str] = []
    count = len(unique)
    for index, pcr in enumerate(unique):
        oligo_row = rows[index] if index < len(rows) else None
        footprint_row = rows[count + index] if count + index < len(rows) else None
        oligo_forward = _number(oligo_row, "tm1")
        oligo_reverse = _number(oligo_row, "tm2")
        footprint_forward = _number(footprint_row, "tm1")
        footprint_reverse = _number(footprint_row, "tm2")
        ta = _number(footprint_row, "ta")
        if None in (
            oligo_forward,
            oligo_reverse,
            footprint_forward,
            footprint_reverse,
            ta,
        ):
            results[pcr] = _local_temperatures(pcr, conditions)
            failures.append(
                f"NEB returned no usable tm1/tm2/ta for part {pcr.part_id!r}; "
                "its Tm and Ta came from the local pydna models instead"
            )
            continue
        results[pcr] = PcrTemperatures(
            oligo_tm_forward=Tm(oligo_forward, ENGINE_NEB),
            oligo_tm_reverse=Tm(oligo_reverse, ENGINE_NEB),
            footprint_tm_forward=Tm(footprint_forward, ENGINE_NEB),
            footprint_tm_reverse=Tm(footprint_reverse, ENGINE_NEB),
            ta_from_footprints=Tm(ta, ENGINE_NEB),
        )
    return results, tuple(failures)


def _local_temperatures(pcr: PcrOligos, conditions: Conditions) -> PcrTemperatures:
    return PcrTemperatures(
        oligo_tm_forward=local_tm(pcr.forward_oligo, conditions),
        oligo_tm_reverse=local_tm(pcr.reverse_oligo, conditions),
        footprint_tm_forward=local_tm(pcr.forward_footprint, conditions),
        footprint_tm_reverse=local_tm(pcr.reverse_footprint, conditions),
        ta_from_footprints=local_ta(pcr.forward_footprint, pcr.reverse_footprint),
    )


def _post_batch(
    payload: Mapping[str, Any], timeout: float
) -> Tuple[list, Optional[str]]:
    degraded = "Tm and Ta came from the local pydna models instead"
    try:
        response = requests.post(NEB_BATCH_ENDPOINT, json=payload, timeout=timeout)
        response.raise_for_status()
        body = response.json()
    except (requests.RequestException, ValueError) as exc:
        return [], f"{NEB_BATCH_ENDPOINT} could not be reached ({exc}); {degraded}"
    if not isinstance(body, Mapping) or not body.get("success"):
        reason = body.get("error") if isinstance(body, Mapping) else "no JSON object"
        return [], f"{NEB_BATCH_ENDPOINT} refused the batch ({reason}); {degraded}"
    rows = body.get("data")
    if not isinstance(rows, list):
        return [], f"{NEB_BATCH_ENDPOINT} replied without data rows; {degraded}"
    return rows, None


def _number(row: Any, key: str) -> Optional[float]:
    """``row[key]`` as a usable temperature, or ``None``.

    A batch HTTP API has many ways to say "not this row" -- a missing key, a
    ``null``, a string placeholder, a per-row error object.  Each becomes
    ``None`` so the caller takes its documented local path instead of raising
    out of the middle of a design.  ``bool`` is excluded because ``True`` is an
    ``int`` and 1.0 C is not a temperature.
    """
    if not isinstance(row, Mapping):
        return None
    raw = row.get(key)
    if isinstance(raw, bool) or not isinstance(raw, (int, float)):
        return None
    value = float(raw)
    return value if math.isfinite(value) else None
