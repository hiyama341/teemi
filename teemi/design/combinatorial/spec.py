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

"""The JSON design specification for :mod:`teemi.design.combinatorial`.

One frozen pydantic tree rooted at :class:`DesignSpec`; there is no schema
version and no migration path.  Every collection is a tuple and no field is a
dict, so a spec is hashable and can key a cache.

A spec carries no URL and no network target of any kind: that is what keeps the
input boundary in :mod:`teemi.design.combinatorial.sources` a path check rather
than a filter to be bypassed.  Sequence content is *not* validated here either
-- alphabet, length and count caps all live in that module, so they are applied
exactly once.

Every bound in :class:`Constraints` names the engine whose number it bounds:
``footprint_*`` bounds come from the footprint Tm engine, the structure ``*_dg``
and ``*_tm`` bounds from primer3.  Nothing in this schema can express a bound on
one engine's number checked against another's, which is the defect class this
package exists to avoid.
"""

from __future__ import annotations

import math
from typing import Annotated, Literal, Mapping, Optional, Tuple, Union

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    StringConstraints,
    model_validator,
)

__all__ = [
    "ID_SEPARATOR",
    "MAX_PARTS_TOTAL",
    "MAX_VARIANTS",
    "Slug",
    "InlineSource",
    "FileSource",
    "PartSource",
    "Part",
    "Pad",
    "Position",
    "Design",
    "SplitSearch",
    "Junction",
    "PrimerDesign",
    "TmWindow",
    "Constraints",
    "Weights",
    "Objective",
    "DesignSpec",
    "load_spec",
    "json_schema",
]

ID_SEPARATOR = "__"
"""Separator for composite ids (``"pMLS1__AKF02530"``).

Composite ids become JSON object keys and filename stems, so they must be
strings rather than tuples and must survive a shell glob.
"""

MAX_PARTS_TOTAL = 512
MAX_VARIANTS = 10_000

Slug = Annotated[str, StringConstraints(pattern=r"^[A-Za-z0-9][A-Za-z0-9_.\-]{0,63}$")]
"""An identifier safe in a report column, a filename and a URL path."""


class _Frozen(BaseModel):
    """``frozen`` makes a spec hashable; ``forbid`` turns a typo'd key into an
    error instead of a design silently computed against the default."""

    model_config = ConfigDict(frozen=True, extra="forbid")


class InlineSource(_Frozen):
    type: Literal["inline"] = "inline"
    sequence: str = Field(min_length=1)
    topology: Literal["linear", "circular"] = "linear"


class FileSource(_Frozen):
    """A sequence file under the resolver's ``base_dir``.

    :attr:`format` is declared rather than sniffed, and :attr:`path` is
    interpreted only by :class:`~teemi.design.combinatorial.sources.SourceResolver`,
    which refuses anything that escapes its ``base_dir``.
    """

    type: Literal["file"] = "file"
    path: str = Field(min_length=1, max_length=512)
    format: Literal["fasta", "genbank"]
    record_id: Optional[str] = Field(default=None, max_length=256)
    topology: Literal["linear", "circular"] = "linear"


PartSource = Annotated[Union[FileSource, InlineSource], Field(discriminator="type")]


class Part(_Frozen):
    id: Slug
    name: Optional[str] = Field(default=None, max_length=200)
    source: PartSource


class Pad(Part):
    """A fixed fragment spliced in next to a position (a linker, a homology arm).

    A pad anchors to a position *id* and a side, never to an index: the legacy
    ``DesignAssembly`` took list indices and silently changed meaning when a
    position was added.  Exactly one of the two anchors is set.
    """

    before_position: Optional[Slug] = None
    after_position: Optional[Slug] = None

    @model_validator(mode="after")
    def _one_anchor(self) -> "Pad":
        if (self.before_position is None) == (self.after_position is None):
            raise ValueError(
                f"pad {self.id!r} needs exactly one of before_position/after_position"
            )
        return self


class Position(_Frozen):
    """One slot of the library; the variants are the product of the positions.

    The six flags are *slot-level* facts, following FLOW's ``construct_part``:
    they describe the place on the DNA, not the sequence that fills it, so
    several alternatives at one position cannot disagree about which strand the
    position sits on or whether its start codon survives the fusion.
    """

    id: Slug
    name: Optional[str] = Field(default=None, max_length=200)
    parts: Tuple[Part, ...] = Field(min_length=1)

    reverse: bool = False
    """This part sits on the other strand; the reverse complement is taken here
    so everything downstream reads one orientation."""

    in_tail: bool = False
    """Build this part into the flanking primer tails instead of amplifying it.
    A request, not a command: a part too long to fit the primer budget comes
    back as its own fragment with a warning saying so."""

    no_tail_fwd: bool = False
    """Omit the homology tail on this part's forward primer, which is the one
    giving homology at its 5' end.  The junction then takes all its homology
    from the upstream part's reverse primer."""

    no_tail_rev: bool = False
    """Omit the homology tail on this part's reverse primer, which lands on its
    3' end.  The junction takes all its homology from the downstream part's
    forward primer."""

    drop_start: bool = False
    """Drop this part's own start codon, as a CDS downstream of a signal
    peptide does to keep one ATG across the fusion."""

    drop_stop: bool = False
    """Drop this part's own stop codon, as a CDS fused to a downstream tag or
    linker does."""


class Design(_Frozen):
    """:attr:`topology` is load-bearing: ``circular`` makes the last fragment
    join the first, and verification then compares rotation-invariantly."""

    name: str = Field(min_length=1, max_length=200)
    topology: Literal["linear", "circular"] = "linear"
    description: Optional[str] = Field(default=None, max_length=4000)


class SplitSearch(_Frozen):
    """How the homology length ``L`` is split into the two tails ``t_f + t_r = L``.

    :attr:`order` enumerates candidate splits, :attr:`accept` decides when to
    stop: ``outward_from_balanced`` starts at ``ceil(L/2)`` and steps +/-1
    outward, ``exhaustive`` tries all ``L + 1`` splits, ``fixed_balanced`` tries
    only ``ceil(L/2)``.  ``first_feasible`` takes the first split that passes
    every constraint; ``argmin`` takes the cheapest under
    :class:`Weights`.
    """

    order: Literal["outward_from_balanced", "exhaustive", "fixed_balanced"] = (
        "outward_from_balanced"
    )
    accept: Literal["first_feasible", "argmin"] = "first_feasible"


class Junction(_Frozen):
    """:attr:`homology_length` is **exactly** ``L`` shared bases.

    pydna's ``assembly_fragments`` applies ``ceil(overlap/2)`` to both sides, so
    its ``overlap=35`` yields 36 bases of homology; here 35 means 35.
    """

    homology_length: int = Field(default=35, ge=12, le=80)

    tune_overhangs: bool = False
    """Search the junction's split, instead of taking pydna's.

    Off by default, and off means pydna: ``pydna.design.assembly_fragments``
    builds the tails, because it is the maintained and validated
    implementation of this manoeuvre and there is no reason to prefer a local
    one.  It offers a single balanced split (``ceil(overlap/2)`` each side),
    which is why turning this on is what buys a search at all.

    FLOW sends this under the same name, and its own description says why it is
    safe to hand to a search: 18/22 and 20/20 leave the two amplicons sharing
    the same 40 bp, so the assembled product is identical either way and only
    the two primer sequences change."""

    split_search: SplitSearch = SplitSearch()
    """How to search, when :attr:`tune_overhangs` says to."""


class PrimerDesign(_Frozen):
    """Footprint design, which happens once per part and is then immutable.

    :attr:`target_tm` and :attr:`footprint_len_min` are pydna's
    ``primer_design(target_tm=..., limit=...)``: pydna extends the footprint from
    ``limit`` until its Tm reaches the target, so the target is the *lower* edge
    of :attr:`Constraints.footprint_tm_window` and the footprint taken is the
    shortest one inside that window.  A window rather than a point target is
    what makes this robust to the 3-4 C offset between the local estimator and
    NEB; a single target plus a hard footprint length cap is not.
    """

    target_tm: float = Field(default=57.0, ge=30.0, le=80.0)
    footprint_len_min: int = Field(default=18, ge=10, le=60)

    polymerase: str = "q5-0"
    """NEB's product code, which selects the Tm and Ta model.

    FLOW sends it under ``assembly.polymerase`` and says choosing the model is
    the backend's job, so this is where that choice lands: it is the ``prodcode``
    the NEB batch is asked for, and NEB returns the annealing temperature its
    own rule gives for that polymerase.  ``q5-0`` and ``phusion-0`` differ by
    2 C in Ta, so it is not cosmetic.

    Offline it changes nothing, because pydna implements only the Phusion rule
    -- see :func:`~teemi.design.combinatorial.tm.local_ta`."""

    junction: Junction = Junction()


class TmWindow(_Frozen):
    min: float = Field(ge=0.0, le=120.0)
    max: float = Field(ge=0.0, le=120.0)

    @model_validator(mode="after")
    def _ordered(self) -> "TmWindow":
        if self.min >= self.max:
            raise ValueError(f"Tm window min {self.min} must be below max {self.max}")
        return self


class Constraints(_Frozen):
    """Hard bounds.  A design that breaches one is reported, never shipped.

    :attr:`footprint_tm_window` and :attr:`pair_tm_delta_max` are about the
    **footprint pair** -- the footprints are what anneal in cycle 1 and what the
    annealing temperature is set from -- and both sides of either comparison come
    from the same Tm engine.  There is deliberately no bound on the whole-oligo
    pair Tm: the tail does not anneal in cycle 1, and from cycle 2 a tail only
    raises Tm, so with Ta set from the footprint an imbalance there has no
    failure mode.

    The ``*_dg`` floors and ``*_tm`` ceilings are primer3 numbers for the
    **complete** primer including its tail, in kcal/mol (primer3 reports ΔG in
    cal/mol; the conversion happens once, at that boundary).  More positive ΔG is
    better.  :attr:`three_prime_dimer_dg_min` is primer3's ``calc_end_stability``
    and is not nested with :attr:`heterodimer_dg_min`: a duplex that occludes the
    3' end is the one a polymerase extends into a primer dimer, and the two are
    scored under different alignment rules.

    :attr:`primer_len_max` is the **orderable** cap and nothing else.  It used
    to double as primer3's limit, which conflated two unrelated facts: a
    synthesiser will happily make a 90-mer, while primer3's thermodynamic
    alignment refuses above 60 nt.  A primer past 60 is therefore ordered as
    designed and its structure numbers describe only its 3'-most 60 bases,
    which the tables say in ``primer3_clamped``.  Setting this above 60 is how
    a 40 bp homology at a 65 C target becomes designable at all.
    """

    footprint_tm_window: TmWindow = TmWindow(min=52.0, max=72.0)
    pair_tm_delta_max: float = Field(default=4.0, ge=0.0, le=40.0)
    hairpin_dg_min: float = Field(default=-6.0, le=0.0)
    homodimer_dg_min: float = Field(default=-6.0, le=0.0)
    heterodimer_dg_min: float = Field(default=-6.0, le=0.0)
    three_prime_dimer_dg_min: float = Field(default=-9.0, le=0.0)
    in_tail_max_bp: int = Field(default=40, ge=0)
    """Longest part that may ride in the primer tails rather than be amplified.

    pydna's ``assembly_fragments`` calls this ``maxlink`` and defaults it to 40,
    which is where the number comes from.  Above it the part comes back as its
    own fragment with a warning."""

    hairpin_tm_max: float = Field(default=45.0, ge=0.0, le=120.0)
    homodimer_tm_max: float = Field(default=45.0, ge=0.0, le=120.0)
    heterodimer_tm_max: float = Field(default=45.0, ge=0.0, le=120.0)
    primer_len_max: int = Field(default=60, ge=20, le=200)
    """What a synthesiser will make, not what primer3 will score.

    The default stays 60 because that is where primer3 stops and so where the
    structure numbers stop being whole-molecule; raising it is a deliberate
    choice to order longer oligos and read their ΔG as describing the 3' end
    only."""

    gc_fraction_min: float = Field(default=0.25, ge=0.0, le=1.0)
    gc_fraction_max: float = Field(default=0.75, ge=0.0, le=1.0)

    @model_validator(mode="after")
    def _gc_ordered(self) -> "Constraints":
        if self.gc_fraction_min >= self.gc_fraction_max:
            raise ValueError(
                f"gc_fraction_min {self.gc_fraction_min} must be below "
                f"gc_fraction_max {self.gc_fraction_max}"
            )
        return self


class Weights(_Frozen):
    """Cost terms ranked over candidate splits when ``accept="argmin"``.

    Only quantities the split actually moves: the footprints are fixed before
    any split is tried, and ``t_f + t_r = L`` fixes the total oligo length, so
    the structure energies of the two tailed primers and the length difference
    between them are what is left to choose on.
    """

    hairpin_dg: float = Field(default=1.0, ge=0.0)
    homodimer_dg: float = Field(default=1.0, ge=0.0)
    heterodimer_dg: float = Field(default=1.0, ge=0.0)
    length_imbalance: float = Field(default=0.05, ge=0.0)


class Objective(_Frozen):
    weights: Weights = Weights()
    price_per_base: float = Field(default=1.8, ge=0.0)
    currency: Annotated[str, StringConstraints(pattern=r"^[A-Z]{3}$")] = "DKK"


class DesignSpec(_Frozen):
    design: Design
    positions: Tuple[Position, ...] = Field(min_length=2)
    pads: Tuple[Pad, ...] = ()
    primer_design: PrimerDesign = PrimerDesign()
    constraints: Constraints = Constraints()
    objective: Objective = Objective()

    @model_validator(mode="after")
    def _coherent(self) -> "DesignSpec":
        self._check_pads()
        self._check_ids()
        self._check_size()
        self._check_flags()
        self._check_feasibility()
        return self

    def _check_ids(self) -> None:
        """Unique ids, and none of them containing the separator.

        Composite ids are built by joining part ids with
        :data:`ID_SEPARATOR`, so a part id containing the separator makes that
        construction ambiguous: parts ``P`` and ``Q__R`` next to each other
        give the same composite id as ``P__Q`` and ``R``, and two different
        molecules, variants or junctions then answer to one name.  Rejecting it
        here is what makes every composite id in the package injective.
        """
        position_ids = [p.id for p in self.positions]
        if len(set(position_ids)) != len(position_ids):
            raise ValueError("position ids must be unique")
        fragment_ids = [part.id for p in self.positions for part in p.parts]
        fragment_ids += [pad.id for pad in self.pads]
        duplicates = sorted({i for i in fragment_ids if fragment_ids.count(i) > 1})
        if duplicates:
            raise ValueError(
                "part and pad ids must be unique across the whole spec; "
                f"repeated: {', '.join(duplicates)}"
            )
        # A pad sits in two namespaces, so it is checked in both.  It is a
        # fragment, hence the check above; but ``ordered_slots`` interleaves
        # pads with positions, so it is also a *slot*, and a pad id equal to a
        # position id gives two slots one name.  The dict keyed on that name
        # keeps whichever came last, and an ``after_position`` pad comes last,
        # so every variant reported the pad as the part chosen at that
        # position.  A part id may still match a position id: parts are never
        # slot keys.
        slots = sorted(set(position_ids) & {pad.id for pad in self.pads})
        if slots:
            raise ValueError(
                "a pad id may not equal a position id; both are slots and the "
                f"two would collapse into one: {', '.join(slots)}"
            )
        composite = sorted(
            {i for i in fragment_ids + position_ids if ID_SEPARATOR in i}
        )
        if composite:
            raise ValueError(
                f"an id may not contain {ID_SEPARATOR!r}, which joins composite "
                f"ids and would make them ambiguous; rename: "
                f"{', '.join(composite)}"
            )

    def _check_flags(self) -> None:
        """What the slot flags need in order to mean anything.

        ``in_tail`` removes a slot from the fragment list and splices it into
        the junction between its neighbours, so it needs neighbours on both
        sides and a single, known sequence: with variants, which junction you
        get would depend on the choice made at a slot that is not a fragment.
        ``no_tail_fwd`` and ``no_tail_rev`` each pin a junction's split to one
        end, so facing each other across one junction they ask for both ends at
        once.
        """
        for position in self.positions:
            if not position.in_tail:
                continue
            if len(position.parts) != 1:
                raise ValueError(
                    f"position {position.id!r} is in_tail with "
                    f"{len(position.parts)} variants; an in_tail slot rides in "
                    "its neighbours' primer tails, so it must have exactly one "
                    "part"
                )
        slots = [
            slot_id
            for slot_id, _ in self.ordered_slots()
            if slot_id not in {p.id for p in self.positions if p.in_tail}
        ]
        if len(slots) < 2:
            raise ValueError(
                "fewer than two fragments are left once the in_tail slots are "
                "built into the tails; there is nothing to join"
            )
        order = [slot_id for slot_id, _ in self.ordered_slots()]
        by_id = {position.id: position for position in self.positions}
        in_tail_ids = {p.id for p in self.positions if p.in_tail}
        if self.design.topology == "linear":
            for edge in (order[0], order[-1]):
                if edge in in_tail_ids:
                    raise ValueError(
                        f"position {edge!r} is in_tail at the end of a linear "
                        "design, where there is no neighbour on both sides to "
                        "carry it"
                    )
        # Two no_tail flags facing each other across one junction.
        fragments = [slot_id for slot_id in order if slot_id not in in_tail_ids]
        edges = list(zip(fragments, fragments[1:]))
        if self.design.topology == "circular":
            edges.append((fragments[-1], fragments[0]))
        for upstream, downstream in edges:
            if (
                getattr(by_id.get(upstream), "no_tail_rev", False)
                and getattr(by_id.get(downstream), "no_tail_fwd", False)
            ):
                raise ValueError(
                    f"the junction {upstream!r} -> {downstream!r} has "
                    "no_tail_rev upstream and no_tail_fwd downstream, so "
                    "neither primer may carry the homology"
                )

    def _check_pads(self) -> None:
        known = {p.id for p in self.positions}
        for pad in self.pads:
            anchor = pad.before_position or pad.after_position
            if anchor not in known:
                raise ValueError(
                    f"pad {pad.id!r} anchors to unknown position {anchor!r}"
                )

    def _check_size(self) -> None:
        n_parts = sum(len(p.parts) for p in self.positions) + len(self.pads)
        if n_parts > MAX_PARTS_TOTAL:
            raise ValueError(f"{n_parts} parts exceeds the cap of {MAX_PARTS_TOTAL}")
        if self.n_variants > MAX_VARIANTS:
            raise ValueError(
                f"{self.n_variants} variants exceeds the cap of {MAX_VARIANTS}"
            )

    def _check_feasibility(self) -> None:
        design = self.primer_design
        window = self.constraints.footprint_tm_window
        if not window.min <= design.target_tm <= window.max:
            raise ValueError(
                f"primer_design.target_tm {design.target_tm} lies outside "
                f"footprint_tm_window [{window.min}, {window.max}], so every "
                "footprint would be rejected"
            )
        # One tail of every junction is at least ceil(L/2), so this bound is the
        # shortest possible longest oligo: above the cap, no split can work.
        shortest_possible = design.footprint_len_min + math.ceil(
            design.junction.homology_length / 2
        )
        if shortest_possible > self.constraints.primer_len_max:
            raise ValueError(
                f"homology_length {design.junction.homology_length} with "
                f"footprint_len_min {design.footprint_len_min} needs at least "
                f"{shortest_possible} nt, over primer_len_max "
                f"{self.constraints.primer_len_max}"
            )
        weights = self.objective.weights.model_dump().values()
        if design.junction.split_search.accept == "argmin" and not any(weights):
            raise ValueError(
                "accept='argmin' with all-zero objective weights would pick an "
                "arbitrary split"
            )

    def ordered_slots(self) -> Tuple[Tuple[str, Tuple[Part, ...]], ...]:
        """Fragments in assembly order: ``(slot_id, candidate parts)`` per slot.

        Pads are resolved to their place here, once, and contribute a single
        candidate; a position contributes its parts.  The product over the slots
        is the library.
        """
        before: dict = {}
        after: dict = {}
        for pad in self.pads:
            anchor = pad.before_position or pad.after_position
            bucket = before if pad.before_position else after
            bucket.setdefault(anchor, []).append(pad)
        slots = []
        for position in self.positions:
            for pad in before.get(position.id, ()):
                slots.append((pad.id, (pad,)))
            slots.append((position.id, position.parts))
            for pad in after.get(position.id, ()):
                slots.append((pad.id, (pad,)))
        return tuple(slots)

    @property
    def n_variants(self) -> int:
        return math.prod(len(p.parts) for p in self.positions)


def load_spec(raw: Mapping) -> DesignSpec:
    """Validate a parsed JSON document.

    Raises :class:`pydantic.ValidationError`, which names the offending field
    path; the caller turns that into a diagnostic.
    """
    return DesignSpec.model_validate(raw)


def json_schema() -> dict:
    return DesignSpec.model_json_schema()
