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

"""FLOW's design document as the input, since it is what drives the app.

:class:`FlowDocument` mirrors ``heterologous-expression.schema.json`` -- every
field, every enum, every rule it states -- so a document the frontend sends is
checked here and nothing downstream has to guess what it meant.

COMPOSITION, NOT REPLACEMENT.  FLOW's schema says of itself:

    It carries no primer-design policy, no reaction conditions and no objective
    weights.  Those belong to the backend's own spec, which owns them and
    defaults every one.  FLOW sends what the user decided and nothing else, so
    an absent field means 'not decided' rather than 'decided to be the default'.

So this module does not replace :mod:`teemi.design.combinatorial.spec`; it
composes with it.  The document supplies the construct, the variants, the slot
flags and the assembly settings the user chose; ``policy`` supplies the bounds,
the weights and the price, which FLOW deliberately never sends.  Both end up in
one :class:`~teemi.design.combinatorial.spec.DesignSpec`, which stays the
internal model.

Mirrored in pydantic rather than validated against the JSON Schema on purpose:
pydantic is already a dependency, the errors name the field, and the result is
one validation boundary instead of two.  The schema remains the authority, and
``tests/test_combinatorial_flow.py`` cross-checks the mirror against it.

A ``library`` source is a ``part_id`` and the backend owns resolving it -- the
schema is explicit that the sequence comes from the id and never from the name
-- so :func:`from_flow` takes the mapping rather than going looking for one.
"""

from __future__ import annotations

from typing import Dict, Literal, Mapping, Optional, Tuple, Union

from pydantic import BaseModel, ConfigDict, Field, StringConstraints, model_validator
from typing_extensions import Annotated

from teemi.design.combinatorial.spec import DesignSpec, load_spec

__all__ = [
    "FlowDocument",
    "TRANSLATED_TYPES",
    "UndecidedPosition",
    "UnknownPart",
    "from_flow",
    "load_flow",
    "plan_flow",
]

#: The part types a start or stop codon means anything for, as FLOW's own
#: ``drop_start`` description lists them.  A ``drop_`` flag on anything else is
#: reported: it is not wrong to send, but nothing will come of it.
TRANSLATED_TYPES = frozenset(
    {
        "selection_marker",
        "fluorescent_protein",
        "signal_peptide",
        "cds",
        "linker_peptide",
        "nls_tag",
    }
)

#: ``Method.circular`` resolved, which the schema's ``allOf`` pins per method.
TOPOLOGY_BY_METHOD = {
    "he-invivo-plasmid": "circular",
    "he-genomic-integration": "linear",
}

PartId = Annotated[str, StringConstraints(pattern=r"^[A-Za-z0-9_]+$")]
NonEmpty = Annotated[str, StringConstraints(min_length=1)]


class _Frozen(BaseModel):
    """``extra="forbid"`` mirrors the schema's ``additionalProperties: false``,
    so a key FLOW does not define is an error rather than a silent no-op."""

    model_config = ConfigDict(frozen=True, extra="forbid", populate_by_name=True)


class Genome(_Frozen):
    """A REFERENCE, never the sequence, as the schema says."""

    accession: Annotated[str, StringConstraints(pattern=r"^GC[AF]_[0-9]{9}\.[0-9]+$")]
    assembly_name: Optional[str] = None
    strain: Optional[str] = None


class Organism(_Frozen):
    key: Literal[
        "a-nidulans",
        "a-oryzae",
        "a-niger",
        "t-reesei",
        "s-cerevisiae",
        "p-pastoris",
        "d-hansenii",
    ]
    genome: Optional[Genome] = None


class LibrarySource(_Frozen):
    type: Literal["library"] = "library"
    part_id: PartId


class InlineSource(_Frozen):
    type: Literal["inline"] = "inline"
    sequence: NonEmpty
    format: Literal["fasta", "genbank"]
    upload_id: Optional[str] = None


Source = Union[LibrarySource, InlineSource]


class Variant(_Frozen):
    """One alternative at a combinatorial position: a CHOICE, not a segment."""

    id: NonEmpty
    name: Optional[str] = None
    source: Optional[Source] = None


class ConstructPart(_Frozen):
    """One position on the DNA.  Slot-level facts live here, content on the
    source or on each variant, so alternatives cannot disagree about the slot."""

    id: NonEmpty
    type: Literal[
        "selection_marker",
        "replicator",
        "fluorescent_protein",
        "platform_flank",
        "promoter",
        "kozak_like",
        "signal_peptide",
        "cds",
        "linker_peptide",
        "nls_tag",
        "terminator",
        "protospacer",
        "sgRNA_scaffold",
        "crRNA_direct_repeat",
        "trna",
        "ribozyme",
    ]
    name: Optional[str] = None
    notes: Optional[str] = None
    source: Optional[Source] = None
    variants: Optional[Tuple[Variant, ...]] = Field(default=None, min_length=1)
    reverse: bool = False
    in_tail: bool = False
    no_tail_fwd: bool = False
    no_tail_rev: bool = False
    drop_start: bool = False
    drop_stop: bool = False

    @model_validator(mode="after")
    def _one_source_of_truth(self) -> "ConstructPart":
        """The schema's ``not: required [source, variants]``.

        Its reason, verbatim: two sources of truth for "what goes here" is a bug
        waiting to be read the wrong way.
        """
        if self.source is not None and self.variants is not None:
            raise ValueError(
                f"position {self.id!r} has both source and variants; when "
                "variants is set it is the sole source of truth for what can "
                "sit here"
            )
        return self


class OligoPart(_Frozen):
    """Something ordered alongside the construct, not sitting in it.

    The slot flags are absent by design: none of them means anything off the
    DNA.  Nothing here is assembled, so these do not reach the plan.
    """

    id: NonEmpty
    type: Literal["crRNA_oligo", "repair_oligo", "ddpcr_probe"]
    name: Optional[str] = None
    notes: Optional[str] = None
    source: Optional[Source] = None


class AssemblySettings(_Frozen):
    """In vivo assembly settings as collected, with no ranges enforced.

    Following FLOW: it does not know what a sensible Tm or overlap is for your
    host and polymerase, and a made-up limit that quietly rejects a valid
    number is worse than no limit.  Every quantity carries its unit in the
    name, because the one real bug in this area is 0.5 uM and 250 nM being the
    same number written two ways.
    """

    homology_bp: Optional[int] = Field(default=None, ge=0)
    target_tm_c: Optional[float] = None
    primer_conc_um: Optional[float] = Field(default=None, gt=0)
    min_footprint_nt: Optional[int] = Field(default=None, ge=1)
    polymerase: Optional[str] = None
    tune_overhangs: Optional[bool] = None


class FlowDocument(_Frozen):
    """A whole FLOW heterologous-expression design, as the frontend sends it."""

    flow_schema_version: Literal["1.0"]
    task: Literal["heterologous-expression"]
    method: Literal["he-invivo-plasmid", "he-genomic-integration"]
    organism: Organism
    #: The wire name is ``construct``; the attribute is ``positions`` because
    #: ``construct`` shadows a pydantic ``BaseModel`` member.  Both names
    #: validate, so a document is read exactly as FLOW writes it.
    positions: Tuple[ConstructPart, ...] = Field(alias="construct", min_length=1)
    name: Optional[str] = Field(default=None, max_length=200)
    topology: Optional[Literal["circular", "linear"]] = None
    oligos: Tuple[OligoPart, ...] = ()
    assembly: AssemblySettings = AssemblySettings()
    design_count: Optional[int] = Field(default=None, ge=1)
    library_version: Optional[str] = None

    @model_validator(mode="after")
    def _topology_matches_method(self) -> "FlowDocument":
        """The schema's two ``allOf`` branches: the method pins the topology."""
        expected = TOPOLOGY_BY_METHOD[self.method]
        if self.topology is not None and self.topology != expected:
            raise ValueError(
                f"method {self.method!r} is {expected} but the document says "
                f"topology {self.topology!r}"
            )
        return self

    @property
    def resolved_topology(self) -> str:
        """The topology to build with, from the method when it is not sent.

        FLOW resolves it from ``Method.circular`` rather than leaving the
        backend to infer, but it is not a required field, so it is derived here
        when absent instead of guessed at later.
        """
        return self.topology or TOPOLOGY_BY_METHOD[self.method]


class UndecidedPosition(ValueError):
    """A position with no part chosen, which FLOW explicitly allows.

    A block drawn amber reading "Choose a part" is a legitimate state of a saved
    design, so this is not a malformed document -- there is simply nothing to
    design against that position yet.  ``fill`` stands a sequence in when a
    preview is wanted before every choice has been made.
    """

    def __init__(self, position_id: str) -> None:
        super().__init__(
            f"position {position_id!r} has no part chosen yet, so there is "
            "nothing to design against it; pass fill={{{0!r}: sequence}} to "
            "stand one in".format(position_id)
        )
        self.position_id = position_id


class UnknownPart(ValueError):
    """A ``part_id`` the library handed in does not have.

    FLOW's schema notes that ``parts.json`` is generated, that ids have been
    retired before, and that its own lookup returns undefined silently -- so an
    unpinned reference can fail quietly.  Here it fails loudly.
    """

    def __init__(self, part_id: str, position_id: str) -> None:
        super().__init__(
            f"position {position_id!r} wants part_id {part_id!r}, which the "
            "library does not have; check library_version"
        )
        self.part_id = part_id


def _sequence(
    source: Source, position_id: str, library: Mapping[str, str]
) -> Tuple[str, str]:
    """``(part id, DNA)`` for one source.

    A library source resolves through the mapping, because the schema is
    explicit that the backend resolves from the id and never from the name.  An
    inline source carries the uploaded file's *text*, so it is parsed -- with
    Biopython over a string buffer rather than pydna's parser, which reads a
    string that looks like a path and would turn browser text into a file read.
    """
    if isinstance(source, LibrarySource):
        if source.part_id not in library:
            raise UnknownPart(source.part_id, position_id)
        return source.part_id, library[source.part_id]

    from io import StringIO

    from Bio import SeqIO

    record = SeqIO.read(StringIO(source.sequence), source.format)
    return "", str(record.seq).upper()


def from_flow(
    document: Union[FlowDocument, Mapping],
    *,
    library: Mapping[str, str],
    policy: Optional[Mapping] = None,
    fill: Optional[Mapping[str, str]] = None,
) -> Tuple[DesignSpec, Tuple[str, ...]]:
    """A FLOW document plus the backend's policy -> one :class:`DesignSpec`.

    ``library`` maps ``part_id`` to DNA.  ``policy`` carries what FLOW says it
    never sends -- ``constraints``, ``objective``, and the ``split_search``
    detail behind ``tune_overhangs`` -- and every one of those has a default, so
    ``policy`` is optional.

    ``fill`` stands DNA in where a document has made no choice yet.  A key is
    either a position id, which replaces that whole position with one part, or
    ``"<position>/<variant>"``, which supplies one alternative at a
    combinatorial position -- the form :class:`UndecidedPosition` names, and the
    only one that can give four different genes to one slot.

    Returns the spec and any notes about the document: a ``drop_`` flag on a
    part type no codon applies to, a ``fill`` that was used, an ``oligos`` list
    that is not part of the construct.  They are notes rather than errors
    because none of them makes the design unbuildable.
    """
    doc = (
        document
        if isinstance(document, FlowDocument)
        else FlowDocument.model_validate(document)
    )
    fill = dict(fill or {})
    policy = dict(policy or {})
    notes: list = []
    unused = set(fill)

    positions = []
    for slot in doc.positions:
        flags = {
            flag: True
            for flag in (
                "reverse",
                "in_tail",
                "no_tail_fwd",
                "no_tail_rev",
                "drop_start",
                "drop_stop",
            )
            if getattr(slot, flag)
        }
        if (slot.drop_start or slot.drop_stop) and slot.type not in TRANSLATED_TYPES:
            notes.append(
                f"{slot.id}: drop_start/drop_stop on a {slot.type}, which is not "
                "translated, so there is no codon of its own to drop"
            )

        if slot.id in fill:
            unused.discard(slot.id)
            notes.append(
                f"{slot.id}: a sequence was stood in for it by the caller, not "
                "chosen in the document"
            )
            parts = [
                {
                    "id": f"{slot.id}.filled",
                    "name": slot.name,
                    "source": {"type": "inline", "sequence": fill[slot.id]},
                }
            ]
        elif slot.variants is not None:
            parts = []
            for variant in slot.variants:
                key = f"{slot.id}/{variant.id}"
                if key in fill:
                    unused.discard(key)
                    notes.append(
                        f"{key}: a sequence was stood in for it by the caller, "
                        "not chosen in the document"
                    )
                    parts.append(
                        {
                            # Not the bare variant id: that is only unique
                            # within its position, and part ids are unique
                            # across the whole spec.
                            "id": f"{slot.id}.{variant.id}",
                            "name": variant.name,
                            "source": {"type": "inline", "sequence": fill[key]},
                        }
                    )
                    continue
                if variant.source is None:
                    raise UndecidedPosition(key)
                part_id, dna = _sequence(variant.source, slot.id, library)
                parts.append(
                    {
                        "id": part_id or f"{slot.id}_{variant.id}",
                        "name": variant.name,
                        "source": {"type": "inline", "sequence": dna},
                    }
                )
        else:
            if slot.source is None:
                raise UndecidedPosition(slot.id)
            part_id, dna = _sequence(slot.source, slot.id, library)
            parts = [
                {
                    "id": part_id or f"{slot.id}_part",
                    "name": slot.name,
                    "source": {"type": "inline", "sequence": dna},
                }
            ]
        positions.append({"id": slot.id, "name": slot.name, "parts": parts, **flags})

    if unused:
        raise KeyError(
            f"fill names {', '.join(sorted(unused))}, which is not a position "
            "or a position/variant in this document; a typo here would "
            "silently leave the real slot undecided"
        )

    if doc.oligos:
        notes.append(
            f"{len(doc.oligos)} oligo(s) sit alongside the construct and are not "
            f"assembled, so they are not in this plan: "
            f"{', '.join(o.id for o in doc.oligos)}"
        )

    settings = doc.assembly
    junction: Dict[str, object] = {}
    if settings.homology_bp is not None:
        junction["homology_length"] = settings.homology_bp
    if settings.tune_overhangs is not None:
        junction["tune_overhangs"] = settings.tune_overhangs
    # The one piece of policy that hangs off a FLOW field: how to search, once
    # tune_overhangs has said to.
    if "split_search" in policy:
        junction["split_search"] = policy["split_search"]

    primer_design: Dict[str, object] = {}
    if settings.target_tm_c is not None:
        primer_design["target_tm"] = settings.target_tm_c
    if settings.min_footprint_nt is not None:
        primer_design["footprint_len_min"] = settings.min_footprint_nt
    if settings.polymerase:
        # "It selects the Tm model; choosing that model is the backend's job."
        primer_design["polymerase"] = settings.polymerase
    if junction:
        primer_design["junction"] = junction

    raw: Dict[str, object] = {
        "design": {
            "name": doc.name or doc.method,
            "topology": doc.resolved_topology,
            "description": f"{doc.task} / {doc.method} / {doc.organism.key}",
        },
        "positions": positions,
    }
    if primer_design:
        raw["primer_design"] = primer_design
    for key in ("constraints", "objective"):
        if key in policy:
            raw[key] = policy[key]

    return load_spec(raw), tuple(notes)


def load_flow(
    path,
    *,
    library: Mapping[str, str],
    policy: Optional[Mapping] = None,
    fill: Optional[Mapping[str, str]] = None,
) -> Tuple[DesignSpec, Tuple[str, ...]]:
    """:func:`from_flow` on a document read from a JSON file."""
    import json
    from pathlib import Path

    return from_flow(
        json.loads(Path(path).read_text()),
        library=library,
        policy=policy,
        fill=fill,
    )


def plan_flow(
    document: Union[FlowDocument, Mapping],
    *,
    library: Mapping[str, str],
    policy: Optional[Mapping] = None,
    fill: Optional[Mapping[str, str]] = None,
    **plan_kwargs,
):
    """Design a FLOW document: :func:`from_flow` then
    :func:`~teemi.design.combinatorial.report.plan`.

    Returns ``(result, notes)``.  ``plan_kwargs`` go straight to ``plan``, so
    ``offline`` and ``include_designs`` work as they do there.  ``sources`` is
    not among them on purpose: every part is resolved to an inline sequence
    before the spec exists, so there is no path for a file source to come from.
    """
    from teemi.design.combinatorial.report import plan
    from teemi.design.combinatorial.sources import SourceResolver

    spec, notes = from_flow(
        document, library=library, policy=policy, fill=fill
    )
    return plan(spec, sources=SourceResolver(), **plan_kwargs), notes
