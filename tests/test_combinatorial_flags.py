#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)

"""The six slot flags, which come from FLOW's ``construct_part``.

They are *slot-level*: they describe the place on the DNA, so every alternative
at one position shares them.  Three change the sequence (``reverse``,
``drop_start``, ``drop_stop``), two pin a junction's split (``no_tail_fwd``,
``no_tail_rev``) and one takes a part out of the fragment list altogether
(``in_tail``).
"""

from __future__ import annotations

import random

import pytest
from Bio.Seq import Seq
from pydantic import ValidationError

import teemi.design.combinatorial as cb
from teemi.design.combinatorial.design import apply_slot_flags, plan_layout
from teemi.design.combinatorial.sources import SourceResolver
from teemi.design.combinatorial.spec import load_spec

from tests.combinatorial_netguard import no_network  # noqa: F401

CDS = "ATG" + "GCTACGGATCATTACGGCAGTTCAGGACTTGACC" * 3 + "TAA"


def dna(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def spec(positions, topology="linear", **primer_design):
    junction = {
        "homology_length": 30,
        # These tests are about the flags, so they pin the path they run on
        # rather than inheriting whichever is the default.
        "tune_overhangs": True,
        "split_search": {"order": "exhaustive", "accept": "argmin"},
    }
    junction.update(primer_design.pop("junction", {}))
    return {
        "design": {"name": "flags", "topology": topology},
        "positions": positions,
        "primer_design": {"target_tm": 57.0, "junction": junction, **primer_design},
    }


def slot(part_id, sequence, **flags):
    return {
        "id": f"at_{part_id}",
        "parts": [{"id": part_id, "source": {"type": "inline", "sequence": sequence}}],
        **flags,
    }


def resolved(positions, **kwargs):
    validated = load_spec(spec(positions, **kwargs))
    records = SourceResolver().resolve_spec(validated)
    return apply_slot_flags(records, validated)


# --------------------------------------------------------------------------- #
# The three that change the sequence
# --------------------------------------------------------------------------- #


def test_drop_start_and_drop_stop_remove_exactly_their_codon():
    records, notes = resolved(
        [
            slot("cds", CDS, drop_start=True, drop_stop=True),
            slot("other", dna(400, 1)),
        ]
    )
    assert str(records["cds"].seq) == CDS[3:-3]
    assert notes == ()


def test_reverse_takes_the_reverse_complement():
    records, notes = resolved(
        [slot("cds", CDS, reverse=True), slot("other", dna(400, 2))]
    )
    assert str(records["cds"].seq) == str(Seq(CDS).reverse_complement())
    assert records["cds"].id == "cds", "the part keeps its id through the flip"
    assert notes == ()


def test_the_codons_go_before_the_flip_not_after():
    """A start codon is a fact about the part's own frame, not the construct.

    Flipping first would look for the ATG at the 5' end of the oriented
    sequence, where a reversed CDS presents the reverse complement of its STOP
    codon instead -- so the flag would miss, and would trim the wrong three
    bases on any part whose rc happens to begin ATG.
    """
    records, notes = resolved(
        [slot("cds", CDS, reverse=True, drop_start=True), slot("other", dna(400, 3))]
    )
    assert str(records["cds"].seq) == str(Seq(CDS[3:]).reverse_complement())
    assert notes == (), "the ATG is found in the part's own reading"


def test_a_codon_that_is_not_there_is_reported_and_nothing_is_cut():
    """Trimming three bases that were not a codon shifts the fusion's frame."""
    plain = dna(400, 4)
    records, notes = resolved(
        [slot("p", plain, drop_start=True), slot("other", dna(400, 5))]
    )
    assert str(records["p"].seq) == plain, "untouched"
    assert [n.code for n in notes] == ["drop_start_codon_absent"]
    assert notes[0].severity == "warning"
    assert notes[0].part_id == "p"


# --------------------------------------------------------------------------- #
# The two that pin a junction's split
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "flags, expected_t_f, bare",
    [
        ({}, None, None),
        ({"upstream": {"no_tail_rev": True}}, 30, "reverse"),
        ({"downstream": {"no_tail_fwd": True}}, 0, "forward"),
    ],
)
def test_a_no_tail_flag_leaves_that_primer_a_bare_footprint(flags, expected_t_f, bare):
    result = cb.plan(
        spec(
            [
                slot("a", dna(400, 6), **flags.get("upstream", {})),
                slot("b", dna(400, 7), **flags.get("downstream", {})),
            ]
        ),
        sources=SourceResolver(),
        offline=True,
    )
    junction = result.junctions[0]
    if expected_t_f is None:
        assert len(junction.profile) > 1, "an unpinned junction searches"
        return
    assert junction.chosen.t_f == expected_t_f
    assert len(junction.profile) == 1, "pinned, so there is one candidate"
    by_id = {p.id: p for p in result.primers}
    tailed = by_id["a__b_fw" if bare == "forward" else "a__b_rv"]
    assert tailed.sequence == tailed.footprint, "no tail at all"


def test_two_no_tail_flags_facing_each_other_are_refused():
    """Neither primer may carry the homology, so there is no junction to make."""
    with pytest.raises(ValidationError, match="neither primer may carry"):
        load_spec(
            spec(
                [
                    slot("a", dna(400, 8), no_tail_rev=True),
                    slot("b", dna(400, 9), no_tail_fwd=True),
                ]
            )
        )


def test_a_pinned_split_the_parts_cannot_supply_says_so():
    short = dna(60, 10)
    result = cb.plan(
        spec(
            [
                slot("a", dna(400, 11)),
                slot("b", short, no_tail_fwd=True),
                slot("c", dna(400, 12)),
            ],
            junction={"homology_length": 80},
        ),
        sources=SourceResolver(),
        offline=True,
    )
    assert not result.feasible
    assert any("no_tail flag pins the split" in d.message for d in result.diagnostics)


# --------------------------------------------------------------------------- #
# in_tail: the one that changes what gets amplified
# --------------------------------------------------------------------------- #


KOZAK = "AAAACAATG"


def in_tail_spec(middle=KOZAK):
    return spec(
        [
            slot("a", dna(400, 13)),
            slot("k", middle, in_tail=True),
            slot("b", dna(400, 14)),
        ]
    )


def test_an_in_tail_part_is_not_amplified_but_is_in_the_product():
    result = cb.plan(
        in_tail_spec(), sources=SourceResolver(), offline=True, include_designs=True
    )
    assert result.feasible, [d.message for d in result.errors]
    variant = result.variants[0]

    assert variant.fragment_ids == ("a", "b"), "one fewer PCR than parts"
    assert len(result.pcrs) == 2
    assert [p.id for p in result.spec.positions] == ["at_a", "at_k", "at_b"]

    product = str(result.designs[0].construct.seq)
    assert len(product) == variant.length_bp == 400 + len(KOZAK) + 400
    assert product[400 : 400 + len(KOZAK)] == KOZAK, "and in the right place"


def test_what_rides_in_the_tails_is_inside_the_homology_not_added_to_it():
    """Asking for 30 bp of homology gives 30, Kozak included.

    This is pydna's accounting in ``assembly_fragments``, where a fragment
    below ``maxlink`` is distributed into the flanking primers and the
    ``overlap`` asked for is still what the amplicons share.  Adding the
    linker's length on top instead would silently hand back a longer homology
    and two longer primers than the spec asked for.
    """
    with_it = cb.plan(in_tail_spec(), sources=SourceResolver(), offline=True)
    junction = with_it.junctions[0]
    assert junction.interposed == KOZAK
    assert len(junction.homology) == 30, "exactly what was asked for"
    assert KOZAK in junction.homology, "and it is in there"
    assert junction.chosen.t_f + junction.chosen.t_r == 30 - len(KOZAK)

    carriers = [
        p
        for p in with_it.primers
        if KOZAK in p.sequence
        or KOZAK in str(Seq(p.sequence).reverse_complement())
    ]
    assert len(carriers) == 2, [p.id for p in carriers]


def test_the_product_matches_what_pydna_assembly_fragments_builds():
    """pydna is the oracle for the geometry, not just for the arithmetic.

    ``assembly_fragments`` is the maintained implementation of exactly this
    manoeuvre -- amplicons with a short ``Dseqrecord`` between them, which it
    distributes into the two flanking primers.  It offers one split
    (``ceil(overlap/2)`` each side) where this module sweeps, so it is not
    called directly; instead the construct it builds is what this module's
    construct has to equal.
    """
    from pydna.design import assembly_fragments, primer_design
    from pydna.dseqrecord import Dseqrecord
    from pydna.assembly import Assembly

    a, b = dna(400, 13), dna(400, 14)
    result = cb.plan(
        in_tail_spec(), sources=SourceResolver(), offline=True, include_designs=True
    )
    assert result.feasible, [d.message for d in result.errors]
    mine = str(result.designs[0].construct.seq)

    fragments = assembly_fragments(
        [
            primer_design(Dseqrecord(a), target_tm=57.0, limit=18),
            Dseqrecord(KOZAK),
            primer_design(Dseqrecord(b), target_tm=57.0, limit=18),
        ],
        overlap=30,
        maxlink=40,
    )
    assert len(fragments) == 2, "pydna folds the short one into its neighbours too"
    products = [
        p
        for p in Assembly(list(fragments), limit=25).assemble_linear()
        if len(p) >= max(len(f) for f in fragments)
    ]
    assert len(products) == 1
    assert str(products[0].seq) == mine == a + KOZAK + b


def test_a_part_too_long_to_ride_in_the_tails_is_amplified_instead():
    """FLOW's schema says as much about teemi's maxlink: a request, not a command."""
    result = cb.plan(
        in_tail_spec(dna(60, 15)), sources=SourceResolver(), offline=True
    )
    assert any(d.code == "in_tail_too_long" for d in result.diagnostics)
    assert result.variants[0].fragment_ids == ("a", "k", "b"), "its own PCR after all"
    assert len(result.pcrs) == 3


def test_an_in_tail_part_needs_no_footprint_of_its_own():
    """A 9 bp Kozak is no template, so nothing anneals to it.

    Asking for a footprint anyway refused the part for being too short to do
    something it was never going to do.
    """
    validated = load_spec(in_tail_spec())
    records = SourceResolver().resolve_spec(validated)
    layout = plan_layout(validated, records)
    assert [slot_id for slot_id, _, tail in layout.slots if tail] == ["at_k"]
    assert len(records["k"]) < validated.primer_design.footprint_len_min

    result = cb.plan(in_tail_spec(), sources=SourceResolver(), offline=True)
    assert "k" not in {f.part_id for f in result.footprints}
    assert result.feasible


def test_in_tail_on_a_combinatorial_position_is_refused():
    """Which junction you get would depend on a choice made at a non-fragment."""
    with pytest.raises(ValidationError, match="must have exactly one part"):
        load_spec(
            spec(
                [
                    slot("a", dna(400, 16)),
                    {
                        "id": "at_k",
                        "in_tail": True,
                        "parts": [
                            {
                                "id": "k1",
                                "source": {"type": "inline", "sequence": KOZAK},
                            },
                            {
                                "id": "k2",
                                "source": {
                                    "type": "inline",
                                    "sequence": "AAAACAGTG",
                                },
                            },
                        ],
                    },
                    slot("b", dna(400, 17)),
                ]
            )
        )


def test_in_tail_at_the_end_of_a_linear_design_is_refused():
    with pytest.raises(ValidationError, match="no neighbour on both sides"):
        load_spec(
            spec(
                [
                    slot("a", dna(400, 18)),
                    slot("b", dna(400, 19)),
                    slot("k", KOZAK, in_tail=True),
                ]
            )
        )


# --------------------------------------------------------------------------- #
# Which implementation builds the tails
# --------------------------------------------------------------------------- #


def plain(tune):
    return spec(
        [slot("a", dna(400, 30)), slot("b", dna(400, 31))],
        junction={"tune_overhangs": tune},
    )


def test_pydna_builds_the_tails_unless_the_spec_asks_for_the_search():
    """The default is the maintained implementation, not this module's.

    ``assembly_fragments`` offers one balanced split, so the default path has a
    single candidate and nothing to choose between -- which is also why
    ``unavoidable`` stays empty on it: one candidate proves nothing.
    """
    by_pydna = cb.plan(plain(False), sources=SourceResolver(), offline=True)
    junction = by_pydna.junctions[0]
    assert junction.designed_by == "pydna"
    assert len(junction.profile) == 1
    assert junction.unavoidable == ()

    searched = cb.plan(plain(True), sources=SourceResolver(), offline=True)
    assert searched.junctions[0].designed_by == "sweep"
    assert len(searched.junctions[0].profile) > 1


def test_the_default_path_is_pydnas_own_output_base_for_base():
    """Not an imitation of ``assembly_fragments`` -- its actual result."""
    from pydna.amplify import pcr
    from pydna.design import assembly_fragments, primer_design
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer

    a, b = dna(400, 30), dna(400, 31)
    result = cb.plan(plain(False), sources=SourceResolver(), offline=True)
    junction = result.junctions[0]

    footprints = {f.part_id: f for f in result.footprints}
    amplicons = [
        pcr(
            Primer(footprints[part].forward),
            Primer(footprints[part].reverse),
            Dseqrecord(sequence),
        )
        for part, sequence in (("a", a), ("b", b))
    ]
    tailed = assembly_fragments(
        amplicons,
        overlap=result.spec.primer_design.junction.homology_length,
        maxlink=result.spec.constraints.in_tail_max_bp,
    )
    assert junction.chosen.forward_primer == str(tailed[-1].forward_primer.seq)
    assert junction.chosen.reverse_primer == str(tailed[0].reverse_primer.seq)


def test_a_no_tail_flag_is_built_here_because_pydna_cannot_say_it():
    """``assembly_fragments`` has no way to leave one side bare.

    So that junction is built by this module even with the search off, and says
    so rather than letting a design believe pydna produced every primer.
    """
    result = cb.plan(
        spec(
            [
                slot("a", dna(400, 32), no_tail_rev=True),
                slot("b", dna(400, 33)),
            ],
            junction={"tune_overhangs": False},
        ),
        sources=SourceResolver(),
        offline=True,
    )
    junction = result.junctions[0]
    assert junction.designed_by == "forced"
    assert junction.chosen.t_r == 0
    by_id = {p.id: p for p in result.primers}
    assert by_id["a__b_rv"].sequence == by_id["a__b_rv"].footprint
