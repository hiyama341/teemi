#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)

"""The slot flags, which come from FLOW's ``construct_part``.

They are *slot-level*: they describe the place on the DNA, so every alternative
at one position shares them.  Three change the sequence (``reverse``,
``drop_start``, ``drop_stop``), two pin a junction's split (``no_tail_fwd``,
``no_tail_rev``) and two decide what is BOUGHT for the slot rather than
amplified -- ``in_tail``, which folds it into the flanking primers, and
``ordered``, which buys it as a synthetic double-stranded fragment.
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


def test_a_part_too_long_to_ride_in_the_tails_is_ordered_instead():
    """FLOW's schema says as much about teemi's maxlink: a request, not a command.

    And the answer to a declined request is a PURCHASE, not an extra PCR.  A
    part asked into the tails is short by construction, so "amplify it after
    all" is the one answer it usually cannot be given -- which is what
    :func:`test_the_case_that_used_to_lose_the_whole_library` below is about.
    """
    result = cb.plan(
        in_tail_spec(dna(60, 15)), sources=SourceResolver(), offline=True
    )
    note = next(d for d in result.diagnostics if d.code == "in_tail_too_long")
    assert note.severity == "warning"
    assert "ordered as a synthetic fragment" in note.message

    variant = result.variants[0]
    assert variant.fragment_ids == ("a", "b"), "still two PCRs, not three"
    assert variant.ordered_ids == ("k",)
    assert len(result.pcrs) == 2
    assert [f.part_id for f in result.ordered] == ["k"]
    assert result.ordered[0].length_bp == 60
    assert result.ordered[0].because == "in_tail_too_long"


def test_a_part_longer_than_the_homology_is_ordered_instead():
    """The second limit on an ``in_tail`` part, and the one that used to raise.

    ``in_tail_max_bp`` is raised to 80 so it is not the limit under test: at 60
    bp the part clears the physical cap and fails only the 30 bp homology,
    because what rides in the tails is spent out of the homology rather than
    added to it.  This used to reach :class:`JunctionError` and come back as
    ``impossible_junction`` with nothing designed -- no variants, no PCRs, a
    library lost to a number the design could simply have declined.
    """
    raw = in_tail_spec(dna(60, 15))
    raw["constraints"] = {"in_tail_max_bp": 80}

    result = cb.plan(raw, sources=SourceResolver(), offline=True)
    codes = [d.code for d in result.diagnostics]
    assert "in_tail_too_long" not in codes, "80 bp of cap, so not that limit"

    note = next(d for d in result.diagnostics if d.code == "in_tail_over_homology")
    assert note.severity == "warning"
    assert note.part_id == "k"
    assert "60 bp" in note.message and "homology_length is only 30" in note.message
    assert "-" not in note.message, "a shortfall is reported as a positive number"

    assert not any(d.code == "impossible_junction" for d in result.errors)
    assert result.variants[0].ordered_ids == ("k",), "bought, not amplified"
    assert len(result.pcrs) == 2


def test_the_case_that_used_to_lose_the_whole_library():
    """A 54 bp signal peptide at 30 bp homology: ANNEAL's own failure.

    This is the defect the ``ordered`` mode exists for.  54 bp is past a 30 bp
    homology, so it cannot ride in the tails; demoting it to its own fragment
    produced primers pydna then refused to amplify --

        No PCR product! Template 54 bp linear limit=13: No forward primers
        anneal... No reverse primers anneal

    -- so the variant came back ``no_product`` and the whole design was
    infeasible.  There was no third answer, and both of the two were
    impossible.  Buying the part is the third answer, and it DESIGNS, VERIFIES
    AND ASSEMBLES: one fewer PCR, every oligo inside the budget, and a product
    equal to the intended construct.

    ``primer_len_max`` is 100 because that is what ANNEAL sends, and a full
    30 bp tail on a footprint is past teemi's default 60.
    """
    middle = dna(54, 21)
    raw = in_tail_spec(middle)
    raw["constraints"] = {"primer_len_max": 100, "in_tail_max_bp": 63}

    result = cb.plan(
        raw, sources=SourceResolver(), offline=True, include_designs=True
    )
    assert result.feasible, [d.message for d in result.errors]
    variant = result.variants[0]
    assert variant.verification.verdict == "verified"
    assert variant.fragment_ids == ("a", "b") and variant.ordered_ids == ("k",)
    assert len(result.pcrs) == 2

    # No primers, no PCR and no footprint for the part that is bought, so it
    # cannot be read out of any of the three tables as if it had them.
    assert "k" not in {p.part_id for p in result.pcrs}
    assert all("k" not in p.part_ids for p in result.primers)
    assert "k" not in {f.part_id for f in result.footprints}

    product = str(result.designs[0].construct.seq)
    assert product == dna(400, 13) + middle + dna(400, 14)


def test_a_fragment_is_ordered_at_any_length_folding_cannot_reach():
    """Not a rescue for one length: the mode has no upper bound of its own.

    Folding is capped by the homology -- it needs ``homology_length >= part``
    -- so nothing from 31 bp up can ride in a 30 bp tail.  A bought fragment
    does not care how long it is, because each neighbour takes a tail of
    exactly ``homology_length`` from it whatever its length, which is the point
    asserted here: the oligos do not grow with the fragment, so 276 bp costs
    the same two tails as 54 bp.
    """
    longest = {}
    for length in (54, 93, 276):
        middle = dna(length, 22)
        raw = in_tail_spec(middle)
        raw["constraints"] = {"primer_len_max": 100, "in_tail_max_bp": 63}

        result = cb.plan(
            raw, sources=SourceResolver(), offline=True, include_designs=True
        )
        assert result.feasible, [d.message for d in result.errors]
        assert result.variants[0].verification.verdict == "verified"
        assert len(result.pcrs) == 2
        assert [f.length_bp for f in result.ordered] == [length]
        assert str(result.designs[0].construct.seq) == (
            dna(400, 13) + middle + dna(400, 14)
        )
        longest[length] = max(len(p.sequence) for p in result.primers)

    assert len(set(longest.values())) == 1, longest


def test_ordering_costs_each_neighbour_a_full_length_tail():
    """The one price of the mode, and it is reported as a number not a mood.

    A folded part splits the homology between two tails; a bought one is a
    molecule of its own, so each neighbour carries the WHOLE homology.  At
    teemi's default ``primer_len_max`` of 60 a 30 bp tail on a 32 nt footprint
    is 62 nt and breaches it -- which is why ANNEAL sends 100.

    The 31 bp part is the case that used to have no answer at all: past a 30 bp
    homology, so it could not fold, and too short for pydna to find a 57 C
    footprint in, so it could not be amplified either.  ``no_footprint`` is now
    gone from it, and what is left is an oligo-length bound with a remedy in
    the spec rather than a dead library.
    """
    raw = in_tail_spec(dna(31, 17))
    at_default = cb.plan(raw, sources=SourceResolver(), offline=True)
    codes = {d.code for d in at_default.diagnostics}
    assert "in_tail_over_homology" in codes
    assert "no_footprint" not in codes, "nothing anneals to a bought fragment"
    assert [d.constraint.constraint for d in at_default.errors] == ["primer_len_max"]
    assert at_default.variants[0].ordered_ids == ("k",)

    raw["constraints"] = {"primer_len_max": 100}
    with_budget = cb.plan(
        raw, sources=SourceResolver(), offline=True, include_designs=True
    )
    assert with_budget.feasible, [d.message for d in with_budget.errors]
    assert with_budget.variants[0].verification.verdict == "verified"


def test_a_part_exactly_the_length_of_the_homology_still_rides_in_the_tails():
    """Budget zero is a real junction, so the guard is strictly greater-than.

    At 30 bp into a 30 bp homology the two neighbours supply none of it: both
    amplicons carry the part in full and the homology *is* the part.  That
    assembles, so declining it would be inventing a floor on how much of the
    overlap the neighbours have to own.

    ``primer_len_max`` is raised because a 30 bp tail plus a footprint is past
    teemi's default 60 -- which is the cap FLOW already sends 90 for, and is a
    fact about oligo length rather than about the budget being zero.
    """
    raw = in_tail_spec(dna(30, 18))
    raw["constraints"] = {"primer_len_max": 90}

    validated = load_spec(raw)
    records = SourceResolver().resolve_spec(validated)
    layout = plan_layout(validated, records)
    assert [s for s, _, mode in layout.slots if mode == "in_tail"] == ["at_k"]
    assert layout.diagnostics == (), "nothing to decline"

    result = cb.plan(raw, sources=SourceResolver(), offline=True)
    assert result.feasible, [d.message for d in result.errors]
    assert result.variants[0].fragment_ids == ("a", "b"), "still one fewer PCR"
    assert len(result.pcrs) == 2


def test_an_in_tail_part_needs_no_footprint_of_its_own():
    """A 9 bp Kozak is no template, so nothing anneals to it.

    Asking for a footprint anyway refused the part for being too short to do
    something it was never going to do.
    """
    validated = load_spec(in_tail_spec())
    records = SourceResolver().resolve_spec(validated)
    layout = plan_layout(validated, records)
    assert [s for s, _, mode in layout.slots if mode == "in_tail"] == ["at_k"]
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
# ordered: the one that changes what is BOUGHT
# --------------------------------------------------------------------------- #


def ordered_spec(middle, **constraints):
    raw = spec(
        [
            slot("a", dna(400, 13)),
            slot("k", middle, ordered=True),
            slot("b", dna(400, 14)),
        ]
    )
    raw["constraints"] = {"primer_len_max": 100, **constraints}
    return raw


def test_an_ordered_part_is_bought_not_amplified_and_is_in_the_product():
    """The mode asked for directly, with no declined request behind it.

    A 54 bp part no primer pair will amplify, written into the spec as
    ``ordered`` rather than arrived at by demotion: one fewer PCR than parts,
    no primers of its own, and the sequence in its place in the product.
    """
    middle = dna(54, 21)
    result = cb.plan(
        ordered_spec(middle),
        sources=SourceResolver(),
        offline=True,
        include_designs=True,
    )
    assert result.feasible, [d.message for d in result.errors]
    assert not [
        d
        for d in result.diagnostics
        if d.code.startswith("in_tail") or d.code == "no_footprint"
    ], "asked for, so there is no declined request and no missing footprint"

    variant = result.variants[0]
    assert variant.fragment_ids == ("a", "b")
    assert variant.ordered_ids == ("k",)
    assert len(result.pcrs) == 2

    fragment = result.ordered[0]
    assert (fragment.part_id, fragment.length_bp, fragment.because) == (
        "k",
        54,
        "asked",
    )
    assert fragment.sequence == middle
    assert fragment.in_constructs == 1

    product = str(result.designs[0].construct.seq)
    assert len(product) == variant.length_bp == 400 + 54 + 400
    assert product[400:454] == middle, "and in the right place"


def test_the_fold_order_boundary_goes_both_ways():
    """One base either side of ``homology_length`` picks a different purchase.

    At exactly the homology the part still folds: the budget is zero, both
    amplicons carry it in full, and the homology IS the part.  One base more
    and there is nothing left to spend it out of, so it is bought instead.
    This is the only place the two modes meet, and the whole bill changes
    across it -- two oligos either side, or two oligos and a fragment.
    """
    at_limit = cb.plan(
        in_tail_spec(dna(30, 18)),
        sources=SourceResolver(),
        offline=True,
    )
    assert at_limit.variants[0].ordered_ids == (), "folded"
    assert at_limit.ordered == ()
    assert at_limit.junctions[0].interposed == dna(30, 18)
    assert at_limit.junctions[0].ordered == ""

    over = cb.plan(in_tail_spec(dna(31, 18)), sources=SourceResolver(), offline=True)
    assert over.variants[0].ordered_ids == ("k",), "bought"
    assert [f.length_bp for f in over.ordered] == [31]
    assert over.junctions[0].interposed == ""
    assert over.junctions[0].ordered == dna(31, 18)


def test_an_ordered_junction_is_pydnas_own_geometry_base_for_base():
    """Not an imitation of ``assembly_fragments``: its actual result.

    pydna documents what happens to an ``Amplicon`` next to a ``Dseqrecord``
    it is not allowed to absorb -- the tail is the whole ``overlap``, because
    the recombining sequence sits entirely on that one primer -- and that is
    the manoeuvre this mode is.  So the fragment is handed to pydna with a
    ``maxlink`` below its own length, and the two primers that come back are
    what teemi's two primers have to equal.
    """
    from pydna.amplify import pcr
    from pydna.design import assembly_fragments
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer

    middle = dna(54, 21)
    a, b = dna(400, 13), dna(400, 14)
    result = cb.plan(ordered_spec(middle), sources=SourceResolver(), offline=True)
    junction = result.junctions[0]
    assert junction.designed_by == "pydna"
    assert junction.ordered == middle
    assert junction.homology == middle, (
        "the two joins are its first and last L bases, so the fragment IS what "
        "there is to report; there is no single shared stretch between the "
        "amplicons to measure instead"
    )

    footprints = {f.part_id: f for f in result.footprints}
    amplicons = [
        pcr(
            Primer(footprints[part].forward),
            Primer(footprints[part].reverse),
            Dseqrecord(sequence),
        )
        for part, sequence in (("a", a), ("b", b))
    ]
    homology = result.spec.primer_design.junction.homology_length
    tailed = assembly_fragments(
        [amplicons[0], Dseqrecord(middle), amplicons[1]],
        overlap=homology,
        maxlink=len(middle) - 1,
    )
    assert len(tailed) == 3, "pydna leaves the fragment standing"
    assert junction.chosen.forward_primer == str(tailed[-1].forward_primer.seq)
    assert junction.chosen.reverse_primer == str(tailed[0].reverse_primer.seq)

    # And the tails are the full homology taken from the fragment, each side.
    assert junction.chosen.forward_primer.startswith(middle[-homology:])
    assert junction.chosen.reverse_primer.startswith(
        str(Seq(middle[:homology]).reverse_complement())
    )
    assert junction.chosen.t_f == junction.chosen.t_r == homology


def test_the_search_is_not_run_on_an_ordered_junction():
    """``tune_overhangs`` has nothing to tune: the fragment pins both tails.

    Sweeping would offer splits of a budget that is not being split, and the
    cheapest of them would be a primer pair that does not assemble.
    """
    raw = ordered_spec(dna(54, 21))
    raw["primer_design"]["junction"]["tune_overhangs"] = True
    result = cb.plan(raw, sources=SourceResolver(), offline=True)
    junction = result.junctions[0]
    assert junction.designed_by == "pydna", "not 'sweep', even though asked"
    assert len(junction.profile) == 1
    assert junction.unavoidable == ()


def test_an_in_tail_slot_next_to_an_ordered_one_joins_the_fragment():
    """Two purchases next to each other are one purchase.

    A Kozak asked into the tails, sitting against a bought signal peptide: the
    two are adjacent in the construct, so the buildable reading is one 63 bp
    fragment rather than a 9 bp tail on top of a 30 bp one.  Two separate
    synthetic pieces would share no sequence and could not be joined at all.
    """
    peptide = dna(54, 21)
    raw = spec(
        [
            slot("a", dna(400, 13)),
            slot("koz", KOZAK, in_tail=True),
            slot("k", peptide, ordered=True),
            slot("b", dna(400, 14)),
        ]
    )
    raw["constraints"] = {"primer_len_max": 100}

    result = cb.plan(
        raw, sources=SourceResolver(), offline=True, include_designs=True
    )
    assert result.feasible, [d.message for d in result.errors]
    note = next(
        d for d in result.diagnostics if d.code == "in_tail_folded_into_ordered"
    )
    assert note.severity == "warning" and note.part_id == "koz"
    assert "k" in note.message

    variant = result.variants[0]
    assert variant.fragment_ids == ("a", "b")
    assert variant.ordered_ids == ("koz", "k"), "both bought, in order"
    assert result.junctions[0].ordered == KOZAK + peptide, "as one piece"
    assert len(result.pcrs) == 2

    product = str(result.designs[0].construct.seq)
    assert product == dna(400, 13) + KOZAK + peptide + dna(400, 14)


def test_an_ordered_fragment_works_in_a_circular_design():
    """The last fragment joins the first, so a bought part has two tails there too.

    Checked separately because the verification is rotation-invariant for a
    circular construct and the junction list wraps, which is where an
    off-by-one in the run grouping would hide.
    """
    middle = dna(54, 21)
    raw = spec(
        [
            slot("a", dna(400, 13)),
            slot("k", middle, ordered=True),
            slot("b", dna(400, 14)),
        ],
        topology="circular",
    )
    raw["constraints"] = {"primer_len_max": 100}

    result = cb.plan(
        raw, sources=SourceResolver(), offline=True, include_designs=True
    )
    assert result.feasible, [d.message for d in result.errors]
    assert result.variants[0].verification.verdict == "verified"
    assert result.variants[0].ordered_ids == ("k",)
    assert len(result.pcrs) == 2
    assert len(str(result.designs[0].construct.seq)) == 400 + 54 + 400


def test_ordered_on_a_combinatorial_position_is_refused():
    """Same reason as in_tail: the junction would depend on a non-fragment."""
    with pytest.raises(ValidationError, match="must have exactly one part"):
        load_spec(
            spec(
                [
                    slot("a", dna(400, 16)),
                    {
                        "id": "at_k",
                        "ordered": True,
                        "parts": [
                            {
                                "id": "k1",
                                "source": {"type": "inline", "sequence": dna(54, 1)},
                            },
                            {
                                "id": "k2",
                                "source": {"type": "inline", "sequence": dna(54, 2)},
                            },
                        ],
                    },
                    slot("b", dna(400, 17)),
                ]
            )
        )


def test_ordered_at_the_end_of_a_linear_design_is_refused():
    with pytest.raises(ValidationError, match="no neighbour on both sides"):
        load_spec(
            spec(
                [
                    slot("a", dna(400, 18)),
                    slot("b", dna(400, 19)),
                    slot("k", dna(54, 1), ordered=True),
                ]
            )
        )


def test_in_tail_and_ordered_together_are_refused():
    """One says fold it into the tails, the other says do not."""
    with pytest.raises(ValidationError, match="cannot both hold"):
        load_spec(
            spec(
                [
                    slot("a", dna(400, 18)),
                    slot("k", dna(54, 1), in_tail=True, ordered=True),
                    slot("b", dna(400, 19)),
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
