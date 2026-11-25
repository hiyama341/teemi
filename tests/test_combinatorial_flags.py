#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)

"""The slot flags, which come from FLOW's ``construct_part``.

They are *slot-level*: they describe the place on the DNA, so every alternative
at one position shares them.  Three change the sequence (``reverse``,
``drop_start``, ``drop_stop``), two pin a junction's split (``no_tail_fwd``,
``no_tail_rev``) and two say how the slot reaches the product without a PCR of
its own -- ``in_tail``, which folds it into the flanking primers, and
``ordered``, which buys it as a synthetic double-stranded fragment.

Neither is ever chosen for the caller.  A fold that will not fit inside
``primer_len_max`` is reported as an error naming the remedies, one of which is
``ordered``; it is not quietly turned into a purchase.
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

    ``t_f`` and ``t_r`` are WHOLE TAILS, so their total is ``L + P`` and not
    ``L``: each one reaches back across its share of the part and on into the
    neighbour beyond it, and between them they span the part once and the
    homology once.  That total is the same for every split -- which is what
    :meth:`DesignSpec.in_tail_limit` is derived from -- so the split moves
    length between the two primers and can never remove any.
    """
    with_it = cb.plan(in_tail_spec(), sources=SourceResolver(), offline=True)
    junction = with_it.junctions[0]
    assert junction.interposed == KOZAK
    assert len(junction.homology) == 30, "exactly what was asked for"
    assert KOZAK in junction.homology, "and it is in there"
    assert junction.chosen.t_f + junction.chosen.t_r == 30 + len(KOZAK)

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


def test_a_part_past_an_explicit_in_tail_cap_is_reported_not_bought():
    """FLOW's schema says as much about teemi's maxlink: a request, not a command.

    But the answer to a declined request is a REPORT, not a purchase.  Buying
    DNA costs money and takes a fortnight, so it is the caller's decision and
    nobody else's: the slot stays folded, the error names the cap that bit and
    the three ways out of it, and the over-length primer it leaves behind is
    reported too.  Choosing ``ordered`` is one of the three, which is how this
    stays a two-line fix for a caller who does want to buy it.

    The homology is 60 here so the derived bound is not what bites: this test
    is about the cap the caller set, and the message names that number because
    raising it is the remedy.
    """
    raw = in_tail_spec(dna(50, 15))
    raw["primer_design"]["junction"]["homology_length"] = 60
    raw["constraints"] = {"primer_len_max": 100, "in_tail_max_bp": 40}

    result = cb.plan(raw, sources=SourceResolver(), offline=True)
    note = next(d for d in result.diagnostics if d.code == "in_tail_too_long")
    assert note.severity == "error", "a refusal to design, not a remark"
    assert "in_tail_max_bp 40" in note.message
    assert "still folded into the tails" in note.message
    assert "set ordered on at_k" in note.message, "the remedy, offered"

    variant = result.variants[0]
    assert variant.fragment_ids == ("a", "b"), "still two PCRs, not three"
    assert variant.ordered_ids == (), "NOTHING is bought behind the caller"
    assert result.ordered == ()
    assert len(result.pcrs) == 2
    assert result.junctions[0].interposed == dna(50, 15), "still in the tails"


def test_a_part_past_what_the_primers_can_carry_names_the_cap_to_raise():
    """The limit a constant cannot express, and the one 40 bp gets wrong twice.

    The two tailed primers have an invariant total length -- ``P + L + f_fwd +
    f_rev``, with no ``t_f`` in it -- so no split makes the longer of them
    shorter than half of that, and half of that has to clear
    ``primer_len_max``.  At a 60 bp homology and the default 60 nt cap that
    leaves 24 bp, so a 50 bp part provably cannot fold however the split falls.

    ``in_tail_max_bp`` is ``None``, which is how a caller says "use the derived
    bound and stop guessing".  The message gives the number, where it came from
    and the ``primer_len_max`` that WOULD carry it, because the remedy is a
    different field from the one the explicit cap points at.
    """
    raw = in_tail_spec(dna(50, 15))
    raw["primer_design"]["junction"]["homology_length"] = 60
    raw["constraints"] = {"in_tail_max_bp": None}
    validated = load_spec(raw)
    assert validated.in_tail_limit() == 2 * 60 - 60 - 2 * 18 == 24

    result = cb.plan(raw, sources=SourceResolver(), offline=True)
    note = next(d for d in result.diagnostics if d.code == "in_tail_too_long")
    assert note.severity == "error"
    assert "24 bp the two tails can carry" in note.message
    assert "in_tail_max_bp" not in note.message, "not the cap that bit"
    assert "invariant total" in note.message
    assert "primer_len_max to 73" in note.message, "the number that would work"
    assert result.variants[0].ordered_ids == ()

    # And the remedy the message names is the remedy: at that cap it folds.
    raw["constraints"] = {"primer_len_max": 73, "in_tail_max_bp": None}
    fixed = cb.plan(raw, sources=SourceResolver(), offline=True)
    assert [d.code for d in fixed.diagnostics if d.code == "in_tail_too_long"] == []


def test_the_derived_bound_overrides_a_cap_set_above_it():
    """A cap above what the oligos carry is not a permission.

    40 bp is the historical default and is *too permissive* at a 50 bp
    homology: the derived bound there is 34, so a 38 bp part used to fold
    quietly and then breach ``primer_len_max`` on a primer nobody asked for.
    It still folds -- that is the caller's request and this function does not
    overrule it -- but now the error arrives at layout time with the arithmetic
    in it, instead of as a mystery oligo three tables later.
    """
    raw = in_tail_spec(dna(38, 15))
    raw["primer_design"]["junction"]["homology_length"] = 50
    raw["constraints"] = {"in_tail_max_bp": 40}
    validated = load_spec(raw)
    assert validated.in_tail_limit() == 34, "the derived bound, not the 40 asked"

    result = cb.plan(raw, sources=SourceResolver(), offline=True)
    note = next(d for d in result.diagnostics if d.code == "in_tail_too_long")
    assert "in_tail_max_bp" not in note.message, "the derived bound is what bit"
    assert "34 bp the two tails can carry" in note.message
    assert result.variants[0].ordered_ids == ()


def test_the_derived_bound_also_allows_a_fold_a_constant_refused():
    """The same constant is *too conservative* in the other direction.

    At a 45 bp homology with ``primer_len_max`` 100 the oligos can carry 119 bp,
    so a 42 bp part folds into two primers well inside the cap.  The default 40
    refused it, which is the cost of a constant where the truth is a function;
    ``in_tail_max_bp: null`` is how to stop paying it.
    """
    raw = in_tail_spec(dna(42, 15))
    raw["primer_design"]["junction"]["homology_length"] = 45
    raw["constraints"] = {"primer_len_max": 100, "in_tail_max_bp": None}
    validated = load_spec(raw)
    assert validated.in_tail_limit() == 2 * 100 - 45 - 2 * 18 == 119

    result = cb.plan(
        raw, sources=SourceResolver(), offline=True, include_designs=True
    )
    assert result.feasible, [d.message for d in result.errors]
    assert result.variants[0].ordered_ids == (), "folded, and nothing bought"
    assert result.ordered == ()
    assert result.variants[0].verification.verdict == "verified"

    refused = dict(raw)
    refused["constraints"] = {"primer_len_max": 100, "in_tail_max_bp": 40}
    note = next(
        d
        for d in cb.plan(refused, sources=SourceResolver(), offline=True).diagnostics
        if d.code == "in_tail_too_long"
    )
    assert "in_tail_max_bp 40" in note.message, "the constant refuses what fits"


def test_the_fold_bound_is_a_floor_no_split_beats_at_any_length():
    """``ceil((P + L + f_fwd + f_rev) / 2)`` and nothing else, at every length.

    The old accounting had a special case at ``L == P``: the split had one
    value there, so the longest oligo was ``P + max(f_fwd, f_rev)`` exactly.
    Under whole tails there is no such point -- the split is free at every
    length -- so one formula covers the lot, and what it gives is a FLOOR.  No
    split can come in under it, because the two tails' total does not contain
    the split at all; a search can only move length from one primer to the
    other.

    ``P`` runs from a third of the homology to four times it, straddling the
    ``L == P`` that used to be a boundary.  This is the inequality
    :meth:`DesignSpec.in_tail_limit` is derived from, so it is measured rather
    than assumed.
    """
    import math

    for part in (9, 30, 54, 68, 93, 120):
        raw = in_tail_spec(dna(part, 15))
        raw["primer_design"]["junction"].update(
            homology_length=30,
            tune_overhangs=True,
            split_search={"order": "exhaustive", "accept": "argmin"},
        )
        raw["constraints"] = {"primer_len_max": 200, "in_tail_max_bp": None}

        result = cb.plan(raw, sources=SourceResolver(), offline=True)
        assert result.variants[0].ordered_ids == (), "folds, at every length"
        footprints = {f.part_id: f for f in result.footprints}
        floor = math.ceil(
            (part + 30 + len(footprints["b"].forward) + len(footprints["a"].reverse))
            / 2
        )
        longest = max(len(p.sequence) for p in result.primers)
        assert longest >= floor, (part, longest, floor)


def test_the_fold_bound_is_reached_exactly_and_one_base_below_it_nothing_fits():
    """The floor is attainable, which is what makes it the right bound.

    A 60 bp part at 30 bp of homology between an 18 nt and a 32 nt footprint
    has a floor of ``ceil((60 + 30 + 18 + 32) / 2) = 70``.  Set
    ``primer_len_max`` to exactly 70 and precisely one split of the 90 bp of
    tail is feasible -- the balancing one -- so feasibility pins it with no
    flag needed, and the design verifies on a 70 nt longest oligo.  Set it to
    69 and nothing fits, which is the bound being tight rather than cautious.

    The balancing split also puts the whole 30 bp homology INSIDE the 60 bp
    part: both amplicons carry the part in full and share 30 bp of it, and
    neither neighbour contributes a base of the overlap.  That is the geometry
    the old ``homology_length`` check denied existed.
    """
    part = dna(60, 15)
    raw = in_tail_spec(part)
    raw["primer_design"]["junction"].update(
        tune_overhangs=True,
        split_search={"order": "exhaustive", "accept": "argmin"},
    )
    raw["constraints"] = {"primer_len_max": 70, "in_tail_max_bp": None}

    at_floor = cb.plan(raw, sources=SourceResolver(), offline=True)
    assert at_floor.feasible, [d.message for d in at_floor.errors]
    assert at_floor.variants[0].verification.verdict == "verified"
    assert max(len(pr.sequence) for pr in at_floor.primers) == 70
    junction = at_floor.junctions[0]
    assert junction.chosen.t_f == 52, "the balancing split, pinned by the cap"
    assert junction.homology in part, "the overlap is 30 bp OF THE PART"
    assert len(junction.homology) == 30

    raw["constraints"] = {"primer_len_max": 69, "in_tail_max_bp": None}
    below = cb.plan(raw, sources=SourceResolver(), offline=True)
    assert not below.feasible
    assert [d.constraint.constraint for d in below.errors] == ["primer_len_max"]


def test_a_part_longer_than_the_homology_folds_with_the_homology_inside_it():
    """THE DEFECT.  A part longer than the overlap folds; it used to be refused.

    The old accounting read ``homology_length`` as a budget the part was spent
    out of, so a 60 bp part at 30 bp of homology left the two neighbours
    ``30 - 60 = -30`` bp to supply and the junction looked impossible.  It is
    not, and the negative number is the tell: **the overlap does not have to
    contain the part.**  Once the part is the longer of the two the homology
    sits inside THE PART -- both amplicons carry the whole part, they share 30
    bp of the part itself, and neither neighbour contributes a base of the
    overlap.  That is what the last assertion here measures.

    pydna's ``assembly_fragments`` has always built exactly this.  The check in
    front of it was the only thing saying no.
    """
    raw = in_tail_spec(dna(60, 15))
    raw["constraints"] = {"primer_len_max": 100, "in_tail_max_bp": None}

    result = cb.plan(
        raw, sources=SourceResolver(), offline=True, include_designs=True
    )
    assert result.feasible, [d.message for d in result.errors]
    codes = [d.code for d in result.diagnostics]
    assert "in_tail_too_long" not in codes, "100 nt of cap carries 60 bp"
    assert "in_tail_over_homology" not in codes, "that check is gone, and was wrong"
    assert "impossible_junction" not in codes, "it used to end here"

    variant = result.variants[0]
    assert variant.verification.verdict == "verified"
    assert variant.fragment_ids == ("a", "b") and variant.ordered_ids == ()
    assert len(result.pcrs) == 2
    assert str(result.designs[0].construct.seq) == (
        dna(400, 13) + dna(60, 15) + dna(400, 14)
    )

    # The 30 bp overlap is a window of the construct that SLIDES with the
    # split -- it is wherever the two amplicons happen to meet -- and the
    # 60 bp part is simply longer than it.  Neither has to contain the other,
    # which is the thing the old budget could not express: it insisted the
    # overlap contain the part and called the shortfall impossible.
    junction = result.junctions[0]
    whole = dna(400, 13) + dna(60, 15) + dna(400, 14)
    start = 400 + 60 - junction.chosen.t_f
    assert junction.homology == whole[start : start + 30]
    assert junction.chosen.t_f + junction.chosen.t_r == 30 + 60
    assert 400 < start + 30 and start < 460, "and it overlaps the part"


def test_the_case_that_used_to_lose_the_whole_library():
    """A 54 bp signal peptide at 30 bp of homology: ANNEAL's own failure.

    54 bp is past a 30 bp homology, so the old accounting would not fold it;
    demoting it to a fragment of its own produced primers pydna then refused to
    amplify --

        No PCR product! Template 54 bp linear limit=13: No forward primers
        anneal... No reverse primers anneal

    -- so the variant came back ``no_product`` and the whole design was
    infeasible.  There was no third answer and both of the two were impossible.

    The fix is not a third answer: it is that the FIRST one always worked.  The
    peptide folds into the two flanking primers, which is where a signal
    peptide belongs -- it is 18 codons, not a part to buy -- and the longest
    oligo it costs is 61 nt.  ``primer_len_max`` is 100 because that is what
    ANNEAL sends.
    """
    middle = dna(54, 21)
    raw = in_tail_spec(middle)
    raw["constraints"] = {"primer_len_max": 100, "in_tail_max_bp": None}

    result = cb.plan(
        raw, sources=SourceResolver(), offline=True, include_designs=True
    )
    assert result.feasible, [d.message for d in result.errors]
    variant = result.variants[0]
    assert variant.verification.verdict == "verified"
    assert variant.fragment_ids == ("a", "b")
    assert variant.ordered_ids == (), "folded into the tails, not bought"
    assert result.ordered == ()
    assert len(result.pcrs) == 2

    # The two tails are the construct either side of the join, base for base,
    # and between them they cover every base of the peptide -- which is what
    # "it is in the primers" has to mean if the product is to contain it.
    a, b = dna(400, 13), dna(400, 14)
    footprints = {f.part_id: f for f in result.footprints}
    t_f, t_r = result.junctions[0].chosen.t_f, result.junctions[0].chosen.t_r
    assert t_f + t_r == 30 + 54, "the part once, the homology once"
    assert result.junctions[0].chosen.forward_primer == (
        (a + middle)[-t_f:] + footprints["b"].forward
    )
    assert result.junctions[0].chosen.reverse_primer == (
        str(Seq((middle + b)[:t_r]).reverse_complement()) + footprints["a"].reverse
    )
    covered = set(range(len(middle) - min(t_f, len(middle)), len(middle)))
    covered |= set(range(min(t_r, len(middle))))
    assert covered == set(range(len(middle))), "every base of it, in an oligo"
    assert max(len(pr.sequence) for pr in result.primers) <= 100

    product = str(result.designs[0].construct.seq)
    assert product == dna(400, 13) + middle + dna(400, 14)




def test_a_part_exactly_the_length_of_the_homology_still_rides_in_the_tails():
    """``L == P`` is not a boundary of any kind, and used to be the last one.

    At 30 bp into a 30 bp homology the two neighbours supply none of it: both
    amplicons carry the part in full and the homology *is* the part.  The old
    accounting called that a budget of zero and allowed it only because the
    guard was strictly greater-than; one base more was refused.  There is no
    budget now and no guard, so this length is unremarkable -- which is the
    point of asserting it still works.

    ``primer_len_max`` is raised because a 30 bp part in both tails plus a
    footprint is past teemi's default 60, and that is a fact about oligo
    length rather than about any boundary here.
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
    assert result.variants[0].ordered_ids == ()
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


def test_an_ordered_fragments_tails_do_not_grow_with_it():
    """The mode has no upper bound of its own, which is its whole appeal.

    Each neighbour takes a tail of exactly ``homology_length`` from the
    fragment whatever its length, so the oligos do not grow with it: 276 bp
    costs the same two tails as 54 bp.  That is the one thing folding cannot
    do -- a folded part is IN the primers, so there it is the part's length
    that drives them -- and it is why ``ordered`` is the right answer for
    something genuinely long, and the wrong one for an 18-codon peptide.
    """
    longest = {}
    for length in (54, 93, 276):
        middle = dna(length, 22)
        result = cb.plan(
            ordered_spec(middle),
            sources=SourceResolver(),
            offline=True,
            include_designs=True,
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

    A folded part splits its two tails; a bought one is a molecule of its own,
    so each neighbour carries the WHOLE homology on top of its footprint.  At
    teemi's default ``primer_len_max`` of 60 a 30 bp tail on a 32 nt footprint
    is 62 nt and breaches it -- which is why ANNEAL sends 100, and why a
    caller who reaches for ``ordered`` to dodge an over-length primer can find
    it waiting for them on the other side.
    """
    raw = ordered_spec(dna(31, 17))
    raw["constraints"] = {}  # teemi's own default cap of 60
    at_default = cb.plan(raw, sources=SourceResolver(), offline=True)
    codes = {d.code for d in at_default.diagnostics}
    assert "no_footprint" not in codes, "nothing anneals to a bought fragment"
    assert [d.constraint.constraint for d in at_default.errors] == ["primer_len_max"]
    assert at_default.variants[0].ordered_ids == ("k",)

    raw["constraints"] = {"primer_len_max": 100}
    with_budget = cb.plan(
        raw, sources=SourceResolver(), offline=True, include_designs=True
    )
    assert with_budget.feasible, [d.message for d in with_budget.errors]
    assert with_budget.variants[0].verification.verdict == "verified"


def test_fold_or_order_is_the_callers_choice_and_not_a_length():
    """Nothing switches mode on its own, so the bill is what the spec asked for.

    The same 54 bp part, written two ways: ``in_tail`` folds it into two
    oligos, ``ordered`` buys it as a fragment and gives each neighbour a full
    homology tail.  Both build.  There is no length at which teemi picks for
    the caller -- that boundary used to sit at ``homology_length`` and was the
    defect -- so the only thing that decides which bill arrives is the flag.
    """
    middle = dna(54, 21)

    folded_raw = in_tail_spec(middle)
    folded_raw["constraints"] = {"primer_len_max": 100, "in_tail_max_bp": None}
    folded = cb.plan(folded_raw, sources=SourceResolver(), offline=True)
    assert folded.variants[0].ordered_ids == ()
    assert folded.ordered == (), "two oligos, and no purchase"
    assert folded.junctions[0].interposed == middle
    assert folded.junctions[0].ordered == ""

    bought = cb.plan(ordered_spec(middle), sources=SourceResolver(), offline=True)
    assert bought.variants[0].ordered_ids == ("k",)
    assert [f.length_bp for f in bought.ordered] == [54]
    assert bought.junctions[0].interposed == ""
    assert bought.junctions[0].ordered == middle

    # Same part, same two PCRs, different purchase -- and the folded one's
    # oligos are the longer, because they contain the part.
    assert len(folded.pcrs) == len(bought.pcrs) == 2
    assert max(len(pr.sequence) for pr in folded.primers) > max(
        len(pr.sequence) for pr in bought.primers
    )


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
        # 0, not ``len(middle) - 1``: both leave the fragment standing, but
        # pydna absorbs every fragment with ``len(f) <= maxlink`` whatever its
        # type, so anything larger than it needs can swallow a short
        # NEIGHBOUR instead.  This is the maxlink the module passes.
        maxlink=0,
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


def test_a_neighbour_shorter_than_the_fragment_is_still_amplified():
    """The fragment's length is not a licence to absorb anything else.

    ``maxlink`` has to be as small as says "leave this fragment standing",
    because ``assembly_fragments`` absorbs every fragment with
    ``len(f) <= maxlink`` whatever its type.  A maxlink read off the fragment
    -- ``len(fragment) - 1``, which reads naturally and was wrong -- is larger
    than a short NEIGHBOUR, so a 276 bp purchase beside a 200 bp part folds
    that part into its own neighbours and leaves it with two bare footprints
    and no PCR.
    """
    fragment = dna(276, 22)
    raw = spec(
        [
            slot("a", dna(400, 13)),
            slot("k", fragment, ordered=True),
            slot("short_b", dna(200, 19)),
            slot("c", dna(400, 14)),
        ]
    )
    raw["constraints"] = {"primer_len_max": 100}

    result = cb.plan(
        raw, sources=SourceResolver(), offline=True, include_designs=True
    )
    assert result.variants[0].verification.verdict == "verified"
    assert result.variants[0].fragment_ids == ("a", "short_b", "c")
    assert result.variants[0].ordered_ids == ("k",)
    assert {pcr.part_id for pcr in result.pcrs} == {"a", "short_b", "c"}
    mine = [pr for pr in result.primers if "short_b" in pr.part_ids]
    assert len(mine) == 2
    assert all(len(pr.sequence) > len(pr.footprint) for pr in mine)
    assert str(result.designs[0].construct.seq) == (
        dna(400, 13) + fragment + dna(200, 19) + dna(400, 14)
    )


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
        maxlink=result.spec.in_tail_limit(),
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
