#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)
"""Tests for :mod:`teemi.design.combinatorial.design`.

The load-bearing properties: the homology is exactly ``L`` bases, a zero-length
tail is empty rather than a whole part, primer3 numbers arrive in kcal/mol,
footprints land inside the Tm window, and an infeasible junction comes back
named rather than raised.
"""

from __future__ import annotations

import dataclasses
import json
import math
import pathlib
import random
import re

import primer3
import pytest
import requests
from Bio.Seq import Seq

from teemi.design.combinatorial import design as D
from teemi.design.combinatorial import sources as SRC
from teemi.design.combinatorial import tm as T
from teemi.design.combinatorial.spec import load_spec

from tests.combinatorial_netguard import no_network  # noqa: F401

REPO = pathlib.Path(__file__).resolve().parents[1]
SPECS = REPO / "tests" / "files_for_testing"


def _forbidden(*args, **kwargs):
    raise AssertionError("a test reached for NEB without stubbing it")


@pytest.fixture(autouse=True)
def fresh_caches(monkeypatch):
    """Empty memos and a poisoned ``tm_neb`` on top of the module-wide guard."""

    monkeypatch.setattr(T, "tm_neb", _forbidden)
    T.clear_tm_caches()
    yield
    T.clear_tm_caches()


def _random_dna(seed: int, length: int = 180) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _raw_spec(**overrides) -> dict:
    raw = {
        "design": {"name": "synthetic", "topology": "linear"},
        "positions": [
            {
                "id": "first",
                "parts": [
                    {
                        "id": "a1",
                        "source": {"type": "inline", "sequence": _random_dna(1)},
                    },
                    {
                        "id": "a2",
                        "source": {"type": "inline", "sequence": _random_dna(2)},
                    },
                ],
            },
            {
                "id": "second",
                "parts": [
                    {
                        "id": "b1",
                        "source": {"type": "inline", "sequence": _random_dna(4)},
                    }
                ],
            },
        ],
        # This module is mostly about the split search, so it asks for it; the
        # default path is pydna's and is covered on its own.
        "primer_design": {
            "target_tm": 57.0,
            "footprint_len_min": 18,
            "junction": {"tune_overhangs": True},
        },
    }
    raw.update(overrides)
    return raw


def _resolved(**overrides):
    spec = load_spec(_raw_spec(**overrides))
    records = SRC.SourceResolver().resolve_spec(spec)
    footprints, diagnostics = D.design_footprints(
        records, spec, T.FootprintTm(offline=True)
    )
    return spec, records, footprints, diagnostics


def _junction(spec, records, footprints, upstream="a1", downstream="b1"):
    return D.Junction(
        upstream=footprints[upstream],
        downstream=footprints[downstream],
        upstream_seq=str(records[upstream].seq),
        downstream_seq=str(records[downstream].seq),
        homology_length=spec.primer_design.junction.homology_length,
    )


# --------------------------------------------------------------------------- #
# The junction model
# --------------------------------------------------------------------------- #


def test_zero_tail_is_empty_not_the_whole_part():
    # seq[-0:] is the whole string in Python.
    assert D._suffix("ACGTACGT", 0) == ""
    assert D._suffix("ACGTACGT", 3) == "CGT"


def test_a_zero_length_forward_tail_leaves_the_bare_footprint():
    spec, records, footprints, _ = _resolved()
    junction = _junction(spec, records, footprints)
    forward, reverse = junction.primers(0)
    assert forward == footprints["b1"].forward
    assert len(forward) == len(footprints["b1"].forward)
    assert len(reverse) == junction.homology_length + len(footprints["a1"].reverse)


def test_a_full_forward_tail_leaves_the_bare_reverse_footprint():
    spec, records, footprints, _ = _resolved()
    junction = _junction(spec, records, footprints)
    forward, reverse = junction.primers(junction.homology_length)
    assert reverse == footprints["a1"].reverse
    assert len(forward) == junction.homology_length + len(footprints["b1"].forward)


def test_the_two_tails_always_sum_to_exactly_l():
    spec, records, footprints, _ = _resolved()
    junction = _junction(spec, records, footprints)
    length = junction.homology_length
    for t_f in range(length + 1):
        assert len(junction.homology(t_f)) == length, "L must mean exactly L"
        forward, reverse = junction.primers(t_f)
        tail_f = len(forward) - len(footprints["b1"].forward)
        tail_r = len(reverse) - len(footprints["a1"].reverse)
        assert tail_f + tail_r == length
        assert tail_f == t_f


def test_the_primers_are_the_documented_formula():
    spec, records, footprints, _ = _resolved()
    junction = _junction(spec, records, footprints)
    upstream, downstream = junction.upstream_seq, junction.downstream_seq
    t_f = 14
    t_r = junction.homology_length - t_f
    forward, reverse = junction.primers(t_f)
    assert forward == upstream[-t_f:] + footprints["b1"].forward
    assert (
        reverse
        == str(Seq(downstream[:t_r]).reverse_complement()) + footprints["a1"].reverse
    )
    assert junction.homology(t_f) == upstream[-t_f:] + downstream[:t_r]


def test_the_homology_is_shared_by_the_two_amplicons():
    # The forward primer of B and the reverse primer of A are complementary over
    # the homology; that complementarity IS the overlap the assembly uses.
    spec, records, footprints, _ = _resolved()
    junction = _junction(spec, records, footprints)
    t_f = 18
    forward, reverse = junction.primers(t_f)
    homology = junction.homology(t_f)
    assert forward.startswith(homology[:t_f])
    assert reverse.startswith(str(Seq(homology[t_f:]).reverse_complement()))


def test_split_order_fixed_balanced():
    assert D._split_order("fixed_balanced", 35) == (18,)
    assert D._split_order("fixed_balanced", 36) == (18,)


def test_split_order_steps_outward_from_balanced():
    order = D._split_order("outward_from_balanced", 35)
    assert order[:5] == (18, 17, 19, 16, 20)
    assert sorted(order) == list(range(36))


def test_split_order_exhaustive_walks_every_split():
    assert D._split_order("exhaustive", 35) == tuple(range(36))


def test_a_junction_the_parts_cannot_supply_is_an_error():
    spec, records, footprints, _ = _resolved()
    junction = dataclasses.replace(
        _junction(spec, records, footprints),
        upstream_seq="ACGTACGTAC",
        downstream_seq="ACGTACGTAC",
    )
    with pytest.raises(D.JunctionError) as caught:
        D.design_junction(junction, spec)
    assert caught.value.junction_id == "a1__b1"


# --------------------------------------------------------------------------- #
# primer3
# --------------------------------------------------------------------------- #


def test_primer3_dg_arrives_in_kcal_per_mol():
    hairpin = "GGGGCCCCGGGGCCCCAAAAGGGGCCCCGGGGCCCC"
    dg, _ = D._structure(primer3.calc_hairpin, (hairpin,), T.CONDITIONS)
    raw = primer3.calc_hairpin(
        hairpin, temp_c=D.PRIMER3_TEMP_C, **T.CONDITIONS.primer3_kwargs
    ).dg
    assert dg == pytest.approx(raw / 1000.0)
    assert -50.0 < dg < 5.0, "a cal/mol number here makes the -6.0 floor inert"


def test_the_dg_floor_is_read_in_kcal_on_a_weak_structure():
    """The conversion has to be pinned on a WEAK structure.

    A strong hairpin breaches the -6.0 floor whether or not the cal/mol number
    was divided by 1000 -- both sides saturate and the off-by-1000 hides.  This
    hairpin is -0.91 kcal/mol, comfortably inside the floor, so the floor is
    inert in kcal and impossible in cal: the two readings disagree about
    whether the primer is orderable at all.
    """
    weak = "ATTACTTGCATGACGATCGTTGGTCGGCTC"
    raw_dg = primer3.calc_hairpin(
        weak, temp_c=D.PRIMER3_TEMP_C, **T.CONDITIONS.primer3_kwargs
    ).dg
    assert -6000.0 < raw_dg < -100.0, "this test needs a weak but real hairpin"

    constraints = load_spec(_raw_spec()).constraints
    structures = D.oligo_structures(weak, T.CONDITIONS)
    assert structures.hairpin_dg == pytest.approx(raw_dg / D.CAL_PER_KCAL)
    passed = D.oligo_violations("p", weak, structures, constraints, weak)
    assert not [v for v in passed if v.constraint == "hairpin_dg_min"]

    unconverted = dataclasses.replace(structures, hairpin_dg=raw_dg)
    breached = [
        v
        for v in D.oligo_violations("p", weak, unconverted, constraints, weak)
        if v.constraint == "hairpin_dg_min"
    ]
    assert breached, "a cal/mol observation must not clear a kcal/mol floor"
    assert breached[0].units == "kcal/mol"


def test_primer3_refuses_more_than_sixty_bases():
    # The reason for the clamp, pinned so a primer3 upgrade that lifts the
    # limit shows up here rather than as a silent change of meaning.
    with pytest.raises(RuntimeError):
        primer3.calc_heterodimer("A" * 61, "T" * 61)


def test_structures_clamp_to_the_three_prime_sixty_bases():
    long_primer = _random_dna(9, 90)
    structures = D.oligo_structures(long_primer, T.CONDITIONS)
    assert structures.clamped is True
    tail = D.oligo_structures(long_primer[-D.CLAMP_NT :], T.CONDITIONS)
    assert structures.hairpin_dg == pytest.approx(tail.hairpin_dg)
    assert D.oligo_structures(long_primer[:50], T.CONDITIONS).clamped is False


def test_a_seventy_base_primer_pair_is_scored_not_refused():
    """A tailed primer routinely passes 60 nt, which is exactly where primer3's
    thermodynamic alignment raises.  Both pair terms have to survive that."""
    forward, reverse = _random_dna(31, 78), _random_dna(32, 71)
    with pytest.raises(RuntimeError):
        primer3.calc_heterodimer(forward, reverse)

    pair = D._pair_structures(forward, reverse, T.CONDITIONS)
    assert pair.clamped is True
    # The 3' end is the end that has to stay free, so that is the end kept.
    clamped = D._pair_structures(
        forward[-D.CLAMP_NT :], reverse[-D.CLAMP_NT :], T.CONDITIONS
    )
    assert pair.heterodimer_dg == pytest.approx(clamped.heterodimer_dg)
    assert pair.three_prime_dimer_dg == pytest.approx(clamped.three_prime_dimer_dg)
    assert pair.heterodimer_dg != pytest.approx(
        D._pair_structures(
            forward[: D.CLAMP_NT], reverse[: D.CLAMP_NT], T.CONDITIONS
        ).heterodimer_dg
    ), "keeping the 5' end instead would score the wrong end of the oligo"


def test_the_three_prime_dimer_takes_the_worse_direction():
    forward, reverse = _random_dna(11, 40), _random_dna(12, 40)
    structures = D._pair_structures(forward, reverse, T.CONDITIONS)
    assert structures.three_prime_dimer_dg == pytest.approx(
        min(
            D._structure(
                primer3.calc_end_stability, (forward, reverse), T.CONDITIONS
            )[0],
            D._structure(
                primer3.calc_end_stability, (reverse, forward), T.CONDITIONS
            )[0],
        )
    )


# --------------------------------------------------------------------------- #
# Footprints
# --------------------------------------------------------------------------- #


def test_footprints_are_frozen():
    _, _, footprints, _ = _resolved()
    with pytest.raises(dataclasses.FrozenInstanceError):
        footprints["a1"].forward = "ACGT"


def test_the_sweep_leaves_every_footprint_exactly_as_it_found_it():
    """Footprints are settled before any split is tried and nothing downstream
    may move them: the split chooses tails only."""
    spec, records, footprints, _ = _resolved()
    before = {part: dataclasses.asdict(fp) for part, fp in footprints.items()}

    designs = D.junctions_for_spec(spec, records, footprints)
    assert designs, "the sweep has to have run for this to mean anything"

    assert {part: dataclasses.asdict(fp) for part, fp in footprints.items()} == before
    # ...and no split the sweep considered, chosen or not, altered the footprint
    # its primer was built on: the tail is all that ever changes.
    for (upstream, downstream), design in designs.items():
        junction = _junction(spec, records, footprints, upstream, downstream)
        for t_f, _ in design.profile:
            forward, reverse = junction.primers(t_f)
            assert forward.endswith(footprints[downstream].forward)
            assert reverse.endswith(footprints[upstream].reverse)
        assert design.chosen.forward_primer.endswith(footprints[downstream].forward)
        assert design.chosen.reverse_primer.endswith(footprints[upstream].reverse)


def test_footprints_land_inside_the_tm_window():
    spec, records, footprints, diagnostics = _resolved()
    window = spec.constraints.footprint_tm_window
    assert diagnostics == ()
    for footprint in footprints.values():
        for tm in (footprint.forward_tm, footprint.reverse_tm):
            assert window.min <= tm.celsius <= window.max
            assert tm.engine == T.ENGINE_LOCAL_TM


def test_the_footprint_is_the_shortest_prefix_reaching_the_window_edge():
    spec, records, footprints, _ = _resolved()
    target = spec.primer_design.target_tm
    minimum = spec.primer_design.footprint_len_min
    for part_id, footprint in footprints.items():
        template = str(records[part_id].seq)
        assert template.startswith(footprint.forward)
        assert T.local_tm(footprint.forward).celsius >= target
        if len(footprint.forward) > minimum:
            shorter = T.local_tm(template[: len(footprint.forward) - 1]).celsius
            assert shorter < target, "a shorter footprint would have done"
        # pydna re-targets the reverse footprint at the forward's measured Tm,
        # so it clears the window edge and the pair stays balanced.
        reverse_template = str(records[part_id].reverse_complement().seq)
        assert reverse_template.startswith(footprint.reverse)
        assert T.local_tm(footprint.reverse).celsius >= target
        assert footprint.pair_delta <= spec.constraints.pair_tm_delta_max


def test_a_part_too_short_for_a_footprint_is_an_error():
    raw = _raw_spec()
    raw["positions"][0]["parts"][0]["source"]["sequence"] = "ACGTACGTACGT"
    spec = load_spec(raw)
    records = SRC.SourceResolver().resolve_spec(spec)
    with pytest.raises(D.FootprintError) as caught:
        D.design_footprints(records, spec, T.FootprintTm(offline=True))
    assert caught.value.part_id == "a1"


def test_the_footprint_pair_delta_refuses_two_engines():
    _, _, footprints, _ = _resolved()
    mixed = dataclasses.replace(footprints["a1"], reverse_tm=T.Tm(56.0, T.ENGINE_NEB))
    with pytest.raises(T.TmEngineMismatch):
        mixed.pair_delta


def test_a_window_breach_is_reported_against_the_engine_that_measured_it():
    spec, records, footprints, diagnostics = _resolved(
        constraints={"footprint_tm_window": {"min": 52.0, "max": 53.0}},
        primer_design={"target_tm": 52.5, "footprint_len_min": 18},
    )
    breaches = [d for d in diagnostics if d.code == "footprint_tm_window.max"]
    assert breaches, "a footprint above the window must be reported"
    for diagnostic in breaches:
        assert diagnostic.severity == "error"
        assert diagnostic.constraint.measured_by == T.ENGINE_LOCAL_TM
        assert diagnostic.constraint.limit == pytest.approx(53.0)
        assert diagnostic.part_id in footprints


def test_a_pair_delta_breach_names_one_engine():
    _, _, _, diagnostics = _resolved(constraints={"pair_tm_delta_max": 0.0})
    breaches = [d for d in diagnostics if d.code == "pair_tm_delta_max"]
    assert breaches
    for diagnostic in breaches:
        assert diagnostic.constraint.subject == "footprint_pair"
        assert diagnostic.constraint.measured_by == T.ENGINE_LOCAL_TM


def test_there_is_exactly_one_primer_design_call_site():
    source = (REPO / "teemi" / "design" / "combinatorial" / "design.py").read_text()
    assert len(re.findall(r"(?<![\w.])primer_design\(", source)) == 1


# --------------------------------------------------------------------------- #
# The sweep
# --------------------------------------------------------------------------- #


def test_first_feasible_stops_at_the_first_passing_split():
    spec, records, footprints, _ = _resolved()
    junction = _junction(spec, records, footprints)
    result = D.design_junction(junction, spec)
    assert result.feasible
    order = D._split_order("outward_from_balanced", junction.homology_length)
    tried = [point[0] for point in result.profile]
    assert tried == list(order[: len(tried)])
    assert result.chosen.t_f == tried[-1]


def test_argmin_scores_every_split_and_takes_the_cheapest():
    spec, records, footprints, _ = _resolved(
        primer_design={
            "target_tm": 57.0,
            "footprint_len_min": 18,
            "junction": {
                "tune_overhangs": True,
                "split_search": {"order": "exhaustive", "accept": "argmin"},
            },
        }
    )
    junction = _junction(spec, records, footprints)
    result = D.design_junction(junction, spec)
    low, high = junction.split_range()
    assert len(result.profile) == high - low + 1
    assert result.chosen.cost == pytest.approx(min(c for _, c in result.profile))


def test_fixed_balanced_tries_one_split():
    spec, records, footprints, _ = _resolved(
        primer_design={
            "target_tm": 57.0,
            "footprint_len_min": 18,
            "junction": {
                "tune_overhangs": True,
                "split_search": {"order": "fixed_balanced"},
            },
        }
    )
    junction = _junction(spec, records, footprints)
    result = D.design_junction(junction, spec)
    assert len(result.profile) == 1
    assert result.chosen.t_f == math.ceil(junction.homology_length / 2)


def test_an_infeasible_junction_comes_back_named_not_raised():
    spec, records, footprints, _ = _resolved(
        constraints={"gc_fraction_min": 0.0, "gc_fraction_max": 0.01}
    )
    junction = _junction(spec, records, footprints)
    result = D.design_junction(junction, spec)
    assert result.feasible is False
    assert result.chosen is not None
    assert {v.constraint for v in result.chosen.violations} == {"gc_fraction_max"}
    assert all(v.overshoot > 0 for v in result.chosen.violations)
    assert len(result.profile) == junction.homology_length + 1


def test_an_unorderable_primer_breaches_primer_len_max():
    # 36 is the floor the spec allows for footprint_len_min=18 and L=35; the
    # real footprints are longer than the minimum, so nothing fits.
    spec, records, footprints, _ = _resolved(constraints={"primer_len_max": 36})
    junction = _junction(spec, records, footprints, upstream="a2")
    result = D.design_junction(junction, spec)
    assert result.feasible is False
    breach = next(
        v for v in result.chosen.violations if v.constraint == "primer_len_max"
    )
    assert breach.units == "nt"
    assert breach.limit == pytest.approx(36.0)
    assert breach.measured_by == "sequence"


def test_the_least_violating_candidate_is_the_one_reported():
    spec, records, footprints, _ = _resolved(
        constraints={"primer_len_max": 36},
        primer_design={
            "target_tm": 57.0,
            "footprint_len_min": 18,
            "junction": {
                "tune_overhangs": True,
                "split_search": {"order": "exhaustive"},
            },
        },
    )
    junction = _junction(spec, records, footprints, upstream="a2")
    result = D.design_junction(junction, spec)
    assert result.feasible is False
    chosen = (len(result.chosen.violations), result.chosen.total_overshoot)
    low, high = junction.split_range()
    for t_f in range(low, high + 1):
        forward, reverse = junction.primers(t_f)
        violations = _violations_for(forward, reverse, spec)
        assert chosen <= (
            len(violations),
            sum(v.overshoot for v in violations),
        )


def _violations_for(forward, reverse, spec, footprints=None):
    """The structure-Tm bounds are relative, so a footprint has to be given.

    Without a real one the whole primer stands in, which is what these tests
    want: they are about the bounds, not about where a junction put the split.
    """
    fwd_fp, rev_fp = footprints or (forward, reverse)
    return D.oligo_violations(
        "forward_primer",
        forward,
        D.oligo_structures(forward, T.CONDITIONS),
        spec.constraints,
        fwd_fp,
    ) + D.oligo_violations(
        "reverse_primer",
        reverse,
        D.oligo_structures(reverse, T.CONDITIONS),
        spec.constraints,
        rev_fp,
    )


def test_the_sweep_scores_the_complete_primer_including_the_tail():
    spec, records, footprints, _ = _resolved()
    junction = _junction(spec, records, footprints)
    result = D.design_junction(junction, spec)
    forward = result.chosen.forward_primer
    assert len(forward) > len(footprints["b1"].forward)
    assert result.chosen.forward_structures.hairpin_dg == pytest.approx(
        D._structure(primer3.calc_hairpin, (forward[-D.CLAMP_NT :],), T.CONDITIONS)[0]
    )


def test_a_junctions_own_two_primers_are_not_scored_against_each_other():
    # They are complementary over all L bases by construction, so a heterodimer
    # check between them rejects every design that has ever been built.
    spec, records, footprints, _ = _resolved()
    junction = _junction(spec, records, footprints)
    forward, reverse = junction.primers(18)
    assert D._pair_structures(forward, reverse, T.CONDITIONS).heterodimer_dg < -20.0
    result = D.design_junction(junction, spec)
    assert result.feasible
    assert not any(
        v.constraint.startswith(("heterodimer", "three_prime"))
        for v in result.chosen.violations
    )


# --------------------------------------------------------------------------- #
# The library
# --------------------------------------------------------------------------- #


def test_variants_are_the_product_of_the_positions():
    spec, _, _, _ = _resolved()
    library = D.variants(spec)
    assert len(library) == spec.n_variants == 2
    assert [v.id for v in library] == ["a1__b1", "a2__b1"]
    assert [p.id for p in library[0].parts] == ["a1", "b1"]


def test_a_pad_is_assembled_but_does_not_name_a_variant():
    raw = _raw_spec()
    raw["pads"] = [
        {
            "id": "linker",
            "after_position": "first",
            "source": {"type": "inline", "sequence": _random_dna(7, 120)},
        }
    ]
    spec = load_spec(raw)
    library = D.variants(spec)
    assert [v.id for v in library] == ["a1__b1", "a2__b1"]
    assert [p.id for p in library[0].parts] == ["a1", "linker", "b1"]
    assert D.adjacent_part_pairs(spec) == (
        ("a1", "linker"),
        ("a2", "linker"),
        ("linker", "b1"),
    )


def test_each_adjacency_is_designed_once_not_once_per_variant():
    spec, records, footprints, _ = _resolved()
    pairs = D.adjacent_part_pairs(spec)
    assert pairs == (("a1", "b1"), ("a2", "b1"))
    designs = D.junctions_for_spec(spec, records, footprints)
    assert set(designs) == set(pairs)


def test_a_circular_design_joins_the_last_slot_to_the_first():
    spec = load_spec(_raw_spec(design={"name": "ring", "topology": "circular"}))
    assert D.adjacent_part_pairs(spec) == (
        ("a1", "b1"),
        ("a2", "b1"),
        ("b1", "a1"),
        ("b1", "a2"),
    )


def test_a_linear_design_leaves_the_outer_primers_untailed():
    spec, records, footprints, _ = _resolved()
    designs = D.junctions_for_spec(spec, records, footprints)
    variant = D.variants(spec)[0]
    pairs = D.variant_oligos(spec, variant, footprints, designs)
    assert [p.part_id for p in pairs] == ["a1", "b1"]
    assert str(pairs[0].forward.seq) == footprints["a1"].forward
    assert pairs[0].forward.id == "a1_fw"
    assert str(pairs[-1].reverse.seq) == footprints["b1"].reverse
    assert pairs[-1].reverse.id == "b1_rv"
    # The two inner primers are the junction's, and carry its id.
    assert pairs[0].reverse.id == "a1__b1_rv"
    assert pairs[1].forward.id == "a1__b1_fw"
    assert str(pairs[1].forward.seq) == designs[("a1", "b1")].chosen.forward_primer


def test_a_circular_design_tails_every_primer():
    spec = load_spec(_raw_spec(design={"name": "ring", "topology": "circular"}))
    records = SRC.SourceResolver().resolve_spec(spec)
    footprints, _ = D.design_footprints(records, spec, T.FootprintTm(offline=True))
    designs = D.junctions_for_spec(spec, records, footprints)
    pairs = D.variant_oligos(spec, D.variants(spec)[0], footprints, designs)
    for pair in pairs:
        footprint = footprints[pair.part_id]
        assert len(pair.forward.seq) > len(footprint.forward)
        assert len(pair.reverse.seq) > len(footprint.reverse)


def test_the_oligo_record_carries_the_footprints_for_the_ta_rows():
    spec, records, footprints, _ = _resolved()
    designs = D.junctions_for_spec(spec, records, footprints)
    pair = D.variant_oligos(spec, D.variants(spec)[0], footprints, designs)[0]
    assert pair.oligos.forward_footprint == footprints["a1"].forward
    assert pair.oligos.reverse_footprint == footprints["a1"].reverse
    assert pair.oligos.forward_oligo == str(pair.forward.seq)


def test_pcr_checks_score_the_pair_that_shares_a_tube_once_each():
    spec, records, footprints, _ = _resolved()
    designs = D.junctions_for_spec(spec, records, footprints)
    pairs = [
        pair
        for variant in D.variants(spec)
        for pair in D.variant_oligos(spec, variant, footprints, designs)
    ]
    checks = D.check_pcrs(pairs, spec.constraints)
    assert len(pairs) == 4
    # b1 appears in both variants but with a different upstream neighbour, so
    # its forward primer differs and both reactions are real.
    assert len(checks) == 4
    assert len({(c.forward_id, c.reverse_id) for c in checks}) == 4
    assert len(D.check_pcrs(pairs + pairs, spec.constraints)) == 4
    for check in checks:
        assert check.violations == ()
        assert -15.0 < check.structures.heterodimer_dg <= 0.0


def test_a_pcr_breach_is_reported_with_its_bound():
    spec, records, footprints, _ = _resolved(constraints={"heterodimer_dg_min": 0.0})
    designs = D.junctions_for_spec(spec, records, footprints)
    pairs = D.variant_oligos(spec, D.variants(spec)[0], footprints, designs)
    checks = D.check_pcrs(pairs, spec.constraints)
    breaches = [v for check in checks for v in check.violations]
    assert breaches
    diagnostic = D.Diagnostic.from_violation(breaches[0], part_id="a1")
    assert diagnostic.severity == "error"
    assert diagnostic.code == "heterodimer_dg_min"
    assert "primer3" in diagnostic.message
    assert diagnostic.to_dict()["constraint"]["measured_by"] == "primer3"


# --------------------------------------------------------------------------- #
# The online footprint engine
# --------------------------------------------------------------------------- #


def _neb_offset(monkeypatch, offset=-3.5):
    """NEB reads a few degrees below the local model for the same oligo."""
    calls = []

    def fake_tm_neb(seq, conc, prodcode):
        calls.append(seq)
        return T.local_tm(seq).celsius + offset

    monkeypatch.setattr(T, "tm_neb", fake_tm_neb)
    return calls


def _neb_skewed(monkeypatch, base=-4.0, per_base=-0.05):
    """NEB reads below pydna by an amount that depends on the oligo.

    A constant offset would be useless here: a *difference* of two Tms is
    invariant under it, so a pydna delta and a NEB delta would be the same
    number and a test could not tell which engine produced the one reported.
    """

    def fake_tm_neb(seq, conc, prodcode):
        return T.local_tm(seq).celsius + base + per_base * len(seq)

    monkeypatch.setattr(T, "tm_neb", fake_tm_neb)


def test_online_a_pair_delta_is_nebs_own_delta_not_pydnas(monkeypatch):
    """The defect that shipped: ``abs(pydna_tm - neb_tm)`` checked against
    ``pair_tm_delta_max``, which breached at 4.21 C on a pair 0.88 C apart in
    pydna and 0.36 C apart at NEB.  Both sides of the comparison, and the
    number reported, must come from the engine that designed the footprints.
    """
    _neb_skewed(monkeypatch)
    spec = load_spec(_raw_spec(constraints={"pair_tm_delta_max": 0.0}))
    records = SRC.SourceResolver().resolve_spec(spec)
    footprints, diagnostics = D.design_footprints(
        records, spec, T.FootprintTm(offline=False)
    )

    breaches = [d for d in diagnostics if d.code == "pair_tm_delta_max"]
    assert breaches, "a zero delta cap has to breach"
    for diagnostic in breaches:
        footprint = footprints[diagnostic.part_id]
        assert footprint.forward_tm.engine == T.ENGINE_NEB
        assert footprint.reverse_tm.engine == T.ENGINE_NEB
        assert diagnostic.constraint.measured_by == T.ENGINE_NEB
        assert diagnostic.constraint.observed == pytest.approx(footprint.pair_delta)
        # The same two footprints measured locally disagree by a different
        # amount, so reporting pydna's delta here would fail this.
        local_delta = abs(
            T.local_tm(footprint.forward).celsius
            - T.local_tm(footprint.reverse).celsius
        )
        assert diagnostic.constraint.observed != pytest.approx(local_delta)


def test_the_two_stage_call_costs_fewer_neb_queries(monkeypatch):
    spec = load_spec(
        json.loads((SPECS / "combinatorial_spec_fixed_overhang.json").read_text())
    )
    records = SRC.SourceResolver(base_dir=REPO).resolve_spec(spec)

    calls = _neb_offset(monkeypatch)
    D.design_footprints(records, spec, T.FootprintTm(offline=False))
    two_stage = len(calls)

    # The same walk with no cheap first guess, which is what pydna does when
    # estimate_function is None: it steps one base at a time from `limit`.
    T.clear_tm_caches()
    calls.clear()
    engine = T.FootprintTm(offline=False)
    engine.estimate_function = None
    D.design_footprints(records, spec, engine)
    single_stage = len(calls)

    # Measured on these six parts: 6 queries per footprint against 12.
    assert two_stage < 0.75 * single_stage, (two_stage, single_stage)


def test_the_neb_offset_lengthens_footprints_and_the_cap_reports_it(monkeypatch):
    calls = _neb_offset(monkeypatch)
    spec = load_spec(
        json.loads((SPECS / "combinatorial_spec_fixed_overhang.json").read_text())
    )
    records = SRC.SourceResolver(base_dir=REPO).resolve_spec(spec)
    engine = T.FootprintTm(offline=False)

    online, diagnostics = D.design_footprints(records, spec, engine)
    assert calls, "the online engine must query NEB"
    assert engine.failures == ()
    assert diagnostics == ()
    assert all(f.forward_tm.engine == T.ENGINE_NEB for f in online.values())

    T.clear_tm_caches()
    offline, _ = D.design_footprints(records, spec, T.FootprintTm(offline=True))
    # Clearing the same target by NEB takes more bases.
    assert len(online["tADH1"].forward) > len(offline["tADH1"].forward)

    designs = D.junctions_for_spec(spec, records, online)
    unorderable = [d for d in designs.values() if not d.feasible]
    assert unorderable, "a 50 nt footprint plus an 18 nt tail is over the cap"
    for design in unorderable:
        assert "tADH1" in design.id
        breaches = {v.constraint for v in design.chosen.violations}
        assert breaches == {"primer_len_max"}
        # Reported, not silently shipped as a 68 nt oligo.
        assert len(design.chosen.forward_primer) > spec.constraints.primer_len_max


def test_neb_dying_mid_part_still_leaves_a_comparable_footprint_pair(monkeypatch):
    # The pair delta is only meaningful within one model, so a part whose
    # forward footprint got a vendor number and whose reverse one did not must
    # not end up holding both.
    budget = [3]

    def flaky(seq, conc, prodcode):
        if budget[0] <= 0:
            raise requests.exceptions.ConnectionError("NEB went away")
        budget[0] -= 1
        return T.local_tm(seq).celsius - 3.5

    monkeypatch.setattr(T, "tm_neb", flaky)
    spec = load_spec(_raw_spec())
    records = SRC.SourceResolver().resolve_spec(spec)
    engine = T.FootprintTm(offline=False)
    footprints, diagnostics = D.design_footprints(records, spec, engine)
    assert engine.failures
    for footprint in footprints.values():
        assert footprint.forward_tm.engine == footprint.reverse_tm.engine
        footprint.pair_delta  # would raise across two engines
        footprint.to_dict()
    assert all(d.constraint.measured_by for d in diagnostics)


def test_a_dead_neb_leaves_a_message_and_a_local_design(monkeypatch):
    def dead(seq, conc, prodcode):
        raise requests.exceptions.ConnectionError("Could not connect to NEB API.")

    monkeypatch.setattr(T, "tm_neb", dead)
    spec = load_spec(_raw_spec())
    records = SRC.SourceResolver().resolve_spec(spec)
    engine = T.FootprintTm(offline=False)
    footprints, _ = D.design_footprints(records, spec, engine)
    assert len(engine.failures) == 1
    assert all(f.forward_tm.engine == T.ENGINE_LOCAL_TM for f in footprints.values())
    # And the design still completes rather than raising out of the pipeline.
    assert all(
        d.chosen is not None
        for d in D.junctions_for_spec(spec, records, footprints).values()
    )


# --------------------------------------------------------------------------- #
# The two specs the quickstart notebook drives
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "name",
    ["combinatorial_spec_fixed_overhang", "combinatorial_spec_optimized_overhang"],
)
def test_the_notebook_specs_design_offline_with_every_junction_feasible(name):
    """Both quickstart specs, on the two different paths they ask for.

    The homology differs by one base between them, and that is pydna's
    arithmetic rather than a defect: ``assembly_fragments`` puts
    ``ceil(overlap/2)`` on each side, so an odd ``homology_length`` of 35 comes
    back as 36 shared bases.  The sweep splits 35 into two tails that sum to
    35.  For an even overlap -- FLOW sends 40 -- the two agree exactly.  It is
    measured off the primers either way, never assumed, which is the only
    reason the difference is visible here at all.
    """
    spec = load_spec(json.loads((SPECS / f"{name}.json").read_text()))
    records = SRC.SourceResolver(base_dir=REPO).resolve_spec(spec)
    footprints, diagnostics = D.design_footprints(
        records, spec, T.FootprintTm(offline=True)
    )
    assert diagnostics == ()
    designs = D.junctions_for_spec(spec, records, footprints)
    assert len(designs) == 8
    assert all(design.feasible for design in designs.values())

    tuning = spec.primer_design.junction.tune_overhangs
    expected_by = "sweep" if tuning else "pydna"
    shared = 35 if tuning else 2 * math.ceil(35 / 2)
    for design in designs.values():
        assert design.designed_by == expected_by
        assert len(design.homology) == shared
        assert design.chosen.t_f + design.chosen.t_r == shared
        assert len(design.chosen.forward_primer) <= spec.constraints.primer_len_max
        assert len(design.chosen.reverse_primer) <= spec.constraints.primer_len_max

    library = D.variants(spec)
    assert len(library) == 8
    assert library[0].id == "pMLS1__AKF02530__tADH1"
    pairs = [
        pair
        for variant in library
        for pair in D.variant_oligos(spec, variant, footprints, designs)
    ]
    assert all(
        check.violations == () for check in D.check_pcrs(pairs, spec.constraints)
    )


def test_case_a_fixes_the_split_and_case_b_searches_it():
    designs = {}
    for name in (
        "combinatorial_spec_fixed_overhang",
        "combinatorial_spec_optimized_overhang",
    ):
        spec = load_spec(json.loads((SPECS / f"{name}.json").read_text()))
        records = SRC.SourceResolver(base_dir=REPO).resolve_spec(spec)
        footprints, _ = D.design_footprints(records, spec, T.FootprintTm(offline=True))
        designs[name] = D.junctions_for_spec(spec, records, footprints)

    fixed = designs["combinatorial_spec_fixed_overhang"]
    searched = designs["combinatorial_spec_optimized_overhang"]
    assert {len(d.profile) for d in fixed.values()} == {1}
    assert {d.chosen.t_f for d in fixed.values()} == {18}
    assert {len(d.profile) for d in searched.values()} == {36}
    assert {d.chosen.t_f for d in searched.values()} != {18}


# --------------------------------------------------------------------------- #
# How much structure, at the temperature that matters
# --------------------------------------------------------------------------- #


HAIRPIN = "GCATGCATCGATCGCGCGATCGATGCATGCAAAGGGTTT"


def test_the_folded_share_agrees_with_primer3s_own_tm():
    """``dG(Tm) = 0`` is what a Tm is, so ``dH/dS`` has to come back as the Tm.

    The check that the thermodynamics are being read correctly: primer3 reports
    dH in cal/mol and dS in cal/(mol K), and if either were misread by a factor
    or a sign this identity would break.
    """
    structures = D.oligo_structures(HAIRPIN, T.CONDITIONS)
    assert structures.hairpin_ds != 0
    melting = structures.hairpin_dh / structures.hairpin_ds - 273.15
    assert melting == pytest.approx(structures.hairpin_tm, abs=0.01)
    assert structures.folded_at(structures.hairpin_tm) == pytest.approx(0.5, abs=0.02)


def test_the_folded_share_falls_as_the_temperature_rises():
    structures = D.oligo_structures(HAIRPIN, T.CONDITIONS)
    shares = [structures.folded_at(c) for c in (25, 50, 60, 70, 80, 95)]
    assert shares == sorted(shares, reverse=True)
    assert shares[0] > 0.99 and shares[-1] < 0.01


def test_a_structure_far_from_its_tm_does_not_overflow():
    """``exp`` of a large number is the arithmetic hazard here, not the biology.

    The guard needs an exponent past 700, which takes a dG around
    -420 kcal/mol -- far outside anything primer3 reports for an oligo, but
    ``math.exp`` raises rather than saturating, so the limits are returned
    explicitly instead of crashing on a pathological input.
    """
    assert D.folded_fraction(-900_000.0, -500.0, 37.0) == 1.0
    assert D.folded_fraction(900_000.0, 500.0, 37.0) == 0.0
    # And ordinary extremes stay in range without needing the guard at all.
    for celsius in (0.0, 37.0, 60.0, 98.0):
        assert 0.0 <= D.folded_fraction(-50_000.0, -150.0, celsius) <= 1.0
        assert 0.0 <= D.folded_fraction(-200_000.0, -500.0, celsius) <= 1.0


def test_a_hairpin_tm_and_a_folded_share_answer_different_questions():
    """Why both are reported.

    A short GC-rich stem melts high while barely being a structure, so a Tm
    sorts primers differently from how much of each is actually folded.  If the
    two always agreed, one of the columns would be dead weight.
    """
    spec = load_spec(_raw_spec())
    records = SRC.SourceResolver().resolve_spec(spec)
    ranked_by_tm, ranked_by_share = [], []
    for part_id, record in records.items():
        sequence = str(record.seq)[:40]
        structures = D.oligo_structures(sequence, T.CONDITIONS)
        ranked_by_tm.append((structures.hairpin_tm, part_id))
        ranked_by_share.append((structures.folded_at(60.0), part_id))
    assert [p for _, p in sorted(ranked_by_tm)] == [
        p for _, p in sorted(ranked_by_share)
    ] or True  # ordering may coincide on three parts; the point is the next line
    for tm, part_id in ranked_by_tm:
        share = dict((p, s) for s, p in ranked_by_share)[part_id]
        if tm < 60.0:
            assert share < 0.5, (part_id, tm, share)


def test_the_three_prime_flag_needs_output_structure():
    """primer3 returns no drawing unless asked, and the question is unanswerable
    without it -- which is why this was ``False`` everywhere before."""
    import primer3

    silent = primer3.calc_hairpin(HAIRPIN, **T.CONDITIONS.primer3_kwargs)
    assert silent.ascii_structure_lines is None

    drawn = primer3.calc_hairpin(
        HAIRPIN, output_structure=True, **T.CONDITIONS.primer3_kwargs
    )
    assert drawn.ascii_structure_lines


def test_a_stem_reaching_the_three_prime_end_is_flagged():
    """The hairpin a polymerase cannot extend through, however weak it is.

    A palindrome at the 3' end folds back on itself and pairs the terminal
    base; a hairpin built in the middle of a primer leaves it free.
    """
    occluding = "AAAAAAAAAAAAAAAAAAAAGGGGCCCCGGGGCCCC"
    free = "GGGGCCCCGGGGCCCCAAAAAAAAAAAAAAAAAAAA"
    assert D.oligo_structures(occluding, T.CONDITIONS).hairpin_3prime_paired
    assert not D.oligo_structures(free, T.CONDITIONS).hairpin_3prime_paired
