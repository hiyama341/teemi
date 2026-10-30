#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)

"""Tests for :mod:`teemi.design.combinatorial.spec`."""

import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from teemi.design.combinatorial.spec import DesignSpec, json_schema, load_spec

from tests.combinatorial_netguard import no_network  # noqa: F401

REPO = Path(__file__).resolve().parents[1]
SPEC_FILES = sorted(
    (REPO / "tests" / "files_for_testing").glob("combinatorial_spec_*overhang.json")
)


def minimal(**overrides) -> dict:
    raw = {
        "design": {"name": "two parts"},
        "positions": [
            {
                "id": "a",
                "parts": [
                    {"id": "a1", "source": {"type": "inline", "sequence": "ACGT"}}
                ],
            },
            {
                "id": "b",
                "parts": [
                    {"id": "b1", "source": {"type": "inline", "sequence": "TTTT"}}
                ],
            },
        ],
    }
    raw.update(overrides)
    return raw


@pytest.mark.parametrize("path", SPEC_FILES, ids=lambda p: p.stem)
def test_notebook_specs_load(path):
    spec = load_spec(json.loads(path.read_text()))
    assert spec.n_variants == 8
    assert spec.primer_design.junction.homology_length == 35
    assert spec.objective.currency == "DKK"


def test_defaults_are_the_documented_ones():
    spec = load_spec(minimal())
    assert spec.primer_design.target_tm == 57.0
    assert spec.primer_design.footprint_len_min == 18
    assert spec.primer_design.junction.split_search.order == "outward_from_balanced"
    assert spec.primer_design.junction.split_search.accept == "first_feasible"
    c = spec.constraints
    assert (c.hairpin_dg_min, c.homodimer_dg_min, c.heterodimer_dg_min) == (
        -6.0,
        -6.0,
        -6.0,
    )
    assert c.pair_tm_delta_max == 4.0
    assert c.primer_len_max == 60
    assert (c.footprint_tm_window.min, c.footprint_tm_window.max) == (52.0, 72.0)


def test_spec_is_hashable_and_frozen():
    spec = load_spec(minimal())
    assert hash(spec) == hash(load_spec(minimal()))
    with pytest.raises(ValidationError):
        spec.positions = ()


def test_unknown_key_is_rejected():
    with pytest.raises(ValidationError, match="hairpin_max_tm"):
        load_spec(minimal(constraints={"hairpin_max_tm": 45.0}))


def test_one_position_is_not_a_junction():
    raw = minimal()
    raw["positions"] = raw["positions"][:1]
    with pytest.raises(ValidationError):
        load_spec(raw)


def test_duplicate_part_ids_across_positions():
    raw = minimal()
    raw["positions"][1]["parts"][0]["id"] = "a1"
    with pytest.raises(ValidationError, match="unique"):
        load_spec(raw)


def test_pad_must_anchor_to_a_known_position():
    pad = {
        "id": "linker",
        "source": {"type": "inline", "sequence": "ACGT"},
        "after_position": "nope",
    }
    with pytest.raises(ValidationError, match="unknown position"):
        load_spec(minimal(pads=[pad]))


def test_pad_needs_exactly_one_anchor():
    pad = {"id": "linker", "source": {"type": "inline", "sequence": "ACGT"}}
    with pytest.raises(ValidationError, match="exactly one"):
        load_spec(minimal(pads=[pad]))
    pad = dict(pad, before_position="a", after_position="b")
    with pytest.raises(ValidationError, match="exactly one"):
        load_spec(minimal(pads=[pad]))


def test_ordered_slots_places_pads_at_their_anchor():
    pads = [
        {
            "id": "up",
            "source": {"type": "inline", "sequence": "ACGT"},
            "before_position": "a",
        },
        {
            "id": "mid",
            "source": {"type": "inline", "sequence": "ACGT"},
            "after_position": "a",
        },
    ]
    spec = load_spec(minimal(pads=pads))
    assert [slot for slot, _ in spec.ordered_slots()] == ["up", "a", "mid", "b"]
    assert [len(parts) for _, parts in spec.ordered_slots()] == [1, 1, 1, 1]
    assert spec.n_variants == 1


def test_target_tm_outside_the_footprint_window():
    with pytest.raises(ValidationError, match="outside"):
        load_spec(
            minimal(
                primer_design={"target_tm": 50.0},
                constraints={"footprint_tm_window": {"min": 52.0, "max": 72.0}},
            )
        )


def test_tm_window_must_be_ordered():
    with pytest.raises(ValidationError, match="below max"):
        load_spec(
            minimal(constraints={"footprint_tm_window": {"min": 72.0, "max": 52.0}})
        )


def test_gc_bounds_must_be_ordered():
    with pytest.raises(ValidationError, match="gc_fraction_min"):
        load_spec(minimal(constraints={"gc_fraction_min": 0.8, "gc_fraction_max": 0.5}))


def test_homology_the_primer_cap_cannot_accommodate():
    # 30 nt footprint floor + ceil(80/2) = 70 nt, over the 60 nt cap.
    with pytest.raises(ValidationError, match="primer_len_max"):
        load_spec(
            minimal(
                primer_design={
                    "footprint_len_min": 30,
                    "junction": {"homology_length": 80},
                }
            )
        )


def test_argmin_with_no_weights_is_refused():
    weights = {
        k: 0.0
        for k in ("hairpin_dg", "homodimer_dg", "heterodimer_dg", "length_imbalance")
    }
    with pytest.raises(ValidationError, match="argmin"):
        load_spec(
            minimal(
                primer_design={"junction": {"split_search": {"accept": "argmin"}}},
                objective={"weights": weights},
            )
        )


def test_json_schema_names_the_blocks():
    assert set(json_schema()["properties"]) == {
        "design",
        "positions",
        "pads",
        "primer_design",
        "constraints",
        "objective",
    }
    assert json_schema()["$defs"]["Constraints"]["additionalProperties"] is False


def test_load_spec_returns_the_root_model():
    assert isinstance(load_spec(minimal()), DesignSpec)


def test_an_id_may_not_contain_the_separator_that_joins_composite_ids():
    """The one check that makes every composite id in the package injective.

    ``Slug`` permits ``_``, so it permits ``__``, which is what primer, variant
    and junction ids are built by joining part ids with.  Parts ``P`` and
    ``Q__R`` side by side produce the same composite id as ``P__Q`` and ``R``:
    two different molecules under one ``primer_id``, two different variants
    under one ``variant.id``.  Both were silently collapsed by the dicts keyed
    on those ids, so a reaction was handed another molecule's numbers and a
    whole variant's oligos went missing from the order sheet.  Rejecting the
    separator here is cheaper and safer than making every consumer defensive.
    """
    for bad in ("A__B", "p__1", "x__", "a__b__c"):
        spec = minimal()
        spec["positions"][0]["parts"][0]["id"] = bad
        with pytest.raises(ValidationError, match="may not contain"):
            load_spec(spec)

    # A single underscore is still a perfectly good part id.
    ok = minimal()
    ok["positions"][0]["parts"][0]["id"] = "my_part.v2-1"
    assert load_spec(ok)


def test_a_position_id_may_not_contain_the_separator_either():
    """Variant ids are joined from part ids, junction ids from position-derived
    ones; the whole id space has to be free of the separator, not just parts."""
    spec = minimal()
    spec["positions"][0]["id"] = "pos__1"
    with pytest.raises(ValidationError, match="may not contain"):
        load_spec(spec)


def test_a_pad_id_may_not_equal_a_position_id():
    """Pads and positions are both slots, so they share that namespace.

    ``ordered_slots`` interleaves them and ``_variant_plans`` keys a dict on the
    result, so two slots with one name keep only the later -- and an
    ``after_position`` pad is the later one, so every variant reported the pad
    as the part chosen at that position and the overview table printed it.
    """
    spec = minimal()
    position_id = spec["positions"][0]["id"]
    spec["pads"] = [
        {
            "id": position_id,
            "source": {"type": "inline", "sequence": "ATGC" * 40},
            "after_position": position_id,
        }
    ]
    with pytest.raises(ValidationError, match="may not equal a position id"):
        load_spec(spec)

    # A distinct pad id is fine, and so is a part id matching a position id --
    # parts are never slot keys.
    spec["pads"][0]["id"] = "linker"
    assert load_spec(spec)
    shared = minimal()
    shared["positions"][0]["parts"][0]["id"] = shared["positions"][0]["id"]
    assert load_spec(shared)
