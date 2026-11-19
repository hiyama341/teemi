#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)

"""Tests for :mod:`teemi.design.combinatorial.idt`, the vendor order shapes."""

import json
from pathlib import Path

import pytest

from teemi.design.combinatorial.idt import (
    DEFAULT_PURIFICATION,
    DEFAULT_SCALE,
    IDT_NAME_MAX,
    IDT_PLATE_COLUMNS,
    IDT_TUBE_COLUMNS,
    oligo_name,
    plate_count,
    plate_frames,
    primers_fasta,
    tube_frame,
    well_positions,
)
from teemi.design.combinatorial.report import plan
from teemi.design.combinatorial.sources import SourceResolver

from tests.combinatorial_netguard import no_network  # noqa: F401

REPO = Path(__file__).resolve().parents[1]
SPEC_A = REPO / "tests" / "files_for_testing" / "combinatorial_spec_fixed_overhang.json"


@pytest.fixture(scope="module")
def linear():
    return plan(
        json.loads(SPEC_A.read_text()),
        sources=SourceResolver(base_dir=REPO),
        offline=True,
        include_designs=True,
    )


class _Oligo:
    """The two attributes :func:`oligo_name` reads, and nothing else."""

    def __init__(self, label, id):
        self.label = label
        self.id = id


# -- the plate geometry ------------------------------------------------------


def test_wells_are_filled_down_the_columns():
    """An 8-channel pipette spans one column, so a column is the unit."""
    wells = well_positions()
    assert wells[:8] == ("A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1")
    # The ninth oligo starts the second column, not the second row.
    assert wells[8] == "A2"
    assert len(wells) == 96
    assert wells[-1] == "H12"


def test_every_well_of_a_plate_is_named_once():
    assert len(set(well_positions())) == 96


def test_an_impossible_plate_is_refused():
    with pytest.raises(ValueError, match="rows"):
        well_positions(rows=0)
    with pytest.raises(ValueError, match="rows"):
        well_positions(rows=27)
    with pytest.raises(ValueError, match="cols"):
        well_positions(cols=0)


def test_plate_count_counts_whole_plates_and_the_remainder():
    assert plate_count(0) == 0
    assert plate_count(1) == 1
    assert plate_count(96) == 1
    assert plate_count(97) == 2


# -- the name ---------------------------------------------------------------


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("pMLS1_fw", "001_pMLS1_fw"),
        ("sgRNA #1", "001_sgRNA_1"),
        ("a (b) c", "001_a_b_c"),
        # ``__`` is teemi's own separator and survives; see _UNSAFE_NAME_CHARS.
        ("sgRNA_#1", "001_sgRNA__1"),
        ("AAA17732+AKF02530__tADH1_rv", "001_AAA17732+AKF02530__tADH1_rv"),
        ("", "001"),
    ],
)
def test_a_name_is_folded_to_what_a_vendor_form_takes(raw, expected):
    assert oligo_name(_Oligo("001", raw)) == expected


def test_a_name_is_capped():
    name = oligo_name(_Oligo("007", "X" * 200))
    assert len(name) <= IDT_NAME_MAX
    assert name.startswith("007_")


def test_two_long_ids_still_get_two_names():
    """The cap cannot collide, because the label leads and labels are unique.

    A shared oligo's id is the join of every junction that asked for it, so ids
    run long and two of them can agree for the first fifty characters.
    """
    long_a = oligo_name(_Oligo("001", "Y" * 200))
    long_b = oligo_name(_Oligo("002", "Y" * 200))
    assert long_a != long_b


def test_a_nameless_oligo_still_gets_a_name():
    assert oligo_name(_Oligo("", "")) == "unnamed_oligo"


# -- the two shapes agree ---------------------------------------------------


def test_the_tube_sheet_has_idts_headers_in_order(linear):
    assert tuple(tube_frame(linear).columns) == IDT_TUBE_COLUMNS


def test_the_plate_has_idts_headers_in_order(linear):
    for frame in plate_frames(linear):
        assert tuple(frame.columns) == IDT_PLATE_COLUMNS


def test_the_plate_and_the_tube_sheet_are_one_order(linear):
    """The same molecules in the same order, so the two files cannot disagree.

    Compared on SEQUENCES rather than on counts: an oligo is identified by its
    sequence and not by its name, so the right number of the wrong oligos is
    exactly the failure a count would pass.
    """
    tubes = tube_frame(linear)
    plated = [row for frame in plate_frames(linear) for row in frame["Sequence"]]
    assert plated == tubes["Sequence"].tolist()
    names = [row for frame in plate_frames(linear) for row in frame["Name"]]
    assert names == tubes["Name"].tolist()


def test_one_row_per_distinct_molecule(linear):
    assert len(tube_frame(linear)) == len(linear.primers)


def test_a_shared_oligo_is_ordered_once(linear):
    """A molecule two junctions need is one tube and one well.

    The fixture has shared oligos -- that is what makes this worth asserting --
    and a sheet that listed them per junction would buy the same DNA twice.
    """
    shared = [p for p in linear.primers if len(p.uses) > 1]
    assert shared, "fixture no longer has a shared oligo; this test is void"
    sequences = tube_frame(linear)["Sequence"].tolist()
    assert len(sequences) == len(set(sequences))


def test_the_first_fragments_pair_lands_stacked(linear):
    """001 in A1 and 002 in B1, one directly below the other.

    Not a guarantee for every pair -- a shared oligo shifts everything after it
    up a well -- but it is what the fill order is for, so the first pair is
    asserted where no dedup can have happened yet.
    """
    first = plate_frames(linear)[0]
    assert first["Well Position"].tolist()[:2] == ["A1", "B1"]
    assert first["Name"].tolist()[0].startswith("001_")
    assert first["Name"].tolist()[1].startswith("002_")


def test_the_scale_and_purification_are_written_on_every_row(linear):
    frame = tube_frame(linear)
    assert set(frame["Scale"]) == {DEFAULT_SCALE}
    assert set(frame["Purification"]) == {DEFAULT_PURIFICATION}


# -- overflow and emptiness -------------------------------------------------


class _Plan:
    def __init__(self, primers):
        self.primers = primers


def _fake(n):
    class P:
        def __init__(self, i):
            self.label = f"{i:03d}"
            self.id = f"oligo{i}"
            self.sequence = "ATGC"
            self.uses = ("x",)

    return _Plan([P(i) for i in range(1, n + 1)])


def test_an_order_past_one_plate_starts_a_second_at_a1():
    plates = plate_frames(_fake(100))
    assert len(plates) == 2
    assert plates[0]["Well Position"].iloc[-1] == "H12"
    assert plates[1]["Well Position"].tolist() == ["A1", "B1", "C1", "D1"]


def test_no_oligo_is_lost_across_a_plate_boundary():
    """Chunking is where a column could silently go missing."""
    for count in (1, 95, 96, 97, 192, 193):
        plates = plate_frames(_fake(count))
        assert sum(len(p) for p in plates) == count
        assert len(plates) == plate_count(count)


def test_an_empty_order_gives_no_plates():
    """Not one empty plate: IDT would refuse it and openpyxl cannot write it."""
    assert plate_frames(_fake(0)) == []
    assert tube_frame(_fake(0)).empty


# -- FASTA ------------------------------------------------------------------


def test_the_fasta_ids_are_the_names_on_the_order(linear):
    """One owner for the oligo name, so a record and a tube are one molecule."""
    text = primers_fasta(linear)
    ids = [line[1:] for line in text.splitlines() if line.startswith(">")]
    assert ids == tube_frame(linear)["Name"].tolist()


def test_the_fasta_sequences_are_unwrapped(linear):
    """An oligo is under 100 bases; a newline in one is a newline in the middle
    of the only thing anybody copies out of this file."""
    lines = primers_fasta(linear).splitlines()
    sequences = [line for line in lines if not line.startswith(">")]
    assert sequences == tube_frame(linear)["Sequence"].tolist()


def test_an_empty_order_gives_empty_fasta():
    assert primers_fasta(_fake(0)) == ""
