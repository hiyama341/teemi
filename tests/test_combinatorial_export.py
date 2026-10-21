#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)

"""Tests for :mod:`teemi.design.combinatorial.export`, the file-writing edge."""

import json
from pathlib import Path

import pytest
from Bio import SeqIO

from teemi.design.combinatorial.export import (
    write_assembly_figures,
    write_genbank_files,
)
from teemi.design.combinatorial.report import VariantPlan, plan
from teemi.design.combinatorial.simulate import Verification
from teemi.design.combinatorial.sources import SourceResolver

from tests.combinatorial_netguard import no_network  # noqa: F401

REPO = Path(__file__).resolve().parents[1]
SPEC_A = REPO / "tests" / "files_for_testing" / "combinatorial_spec_fixed_overhang.json"
CIRCULAR = REPO / "tests" / "files_for_testing" / "combinatorial_spec_example.json"


@pytest.fixture(scope="module")
def linear():
    return plan(
        json.loads(SPEC_A.read_text()),
        sources=SourceResolver(base_dir=REPO),
        offline=True,
        include_designs=True,
    )


@pytest.fixture(scope="module")
def circular():
    return plan(
        json.loads(CIRCULAR.read_text()),
        sources=SourceResolver(base_dir=REPO),
        offline=True,
        include_designs=True,
    )


def test_genbank_files_are_named_after_their_variant(linear, tmp_path):
    designs = list(linear.designs)[:2]
    written = write_genbank_files(designs, tmp_path)
    assert [p.name for p in written] == [f"{d.id}.gb" for d in designs]
    assert all(p.parent == tmp_path.resolve() for p in written)


def test_assembly_figures_are_named_after_their_variant(linear, tmp_path):
    designs = list(linear.designs)[:1]
    written = write_assembly_figures(designs, tmp_path)
    assert [p.name for p in written] == [f"{designs[0].id}.txt"]
    text = written[0].read_text()
    assert text.startswith(f"{designs[0].id}  {designs[0].length_bp} bp  linear")
    for part_id in designs[0].fragment_ids:
        assert part_id in text


def test_the_genbank_file_carries_the_part_and_primer_annotations(linear, tmp_path):
    """pydna propagates the fragments' own features into the contig, so no
    feature is constructed on the assembled construct here."""
    design = list(linear.designs)[0]
    path = write_genbank_files([design], tmp_path)[0]
    record = SeqIO.read(path, "genbank")
    assert len(record.seq) == design.length_bp
    kinds = [f.type for f in record.features]
    assert kinds.count("misc_feature") == len(design.fragment_ids)
    assert kinds.count("primer_bind") == 2 * len(design.fragment_ids)
    labels = {f.qualifiers["label"][0] for f in record.features}
    assert set(design.fragment_ids) <= labels


def test_a_circular_construct_is_written_as_circular(circular, tmp_path):
    design = list(circular.designs)[0]
    path = write_genbank_files([design], tmp_path)[0]
    record = SeqIO.read(path, "genbank")
    assert record.annotations["topology"] == "circular"
    figure = write_assembly_figures([design], tmp_path)[0].read_text()
    assert "circular" in figure.splitlines()[0]


def test_a_variant_with_no_construct_writes_nothing(linear, tmp_path):
    """``plan`` without ``include_designs`` keeps no records, so there is
    nothing to write and the returned paths say so."""
    bare = plan(
        json.loads(SPEC_A.read_text()),
        sources=SourceResolver(base_dir=REPO),
        offline=True,
    )
    assert bare.designs == ()
    assert write_genbank_files(bare.variants, tmp_path) == ()
    assert write_assembly_figures(bare.variants, tmp_path) == ()
    assert list(tmp_path.iterdir()) == []


def test_the_output_directory_is_created(linear, tmp_path):
    target = tmp_path / "deep" / "nested"
    written = write_genbank_files(list(linear.designs)[:1], target)
    assert written[0].is_file()


def test_a_traversing_variant_id_is_refused(linear, tmp_path):
    design = list(linear.designs)[0]
    escaping = VariantPlan(
        id="../escaped",
        choices=design.choices,
        fragment_ids=design.fragment_ids,
        length_bp=design.length_bp,
        verification=design.verification,
        construct=design.construct,
    )
    with pytest.raises(ValueError, match="outside"):
        write_genbank_files([escaping], tmp_path / "out")
    with pytest.raises(ValueError, match="outside"):
        write_assembly_figures([escaping], tmp_path / "out")
    assert not (tmp_path / "escaped.gb").exists()


def test_writing_the_whole_library_gives_one_file_per_variant(circular, tmp_path):
    genbank = write_genbank_files(circular.designs, tmp_path)
    figures = write_assembly_figures(circular.designs, tmp_path)
    assert len(genbank) == len(figures) == len(circular.variants) == 12
    assert len({p.name for p in genbank}) == 12


def test_a_hand_built_plan_with_no_construct_is_skipped_silently(tmp_path):
    plan_without_record = VariantPlan(
        id="v1",
        choices=(("a", "a1"),),
        fragment_ids=("a1",),
        length_bp=0,
        verification=Verification(
            verdict="no_product",
            message="nothing assembled",
            n_products=0,
            expected_bp=0,
            observed_bp=None,
        ),
    )
    assert write_genbank_files([plan_without_record], tmp_path) == ()
    assert write_assembly_figures([plan_without_record], tmp_path) == ()


def test_every_primer_label_in_a_genbank_is_on_the_order_sheet(linear, tmp_path):
    """A label in the map has to be findable in the tube rack.

    ``pydna.amplify.pcr`` names a ``primer_bind`` feature after the junction
    that produced it, which for a molecule several junctions share is only one
    of its names.  Left alone, a GenBank would point at
    ``AKF02530__tADH1_rv`` while the order sheet listed
    ``AAA17732+AKF02530__tADH1_rv`` and nothing connected the two.
    """
    order_sheet = {primer.id for primer in linear.primers}
    shared = [primer for primer in linear.primers if len(primer.uses) > 1]
    assert shared, "this library has to contain a shared oligo to be a real test"

    written = write_genbank_files(linear.designs, tmp_path)
    assert written
    labelled = 0
    for path in written:
        for feature in SeqIO.read(path, "genbank").features:
            if feature.type != "primer_bind":
                continue
            for label in feature.qualifiers.get("label") or []:
                assert label in order_sheet, (label, path.name)
                labelled += 1
    assert labelled, "the constructs have to carry primer_bind features at all"
    # And the shared molecule really does appear under its shared name.
    assert any(
        label in {p.id for p in shared}
        for path in written
        for feature in SeqIO.read(path, "genbank").features
        if feature.type == "primer_bind"
        for label in feature.qualifiers.get("label") or []
    )
