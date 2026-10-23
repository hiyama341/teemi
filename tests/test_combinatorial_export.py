#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)

"""Tests for :mod:`teemi.design.combinatorial.export`, the file-writing edge."""

import json
import zipfile
from pathlib import Path

import pytest
from Bio import SeqIO

from teemi.design.combinatorial.export import (
    write_assembly_figures,
    write_bundle,
    write_bundle_zip,
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
    order_sheet = {primer.label: primer for primer in linear.primers}
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
            # The name is still on the map, as prose beside the number.
            note = (feature.qualifiers.get("note") or [""])[0]
            name = order_sheet[feature.qualifiers["label"][0]].id
            assert name.replace("+", " + ") in note
    assert labelled, "the constructs have to carry primer_bind features at all"
    # And the shared molecule really does appear, under its shared name.
    assert any(
        primer.id.replace("+", " + ") in note and f"used {len(primer.uses)}x" in note
        for primer in shared
        for path in written
        for feature in SeqIO.read(path, "genbank").features
        if feature.type == "primer_bind"
        for note in feature.qualifiers.get("note") or []
    )


def test_a_primer_label_survives_the_genbank_line_wrap(tmp_path):
    """The label is the number because a long name does not survive the file.

    GenBank wraps a qualifier at 58 characters and Biopython rejoins the halves
    with a space, so labelling by name sent
    ``gene_gene-1+gene_gene-2+gene_gene-3+gene_gene-4__TER_ANID_TEF1_fw`` out
    and read ``...__TE R_ANID_TEF1_fw`` back -- a map naming a primer that is
    on no order sheet, the exact failure the labelling exists to prevent.

    No small fixture reaches 58 characters on its own, so this one lengthens
    the part ids until the junction names do.  Part ids are the user's to
    choose and nothing bounds their length, so this is ordinary input.

    Only the number is asserted here.  The note beside it is prose, and a name
    whose length comes from one long part id is a single token that GenBank
    cannot hold on a line however it is spaced -- which is the whole reason the
    label is a number.  That the note carries the name in the ordinary case is
    :func:`test_every_primer_label_in_a_genbank_is_on_the_order_sheet`.
    """
    spec = json.loads(SPEC_A.read_text())
    for position in spec["positions"]:
        for part in position["parts"]:
            part["id"] = f"{part['id']}-{'x' * 20}"
    result = plan(
        spec,
        sources=SourceResolver(base_dir=REPO),
        offline=True,
        include_designs=True,
    )
    long_named = [p for p in result.primers if len(p.id) > 58]
    assert long_named, "the lengthened ids have to produce a name that wraps"

    written = write_genbank_files(result.designs, tmp_path)
    by_label = {p.label: p for p in result.primers}
    found = set()
    for path in written:
        for feature in SeqIO.read(path, "genbank").features:
            if feature.type != "primer_bind":
                continue
            label = feature.qualifiers["label"][0]
            assert " " not in label and label.isdigit(), label
            assert label in by_label, (label, path.name)
            found.add(label)
    # And a name too long for a label is on a map, found by its number.
    assert found & {p.label for p in long_named}


# --------------------------------------------------------------------------- #
# The bundle
# --------------------------------------------------------------------------- #


PROMISED = (
    "LAYOUT.txt",
    "SUMMARY.txt",
    "design.json",
    "primers/primers.csv",
    "primers/analysis.csv",
    "pcrs/pcrs.csv",
    "parts/parts.csv",
    "constructs/assemblies.txt",
)


def test_the_bundle_writes_everything_its_own_layout_promises(linear, tmp_path):
    """The map in the folder has to describe the folder.

    ``LAYOUT.txt`` is what a person reads to find their way around, so a file
    it names and the bundle does not write is worse than no map at all.
    """
    folder = write_bundle(linear, tmp_path / "build")
    for relative in PROMISED:
        assert (folder / relative).is_file(), relative

    layout = (folder / "LAYOUT.txt").read_text()
    for relative in PROMISED:
        if relative != "LAYOUT.txt":
            named = relative.split("/")[-1]
            assert named in layout or relative.replace("/", "/*") in layout, relative

    written = {
        str(path.relative_to(folder)) for path in folder.rglob("*") if path.is_file()
    }
    maps = {name for name in written if name.endswith(".gb")}
    assert len(maps) == len(list(linear.designs))
    assert written == set(PROMISED) | maps


def test_the_zip_holds_the_same_files_under_one_root(linear, tmp_path):
    """What a backend hands over: one file, the folder intact inside it."""
    folder = write_bundle(linear, tmp_path / "loose")
    archive = write_bundle_zip(linear, tmp_path / "build.zip")
    assert archive.is_file() and archive.suffix == ".zip"

    with zipfile.ZipFile(archive) as opened:
        assert opened.testzip() is None, "a corrupt member"
        names = opened.namelist()
    roots = {name.split("/")[0] for name in names}
    assert roots == {"build"}, roots
    assert {name.split("/", 1)[1] for name in names} == {
        str(path.relative_to(folder)) for path in folder.rglob("*") if path.is_file()
    }


def test_every_construct_map_annotates_its_primers_over_the_whole_primer(
    circular, tmp_path
):
    """Including the tails, which is what the maps are for.

    ``pydna.amplify.pcr`` annotates only where a primer anneals to its
    template, so the feature stopped at the footprint and under-reported a
    49 nt primer as 35 bp.  The circular design is the harder case: its
    origin-crossing primer needs a two-part location rather than a clamp.
    """
    folder = write_bundle(circular, tmp_path / "build")
    by_id = {p.label: p for p in circular.primers}
    checked = 0
    for path in sorted((folder / "constructs").glob("*.gb")):
        record = SeqIO.read(path, "genbank")
        for feature in record.features:
            if feature.type != "primer_bind":
                continue
            label = (feature.qualifiers.get("label") or [""])[0]
            primer = by_id.get(label)
            assert primer is not None, (label, path.name)
            span = sum(len(part) for part in feature.location.parts)
            assert span == len(primer.sequence), (label, span, len(primer.sequence))
            checked += 1
    assert checked, "the maps have to carry primer_bind features at all"
    # And the parts are annotated too, which is the other half of a map.
    record = SeqIO.read(next((folder / "constructs").glob("*.gb")), "genbank")
    parts = [f for f in record.features if f.type == "misc_feature"]
    assert parts
    for feature in parts:
        label = feature.qualifiers["label"]
        name = label[0] if isinstance(label, list) else label
        assert len(name) > 1, "a part label must not be a single character"


def test_the_parts_table_counts_reactions_and_constructs_separately(linear, tmp_path):
    """A part at a combinatorial position is amplified once per neighbour pair,
    not once per construct, and sourcing cares about the first number while
    risk cares about the second."""
    import pandas as pd

    folder = write_bundle(linear, tmp_path / "build")
    parts = pd.read_csv(folder / "parts" / "parts.csv")
    assert set(parts["part_id"]) == {f.part_id for f in linear.footprints}
    assert (parts["reactions"] >= 1).all()
    assert (
        parts["in_constructs"] >= parts["reactions"]
    ).any(), "this library has a part reused across constructs by one reaction"
    assert "position" in parts.columns


def test_part_provenance_comes_from_the_caller(linear, tmp_path):
    """Where a part can be obtained is not something the design knows.

    The spec carries sequences and ids; an Addgene number lives in whatever
    library the caller resolved the part from, so it is passed in and any keys
    become columns.
    """
    import pandas as pd

    first = sorted(f.part_id for f in linear.footprints)[0]
    folder = write_bundle(
        linear,
        tmp_path / "build",
        part_info={first: {"addgene": "#12345", "doi": "10.1000/x"}},
    )
    parts = pd.read_csv(folder / "parts" / "parts.csv")
    assert "addgene" in parts.columns and "doi" in parts.columns
    row = parts[parts["part_id"] == first].iloc[0]
    assert row["addgene"] == "#12345"

    # And the summary says which parts it was given nothing for, rather than
    # leaving a blank column to be noticed.
    summary = (folder / "SUMMARY.txt").read_text()
    assert "no provenance given for" in summary


def test_the_summary_says_what_to_order_run_source_and_check(linear, tmp_path):
    folder = write_bundle(linear, tmp_path / "build")
    summary = (folder / "SUMMARY.txt").read_text()
    for heading in (
        "WHAT TO ORDER",
        "WHAT TO RUN",
        "WHAT TO SOURCE",
        "WHAT WAS CHECKED",
    ):
        assert heading in summary, heading
    assert linear.spec.design.name in summary
    assert f"{len(linear.primers)} distinct oligos" in summary
    assert f"{len(linear.pcrs)} reactions" in summary
    assert f"{len(linear.variants)} constructs simulated" in summary
    # A flag is explained as a forecast, so nobody reads it as a refusal.
    flagged = [p for p in linear.to_dict()["pcrs"] if p["flagged_by"]]
    if flagged:
        assert "prediction, not a refusal" in summary


def test_assemblies_are_one_file_with_one_figure_per_construct(linear, tmp_path):
    folder = write_bundle(linear, tmp_path / "build")
    text = (folder / "constructs" / "assemblies.txt").read_text()
    for design in linear.designs:
        assert design.id in text, design.id
    assert text.count("-" * 72) == max(len(list(linear.designs)) - 1, 0)
    # The per-figure files are folded into the one file, not left beside it.
    assert not list((folder / "constructs").glob("*.txt")) == []
    assert sorted(p.name for p in (folder / "constructs").glob("*.txt")) == [
        "assemblies.txt"
    ]
