#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)

"""Tests for :mod:`teemi.design.combinatorial.sources`, the input boundary."""

import json
from pathlib import Path

import pytest
from pydna.dseqrecord import Dseqrecord

from teemi.design.combinatorial.sources import SourceError, SourceResolver
from teemi.design.combinatorial.spec import Part, load_spec

from tests.combinatorial_netguard import no_network  # noqa: F401

REPO = Path(__file__).resolve().parents[1]
FASTA = "tests/files_for_testing/promoter_test_data.fasta"


def file_part(part_id="p", **source) -> Part:
    return Part.model_validate(
        {
            "id": part_id,
            "source": {"type": "file", "format": "fasta", "path": FASTA, **source},
        }
    )


def inline_part(sequence, part_id="p", topology="linear") -> Part:
    return Part.model_validate(
        {
            "id": part_id,
            "source": {"type": "inline", "sequence": sequence, "topology": topology},
        }
    )


def test_resolves_the_notebook_spec():
    raw = json.loads((
        REPO / "tests" / "files_for_testing" / "combinatorial_spec_fixed_overhang.json"
    ).read_text())
    records = SourceResolver(base_dir=REPO).resolve_spec(load_spec(raw))
    assert len(records) == 6
    assert all(isinstance(r, Dseqrecord) for r in records.values())
    assert records["pMLS1"].id == "pMLS1"
    assert len(records["pMLS1"]) == 1000


def test_record_is_annotated_with_the_part_id():
    record = SourceResolver(base_dir=REPO).resolve(file_part(record_id="MLS1"))
    (feature,) = record.features
    assert feature.type == "misc_feature"
    assert feature.qualifiers["label"] == "p"
    assert int(feature.location.start) == 0
    assert int(feature.location.end) == len(record)


def test_inline_sequence_is_uppercased_and_topology_honoured():
    record = SourceResolver().resolve(inline_part("acgtacgt", topology="circular"))
    assert str(record.seq) == "ACGTACGT"
    assert record.circular is True


@pytest.mark.parametrize("sequence", ["ACGT-ACGT", "ACG TACGT", "ACGU"])
def test_what_is_not_dna_at_all_is_refused(sequence):
    with pytest.raises(SourceError, match="not DNA"):
        SourceResolver().resolve(inline_part(sequence))


@pytest.mark.parametrize("sequence", ["ACGTN", "NACGT", "ACGTRACGT"])
def test_a_degenerate_code_near_an_end_is_refused(sequence):
    """Where a primer is taken from, an ambiguous base has no Tm and no primer."""
    with pytest.raises(SourceError, match="degenerate code"):
        SourceResolver().resolve(inline_part(sequence))


def test_a_degenerate_code_deep_inside_is_carried_through():
    """Real library entries have them, and no primer ever reads that far in.

    AMA1 carries six, the nearest 1,252 bp from an end.  Refusing the part
    whole would turn away a construct that is perfectly designable: the bases a
    footprint and the homology are taken from are unambiguous, and the
    degenerate ones are only ever copied into the product.
    """
    interior = "ACGT" * 50 + "N" + "ACGT" * 50  # 201 bp, the N at 200
    record = SourceResolver().resolve(inline_part(interior))
    assert str(record.seq) == interior
    assert str(record.seq)[200] == "N"

    # One base closer to the end than a primer can reach is still refused.
    reaches = "ACGT" * 19 + "N" + "ACGT" * 50
    with pytest.raises(SourceError, match="degenerate code"):
        SourceResolver().resolve(inline_part(reaches))


def test_length_cap():
    with pytest.raises(SourceError, match="over the"):
        SourceResolver(max_part_bp=8).resolve(inline_part("ACGT" * 4))


def test_path_escaping_base_dir_is_refused(tmp_path):
    outside = tmp_path / "secret.fasta"
    outside.write_text(">x\nACGT\n")
    base = tmp_path / "base"
    base.mkdir()
    resolver = SourceResolver(base_dir=base)
    with pytest.raises(SourceError, match="outside base_dir"):
        resolver.resolve(file_part(path="../secret.fasta"))


def test_symlink_out_of_base_dir_is_refused(tmp_path):
    outside = tmp_path / "secret.fasta"
    outside.write_text(">x\nACGT\n")
    base = tmp_path / "base"
    base.mkdir()
    (base / "link.fasta").symlink_to(outside)
    with pytest.raises(SourceError, match="outside base_dir"):
        SourceResolver(base_dir=base).resolve(file_part(path="link.fasta"))


def test_absolute_path_is_refused():
    with pytest.raises(SourceError, match="absolute"):
        SourceResolver(base_dir=REPO).resolve(file_part(path=str(REPO / FASTA)))


def test_file_source_without_base_dir_is_refused():
    with pytest.raises(SourceError, match="base_dir"):
        SourceResolver().resolve(file_part(record_id="MLS1"))


def test_missing_file():
    with pytest.raises(SourceError, match="no such file"):
        SourceResolver(base_dir=REPO).resolve(file_part(path="tests/nope.fasta"))


def test_file_size_cap():
    with pytest.raises(SourceError, match="byte cap"):
        SourceResolver(base_dir=REPO, max_file_bytes=10).resolve(
            file_part(record_id="MLS1")
        )


def test_multi_record_file_needs_a_record_id():
    with pytest.raises(SourceError, match="set record_id"):
        SourceResolver(base_dir=REPO).resolve(file_part())


def test_unknown_record_id_lists_what_is_there():
    with pytest.raises(SourceError, match="MLS1"):
        SourceResolver(base_dir=REPO).resolve(file_part(record_id="absent"))


def test_single_record_file_needs_no_record_id(tmp_path):
    (tmp_path / "one.fasta").write_text(">only\nACGTACGT\n")
    record = SourceResolver(base_dir=tmp_path).resolve(file_part(path="one.fasta"))
    assert str(record.seq) == "ACGTACGT"


def test_genbank_features_are_kept(tmp_path):
    (tmp_path / "g.gb").write_text(
        Dseqrecord("ACGTACGTAA", name="gene").format("genbank")
    )
    part = Part.model_validate(
        {"id": "g", "source": {"type": "file", "format": "genbank", "path": "g.gb"}}
    )
    record = SourceResolver(base_dir=tmp_path).resolve(part)
    assert str(record.seq) == "ACGTACGTAA"
    # the file's own features survive alongside the part-id feature
    assert [f.type for f in record.features].count("misc_feature") == 1


def test_error_attributes_the_part():
    with pytest.raises(SourceError) as excinfo:
        SourceResolver().resolve(inline_part("ACGTN", part_id="pURA1"))
    assert excinfo.value.part_id == "pURA1"
    assert "pURA1" in str(excinfo.value)


def test_records_are_not_shared_between_parts():
    resolver = SourceResolver(base_dir=REPO)
    first = resolver.resolve(file_part("one", record_id="MLS1"))
    second = resolver.resolve(file_part("two", record_id="MLS1"))
    assert first is not second
    assert str(first.seq) == str(second.seq)
    assert (first.id, second.id) == ("one", "two")
    first.add_feature(0, 10, type_="misc_feature", label="scribble")
    assert len(second.features) == 1
