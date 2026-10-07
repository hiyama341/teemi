#!/usr/bin/env python

# Test benchlingapi module
#
# All Benchling access is mocked: the module-level ``session`` is replaced by a
# MagicMock, and the module is imported with dotenv/Session patched so that no
# .env file is read and no real Benchling session is created.

import datetime as _real_datetime
import importlib
import types
from unittest import mock

import pandas as pd
import pytest

# from_benchling relies on these sub-modules being importable as attributes of
# ``Bio`` / ``pydna`` (the module itself only does ``import Bio`` / ``import pydna``).
import Bio.SeqRecord  # noqa: F401
import pydna.primer

# Importing the module we are testing (without touching .env or the network)
with mock.patch("dotenv.load_dotenv"), mock.patch(
    "dotenv.find_dotenv", return_value=""
), mock.patch("benchlingapi.Session"):
    import teemi.legacy.lims.benchling_api as benchling_api


INVENTORY_COLUMNS = [
    "batchEntId",
    "parentBoxPlateName",
    "parentBoxPlatePos",
    "volume",
    "Concentration (ng/ul)",
]


def _write_inventory_csv(path):
    pd.DataFrame(
        [
            ["seq_abc", "BoxA", "A1", 50.0, 120.0],
            ["seq_other", "BoxA", "A2", 10.0, 5.0],
            ["seq_abc", "BoxB", "H12", 20.0, 80.0],
        ],
        columns=INVENTORY_COLUMNS,
    ).to_csv(path, index=False)
    return path


def _fake_benchling_dump(features=None):
    return {
        "id": "seq_abc",
        "name": "pTEST_benchling_name",
        "bases": "ATGCATGCATGCATGCATGC",
        "annotations": features if features is not None else [],
        "fields": {"Resistance": {"value": "AmpR"}},
        "customFields": {"comment": "made for testing", "Owner": "LL"},
        "isCircular": False,
        "folderId": "lib_123",
        "aliases": [],
    }


@pytest.fixture
def fake_session(monkeypatch):
    session = mock.MagicMock(name="benchling_session")
    monkeypatch.setattr(benchling_api, "session", session)
    return session


@pytest.fixture
def frozen_today(monkeypatch):
    class _FrozenDatetime(_real_datetime.datetime):
        @classmethod
        def today(cls):
            return cls(2024, 3, 5, 12, 0, 0)

    monkeypatch.setattr(
        benchling_api, "datetime", types.SimpleNamespace(datetime=_FrozenDatetime)
    )


@pytest.fixture
def inventory_csv(tmp_path, monkeypatch):
    """Make from_benchling read batches from a temporary inventory CSV."""
    csv_path = _write_inventory_csv(tmp_path / "inventory.csv")
    real_update = benchling_api.update_loc_vol_conc
    monkeypatch.setattr(
        benchling_api,
        "update_loc_vol_conc",
        lambda seqRecord: real_update(seqRecord, DBpath=str(csv_path)),
    )
    return csv_path


def test_module_initialises_session_from_environment(monkeypatch):
    monkeypatch.setenv("API_KEY", "test-api-key")
    monkeypatch.setenv("HOME_url", "https://example.test/api/v2")
    fake_session_class = mock.MagicMock(name="Session")

    with mock.patch(
        "dotenv.find_dotenv", return_value="/nonexistent/.env"
    ) as find_dotenv, mock.patch("dotenv.load_dotenv") as load_dotenv, mock.patch(
        "benchlingapi.Session", fake_session_class
    ):
        module = importlib.reload(benchling_api)

    find_dotenv.assert_called_once_with()
    load_dotenv.assert_called_once_with("/nonexistent/.env")
    fake_session_class.assert_called_once_with(
        api_key="test-api-key", home="https://example.test/api/v2"
    )
    assert module.api_key == "test-api-key"
    assert module.home_url == "https://example.test/api/v2"
    assert module.session is fake_session_class.return_value


def test_sequence_to_benchling(fake_session):
    fake_session.Folder.find_by_name.return_value.id = "lib_123"

    result = benchling_api.sequence_to_benchling(
        "Primers folder", "P001_fw", "ATGCATGCAT", "Primer"
    )

    assert result is None
    fake_session.Folder.find_by_name.assert_called_once_with("Primers folder")
    fake_session.DNASequence.assert_called_once_with(
        name="P001_fw", bases="ATGCATGCAT", folder_id="lib_123", is_circular=False
    )
    dna = fake_session.DNASequence.return_value
    assert dna.method_calls == [
        mock.call.save(),
        mock.call.set_schema("Primer"),
        mock.call.register(),
    ]


def test_sequence_to_benchling_unknown_schema_is_not_set(fake_session):
    fake_session.Folder.find_by_name.return_value.id = "lib_456"

    benchling_api.sequence_to_benchling("Misc", "seq1", "GGGCCC", "Not a schema")

    fake_session.DNASequence.assert_called_once_with(
        name="seq1", bases="GGGCCC", folder_id="lib_456", is_circular=False
    )
    dna = fake_session.DNASequence.return_value
    assert dna.method_calls == [mock.call.save(), mock.call.register()]


def test_update_loc_vol_conc(tmp_path):
    csv_path = _write_inventory_csv(tmp_path / "inventory.csv")
    record = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATGC"), id="seq_abc")
    record.annotations["batches"] = [{"location": "stale"}]

    result = benchling_api.update_loc_vol_conc(record, DBpath=str(csv_path))

    assert result is record
    assert record.annotations["batches"] == [
        {"box": "BoxA", "position": "A1", "volume": 50, "concentration": 120, "location": "BoxA_A1"},
        {"box": "BoxB", "position": "H12", "volume": 20, "concentration": 80, "location": "BoxB_H12"},
    ]
    assert all(type(b["volume"]) is int for b in record.annotations["batches"])


def test_update_loc_vol_conc_no_matching_batches(tmp_path):
    csv_path = _write_inventory_csv(tmp_path / "inventory.csv")
    record = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATGC"), id="seq_not_in_inventory")

    benchling_api.update_loc_vol_conc(record, DBpath=str(csv_path))

    assert record.annotations["batches"] == []


def test_from_benchling(fake_session, frozen_today, inventory_csv):
    fake_session.DNASequence.find_by_name.return_value.dump.return_value = (
        _fake_benchling_dump()
    )

    record = benchling_api.from_benchling("pTEST")

    fake_session.DNASequence.find_by_name.assert_called_once_with("pTEST")
    assert type(record) is Bio.SeqRecord.SeqRecord
    assert str(record.seq) == "ATGCATGCATGCATGCATGC"
    assert record.id == "seq_abc"
    assert record.name == "pTEST"
    assert record.features == []

    annotations = record.annotations
    # Benchling "fields" + "customFields" end up in the annotations ...
    assert annotations["Resistance"] == {"value": "AmpR"}
    assert annotations["Owner"] == "LL"
    # ... with "comment" renamed to "commentary" (Bio.SeqIO cannot write "comment")
    assert "comment" not in annotations
    assert annotations["commentary"] == "made for testing"
    assert annotations["data_file_division"] == "PLN"
    assert annotations["date"] == "05-MAR-2024"
    assert annotations["molecule_type"] == "DNA"
    assert annotations["location"] == "unknown"
    assert "topology" in annotations
    # Benchling-only keys are not carried over
    assert "folderId" not in annotations
    # batches come from the inventory CSV
    assert annotations["batches"] == [
        {"box": "BoxA", "position": "A1", "volume": 50, "concentration": 120, "location": "BoxA_A1"},
        {"box": "BoxB", "position": "H12", "volume": 20, "concentration": 80, "location": "BoxB_H12"},
    ]


def test_from_benchling_primer_schema(fake_session, frozen_today, inventory_csv):
    dump = _fake_benchling_dump()
    dump["customFields"] = {}
    fake_session.DNASequence.find_by_name.return_value.dump.return_value = dump

    primer = benchling_api.from_benchling("P001_fw", schema="Primer")

    assert isinstance(primer, pydna.primer.Primer)
    assert str(primer.seq) == "ATGCATGCATGCATGCATGC"
    assert primer.name == "P001_fw"
    assert primer.annotations["commentary"] is None
    assert [b["location"] for b in primer.annotations["batches"]] == ["BoxA_A1", "BoxB_H12"]


@pytest.mark.xfail(
    strict=True,
    raises=TypeError,
    reason=(
        "Bug: from_benchling passes strand= to Bio.SeqFeature.SeqFeature, "
        "which Biopython >= 1.82 no longer accepts"
    ),
)
def test_from_benchling_converts_features(fake_session, frozen_today, inventory_csv):
    features = [
        {"start": 2, "end": 8, "strand": 1, "type": "CDS", "name": "geneA", "color": "#ff0000"},
        # wraps around the origin -> compound location
        {"start": 15, "end": 3, "strand": -1, "type": "misc_feature", "name": "wrap", "color": "#00ff00"},
    ]
    fake_session.DNASequence.find_by_name.return_value.dump.return_value = (
        _fake_benchling_dump(features=features)
    )

    record = benchling_api.from_benchling("pTEST")

    gene, wrap = record.features
    assert gene.type == "CDS"
    assert (int(gene.location.start), int(gene.location.end), gene.location.strand) == (2, 8, 1)
    assert gene.qualifiers == {"name": "geneA", "color": "#ff0000", "label": "geneA"}
    assert wrap.type == "misc_feature"
    assert len(wrap.location.parts) == 2
    assert wrap.qualifiers["label"] == "wrap"
