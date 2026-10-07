#!/usr/bin/env python

# The LIMS helpers moved to teemi.legacy.lims; the old teemi.lims import paths
# are kept (with a DeprecationWarning) because the published notebooks use them.

import importlib
import sys
from unittest import mock

import pytest


def _fresh_import(monkeypatch, name):
    monkeypatch.delitem(sys.modules, name, raising=False)
    with mock.patch("dotenv.load_dotenv"), mock.patch(
        "dotenv.find_dotenv", return_value=""
    ), mock.patch("benchlingapi.Session"):
        with pytest.warns(DeprecationWarning, match="teemi.legacy.lims"):
            return importlib.import_module(name)


def test_csv_database_shim_reexports_legacy_module(monkeypatch):
    from teemi.legacy.lims import csv_database as legacy

    shim = _fresh_import(monkeypatch, "teemi.lims.csv_database")

    for name in [
        "add_annotations",
        "add_sequences_to_dataframe",
        "add_unique_ids",
        "get_box",
        "get_database",
        "get_dna_from_box_name",
        "get_dna_from_plate_name",
        "get_unique_id",
        "update_database",
    ]:
        assert getattr(shim, name) is getattr(legacy, name)


def test_benchling_api_shim_reexports_legacy_module(monkeypatch):
    from teemi.legacy.lims import benchling_api as legacy

    shim = _fresh_import(monkeypatch, "teemi.lims.benchling_api")

    for name in ["sequence_to_benchling", "from_benchling", "update_loc_vol_conc"]:
        assert getattr(shim, name) is getattr(legacy, name)
