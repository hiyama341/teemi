#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)

"""The package's public surface, and the quickstart notebook's contract.

Everything here runs the calls ``notebooks/combinatorial_quickstart.ipynb``
makes, in the order it makes them, so a change that breaks the notebook breaks a
test instead.
"""

import ast
import json
import re
from pathlib import Path

import pytest
from Bio import SeqIO

import teemi.design.combinatorial as cb
from teemi.design.combinatorial.sources import SourceResolver

from tests.combinatorial_netguard import no_network  # noqa: F401

REPO = Path(__file__).resolve().parents[1]
PACKAGE = REPO / "teemi" / "design" / "combinatorial"
NOTEBOOK = REPO / "notebooks" / "combinatorial_quickstart.ipynb"
SPEC_FILES = sorted(
    (REPO / "tests" / "files_for_testing").glob("combinatorial_spec_*overhang.json")
)


#: The notebook lives under gitignored scratch, so a fresh clone does not have
#: it and the suite must not depend on its being there.
needs_notebook = pytest.mark.skipif(
    not NOTEBOOK.is_file(),
    reason=f"{NOTEBOOK.relative_to(REPO)} is gitignored scratch and is not here",
)


def notebook_source() -> str:
    cells = json.loads(NOTEBOOK.read_text())["cells"]
    return "\n".join(
        "".join(cell["source"]) for cell in cells if cell["cell_type"] == "code"
    )


# --------------------------------------------------------------------------- #
# The surface
# --------------------------------------------------------------------------- #


def test_the_public_surface_is_small_and_complete():
    assert sorted(cb.__all__) == sorted(
        [
            "plan",
            "DesignResult",
            "InfeasibleDesignError",
            "variant_overview_dataframe",
            "primer_dataframe",
            "pcr_dataframe",
            "write_genbank_files",
            "write_assembly_figures",
            "DesignSpec",
            "load_spec",
            "json_schema",
            "SourceResolver",
            "SourceError",
            "Diagnostic",
            "Violation",
            "Verification",
        ]
    )
    for name in cb.__all__:
        assert hasattr(cb, name), name


@needs_notebook
def test_the_notebook_uses_only_the_public_surface():
    used = set(re.findall(r"\bcb\.([A-Za-z_][A-Za-z0-9_]*)", notebook_source()))
    assert used, "the notebook should exercise the package"
    assert used <= set(cb.__all__), used - set(cb.__all__)


def test_the_package_is_the_seven_modules_the_map_names():
    assert sorted(p.name for p in PACKAGE.glob("*.py")) == [
        "__init__.py",
        "design.py",
        "export.py",
        "report.py",
        "simulate.py",
        "sources.py",
        "spec.py",
        "tm.py",
    ]


def test_there_is_one_schema_with_no_version_and_no_migrations():
    schema = cb.json_schema()
    assert schema["title"] == "DesignSpec"
    assert "schema_version" not in schema["properties"]
    assert set(schema["required"]) == {"design", "positions"}
    assert not [name for name in dir(cb) if "migrat" in name.lower()]


# --------------------------------------------------------------------------- #
# The notebook, cell by cell
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("spec_path", SPEC_FILES, ids=lambda p: p.stem)
def test_the_notebook_runs_against_both_specs(spec_path, tmp_path):
    SPEC = json.loads(spec_path.read_text())

    result = cb.plan(
        SPEC,
        sources=SourceResolver(base_dir=REPO),
        offline=True,
        include_designs=True,
    )
    R = result.to_dict()

    assert isinstance(result.summary(), str) and result.summary()
    for d in R["diagnostics"]:
        assert d["severity"] in ("info", "warning", "error", "critical")
        assert isinstance(d["message"], str)

    overview = cb.variant_overview_dataframe(result)
    primers = cb.primer_dataframe(result)
    pcrs = cb.pcr_dataframe(result)
    assert len(overview) == 8
    assert not primers.empty and not pcrs.empty

    designs = list(result.designs)[:1]
    figures = cb.write_assembly_figures(designs, tmp_path)
    genbank = cb.write_genbank_files(designs, tmp_path)
    assert len(figures) == len(genbank) == 1

    vid = R["variants"][0]["id"]
    assert next(tmp_path.glob(f"{vid}*.txt")).read_text()
    record = SeqIO.read(next(tmp_path.glob(f"{vid}*.gb")), "genbank")
    assert len(record.features) == 9


def test_both_specs_describe_the_same_library_with_different_split_searches():
    """Case A fixes the split at the midpoint; case B searches all 36.  Same
    parts, same homology, different primers."""
    plans = {
        path.stem: cb.plan(
            json.loads(path.read_text()),
            sources=SourceResolver(base_dir=REPO),
            offline=True,
        )
        for path in SPEC_FILES
    }
    a = plans["combinatorial_spec_fixed_overhang"]
    b = plans["combinatorial_spec_optimized_overhang"]
    assert [v.id for v in a.variants] == [v.id for v in b.variants]
    assert {j.chosen.t_f for j in a.junctions} == {18}
    assert {j.chosen.t_f for j in b.junctions} != {18}
    assert all(len(j.profile) == 1 for j in a.junctions)
    assert all(len(j.profile) == 36 for j in b.junctions)


def test_footprints_come_from_the_one_pydna_primer_design_call():
    """One call site, so no second code path can produce a footprint under
    different settings.  Counted over the syntax tree, since the docstrings
    name the call too."""
    calls = []
    for path in PACKAGE.glob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Call):
                func = node.func
                name = getattr(func, "id", None) or getattr(func, "attr", None)
                if name == "primer_design":
                    calls.append(path.name)
    assert calls == ["design.py"], calls
