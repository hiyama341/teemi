#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)

"""Tests for :mod:`teemi.design.combinatorial.report`: :func:`plan`, the
result, and the three tables."""

import copy
import json
import random
from pathlib import Path

import pytest

from teemi.design.combinatorial.report import (
    _shared_id,
    COLUMN_LABELS,
    PCR_COLUMNS,
    PRIMER_COLUMNS,
    InfeasibleDesignError,
    pcr_dataframe,
    plan,
    primer_dataframe,
    variant_overview_dataframe,
)
from teemi.design.combinatorial.design import STRUCTURE_BOUNDS
from teemi.design.combinatorial.spec import ID_SEPARATOR
from teemi.design.combinatorial.sources import SourceResolver
from teemi.design.combinatorial.spec import load_spec
from teemi.design.combinatorial.tm import (
    ENGINE_LOCAL_TA,
    ENGINE_LOCAL_TM,
    ENGINE_NEB,
    local_ta,
    local_tm,
)

from tests.combinatorial_netguard import no_network  # noqa: F401

REPO = Path(__file__).resolve().parents[1]
# The suite owns its inputs: ``notebooks/`` is gitignored scratch the user
# edits freely, so a test that read it would break on an unrelated edit and
# would not run at all in a fresh clone.
FIXTURES = REPO / "tests" / "files_for_testing"
SPEC_B = FIXTURES / "combinatorial_spec_optimized_overhang.json"
SPEC_A = FIXTURES / "combinatorial_spec_fixed_overhang.json"
CIRCULAR = REPO / "tests" / "files_for_testing" / "combinatorial_spec_example.json"


def raw(path: Path) -> dict:
    return json.loads(path.read_text())


@pytest.fixture(scope="module")
def result():
    return plan(
        raw(SPEC_B),
        sources=SourceResolver(base_dir=REPO),
        offline=True,
        include_designs=True,
    )


@pytest.fixture(scope="module")
def circular_result():
    return plan(raw(CIRCULAR), sources=SourceResolver(base_dir=REPO), offline=True)


# --------------------------------------------------------------------------- #
# The result
# --------------------------------------------------------------------------- #


def test_plan_designs_the_whole_library(result):
    assert len(result.variants) == 8
    assert all(v.verification.ok for v in result.variants)
    assert len(result.footprints) == 6
    assert len(result.junctions) == 8
    assert all(j.feasible for j in result.junctions)
    assert len(result.primers) == 16
    assert len(result.pcrs) == 14


def test_an_oligo_is_identified_by_its_sequence_not_its_name(result):
    """Reuse is found, not assumed: one row, one tube, one price per molecule.

    Two of the CDS parts in this library share their 3' end, so the junctions
    either side of them need the same molecules.  Each such molecule has to
    appear once, carrying every name that asked for it, and the PCR table has
    to reference the oligo that serves it rather than a name of its own.
    """
    sequences = [p.sequence for p in result.primers]
    assert len(set(sequences)) == len(sequences)

    reused = [p for p in result.primers if len(p.uses) > 1]
    assert reused, "this library does contain reusable oligos"
    assert sum(len(p.uses) for p in result.primers) == 20, "every name accounted for"

    # A molecule one junction needs keeps that junction's name; one that several
    # need is named for all of them, so no reaction is sent to a tube labelled
    # with a part it does not contain.
    for primer in result.primers:
        if len(primer.uses) == 1:
            assert primer.id == primer.uses[0]
        else:
            for use in primer.uses:
                part = use.split(ID_SEPARATOR)[0]
                assert part in primer.id, (primer.id, use)
    assert any("+" in p.id for p in reused)

    known = {p.id for p in result.primers}
    for pcr in result.pcrs:
        assert pcr.forward_id in known and pcr.reverse_id in known
    assert len({(p.part_id, p.forward_id, p.reverse_id) for p in result.pcrs}) == len(
        result.pcrs
    )


def test_the_price_counts_each_molecule_once(result):
    """What you pay is what you order, so a reused oligo is paid for once."""
    expected = sum(
        len(sequence) * result.spec.objective.price_per_base
        for sequence in {p.sequence for p in result.primers}
    )
    assert result.price == pytest.approx(expected)
    assert result.price < sum(
        len(p.sequence) * result.spec.objective.price_per_base
        for p in result.primers
        for _ in p.uses
    )


def test_a_reaction_shows_the_breaches_of_its_own_primers(result):
    """A reaction carrying a flagged oligo cannot read as unflagged.

    The two tables measure different things -- a PCR row has the *pair*
    (heterodimer, 3' dimer), a primer row has the single molecule (hairpin,
    homodimer) -- and a reader comparing them concluded, reasonably, that they
    disagreed.  They do not: the PCR row now carries the worse of its two
    primers' homodimers and is flagged if either primer is.
    """
    by_id = {p.id: p for p in result.primers}
    for pcr in result.pcrs:
        pair = [by_id[pcr.forward_id], by_id[pcr.reverse_id]]
        row = pcr.to_dict()
        assert row["homodimer_dg_worst"] == min(
            p.structures.homodimer_dg for p in pair
        )
        if any(p.to_dict()["dg_flagged"] for p in pair):
            assert row["dg_flagged"], pcr.part_id


def test_homology_is_exactly_the_requested_length(result):
    requested = result.spec.primer_design.junction.homology_length
    for junction in result.junctions:
        assert len(junction.homology) == requested
        assert junction.chosen.t_f + junction.chosen.t_r == requested


def test_no_primer_exceeds_the_orderable_cap(result):
    cap = result.spec.constraints.primer_len_max
    assert max(len(p.sequence) for p in result.primers) <= cap


def test_to_dict_serialises_under_a_bare_json_dumps(result):
    text = json.dumps(result.to_dict())
    assert json.loads(text)["design"]["topology"] == "linear"


def test_to_dict_carries_what_a_caller_reads(result):
    payload = result.to_dict()
    assert set(payload) == {
        "design",
        "feasible",
        "cost",
        "variants",
        "footprints",
        "junctions",
        "primers",
        "pcrs",
        "diagnostics",
    }
    assert payload["variants"][0]["id"] == result.variants[0].id
    for diagnostic in payload["diagnostics"]:
        assert diagnostic["severity"] in ("info", "warning", "error", "critical")
        assert diagnostic["message"]
        assert "code" in diagnostic and "constraint" in diagnostic


def test_plan_is_deterministic_and_include_designs_does_not_change_the_payload():
    """``hash(spec)`` is a usable cache key only if the result is a function of
    the spec alone."""
    spec = load_spec(raw(SPEC_B))
    first = plan(spec, sources=SourceResolver(base_dir=REPO), offline=True)
    second = plan(
        spec,
        sources=SourceResolver(base_dir=REPO),
        offline=True,
        include_designs=True,
    )
    assert json.dumps(first.to_dict(), sort_keys=True) == json.dumps(
        second.to_dict(), sort_keys=True
    )
    assert hash(spec) == hash(load_spec(raw(SPEC_B)))


def test_every_variant_is_verified_even_without_include_designs(result):
    without = plan(raw(SPEC_B), sources=SourceResolver(base_dir=REPO), offline=True)
    assert [v.verification.verdict for v in without.variants] == [
        v.verification.verdict for v in result.variants
    ]
    assert without.designs == ()
    assert len(result.designs) == 8


def test_feasible_errors_and_raise_for_status_agree(circular_result):
    assert circular_result.errors
    assert circular_result.feasible is False
    with pytest.raises(InfeasibleDesignError) as excinfo:
        circular_result.raise_for_status()
    assert excinfo.value.diagnostics == circular_result.errors


def test_summary_names_the_design_and_the_tallies(result):
    text = result.summary()
    assert result.spec.design.name in text
    assert "8 variants" in text
    assert "8 verified" in text
    assert result.spec.objective.currency in text


def test_price_is_the_ordered_bases_times_the_spec_price(result):
    bases = sum(len(p.sequence) for p in result.primers)
    assert result.price == pytest.approx(bases * result.spec.objective.price_per_base)
    assert result.to_dict()["cost"]["bases"] == bases


# --------------------------------------------------------------------------- #
# The one input boundary
# --------------------------------------------------------------------------- #


def test_an_invalid_spec_becomes_a_critical_diagnostic_not_an_exception():
    result = plan({"design": {"name": "x"}})
    assert result.spec is None
    assert result.feasible is False
    assert [d.code for d in result.diagnostics] == ["invalid_spec"]
    assert "positions" in result.diagnostics[0].message
    assert result.summary().startswith("no design")
    assert result.to_dict()["design"] is None
    assert variant_overview_dataframe(result).empty


def test_a_path_escaping_base_dir_becomes_a_critical_diagnostic():
    escaping = copy.deepcopy(raw(SPEC_B))
    escaping["positions"][0]["parts"][0]["source"]["path"] = "../../etc/passwd"
    result = plan(escaping, sources=SourceResolver(base_dir=REPO), offline=True)
    assert [d.code for d in result.diagnostics] == ["unresolved_source"]
    assert result.diagnostics[0].part_id == "pMLS1"
    assert result.variants == ()


def test_a_part_too_short_for_a_footprint_becomes_a_critical_diagnostic():
    short = {
        "design": {"name": "short"},
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
                    {
                        "id": "b1",
                        "source": {"type": "inline", "sequence": "ACGT" * 100},
                    }
                ],
            },
        ],
    }
    result = plan(short, offline=True)
    assert [d.code for d in result.diagnostics] == ["no_footprint"]
    assert result.diagnostics[0].part_id == "a1"


def test_a_homology_longer_than_the_parts_becomes_a_critical_diagnostic():
    """Two parts that between them are shorter than the requested overlap stop
    the plan at the one boundary, named, rather than raising out of it.

    39 + 39 bp of template cannot supply ``t_f + t_r == 80``, whatever the
    split; both parts still carry footprints (13-17 nt here), so what is
    impossible is the junction and not the footprint.
    """
    result = plan(
        {
            "design": {"name": "homology past the parts"},
            "positions": [
                {
                    "id": slot,
                    "parts": [
                        {"id": slot, "source": {"type": "inline", "sequence": seq}}
                    ],
                }
                for slot, seq in (
                    ("a", "CCGGACCGCGACACGCGAGATACGGTGTGCAGGCGATCA"),
                    ("b", "CCAGGAGGTCTTTGACGCTGTGCACCCACGGGACGGAGT"),
                )
            ],
            "primer_design": {
                "target_tm": 57.0,
                "footprint_len_min": 10,
                "junction": {"homology_length": 80},
            },
            "constraints": {"footprint_tm_window": {"min": 52.0, "max": 85.0}},
        },
        offline=True,
    )
    assert [d.code for d in result.diagnostics] == ["impossible_junction"]
    assert result.diagnostics[0].severity == "critical"
    assert result.diagnostics[0].junction_id == "a__b"
    assert "80" in result.diagnostics[0].message
    assert result.feasible is False
    assert result.variants == ()


def test_file_sources_need_a_resolver_with_a_base_dir():
    result = plan(raw(SPEC_B))
    assert [d.code for d in result.diagnostics] == ["unresolved_source"]
    assert "base_dir" in result.diagnostics[0].message


def test_offline_makes_no_network_request(monkeypatch):
    import requests

    def refuse(*args, **kwargs):
        raise AssertionError("offline=True must issue no request")

    monkeypatch.setattr(requests, "get", refuse)
    monkeypatch.setattr(requests, "post", refuse)
    monkeypatch.setattr(requests.Session, "request", refuse)
    result = plan(raw(SPEC_A), sources=SourceResolver(base_dir=REPO), offline=True)
    assert len(result.variants) == 8
    assert all(p.footprint_tm.engine == ENGINE_LOCAL_TM for p in result.primers)


# --------------------------------------------------------------------------- #
# Constraint reporting
# --------------------------------------------------------------------------- #


def test_a_breach_is_reported_against_the_junction_whose_split_caused_it(
    circular_result,
):
    junction_breaches = [d for d in circular_result.diagnostics if d.junction_id]
    assert junction_breaches
    for diagnostic in junction_breaches:
        assert diagnostic.constraint is not None
        assert diagnostic.constraint.measured_by
        assert diagnostic.junction_id in {j.id for j in circular_result.junctions}


def test_an_untailed_terminal_primer_is_checked_too(result):
    """A linear design's two end primers are bare footprints that no junction
    produced, so nothing else would check their length or structure."""
    flagged = {
        d.message
        for d in result.diagnostics
        if d.code == "homodimer_dg_min" and d.part_id == "tADH1"
    }
    assert flagged, "tADH1's untailed reverse primer breaches the -6.0 floor"
    terminal = next(p for p in result.primers if p.id == "tADH1_rv")
    assert terminal.violations
    assert terminal.sequence == terminal.footprint


def test_a_breach_is_reported_rather_than_an_unorderable_primer_shipped():
    """A footprint forced long by a high target Tm pushes the tailed oligo past
    ``primer_len_max``; the design comes back saying so."""
    long_footprints = copy.deepcopy(raw(SPEC_B))
    long_footprints["primer_design"]["target_tm"] = 72.0
    long_footprints["constraints"]["footprint_tm_window"] = {
        "min": 52.0,
        "max": 80.0,
    }
    result = plan(long_footprints, sources=SourceResolver(base_dir=REPO), offline=True)
    breaches = [d for d in result.diagnostics if d.code == "primer_len_max"]
    assert breaches
    assert all(
        d.constraint.limit == 60.0 and d.constraint.observed > 60.0 for d in breaches
    )


def test_a_verification_failure_is_an_error_naming_its_variant():
    """Two parts with the same sequence make the homology ambiguous, so the
    library would scramble; the plan says which variant and why."""
    rng = random.Random(5)
    unique = "".join(rng.choice("ACGTACGTGC") for _ in range(400))
    repeated = "".join(rng.choice("ACGTACGTGC") for _ in range(400))
    result = plan(
        {
            "design": {"name": "ambiguous homology"},
            "positions": [
                {
                    "id": slot,
                    "parts": [
                        {
                            "id": slot,
                            "source": {"type": "inline", "sequence": sequence},
                        }
                    ],
                }
                for slot, sequence in (
                    ("a", unique),
                    ("b", repeated),
                    ("c", repeated),
                )
            ],
        },
        offline=True,
    )
    failures = [d for d in result.diagnostics if d.code.startswith("assembly_")]
    assert [d.code for d in failures] == ["assembly_multiple_products"]
    assert failures[0].variant_id == result.variants[0].id
    assert result.feasible is False


# --------------------------------------------------------------------------- #
# The three tables
# --------------------------------------------------------------------------- #


def test_variant_overview_has_one_column_per_position_headed_by_its_name(result):
    frame = variant_overview_dataframe(result)
    assert list(frame.columns) == [
        "variant_id",
        "Promoter",
        "CDS",
        "Terminator",
        "length_bp",
        "assembly",
    ]
    assert len(frame) == 8
    assert set(frame["Promoter"]) == {"pMLS1", "pURA1"}
    assert set(frame["assembly"]) == {"verified"}


def test_variant_overview_headings_follow_the_spec(result):
    renamed = copy.deepcopy(raw(SPEC_B))
    renamed["positions"][0]["name"] = "Upstream regulator"
    other = plan(renamed, sources=SourceResolver(base_dir=REPO), offline=True)
    assert "Upstream regulator" in variant_overview_dataframe(other).columns


def headings(columns):
    return [COLUMN_LABELS.get(column, column) for column in columns]


def test_the_primer_table_is_readable_by_default(result):
    frame = primer_dataframe(result)
    assert list(frame.columns) == headings(PRIMER_COLUMNS)
    assert len(frame.columns) <= 16
    assert len(frame) == len(result.primers)


def test_the_pcr_table_is_readable_by_default(result):
    frame = pcr_dataframe(result)
    assert list(frame.columns) == headings(PCR_COLUMNS)
    assert len(frame.columns) <= 16
    assert len(frame) == len(result.pcrs)


def _quantity_keys(columns):
    """The wire keys that carry a dG or a temperature, engines excluded."""
    return [
        key
        for key in columns
        if ("_dg" in key or "_tm" in key or key.startswith("ta_"))
        and not key.endswith(("_engine", "_flagged"))
    ]


@pytest.mark.parametrize("columns", [PRIMER_COLUMNS, PCR_COLUMNS])
def test_a_dg_or_tm_column_cannot_reach_a_table_without_its_unit(columns):
    """A bare number in a table is a number read in the wrong scale sooner or
    later: primer3 answers in cal/mol, NEB in C, and the floors are kcal/mol.

    Driven from the WIRE KEYS, not from the headings.  Checking the headings
    only ever asserted that labels which exist carry a unit -- a quantity whose
    label was forgotten keeps its bare lowercase key, matches no heading
    pattern, and would sail through.  This way a new quantity column is a
    failing test until it is given a unit.
    """
    quantities = _quantity_keys(columns)
    assert quantities, "the projection has to contain some quantity at all"
    for key in quantities:
        assert key in COLUMN_LABELS, f"{key} reaches a table with no unit"
        heading = COLUMN_LABELS[key]
        expected = "(kcal/mol)" if "_dg" in key else "(C)"
        assert heading.endswith(expected), (key, heading)


def test_a_dg_or_tm_heading_carries_its_unit(result):
    """The same guarantee seen from the rendered frame."""
    for frame in (
        primer_dataframe(result),
        pcr_dataframe(result),
    ):
        for column in frame.columns:
            if "deltaG" in column:
                assert column.endswith("(kcal/mol)"), column
            elif "Tm" in column or column.startswith("Ta"):
                assert column.endswith("(C)"), column


def test_no_two_columns_share_a_heading():
    """A label collision would silently merge two quantities into one column."""
    headings = [COLUMN_LABELS.get(c, c) for c in (*PRIMER_COLUMNS, *PCR_COLUMNS)]
    for group in (PRIMER_COLUMNS, PCR_COLUMNS):
        rendered = [COLUMN_LABELS.get(c, c) for c in group]
        assert len(set(rendered)) == len(rendered), rendered
    # And no label is the bare key of some other column.
    assert not set(COLUMN_LABELS.values()) & set(headings) - set(
        COLUMN_LABELS.values()
    )


def test_full_tables_are_a_superset_of_the_readable_ones(result):
    for narrow, wide in (
        (primer_dataframe(result), primer_dataframe(result, full=True)),
        (pcr_dataframe(result), pcr_dataframe(result, full=True)),
    ):
        assert set(narrow.columns) < set(wide.columns)


@pytest.mark.parametrize("full", [False, True])
def test_no_table_column_conflates_the_two_tm_quantities(result, full):
    """The bug class this package exists to prevent: a footprint Tm and a whole
    oligo Tm printed as one quantity, or a Ta whose name does not say which pair
    it came from."""
    for frame in (
        primer_dataframe(result, full=full),
        pcr_dataframe(result, full=full),
    ):
        for column in frame.columns:
            lower = column.lower()
            assert lower not in ("tm", "ta", "tm_c", "ta_c", "tm (c)", "ta (c)")
            # Case-insensitive on purpose: the headings spell it ``Tm`` and the
            # wire keys ``tm``, and a check that silently stopped matching one
            # of them would stop guarding this bug class at all.
            if "tm" in lower and "dg" not in lower and "deltag" not in lower:
                assert lower.startswith(
                    (
                        "footprint_tm",
                        "oligo_tm",
                        "hairpin_tm",
                        "homodimer_tm",
                        "heterodimer_tm",
                        "heterodimer_fwd_rev_tm",
                        "tm_engine",
                    )
                ), column


def test_every_tm_column_is_accompanied_by_its_engine(result):
    """One ``tm_engine`` heading per row, covering both Tm columns.

    Both Tm of a row are measured by the same engine, so a second identical
    column said nothing; the engine is still on every record for provenance,
    and Ta keeps its own because its model is a different one.
    """
    primers = primer_dataframe(result)
    assert set(primers["tm_engine"]) == {ENGINE_LOCAL_TM}
    assert "oligo_tm_engine" not in primers.columns
    pcrs = pcr_dataframe(result)
    assert set(pcrs["ta_engine"]) == {ENGINE_LOCAL_TA}
    assert set(pcrs["tm_engine"]) == {ENGINE_LOCAL_TM}
    # Dropped from the readable table, never from the record.
    for row in (result.primers[0].to_dict(), result.pcrs[0].to_dict()):
        assert row["oligo_tm_engine"] == ENGINE_LOCAL_TM


def test_ta_comes_from_the_footprint_pair_not_the_tailed_oligos(result):
    """NEB's recommendation for a 50-60 nt tailed oligo runs well above the
    footprints it is labelled as deriving from, and a PCR run at it fails."""
    by_id = {p.id: p for p in result.primers}
    tailed_differs = 0
    for pcr in result.pcrs:
        forward, reverse = by_id[pcr.forward_id], by_id[pcr.reverse_id]
        from_footprints = local_ta(forward.footprint, reverse.footprint)
        from_oligos = local_ta(forward.sequence, reverse.sequence)
        assert pcr.temperatures.ta_from_footprints.celsius == pytest.approx(
            from_footprints.celsius
        )
        if abs(from_oligos.celsius - from_footprints.celsius) > 1.0:
            tailed_differs += 1
    assert tailed_differs, "the tailed oligos would have given a different Ta"


def test_the_pcr_amplicon_length_is_the_part_plus_the_two_tails(result):
    for pcr in result.pcrs:
        forward = next(p for p in result.primers if p.id == pcr.forward_id)
        reverse = next(p for p in result.primers if p.id == pcr.reverse_id)
        part = next(f for f in result.footprints if f.part_id == pcr.part_id)
        tails = (len(forward.sequence) - len(forward.footprint)) + (
            len(reverse.sequence) - len(reverse.footprint)
        )
        assert pcr.amplicon_bp - tails > 0
        assert part.part_id == pcr.part_id


def test_tables_are_empty_rather_than_broken_when_there_is_no_design():
    result = plan({"design": {"name": "x"}})
    assert primer_dataframe(result).empty
    assert pcr_dataframe(result).empty
    assert variant_overview_dataframe(result).empty


# --------------------------------------------------------------------------- #
# The NEB path
# --------------------------------------------------------------------------- #


class _Reply:
    """A requests-shaped reply whose numbers are functions of the sequences.

    Deriving every number from the sequence itself, rather than from the row's
    position, is what lets the tests below tell a footprint row from an oligo
    row without assuming the order the request was built in.
    """

    status_code = 200

    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return self._payload


def _stub_tm(seq: str) -> float:
    return 40.0 + len(seq) / 10.0


def _stub_ta(first: str, second: str) -> float:
    return 20.0 + len(first) + len(second)


@pytest.fixture
def neb(monkeypatch):
    """NEB answered by a stub; the counts are what went out."""
    import requests

    from teemi.design.combinatorial.tm import clear_tm_caches

    clear_tm_caches()
    calls = {"get": 0, "post": 0}

    def get(url, params=None, headers=None, **kwargs):
        calls["get"] += 1
        return _Reply({"success": True, "data": {"tm1": _stub_tm(params["seq1"])}})

    def post(url, json=None, timeout=None, **kwargs):
        calls["post"] += 1
        return _Reply(
            {
                "success": True,
                "data": [
                    {
                        "tm1": _stub_tm(first),
                        "tm2": _stub_tm(second),
                        "ta": _stub_ta(first, second),
                    }
                    for first, second in json["seqpairs"]
                ],
            }
        )

    monkeypatch.setattr(requests, "get", get)
    monkeypatch.setattr(requests, "post", post)
    yield calls
    clear_tm_caches()


def test_online_the_authoritative_numbers_cost_one_batched_post(neb):
    result = plan(raw(SPEC_A), sources=SourceResolver(base_dir=REPO), offline=False)
    assert neb["post"] == 1, "every PCR's Tm and Ta travel in one request"
    assert neb["get"] > 0
    assert {p.footprint_tm.engine for p in result.primers} == {ENGINE_NEB}
    assert {p.oligo_tm.engine for p in result.primers} == {ENGINE_NEB}


def test_online_ta_is_read_off_the_footprint_rows(neb):
    """The bug this guards: Ta taken from the tailed-oligo rows, which read well
    above the footprints the column says it derives from."""
    result = plan(raw(SPEC_A), sources=SourceResolver(base_dir=REPO), offline=False)
    by_id = {p.id: p for p in result.primers}
    for pcr in result.pcrs:
        forward, reverse = by_id[pcr.forward_id], by_id[pcr.reverse_id]
        from_footprints = _stub_ta(forward.footprint, reverse.footprint)
        from_oligos = _stub_ta(forward.sequence, reverse.sequence)
        assert pcr.temperatures.ta_from_footprints.celsius == pytest.approx(
            from_footprints
        )
        if from_oligos != from_footprints:
            assert pcr.temperatures.ta_from_footprints.celsius != pytest.approx(
                from_oligos
            )


def test_online_the_two_tm_quantities_are_read_off_their_own_rows(neb):
    result = plan(raw(SPEC_A), sources=SourceResolver(base_dir=REPO), offline=False)
    by_id = {p.id: p for p in result.primers}
    for pcr in result.pcrs:
        forward, reverse = by_id[pcr.forward_id], by_id[pcr.reverse_id]
        temps = pcr.temperatures
        assert temps.footprint_tm_forward.celsius == pytest.approx(
            _stub_tm(forward.footprint)
        )
        assert temps.oligo_tm_forward.celsius == pytest.approx(
            _stub_tm(forward.sequence)
        )
        assert temps.footprint_tm_reverse.celsius == pytest.approx(
            _stub_tm(reverse.footprint)
        )
        assert temps.oligo_tm_reverse.celsius == pytest.approx(
            _stub_tm(reverse.sequence)
        )


def test_an_unreachable_neb_is_a_warning_and_local_numbers_not_an_exception(
    monkeypatch,
):
    import requests

    from teemi.design.combinatorial.tm import clear_tm_caches

    clear_tm_caches()

    def refuse(*args, **kwargs):
        raise requests.exceptions.ConnectionError("no route to host")

    monkeypatch.setattr(requests, "get", refuse)
    monkeypatch.setattr(requests, "post", refuse)
    result = plan(raw(SPEC_A), sources=SourceResolver(base_dir=REPO), offline=False)
    clear_tm_caches()

    warnings = [d for d in result.diagnostics if d.code == "neb_unreachable"]
    assert warnings
    assert all(d.severity == "warning" for d in warnings)
    assert len(result.variants) == 8
    assert {p.footprint_tm.engine for p in result.primers} == {ENGINE_LOCAL_TM}


# --------------------------------------------------------------------------- #
# One engine per comparison, over a whole plan
# --------------------------------------------------------------------------- #

# Which source is allowed to have measured BOTH sides of each bound.  A Tm
# bound is the interesting case: its limit is a spec number, but the
# observation is only meaningful against it when the engine that produced the
# footprints produced it too.  Anything else is a property of the molecule
# (counted off the string) or of primer3's structure model.
TM_ENGINE_BOUNDS = frozenset(
    {"footprint_tm_window.min", "footprint_tm_window.max", "pair_tm_delta_max"}
)
PRIMER3_BOUNDS = frozenset(
    {
        "hairpin_dg_min",
        "hairpin_tm_max",
        "homodimer_dg_min",
        "homodimer_tm_max",
        "heterodimer_dg_min",
        "heterodimer_tm_max",
        "three_prime_dimer_dg_min",
    }
)
SEQUENCE_BOUNDS = frozenset({"primer_len_max", "gc_fraction_min", "gc_fraction_max"})


def _all_violations(result):
    """Every breach the plan recorded, wherever it was recorded."""
    found = [d.constraint for d in result.diagnostics if d.constraint is not None]
    for junction in result.junctions:
        found.extend(junction.chosen.violations)
    for primer in result.primers:
        found.extend(primer.violations)
    for pcr in result.pcrs:
        found.extend(pcr.violations)
    return found


def _neb_offset_by_oligo(monkeypatch):
    """NEB reading below pydna by an oligo-dependent amount.

    The offset has to depend on the sequence: a *difference* of two Tms is
    invariant under a constant shift, so with a flat offset a pydna delta and a
    NEB delta would be the same number and no test could tell which engine
    produced the one that got reported.
    """
    import requests

    def get(url, params=None, headers=None, **kwargs):
        seq = params["seq1"]
        celsius = local_tm(seq).celsius - 4.0 - 0.05 * len(seq)
        return _Reply({"success": True, "data": {"tm1": celsius}})

    def post(url, json=None, timeout=None, **kwargs):
        rows = []
        for first, second in json["seqpairs"]:
            rows.append(
                {
                    "tm1": local_tm(first).celsius - 4.0 - 0.05 * len(first),
                    "tm2": local_tm(second).celsius - 4.0 - 0.05 * len(second),
                    "ta": _stub_ta(first, second),
                }
            )
        return _Reply({"success": True, "data": rows})

    monkeypatch.setattr(requests, "get", get)
    monkeypatch.setattr(requests, "post", post)


@pytest.mark.parametrize("offline", [True, False])
def test_no_bound_is_ever_checked_against_another_engines_number(offline, monkeypatch):
    """The defect that shipped three times, caught over a whole library.

    Every breach names the one source behind both its observation and its
    limit, and for a Tm bound that source must be the engine that actually
    designed this plan's footprints -- pydna offline, NEB online.  The final
    ``else`` is the point of the test: a bound added later cannot be recorded
    without being classified here, so nobody can quietly introduce a
    comparison whose two sides come from different models.
    """
    tight = copy.deepcopy(raw(SPEC_B))
    target = tight["primer_design"]["target_tm"]
    tight["constraints"] = {
        # Straddling target_tm, and derived from it rather than written out: the
        # footprint walk stops at the first prefix to REACH the target, so it
        # overshoots a ceiling set this close and the window bound is exercised
        # under either engine.  A literal here would silently stop exercising
        # it the moment the spec's target moved.
        "footprint_tm_window": {"min": target - 1.0, "max": target + 0.2},
        "pair_tm_delta_max": 0.2,
        "primer_len_max": 45,
        "gc_fraction_min": 0.45,
        "gc_fraction_max": 0.55,
        "hairpin_dg_min": -0.5,
        "homodimer_dg_min": -0.5,
        "heterodimer_dg_min": -0.5,
        "three_prime_dimer_dg_min": -0.5,
        "hairpin_tm_max": 10.0,
        "homodimer_tm_max": 10.0,
        "heterodimer_tm_max": 10.0,
    }
    from teemi.design.combinatorial.tm import clear_tm_caches

    clear_tm_caches()
    if not offline:
        _neb_offset_by_oligo(monkeypatch)
    result = plan(tight, sources=SourceResolver(base_dir=REPO), offline=offline)
    clear_tm_caches()

    engines = {p.footprint_tm.engine for p in result.primers}
    assert engines == {ENGINE_LOCAL_TM if offline else ENGINE_NEB}
    footprint_engine = engines.pop()

    violations = _all_violations(result)
    assert violations, "the tightened spec has to breach something"
    seen = set()
    for violation in violations:
        seen.add(violation.constraint)
        assert violation.measured_by, violation
        if violation.constraint in TM_ENGINE_BOUNDS:
            assert violation.measured_by == footprint_engine, violation
        elif violation.constraint in PRIMER3_BOUNDS:
            assert violation.measured_by == "primer3", violation
        elif violation.constraint in SEQUENCE_BOUNDS:
            assert violation.measured_by == "sequence", violation
        else:
            raise AssertionError(
                f"{violation.constraint!r} is checked but no test says which "
                "engine measured both sides of it"
            )
    # Not a vacuous pass: a mislabelled bound is only caught if it breached,
    # so every source, and both Tm-engine bounds, have to be exercised here.
    assert "pair_tm_delta_max" in seen
    assert seen & {"footprint_tm_window.min", "footprint_tm_window.max"}
    assert seen & PRIMER3_BOUNDS and seen & SEQUENCE_BOUNDS


# --------------------------------------------------------------------------- #
# A breach nothing could have avoided is a warning, not a blocker
# --------------------------------------------------------------------------- #


def test_a_breach_no_option_avoids_is_a_warning_not_an_error(result):
    """An optimal design must not report itself infeasible.

    ``tADH1``'s reverse primer is the terminator's own 3' end: ``GCATGC`` is
    self-complementary and sits inside every footprint length, there is no 5'
    tail at a terminal end to re-split, and trimming a genetic part is not
    allowed.  The -6.0 floor is still breached and still shown; what changed is
    that it no longer blocks, because no input could have done better.
    """
    assert result.feasible is True
    assert result.errors == ()

    warnings = [d for d in result.diagnostics if d.severity == "warning"]
    assert warnings, "the two terminal homodimers still have to be reported"
    for diagnostic in warnings:
        assert diagnostic.code in STRUCTURE_BOUNDS
        assert "no option avoids it" in diagnostic.message
        assert diagnostic.constraint is not None
        # The floor itself is untouched.
        assert diagnostic.constraint.limit == -6.0
        assert diagnostic.constraint.observed < -6.0

    # Still visible in both tables, and marked as not worth acting on.
    for frame in (primer_dataframe(result), pcr_dataframe(result)):
        flagged = frame[frame["flagged_by"] != ""]
        assert not flagged.empty
        assert set(flagged["flag_severity"]) == {"warning"}


def test_a_breach_some_other_split_avoids_stays_an_error(circular_result):
    """Exhaustion is proved per bound, never assumed for a whole junction.

    In the circular library every primer is tailed, so the sweep has 36 splits
    to choose from.  Where some candidate avoided a bound and the chosen one
    still breaches it, the optimizer made a trade-off a caller can influence,
    and that stays blocking.
    """
    assert circular_result.feasible is False
    severities = {d.severity for d in circular_result.diagnostics}
    assert severities >= {"warning", "error"}, severities

    by_junction = {j.id: j for j in circular_result.junctions}
    errors = [
        d
        for d in circular_result.diagnostics
        if d.severity == "error" and d.junction_id and d.code in STRUCTURE_BOUNDS
    ]
    assert errors, "this fixture has to contain an avoidable structure breach"
    for diagnostic in errors:
        junction = by_junction[diagnostic.junction_id]
        assert diagnostic.code not in junction.unavoidable
        assert "no option avoids it" not in diagnostic.message
        # Some candidate really did avoid it, so there was a choice.
        assert len(junction.profile) > 1


def test_unavoidable_only_ever_names_a_bound_the_sweep_owns(circular_result):
    """A length or GC breach is not something a split was going to fix."""
    for junction in circular_result.junctions:
        assert set(junction.unavoidable) <= STRUCTURE_BOUNDS
        breached = {v.constraint for v in junction.chosen.violations}
        assert set(junction.unavoidable) <= breached


def test_a_non_structure_breach_always_blocks():
    """The escape hatch is for secondary structure and nothing else.

    An unorderable primer is unorderable however hard the search tried, so a
    ``primer_len_max`` breach must block whatever the sweep found.
    """
    tight = copy.deepcopy(raw(SPEC_B))
    tight.setdefault("constraints", {})["primer_len_max"] = 40
    result = plan(tight, sources=SourceResolver(base_dir=REPO), offline=True)

    length_breaches = [d for d in result.diagnostics if d.code == "primer_len_max"]
    assert length_breaches, "a 40 nt cap has to be breached by this library"
    assert all(d.severity == "error" for d in length_breaches)
    assert result.feasible is False


@pytest.mark.parametrize("spec_path", [SPEC_A, SPEC_B, CIRCULAR])
def test_only_a_structure_bound_is_ever_downgraded(spec_path):
    """The invariant behind the escape hatch, over every fixture.

    This is the test that keeps the hatch narrow: whatever a future bound is
    called, it cannot be quietly downgraded into a warning unless it is one the
    split sweep is actually responsible for.
    """
    result = plan(raw(spec_path), sources=SourceResolver(base_dir=REPO), offline=True)
    downgraded = [d for d in result.diagnostics if d.severity == "warning" and d.code]
    for diagnostic in downgraded:
        # Warnings that are not breaches at all (an unreachable NEB) carry no
        # constraint; a breach that was downgraded has to be a structure bound.
        if diagnostic.constraint is not None:
            assert diagnostic.code in STRUCTURE_BOUNDS, diagnostic.message


def test_shared_id_never_names_a_junction_the_design_does_not_have():
    """``_shared_id`` merges one token column, or none at all.

    Joining each differing column independently names the cross product: two
    uses differing in both tokens would produce ``A+B__X+Y_rv``, asserting the
    junctions ``(A, Y)`` and ``(B, X)`` that the design need not contain.  A
    differing last column is worse -- the direction suffix rides on it, so the
    name comes out as ``...__X_rv+Y_rv`` and no longer parses at all.  In both
    cases the id falls back to one real name and ``used_as`` carries the truth.
    """
    assert _shared_id(["A__X_rv"]) == "A__X_rv"
    # One column differs: safe to merge, and sorted so it is deterministic.
    assert _shared_id(["B__X_rv", "A__X_rv"]) == "A+B__X_rv"
    # The last column carries the direction, so merging it would give
    # "A__X_rv+Y_rv" -- a name that no longer parses.  Refused.
    assert _shared_id(["A__Y_rv", "A__X_rv"]) == "A__Y_rv"
    # Two columns differ: no merged name could be truthful.
    assert _shared_id(["A__X_rv", "B__Y_rv"]) == "A__X_rv"
    # Different shapes (a bare footprint beside a tailed primer): no splicing.
    assert _shared_id(["X_fw", "A__X_fw"]) == "X_fw"


def test_the_merged_name_is_independent_of_enumeration_order(result):
    """Two runs must put the same oligo under the same name."""
    again = plan(raw(SPEC_B), sources=SourceResolver(base_dir=REPO), offline=True)
    assert [p.id for p in again.primers] == [p.id for p in result.primers]
    for primer in result.primers:
        # Only a merged name is order-independent by construction: the fallback
        # deliberately keeps one real name, and which one is whichever the
        # enumeration reached first -- deterministic for a spec, as the line
        # above checks, but not invariant under reordering ``uses``.
        if "+" in primer.id:
            assert primer.id == _shared_id(list(reversed(primer.uses)))


def test_a_tailed_primer_shows_its_own_breaches_in_its_own_row(circular_result):
    """A row cannot print a homodimer past the floor and call itself clean.

    Tailed primers used to be given ``violations=()`` on the grounds that their
    breaches are reported against their junction.  The diagnostic was there, but
    both tables read every tailed primer as unflagged -- and since the PCR row's
    flag is fed from those same tuples, the whole circular library showed zero
    flags against twelve diagnostics.  Reporting and display are now separate:
    the row always carries its own breaches, the diagnostic is still raised once.
    """
    tailed = [p for p in circular_result.primers if p.tailed]
    assert tailed, "every primer in a circular design is tailed"
    flagged = [p for p in tailed if p.violations]
    assert flagged, "this fixture has tailed primers that breach a bound"

    frame = primer_dataframe(circular_result)
    shown = frame[frame["flagged_by"] != ""]
    assert len(shown) == len(flagged)
    assert set(shown["flag_severity"]) <= {"error", "warning"}

    # Still reported exactly once, against the junction that chose the split.
    for primer in flagged:
        codes = [
            d
            for d in circular_result.diagnostics
            if d.code in {v.constraint for v in primer.violations}
            and d.constraint is not None
            and any(
                abs(d.constraint.observed - v.observed) < 1e-9
                for v in primer.violations
            )
        ]
        assert codes, primer.id
        assert all(d.junction_id for d in codes), [d.to_dict() for d in codes]


def test_flagged_by_always_has_its_detail_in_violations(result, circular_result):
    """Reading ``flagged_by`` then ``violations`` must not come up empty.

    ``PcrRecord`` derives its flag from the pair bounds *and* its two primers'
    own bounds, but only ever serialised the pair ones, so a consumer that saw
    ``flagged_by='homodimer_dg_min'`` found an empty ``violations`` array and no
    way to learn the observed value or the limit.
    """
    for design in (result, circular_result):
        for row in design.to_dict()["pcrs"] + design.to_dict()["primers"]:
            named = {c for c in row["flagged_by"].split(",") if c}
            detailed = {v["constraint"] for v in row["violations"]}
            assert named == detailed, row
            for violation in row["violations"]:
                assert violation["measured_by"]
                assert violation["units"]


def test_a_shared_oligo_is_attributed_to_every_part_it_amplifies(circular_result):
    """``part_id`` has to agree with the row's own id and ``used_as``.

    Taking it from whichever pair happened to introduce the sequence first
    attributed a molecule serving three templates to one arbitrary part, so
    filtering the order sheet by part silently returned an incomplete set.
    """
    shared = [p for p in circular_result.primers if len(p.uses) > 1]
    assert shared
    by_part = {}
    for pair_list in circular_result.pcrs:
        by_part.setdefault(pair_list.forward_id, set()).add(pair_list.part_id)
        by_part.setdefault(pair_list.reverse_id, set()).add(pair_list.part_id)
    for primer in shared:
        named = set(primer.part_id.split("+"))
        assert named == by_part[primer.id], (primer.id, named, by_part[primer.id])
