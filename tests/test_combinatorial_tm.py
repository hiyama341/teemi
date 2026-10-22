#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)
"""Tests for :mod:`teemi.design.combinatorial.tm`.

Two properties carry most of the weight: a number from one model can never be
compared with a number from another, and the annealing temperature comes off
the footprint rows of the one batched request, not the oligo rows.
"""

from __future__ import annotations

import pytest
import requests
from pydna.tm import ta_dbd, tm_default

from teemi.design.combinatorial import tm as T

from tests.combinatorial_netguard import no_network  # noqa: F401

FORWARD = "ATGCGTACGTTAGCACCGTAGCATCGATCGGATCCTAGCA"
REVERSE = "TGCTAGGATCCGATCGATGCTACGGTGCTAACGTACGCAT"
FOOTPRINT_F = "ATGCGTACGTTAGCACCG"
FOOTPRINT_R = "TGCTAGGATCCGATCGAT"


class FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return self._payload


def _forbidden(*args, **kwargs):
    raise AssertionError("a test reached for NEB without stubbing it")


@pytest.fixture(autouse=True)
def fresh_caches(monkeypatch):
    """Empty memos and a poisoned ``tm_neb`` on top of the module-wide guard.

    The shared guard blocks the verbs; this also blocks the pydna function that
    would reach for one, so a test that forgets to stub NEB names the call it
    missed rather than the socket.
    """

    monkeypatch.setattr(T, "tm_neb", _forbidden)
    T.clear_tm_caches()
    yield
    T.clear_tm_caches()


def _pcr(part_id="p", suffix=""):
    return T.PcrOligos(
        part_id=part_id,
        forward_oligo=FORWARD + suffix,
        reverse_oligo=REVERSE + suffix,
        forward_footprint=FOOTPRINT_F,
        reverse_footprint=FOOTPRINT_R,
    )


# --------------------------------------------------------------------------- #
# One engine per comparison
# --------------------------------------------------------------------------- #


def test_tm_cannot_be_subtracted_or_ordered():
    a, b = T.Tm(60.0, T.ENGINE_NEB), T.Tm(58.0, T.ENGINE_NEB)
    with pytest.raises(TypeError):
        a - b
    with pytest.raises(TypeError):
        a < b


def test_delta_within_one_engine():
    a, b = T.Tm(60.0, T.ENGINE_NEB), T.Tm(58.5, T.ENGINE_NEB)
    assert a.delta(b) == pytest.approx(1.5)
    assert b.delta(a) == pytest.approx(1.5)


def test_delta_refuses_two_engines():
    local = T.Tm(60.0, T.ENGINE_LOCAL_TM)
    neb = T.Tm(56.0, T.ENGINE_NEB)
    with pytest.raises(T.TmEngineMismatch):
        local.delta(neb)


def test_require_same_engine():
    assert T.require_same_engine(T.Tm(1.0, "x"), T.Tm(2.0, "x")) == "x"
    with pytest.raises(T.TmEngineMismatch):
        T.require_same_engine(T.Tm(1.0, "x"), T.Tm(2.0, "y"))


# --------------------------------------------------------------------------- #
# Conditions
# --------------------------------------------------------------------------- #


def test_monovalent_total_is_na_plus_k_plus_half_tris():
    assert T.CONDITIONS.monovalent_mm == pytest.approx(77.5)


def test_neb_conc_is_the_same_primer_concentration_in_micromolar():
    assert T.CONDITIONS.neb_conc_um == pytest.approx(0.25)
    assert T.CONDITIONS.pydna_tm_kwargs["dnac1"] == pytest.approx(250.0)


def test_pydna_kwargs_are_pydna_s_own_defaults():
    # An offline run must reproduce pydna's numbers exactly.
    assert tm_default(FORWARD, **T.CONDITIONS.pydna_tm_kwargs) == pytest.approx(
        tm_default(FORWARD)
    )


def test_primer3_kwargs_carry_the_same_salts():
    kwargs = T.CONDITIONS.primer3_kwargs
    assert kwargs["mv_conc"] == pytest.approx(T.CONDITIONS.monovalent_mm)
    assert kwargs["dv_conc"] == pytest.approx(T.CONDITIONS.mg_mm)
    assert kwargs["dna_conc"] == pytest.approx(T.CONDITIONS.primer_nm)


# --------------------------------------------------------------------------- #
# Local models
# --------------------------------------------------------------------------- #


def test_local_tm_delegates_to_pydna():
    value = T.local_tm(FORWARD)
    assert value.engine == T.ENGINE_LOCAL_TM
    assert value.celsius == pytest.approx(
        tm_default(FORWARD, **T.CONDITIONS.pydna_tm_kwargs)
    )


def test_local_ta_delegates_to_pydna_ta_dbd():
    value = T.local_ta(FOOTPRINT_F, FOOTPRINT_R)
    assert value.engine == T.ENGINE_LOCAL_TA
    assert value.celsius == pytest.approx(ta_dbd(FOOTPRINT_F, FOOTPRINT_R, ""))


# --------------------------------------------------------------------------- #
# FootprintTm
# --------------------------------------------------------------------------- #


def test_offline_engine_makes_no_request():
    engine = T.FootprintTm(offline=True)
    # primer3 is the cheap walk either way; what offline changes is the
    # authoritative pass, which becomes pydna's own model instead of NEB.
    assert engine.estimate_function(FORWARD) == pytest.approx(
        T.primer3_tm(FORWARD)
    )
    assert engine.tm(FORWARD).engine == T.ENGINE_LOCAL_TM
    assert engine.failures == ()
    assert engine.engine == T.ENGINE_LOCAL_TM


def test_online_engine_tags_neb_values_and_passes_a_two_stage_pair(monkeypatch):
    calls = []

    def fake_tm_neb(seq, conc, prodcode):
        calls.append((seq, conc, prodcode))
        return 54.5

    monkeypatch.setattr(T, "tm_neb", fake_tm_neb)
    engine = T.FootprintTm(offline=False)
    assert engine.estimate_function is not None
    value = engine.tm(FORWARD)
    assert value == T.Tm(54.5, T.ENGINE_NEB)
    assert engine.tm_func(FORWARD) == pytest.approx(54.5)
    # One round trip for one distinct sequence; the second read is a cache hit.
    assert len(calls) == 1
    assert calls[0][1] == pytest.approx(T.CONDITIONS.neb_conc_um)
    assert calls[0][2] == T.NEB_PRODCODE
    # The estimate function is primer3, which is free, and is the same model
    # the structure bounds are read against -- not NEB, and not pydna's.
    assert engine.estimate_function(FORWARD) == pytest.approx(
        T.primer3_tm(FORWARD)
    )
    assert engine.estimate_function(FORWARD) != pytest.approx(
        T.local_tm(FORWARD).celsius
    )


def test_a_dead_endpoint_degrades_once_and_is_not_asked_again(monkeypatch):
    calls = []

    def dead(seq, conc, prodcode):
        calls.append(seq)
        raise requests.exceptions.ConnectionError("Could not connect to NEB API.")

    monkeypatch.setattr(T, "tm_neb", dead)
    engine = T.FootprintTm(offline=False)
    first = engine.tm(FORWARD)
    second = engine.tm(REVERSE)
    assert first.engine == T.ENGINE_LOCAL_TM
    assert second.engine == T.ENGINE_LOCAL_TM
    assert first.celsius == pytest.approx(T.local_tm(FORWARD).celsius)
    # tm_neb takes no timeout, so one hang must not become one hang per oligo.
    assert len(calls) == 1
    assert len(engine.failures) == 1
    assert "tm_neb" in engine.failures[0]


# --------------------------------------------------------------------------- #
# The one batched POST
# --------------------------------------------------------------------------- #


def test_offline_batch_makes_no_request_and_uses_the_local_ta_model():
    pcr = _pcr()
    results, failures = T.neb_pair_temperatures([pcr], offline=True)
    assert failures == ()
    temperatures = results[pcr]
    assert temperatures.ta_from_footprints.engine == T.ENGINE_LOCAL_TA
    assert temperatures.footprint_tm_forward.engine == T.ENGINE_LOCAL_TM
    assert temperatures.footprint_tm_forward.celsius == pytest.approx(
        T.local_tm(FOOTPRINT_F).celsius
    )


def test_one_post_carries_both_row_sets_and_ta_comes_from_the_footprint_rows(
    monkeypatch,
):
    posts = []

    def fake_post(url, json=None, timeout=None):
        posts.append((url, json, timeout))
        count = len(json["seqpairs"]) // 2
        oligo_rows = [{"tm1": 70.0, "tm2": 71.0, "ta": 99.0} for _ in range(count)]
        footprint_rows = [{"tm1": 56.0, "tm2": 57.0, "ta": 55.0} for _ in range(count)]
        return FakeResponse({"success": True, "data": oligo_rows + footprint_rows})

    monkeypatch.setattr(requests, "post", fake_post)
    first, second = _pcr("one"), _pcr("two", suffix="GGCC")
    results, failures = T.neb_pair_temperatures([first, second], offline=False)

    assert failures == ()
    assert len(posts) == 1, "N PCRs must cost one request, not N"
    url, payload, timeout = posts[0]
    assert url == T.NEB_BATCH_ENDPOINT
    assert timeout == T.NEB_TIMEOUT_S
    assert payload["conc"] == pytest.approx(T.CONDITIONS.neb_conc_um)
    assert payload["seqpairs"] == [
        [first.forward_oligo, first.reverse_oligo],
        [second.forward_oligo, second.reverse_oligo],
        [first.forward_footprint, first.reverse_footprint],
        [second.forward_footprint, second.reverse_footprint],
    ]

    temperatures = results[first]
    assert temperatures.oligo_tm_forward == T.Tm(70.0, T.ENGINE_NEB)
    assert temperatures.oligo_tm_reverse == T.Tm(71.0, T.ENGINE_NEB)
    assert temperatures.footprint_tm_forward == T.Tm(56.0, T.ENGINE_NEB)
    assert temperatures.footprint_tm_reverse == T.Tm(57.0, T.ENGINE_NEB)
    # 99.0 is NEB's recommendation for the tailed oligos; the footprints are
    # what anneal, so 55.0 is the only correct answer here.
    assert temperatures.ta_from_footprints == T.Tm(55.0, T.ENGINE_NEB)


def test_batch_deduplicates_identical_pcrs(monkeypatch):
    posts = []

    def fake_post(url, json=None, timeout=None):
        posts.append(json)
        count = len(json["seqpairs"]) // 2
        rows = [{"tm1": 60.0, "tm2": 61.0, "ta": 58.0} for _ in range(2 * count)]
        return FakeResponse({"success": True, "data": rows})

    monkeypatch.setattr(requests, "post", fake_post)
    pcr = _pcr()
    results, _ = T.neb_pair_temperatures([pcr, pcr, pcr], offline=False)
    assert len(results) == 1
    assert len(posts[0]["seqpairs"]) == 2


def test_an_unreachable_endpoint_degrades_to_local_values_and_reports(monkeypatch):
    def dead(url, json=None, timeout=None):
        raise requests.exceptions.ConnectionError("no route to host")

    monkeypatch.setattr(requests, "post", dead)
    pcr = _pcr()
    results, failures = T.neb_pair_temperatures([pcr], offline=False)
    assert len(failures) == 1
    assert T.NEB_BATCH_ENDPOINT in failures[0]
    assert results[pcr].ta_from_footprints.engine == T.ENGINE_LOCAL_TA
    assert results[pcr].oligo_tm_forward.engine == T.ENGINE_LOCAL_TM


@pytest.mark.parametrize(
    "footprint_row",
    [
        {"tm1": 56.0, "tm2": 57.0},  # no ta
        {"tm1": 56.0, "tm2": 57.0, "ta": None},
        {"tm1": 56.0, "tm2": 57.0, "ta": "n/a"},
        {"tm1": 56.0, "tm2": 57.0, "ta": float("nan")},
        {"error": "bad sequence"},
    ],
)
def test_a_row_without_a_usable_number_degrades_just_that_pcr(
    monkeypatch, footprint_row
):
    def fake_post(url, json=None, timeout=None):
        return FakeResponse(
            {
                "success": True,
                "data": [{"tm1": 70.0, "tm2": 71.0, "ta": 99.0}, footprint_row],
            }
        )

    monkeypatch.setattr(requests, "post", fake_post)
    pcr = _pcr("tadh1")
    results, failures = T.neb_pair_temperatures([pcr], offline=False)
    assert len(failures) == 1
    assert "tadh1" in failures[0]
    assert results[pcr].ta_from_footprints.engine == T.ENGINE_LOCAL_TA


def test_a_refused_batch_degrades(monkeypatch):
    monkeypatch.setattr(
        requests,
        "post",
        lambda url, json=None, timeout=None: FakeResponse(
            {"success": False, "error": "prodcode unknown"}
        ),
    )
    pcr = _pcr()
    results, failures = T.neb_pair_temperatures([pcr], offline=False)
    assert "prodcode unknown" in failures[0]
    assert results[pcr].footprint_tm_forward.engine == T.ENGINE_LOCAL_TM


def test_an_empty_batch_is_not_a_request():
    assert T.neb_pair_temperatures([], offline=False) == ({}, ())


def test_the_report_columns_keep_the_two_tm_quantities_apart(monkeypatch):
    def fake_post(url, json=None, timeout=None):
        return FakeResponse(
            {
                "success": True,
                "data": [
                    {"tm1": 70.0, "tm2": 71.0, "ta": 99.0},
                    {"tm1": 56.0, "tm2": 57.0, "ta": 55.0},
                ],
            }
        )

    monkeypatch.setattr(requests, "post", fake_post)
    pcr = _pcr()
    results, _ = T.neb_pair_temperatures([pcr], offline=False)
    row = results[pcr].to_dict()
    assert row == {
        "oligo_tm_forward": 70.0,
        "oligo_tm_reverse": 71.0,
        "oligo_tm_engine": T.ENGINE_NEB,
        "footprint_tm_forward": 56.0,
        "footprint_tm_reverse": 57.0,
        "footprint_tm_engine": T.ENGINE_NEB,
        "ta_from_footprints": 55.0,
        "ta_engine": T.ENGINE_NEB,
    }
    # No column holds "the Tm" without saying which, and the Ta column names
    # the pair it derives from.
    assert "tm" not in row
    assert "ta" not in row


def test_the_cheap_walk_is_primer3_and_the_reported_value_is_the_engines(monkeypatch):
    """primer3 chooses where to start looking; the engine says what it is.

    The walk used to ask the authoritative engine at every candidate length,
    which on a 67-gene library was 1,481 NEB requests.  primer3 now does that
    climb for free and the engine takes over at the end -- about 3.5 requests
    per footprint instead of ten -- and the Tm that reaches the window check
    and the tables is still the engine's, never primer3's.
    """
    from pydna.tm import tm_default

    import teemi.design.combinatorial.design as D

    asked = []

    def fake_tm_neb(seq, conc, prodcode):
        asked.append(seq)
        return tm_default(seq) - 1.8

    T.clear_tm_caches()
    monkeypatch.setattr(T, "tm_neb", fake_tm_neb)
    engine = T.FootprintTm(offline=False)
    template = FORWARD + "GCATGCATCGATCGGCATTACGGATCCAAGGTTACCGGAT"

    footprint, tm = D._shortest_footprint(template, template[:18], 60.0, engine)
    assert tm.engine == T.ENGINE_NEB, "the reported value is the engine's"
    assert tm.celsius >= 60.0
    assert T.primer3_tm(footprint) >= 60.0, "primer3 cleared it too"

    # The engine was asked only near the end, not from the seed up.
    assert len(asked) < len(footprint) - 18 + 1 or len(asked) <= 4, len(asked)
    assert all(len(seq) >= 18 for seq in asked)
    T.clear_tm_caches()
