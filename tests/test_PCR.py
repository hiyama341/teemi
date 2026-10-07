#!/usr/bin/env python

# Test PCR module

from pydna.dseqrecord import Dseqrecord
from pydna.amplify import pcr
from pydna.primer import Primer
from Bio.SeqRecord import SeqRecord


# Importing the module we are  testing
from teemi.build.PCR import *

def test_amplicon_by_name(): 
    
    # initialize
    middle = 'a'*2000
    template = Dseqrecord("tacactcaccgtctatcattatcagcgacgaagcgagcgcgaccgcgagcgcgagcgca"+middle+"caggagcgagacacggcgacgcgagcgagcgagcgatactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatcattatc")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon = pcr(p1, p2, template)
    amplicon.name = 'AMPICON_FOR_TESTING_amplicon_byname_function'

    amplicon_list = [amplicon]
    my_amplicon = amplicon_by_name('AMPICON_FOR_TESTING_amplicon_byname_function', amplicon_list)

    assert my_amplicon.name == 'AMPICON_FOR_TESTING_amplicon_byname_function'
    assert len(my_amplicon) == 2123
    
def test_calculate_processing_speed(): 
    
    # initialize
    template = Dseqrecord("tacactcaccgtctatcattatctactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatcattatc")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon = pcr(p1, p2, template)
    
    # the tests
    amplicon.annotations['polymerase'] = "OneTaq Hot Start"
    calculate_processing_speed(amplicon)
    assert amplicon.annotations['proc_speed'] == 60
    
    amplicon.annotations['polymerase'] = "Q5 Hot Start"
    calculate_processing_speed(amplicon)
    assert amplicon.annotations['proc_speed'] == 30
    

    amplicon.annotations['polymerase'] = "Phusion"
    calculate_processing_speed(amplicon)
    assert amplicon.annotations['proc_speed'] == 30
    

def test_calculate_elongation_time():
    
    # initialize
    template = Dseqrecord("tacactcaccgtctatcattatctactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatcattatc")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon = pcr(p1, p2, template)
    amplicon.annotations['polymerase'] = "OneTaq Hot Start"
    calculate_processing_speed(amplicon)
    #the test
    calculate_elongation_time(amplicon)
    
    assert amplicon.annotations['elongation_time'] == 4
    
    
    # initialize 2
    middle = 'a'*2000
    template = Dseqrecord("tacactcaccgtctatcattatcagcgacgaagcgagcgcgaccgcgagcgcgagcgca"+middle+"caggagcgagacacggcgacgcgagcgagcgagcgatactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatcattatc")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon = pcr(p1, p2, template)
    amplicon.annotations['polymerase'] = "OneTaq Hot Start"
    calculate_processing_speed(amplicon)
    calculate_elongation_time(amplicon)

    #tests
    assert amplicon.annotations['elongation_time'] == 128
    
    
    # initialize 2
    middle = 'a'*3000
    template = Dseqrecord("tacactcaccgtctatcattatcagcgacgaagcgagcgcgaccgcgagcgcgagcgca"+middle+"caggagcgagacacggcgacgcgagcgagcgagcgatactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatcattatc")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon = pcr(p1, p2, template)
    amplicon.annotations['polymerase'] = "OneTaq Hot Start"
    calculate_processing_speed(amplicon)
    calculate_elongation_time(amplicon)

    #test
    assert amplicon.annotations['elongation_time'] == 188

    

def test_PCR_program(): 
    # intialize an amplicon 
    middle = 'a'*2000
    template = Dseqrecord("tacactcaccgtctatcattatcagcgacgaagcgagcgcgaccgcgagcgcgagcgca"+middle+"caggagcgagacacggcgacgcgagcgagcgagcgatactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatcattatc")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon = pcr(p1, p2, template)
    amplicon.name = 'AMPICON_FOR_TESTING_PCR_program_function'
    amplicon.annotations['polymerase'] = "OneTaq Hot Start"
    calculate_processing_speed(amplicon)

    # initialize a string object
    program = Q5_NEB_PCR_program(amplicon)

    assert program[1:5] == '98°C'
    assert program[35:39] == '61.0'
    assert program[61:65] == '72°C'


def test_grouper():
    elong_times = [60,60, 46, 60, 45, 30, 200, 100]
    elong_times.sort()
    elong_time_max_diff = 10
    groups = dict(enumerate(grouper(elong_times,elong_time_max_diff), 1))

    assert groups == {1: [30], 2: [45, 46], 3: [60, 60, 60], 4: [100], 5: [200]}



def test_calculate_required_thermal_cyclers():
    # amplicon1
    middle = 'c'*100
    template = Dseqrecord("tacactcaccgtctatcattatcagcgacgaagcgagcgcgaccgcgagcgcgagcgca"+middle+"caggagcgagacacggcgacgcgagcgagcgagcgatactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatca")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon1 = pcr(p1, p2, template)
    amplicon1.name = 'AMPLICON1'
    amplicon1.annotations['polymerase'] = "OneTaq Hot Start"
    calculate_processing_speed(amplicon1)
    calculate_elongation_time(amplicon1)

    # amplicon2
    middle = 'a'*200
    template = Dseqrecord("tacactcaccgtctatcattatcagcgacgaagcgagcgcgaccgcgagcgcgagcgca"+middle+"caggagcgagacacggcgacgcgagcgagcgagcgatactatcgactgtatcatctgatagcac")
    p1 = Primer("tacactcaccgtctatca")
    p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    amplicon2 = pcr(p1, p2, template)
    amplicon2.name = 'AMPLICON2'
    amplicon2.annotations['polymerase'] = "OneTaq Hot Start"
    calculate_processing_speed(amplicon2)
    calculate_elongation_time(amplicon2)

    # adding extra annotations
    amplicon1_program = Q5_NEB_PCR_program(amplicon1)
    amplicon2_program = Q5_NEB_PCR_program(amplicon2)

    amplicons = [amplicon1, amplicon2]
    # running the function
    thermal_cyclers = calculate_required_thermal_cyclers(amplicons, polymerase='Q5 Hot Start') # pol

    assert thermal_cyclers.iloc[0]['tas'] == 59
    assert thermal_cyclers.iloc[0]['elong_times'] == 20
    assert thermal_cyclers.iloc[0]['amplicons'] == 'AMPLICON1, AMPLICON2'


def test_pcr_locations():
    from pydna.dseqrecord import Dseqrecord
    from pydna.amplify import pcr


    dna = Dseqrecord('ATGATATATGGCTCGACTGCAGGGGGATTTTTCCGGATCGCGGTCGATGACTGATACTACTACGACTACTAG')
    primer1 = Dseqrecord('ATGATATATGGCTCGAC')
    primer2 = Dseqrecord('TACTACGACTACTAG').reverse_complement()

    #names 
    dna.name = 'dna'
    primer1.name = 'primer1' 
    primer2.name  = 'primer2' 

    #annott
    dna.annotations['batches'] = [{'location':'Freezer1', 'concentration':123}]
    primer1.annotations['batches']  = [{'location':'Freezer2','concentration':274}]
    primer2.annotations['batches']   = [{'location':'Freezer3','concentration':124}]

    # make a pcr_prod
    gRNA1_pcr_prod = pcr(primer1,primer2, dna)

    # run the function
    pcr_locations_df = pcr_locations([gRNA1_pcr_prod])

    assert pcr_locations_df.iloc[0]['location'] == 'Freezer1'
    assert pcr_locations_df.iloc[0]['name'] == '72bp_PCR_prod'
    assert pcr_locations_df.iloc[0]['template'] == 'Freezer1'
    assert pcr_locations_df.iloc[0]['fw'] == 'Freezer2'
    assert pcr_locations_df.iloc[0]['rv'] == 'Freezer3'


def test_nanophotometer_concentrations(): 
    list_of_conc = nanophotometer_concentrations(path = '../teemi/tests/files_for_testing/2021-03-29_G8H_CPR_library_part_concentrations.tsv')

    assert list_of_conc[0] ==142.8
    assert list_of_conc[1] ==134.5

    assert list_of_conc[-2] ==17.5
    assert list_of_conc[-1] ==39.9


def test_calculate_volumes(): 

    calculate_volumes_df = calculate_volumes(vol_p_reac = 20, 
            no_of_reactions = 3,
            standard_reagents = ["Template", "Primer 1", "Primer 2", "H20", "Pol"],
            standard_volumes = [1, 2.5, 2.5, 19, 25])
    
    # template
    assert calculate_volumes_df.iloc[0]['vol_p_reac'] == 0.4
    assert round(calculate_volumes_df.iloc[0]['vol_p_3_reac'],2) == 1.2

    # h2o
    assert calculate_volumes_df.iloc[3]['vol_p_reac'] == 7.6
    assert round(calculate_volumes_df.iloc[3]['vol_p_3_reac'],2) == 22.8

    # total
    assert calculate_volumes_df.iloc[5]['vol_p_reac'] == 20.0
    assert round(calculate_volumes_df.iloc[5]['vol_p_3_reac'],2) == 60.0


#########################################################################
# Offline tests: the NEB Tm API is always mocked below (no network access)
#########################################################################
import json as _json

import pandas as pd
import pytest
import requests
from Bio.SeqUtils import MeltingTemp as _BioMeltingTemp

import teemi.build.PCR as PCR_module


class _FakeNEBResponse:
    def __init__(self, payload=None, http_error=None, json_error=None):
        self._payload = payload
        self._http_error = http_error
        self._json_error = json_error

    def raise_for_status(self):
        if self._http_error is not None:
            raise self._http_error

    def json(self):
        if self._json_error is not None:
            raise self._json_error
        return self._payload


def _make_amplicon(template_seq, fw_seq, rv_seq, name=None, template_name=None):
    template = Dseqrecord(template_seq)
    if template_name is not None:
        template.name = template_name
    fw = Primer(fw_seq, id="fw_" + (name or "x"), name="fw_" + (name or "x"))
    rv = Primer(rv_seq, id="rv_" + (name or "x"), name="rv_" + (name or "x"))
    amplicon = pcr(fw, rv, template)
    if name is not None:
        amplicon.name = name
    return amplicon


def test_post_neb_tm_api_posts_json_payload_and_returns_decoded_json(monkeypatch):
    calls = []

    def fake_post(url, data=None, headers=None, timeout=None):
        calls.append({"url": url, "data": data, "headers": headers, "timeout": timeout})
        return _FakeNEBResponse(payload={"success": True, "data": [{"tm1": 60.1}]})

    monkeypatch.setattr(PCR_module.requests, "post", fake_post)
    payload = {"seqpairs": [["ACGTACGTACGT"]], "conc": 0.5, "prodcode": "q5-0"}

    result = PCR_module._post_neb_tm_api(payload)

    assert result == {"success": True, "data": [{"tm1": 60.1}]}
    assert len(calls) == 1
    assert calls[0]["url"] == "https://tmapi.neb.com/tm/batch"
    assert calls[0]["headers"] == {"content-type": "application/json"}
    assert calls[0]["timeout"] == 15
    assert _json.loads(calls[0]["data"]) == payload


@pytest.mark.parametrize(
    "fake_post",
    [
        # network failure
        lambda *a, **k: (_ for _ in ()).throw(requests.ConnectionError("offline")),
        # HTTP error status
        lambda *a, **k: _FakeNEBResponse(http_error=requests.HTTPError("500")),
        # response body is not valid JSON
        lambda *a, **k: _FakeNEBResponse(json_error=ValueError("no json")),
    ],
    ids=["connection_error", "http_error", "invalid_json"],
)
def test_post_neb_tm_api_returns_none_when_api_unavailable(monkeypatch, fake_post):
    monkeypatch.setattr(PCR_module.requests, "post", fake_post)

    assert PCR_module._post_neb_tm_api({"seqpairs": [["ACGT"]]}) is None


def test_primer_tm_neb_uses_api_tm(monkeypatch):
    payloads = []

    def fake_api(payload):
        payloads.append(payload)
        return {"success": True, "data": [{"tm1": 63.4, "tm2": None, "ta": None}]}

    monkeypatch.setattr(PCR_module, "_post_neb_tm_api", fake_api)

    assert primer_tm_neb("ACGTGCTAGCTAGCTAGCAT", conc=0.25, prodcode="phusion-0") == 63.4
    assert payloads == [
        {"seqpairs": [["ACGTGCTAGCTAGCTAGCAT"]], "conc": 0.25, "prodcode": "phusion-0"}
    ]


@pytest.mark.parametrize("api_response", [None, {"success": False, "error": ["bad"]}])
def test_primer_tm_neb_falls_back_to_local_nearest_neighbour_tm(monkeypatch, api_response):
    monkeypatch.setattr(PCR_module, "_post_neb_tm_api", lambda payload: api_response)
    primer = "tacactcaccgtctatcattatc"

    tm = primer_tm_neb(primer)

    assert tm == round(_BioMeltingTemp.Tm_NN(primer), 1)
    assert isinstance(tm, float)


def test_fallback_primer_tm_uses_wallace_rule_without_biopython(monkeypatch):
    monkeypatch.setattr(PCR_module, "_MeltingTemp", None)

    # Wallace rule: Tm = 2 * (A + T) + 4 * (G + C)
    assert PCR_module._fallback_primer_tm("AATTGGCC") == 24.0
    assert PCR_module._fallback_primer_tm("aattgc") == 16.0
    # non-str sequence objects are converted with str()
    assert PCR_module._fallback_primer_tm(Primer("GGGGAAAA").seq) == 24.0


def test_primer_ta_neb_uses_api_ta(monkeypatch):
    payloads = []

    def fake_api(payload):
        payloads.append(payload)
        return {"success": True, "data": [{"tm1": 62.0, "tm2": 64.0, "ta": 63}]}

    monkeypatch.setattr(PCR_module, "_post_neb_tm_api", fake_api)

    assert primer_ta_neb("AAAACCCCGGGG", "TTTTGGGGCCCC") == 63
    assert payloads == [
        {"seqpairs": [["AAAACCCCGGGG", "TTTTGGGGCCCC"]], "conc": 0.5, "prodcode": "q5-0"}
    ]


def test_primer_ta_neb_falls_back_to_lowest_primer_tm(monkeypatch):
    monkeypatch.setattr(PCR_module, "_post_neb_tm_api", lambda payload: None)
    monkeypatch.setattr(PCR_module, "_MeltingTemp", None)

    # Wallace Tm: AAAAAAAAAA -> 20, GGGGGGGGGG -> 40, so Ta = min = 20
    assert primer_ta_neb("AAAAAAAAAA", "GGGGGGGGGG") == 20
    # rounding of the minimum Tm
    monkeypatch.setattr(PCR_module, "_fallback_primer_tm", lambda p: {"P1": 57.6, "P2": 61.2}[p])
    assert primer_ta_neb("P1", "P2") == 58


def test_Q5_NEB_PCR_program_offline(monkeypatch):
    def fake_api(payload):
        pair = payload["seqpairs"][0]
        if len(pair) == 2:
            return {"success": True, "data": [{"ta": 62}]}
        tms = {"tacactcaccgtctatcattatc": 59.5, "gtgctatcagatgatacagtcg": 61.27}
        return {"success": True, "data": [{"tm1": tms[pair[0]]}]}

    monkeypatch.setattr(PCR_module, "_post_neb_tm_api", fake_api)

    middle = "a" * 2000
    amplicon = _make_amplicon(
        "tacactcaccgtctatcattatcagcgacgaagcgagcgcgaccgcgagcgcgagcgca"
        + middle
        + "caggagcgagacacggcgacgcgagcgagcgagcgatactatcgactgtatcatctgatagcac",
        "tacactcaccgtctatcattatc",
        "gtgctatcagatgatacagtcg",
        name="AMP",
    )
    amplicon.annotations["polymerase"] = "Q5 Hot Start"
    calculate_processing_speed(amplicon)

    program = Q5_NEB_PCR_program(amplicon)
    lines = program.splitlines()

    # 30 s/kb * 2123 bp = 63.69 s -> 64 s -> 1 min 4 s
    assert amplicon.annotations["elongation_time"] == 64
    assert amplicon.annotations["ta Q5 Hot Start"] == 62
    assert amplicon.forward_primer.annotations["tm Q5 Hot Start"] == 59.5
    assert amplicon.reverse_primer.annotations["tm Q5 Hot Start"] == 61.27
    assert lines[0].endswith("|tmf:59.5")
    assert lines[1].endswith("|tmr:61.3")
    assert "\\ 62.0°C" in lines[2]
    assert lines[2].endswith("|30s/kb")
    assert " 1: 4|2min|" in lines[3]
    assert lines[4].endswith("|2123bp")


def test_calculate_processing_speed_keeps_preset_proc_speed(capsys):
    amplicon = _make_amplicon(
        "tacactcaccgtctatcattatctactatcgactgtatcatctgatagcac",
        "tacactcaccgtctatcattatc",
        "gtgctatcagatgatacagtcg",
    )
    amplicon.annotations["polymerase"] = "Q5 Hot Start"
    amplicon.annotations["proc_speed"] = 45
    amplicon.forward_primer.annotations["proc_speed"] = 45

    result = calculate_processing_speed(amplicon)

    assert result is amplicon
    assert amplicon.annotations["proc_speed"] == 45
    assert "proc_speed already set" in capsys.readouterr().out


def test_calculate_elongation_time_keeps_preset_elongation_time(capsys):
    amplicon = _make_amplicon(
        "tacactcaccgtctatcattatctactatcgactgtatcatctgatagcac",
        "tacactcaccgtctatcattatc",
        "gtgctatcagatgatacagtcg",
    )
    amplicon.annotations["proc_speed"] = 30
    amplicon.annotations["elongation_time"] = 99
    amplicon.forward_primer.annotations["elongation_time"] = 99

    result = calculate_elongation_time(amplicon)

    assert result is amplicon
    assert amplicon.annotations["elongation_time"] == 99
    assert "elongation_time already set" in capsys.readouterr().out


def test_pcr_locations_falls_back_to_amplicon_batches_and_empty(capsys):
    template_seq = "ATGATATATGGCTCGACTGCAGGGGGATTTTTCCGGATCGCGGTCGATGACTGATACTACTACGACTACTAG"

    # template has an empty batch list, the amplicon itself carries the location
    amp1 = _make_amplicon(template_seq, "ATGATATATGGCTCGAC", "CTAGTAGTCGTAGTA", name="AMP1")
    amp1.template.annotations["batches"] = []
    amp1.annotations["batches"] = [{"location": "PCRplate_A1"}]
    amp1.forward_primer.annotations["batches"] = [{"location": "Primers_B1"}]
    amp1.reverse_primer.annotations["batches"] = [{"location": "Primers_B2"}]

    # nothing has a location
    amp2 = _make_amplicon(template_seq, "ATGATATATGGCTCGAC", "CTAGTAGTCGTAGTA", name="AMP2")
    amp2.forward_primer.annotations["batches"] = []

    df = pcr_locations([amp1, amp2])

    assert list(df.columns) == ["location", "name", "template", "fw", "rv"]
    assert df.iloc[0].tolist() == ["PCRplate_A1", "AMP1", "PCRplate_A1", "Primers_B1", "Primers_B2"]
    assert df.iloc[1].tolist() == ["Empty", "AMP2", "Empty", "Empty", "Empty"]

    out = capsys.readouterr().out
    assert "No batches were found for AMP2. Please check the object." in out
    assert "AMP2: Foward primer location was not found" in out
    assert "AMP2: Reverse primer location was not found" in out
    assert "AMP1" not in out


def _located_amplicons():
    template_seq = "ATGATATATGGCTCGACTGCAGGGGGATTTTTCCGGATCGCGGTCGATGACTGATACTACTACGACTACTAG"
    amplicons = []
    for i, name in enumerate(["AMP_A", "AMP_B"], start=1):
        amp = _make_amplicon(template_seq, "ATGATATATGGCTCGAC", "CTAGTAGTCGTAGTA", name=name)
        amp.annotations["batches"] = [{"location": f"pcr_{i}"}]
        amp.annotations["template_name"] = f"template_{i}"
        amp.template.annotations["batches"] = [{"location": f"tmpl_{i}"}]
        amp.forward_primer.annotations["batches"] = [{"location": f"fw_loc_{i}"}]
        amp.reverse_primer.annotations["batches"] = [{"location": f"rv_loc_{i}"}]
        amplicons.append(amp)
    return amplicons


def test_set_plate_locations():
    amplicons = _located_amplicons()

    df = set_plate_locations(amplicons)

    assert df.index.name == "name"
    assert list(df.index) == ["AMP_A", "AMP_B"]
    assert list(df.columns) == [
        "location",
        "template_name",
        "template_location",
        "fw_name",
        "fw_location",
        "rv_name",
        "rv_location",
    ]
    assert df.loc["AMP_A"].tolist() == [
        "pcr_1", "template_1", "tmpl_1", "fw_AMP_A", "fw_loc_1", "rv_AMP_A", "rv_loc_1"
    ]
    assert df.loc["AMP_B"].tolist() == [
        "pcr_2", "template_2", "tmpl_2", "fw_AMP_B", "fw_loc_2", "rv_AMP_B", "rv_loc_2"
    ]


def test_update_amplicon_annotations():
    amplicons = _located_amplicons()

    result = update_amplicon_annotations(
        ["AMP_B", "AMP_A"], amplicons, ["H1", "H2"], [12.5, 30.0], [40, 50]
    )

    assert result is None
    assert amplicons[0].annotations["batches"][0] == {
        "location": "H2", "concentration": 30.0, "volume": 50
    }
    assert amplicons[1].annotations["batches"][0] == {
        "location": "H1", "concentration": 12.5, "volume": 40
    }


def test_get_amplicons_by_row_and_column():
    template_seq = "ATGATATATGGCTCGACTGCAGGGGGATTTTTCCGGATCGCGGTCGATGACTGATACTACTACGACTACTAG"
    amps = [
        _make_amplicon(template_seq, "ATGATATATGGCTCGAC", "CTAGTAGTCGTAGTA", name=n)
        for n in ["A1_amp", "A2_amp", "B1_amp"]
    ]
    amplicon_df = pd.DataFrame(
        {
            "name": ["A1_amp", "A2_amp", "B1_amp", "not_in_list"],
            "prow": ["A", "A", "B", "A"],
            "pcol": [1, 2, 1, 1],
        }
    )

    row_a = get_amplicons_by_row("A", amplicon_df, amps)
    row_c = get_amplicons_by_row("C", amplicon_df, amps)
    col_1 = get_amplicons_by_column(1, amplicon_df, amps)

    assert [len(group) for group in row_a] == [1, 1]
    assert row_a[0][0] is amps[0] and row_a[1][0] is amps[1]
    assert row_c == []
    assert [len(group) for group in col_1] == [1, 1]
    assert col_1[0][0] is amps[0] and col_1[1][0] is amps[2]
