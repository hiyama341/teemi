#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)

"""FLOW's design document as the input boundary.

The schema is the authority and :class:`FlowDocument` is a mirror of it, so the
mirror is cross-checked against it here rather than trusted: for every case
below, both have to agree on accept or reject.  Drift between them would
otherwise be invisible until a real document was refused or, worse, accepted.

The fixtures are the suite's own copies.  ``notebooks/specs/`` holds the ones
the notebooks iterate on and is gitignored, so a test reading those would break
on an unrelated edit and would not run at all in a fresh clone.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest
from pydantic import ValidationError

import teemi.design.combinatorial as cb
from teemi.design.combinatorial.flow import (
    TOPOLOGY_BY_METHOD,
    TRANSLATED_TYPES,
    FlowDocument,
    UndecidedPosition,
    UnknownPart,
    from_flow,
    plan_flow,
)

from tests.combinatorial_netguard import no_network  # noqa: F401

FIXTURES = Path(__file__).resolve().parents[1] / "tests" / "files_for_testing"
SCHEMA = FIXTURES / "flow_heterologous_expression.schema.json"
PLASMID = FIXTURES / "flow_example_invivo_plasmid.json"
INTEGRATION = FIXTURES / "flow_example_genomic_integration.json"
PARTS = FIXTURES / "flow_parts_used.json"


def document(path: Path) -> dict:
    return json.loads(path.read_text())


@pytest.fixture(scope="module")
def library() -> dict:
    saved = json.loads(PARTS.read_text())["parts"]
    return {part_id: row["sequence"] for part_id, row in saved.items()}


# --------------------------------------------------------------------------- #
# The mirror against the schema
# --------------------------------------------------------------------------- #


def _first_library_source(raw: dict) -> dict:
    """The first ``library`` source in the document, wherever it sits.

    The two documents differ in shape -- the plasmid's first position carries
    ``variants`` and the integration's carries ``source`` -- so a mutation that
    wants a library source has to go looking for one.
    """
    for slot in raw["construct"]:
        if isinstance(slot.get("source"), dict):
            return slot["source"]
        for variant in slot.get("variants") or ():
            if isinstance(variant.get("source"), dict):
                return variant["source"]
    raise AssertionError("no library source in this document")


#: ``(label, mutation)``.  ``None`` means the document unchanged.
CASES = [
    ("unchanged", None),
    ("both source and variants", lambda d: d["construct"][0].update(
        {"variants": [{"id": "v1", "source": {"type": "library", "part_id": "X"}}]}
    )),
    ("topology against method", lambda d: d.update({"topology": "circular"})),
    ("unknown top-level key", lambda d: d.update({"surprise": 1})),
    ("unknown part key", lambda d: d["construct"][0].update({"surprise": 1})),
    ("part type not in the vocabulary", lambda d: d["construct"][0].update(
        {"type": "not_a_type"}
    )),
    ("malformed accession", lambda d: d["organism"]["genome"].update(
        {"accession": "nope"}
    )),
    ("organism not in the list", lambda d: d.update({"organism": {"key": "e-coli"}})),
    ("part_id with a hyphen", lambda d: _first_library_source(d).update(
        {"part_id": "NOT-ALLOWED"}
    )),
    ("empty construct", lambda d: d.update({"construct": []})),
    ("empty variants list", lambda d: d["construct"][0].update(
        {"source": None, "variants": []}
    )),
    ("inline without a format", lambda d: d["construct"][0].update(
        {"source": {"type": "inline", "sequence": "ACGT"}}
    )),
    ("wrong schema version", lambda d: d.update({"flow_schema_version": "2.0"})),
    ("missing task", lambda d: d.pop("task")),
    ("design_count of zero", lambda d: d.update({"design_count": 0})),
    ("primer_conc_um of zero", lambda d: d.setdefault("assembly", {}).update(
        {"primer_conc_um": 0}
    )),
    ("oligo type not in the vocabulary", lambda d: d.update(
        {"oligos": [{"id": "o1", "type": "not_an_oligo"}]}
    )),
]


def mutated(path: Path, mutate) -> dict:
    raw = copy.deepcopy(document(path))
    if mutate is not None:
        mutate(raw)
    return raw


def accepted_by_mirror(raw: dict) -> bool:
    try:
        FlowDocument.model_validate(raw)
    except ValidationError:
        return False
    return True


def accepted_by_schema(raw: dict) -> bool:
    jsonschema = pytest.importorskip("jsonschema")
    validator = jsonschema.Draft202012Validator(json.loads(SCHEMA.read_text()))
    return not list(validator.iter_errors(raw))


@pytest.mark.parametrize("path", [PLASMID, INTEGRATION], ids=["plasmid", "integration"])
@pytest.mark.parametrize("label,mutate", CASES, ids=[c[0] for c in CASES])
def test_the_mirror_and_the_schema_agree(path, label, mutate):
    """Neither may accept what the other rejects."""
    raw = mutated(path, mutate)
    by_schema = accepted_by_schema(raw)
    by_mirror = accepted_by_mirror(raw)
    assert by_mirror == by_schema, (
        f"{label}: schema {'accepts' if by_schema else 'rejects'} but the "
        f"mirror {'accepts' if by_mirror else 'rejects'}"
    )


def test_both_real_documents_parse():
    """The point of the whole exercise: FLOW's own files, unmodified."""
    for path, positions, oligos in ((PLASMID, 6, 1), (INTEGRATION, 8, 0)):
        doc = FlowDocument.model_validate(document(path))
        assert len(doc.positions) == positions
        assert len(doc.oligos) == oligos
        assert doc.resolved_topology == TOPOLOGY_BY_METHOD[doc.method]


def test_the_wire_name_is_construct_and_the_attribute_is_positions():
    """``construct`` shadows a pydantic member, so the attribute is renamed.

    Both names have to validate, or a document could only be read by going
    through the alias and the mirror would stop being a mirror.
    """
    raw = document(INTEGRATION)
    assert "construct" in raw and "positions" not in raw
    by_wire = FlowDocument.model_validate(raw)
    by_name = FlowDocument.model_validate(
        {**{k: v for k, v in raw.items() if k != "construct"},
         "positions": raw["construct"]}
    )
    assert by_wire == by_name


def test_topology_is_derived_when_the_document_leaves_it_out():
    """It is not a required field, so it is resolved rather than guessed later."""
    raw = document(INTEGRATION)
    raw.pop("topology")
    assert FlowDocument.model_validate(raw).resolved_topology == "linear"
    raw = document(PLASMID)
    raw.pop("topology")
    assert FlowDocument.model_validate(raw).resolved_topology == "circular"


# --------------------------------------------------------------------------- #
# Document plus policy -> spec
# --------------------------------------------------------------------------- #


def test_the_document_supplies_the_design_and_policy_the_bounds(library):
    """The division of labour the schema itself sets out.

    It carries no primer-design policy, no reaction conditions and no objective
    weights; those belong to the backend's spec, which defaults every one.  So
    a document with no policy at all still makes a complete spec.
    """
    gene = library["FP_MRFP1"]
    spec, _ = from_flow(
        document(INTEGRATION), library=library, fill={"pos4": gene}
    )

    # From the document.
    assert spec.design.topology == "linear"
    assert spec.primer_design.junction.homology_length == 40
    assert spec.primer_design.target_tm == 65
    assert spec.primer_design.footprint_len_min == 15
    assert [p.id for p in spec.positions][:2] == ["pos1", "pos2"]

    # Defaulted by the backend, because FLOW never sends them.
    assert spec.constraints.homodimer_dg_min == -6.0
    assert spec.constraints.primer_len_max == 60
    assert spec.objective.price_per_base > 0

    # And overridable as policy.
    tightened, _ = from_flow(
        document(INTEGRATION),
        library=library,
        policy={"constraints": {"primer_len_max": 90}},
        fill={"pos4": gene},
    )
    assert tightened.constraints.primer_len_max == 90


def test_tune_overhangs_comes_across_under_its_own_name(library):
    gene = library["FP_MRFP1"]
    plasmid, _ = from_flow(document(PLASMID), library=library, fill={"pos3": gene})
    assert plasmid.primer_design.junction.tune_overhangs is True

    integration, _ = from_flow(
        document(INTEGRATION), library=library, fill={"pos4": gene}
    )
    assert integration.primer_design.junction.tune_overhangs is False


def test_every_slot_flag_survives_the_translation(library):
    gene = library["FP_MRFP1"]
    spec, _ = from_flow(document(INTEGRATION), library=library, fill={"pos4": gene})
    by_id = {p.id: p for p in spec.positions}
    assert by_id["pos3"].in_tail is True
    assert by_id["pos5"].drop_start is True and by_id["pos5"].drop_stop is True
    assert by_id["pos6"].drop_start is True
    assert by_id["pos7"].no_tail_rev is True
    assert by_id["pos8"].reverse is True
    assert by_id["pos1"].reverse is False


def test_a_position_with_no_choice_yet_is_named_not_guessed(library):
    """FLOW allows it; there is simply nothing to design against it."""
    with pytest.raises(UndecidedPosition) as caught:
        from_flow(document(INTEGRATION), library=library)
    assert caught.value.position_id == "pos4"
    assert "fill=" in str(caught.value)


def test_a_retired_part_id_fails_loudly(library):
    """FLOW's own schema warns that its lookup returns undefined silently."""
    raw = document(INTEGRATION)
    raw["construct"][0]["source"]["part_id"] = "PRO_RETIRED_LONG_AGO"
    with pytest.raises(UnknownPart) as caught:
        from_flow(raw, library=library, fill={"pos4": library["FP_MRFP1"]})
    assert caught.value.part_id == "PRO_RETIRED_LONG_AGO"
    assert "library_version" in str(caught.value)


def test_a_drop_flag_on_an_untranslated_type_is_reported(library):
    """A note, not an error: it is not wrong to send, but nothing comes of it."""
    raw = document(INTEGRATION)
    assert raw["construct"][1]["type"] == "promoter"
    assert raw["construct"][1]["type"] not in TRANSLATED_TYPES
    raw["construct"][1]["drop_start"] = True
    _, notes = from_flow(raw, library=library, fill={"pos4": library["FP_MRFP1"]})
    assert any("not" in n and "translated" in n for n in notes), notes


def test_the_oligos_are_reported_as_not_assembled(library):
    """They sit alongside the construct, so they are not in the plan."""
    _, notes = from_flow(
        document(PLASMID), library=library, fill={"pos3": library["FP_MRFP1"]}
    )
    assert any("not assembled" in n and "o1" in n for n in notes), notes


def test_an_inline_source_is_parsed_from_the_file_text(library):
    """FLOW sends the uploaded file's text, not bare DNA."""
    raw = document(PLASMID)
    inline = raw["construct"][2]["source"]
    assert inline["type"] == "inline" and inline["format"] == "genbank"
    spec, _ = from_flow(raw, library=library)
    sequence = spec.positions[2].parts[0].source.sequence
    assert sequence == "ATGAGCAAAGGTGAAGAATTG"


# --------------------------------------------------------------------------- #
# End to end
# --------------------------------------------------------------------------- #


def test_plan_flow_designs_a_real_document(library):
    result, notes = plan_flow(
        document(INTEGRATION),
        library=library,
        policy={"constraints": {"primer_len_max": 90}},
        fill={"pos4": library["FP_MRFP1"]},
        offline=True,
        include_designs=True,
    )
    assert result.spec is not None
    assert len(result.variants) == 1
    variant = result.variants[0]
    assert variant.verification.verdict == "verified"

    # pos3 is in_tail at 6 bp, so it is in the construct but not amplified.
    assert "KOZ_SC_TIS8" not in variant.fragment_ids
    assert "KOZ_SC_TIS8" in variant.id
    assert len(result.pcrs) == len(variant.fragment_ids)

    # The no_tail_rev position's junction cannot go through pydna.
    built_by = {j.designed_by for j in result.junctions}
    assert built_by == {"pydna", "forced"}

    json.dumps(result.to_dict())
    assert result.designs and len(result.designs[0].construct) == variant.length_bp
