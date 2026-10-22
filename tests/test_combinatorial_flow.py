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


# --------------------------------------------------------------------------- #
# A combinatorial position whose variants carry no sequences
# --------------------------------------------------------------------------- #


CDS_LIBRARY = FIXTURES / "flow_example_cds_library.json"

#: Four real coding sequences, standing in for four genes of the user's own.
GENE_PARTS = ("FP_GP_EGFP", "FP_MCITRINE", "FP_MRFP1", "MRK_GP_KANMX")


@pytest.fixture()
def genes(library) -> dict:
    return {
        f"gene/gene-{n}": library[part]
        for n, part in enumerate(GENE_PARTS, start=1)
    }


def test_fill_can_target_one_variant_rather_than_a_whole_position(library, genes):
    """A library of four genes needs four sequences at one position.

    Keyed by position alone, ``fill`` can only replace the whole slot with a
    single part, which would collapse the library to one member.  The key is
    the ``position/variant`` form :class:`UndecidedPosition` names.
    """
    spec, notes = from_flow(sourceless_genes(), library=library, fill=genes)

    gene_slot = next(p for p in spec.positions if p.id == "gene")
    assert len(gene_slot.parts) == 4
    assert [p.id for p in gene_slot.parts] == [
        f"gene.gene-{n}" for n in (1, 2, 3, 4)
    ]
    assert {p.source.sequence for p in gene_slot.parts} == {
        library[part] for part in GENE_PARTS
    }
    assert spec.n_variants == 4 == document(CDS_LIBRARY)["design_count"]
    assert sum("stood in" in n for n in notes) == 4


def sourceless_genes() -> dict:
    """The CDS-library document with its four sequences taken back out.

    Built here rather than kept as a second fixture: it is the same document in
    the state FLOW saves it in before any gene has been chosen, and deriving it
    cannot drift from the real one.
    """
    raw = document(CDS_LIBRARY)
    for slot in raw["construct"]:
        if slot["id"] == "gene":
            for variant in slot["variants"]:
                variant.pop("source", None)
    return raw


def test_a_variant_with_no_sequence_and_no_fill_is_named(library):
    """Named by position AND variant, so it says which of the four is missing."""
    partial = {f"gene/gene-{n}": library["FP_MRFP1"] for n in (1, 2, 3)}
    with pytest.raises(UndecidedPosition) as caught:
        from_flow(sourceless_genes(), library=library, fill=partial)
    assert caught.value.position_id == "gene/gene-4"

    # And with nothing supplied at all, the first one is named.
    with pytest.raises(UndecidedPosition) as caught:
        from_flow(sourceless_genes(), library=library)
    assert caught.value.position_id == "gene/gene-1"


def test_a_fill_key_that_matches_nothing_is_refused(library, genes):
    """A typo would otherwise leave the real slot undecided and say nothing."""
    with pytest.raises(KeyError, match="not a position"):
        from_flow(
            document(CDS_LIBRARY),
            library=library,
            fill={**genes, "gene/gene-5": library["FP_MRFP1"]},
        )
    with pytest.raises(KeyError, match="not a position"):
        from_flow(
            document(CDS_LIBRARY),
            library=library,
            fill={**genes, "promotor": library["FP_MRFP1"]},
        )


def test_no_tail_on_the_backbone_shares_one_primer_across_the_library(
    library, genes
):
    """The point of this document, and the thing worth asserting about it.

    ``no_tail_rev`` on the promoter and ``no_tail_fwd`` on the terminator leave
    those two primers as bare footprints with no gene-specific tail, so one
    molecule serves every member and adding a fifth gene does not change them.
    """
    # The document supplies its own four sequences, so nothing is stood in.
    result, notes = plan_flow(
        document(CDS_LIBRARY),
        library=library,
        policy={"constraints": {"primer_len_max": 90}},
        offline=True,
    )
    assert notes == ()
    assert len(result.variants) == 4
    assert all(v.verification.verdict == "verified" for v in result.variants)

    shared = [p for p in result.primers if len(p.uses) == 4]
    assert len(shared) == 2, [p.id for p in result.primers]
    assert {p.part_id for p in shared} == {"PRO_ANID_TEF1", "TER_ANID_TEF1"}
    for primer in shared:
        assert primer.sequence == primer.footprint, "a bare footprint, no tail"

    # Fewer oligos than primer slots, and the price counts each molecule once.
    slots = sum(len(p.uses) for p in result.primers)
    assert slots > len(result.primers)

    # Every CDS seam is pinned by a flag, so pydna is not asked for those.
    by_who = {j.designed_by for j in result.junctions}
    assert by_who == {"forced", "pydna"}
    forced = [j for j in result.junctions if j.designed_by == "forced"]
    assert all(len(j.profile) == 1 for j in forced)
    assert all(j.unavoidable == () for j in forced), (
        "one candidate proves nothing, so a pinned junction may not downgrade"
    )


def test_the_polymerase_the_document_names_reaches_the_neb_request(
    library, monkeypatch
):
    """FLOW sends it; the schema says choosing the Tm model is the backend's job.

    It is not cosmetic: pydna's own note is that Q5 takes the lower of the two
    Tms and adds 1 C while Phusion adds 3, so a document asking for Phusion and
    getting Q5's annealing temperature would be 2 C out.  Both NEB paths are
    checked -- the per-sequence walk that designs the footprints and the batch
    that returns the authoritative pair -- because the product code has to
    reach both or the two would disagree.
    """
    import requests

    import teemi.design.combinatorial.tm as tm_module

    spec, _ = from_flow(document(CDS_LIBRARY), library=library)
    assert spec.primer_design.polymerase == "phusion-0"

    walked, batched = [], []

    def fake_tm_neb(seq, conc=0.5, prodcode="q5-0"):
        walked.append(prodcode)
        return 60.0

    class FakeReply:
        status_code = 200

        def __init__(self, rows):
            self._rows = rows

        def raise_for_status(self):
            return None

        def json(self):
            return {"data": self._rows}

    def fake_post(url, json=None, **kwargs):
        batched.append(json.get("prodcode"))
        return FakeReply(
            [{"tm1": 60.0, "tm2": 60.0, "ta": 63.0} for _ in json["seqpairs"]]
        )

    tm_module.clear_tm_caches()
    monkeypatch.setattr(tm_module, "tm_neb", fake_tm_neb)
    monkeypatch.setattr(requests, "post", fake_post)
    try:
        result, _ = plan_flow(
            document(CDS_LIBRARY),
            library=library,
            policy={"constraints": {"primer_len_max": 90}},
            offline=False,
        )
    finally:
        tm_module.clear_tm_caches()

    assert walked and set(walked) == {"phusion-0"}, set(walked)
    assert batched and set(batched) == {"phusion-0"}, set(batched)
    assert result.spec.primer_design.polymerase == "phusion-0"


def test_a_document_that_names_no_polymerase_gets_the_default(library):
    spec, _ = from_flow(
        document(INTEGRATION), library=library, fill={"pos4": library["FP_MRFP1"]}
    )
    assert spec.primer_design.polymerase == "q5-0"


def test_a_hairpin_that_melts_before_the_anneal_is_not_a_breach(library):
    """The bound is relative to where the primer has to work.

    An absolute 45 C cap asked the wrong question: a hairpin melting at 50 C is
    gone by the time a 63 C anneal happens, so flagging it said nothing about
    the reaction.  These are also ΔG -0.6 to -1.9 kcal/mol against a -6.0
    floor, which is the criterion that says whether there is much structure at
    all -- so the cap was firing on structures its own ΔG bound calls trivial.
    """
    result, _ = plan_flow(
        document(CDS_LIBRARY),
        library=library,
        policy={"constraints": {"primer_len_max": 90}},
        offline=True,
    )
    breaches = [d.code for d in result.diagnostics if d.constraint is not None]
    assert "hairpin_tm_vs_anneal" not in breaches, [
        d.message for d in result.diagnostics if d.code == "hairpin_tm_vs_anneal"
    ]

    # The hairpins are still there and still reported; they just do not block.
    hairpins = [p.structures.hairpin_tm for p in result.primers]
    assert max(hairpins) > 45.0, "this design does have hairpins above the old cap"
    assert all(p.structures.hairpin_dg > -6.0 for p in result.primers)


def test_the_margin_moves_the_ceiling_both_ways(library):
    """Zero means "a breach once it survives to the annealing temperature"."""
    import copy

    def breaches(margin):
        raw = copy.deepcopy(document(CDS_LIBRARY))
        result, _ = plan_flow(
            raw,
            library=library,
            policy={
                "constraints": {
                    "primer_len_max": 90,
                    "structure_tm_margin_c": margin,
                }
            },
            offline=True,
        )
        return [d.code for d in result.diagnostics if d.code.endswith("_vs_anneal")]

    assert breaches(0.0) == []
    assert breaches(25.0), "a wide margin pulls the ceiling below the hairpins"
    # Negative pushes the ceiling above the annealing Tm, so nothing can breach.
    assert breaches(-20.0) == []


def test_the_structure_tm_bound_never_crosses_two_engines(library):
    """Both sides of the comparison are primer3's, and the gap is real.

    primer3's Tm for these footprints runs about 2 C above NEB's, so a primer3
    structure Tm checked against the authoritative NEB footprint Tm would be
    the cross-model comparison this package exists to avoid.
    """
    import primer3

    from teemi.design.combinatorial.design import CONDITIONS, anneal_tm

    result, _ = plan_flow(
        document(CDS_LIBRARY),
        library=library,
        policy={"constraints": {"primer_len_max": 90}},
        offline=True,
    )
    gaps = []
    for primer in result.primers:
        own = anneal_tm(primer.footprint, CONDITIONS)
        assert own == pytest.approx(
            primer3.calc_tm(primer.footprint, **CONDITIONS.primer3_kwargs)
        )
        gaps.append(own - primer.footprint_tm.celsius)
    assert max(abs(g) for g in gaps) > 0.5, (
        "if the two engines agreed, keeping them apart would not matter"
    )
