#!/usr/bin/env python

"""Tests for teemi.design.combinatorial_design."""

import math
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from pandas.testing import assert_frame_equal
from pydna.assembly import Assembly
from pydna.dseqrecord import Dseqrecord
from pydna.tm import tm_default

from teemi.design.combinatorial_design import DesignAssembly, count_unique_parts

PARTS = Path(__file__).parent / "files_for_testing"


def load_parts(filename, count=2):
    records = SeqIO.parse(PARTS / filename, "fasta")
    return [Dseqrecord(record) for record in list(records)[:count]]


@pytest.fixture(scope="module")
def promoters():
    return load_parts("promoter_test_data.fasta")


@pytest.fixture(scope="module")
def cds():
    return load_parts("CDS_test_data.fasta")


@pytest.fixture(scope="module")
def terminators():
    return load_parts("terminator_test_data.fasta")


@pytest.fixture(scope="module")
def library(promoters, cds, terminators):
    """A 2 x 2 x 2 promoter-CDS-terminator library."""
    return DesignAssembly([promoters, cds, terminators])


def homology(upstream, downstream):
    """Length of the sequence the two amplicons share at their junction."""
    upstream, downstream = str(upstream.seq), str(downstream.seq)
    for length in range(min(len(upstream), len(downstream)), 0, -1):
        if upstream[-length:] == downstream[:length]:
            return length
    return 0


def test_library_covers_every_combination(library):
    assert library.combinatorial_list_of_names == [
        ("MLS1", "AKF02530", "ADH1"),
        ("MLS1", "AKF02530", "CYC1"),
        ("MLS1", "AAA17732", "ADH1"),
        ("MLS1", "AAA17732", "CYC1"),
        ("URA1", "AKF02530", "ADH1"),
        ("URA1", "AKF02530", "CYC1"),
        ("URA1", "AAA17732", "ADH1"),
        ("URA1", "AAA17732", "CYC1"),
    ]
    assert library.systematic_names == [
        (1, 1, 1),
        (1, 1, 2),
        (1, 2, 1),
        (1, 2, 2),
        (2, 1, 1),
        (2, 1, 2),
        (2, 2, 1),
        (2, 2, 2),
    ]
    assert [len(assembly) for assembly in library.list_of_assemblies] == [3] * 8


def test_variants_dataframe(library):
    assert_frame_equal(
        library.show_variants_lib_df(),
        pd.DataFrame(
            {
                0: ["MLS1", "MLS1", "MLS1", "MLS1", "URA1", "URA1", "URA1", "URA1"],
                1: ["AKF02530", "AKF02530", "AAA17732", "AAA17732"] * 2,
                2: ["ADH1", "CYC1"] * 4,
                "Systematic_name": [
                    (1, 1, 1),
                    (1, 1, 2),
                    (1, 2, 1),
                    (1, 2, 2),
                    (2, 1, 1),
                    (2, 1, 2),
                    (2, 2, 1),
                    (2, 2, 2),
                ],
                "Variant": range(8),
            }
        ),
    )


def test_amplicons_assemble_into_one_contig_of_the_expected_length(
    library, promoters, cds, terminators
):
    """Overlapping primer tails must fuse the parts without losing or doubling sequence."""
    for variant, assembly in enumerate(library.list_of_assemblies):
        parts = [
            promoters[library.systematic_names[variant][0] - 1],
            cds[library.systematic_names[variant][1] - 1],
            terminators[library.systematic_names[variant][2] - 1],
        ]
        contigs = Assembly(assembly, limit=30).assemble_linear()
        assert len(contigs) == 1
        assert len(contigs[0]) == sum(len(part) for part in parts)
        assert str(contigs[0].seq) == "".join(str(part.seq) for part in parts)


def test_neighbouring_amplicons_share_the_requested_overlap(
    promoters, cds, terminators
):
    library = DesignAssembly([promoters, cds, terminators], overlap=60)
    for assembly in library.list_of_assemblies:
        for upstream, downstream in zip(assembly, assembly[1:]):
            assert homology(upstream, downstream) >= 60


def test_primers_to_order(library):
    """The orderable content of the primer table: identity, template and sequence."""
    primers = library.primer_list_to_dataframe()
    assert list(primers.columns) == [
        "id",
        "anneals to",
        "sequence",
        "annealing temperature",
        "length",
        "price(DKK)",
        "description",
        "footprint",
        "len_footprint",
    ]
    assert list(
        primers[["id", "anneals to", "sequence"]].itertuples(index=False, name=None)
    ) == [
        ("P001", "MLS1", "TTTAATCTTTAGGGAGGGTAAAG"),
        ("P002", "MLS1", "TTCCATTTCATTATCCATTTTCTTAATTCTTTTATGTGCTTTTACT"),
        ("P003", "AKF02530", "CATAAAAGAATTAAGAAAATGGATAATGAAATGGAAACTATG"),
        ("P004", "AKF02530", "AAATCATAAGAAATTCGCATAAGAATCTGGATTATTTTACATAACT"),
        ("P005", "ADH1", "AATAATCCAGATTCTTATGCGAATTTCTTATGATTTATGATTTT"),
        ("P006", "ADH1", "CGTAAAAAAAGCATGCACG"),
        ("P007", "AKF02530", "AAAGGAAAAGGGGCCTGTATAAGAATCTGGATTATTTTACATAACT"),
        ("P008", "CYC1", "AATAATCCAGATTCTTATACAGGCCCCTTTTCC"),
        ("P009", "CYC1", "GTCGACAACTAAACTGGAATG"),
        ("P010", "MLS1", "TTACATTTCCATATCCATTTTCTTAATTCTTTTATGTGCTTTTACT"),
        ("P011", "AAA17732", "CATAAAAGAATTAAGAAAATGGATATGGAAATGTAAACTATG"),
        ("P012", "URA1", "GTTGTATTAATTTTCTCGAAGGG"),
        ("P013", "URA1", "TTCCATTTCATTATCCATGTTTGGTACGGAAGTTCAA"),
        ("P014", "AKF02530", "TGAACTTCCGTACCAAACATGGATAATGAAATGGAAACTATG"),
        ("P015", "URA1", "TTACATTTCCATATCCATGTTTGGTACGGAAGTTCAA"),
        ("P016", "AAA17732", "TGAACTTCCGTACCAAACATGGATATGGAAATGTAAACTATG"),
    ]


def test_primer_table_rows_are_consistent(library):
    for _, row in library.primer_list_to_dataframe().iterrows():
        assert row["length"] == len(row["sequence"])
        assert row["price(DKK)"] == pytest.approx(row["length"] * 1.8)
        # The 5' tail carries the overlap; only the 3' footprint anneals.
        assert row["sequence"].endswith(row["footprint"])
        assert row["len_footprint"] == len(row["footprint"])
        assert row["annealing temperature"] == pytest.approx(
            round(tm_default(row["footprint"]), 2)
        )


def test_primer_descriptions_name_the_parts_they_join(library):
    descriptions = set(library.primer_list_to_dataframe()["description"])
    assert "Anneals to MLS1" in descriptions
    assert "Anneals to MLS1, overlaps to AKF02530" in descriptions
    assert "Anneals to AKF02530, overlaps to MLS1" in descriptions
    # pydna names unlabelled PCR products "<n>bp_PCR_prod"; real part names belong here.
    assert not [text for text in descriptions if "PCR_prod" in text]


def test_pcrs_to_run(library):
    assert_frame_equal(
        library.pcr_list_to_dataframe(),
        pd.DataFrame(
            {
                "pcr_number": [f"PCR{number}" for number in range(1, 15)],
                "template": [
                    "MLS1",
                    "AKF02530",
                    "ADH1",
                    "AKF02530",
                    "CYC1",
                    "MLS1",
                    "AAA17732",
                    "AAA17732",
                    "URA1",
                    "AKF02530",
                    "AKF02530",
                    "URA1",
                    "AAA17732",
                    "AAA17732",
                ],
                "forward_primer": [
                    "P001",
                    "P003",
                    "P005",
                    "P003",
                    "P008",
                    "P001",
                    "P011",
                    "P011",
                    "P012",
                    "P014",
                    "P014",
                    "P012",
                    "P016",
                    "P016",
                ],
                "reverse_primer": [
                    "P002",
                    "P004",
                    "P006",
                    "P007",
                    "P009",
                    "P010",
                    "P004",
                    "P007",
                    "P013",
                    "P004",
                    "P007",
                    "P015",
                    "P004",
                    "P007",
                ],
                "f_tm": [
                    55.59,
                    54.77,
                    54.72,
                    54.77,
                    56.44,
                    55.59,
                    54.29,
                    54.29,
                    55.52,
                    54.77,
                    54.77,
                    55.52,
                    54.29,
                    54.29,
                ],
                "r_tm": [
                    56.13,
                    54.64,
                    56.27,
                    54.64,
                    57.12,
                    56.13,
                    54.64,
                    54.64,
                    55.64,
                    54.64,
                    54.64,
                    55.64,
                    54.64,
                    54.64,
                ],
            }
        ),
    )


def test_variants_reuse_primers_and_pcrs_instead_of_duplicating_them(library):
    """8 variants need 24 amplicons, but only the distinct ones have to be made."""
    assert len(library.list_of_assemblies) * 3 == 24
    assert len(library.amplicons) == 14
    assert len(library.primers) == 16
    assert len({primer.id for primer in library.primers}) == 16
    # Both CDS parts end in the same footprint, so one reverse primer serves several PCRs.
    pcrs = library.pcr_list_to_dataframe()
    assert sorted(pcrs.loc[pcrs["reverse_primer"] == "P004", "pcr_number"]) == [
        "PCR10",
        "PCR13",
        "PCR2",
        "PCR7",
    ]


def test_tm_func_is_used_for_primer_design_and_reporting(promoters, cds):
    def hot(sequence):
        return tm_default(sequence) + 10

    library = DesignAssembly([promoters, cds], target_tm=65.0, tm_func=hot)
    primers = library.primer_list_to_dataframe()
    for _, row in primers.iterrows():
        assert row["annealing temperature"] == pytest.approx(
            round(hot(row["footprint"]), 2)
        )
    assert primers["annealing temperature"].between(63, 67).all()


class TestPads:
    def test_pad_is_carried_by_the_primers_into_every_construct(
        self, library, promoters, cds, terminators
    ):
        kozak = Dseqrecord("TCGGTC", name="kozak")
        padded = DesignAssembly([promoters, cds, terminators], [kozak], [1])

        # A pad rides in the primer tails, so it adds no PCR of its own.
        assert [len(assembly) for assembly in padded.list_of_assemblies] == [3] * 8
        assert len(padded.amplicons) == len(library.amplicons)

        for variant, assembly in enumerate(padded.list_of_assemblies):
            contig = Assembly(assembly, limit=30).assemble_linear()[0]
            unpadded = Assembly(
                library.list_of_assemblies[variant], limit=30
            ).assemble_linear()[0]
            assert len(contig) == len(unpadded) + len(kozak)
            assert str(contig.seq).count(str(kozak.seq)) == 1
            promoter = promoters[padded.systematic_names[variant][0] - 1]
            assert str(contig.seq).startswith(str(promoter.seq) + str(kozak.seq))

    @pytest.mark.parametrize("position", [0, 1, 2, 3])
    def test_pad_can_sit_at_any_position(self, promoters, cds, terminators, position):
        kozak = Dseqrecord("TCGGTC", name="kozak")
        padded = DesignAssembly([promoters, cds, terminators], [kozak], [position])
        contig = Assembly(padded.list_of_assemblies[0], limit=30).assemble_linear()[0]
        expected = [str(promoters[0].seq), str(cds[0].seq), str(terminators[0].seq)]
        expected.insert(position, str(kozak.seq))
        assert str(contig.seq) == "".join(expected)

    def test_several_pads(self, promoters, cds, terminators):
        pads = [Dseqrecord("TCGGTC"), Dseqrecord("AAGCTT")]
        padded = DesignAssembly([promoters, cds, terminators], pads, [1, 2])
        contig = Assembly(padded.list_of_assemblies[0], limit=30).assemble_linear()[0]
        assert str(contig.seq) == (
            str(promoters[0].seq)
            + "TCGGTC"
            + str(cds[0].seq)
            + "AAGCTT"
            + str(terminators[0].seq)
        )


class TestInvalidInput:
    def test_no_parts(self):
        with pytest.raises(ValueError, match="one non-empty list of parts"):
            DesignAssembly([])

    def test_position_without_parts(self, promoters):
        with pytest.raises(ValueError, match="one non-empty list of parts"):
            DesignAssembly([promoters, []])

    def test_biopython_records_are_rejected(self):
        records = list(SeqIO.parse(PARTS / "promoter_test_data.fasta", "fasta"))[:1]
        with pytest.raises(TypeError, match="must be Dseqrecord"):
            DesignAssembly([records])

    def test_pads_without_positions(self, promoters, cds):
        with pytest.raises(ValueError, match="same length"):
            DesignAssembly([promoters, cds], [Dseqrecord("TCGGTC")], [])

    def test_pad_too_long_to_fit_in_a_primer(self, promoters, cds):
        with pytest.raises(ValueError, match="do not fit in the primer tails"):
            DesignAssembly([promoters, cds], [Dseqrecord("A" * 41)], [1])

    def test_pad_position_outside_the_construct(self, promoters, cds):
        with pytest.raises(ValueError, match="outside the construct"):
            DesignAssembly([promoters, cds], [Dseqrecord("TCGGTC")], [3])


class TestCountUniqueParts:
    @pytest.fixture
    def predictions(self):
        """Ranked predictions, best first."""
        return pd.DataFrame(
            {
                "G8H": ["Rsep", "Smus", "Rsep", "Vmin", "Smus"],
                "pG8H": ["ENO2", "ENO2", "ENO2", "ENO2", "ENO2"],
                "CPR": ["Clo", "Cac", "Ara", "Ara", "Aan"],
                "pCPR": ["TPI1", "TPI1", "TPI1", "TPI1", "CCW12"],
            }
        )

    def test_stops_as_soon_as_the_parts_span_the_target(self, predictions):
        """Reading a 4th prediction would add a 3rd G8H and overshoot 6 combinations."""
        assert count_unique_parts(predictions, 6) == {
            "G8H": ["Rsep", "Smus"],
            "pG8H": ["ENO2"],
            "CPR": ["Clo", "Cac", "Ara"],
            "pCPR": ["TPI1"],
            "sum_of_parts": "6",
            "prediction_number": "3",
        }

    def test_reported_parts_span_the_reported_number_of_combinations(self, predictions):
        parts = count_unique_parts(predictions, 6)
        collected = [v for k, v in parts.items() if isinstance(v, list)]
        assert math.prod(len(part) for part in collected) == int(parts["sum_of_parts"])

    def test_first_prediction_is_enough_for_a_single_combination(self, predictions):
        assert count_unique_parts(predictions, 1) == {
            "G8H": ["Rsep"],
            "pG8H": ["ENO2"],
            "CPR": ["Clo"],
            "pCPR": ["TPI1"],
            "sum_of_parts": "1",
            "prediction_number": "1",
        }

    def test_any_columns_are_accepted(self):
        predictions = pd.DataFrame({"promoter": ["pA", "pB"], "gene": ["g1", "g1"]})
        assert count_unique_parts(predictions, 2) == {
            "promoter": ["pA", "pB"],
            "gene": ["g1"],
            "sum_of_parts": "2",
            "prediction_number": "2",
        }

    def test_too_few_predictions_to_reach_the_target(self, predictions):
        with pytest.raises(ValueError, match="span only 24 combinations"):
            count_unique_parts(predictions, 1000)
