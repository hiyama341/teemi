#!/usr/bin/env python

# Test teselagen_helpers module

import random

import pandas as pd

# teselagen_helpers only does ``import pydna`` but uses pydna.primer and
# pydna.amplify, so those submodules have to be imported somewhere - we need
# them here anyway to build the inputs.
from pydna.amplify import pcr  # noqa: F401
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer  # noqa: F401

from teemi.design.teselagen_helpers import (
    amplicon_matrix_teselagen,
    primer_matrix_teselagen,
)

NO_FRAGS = 5
FRAG_LENGTH = 90
PRIMER_LENGTH = 20

PARTS_COLUMNS = [
    "Name",
    "Sequence",
    "Overlaps_to",
    "Forward Oligo Name",
    "Reverse Oligo Name",
    "Location",
    "volume",
    "concentration",
    "Q5_fw_tm",
    "Q5_rv_tm",
    "Q5_ta",
]


def _random_dna(length, rng):
    return "".join(rng.choice("ATGC") for _ in range(length))


def _reverse_complement(sequence):
    return str(Dseqrecord(sequence).seq.reverse_complement())


def _make_inputs(no_combs=1, seed=11):
    """Builds fragments plus the teselagen-style primer/parts dataframes.

    Every combination holds ``NO_FRAGS`` fragments. The forward primer is the
    first 20 bp of a fragment and the reverse primer the reverse complement of
    its last 20 bp, so a PCR regenerates the fragment exactly.
    """
    rng = random.Random(seed)
    combinations_matrix = []
    primer_rows = []
    part_rows = []

    for comb_no in range(no_combs):
        fragments = []
        for frag_no in range(NO_FRAGS):
            fragment = Dseqrecord(_random_dna(FRAG_LENGTH, rng))
            fragment.name = "c%df%d" % (comb_no, frag_no)
            fragments.append(fragment)
        combinations_matrix.append(fragments)

        for frag_no, fragment in enumerate(fragments):
            sequence = str(fragment.seq)
            forward_name = "PR_c%df%d_fw" % (comb_no, frag_no)
            reverse_name = "PR_c%df%d_rv" % (comb_no, frag_no)
            primer_rows.append(
                {
                    "Name": forward_name,
                    "Sequence": sequence[:PRIMER_LENGTH],
                    "Location": "o1_c%df%d_fw" % (comb_no, frag_no),
                }
            )
            primer_rows.append(
                {
                    "Name": reverse_name,
                    "Sequence": _reverse_complement(sequence[-PRIMER_LENGTH:]),
                    "Location": "o1_c%df%d_rv" % (comb_no, frag_no),
                }
            )
            # the function looks up the overlapping part for fragment 1 and 4
            overlaps_to = {1: fragments[2].name, 4: fragments[3].name}.get(
                frag_no, "no_overlap"
            )
            part_rows.append(
                {
                    "Name": "PCR_c%df%d" % (comb_no, frag_no),
                    "Sequence": sequence,
                    "Overlaps_to": overlaps_to,
                    "Forward Oligo Name": forward_name,
                    "Reverse Oligo Name": reverse_name,
                    "Location": "l5_c%df%d" % (comb_no, frag_no),
                    "volume": 100.0,
                    "concentration": 10.0,
                    "Q5_fw_tm": 60.0 + frag_no,
                    "Q5_rv_tm": 61.0 + frag_no,
                    "Q5_ta": 63,
                }
            )

    primers = pd.DataFrame(primer_rows)
    parts = pd.DataFrame(part_rows, columns=PARTS_COLUMNS)
    return combinations_matrix, primers, parts


def test_primer_matrix_teselagen_returns_primer_pairs():
    combinations_matrix, primers, parts = _make_inputs()

    primer_matrix = primer_matrix_teselagen(
        primers, parts, combinations_matrix, no_combs=1, no_frags=NO_FRAGS
    )

    assert len(primer_matrix) == 1
    assert len(primer_matrix[0]) == NO_FRAGS

    for frag_no, (forward, reverse) in enumerate(primer_matrix[0]):
        fragment_seq = str(combinations_matrix[0][frag_no].seq)
        assert forward.id == "PR_c0f%d_fw" % frag_no
        assert reverse.id == "PR_c0f%d_rv" % frag_no
        assert str(forward.seq) == fragment_seq[:PRIMER_LENGTH]
        assert str(reverse.seq) == _reverse_complement(fragment_seq[-PRIMER_LENGTH:])
        assert forward.annotations["batches"] == [
            {"location": "o1_c0f%d_fw" % frag_no, "volume": 100, "concentration": 10}
        ]
        assert reverse.annotations["batches"] == [
            {"location": "o1_c0f%d_rv" % frag_no, "volume": 100, "concentration": 10}
        ]


def test_primer_matrix_teselagen_handles_several_combinations():
    combinations_matrix, primers, parts = _make_inputs(no_combs=2)

    primer_matrix = primer_matrix_teselagen(
        primers, parts, combinations_matrix, no_combs=2, no_frags=NO_FRAGS
    )

    assert [len(row) for row in primer_matrix] == [NO_FRAGS, NO_FRAGS]
    assert primer_matrix[0][0][0].id == "PR_c0f0_fw"
    assert primer_matrix[1][0][0].id == "PR_c1f0_fw"
    # the two combinations hold different fragments and therefore differ
    assert str(primer_matrix[0][2][0].seq) != str(primer_matrix[1][2][0].seq)


def test_primer_matrix_teselagen_uses_overlaps_to_for_fragment_1_and_4():
    """Fragments 1 and 4 are disambiguated by the part they overlap with."""
    combinations_matrix, primers, parts = _make_inputs()

    # add decoy parts with the very same sequence but a different overlap -
    # they come first in the dataframe, so they would win without the filter
    decoys = []
    for frag_no in (1, 4):
        sequence = str(combinations_matrix[0][frag_no].seq)
        decoy_fw = "PR_decoy%d_fw" % frag_no
        decoy_rv = "PR_decoy%d_rv" % frag_no
        primers = pd.concat(
            [
                pd.DataFrame(
                    [
                        {
                            "Name": decoy_fw,
                            "Sequence": sequence[:PRIMER_LENGTH],
                            "Location": "decoy_fw",
                        },
                        {
                            "Name": decoy_rv,
                            "Sequence": _reverse_complement(sequence[-PRIMER_LENGTH:]),
                            "Location": "decoy_rv",
                        },
                    ]
                ),
                primers,
            ],
            ignore_index=True,
        )
        decoys.append(
            {
                "Name": "PCR_decoy%d" % frag_no,
                "Sequence": sequence,
                "Overlaps_to": "some_other_part",
                "Forward Oligo Name": decoy_fw,
                "Reverse Oligo Name": decoy_rv,
                "Location": "decoy_location",
                "volume": 100.0,
                "concentration": 10.0,
                "Q5_fw_tm": 50.0,
                "Q5_rv_tm": 50.0,
                "Q5_ta": 50,
            }
        )
    parts = pd.concat(
        [pd.DataFrame(decoys, columns=PARTS_COLUMNS), parts], ignore_index=True
    )

    primer_matrix = primer_matrix_teselagen(
        primers, parts, combinations_matrix, no_combs=1, no_frags=NO_FRAGS
    )

    assert primer_matrix[0][1][0].id == "PR_c0f1_fw"
    assert primer_matrix[0][4][0].id == "PR_c0f4_fw"
    assert primer_matrix[0][1][1].id == "PR_c0f1_rv"
    assert primer_matrix[0][4][1].id == "PR_c0f4_rv"
    # fragment 0 does not use the Overlaps_to column, the decoys do not apply
    assert primer_matrix[0][0][0].id == "PR_c0f0_fw"


def test_amplicon_matrix_teselagen():
    combinations_matrix, primers, parts = _make_inputs()
    primer_matrix = primer_matrix_teselagen(
        primers, parts, combinations_matrix, no_combs=1, no_frags=NO_FRAGS
    )

    amplicon_matrix = amplicon_matrix_teselagen(
        parts, primer_matrix, combinations_matrix, no_combs=1, no_frags=NO_FRAGS
    )

    assert len(amplicon_matrix) == 1
    assert len(amplicon_matrix[0]) == NO_FRAGS

    for frag_no, amplicon in enumerate(amplicon_matrix[0]):
        fragment = combinations_matrix[0][frag_no]
        # the primers amplify the whole fragment
        assert str(amplicon.seq) == str(fragment.seq)
        assert amplicon.name == "PCR_c0f%d" % frag_no
        assert amplicon.annotations["template_name"] == "c0f%d" % frag_no
        assert amplicon.annotations["batches"] == [
            {
                "location": "l5_c0f%d" % frag_no,
                "volume": 100.0,
                "concentration": 10.0,
            }
        ]
        assert amplicon.annotations["ta Q5 Hot Start"] == 63
        assert amplicon.forward_primer.annotations["tm Q5 Hot Start"] == 60.0 + frag_no
        assert amplicon.reverse_primer.annotations["tm Q5 Hot Start"] == 61.0 + frag_no


def test_amplicon_matrix_teselagen_annotates_the_amplicon():
    combinations_matrix, primers, parts = _make_inputs()
    primer_matrix = primer_matrix_teselagen(
        primers, parts, combinations_matrix, no_combs=1, no_frags=NO_FRAGS
    )

    amplicon = amplicon_matrix_teselagen(
        parts, primer_matrix, combinations_matrix, no_combs=1, no_frags=NO_FRAGS
    )[0][0]

    # exactly one PCR_product feature and the two primer_bind features that the
    # function adds itself (the ones added by pydna.amplify.pcr are removed)
    features = [
        (feature.type, int(feature.location.start), int(feature.location.end))
        for feature in amplicon.features
    ]
    assert features == [
        ("PCR_product", 0, FRAG_LENGTH),
        ("primer_bind", 0, PRIMER_LENGTH),
        ("primer_bind", FRAG_LENGTH - PRIMER_LENGTH, FRAG_LENGTH),
    ]
    labels = [feature.qualifiers["label"] for feature in amplicon.features]
    assert labels[1:] == ["PR_c0f0_fw", "PR_c0f0_rv"]


def test_amplicon_matrix_teselagen_several_combinations():
    combinations_matrix, primers, parts = _make_inputs(no_combs=2)
    primer_matrix = primer_matrix_teselagen(
        primers, parts, combinations_matrix, no_combs=2, no_frags=NO_FRAGS
    )

    amplicon_matrix = amplicon_matrix_teselagen(
        parts, primer_matrix, combinations_matrix, no_combs=2, no_frags=NO_FRAGS
    )

    assert [len(row) for row in amplicon_matrix] == [NO_FRAGS, NO_FRAGS]
    assert [amplicon.name for amplicon in amplicon_matrix[1]] == [
        "PCR_c1f%d" % frag_no for frag_no in range(NO_FRAGS)
    ]
    assert str(amplicon_matrix[1][3].seq) == str(combinations_matrix[1][3].seq)
