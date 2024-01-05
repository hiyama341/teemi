#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright d notice and this permission notice shall be included in all
# copies or substantial portions of the Software.


from dataclasses import dataclass
from typing import Tuple
from math import fabs
from Bio.SeqRecord import SeqRecord
from teemi.design.cloning import (
    find_all_occurrences_of_a_sequence,
    find_sequence_location,
    add_feature_annotation_to_seqrecord,
)


@dataclass
class CRISPRSequenceCutter:
    """Cuts DNA through CRISPR-cas9 double-stranded break.

    Parameters
    ----------
    sequence : str
        This can be any search query
    chromosome : Bio.SeqRecord
        This could be any sequence to search in.
    """

    sequence: str
    chromosome: SeqRecord

    def __post_init__(self) -> None:
        if len(self.sequence) != 20:
            raise ValueError(
                f"sgRNA must be 20 base pairs and not {len(self.sequence)}"
            )

        # Check that we don't have multiple sgRNA sites
        self.find_all_occurrences_of_a_sequence = find_all_occurrences_of_a_sequence(
            self.sequence, self.chromosome
        )
        if self.find_all_occurrences_of_a_sequence != 1:
            raise ValueError(
                f"sgRNA sites occurring in the chromosome should be 1 and not {len(self.find_all_occurrences_of_a_sequence)}"
            )

        # Find searched sequence
        self.start_end_strand_index = find_sequence_location(
            self.sequence, self.chromosome
        )

        # Add attributes
        (
            self.start_location,
            self.end_location,
            self.strand,
        ) = self.start_end_strand_index

        # CRISPR location
        self.crispr_base17_index = self.get_crispr_db_break_location()
        self.crispr_seq = SeqRecord(
            self.chromosome.seq[self.start_location : self.end_location],
            name=f"sgRNA_{self.sequence.name}",
            id="",
        )

        # Get upstream and downstream sequences - get absolute values and add annotations
        (
            self.upstream_crispr_sequence,
            self.downstream_crispr_sequence,
        ) = self.get_upstream_and_downstream_sequences()
        self.upstream_crispr_sequence.name, self.downstream_crispr_sequence.name = (
            f"Upstream seq of DB",
            f"Downstream seq of DB",
        )

        # Add annotations
        add_feature_annotation_to_seqrecord(
            self.upstream_crispr_sequence,
            label=f"{self.chromosome.name}_upstream_sequence_of_CRISPR_cut",
            strand=self.strand,
        )
        add_feature_annotation_to_seqrecord(
            self.downstream_crispr_sequence,
            label=f"{self.chromosome.name}_downstream_sequence_of_CRISPR_cut",
            strand=self.strand,
        )

    def get_crispr_db_break_location(self) -> int:
        return (
            fabs(self.start_location) + 17
            if self.strand == 1
            else fabs(self.end_location) - 3
        )

    def display_cut(self) -> str:
        break_sequence = f"{str(self.upstream_crispr_sequence.seq)[-20:]}---x----{str(self.downstream_crispr_sequence.seq[:20])}"
        return break_sequence

    def get_upstream_and_downstream_sequences(self) -> Tuple[SeqRecord, SeqRecord]:
        upstream = self.chromosome[: int(fabs(self.crispr_base17_index))]
        downstream = self.chromosome[int(fabs(self.crispr_base17_index)) :]

        return upstream, downstream

    # TODO: add method for calculating x
