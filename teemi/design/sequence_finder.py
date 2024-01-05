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
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.


from typing import Tuple, List
from Bio.SeqRecord import SeqRecord
from math import fabs
from dataclasses import dataclass, field
from teemi.design.cloning import (
    find_sequence_location,
    find_all_occurrences_of_a_sequence,
    add_feature_annotation_to_seqrecord,
)


@dataclass
class SequenceFinder:
    """Finds upstream and downstream sequences from a sequence input,
    annotates them and saves them.

    Parameters
    ----------
    sequence : Bio.SeqRecord
        This can be any search query
    chromosome : Bio.SeqRecord
        This could be any sequence to search in.
    """

    sequence: SeqRecord
    chromosome: SeqRecord
    start_location: int = field(init=False)
    end_location: int = field(init=False)
    strand: int = field(init=False)
    upstream_sequence: SeqRecord = field(init=False)
    downstream_sequence: SeqRecord = field(init=False)

    def __post_init__(self) -> None:
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

        # Get upstream and downstream sequences
        if self.strand > 0:
            self.upstream_sequence = self.chromosome[: int(fabs(self.start_location))]
            self.downstream_sequence = self.chromosome[int(fabs(self.end_location)) :]
        else:
            self.upstream_sequence = self.chromosome[: int(fabs(self.end_location))]
            self.downstream_sequence = self.chromosome[int(fabs(self.start_location)) :]

        # Add annotations
        self.upstream_sequence = self._add_feature_annotation_to_seqrecord(
            self.upstream_sequence,
            f"{self.chromosome.name}_upstream_sequence",
            self.strand,
        )
        self.downstream_sequence = self._add_feature_annotation_to_seqrecord(
            self.downstream_sequence,
            f"{self.chromosome.name}_downstream_sequence",
            self.strand,
        )

    def get_upstream_and_downstream_sequences(self) -> Tuple[SeqRecord, SeqRecord]:
        self.upstream_sequence.id = f"Upstream_sequence_{self.chromosome.name}"
        self.upstream_sequence.name = f"Upstream_sequence_{self.chromosome.name}"
        self.downstream_sequence.id = f"Downstream_sequence_{self.chromosome.name}"
        self.downstream_sequence.name = f"Downstream_sequence_{self.chromosome.id}"
        return self.upstream_sequence, self.downstream_sequence

    def find_occurrence_of_sequence_in_target_sequence(self) -> List[int]:
        return find_all_occurrences_of_a_sequence(self.sequence, self.chromosome)

    def _add_feature_annotation_to_seqrecord(
        self, seq: SeqRecord, label: str, strand: int
    ) -> SeqRecord:
        add_feature_annotation_to_seqrecord(seq, label=label, strand=strand)
        return seq
