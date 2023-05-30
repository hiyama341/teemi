#!/usr/bin/env python
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from dataclasses import asdict
from typing import List
import pytest
from teemi.design.sequence_finder import SequenceFinder


def test_sequence_finder_init():
    seq = SeqRecord(Seq("TAGCCAG"), id="seq1", name = 'seq1')
    chrom = SeqRecord(Seq("ATGCTAGCCAGTCGTA"), id="chrom1", name= "chrom1")
    sf = SequenceFinder(seq, chrom)

    assert sf.sequence.name == "seq1"
    assert sf.chromosome.name == "chrom1"
    assert sf.start_location == 4
    assert sf.end_location == 11
    assert sf.strand == 1
    assert str(sf.upstream_sequence.seq) == "ATGC"
    assert str(sf.downstream_sequence.seq) == "TCGTA"


def test_sequence_finder_get_upstream_and_downstream_sequences():
    seq = SeqRecord(Seq("AGC"), id="seq1", name = 'seq1')
    chrom = SeqRecord(Seq("ATGCTAGCCAGTCGTA"), id="chrom1", name = 'chrom1')
    sf = SequenceFinder(seq, chrom)
    upstream, downstream = sf.get_upstream_and_downstream_sequences()

    assert upstream.id == "Upstream_sequence_chrom1"
    assert upstream.name == "Upstream_sequence_chrom1"
    assert str(upstream.seq) == "ATGCT"
    assert downstream.id == "Downstream_sequence_chrom1"
    assert downstream.name == "Downstream_sequence_chrom1"
    assert str(downstream.seq) == "CAGTCGTA"


def test_sequence_finder_find_occurrence_of_sequence_in_target_sequence():
    seq = SeqRecord(Seq("AGC"), id="seq1")
    chrom = SeqRecord(Seq("ATGCTAGCCAGTCGTA"), id="chrom1")
    sf = SequenceFinder(seq, chrom)
    assert sf.find_occurrence_of_sequence_in_target_sequence() == 2
