#!/usr/bin/env python

# Test retrieve_gene_homologs module

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from teemi.design.retrieve_gene_homologs import (
    alignment_identity,
    filter_blast_results,
    codon_optimize_with_dnachisel,
    find_all_starts,
    find_first_in_register_stop,
    all_orfs,
    longest_orf,
)


# --------------------------------------------------------------------------
# alignment_identity
# --------------------------------------------------------------------------
def test_alignment_identity_identical_and_mismatched():
    reference = SeqRecord(Seq("ATGGCCATTGTAATG"), id="ref")  # 15 bp
    identical = SeqRecord(Seq("ATGGCCATTGTAATG"), id="identical")
    one_mismatch = SeqRecord(Seq("ATGGCCATTGTAATT"), id="one_mismatch")

    scores = alignment_identity([identical, one_mismatch], reference)

    # globalmx scores 1 per match, 0 per mismatch/gap, normalized by ref length
    assert scores == [1.0, 14 / 15]


def test_alignment_identity_truncated_query():
    reference = SeqRecord(Seq("ATGGCCATTGTAATG"), id="ref")  # 15 bp
    truncated = SeqRecord(Seq("ATGGCCATTG"), id="truncated")  # first 10 bp

    scores = alignment_identity([truncated], reference)

    # 10 matching bases out of 15 reference bases (gaps are not penalized)
    assert scores == [10 / 15]


def test_alignment_identity_empty_query_list():
    reference = SeqRecord(Seq("ATGGCCATTGTAATG"), id="ref")
    assert alignment_identity([], reference) == []


# --------------------------------------------------------------------------
# filter_blast_results
# --------------------------------------------------------------------------
class FakeHSP:
    """Minimal stand-in for Bio.Blast.Record.HSP."""

    def __init__(self, identities, align_length, expect):
        self.identities = identities
        self.align_length = align_length
        self.expect = expect
        self.query = "M" * align_length
        self.match = "|" * align_length
        self.sbjct = "M" * align_length


class FakeAlignment:
    """Minimal stand-in for Bio.Blast.Record.Alignment."""

    def __init__(self, hit_def, accession, length, hsps):
        self.hit_def = hit_def
        self.title = "gi|000| " + hit_def
        self.accession = accession
        self.length = length
        self.hsps = hsps

    def __str__(self):
        return "FakeAlignment(%s)" % self.hit_def


class FakeBlastRecord:
    """Minimal stand-in for a parsed NCBIXML blast record (no network used)."""

    def __init__(self, alignments):
        self.alignments = alignments


def _blast_record():
    return FakeBlastRecord(
        [
            # passes the default thresholds: e-value < 0.4, 0.1 < identity < 1
            FakeAlignment("Good homolog", "ACC_GOOD", 300, [FakeHSP(80, 100, 1e-20)]),
            # e-value too high
            FakeAlignment("High evalue", "ACC_EVAL", 310, [FakeHSP(70, 100, 0.5)]),
            # identity below the lower threshold
            FakeAlignment("Low identity", "ACC_LOW", 320, [FakeHSP(5, 100, 1e-10)]),
            # identity of exactly 1 is excluded (likely the query itself)
            FakeAlignment("Identical hit", "ACC_SELF", 330, [FakeHSP(100, 100, 0.0)]),
        ]
    )


def test_filter_blast_results_keeps_only_hits_within_thresholds():
    df = filter_blast_results(_blast_record())

    assert list(df.columns) == ["Name", "Identity", "E_value", "Length", "ACC_number"]
    assert len(df) == 1
    row = df.iloc[0]
    assert row["Name"] == "Good homolog"
    assert row["Identity"] == 0.8
    assert row["E_value"] == 1e-20
    assert row["Length"] == 300
    assert row["ACC_number"] == "ACC_GOOD"


def test_filter_blast_results_custom_thresholds_keep_everything():
    df = filter_blast_results(
        _blast_record(),
        E_VALUE_THRESH=1.0,
        LOWER_PROTEIN_IDENTITY_THRESH=0.0,
        UPPER__PROTEIN_IDENTITY_THRESH=1.1,
    )

    assert len(df) == 4
    assert list(df["ACC_number"]) == ["ACC_GOOD", "ACC_EVAL", "ACC_LOW", "ACC_SELF"]
    assert list(df["Identity"]) == [0.8, 0.7, 0.05, 1.0]
    assert list(df["Length"]) == [300, 310, 320, 330]


def test_filter_blast_results_handles_several_hsps_per_alignment():
    record = FakeBlastRecord(
        [
            FakeAlignment(
                "Two hsps",
                "ACC_TWO",
                400,
                [FakeHSP(90, 100, 1e-30), FakeHSP(40, 100, 1e-5)],
            )
        ]
    )

    df = filter_blast_results(record)

    assert len(df) == 2
    assert list(df["Identity"]) == [0.9, 0.4]
    assert set(df["Name"]) == {"Two hsps"}
    assert set(df["ACC_number"]) == {"ACC_TWO"}


def test_filter_blast_results_returns_empty_dataframe_when_nothing_passes():
    record = FakeBlastRecord(
        [FakeAlignment("No good", "ACC_NONE", 100, [FakeHSP(1, 100, 10.0)])]
    )

    df = filter_blast_results(record)

    assert df.empty
    assert list(df.columns) == ["Name", "Identity", "E_value", "Length", "ACC_number"]


def test_filter_blast_results_show_alignment_prints_summary(capsys):
    df = filter_blast_results(_blast_record(), show_alignment=True)
    printed = capsys.readouterr().out

    assert len(df) == 1
    assert "Alignment# 1" in printed
    assert "Name: Good homolog" in printed
    assert "Title: gi|000| Good homolog" in printed
    assert "E value: 1e-20" in printed
    assert "Identitiy 0.80" in printed
    assert "FakeAlignment(Good homolog)" in printed
    assert "TOTAL HOMOLOGS 1" in printed
    # the discarded hits are never printed
    assert "Low identity" not in printed


# --------------------------------------------------------------------------
# codon_optimize_with_dnachisel
# --------------------------------------------------------------------------
# NOTE: the function hands ``seq.seq`` straight to dnachisel, which only accepts
# plain strings. With modern Biopython (Seq is no longer a str subclass) a real
# SeqRecord therefore crashes, so the tests below use a tiny record-like object
# whose ``.seq`` is a string. See the report for the suggested source fix.
class StrRecord:
    """Record-like object exposing the attributes the function reads."""

    def __init__(self, seq, id, name, description):
        self.seq = seq
        self.id = id
        self.name = name
        self.description = description


# 63 bp (21 codons), GC ~40%
GENE = "ATGGCTAGCAAAGGTGAAGAATTATTCACTGGTGTTGTCCCAATTTTGGTTGAATTAGATTAA"


def _gc_fraction(sequence):
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def _window_gc_within(sequence, window, lower, upper):
    return all(
        lower <= _gc_fraction(sequence[i : i + window]) <= upper
        for i in range(0, len(sequence) - window + 1)
    )


def test_codon_optimize_with_dnachisel_requires_species_or_table():
    with pytest.raises(ValueError):
        codon_optimize_with_dnachisel([StrRecord(GENE, "g1", "gene1", "a gene")])


def test_codon_optimize_with_dnachisel_species(capsys):
    record = StrRecord(GENE, "g1", "gene1", "a nice gene")

    optimized = codon_optimize_with_dnachisel([record], species="s_cerevisiae", window=20)

    assert len(optimized) == 1
    result = optimized[0]
    # metadata is carried over from the input record
    assert (result.id, result.name, result.description) == ("g1", "gene1", "a nice gene")
    # the sequence keeps its length and only contains DNA letters
    sequence = str(result.seq)
    assert len(sequence) == len(GENE)
    assert set(sequence) <= set("ATGC")
    # the GC constraint that was enforced is actually satisfied
    assert _window_gc_within(sequence, 20, 0.3, 0.7)
    # dnachisel annotates the record with the constraint/objective it used
    labels = [f.qualifiers["label"] for f in result.features]
    assert "@30-70% GC/20bp" in labels
    assert "~best-codon-optimize (s_cerevisiae)" in labels
    # the summaries are printed by the function
    printed = capsys.readouterr().out
    assert "EnforceGCContent" in printed
    assert "MaximizeCAI" in printed


def test_codon_optimize_with_dnachisel_custom_codon_table():
    import python_codon_tables

    table = python_codon_tables.get_codons_table("s_cerevisiae")
    record = StrRecord(GENE, "g2", "gene2", "another gene")

    optimized = codon_optimize_with_dnachisel(
        [record],
        codon_usage_table=table,
        lower_GC=0.35,
        upper_GC=0.65,
        window=20,
    )

    assert len(optimized) == 1
    result = optimized[0]
    assert (result.id, result.name, result.description) == ("g2", "gene2", "another gene")
    sequence = str(result.seq)
    assert len(sequence) == len(GENE)
    assert _window_gc_within(sequence, 20, 0.35, 0.65)
    labels = [f.qualifiers["label"] for f in result.features]
    assert "@35-65% GC/20bp" in labels
    assert "~best-codon-optimize" in labels


def test_codon_optimize_with_dnachisel_several_sequences():
    records = [
        StrRecord(GENE, "g1", "gene1", "first"),
        StrRecord(GENE[:30], "g2", "gene2", "second"),
    ]

    optimized = codon_optimize_with_dnachisel(records, species="e_coli", window=20)

    assert [r.id for r in optimized] == ["g1", "g2"]
    assert [len(r.seq) for r in optimized] == [len(GENE), 30]
    for result in optimized:
        assert _window_gc_within(str(result.seq), 20, 0.3, 0.7)


# --------------------------------------------------------------------------
# find_all_starts / find_first_in_register_stop / all_orfs / longest_orf
# --------------------------------------------------------------------------
def test_find_all_starts():
    # start codons at index 2, 11 and 14 (the last two overlap)
    assert find_all_starts("ccatgcccTGA".lower()) == (2,)
    assert find_all_starts("ccatgccctgaatgatgttttaa") == (2, 11, 14)


def test_find_all_starts_without_any_start_codon():
    assert find_all_starts("cccaaatttggg") == ()


def test_find_first_in_register_stop():
    # "taa" starts at index 3 and ends at 6, which is divisible by three
    assert find_first_in_register_stop("ccctaa") == 6


def test_find_first_in_register_stop_skips_out_of_register_stops():
    # "tga" ends at index 4 (out of register), "taa" ends at 9 (in register)
    assert find_first_in_register_stop("ctgaattaa") == 9


def test_find_first_in_register_stop_without_in_register_stop():
    # the only stop codon ends at index 4
    assert find_first_in_register_stop("ataa") == -1


# two ORFs: "atgccctga" at 2-11 and "atgatgttttaa" at 11-23. The start codon at
# index 14 shares the stop at 23 and is therefore discarded as the shorter ORF.
ORF_SEQ = "ccatgccctgaatgatgttttaa"


def test_all_orfs_sorted_by_decreasing_length():
    assert all_orfs(ORF_SEQ) == ((11, 23), (2, 11))


def test_all_orfs_is_case_insensitive():
    assert all_orfs(ORF_SEQ.upper()) == all_orfs(ORF_SEQ)


def test_all_orfs_without_stop_codon_in_register():
    assert all_orfs("atgaaacccggg") == ()


def test_all_orfs_without_start_codon():
    assert all_orfs("cccaaatttggg") == ()


def test_longest_orf():
    assert longest_orf(ORF_SEQ) == "atgatgttttaa"


def test_longest_orf_keeps_the_case_of_the_input():
    assert longest_orf(ORF_SEQ.upper()) == "ATGATGTTTTAA"


def test_longest_orf_several_orfs():
    assert longest_orf(ORF_SEQ, n=2) == ("atgatgttttaa", "atgccctga")
    # n larger than the number of ORFs simply returns all of them
    assert longest_orf(ORF_SEQ, n=5) == ("atgatgttttaa", "atgccctga")


def test_longest_orf_with_a_single_orf_returns_a_string():
    single = "ccatgaaataa"  # one ORF only: atgaaataa
    assert all_orfs(single) == ((2, 11),)
    assert longest_orf(single, n=3) == "atgaaataa"


def test_longest_orf_without_orfs():
    assert longest_orf("cccaaatttggg") == ""
