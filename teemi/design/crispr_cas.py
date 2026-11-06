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

from Bio import SeqIO
import pandas as pd
from typing import Counter
import re
from collections import Counter
from Bio.Seq import Seq
from typing import Callable, Dict, NamedTuple, Optional, List, Tuple
from pydna.dseqrecord import Dseqrecord
from Bio.SeqRecord import SeqRecord
import pandas as pd
from typing import Tuple, List
from Bio import SeqIO


class SgRNAargs:
    """
    A class used to store the parameters for sgRNA generation and filtering.

    Attributes
    ----------
    dseqrecord : Dseqrecord
        Dseqrecord containing the genomic data.
    locus_tag : list of str
        List of locus tags (i.e., genes) to be analyzed.
    cas_type : str
        The type of the CRISPR-Cas system used ('cas9', 'dcas9', 'cas12a' or
        'cas3'). 'dcas9' is treated as 'cas9'.
    step : list of str, optional
        Pipeline components to execute (default is ['find', 'filter']).
    downstream : int, optional
        Number of downstream nucleotides to consider when filtering sgRNAs (default is 15).
    off_target_seed : int, optional
        How many PAM-adjacent bases to use as seed-sequence for off-target hits
        (default is 13). Must be between 1 and ``protospacer_len``, because the
        finders look the seed up as a slice of the protospacer.
    pam_remove : list of str, optional
        List of PAM sequences to filter (default is None).
    downstream_remove : list of str, optional
        List of downstream sequences to filter (default is None).
    sgrna_remove : list of str, optional
        List of specific sgRNA targets to filter (default is None).
    gc_upper : float, optional
        GC cutoff filter (maximum-value) (default is 1).
    gc_lower : float, optional
        GC cutoff filter (minimum-value) (default is 0).
    off_target_upper : int, optional
        Maximum number of seed nucleotides allowed for off-target (default is 10).
    extension_to_promoter_region : int, optional
        Additional parameter to extend the promoter region (default is 0).
    upstream_tss : int, optional
        Number of upstream nucleotides to consider from the TSS (default is 100).
    dwstream_tss : int, optional
        Number of downstream nucleotides to consider from the TSS (default is 100).
    target_non_template_strand : bool, optional
        Whether to target the non-template strand (default is True).
    strain_name : str
        Strain name extracted from the Dseqrecord.
    """

    def __init__(
        self,
        dseqrecord: Dseqrecord,
        locus_tag: List[str],
        cas_type: str = "cas9",
        step: List[str] = ["find", "filter"],
        downstream: int = 15,
        off_target_seed: int = 13,
        pam_remove: Optional[List[str]] = None,
        downstream_remove: Optional[List[str]] = None,
        sgrna_remove: Optional[List[str]] = None,
        gc_upper: float = 1,
        gc_lower: float = 0,
        off_target_upper: int = 10,
        extension_to_promoter_region: int = 0,
        upstream_tss: int = 100,
        dwstream_tss: int = 100,
        target_non_template_strand=True,
        protospacer_len=20,
    ):
        if not isinstance(dseqrecord, Dseqrecord):
            raise ValueError("Input must be an instance of Dseqrecord.")
        if not 1 <= off_target_seed <= protospacer_len:
            # The finders look the seed up as a slice of the protospacer, so a
            # seed longer than the protospacer indexes from the wrong end and
            # silently reports no off-targets at all.
            raise ValueError(
                f"off_target_seed must be between 1 and protospacer_len "
                f"({protospacer_len}), got {off_target_seed}. The Cas3 finder "
                f"uses a 34 nt protospacer, so pass protospacer_len=34 for it."
            )

        self.dseqrecord = dseqrecord
        self.locus_tag = locus_tag if isinstance(locus_tag, list) else [locus_tag]
        self.step = step
        self.downstream = downstream
        self.off_target_seed = off_target_seed
        self.pam_remove = pam_remove
        self.downstream_remove = downstream_remove
        self.sgrna_remove = sgrna_remove
        self.gc_upper = gc_upper
        self.gc_lower = gc_lower
        self.off_target_upper = off_target_upper
        self.cas_type = cas_type
        self.extension_to_promoter_region = extension_to_promoter_region
        self.upstream_tss = upstream_tss
        self.dwstream_tss = dwstream_tss
        self.target_non_template_strand = target_non_template_strand
        self.protospacer_len = protospacer_len

        # Extract strain name from the Dseqrecord
        self.strain_name = self.dseqrecord.id


# TODO incorporate this into the code?
def revcomp(x: str) -> str:
    """
    Compute the reverse complement of a DNA sequence.

    Parameters
    ----------
    x : str
        A string representing a DNA sequence.

    Returns
    -------
    str
        The reverse complement of the input DNA sequence.
    """
    return str(Seq(x).reverse_complement())


def _extract_feature_sequence_as_str(feature, record: Dseqrecord) -> str:
    """
    Extract a feature sequence using a plain Biopython Seq parent sequence.

    This avoids Biopython/pydna incompatibilities for reverse-complement
    operations on compound locations.
    """
    return str(feature.extract(Seq(str(record.seq))))


def _scan_breaks(feature, record: Dseqrecord, coding_sequence: str) -> dict:
    """Offsets a guide window must not cross, per scanned strand.

    ``feature.extract()`` concatenates the parts of a compound location, so a
    window spanning two parts reads a sequence the molecule does not carry, and
    a protospacer or PAM assembled across an intron cannot be targeted. Parts
    that *are* adjacent on the DNA are not boundaries -- a CDS running across
    the origin of a circular record is a compound location whose concatenation
    is the real sequence -- so adjacency is tested, not assumed.

    Keyed by the flag the finders scan under: ``-1`` for the coding sequence and
    ``1`` for its reverse complement, which is the whole sequence reversed, so
    an offset in it is measured from the other end. Mapping the offsets rather
    than slicing each run keeps the finders scanning, and the position formulas
    reading, exactly the strings they always have. Both lists are empty for a
    one-part location, which is every CDS in a bacterial genome.
    """
    parts = list(feature.location.parts)
    record_len = len(record.seq)
    circular = bool(getattr(record, "circular", False))

    offsets, offset = [], 0
    for previous, part in zip(parts, parts[1:]):
        offset += len(previous)
        if part.strand != previous.strand:
            offsets.append(offset)
            continue
        if part.strand == -1:
            # The parts of a minus-strand feature run from high coordinates down.
            adjacent = int(part.end) == int(previous.start)
            wraps = int(previous.start) == 0 and int(part.end) == record_len
        else:
            adjacent = int(part.start) == int(previous.end)
            wraps = int(previous.end) == record_len and int(part.start) == 0
        if not (adjacent or (circular and wraps)):
            offsets.append(offset)

    return {
        -1: offsets,
        1: [len(coding_sequence) - offset for offset in offsets],
    }


def parse_genbank_record(dseqrecord: SeqRecord) -> Tuple[List[str], pd.DataFrame]:
    """
    Parse a Dseqrecord and extract gene information and sequences.

    Parameters
    ----------
    dseqrecord : Dseqrecord
        The Dseqrecord to parse.

    Returns
    -------
    sequences : List[str]
        A list of sequences of positive and negative strands (reverse complement).
    gene_df : pd.DataFrame
        A DataFrame of gene information containing locus tags, gene names, strands, starts, and ends.
    """
    sequences = list()

    # Extract gene information and sequences from the Dseqrecord
    sequences.append(str(dseqrecord.seq))  # Store sequence of positive strand
    sequences.append(
        str(dseqrecord.reverse_complement().seq)
    )  # Store sequence of negative strand (reverse complement)

    return sequences


def find_off_target_hits(
    sequences: List[str], off_target_seed: int, cas_type: str
) -> Counter:
    """
    Find all potential off-target hits in both strands and count their frequencies.

    Parameters
    ----------
    sequences : List[str]
        A list of sequences in which to find off-target hits.
    off_target_seed : int
        The length of the off-target seed sequence to match.
    cas_type : str
        The type of the CRISPR-Cas system used ('cas9', 'dcas9', 'cas12a' or
        'cas3'). 'dcas9' is treated as 'cas9'.

    Returns
    -------
    off_target_counter : Counter
        A Counter object containing the frequency of each off-target hit.

    """
    # Initialize list to store off-target hits
    off_target_hits = list()
    # Define PAM pattern based on cas_type
    pam_patterns = {"cas3": r"TTC", "cas9": r"(?=GG)", "cas12a": r"TTT[ACG]"}

    # dCas9 is catalytically dead Cas9: same PAM, same protospacer geometry.
    cas_type = "cas9" if cas_type == "dcas9" else cas_type
    if cas_type not in pam_patterns:
        raise ValueError(
            f"unknown cas_type {cas_type!r}; expected one of "
            f"{sorted(pam_patterns)}, or 'dcas9' for CRISPRi"
        )

    # Find all potential off-target hits in both strands
    for contig in sequences:
        for match in re.finditer(pam_patterns[cas_type], contig):
            if cas_type == "cas12a":
                # For Cas12a, the seed sequence is immediately downstream of the PAM
                seed_sequence = contig[match.end() : (match.end() + off_target_seed)]
            elif cas_type == "cas3":
                # For Cas3, the seed sequence is immediately downstream of the PAM
                seed_sequence = contig[match.end() : (match.end() + off_target_seed)]
            else:  # CAS9
                # (?=GG) matches at the first G of NGG, so the PAM's N sits at
                # match.start() - 1 and the protospacer's last base just before
                # it. The seed is the off_target_seed bases ending there, which
                # is the register the finder looks up as
                # sgrna[protospacer_len - off_target_seed : protospacer_len].
                pam_n = match.start() - 1
                seed_start = pam_n - off_target_seed
                if seed_start < 0:
                    # Not enough sequence 5' of the PAM. Clamping to 0 here would
                    # invent a seed that runs into the PAM itself.
                    continue
                seed_sequence = contig[seed_start:pam_n]
            if len(seed_sequence) != off_target_seed:
                continue
            # Count the frequency of each off-target hit
            off_target_hits.append(seed_sequence)

    # Count the frequency of each off-target hit
    off_target_counter = Counter(off_target_hits)

    return off_target_counter


class _Nuclease(NamedTuple):
    """Everything that differs between the nucleases when scanning for guides.

    The three finders were near-copies of one another -- cas12a and cas3 were
    96% byte-identical -- and every difference between them reduces to one of
    these fields. Keeping them in one place is what lets a single scan serve
    all three, so a fix lands once rather than three times.

    ``pam_is_3_prime`` says which side of the protospacer the PAM sits on, which
    decides both how the matched window is read and which end of the guide the
    off-target seed is taken from: the seed is always the ``off_target_seed``
    bases nearest the PAM.

    ``strand_sign`` is separate because it does not follow from the PAM side but
    from how the pattern is written. Cas9's ``(?=CC)`` is the reverse complement
    of its NGG PAM, so a match means the guide lies on the opposite strand from
    the one being scanned; cas12a and cas3 match their PAM directly, so a match
    means the guide lies on the scanned strand.

    The messages are preserved verbatim, including which nucleases stay silent:
    two tests assert them on stdout. Printing from a library is a defect in its
    own right and wants replacing with ``warnings``, but not in a change whose
    purpose is to leave behaviour untouched.
    """

    pam_pattern: str
    protospacer_len: int
    pam_len: int
    pam_is_3_prime: bool
    strand_sign: int
    empty_sgrna_message: Optional[str] = None
    short_sgrna_message: Optional[str] = None
    short_pam_message: Optional[str] = None


_CAS9_SHORT_SGRNA_MESSAGE = (
    "sgRNA generated were outside the designated border in {locus_tag}. "
    "To incorporate this extent borders. Skipping to next locus tag."
)

_NUCLEASES: Dict[str, _Nuclease] = {
    "cas9": _Nuclease(
        pam_pattern=r"(?=CC)",
        protospacer_len=20,
        pam_len=3,
        pam_is_3_prime=True,
        strand_sign=1,
        # Cas9 has no empty-protospacer branch of its own: an empty guide
        # falls into the length check, so both cases print the one message.
        empty_sgrna_message=_CAS9_SHORT_SGRNA_MESSAGE,
        short_sgrna_message=_CAS9_SHORT_SGRNA_MESSAGE,
        short_pam_message=(
            "Pam was found outside designated locus_tag: {locus_tag}. "
            "To incorporate this extent borders. Skipping to next locus tag."
        ),
    ),
    "cas12a": _Nuclease(
        pam_pattern=r"TTT[ACG]",
        protospacer_len=23,
        pam_len=4,
        pam_is_3_prime=False,
        strand_sign=-1,
    ),
    "cas3": _Nuclease(
        pam_pattern=r"TTC",
        protospacer_len=34,
        pam_len=3,
        pam_is_3_prime=False,
        strand_sign=-1,
        empty_sgrna_message=(
            "No sgRNA found for locus tag {locus_tag}. Skipping to next locus tag."
        ),
    ),
}

_HIT_COLUMNS = [
    "strain_name",
    "locus_tag",
    "gene_loc",
    "gene_strand",
    "sgrna_strand",
    "sgrna_loc",
    "gc",
    "pam",
    "sgrna",
    "sgrna_seed_sequence",
    "off_target_count",
]


def _find_sgrna_hits(
    record: Dseqrecord,
    strain_name: str,
    locus_tags: List[str],
    off_target_counter: Counter,
    off_target_seed: int,
    revcomp: Callable,
    nuclease: _Nuclease,
) -> pd.DataFrame:
    """Scan the requested CDSs for guides. Shared by all three finders.

    Each CDS is read 5'->3' along the gene and scanned twice, once in that
    orientation and once in its reverse complement. The flag paired with each
    string records which of the two is being scanned -- ``-1`` for the coding
    sequence, ``1`` for its reverse complement -- and not which strand the gene
    is on; the position formulas depend on the former.
    """
    sgrna_hits = []
    protospacer_len = nuclease.protospacer_len
    pam_len = nuclease.pam_len

    for feature in record.features:
        if feature.type != "CDS":
            continue
        locus_tag = feature.qualifiers.get("locus_tag", ["NA"])[0]
        if locus_tag not in locus_tags and "all" not in locus_tags:
            continue

        coding_sequence = _extract_feature_sequence_as_str(feature, record)
        breaks = _scan_breaks(feature, record, coding_sequence)
        location = feature.location
        gene_strand = feature.location.strand

        for flag, scanned in (
            (-1, coding_sequence),
            (1, revcomp(coding_sequence)),
        ):
            for match in re.finditer(nuclease.pam_pattern, scanned):
                # A window spanning two genomically disjoint runs is not a
                # targetable site, so it is skipped before anything is read.
                window_end = match.start() + protospacer_len + pam_len
                if any(match.start() < offset < window_end for offset in breaks[flag]):
                    continue

                window = str(scanned[match.start() : window_end])
                if nuclease.pam_is_3_prime:
                    # The pattern matched the complement of the PAM, so the
                    # guide reads along the other strand.
                    guide = revcomp(window)
                    sgrna = guide[0:protospacer_len]
                    pam = guide[protospacer_len : protospacer_len + pam_len]
                else:
                    sgrna = window[pam_len:]
                    pam = window[:pam_len]

                if not sgrna:
                    if nuclease.empty_sgrna_message:
                        print(nuclease.empty_sgrna_message.format(locus_tag=locus_tag))
                    continue
                if len(sgrna) != protospacer_len:
                    if nuclease.short_sgrna_message:
                        print(nuclease.short_sgrna_message.format(locus_tag=locus_tag))
                    continue
                if len(pam) != pam_len:
                    if nuclease.short_pam_message:
                        print(nuclease.short_pam_message.format(locus_tag=locus_tag))
                    continue

                gc_content = len([b for b in sgrna if b in ("C", "G")]) / len(sgrna)

                # One past the protospacer's end, in coordinates of the CDS
                # as extracted 5'->3': the guide occupies
                # ``coding_sequence[position_sgrna - protospacer_len :
                # position_sgrna]``, read forwards or reverse complemented
                # according to which strand it lies on. The PAM sits at the
                # start of the matched window on the scanned strand whichever
                # side of the protospacer it is on, so the protospacer starts
                # ``pam_len`` into that window either way.
                if flag == 1:
                    position_sgrna = len(scanned) - match.start() - pam_len
                else:
                    position_sgrna = match.start() + pam_len + protospacer_len

                strand_sgrna = nuclease.strand_sign * flag * gene_strand

                # The seed is the bases nearest the PAM, which is the 3' end of
                # the guide for a 3' PAM and the 5' end for a 5' PAM.
                if nuclease.pam_is_3_prime:
                    sgrna_seed = sgrna[
                        (protospacer_len - off_target_seed) : protospacer_len
                    ]
                else:
                    sgrna_seed = sgrna[:off_target_seed]

                sgrna_hits.append(
                    (
                        strain_name,
                        locus_tag,
                        int(location.start) + 1,
                        gene_strand,
                        strand_sgrna,
                        position_sgrna,
                        gc_content,
                        pam,
                        sgrna,
                        sgrna_seed,
                        max(off_target_counter[sgrna_seed] - 1, 0),
                    )
                )

    return pd.DataFrame(sgrna_hits, columns=_HIT_COLUMNS)


def find_sgrna_hits_cas9(
    record: Dseqrecord,
    strain_name: str,
    locus_tags: List[str],
    off_target_counter: Counter,
    off_target_seed: int,
    revcomp: Callable,
) -> pd.DataFrame:
    """Find SpCas9 guides (20 nt protospacer, 3' NGG PAM) in the given CDSs.

    Parameters
    ----------
    record : Dseqrecord
        The record to scan.
    strain_name : str
        Value for the ``strain_name`` column.
    locus_tags : List[str]
        Locus tags to scan, or ``["all"]`` for every CDS.
    off_target_counter : Counter
        Seed frequencies from :func:`find_off_target_hits`.
    off_target_seed : int
        Length of the PAM-proximal seed used for the off-target lookup.
    revcomp : callable
        Reverse-complement function.

    Returns
    -------
    pd.DataFrame
        One row per guide, with the columns listed in ``_HIT_COLUMNS``.
        ``sgrna_loc`` is one past the protospacer's end, in coordinates of the
        CDS as extracted 5'->3'.
    """
    return _find_sgrna_hits(
        record,
        strain_name,
        locus_tags,
        off_target_counter,
        off_target_seed,
        revcomp,
        _NUCLEASES["cas9"],
    )


def find_sgrna_hits_cas12a(
    record: Dseqrecord,
    strain_name: str,
    locus_tags: List[str],
    off_target_counter: Counter,
    off_target_seed: int,
    revcomp: Callable,
    protospacer_len=23,
) -> pd.DataFrame:
    """Find Cas12a guides (5' TTTV PAM) in the given CDSs.

    Takes the same arguments as :func:`find_sgrna_hits_cas9`, plus
    ``protospacer_len``, which is the one nuclease whose protospacer length
    ``extract_sgRNAs`` forwards from ``SgRNAargs`` (whose default is 20, not the
    23 in this signature).
    """
    return _find_sgrna_hits(
        record,
        strain_name,
        locus_tags,
        off_target_counter,
        off_target_seed,
        revcomp,
        _NUCLEASES["cas12a"]._replace(protospacer_len=protospacer_len),
    )


def find_sgrna_hits_cas3(
    record: Dseqrecord,
    strain_name: str,
    locus_tags: List[str],
    off_target_counter: Counter,
    off_target_seed: int,
    revcomp: Callable,
) -> pd.DataFrame:
    """Find Cas3 guides (34 nt protospacer, 5' TTC PAM) in the given CDSs.

    Takes the same arguments as :func:`find_sgrna_hits_cas9`.
    """
    return _find_sgrna_hits(
        record,
        strain_name,
        locus_tags,
        off_target_counter,
        off_target_seed,
        revcomp,
        _NUCLEASES["cas3"],
    )


def filter_guides(args: SgRNAargs, hitframe: pd.DataFrame) -> pd.DataFrame:
    """
    Filter the sgRNAs based on the specified criteria.

    Parameters
    ----------
    args : SgRNAargs
        An object of class SgRNAargs containing the parameters for sgRNA filtering.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the filtered sgRNAs.
    """

    def exclude_rows_based_on_patterns(
        frame: pd.DataFrame, column: str, patterns: List[str]
    ) -> pd.DataFrame:
        """Exclude rows from frame where frame[column] contains any of the provided patterns"""
        if not patterns:
            return frame
        if column not in frame.columns:
            raise ValueError(
                f"cannot filter on {column!r}: the sgRNA hit finders do not "
                f"produce a {column!r} column"
            )
        mask = frame[column].apply(
            lambda x: not any(pattern in x for pattern in patterns)
        )
        return frame[mask]

    # Each mask is built from the frame being filtered rather than from
    # ``hitframe``, so the chain still works once earlier steps have dropped
    # rows, and on the duplicate index that pd.concat gives a multi-contig
    # guide table.
    return (
        hitframe.pipe(exclude_rows_based_on_patterns, "pam", args.pam_remove)
        .pipe(exclude_rows_based_on_patterns, "sgrna", args.sgrna_remove)
        .pipe(exclude_rows_based_on_patterns, "downstream", args.downstream_remove)
        .loc[lambda frame: frame.gc >= args.gc_lower]
        .loc[lambda frame: frame.gc <= args.gc_upper]
        .loc[lambda frame: frame.off_target_count <= args.off_target_upper]
    )


def extract_sgRNAs(args: SgRNAargs) -> Tuple[pd.DataFrame, Counter, pd.DataFrame]:
    """
    Execute all three functions together to extract gene information, off-target hits,
    and sgRNA hits from a given genbank file.

    Parameters
    ----------
    args : SgRNAargs
        An instance of the SgRNAargs class.

    Returns
    -------
    gene_df : pd.DataFrame
        A DataFrame of gene information containing locus tag, gene name, strand, start, end.
    off_target_counter : Counter
        Counter object containing the frequency of each off-target hit.
    sgrna_df : pd.DataFrame
        A DataFrame of sgRNA hits information containing genbank file path, locus tag,
        gene name, strand, offset, position offset, GC content, sgrna, PAM,
        downstream sequence, sgrna_pam_downstream, seed, off-target count.
    """
    # Extract gene information
    sequences = parse_genbank_record(args.dseqrecord)

    # If 'find' is in the steps, execute off-target and sgRNA finding
    if "find" in args.step:
        # Find all potential off-target hits
        off_target_counter = find_off_target_hits(
            sequences, args.off_target_seed, cas_type=args.cas_type
        )

        if "cas3" in args.cas_type:
            # Find all potential sgRNA hits
            sgrna_df = find_sgrna_hits_cas3(
                args.dseqrecord,
                args.strain_name,
                args.locus_tag,
                off_target_counter,
                args.off_target_seed,
                revcomp=revcomp,
            )

        if "cas9" in args.cas_type:
            # Find all potential sgRNA hits
            sgrna_df = find_sgrna_hits_cas9(
                args.dseqrecord,
                args.strain_name,
                args.locus_tag,
                off_target_counter,
                args.off_target_seed,
                revcomp=revcomp,
            )

        if "cas12a" in args.cas_type:
            # Find all potential sgRNA hits
            sgrna_df = find_sgrna_hits_cas12a(
                args.dseqrecord,
                args.strain_name,
                args.locus_tag,
                off_target_counter,
                args.off_target_seed,
                revcomp,
                protospacer_len=args.protospacer_len,
            )

        # Sort sgrna_df by 'off-targets' in ascending order
        sgrna_df.sort_values(by="off_target_count", ascending=True, inplace=True)

    # Filter guides if 'filter' is in the steps
    if "filter" in args.step:
        sgrna_df = filter_guides(args, sgrna_df)

    # use ML models to predict efficient sgRNAs
    if "predict" in args.step and args.cas_type == "cas9":
        print("We will try to implement this part later")

    return sgrna_df
