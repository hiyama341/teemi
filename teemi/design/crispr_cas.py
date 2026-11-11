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
        'cas3'). 'dcas9' is treated as 'cas9'. Each names one enzyme rather than
        a family, since a PAM and a spacer length only mean anything together
        with the enzyme they were measured for: 'cas9' and 'dcas9' are SpCas9
        (3' NGG, 20 nt), 'cas12a' is As/LbCas12a (5' TTTV, 23 nt) and 'cas3' is
        the Type I-C CASCADE-Cas3 of pCRISPR-Cas3 (5' TTC, 34 nt). The sources
        for each are given with the parameter table in this module.
    step : list of str, optional
        Pipeline components to execute (default is ['find', 'filter']).
    downstream : int, optional
        Number of downstream nucleotides to consider when filtering sgRNAs (default is 15).
    off_target_seed : int, optional
        How many PAM-adjacent bases to use as seed-sequence for off-target hits
        (default is 13). Must be between 1 and the protospacer length that
        ``cas_type`` resolves to, because the finders look the seed up as a
        slice of the protospacer. The 13 is CRISPy-web's window and applies to
        all three nucleases; the resulting count is of exact matches only, where
        CRISPy-web counts near-misses too. Both are discussed with the parameter
        table in this module.
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
        Maximum ``off_target_count`` a guide may have and still be kept by
        ``filter_guides`` (default is 10). It is a number of off-target sites,
        not a number of nucleotides.
    extension_to_promoter_region : int, optional
        Additional parameter to extend the promoter region (default is 0).
    upstream_tss : int, optional
        Number of upstream nucleotides to consider from the TSS (default is 100).
    dwstream_tss : int, optional
        Number of downstream nucleotides to consider from the TSS (default is 100).
    target_non_template_strand : bool, optional
        Whether to target the non-template strand (default is True).
    protospacer_len : int, optional
        Spacer length to scan for. ``None``, the default, means the length the
        chosen ``cas_type`` was measured for -- 20 nt for 'cas9' and 'dcas9',
        23 for 'cas12a', 34 for 'cas3'. Pass a number only to override that.
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
        protospacer_len: Optional[int] = None,
    ):
        if not isinstance(dseqrecord, Dseqrecord):
            raise ValueError("Input must be an instance of Dseqrecord.")
        # Resolving the length here is also what validates ``cas_type``, so an
        # unknown one is rejected on construction rather than part-way through
        # ``extract_sgRNAs``.
        nuclease = _nuclease_for(cas_type)
        if protospacer_len is None:
            protospacer_len = nuclease.protospacer_len
        if not 1 <= off_target_seed <= protospacer_len:
            # The finders look the seed up as a slice of the protospacer, so a
            # seed longer than the protospacer indexes from the wrong end and
            # silently reports no off-targets at all.
            raise ValueError(
                f"off_target_seed must be between 1 and protospacer_len "
                f"({protospacer_len}), got {off_target_seed}. That ceiling is "
                f"the protospacer length for cas_type={cas_type!r}; pass "
                f"protospacer_len to scan for a different length."
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
        'cas3'). 'dcas9' is treated as 'cas9'. Each names one enzyme rather than
        a family, since a PAM and a spacer length only mean anything together
        with the enzyme they were measured for: 'cas9' and 'dcas9' are SpCas9
        (3' NGG, 20 nt), 'cas12a' is As/LbCas12a (5' TTTV, 23 nt) and 'cas3' is
        the Type I-C CASCADE-Cas3 of pCRISPR-Cas3 (5' TTC, 34 nt). The sources
        for each are given with the parameter table in this module.

    Returns
    -------
    off_target_counter : Counter
        A Counter object containing the frequency of each off-target hit.

    """
    # Initialize list to store off-target hits
    off_target_hits = list()
    # The PAM and its side both come from the one parameter table, so this scan
    # and the finders cannot drift apart. ``dcas9`` resolves to ``cas9`` there.
    nuclease = _nuclease_for(cas_type)

    # Find all potential off-target hits in both strands
    for contig in sequences:
        for match in re.finditer(nuclease.off_target_pam_pattern, contig):
            if not nuclease.pam_is_3_prime:
                # A 5' PAM puts the seed immediately downstream of it, which is
                # the register the finder looks up as sgrna[:off_target_seed].
                seed_sequence = contig[match.end() : (match.end() + off_target_seed)]
            else:
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

    ``off_target_pam_pattern`` is the same PAM written for the other scan.
    :func:`find_off_target_hits` reads its seed off the strand it is scanning, so
    it needs Cas9's PAM in the same sense as that seed -- ``(?=GG)`` rather than
    ``(?=CC)``. The two differ for Cas9 alone, and because both scans read both
    strands they collect the same seeds either way. The field exists so that one
    table holds every per-nuclease constant: the two were separate dictionaries
    before, with Cas9's entry written in opposite orientations and nothing
    saying why.

    The messages are preserved verbatim, including which nucleases stay silent:
    two tests assert them on stdout. Printing from a library is a defect in its
    own right and wants replacing with ``warnings``, but not in a change whose
    purpose is to leave behaviour untouched.
    """

    pam_pattern: str
    off_target_pam_pattern: str
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

# Where these numbers come from.
#
# Each entry is one named enzyme, not a family: a PAM and a spacer length mean
# nothing apart from the ortholog or subtype they were measured for. Each is
# cited to a source that states both constants for that enzyme, so a reader can
# check a number without taking this table's word for it.
#
# ``cas9`` is SpCas9, Streptococcus pyogenes Cas9 (Type II-A). Jinek et al. 2012
# (Science 337:816, doi:10.1126/science.1225829) gives both constants: the PAM
# "conforms to an NGG consensus sequence", and the 5'-terminal 20 nt of the
# crRNA are available for target binding. Anders et al. 2014 (Nature 513:569,
# doi:10.1038/nature13579) is the structural read-out of those two G:C pairs.
# NGG is not unique to SpCas9 -- several orthologs recognise NGG or a PAM that
# subsumes it -- so what fixes this entry is NGG paired with a 20 nt spacer,
# plus the fact that the Streptomyces vectors teemi targets carry a
# codon-optimised spCas9 (Tong et al. 2015, ACS Synth Biol 4:1020,
# doi:10.1021/acssynbio.5b00038). ``dcas9`` is the same protein with both
# nuclease domains mutated, for CRISPRi: same PAM, same spacer, no cleavage, and
# the PAM still "indispensible for silencing" (Qi et al. 2013, Cell 152:1173,
# doi:10.1016/j.cell.2013.02.022). Cite Tong et al. 2015, not Qi et al., for the
# D10A/H840A pair -- a 2021 correction records Qi's own second mutation as H841A.
#
# ``cas12a`` is AsCas12a / LbCas12a (Type V-A), from Acidaminococcus sp. BV3L6
# and Lachnospiraceae bacterium ND2006. The 5' TTTV PAM is Toth et al. 2018
# (NAR, doi:10.1093/nar/gky815): "the consensus sequences of the PAM requirement
# of As- and LbCpf1 is TTTV (V: all but T)". Excluding TTTT is deliberate, which
# is why the pattern is ``TTT[ACG]`` and not ``TTT.``. The 23 nt spacer is a
# convention, not a measured optimum: Zetsche et al. 2015 measured the native
# Francisella spacer at 23-25 nt and Tu et al. 2017 call 23-25 nt "commonly used
# as spacer sequence for Cpf1". The one published titration preferred 21 nt, but
# in FnCpf1, not in either ortholog above.
#
# ``cas3`` is the minimal Type I-C CASCADE-Cas3 system, four genes from
# Pseudomonas aeruginosa codon-optimised for streptomycetes as pCRISPR-Cas3
# (Addgene #209427) by Whitford et al. 2025 (NAR 53:gkaf214,
# doi:10.1093/nar/gkaf214) -- the same group as CRISPy-web and teemi itself, and
# the system both constants are taken from. That paper states the PAM ("The
# system recognizes a T-rich 5'-TTC-3' PAM sequence, making it an attractive
# system to use in high GC organisms"), searches spacers as ``5'-TTC-N1-N8-3'``,
# and specifies "the desired 34 nt spacer sequence". TTC does not by itself
# identify I-C, since Type I-B uses it too, so the subtype here is the one being
# modelled rather than one inferred from the PAM. For anyone extending this:
# I-C cuts bi-directionally from the target site (Csorgo et al. 2020), not
# unidirectionally.
#
# ``off_target_seed`` stays one number for all three nucleases, and 13 is a
# CRISPy-web figure rather than a measured property of any enzyme. Blin et al.
# 2020 (Synth Syst Biotechnol 5:99, doi:10.1016/j.synbio.2020.05.005) give the
# window: candidate sgRNAs are sorted by off-target hits "with 0, 1, or 2
# mismatches in the 13 bp seed region right upstream of the PAM". CRISPy-web
# holds that 13 fixed and varies only the PAM when it changes nuclease, so one
# seed length for every enzyme is upstream's design and not an omission here.
# The number predates that paper: CRISPy-web 2016 (Blin et al.,
# doi:10.1016/j.synbio.2016.01.003) counts "direct off-target hits for sequences
# with the same 13 bp core sequence upstream of the PAM", and attributes the
# algorithm to Ronda et al. 2014 (Biotechnol Bioeng 111:1604) -- the origin as
# the 2016 paper gives it, not as a paper read here.
#
# Biophysically there is no one number to pin it to. Hsu et al. 2013
# (doi:10.1038/nbt.2647) report the opposite of a fixed seed: "most bases within
# the 20-bp target site provide varying degrees of specificity". 13 should not be
# defended as a measured Cas9 seed, and it does not need to be -- it is the
# search window of the tool this module's Cas9 output is checked against.
#
# Whitford et al. 2025 do state a Cas3 figure -- "the seed off-target search
# considers the PAM and the PAM proximal first 8 nt", counted "only if a PAM
# site is present", which is exactly what the scan below does, so
# ``off_target_seed=8`` reproduces their metric (8 and not 11: the TTC is
# matched separately, and counting it again would specify 14 bases). It is
# deliberately not the default. Over the first 500 CDSs of NC_003888.3, Cas3 at
# 13 keeps all 8712 candidate guides under the default ``off_target_upper=10``
# with a median count of 0; at 8 it keeps 54.2%, the median count is 9, the 90th
# percentile is 32, and three genes are left with no guide at all. The seed and
# that cutoff are a coupled pair -- Cas9 at seed 10 keeps 51.6%, with a median
# of exactly 10 -- and Whitford's script applies no numeric cutoff, so moving
# one of the two alone does not reproduce their method. Their paper calls the
# seed length an adjustable parameter, which is what it stays here.
#
# The larger gap is in what the count means. CRISPy-web searches its 13 bp
# window allowing up to two or three mismatches; this module counts exact
# 13-mers. ``tests/files_for_testing/crispy_web/cas9_SCO5087_0_mismatches.csv``
# carries the difference in its own columns: guide CY00000196 has 0 off-targets
# at 0 mismatches, 30 at one and 591 at two. ``off_target_count`` is the first
# of those three numbers, so at the shipped defaults ``off_target_upper`` rejects
# almost nothing, and a guide this module reports as unique may not be.
_NUCLEASES: Dict[str, _Nuclease] = {
    "cas9": _Nuclease(
        # SpCas9: 3' NGG, 20 nt spacer. Jinek 2012, 10.1126/science.1225829.
        pam_pattern=r"(?=CC)",
        off_target_pam_pattern=r"(?=GG)",
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
        # As/LbCas12a: 5' TTTV. Toth 2018, 10.1093/nar/gky815. 23 nt spacer is
        # the Cpf1 convention (Zetsche 2015, Tu 2017), not a measured optimum.
        pam_pattern=r"TTT[ACG]",
        off_target_pam_pattern=r"TTT[ACG]",
        protospacer_len=23,
        pam_len=4,
        pam_is_3_prime=False,
        strand_sign=-1,
    ),
    "cas3": _Nuclease(
        # Type I-C CASCADE-Cas3: 5' TTC, 34 nt spacer, both from the
        # pCRISPR-Cas3 paper. Whitford 2025, 10.1093/nar/gkaf214.
        pam_pattern=r"TTC",
        off_target_pam_pattern=r"TTC",
        protospacer_len=34,
        pam_len=3,
        pam_is_3_prime=False,
        strand_sign=-1,
        empty_sgrna_message=(
            "No sgRNA found for locus tag {locus_tag}. Skipping to next locus tag."
        ),
    ),
}


def _nuclease_for(cas_type: str) -> _Nuclease:
    """Look up the parameter set for a ``cas_type``.

    ``dcas9`` is an alias for ``cas9``: catalytically dead Cas9 has the same PAM
    and the same protospacer geometry, it just does not cut.
    """
    key = "cas9" if cas_type == "dcas9" else cas_type
    if key not in _NUCLEASES:
        raise ValueError(
            f"unknown cas_type {cas_type!r}; expected one of "
            f"{sorted(_NUCLEASES)}, or 'dcas9' for CRISPRi"
        )
    return _NUCLEASES[key]


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
    ``protospacer_len``, kept because callers of this function may pass it. The
    23 here is the measured As/LbCas12a convention and matches what a default
    ``SgRNAargs`` now resolves to.
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

        # One lookup, where there were three independent substring guards.
        # Those relied on ``"cas9" in "dcas9"`` for the CRISPRi alias, matched a
        # nuclease anywhere in the string, and left ``sgrna_df`` unbound for
        # anything they did not recognise. They also forwarded
        # ``protospacer_len`` to the cas12a finder alone, which is why a default
        # ``SgRNAargs`` scanned for 23 nt Cas12a spacers everywhere except here,
        # where it scanned for 20.
        sgrna_df = _find_sgrna_hits(
            args.dseqrecord,
            args.strain_name,
            args.locus_tag,
            off_target_counter,
            args.off_target_seed,
            revcomp,
            _nuclease_for(args.cas_type)._replace(protospacer_len=args.protospacer_len),
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
