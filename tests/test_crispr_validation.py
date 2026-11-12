#!/usr/bin/env python
"""Validate the CRISPR modules against references produced outside teemi.

The unit tests elsewhere check that the code does what it is written to do.
These check that what it is written to do is right, against references produced
outside this code base:

* **CRISPy-web** output for *Streptomyces coelicolor* A3(2) SCO5087, exported
  from the web tool (``tests/files_for_testing/crispy_web/``). CRISPy-web is a
  sound oracle for Cas9 guide lists and, on this plus-strand gene, for
  base-editing mutation calls.
* **A hand curation of SCO7048**, which is on the minus strand. CRISPy-web is
  not reliable for base editing there, so those rows were worked out by hand.
  This is the reference that matters most: placing an edit on a minus-strand
  gene is what the ``sgrna_strand == gene_strand`` comparison in
  ``process_base_editing`` exists to get right, and reading ``sgrna_strand``
  on its own puts every one of these mutations on the wrong residue.

* **GenScript's gRNA design tool** run with LbCas12a over two *S. cerevisiae*
  genes (``tests/files_for_testing/genscript/``). Cas12a had no external
  reference at all: CRISPy-web's export is Cas9, and *S. coelicolor* is so GC
  rich that SCO5087, the gene the Cas9 reference covers, carries no TTTV site
  on either strand. These two genes were picked instead because they are on
  opposite strands, so between them they exercise both position formulas.

* **The spacer-finding notebook published with pCRISPR-Cas3**
  (``tests/files_for_testing/type1_protospacers/``). Whitford et al. 2025 wrote
  their own script for this, because "the majority of existing computational
  tools do not support type I/type I-C CRISPR systems", and released it as
  `NBChub/Type1_Protospacers <https://github.com/NBChub/Type1_Protospacers>`_.
  It is the reference implementation of the system teemi models as ``cas3``,
  and no vendor tool offers Type I-C at all. Unlike the three above it is the
  same laboratory as teemi, so what it settles is the mechanics rather than the
  biology. Its rows here were transcribed from the notebook's search loop
  (commit ``d1c2aab``) rather than exported from a service, which is why each
  one is re-derived from the chromosome before it is compared.

The genome file is the same one StreptoCAD validates against, byte for byte.
"""

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from pydna.dseqrecord import Dseqrecord

from teemi.design.crispr_best import (
    filter_sgrnas_for_base_editing,
    identify_base_editing_sites,
    process_base_editing,
)
from teemi.design.crispr_cas import (
    _NUCLEASES,
    SgRNAargs,
    _scan_breaks,
    extract_sgRNAs,
)

GENOME = "tests/files_for_testing/gibson/Streptomyces_coelicolor_A3_chromosome.gb"
CRISPY_DIR = "tests/files_for_testing/crispy_web"
GENSCRIPT_DIR = "tests/files_for_testing/genscript"

# GenScript designs against whole chromosomes and reports chromosome
# coordinates, so each record here is that gene's own slice of its chromosome
# and ``offset`` converts back: chromosome position = position in record + offset.
GENSCRIPT_CAS12A = [
    # gene, locus tag, record, export, chromosome offset, teemi guides
    ("CAN1", "YEL063C", "CAN1_YEL063C.gb", "cas12a_CAN1_lbcas12a.csv", 31599, 81),
    ("LYP1", "YNL268W", "LYP1_YNL268W.gb", "cas12a_LYP1_lbcas12a.csv", 138449, 67),
]

TYPE1_DIR = "tests/files_for_testing/type1_protospacers"
TYPE1_CAS3 = [
    # locus tag, reference extract, spacers
    ("SCO5087", "cas3_SCO5087.csv", 31),
    ("SCO7048", "cas3_SCO7048.csv", 31),
]
# The seed length the notebook ships, and the one its off-target column counts.
TYPE1_SEED = 8
# The spacers the paper reports cloning, as opposed to ones a scan proposes.
TYPE1_CLONED = "cas3_cloned_spacers.csv"

# Hand-curated base-editing results for SCO7048, a minus-strand CDS: the guide
# and the amino-acid changes its C-to-T edits produce.
CURATED_SCO7048 = [
    ("TCCCGGCGAATAGCCGGCGC", "R26W, R27*"),
    ("CTCCCGGCGAATAGCCGGCG", "R26W, R27*"),
    ("CCCGGCGAATAGCCGGCGCG", "R27*"),
    ("CTGGCTCCAGCTGTCCCCCA", "W66*, S67N"),
    ("GCTGGCTCCAGCTGTCCCCC", "W66*, S67N"),
    ("CTGGAGCCAGCTGTTCGTCG", "Q68*"),
    ("ACCGGCAGGTCTTCACGCGC", "Q93*"),
    ("GCGACCACACCAGGATGTCG", "V146I, W147*"),
    ("CGACCACACCAGGATGTCGG", "V146I, W147*"),
    ("CGCGACCACACCAGGATGTC", "W147*"),
    ("CCGCGACCACACCAGGATGT", "W147*"),
    ("GGGCAGGCGAGCCGGGCGAC", "Q151*, A152V"),
    ("CCGTCCACCAGGGACTGTTC", "H231Y, Q232*"),
    ("GCCGTCCACCAGGGACTGTT", "H231Y, Q232*"),
    ("CCGGGCCCAGGCCATCACCG", "A305V, Q306*"),
    ("CTACCACCAGACGGGGCTGC", "H321Y, Q322*"),
]


def _split(mutations):
    return [m.strip() for m in str(mutations).split(",") if m.strip()]


@pytest.fixture(scope="module")
def coelicolor():
    return Dseqrecord(SeqIO.read(GENOME, "genbank"))


def _gene_sequences(record, locus_tags):
    """Map each locus tag to its CDS read 5'->3'."""
    return {
        feature.qualifiers["locus_tag"][0]: str(feature.extract(Seq(str(record.seq))))
        for feature in record.features
        if feature.type == "CDS"
        and feature.qualifiers.get("locus_tag", [""])[0] in locus_tags
    }


def _guides(record, locus_tag):
    return extract_sgRNAs(
        SgRNAargs(
            dseqrecord=record,
            locus_tag=[locus_tag],
            cas_type="cas9",
            gc_upper=0.9999,
            gc_lower=0.0001,
            off_target_seed=13,
            off_target_upper=10,
            step=["find", "filter"],
        )
    )


def _stop_codon_guides(record, locus_tag):
    guides = identify_base_editing_sites(_guides(record, locus_tag))
    editable = filter_sgrnas_for_base_editing(guides)
    return process_base_editing(
        editable,
        gene_sequences=_gene_sequences(record, [locus_tag]),
        only_stop_codons=True,
        editing_context=False,
    ).reset_index(drop=True)


@pytest.fixture(scope="module")
def sco5087_stops(coelicolor):
    return _stop_codon_guides(coelicolor, "SCO5087")


@pytest.fixture(scope="module")
def sco7048_stops(coelicolor):
    return _stop_codon_guides(coelicolor, "SCO7048")


# ---------------------------------------------------------------------------
# Base editing against the hand curation of a MINUS-strand gene
# ---------------------------------------------------------------------------


def test_base_editing_matches_hand_curation_on_a_minus_strand_gene(sco7048_stops):
    curated = pd.DataFrame(CURATED_SCO7048, columns=["sgrna", "mutations"])

    merged = sco7048_stops[["sgrna", "mutations"]].merge(
        curated, on="sgrna", suffixes=("_teemi", "_curated")
    )

    # The curation is exhaustive -- it is every stop-codon guide for SCO7048 --
    # so assert set equality, not containment. An inner join on its own passes
    # silently both when teemi stops producing these guides and when it starts
    # producing extra ones.
    assert set(sco7048_stops["sgrna"]) == set(curated["sgrna"]), (
        f"missing: {sorted(set(curated['sgrna']) - set(sco7048_stops['sgrna']))}; "
        f"unexpected: {sorted(set(sco7048_stops['sgrna']) - set(curated['sgrna']))}"
    )
    assert len(merged) == len(curated)

    wrong = {
        row.sgrna: (_split(row.mutations_teemi), _split(row.mutations_curated))
        for row in merged.itertuples()
        if _split(row.mutations_teemi) != _split(row.mutations_curated)
    }
    assert not wrong, f"mutation calls disagree with the curation: {wrong}"


# ---------------------------------------------------------------------------
# Base editing against CRISPy-web on a PLUS-strand gene
# ---------------------------------------------------------------------------


def test_base_editing_matches_crispy_web_on_a_plus_strand_gene(sco5087_stops):
    expected = (
        pd.read_csv(f"{CRISPY_DIR}/BEST_SCO5087_only_stop_codons.csv", index_col=0)[
            ["Sequence", "C to T mutations"]
        ]
        .rename(columns={"Sequence": "sgrna", "C to T mutations": "mutations"})
        .sort_values("sgrna")
        .reset_index(drop=True)
    )
    actual = (
        sco5087_stops[["sgrna", "mutations"]]
        .sort_values("sgrna")
        .reset_index(drop=True)
    )

    assert list(actual["sgrna"]) == list(expected["sgrna"])
    for sgrna, got, want in zip(
        actual["sgrna"], actual["mutations"], expected["mutations"]
    ):
        assert _split(got) == _split(want), f"guide {sgrna}"


# ---------------------------------------------------------------------------
# Cas9 guide finding against CRISPy-web
# ---------------------------------------------------------------------------


def test_cas9_finds_every_guide_crispy_web_reports(coelicolor):
    """Every CRISPy-web guide must be found, with the same PAM.

    teemi reports a superset: 289 guides against CRISPy-web's 206 for this
    gene, and the 83 extras are not explained by GC, by position in the CDS or
    by duplicates, so they are most likely filters CRISPy-web applies by
    default and teemi does not. That is why this asserts containment rather
    than equality -- it still fails if the finder loses guides or moves a PAM.
    """
    expected = pd.read_csv(f"{CRISPY_DIR}/cas9_SCO5087_0_mismatches.csv", index_col=0)
    actual = _guides(coelicolor, "SCO5087")

    assert actual["sgrna"].is_unique
    assert not set(expected["Sequence"]) - set(actual["sgrna"])

    merged = actual[["sgrna", "pam"]].merge(
        expected[["Sequence", "PAM"]].rename(columns={"Sequence": "sgrna"}), on="sgrna"
    )
    assert len(merged) == len(expected)
    mismatched = merged[merged["pam"] != merged["PAM"]]
    assert mismatched.empty, f"PAM differs for {list(mismatched['sgrna'])}"


@pytest.mark.parametrize(
    "gene,locus_tag,record,export,offset,total",
    GENSCRIPT_CAS12A,
    ids=[row[0] for row in GENSCRIPT_CAS12A],
)
def test_cas12a_finds_every_guide_genscript_reports(
    gene, locus_tag, record, export, offset, total
):
    """Every GenScript guide must be found, with the same PAM, in the same place.

    This is the Cas12a counterpart of the CRISPy-web test above, and the only
    reference for the nuclease that teemi gets wholly from outside. It is worth
    more than a second Cas9 reference would be, because GenScript encodes the
    PAM and the spacer length itself: agreeing with it checks the two constants
    in ``_NUCLEASES["cas12a"]``, which no test written against teemi's own table
    can do.

    As with CRISPy-web, teemi reports a superset -- 81 guides to GenScript's 26
    here, 67 to 29 for LYP1 -- because GenScript ranks and truncates by its own
    scores, so this asserts containment rather than equality.

    The position assertion is the part the sequence comparison cannot stand in
    for. ``sgrna_loc`` is in coordinates of the CDS as extracted 5'->3'; mapping
    it back out to the chromosome and landing on the base GenScript names ties
    teemi's internal convention to an externally stated genomic coordinate.
    Both position formulas are covered because CAN1 is on the minus strand and
    LYP1 on the plus.
    """
    rec = SeqIO.read(f"{GENSCRIPT_DIR}/{record}", "genbank")
    cds = next(f for f in rec.features if f.type == "CDS")
    assert cds.qualifiers["locus_tag"] == [locus_tag]

    molecule = Dseqrecord(rec)
    molecule.id = gene
    molecule.features = rec.features
    guides = extract_sgRNAs(
        SgRNAargs(
            dseqrecord=molecule,
            locus_tag=[locus_tag],
            cas_type="cas12a",
            step=["find"],
            off_target_seed=13,
        )
    )
    expected = pd.read_csv(f"{GENSCRIPT_DIR}/{export}")

    assert len(guides) == total
    assert guides["sgrna"].is_unique
    found = {row.sgrna: row for row in guides.itertuples()}
    missing = sorted(set(expected["Sequence"]) - set(found))
    assert (
        not missing
    ), f"{gene}: teemi misses {len(missing)} guides, e.g. {missing[:3]}"

    # The spacer length is one of the two constants this reference checks.
    assert set(guides["sgrna"].str.len()) == {_NUCLEASES["cas12a"].protospacer_len}

    start, end = int(cds.location.start), int(cds.location.end)
    on_plus = cds.location.strand == 1
    for row in expected.itertuples():
        hit = found[row.Sequence]
        assert hit.pam == row.PAM, f"{gene} {row.Sequence}: {hit.pam} != {row.PAM}"

        # GenScript's span runs PAM then protospacer along its own strand, so
        # the protospacer begins 4 nt in on the plus strand and at the span's
        # start on the minus.
        span_start = int(row.Location.split("..")[0])
        want = (
            span_start + _NUCLEASES["cas12a"].pam_len
            if row.Strand == "+"
            else span_start
        )

        begin = hit.sgrna_loc - len(hit.sgrna)
        in_record = start + begin if on_plus else end - hit.sgrna_loc
        assert in_record + offset + 1 == want, (
            f"{gene} {row.Sequence}: teemi puts it at "
            f"{in_record + offset + 1}, GenScript at {want}"
        )


# ---------------------------------------------------------------------------
# Every reported guide must exist on the molecule
# ---------------------------------------------------------------------------


def _full_site(row, cas_type):
    """The protospacer and PAM in the order they occur on the DNA."""
    return row.sgrna + row.pam if cas_type == "cas9" else row.pam + row.sgrna


@pytest.mark.parametrize("cas_type", ["cas9", "cas12a", "cas3"])
def test_guides_on_a_compound_cds_exist_on_the_molecule(coelicolor, cas_type):
    """A guide the molecule does not carry cannot be targeted.

    ``feature.extract()`` concatenates the parts of a compound location, so a
    window spanning two of them reads a sequence that exists only in the
    concatenation. SCO7048 and SCO5087 cannot catch this -- both are one-part
    CDSs -- but SCO0022 is ``join(22156..22299, 22290..22430, 22624..>22761)``,
    whose first two parts overlap by 10 bp and whose second and third are 193 bp
    apart. It is the one compound CDS in this genome, and before the finders
    skipped these windows it yielded 10 cas9 and 2 cas3 guides that occur
    nowhere in the chromosome on either strand.

    This is the assertion no coordinate convention can substitute for: it reads
    the sequence back out of the record and is therefore independent of how the
    annotation describes the gene.
    """
    guides = extract_sgRNAs(
        SgRNAargs(
            dseqrecord=coelicolor,
            locus_tag=["SCO0022"],
            cas_type=cas_type,
            step=["find"],
            off_target_seed=13,
        )
    )
    assert not guides.empty, f"no {cas_type} guides found for SCO0022"

    genome = str(coelicolor.seq)
    reverse = str(Seq(genome).reverse_complement())
    absent = [
        _full_site(row, cas_type)
        for row in guides.itertuples()
        if _full_site(row, cas_type) not in genome
        and _full_site(row, cas_type) not in reverse
    ]
    assert not absent, (
        f"{len(absent)} of {len(guides)} {cas_type} guides are absent from the "
        f"chromosome on both strands: {absent}"
    )


def test_a_cds_crossing_the_origin_of_a_circular_record_keeps_its_guides():
    """Adjacent parts are not an intron, even though the location is compound.

    A CDS running across the origin of a circular molecule is a two-part
    location whose parts *are* contiguous on the DNA, so ``extract()``'s
    concatenation is the real sequence and every window in it is targetable.
    Skipping windows per part rather than per contiguous run would silently
    drop them, which is why ``_scan_breaks`` tests for adjacency instead of
    assuming a compound location means splicing.
    """
    length = 400
    bases = list("".join("ACGT"[(i * 7 + i // 4) % 4] for i in range(length)))
    for position in (0, 1, 4, 5, 10, 11):
        bases[position] = "G"  # put NGG PAMs just past the origin
    sequence = "".join(bases)

    feature = SeqFeature(
        CompoundLocation([SimpleLocation(379, length, 1), SimpleLocation(0, 40, 1)]),
        type="CDS",
        qualifiers={"locus_tag": ["WRAP1"]},
    )

    def guides_for(circular):
        record = Dseqrecord(sequence, circular=circular)
        record.features = [feature]
        return extract_sgRNAs(
            SgRNAargs(
                dseqrecord=record,
                locus_tag=["WRAP1"],
                cas_type="cas9",
                step=["find"],
                off_target_seed=13,
            )
        )

    circular = guides_for(True)
    twice = sequence + sequence
    for row in circular.itertuples():
        site = row.sgrna + row.pam
        assert (
            site in twice or str(Seq(site).reverse_complement()) in twice
        ), f"{site} was reported but does not occur on the circular molecule"

    # The same parts on a linear record genuinely are not contiguous, so almost
    # every window really is unreadable there -- the guard is reading topology,
    # not just counting parts.
    assert len(circular) > len(guides_for(False))


def _spliced_record():
    """A two-exon CDS whose only cas12a PAM sits where the window must cross.

    The filler is a fixed-seed LCG rather than a literal so the sequence is both
    reproducible and non-repetitive: a periodic filler makes the chimeric window
    occur elsewhere in the record by chance, which would make the assertion below
    pass for the wrong reason.
    """
    exon, intron = 40, 60

    bases, state = [], 12345
    for _ in range(exon + intron + exon):
        state = (1103515245 * state + 12345) % (2**31)
        bases.append("ACGT"[(state >> 16) & 3])
    # Scrub accidental TTT runs so the only forward cas12a PAM is the placed one.
    for i in range(len(bases) - 2):
        if bases[i] == bases[i + 1] == bases[i + 2] == "T":
            bases[i + 2] = "A"
    # A PAM at CDS offset 30 puts the 27 nt window at [30, 57), across the
    # junction at 40. Cas12a reads its protospacer 3' of a 4 nt PAM, and a
    # default SgRNAargs now resolves protospacer_len to the As/LbCas12a 23.
    bases[30:34] = list("TTTA")
    sequence = "".join(bases)

    feature = SeqFeature(
        CompoundLocation(
            [
                SimpleLocation(0, exon, 1),
                SimpleLocation(exon + intron, exon + intron + exon, 1),
            ]
        ),
        type="CDS",
        qualifiers={"locus_tag": ["SPLICED1"]},
    )
    record = Dseqrecord(sequence, circular=False)
    record.features = [feature]
    return record, feature, sequence, exon


def test_a_cas12a_window_spanning_an_intron_is_not_reported():
    """The Cas12a guard needs its own case: SCO0022 never exercises it.

    SCO0022 is the only compound CDS in the Streptomyces fixture and no Cas12a
    window happens to cross one of its boundaries, so the parametrised test
    above passes without the Cas12a skip ever running. A wrong window length
    there would go unnoticed, hence this record, built so that a crossing window
    must exist.
    """
    record, feature, sequence, exon = _spliced_record()
    coding = str(feature.extract(Seq(sequence)))
    reverse = str(Seq(sequence).reverse_complement())

    # The window the finder must refuse, and the premise that makes it a test:
    # it is a chimera, present nowhere on the molecule. Its length has to track
    # the protospacer the finder actually scans for, or the assertion below
    # compares a 20-mer against 23-mers and passes without testing anything.
    window = 4 + _NUCLEASES["cas12a"].protospacer_len
    chimera = coding[30 : 30 + window]
    assert len(chimera) == 27
    assert chimera.startswith("TTTA")
    assert chimera not in sequence and chimera not in reverse

    guides = extract_sgRNAs(
        SgRNAargs(
            dseqrecord=record,
            locus_tag=["SPLICED1"],
            cas_type="cas12a",
            step=["find"],
            off_target_seed=13,
        )
    )
    assert chimera[4:] not in set(guides["sgrna"])
    for row in guides.itertuples():
        site = row.pam + row.sgrna
        assert site in sequence or site in reverse, f"{site} is not on the molecule"


def test_parts_on_opposite_strands_are_a_break():
    """A strand change between parts is a discontinuity, whatever the coordinates.

    ``extract()`` reverse-complements each part according to its own strand, so
    the join of a plus part and a minus part is not a walk along either strand
    and no window may span it. The branch is unreachable from the repository's
    fixtures, so it is asserted directly.
    """
    record = Dseqrecord("A" * 200, circular=False)
    feature = SeqFeature(
        CompoundLocation([SimpleLocation(0, 40, 1), SimpleLocation(40, 80, -1)]),
        type="CDS",
        qualifiers={"locus_tag": ["TRANS1"]},
    )
    coding = str(feature.extract(Seq(str(record.seq))))

    breaks = _scan_breaks(feature, record, coding)
    assert breaks[-1] == [40]
    assert breaks[1] == [len(coding) - 40]


# Guides whose window sits flush against one of SCO0022's genomic breaks: the
# first begins on the first base of a part, the other two end on the last base
# of one. All three are wholly inside a single part and present on the
# chromosome, so all three must survive the skip.
SCO0022_BOUNDARY_GUIDES = [
    ("GAGCGAACGCCGCACTGCTC", "GGG"),  # window starts at a break
    ("AACCTGGCCCCACTCGCTGA", "AGG"),  # window ends at a break
    ("GTTCTCGTGGCCGACGAGGT", "CGG"),  # window ends at a break
]


def test_the_skip_does_not_discard_guides_that_merely_touch_a_break(coelicolor):
    """Pin both ends of the crossing test, which are otherwise unpinned.

    ``test_guides_on_a_compound_cds_exist_on_the_molecule`` is one-sided: it
    asserts that every guide reported is real, never that a real guide was not
    thrown away. Loosening either bound of ``match.start() < offset <
    window_end`` to ``<=`` therefore passes the whole suite while silently
    dropping targetable guides -- 68 to 67 for the lower bound and 68 to 66 for
    the upper. A window that merely touches a break at one end lies entirely
    within one part, and only a strict comparison keeps it.
    """
    guides = extract_sgRNAs(
        SgRNAargs(
            dseqrecord=coelicolor,
            locus_tag=["SCO0022"],
            cas_type="cas9",
            step=["find"],
            off_target_seed=13,
        )
    )

    found = set(zip(guides["sgrna"], guides["pam"]))
    missing = [g for g in SCO0022_BOUNDARY_GUIDES if g not in found]
    assert not missing, (
        f"guides flush against a break were discarded: {missing}; the crossing "
        "test must compare strictly at both ends"
    )
    # The count pins the other direction: a skip that fires too rarely would
    # readmit the chimeras the strand-aware test above rejects.
    assert len(guides) == 68, f"expected 68 cas9 guides for SCO0022, got {len(guides)}"


POSITION_LOCI = ["SCO5087", "SCO7048", "SCO0022"]


@pytest.mark.parametrize("cas_type", ["cas9", "cas12a", "cas3"])
def test_sgrna_loc_locates_the_guide_in_the_coding_sequence(coelicolor, cas_type):
    """``sgrna_loc`` has to say where the guide actually is.

    The column is one past the protospacer's end in coordinates of the CDS as
    extracted 5'->3', which is how ``crispr_best.mutate_sequence`` finds a
    guide: it slices ``[sgrna_loc - len(sgrna) : sgrna_loc]``. A value out by
    even one base therefore edits the wrong residue.

    Reading that slice back out of the CDS checks the number against the
    sequence rather than against the formula that produced it. Nothing checked
    the cas12a or cas3 positions that way before, and both were wrong while
    every cas9 position was right.

    The two scans are covered by requiring both readings to occur: a guide
    found on the coding scan of a 5' PAM lies along the CDS, and one found on
    the reverse-complement scan lies against it, so a test that saw only one
    orientation would be exercising only one of the two position formulas.
    """
    guides = extract_sgRNAs(
        SgRNAargs(
            dseqrecord=coelicolor,
            locus_tag=POSITION_LOCI,
            cas_type=cas_type,
            step=["find"],
            off_target_seed=13,
        )
    )
    assert not guides.empty, f"no {cas_type} guides found"

    coding = _gene_sequences(coelicolor, POSITION_LOCI)
    along = against = 0
    misplaced = []
    for row in guides.itertuples():
        cds = coding[row.locus_tag]
        start = row.sgrna_loc - len(row.sgrna)
        found = cds[start : row.sgrna_loc] if start >= 0 else ""
        if row.sgrna == found:
            along += 1
        elif row.sgrna == str(Seq(found).reverse_complement()):
            against += 1
        else:
            misplaced.append((row.locus_tag, row.sgrna_loc, row.sgrna))

    assert not misplaced, (
        f"{len(misplaced)} of {len(guides)} {cas_type} guides are not at the "
        f"sgrna_loc they report, e.g. {misplaced[:3]}"
    )
    assert along and against, (
        f"{cas_type} guides were all read in one orientation "
        f"(along={along}, against={against}), so only one of the two "
        f"position formulas was exercised"
    )


# ---------------------------------------------------------------------------
# Cas3 against the notebook published with pCRISPR-Cas3
# ---------------------------------------------------------------------------


def _type1_footprint(row):
    """0-based half-open span of PAM+protospacer for one reference row.

    The notebook writes its coordinates 1-based on the forward strand
    (``hits["Start"].append(match.span()[1]+1)``) and 0-based on the reverse
    (``hits["Start"].append(match.span()[0]-spacer_length)``). That asymmetry
    is theirs; normalising it here is what lets one number cover both strands.
    """
    pam_len = _NUCLEASES["cas3"].pam_len
    if row.Strand == "Forward":
        return row.Start - 1 - pam_len, row.End
    return row.Start, row.End + pam_len


def _type1_expected(expected, genome, lo, hi, locus_tag):
    """Read every reference row back out of the chromosome.

    Returns the set of ``(0-based protospacer start, spacer, strand)`` the
    notebook reports, with its ``Forward``/``Reverse`` read as the chromosome's
    plus/minus. The fixture was produced by transcribing the notebook's search
    loop rather than by executing the notebook, so each row is re-derived here
    from the genome: a mis-transcribed offset has to show up as a window that
    does not carry TTC, or whose slice is not the spacer, instead of quietly
    agreeing with teemi.
    """
    pam_len = _NUCLEASES["cas3"].pam_len
    want = set()
    for row in expected.itertuples():
        start, stop = _type1_footprint(row)
        assert (
            lo <= start and stop <= hi
        ), f"{locus_tag}: reference window {start}-{stop} leaves the CDS"
        site = genome[start:stop]
        if row.Strand == "Forward":
            assert (
                site[:pam_len] == "TTC"
            ), f"{locus_tag}: {site[:pam_len]} is not a PAM"
            assert site[pam_len:] == row.Spacer
            want.add((start + pam_len, row.Spacer, 1))
        else:
            assert (
                site[-pam_len:] == "GAA"
            ), f"{locus_tag}: {site[-pam_len:]} is not a PAM"
            assert str(Seq(site[:-pam_len]).reverse_complement()) == row.Spacer
            want.add((start, row.Spacer, -1))
    return want


def _type1_cds(record, locus_tag):
    return next(
        f
        for f in record.features
        if f.type == "CDS" and f.qualifiers.get("locus_tag", [""])[0] == locus_tag
    )


@pytest.mark.parametrize("locus_tag, export, total", TYPE1_CAS3)
def test_cas3_finds_exactly_the_spacers_the_reference_notebook_reports(
    coelicolor, locus_tag, export, total
):
    """Cas3 had no external reference at all until this one.

    CRISPy-web is Cas9 and, as of version 3, Cas9/CRISPRi/TnpB -- still no
    Type I. GenScript has no Type I-C either. The only implementation of this
    system is the one its authors published alongside the plasmid, and it is
    the same source teemi takes ``TTC`` and ``34`` from, so agreeing with it
    checks both constants in ``_NUCLEASES["cas3"]`` against the code that
    designed the spacers the paper actually cloned.

    What it does not do is make the biology independent: it is the same lab as
    teemi and the same lab as CRISPy-web. What it pins independently is the
    mechanics -- which side of the protospacer the PAM sits on, how far the
    spacer runs, how the reverse strand is reached, and where on the
    chromosome the result lands. Those are exactly where the cas12a and cas3
    positions were wrong before.

    Equality holds here rather than containment: the notebook ranks but does
    not truncate, so for a given region it reports every site, as teemi does.
    Both genes are covered because SCO5087 is on the plus strand and SCO7048
    on the minus, and ``sgrna_strand`` is read against the chromosome rather
    than against the gene -- for SCO7048 the 20 guides teemi puts on the minus
    strand are the notebook's 20 ``Reverse`` rows, not its 11 ``Forward`` ones.
    The strand is part of the compared key for that reason: with only the
    counts checked, flipping ``_NUCLEASES["cas3"].strand_sign`` still leaves
    both signs present and the test passes while every guide is mislabelled.
    """
    expected = pd.read_csv(f"{TYPE1_DIR}/{export}")
    cds = _type1_cds(coelicolor, locus_tag)
    lo, hi = int(cds.location.start), int(cds.location.end)
    want = _type1_expected(expected, str(coelicolor.seq).upper(), lo, hi, locus_tag)

    guides = extract_sgRNAs(
        SgRNAargs(
            dseqrecord=coelicolor,
            locus_tag=[locus_tag],
            cas_type="cas3",
            step=["find"],
            off_target_seed=13,
        )
    )
    assert len(guides) == total
    assert guides["sgrna"].is_unique

    on_plus = cds.location.strand == 1
    got = {
        (
            lo + row.sgrna_loc - len(row.sgrna) if on_plus else hi - row.sgrna_loc,
            row.sgrna,
            int(row.sgrna_strand),
        )
        for row in guides.itertuples()
    }
    assert got == want, (
        f"{locus_tag}: {len(want - got)} spacers only in the reference, "
        f"{len(got - want)} only in teemi; e.g. {sorted(want - got)[:2]} "
        f"vs {sorted(got - want)[:2]}"
    )

    # The two constants this reference is here to check.
    assert set(guides["pam"]) == {"TTC"}
    assert set(guides["sgrna"].str.len()) == {_NUCLEASES["cas3"].protospacer_len}
    # Both formulas have to be exercised, or half of the comparison above is
    # vacuous. Which guide is on which strand is pinned by ``got == want``.
    assert set(guides["sgrna_strand"]) == {1, -1}


@pytest.mark.parametrize("locus_tag, export, total", TYPE1_CAS3)
def test_cas3_off_target_count_reproduces_the_reference_seed_metric(
    coelicolor, locus_tag, export, total
):
    """``off_target_seed=8`` is the notebook's own metric, to the hit.

    The module's notes claim that 8 reproduces what the paper describes as
    "the seed off-target search considers the PAM and the PAM proximal first
    8 nt", counted "only if a PAM site is present". This is that claim as a
    test rather than as prose. The notebook builds one string of
    ``seed_length+3`` bases -- PAM included -- and counts it with
    ``count_overlap`` on both strands; teemi matches the PAM separately and
    counts the 8 bases after it, over a sequence list that already holds both
    strands, which comes to the same set of sites.

    The one difference is the on-target site. teemi subtracts it
    (``max(off_target_counter[sgrna_seed] - 1, 0)``) while the notebook leaves
    it in, which is why its README says "a match count of 1 indicates that no
    off-target sites were found". So the reference count is teemi's plus one,
    exactly, and that offset is the assertion: were teemi counting one strand,
    or anchoring the seed in the wrong register, the two columns would diverge
    by far more than a constant.
    """
    expected = pd.read_csv(f"{TYPE1_DIR}/{export}")
    assert expected["Spacer"].is_unique
    reference = dict(zip(expected["Spacer"], expected["Seed 0 nt matches"]))

    guides = extract_sgRNAs(
        SgRNAargs(
            dseqrecord=coelicolor,
            locus_tag=[locus_tag],
            cas_type="cas3",
            step=["find"],
            off_target_seed=TYPE1_SEED,
        )
    )
    assert len(guides) == total
    assert set(guides["sgrna_seed_sequence"].str.len()) == {TYPE1_SEED}

    wrong = [
        (row.sgrna, row.off_target_count, reference[row.sgrna])
        for row in guides.itertuples()
        if row.off_target_count != reference[row.sgrna] - 1
    ]
    assert not wrong, (
        f"{locus_tag}: {len(wrong)} of {len(guides)} seed counts are not the "
        f"reference count minus one, e.g. {wrong[:3]}"
    )


def test_cas3_finds_the_spacers_pcrispr_cas3_was_cloned_with(coelicolor):
    """The guides that actually cut, rather than guides a scan proposes.

    Every other reference here is a list of candidates. These ten are the
    protospacers the pCRISPR-Cas3 paper reports building into plasmids and
    conjugating into *Streptomyces*, taken from the supplementary overview
    workbook of its data record (DTU Data, doi:10.11583/DTU.22786157). They
    are not merely proposed, then: each was ordered and cloned. Not each
    produced a clean deletion -- the paper calls one of the region 21-23
    outcomes "fully unspecific" -- so what is asserted here is that teemi
    finds the guides that were actually made, not that every one worked.

    Three of the workbook's plasmid identifiers appear in the paper's own text
    and tie to the same targets ("Targeting of regions 1-3 using p245",
    "targeting regions 21-23 with plasmid p246", and p253 for the
    undecylprodigiosin BGC), which is what corroborates the workbook against
    the publication. Each row's spacer and coordinate are re-derived from this
    chromosome below regardless.

    It cannot validate ``TTC`` or ``34`` and is not written as though it does.
    The paper chose these spacers by searching for that PAM at that length, so
    conformance is its design rule restated, not evidence for it; the
    independent source for the PAM is Leon et al. 2021, cited beside the
    parameter table. What it does validate is the scanner, over a far wider
    stretch of the chromosome than two genes reach -- 92 kb to 7.4 Mb, six
    genes, four protospacers on the minus strand and three on the plus.

    Seven of the ten lie inside a single-part CDS and are therefore reachable
    by a locus tag. The other three are intergenic: Cas3 deletes processively
    from its target, so the paper also aims it between genes, which
    ``extract_sgRNAs`` cannot do because it takes locus tags and nothing else.
    That split is asserted rather than described, so re-annotation cannot
    quietly shrink what this test covers.
    """
    expected = pd.read_csv(f"{TYPE1_DIR}/{TYPE1_CLONED}", keep_default_na=False)
    genome = str(coelicolor.seq).upper()
    pam_len = _NUCLEASES["cas3"].pam_len
    spacer_len = _NUCLEASES["cas3"].protospacer_len

    # Re-derive every published spacer from this chromosome before comparing
    # anything, so a wrong coordinate in the fixture cannot pass as a match.
    for row in expected.itertuples():
        begin = row.Start - 1
        end = begin + spacer_len
        assert len(row.Spacer) == spacer_len, f"{row.Name} is not {spacer_len} nt"
        assert row.End == end, f"{row.Name}: End disagrees with Start + length"
        if row.Strand == 1:
            assert genome[begin - pam_len : begin] == "TTC", f"{row.Name}: no PAM"
            assert genome[begin:end] == row.Spacer, f"{row.Name}: not at Start"
        else:
            assert genome[end : end + pam_len] == "GAA", f"{row.Name}: no PAM"
            assert (
                str(Seq(genome[begin:end]).reverse_complement()) == row.Spacer
            ), f"{row.Name}: not at Start"

    assert len(expected) == 10
    reachable = expected[expected.LocusTag != ""]
    assert len(reachable) == 7, "the CDS/intergenic split has moved"

    spans = {
        f.qualifiers["locus_tag"][0]: f
        for f in coelicolor.features
        if f.type == "CDS"
        and f.qualifiers.get("locus_tag", [""])[0] in set(reachable.LocusTag)
    }
    guides = extract_sgRNAs(
        SgRNAargs(
            dseqrecord=coelicolor,
            locus_tag=sorted(spans),
            cas_type="cas3",
            step=["find"],
            off_target_seed=13,
        )
    )
    found = set()
    for row in guides.itertuples():
        cds = spans[row.locus_tag]
        lo, hi = int(cds.location.start), int(cds.location.end)
        begin = (
            lo + row.sgrna_loc - len(row.sgrna)
            if cds.location.strand == 1
            else hi - row.sgrna_loc
        )
        found.add((begin + 1, row.sgrna, int(row.sgrna_strand)))

    missing = [
        row.Name
        for row in reachable.itertuples()
        if (row.Start, row.Spacer, row.Strand) not in found
    ]
    assert not missing, (
        f"teemi does not report {len(missing)} of {len(reachable)} published "
        f"pCRISPR-Cas3 spacers at the coordinate and strand the paper gives: "
        f"{missing}"
    )
