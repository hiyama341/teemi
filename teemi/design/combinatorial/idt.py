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

"""The order sheet as something a vendor's form will take.

:func:`~teemi.design.combinatorial.report.primer_dataframe` is an order
*sheet*: seventeen columns of Tm, dG and ``flagged_by``, which is what decides
whether to order.  It is not an order *form*.  Turning one into the other by
hand means copying two columns into a web form, and that is the single step in
the whole chain where a typo becomes a tube of the wrong DNA.

IDT takes bulk oligos in two shapes and this module builds both:

* **tubes** -- ``Name, Sequence, Scale, Purification``, pasted into the bulk
  entry form.  The scale and purification codes are IDT's own, and the scale
  is PER ROW because IDT's scales are length-limited: see :data:`IDT_SCALES`.
* **a plate** -- ``Well Position, Name, Sequence``, uploaded as a workbook with
  one worksheet per plate.  Scale and purification are absent because they are
  chosen once for the whole plate in IDT's form, not per well.

Everything here is PURE.  :mod:`teemi.design.combinatorial.export` is the only
file-writing edge of the package and it stays that way: these return frames and
strings, and ``export`` writes them.

Both shapes read ``result.primers`` in its own order, so neither is a second
computation and the two files cannot disagree about what is in the order.  That
order matters more than it looks: see :func:`plate_frames`.
"""

from __future__ import annotations

import re
from collections import Counter
from typing import TYPE_CHECKING, List, Optional, Sequence, Tuple

if TYPE_CHECKING:  # pragma: no cover - import cycle at runtime, types only
    import pandas as pd

__all__ = [
    "IDT_PLATE_COLUMNS",
    "IDT_TUBE_COLUMNS",
    "IDT_PLATE_MINIMUM",
    "IDT_SCALES",
    "IDT_STANDARD_MAX",
    "DEFAULT_SCALE",
    "DEFAULT_PURIFICATION",
    "well_positions",
    "oligo_name",
    "scale_for_length",
    "scale_mix",
    "tube_frame",
    "plate_frames",
    "primers_fasta",
    # Read from export.py, so it is part of the surface whether or not it
    # looks like a helper.
    "plate_count",
]

#: IDT's plate upload expects exactly these three headers, in this order.
#:
#: FLAT ``Sequence``, NOT ``sequence(5'-3')`` like every other column of DNA in
#: a bundle.  IDT reads this header, so renaming it to say which end is which
#: would be a rejected upload.  The orientation is stated in the bundle's
#: ``00_LAYOUT.txt`` line instead, which is the one place that can carry it.
IDT_PLATE_COLUMNS = ("Well Position", "Name", "Sequence")

#: IDT's bulk tube entry.  ``Scale`` and ``Purification`` take IDT's own codes.
IDT_TUBE_COLUMNS = ("Name", "Sequence", "Scale", "Purification")

#: A 96-well plate, lettered rows and numbered columns.
PLATE_ROWS = 8
PLATE_COLS = 12

#: Oligos IDT requires before it will make a 96-well plate at all.
#:
#: Below this the order has to go as tubes, so :func:`plate_frames` still lays
#: the plate out -- the layout is not wrong, it is merely unorderable -- and
#: ``build_summary`` says so, with the count beside the limit.
#:
#: NOT CITED, like :data:`IDT_NAME_MAX` below it, and for the same reason: this
#: is read off IDT's plate-ordering page rather than a document with a version
#: on it, so it is a number that can move without anything here noticing.  It
#: is a floor, so being wrong about it is one-sided -- too high turns an
#: orderable plate into a warning, too low lets an unorderable one through
#: silently -- and 24 is the figure the order form itself enforces.  The 96
#: and 48 of 384-well plates and HOTPlates were quoted here as though equally
#: certain and are now left out, since neither shape is offered.
IDT_PLATE_MINIMUM = 24

#: Longest oligo name written into an order.
#:
#: A GUESS, AND KNOWN TO BE ONE.  IDT publishes no limit on the length or the
#: characters of an oligo name; 50 is what a third-party format converter
#: claims for their pooled-oligo format, and it is short enough that no
#: reasonable form will refuse it.  A real rejected order would be better
#: evidence than this constant.  :func:`oligo_name` is built so that truncating
#: here cannot collide.
IDT_NAME_MAX = 50

#: IDT's standard DNA oligo scales: IDT's own ordering code, then the shortest
#: and the longest oligo that scale will synthesise, smallest scale first.
#:
#: READ OFF IDT'S OWN ORDER PAGE on 2026-11-25 --
#: https://eu.idtdna.com/pages/products/custom-dna-rna/dna-oligos/custom-dna-oligos
#: -- where the same three ranges are given twice, once in the tube price table
#: and again in the yield guarantees below it.
#:
#: CITED, unlike :data:`IDT_PLATE_MINIMUM` and :data:`IDT_NAME_MAX` above, and
#: that is not inconsistency: a plate minimum is a floor, so being wrong about
#: it is one-sided and cheap.  These are not.  Too high an upper bound writes a
#: scale IDT will not fill and the order comes back as a question days later;
#: too low buys a dearer scale than the oligo needed.  A 25 nmol row for a
#: 65-mer is exactly that first failure, and it is what this table exists to
#: stop -- the three shipped fixtures all top out at 59 nt, so nothing here
#: noticed for as long as that was the only DNA anyone tried.
#:
#: THE MINIMUM IS PART OF THE RANGE, not decoration.  25 nmol starts at 15
#: bases, so a 12-mer is unorderable at the smallest scale despite sitting far
#: under its 60-base cap, and the scale that takes it is the next one up.
#:
#: 1, 2, 5 and 10 umole are LEFT OUT.  They synthesise the same 5-100 bases as
#: 250 nmol at several times the price, so no length can reach for one and an
#: entry each would be a row :func:`scale_for_length` could never choose.  They
#: buy more DNA, not longer DNA, and how much DNA is the caller's to ask for:
#: :func:`tube_frame` takes a ``scale`` for that.
IDT_SCALES = (
    ("25nm", 15, 60),
    ("100nm", 10, 90),
    ("250nm", 5, 100),
)

#: The longest oligo any standard scale takes.  Past this is Ultramer
#: territory, which this module does not order -- see :func:`scale_for_length`.
IDT_STANDARD_MAX = max(high for _, _, high in IDT_SCALES)

#: The shortest oligo any standard scale takes.
IDT_STANDARD_MIN = min(low for _, low, _ in IDT_SCALES)

#: IDT's codes for the smallest synthesis scale and for no purification.
#:
#: The smallest scale is no longer what goes on every row -- that was the bug;
#: :func:`scale_for_length` decides per oligo -- so this is the foot of the
#: ladder, read off it rather than written again so the two cannot drift.
DEFAULT_SCALE = IDT_SCALES[0][0]
DEFAULT_PURIFICATION = "STD"

#: Characters IDT's own sample templates use.  Anything else is folded to an
#: underscore rather than risking an order refused for a reason the person who
#: placed it cannot see -- a primer id can pick up a ``+`` from a shared
#: molecule's name, and an uploaded part's name can carry anything at all.
#:
#: A RUN of unsafe characters becomes ONE underscore, but runs of underscores
#: are then left alone, which is the opposite of what a sanitiser usually does.
#: ``ID_SEPARATOR`` is ``__``: a primer id is ``<left>__<right>_<direction>``,
#: so collapsing doubles would turn two different junction names into one and
#: make ``A__B`` and ``A_B`` the same oligo on an order.  The cost is that
#: ``sgRNA_#1`` comes out as ``sgRNA__1`` rather than ``sgRNA_1``, which no
#: form minds.
_UNSAFE_NAME_CHARS = re.compile(r"[^A-Za-z0-9_.+-]+")


def well_positions(rows: int = PLATE_ROWS, cols: int = PLATE_COLS) -> Tuple[str, ...]:
    """One plate's wells in fill order: ``A1, B1 ... H1, A2 ... H12``.

    COLUMN-MAJOR, down the columns rather than across the rows, because an
    8-channel pipette spans exactly one column.  Filled that way a part-full
    plate's oligos sit in whole columns and the empties are whole columns too,
    so every column that holds anything can be worked in one draw.  Filled
    across the rows the same oligos are smeared over every row of the plate and
    not one column is complete.

    IDT reads the ``Well Position`` column, so this decides only where the
    oligos physically land in the plate that arrives.

    ``teemi.legacy.build.robot_assembly.well_keys_96(row=False)`` computes the
    same list.  It is not reused: it is fixed at 96 wells, it is in the legacy
    tree, and importing from there would tie the newest output in the package
    to the oldest module in it.  The two must not disagree, which is cheap to
    keep true -- a plate is lettered from A and numbered from 1.
    """
    if not 1 <= rows <= 26:
        raise ValueError(f"rows must be between 1 and 26, got {rows}")
    if cols < 1:
        raise ValueError(f"cols must be at least 1, got {cols}")
    return tuple(
        f"{chr(ord('A') + row)}{col}"
        for col in range(1, cols + 1)
        for row in range(rows)
    )


def oligo_name(primer) -> str:
    """What to call one oligo on an order: its label, then its id.

    THE LABEL COMES FIRST, and that is what makes the length cap safe.  A
    primer id is as long as the junction it serves and a shared molecule's id
    is the join of every name that asked for it, so ids routinely run past any
    vendor's field -- and two of them truncated to the same width can truncate
    to the same string.  A label is unique within one plan by construction, so
    leading with it keeps every name in an order distinct however hard the tail
    is cut.

    The label alone would not do.  It is positional -- insert a gene and
    everything after it renumbers -- so a tube labelled only ``007`` cannot be
    checked against anything once it is out of the box the order came in.  The
    id says what the molecule is for.
    """
    label = getattr(primer, "label", "") or ""
    tail = _UNSAFE_NAME_CHARS.sub("_", str(getattr(primer, "id", "") or "")).strip("_")
    name = f"{label}_{tail}" if label and tail else (label or tail)
    return name[:IDT_NAME_MAX].rstrip("_") or "unnamed_oligo"


def _rows(result) -> List[Tuple[str, str]]:
    """``(name, sequence)`` per distinct molecule, in order-sheet order.

    The one place either shape reads the result, so the tube sheet and the
    plate cannot come to different conclusions about what is in the order.
    """
    return [(oligo_name(p), str(p.sequence).upper()) for p in result.primers]


def scale_for_length(length: int) -> str:
    """The cheapest standard scale IDT will make a ``length``-mer at.

    SMALLEST FIRST, and the first range that *contains* the length wins rather
    than the first whose cap clears it.  The scales overlap -- 100 nmol starts
    below 25 nmol and reaches further -- and the cheapest is the one with the
    narrowest range, so walking up the ladder stops at the least DNA that can
    be bought for the length asked.  Checking only the cap would put a 12-mer
    at 25 nmol, which starts at 15.

    Raises for a length no standard scale takes, and that is the point rather
    than an inconvenience: the alternative is a ``Scale`` cell IDT cannot fill,
    which fails at the vendor days later instead of here.  Past
    :data:`IDT_STANDARD_MAX` an oligo needs an Ultramer, a different product
    this module does not order.
    """
    for code, low, high in IDT_SCALES:
        if low <= length <= high:
            return code
    if length > IDT_STANDARD_MAX:
        raise ValueError(
            f"{length} nt is past {IDT_STANDARD_MAX}, the longest oligo IDT's "
            "standard scales take; it needs an Ultramer, which this order "
            "sheet cannot buy"
        )
    raise ValueError(
        f"{length} nt is under {IDT_STANDARD_MIN}, the shortest oligo IDT's "
        "standard scales take"
    )


def _scales(result, scale: Optional[str]) -> List[str]:
    """The ``Scale`` column: one code per row, in order-sheet order.

    ``None`` lets each oligo's own length choose, which is the whole point of
    the ladder.  An explicit code is written on every row and CHECKED against
    the ladder where it is a code this module knows -- the failure being fixed
    here is a scale written on a row too long to hold it, and a caller can
    write one just as a hard-coded default could.

    A code the ladder does not list goes through UNCHECKED.  Ultramers, a
    modified oligo's scale, whatever IDT adds next: those are real orders, and
    a module that refused every code it had not been taught would be worse than
    one that trusts a caller who knows something it does not.
    """
    rows = _rows(result)
    if scale is None:
        out: List[str] = []
        for name, sequence in rows:
            try:
                out.append(scale_for_length(len(sequence)))
            except ValueError as exc:
                raise ValueError(f"oligo {name}: {exc}") from exc
        return out
    known = {code: (low, high) for code, low, high in IDT_SCALES}
    if scale in known:
        low, high = known[scale]
        over = sorted(
            {len(sequence) for _, sequence in rows if not low <= len(sequence) <= high}
        )
        if over:
            raise ValueError(
                f"scale {scale!r} takes {low}-{high} nt, but this order holds "
                f"oligos of {', '.join(f'{n} nt' for n in over)}; pass "
                "scale=None to let each oligo's length choose its scale"
            )
    return [scale] * len(rows)


def scale_mix(result, *, scale: Optional[str] = None) -> Tuple[Tuple[str, int], ...]:
    """``(scale, count)`` over the whole order, smallest scale first.

    What the summary says out loud, so that an order spanning two scales is
    visible as one before anybody pastes it into a form and reads a single
    number off the total.  Built from the same :func:`_scales` the sheet's own
    column is, so the count here and the cells there are one decision made
    once rather than two that can disagree.

    A scale the ladder does not list sorts last; there is no price order to put
    it in.
    """
    counts = Counter(_scales(result, scale))
    ladder = [code for code, _, _ in IDT_SCALES]
    return tuple(
        sorted(
            counts.items(),
            key=lambda kv: (
                ladder.index(kv[0]) if kv[0] in ladder else len(ladder),
                kv[0],
            ),
        )
    )


def tube_frame(
    result,
    *,
    scale: Optional[str] = None,
    purification: str = DEFAULT_PURIFICATION,
) -> "pd.DataFrame":
    """The order as tubes, for IDT's bulk entry form.

    One row per distinct molecule, which is one row per tube and one line on
    the invoice: an oligo several junctions share was already a single record
    by the time it got here, because the design identifies an oligo by its
    sequence and not by its name.

    ``scale`` is IDT's code, and LEFT ALONE IT IS PER ROW: each oligo at the
    cheapest scale its own length can be made at.  One code on every row was
    the earlier behaviour and it was wrong, because IDT's scales are
    length-limited and an order holds oligos of different lengths -- 25 nmol
    stops at 60 bases, a tailed junction primer runs past that routinely, and
    such a row is an order IDT will not fill.  Pass a code to put that code on
    every row instead, checked against the ladder: see :func:`_scales`.

    ``purification`` is written on every row, since unlike scale it is not
    length-bound.  The default is none, which is what a PCR primer needs, and
    it cannot usefully be raised on its own: IDT offers no purification at all
    at 25 nmol -- every column of its purification table is N/A on that row --
    so asking for one also means naming a ``scale`` of at least ``100nm``.
    """
    import pandas as pd

    rows = _rows(result)
    return pd.DataFrame(
        {
            IDT_TUBE_COLUMNS[0]: [name for name, _ in rows],
            IDT_TUBE_COLUMNS[1]: [sequence for _, sequence in rows],
            IDT_TUBE_COLUMNS[2]: _scales(result, scale),
            IDT_TUBE_COLUMNS[3]: [purification] * len(rows),
        },
        columns=list(IDT_TUBE_COLUMNS),
    )


def plate_frames(
    result, *, rows: int = PLATE_ROWS, cols: int = PLATE_COLS
) -> List["pd.DataFrame"]:
    """The same order laid out on plates, one frame each.

    Wells are filled in the order sheet's own order: ``001`` into ``A1``,
    ``002`` into ``B1`` below it, and an order past one plate starts the next
    at ``A1`` again.

    THE PAIRS ARE USUALLY NEIGHBOURS, AND THAT IS NOT A GUARANTEE.  Primers are
    numbered in the order the design's walk first met each molecule, so ``001``
    is the first fragment's forward primer and ``002`` its reverse -- the sheet
    reads forward, reverse, forward, reverse, and a fragment's pair lands in
    ``A1``/``B1``, one directly below the other.  But a molecule two junctions
    share is bought once, so where that happens one row is missing from the
    sheet and every later pair shifts up a well.

    Re-pairing them would mean a row per reaction instead of a row per
    molecule, which buys the same tube twice -- so the plate follows the sheet,
    and ``02_primers/uses.csv`` is where which junctions a well serves is
    written down.  What the 8-channel pipette needs is unaffected either way:
    wells fill down whole columns with no gaps, which is decided by
    :func:`well_positions` and the oligo count, not by which oligo is where.

    An order with no oligos gives NO plates, rather than one empty plate that
    IDT would refuse.
    """
    import pandas as pd

    wells = well_positions(rows, cols)
    capacity = len(wells)
    laid = _rows(result)

    plates: List["pd.DataFrame"] = []
    for start in range(0, len(laid), capacity):
        chunk = laid[start : start + capacity]
        plates.append(
            pd.DataFrame(
                {
                    IDT_PLATE_COLUMNS[0]: list(wells[: len(chunk)]),
                    IDT_PLATE_COLUMNS[1]: [name for name, _ in chunk],
                    IDT_PLATE_COLUMNS[2]: [sequence for _, sequence in chunk],
                },
                columns=list(IDT_PLATE_COLUMNS),
            )
        )
    return plates


def primers_fasta(result) -> str:
    """Every oligo as FASTA, for whatever else has to read them.

    A vendor that is not IDT, an alignment, a BLAST, a LIMS import: all of them
    read FASTA and none of them read a plate.  The header is the same name the
    order carries, so a record here and a tube there are the same molecule by
    name as well as by sequence.

    Unwrapped.  An oligo is under 100 bases and wrapping it at 60 would put a
    newline in the middle of the one thing a reader wants to copy.
    """
    return "".join(
        f">{name}\n{sequence}\n" for name, sequence in _rows(result)
    )


def plate_count(n_oligos: int, rows: int = PLATE_ROWS, cols: int = PLATE_COLS) -> int:
    """How many plates ``n_oligos`` fill.  Zero oligos fill no plates."""
    capacity = rows * cols
    return -(-n_oligos // capacity) if n_oligos > 0 else 0
