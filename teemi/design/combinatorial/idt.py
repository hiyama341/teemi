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
  entry form.  The scale and purification codes are IDT's own (``25nm``,
  ``STD``).
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
from typing import TYPE_CHECKING, List, Sequence, Tuple

if TYPE_CHECKING:  # pragma: no cover - import cycle at runtime, types only
    import pandas as pd

__all__ = [
    "IDT_PLATE_COLUMNS",
    "IDT_TUBE_COLUMNS",
    "IDT_PLATE_MINIMUM",
    "DEFAULT_SCALE",
    "DEFAULT_PURIFICATION",
    "well_positions",
    "oligo_name",
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

#: IDT's codes for the smallest synthesis scale and for no purification.
DEFAULT_SCALE = "25nm"
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


def tube_frame(
    result,
    *,
    scale: str = DEFAULT_SCALE,
    purification: str = DEFAULT_PURIFICATION,
) -> "pd.DataFrame":
    """The order as tubes, for IDT's bulk entry form.

    One row per distinct molecule, which is one row per tube and one line on
    the invoice: an oligo several junctions share was already a single record
    by the time it got here, because the design identifies an oligo by its
    sequence and not by its name.

    ``scale`` and ``purification`` are IDT's codes and are written on every
    row, since the form takes them per oligo.  The defaults are the smallest
    scale and no purification, which is what a PCR primer needs; purification
    at IDT requires a scale of at least 100 nmol, so the two defaults cannot be
    changed independently without checking that pairing.
    """
    import pandas as pd

    rows = _rows(result)
    return pd.DataFrame(
        {
            IDT_TUBE_COLUMNS[0]: [name for name, _ in rows],
            IDT_TUBE_COLUMNS[1]: [sequence for _, sequence in rows],
            IDT_TUBE_COLUMNS[2]: [scale] * len(rows),
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
