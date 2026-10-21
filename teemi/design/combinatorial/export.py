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

"""Writing the assembled constructs out.

The only file-writing edge of the package.  :func:`plan` is pure; the caller
hands ``result.designs`` to one of these.

Nothing here builds annotations.  Each part record is annotated once in
:mod:`teemi.design.combinatorial.sources`, ``pydna.amplify.pcr`` adds a
``primer_bind`` feature per primer with that primer's own id, and pydna's
assembly propagates both into the product -- so the GenBank file is annotated
because the fragments were, and constructing features on the finished construct
would only risk disagreeing with it.

A variant id becomes a filename, so it is checked against the output directory
the way :mod:`.sources` checks an input path against its ``base_dir``: the
resolved target must stay inside ``out_dir``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Tuple, Union

from teemi.design.combinatorial.report import VariantPlan

__all__ = ["write_genbank_files", "write_assembly_figures"]


def write_genbank_files(
    designs: Iterable[VariantPlan], out_dir: Union[str, Path]
) -> Tuple[Path, ...]:
    """Write ``<variant_id>.gb`` per design and return the paths written.

    A variant with no assembled construct -- verification found no product or
    more than one -- is not written, so the returned paths say what exists
    rather than a count having to be trusted.
    """
    written: List[Path] = []
    for design in designs:
        if design.construct is None:
            continue
        path = _target(out_dir, design.id, ".gb")
        path.write_text(design.construct.format("genbank"))
        written.append(path)
    return tuple(written)


def write_assembly_figures(
    designs: Iterable[VariantPlan], out_dir: Union[str, Path]
) -> Tuple[Path, ...]:
    """Write ``<variant_id>.txt`` with pydna's own assembly figure per design."""
    written: List[Path] = []
    for design in designs:
        if design.construct is None:
            continue
        path = _target(out_dir, design.id, ".txt")
        path.write_text(
            f"{design.id}  {design.length_bp} bp  "
            f"{'circular' if design.construct.circular else 'linear'}\n\n"
            f"{design.construct.figure()}\n"
        )
        written.append(path)
    return tuple(written)


def _target(out_dir: Union[str, Path], stem: str, suffix: str) -> Path:
    directory = Path(out_dir).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    path = (directory / f"{stem}{suffix}").resolve()
    if path.parent != directory:
        raise ValueError(
            f"variant id {stem!r} would write outside {directory}; refusing"
        )
    return path
