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
resolved target must stay inside ``out_dir``.  It is also made *writable* and
*unique as a filename*, which a legal id is not automatically: see
:func:`_stems`.
"""

from __future__ import annotations

import hashlib
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple, Union

from teemi.design.combinatorial.report import VariantPlan

__all__ = ["write_genbank_files", "write_assembly_figures"]


def write_genbank_files(
    designs: Iterable[VariantPlan], out_dir: Union[str, Path]
) -> Tuple[Path, ...]:
    """Write one GenBank file per design and return the paths written.

    A variant with no assembled construct -- verification found no product or
    more than one -- is not written, so the returned paths say what exists
    rather than a count having to be trusted.  The file is named for the
    variant where the id can carry it; see :func:`_stems`.
    """
    plans = [design for design in designs if design.construct is not None]
    stems = _stems([design.id for design in plans], ".gb")
    written: List[Path] = []
    for design in plans:
        path = _target(out_dir, stems[design.id], ".gb")
        path.write_text(design.construct.format("genbank"))
        written.append(path)
    return tuple(written)


def write_assembly_figures(
    designs: Iterable[VariantPlan], out_dir: Union[str, Path]
) -> Tuple[Path, ...]:
    """Write pydna's own assembly figure per design, one text file each."""
    plans = [design for design in designs if design.construct is not None]
    stems = _stems([design.id for design in plans], ".txt")
    written: List[Path] = []
    for design in plans:
        path = _target(out_dir, stems[design.id], ".txt")
        path.write_text(
            f"{design.id}  {design.length_bp} bp  "
            f"{'circular' if design.construct.circular else 'linear'}\n\n"
            f"{design.construct.figure()}\n"
        )
        written.append(path)
    return tuple(written)


#: Longest filename component the common filesystems accept, in bytes.
MAX_COMPONENT_BYTES = 255

#: Characters of digest appended to a stem that cannot be the id itself.
_DIGEST_CHARS = 8


def _stems(ids: Sequence[str], suffix: str) -> Dict[str, str]:
    """A writable, collision-free filename stem for each variant id.

    Two ids can both be legal and distinct and still name one file.  A variant
    id is up to ``n_positions * 64`` characters, so from four positions onward
    it runs past the filesystem's limit on a single name component and the write
    fails outright.  And ids differing only in case -- ``cdsA`` and ``CDSA`` are
    two valid, different parts -- are one file on a case-insensitive volume,
    which is the default on macOS: the second write silently replaced the first
    while the returned tuple still claimed both existed.

    So an id is used as-is only when it fits and is unique after case folding.
    Otherwise the stem is a truncation plus a digest **of the id**, which keeps
    it deterministic: the same design lands in the same file on every run, and
    the choice depends on the set of ids rather than their order.
    """
    folded = Counter(variant_id.casefold() for variant_id in ids)
    budget = MAX_COMPONENT_BYTES - len(suffix.encode()) - _DIGEST_CHARS - 1
    stems: Dict[str, str] = {}
    for variant_id in ids:
        fits = len(f"{variant_id}{suffix}".encode()) <= MAX_COMPONENT_BYTES
        unique = folded[variant_id.casefold()] == 1
        if fits and unique:
            stems[variant_id] = variant_id
            continue
        digest = hashlib.sha256(variant_id.encode()).hexdigest()[:_DIGEST_CHARS]
        stems[variant_id] = f"{variant_id.encode()[:budget].decode()}-{digest}"
    return stems


def _target(out_dir: Union[str, Path], stem: str, suffix: str) -> Path:
    directory = Path(out_dir).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    path = (directory / f"{stem}{suffix}").resolve()
    if path.parent != directory:
        raise ValueError(
            f"variant id {stem!r} would write outside {directory}; refusing"
        )
    return path
