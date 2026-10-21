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

"""Resolving a spec's part sources to sequence records.

This is the untrusted-input boundary of the package.  Everything downstream may
assume: uppercase ACGT only, at most :data:`MAX_PART_BP` bases, a known
topology, and a :class:`~pydna.dseqrecord.Dseqrecord` whose id is the part id.

Two properties are structural rather than checked:

* **No network.**  A spec has no URL field, and this module makes no request, so
  there is nothing to point at an internal address.
* **No path escape.**  A ``file`` source is relative, is resolved under an
  explicit ``base_dir``, and is refused if the real path lands outside it --
  which also refuses a symlink out, since :meth:`pathlib.Path.resolve` follows
  links before the check.

Files are read with Biopython in the format the spec declares.  Note that
``pydna.parsers.parse`` reads a *file* when handed a string that happens to be a
path, which in this module would be an arbitrary-file-read vector, so it is not
used here.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pydna.dseqrecord import Dseqrecord

from teemi.design.combinatorial.spec import (
    DesignSpec,
    FileSource,
    InlineSource,
    Part,
)

__all__ = [
    "DNA_ALPHABET",
    "MAX_PART_BP",
    "MAX_FILE_BYTES",
    "SourceError",
    "SourceResolver",
]

DNA_ALPHABET = frozenset("ACGT")
"""Unambiguous DNA only.  A degenerate code would make Tm, the primer footprint
and the assembled construct all ambiguous, so it is rejected rather than
guessed at."""

MAX_PART_BP = 50_000
MAX_FILE_BYTES = 10_000_000


class SourceError(ValueError):
    """A part source that cannot be resolved into a valid record.

    :attr:`part_id` lets the caller attribute the failure to a part in a
    diagnostic without parsing the message.
    """

    def __init__(self, message: str, *, part_id: Optional[str] = None) -> None:
        super().__init__(message if part_id is None else f"part {part_id!r}: {message}")
        self.part_id = part_id


class SourceResolver:
    """Resolves part sources to records, under a fixed filesystem boundary.

    ::

        resolver = SourceResolver(base_dir=repo_root)
        records = resolver.resolve_spec(spec)

    ``base_dir`` is the only directory any ``file`` source can read from; left
    at ``None``, file sources are refused outright and only ``inline`` works.
    """

    def __init__(
        self,
        base_dir: Optional[Union[str, Path]] = None,
        *,
        max_part_bp: int = MAX_PART_BP,
        max_file_bytes: int = MAX_FILE_BYTES,
    ) -> None:
        self.base_dir = None if base_dir is None else Path(base_dir).resolve()
        self.max_part_bp = max_part_bp
        self.max_file_bytes = max_file_bytes
        self._files: Dict[Path, List[SeqRecord]] = {}

    def resolve_spec(self, spec: DesignSpec) -> Dict[str, Dseqrecord]:
        """Every part and pad of the spec, keyed by part id (unique per spec)."""
        return {
            part.id: self.resolve(part)
            for _, parts in spec.ordered_slots()
            for part in parts
        }

    def resolve(self, part: Part) -> Dseqrecord:
        """One part (or pad) as a validated record.

        The record carries a whole-length ``misc_feature`` labelled with the part
        id, plus any features the source file had.  pydna propagates features
        from the assembly fragments into the contig, so annotating here is what
        makes the exported construct annotated.
        """
        source = part.source
        if isinstance(source, InlineSource):
            sequence = self._validated_dna(source.sequence, part.id, "inline sequence")
            features: List = []
        else:
            record = self._file_record(source, part.id)
            sequence = self._validated_dna(
                str(record.seq), part.id, f"{source.path}:{record.id}"
            )
            features = list(record.features)

        bio = SeqRecord(
            Seq(sequence),
            id=part.id,
            name=part.id[:16],  # GenBank LOCUS field
            description=part.name or "",
            features=features,
            annotations={"molecule_type": "DNA"},
        )
        resolved = Dseqrecord(bio, circular=source.topology == "circular")
        resolved.add_feature(0, len(resolved), type_="misc_feature", label=part.id)
        return resolved

    def _validated_dna(self, text: str, part_id: str, origin: str) -> str:
        sequence = text.strip().upper()
        if not sequence:
            raise SourceError(f"{origin} is empty", part_id=part_id)
        unexpected = sorted(set(sequence) - DNA_ALPHABET)
        if unexpected:
            raise SourceError(
                f"{origin} has characters outside ACGT: " f"{''.join(unexpected[:8])}",
                part_id=part_id,
            )
        if len(sequence) > self.max_part_bp:
            raise SourceError(
                f"{origin} is {len(sequence)} bp, over the "
                f"{self.max_part_bp} bp cap",
                part_id=part_id,
            )
        return sequence

    def _file_record(self, source: FileSource, part_id: str) -> SeqRecord:
        records = self._parse(self._resolve_path(source, part_id), source.format)
        if source.record_id is None:
            if len(records) != 1:
                raise SourceError(
                    f"{source.path} holds {len(records)} records; "
                    "set record_id to say which",
                    part_id=part_id,
                )
            return records[0]
        matches = [r for r in records if source.record_id in (r.id, r.name)]
        if not matches:
            available = ", ".join(r.id for r in records[:10])
            raise SourceError(
                f"no record {source.record_id!r} in {source.path} "
                f"(found: {available})",
                part_id=part_id,
            )
        if len(matches) > 1:
            raise SourceError(
                f"record {source.record_id!r} is not unique in {source.path}",
                part_id=part_id,
            )
        return matches[0]

    def _resolve_path(self, source: FileSource, part_id: str) -> Path:
        if self.base_dir is None:
            raise SourceError(
                "file sources need SourceResolver(base_dir=...)", part_id=part_id
            )
        if Path(source.path).is_absolute():
            raise SourceError(
                f"{source.path} is absolute; paths are relative to base_dir",
                part_id=part_id,
            )
        try:
            path = (self.base_dir / source.path).resolve()
            inside = path.is_relative_to(self.base_dir)
        except (OSError, ValueError) as exc:
            raise SourceError(
                f"{source.path} is not a usable path: {exc}", part_id=part_id
            )
        if not inside:
            raise SourceError(
                f"{source.path} resolves outside base_dir {self.base_dir}",
                part_id=part_id,
            )
        if not path.is_file():
            raise SourceError(f"no such file: {source.path}", part_id=part_id)
        if path.stat().st_size > self.max_file_bytes:
            raise SourceError(
                f"{source.path} is larger than the {self.max_file_bytes} byte cap",
                part_id=part_id,
            )
        return path

    def _parse(self, path: Path, file_format: str) -> List[SeqRecord]:
        if path not in self._files:
            try:
                records = list(SeqIO.parse(path, file_format))
            except (OSError, ValueError) as exc:
                raise SourceError(f"cannot read {path.name} as {file_format}: {exc}")
            if not records:
                raise SourceError(f"{path.name} holds no {file_format} record")
            self._files[path] = records
        return self._files[path]
