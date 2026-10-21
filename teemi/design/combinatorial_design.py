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

"""Design combinatorial libraries of DNA parts.

The library is the cartesian product of the parts offered for each position of a
construct. Every combination becomes a list of PCR products that overlap their
neighbours, ready to be fused by Gibson assembly or in-vivo homologous
recombination.
"""

import itertools
import math
from typing import Callable, Sequence

import pandas as pd
from pydna.amplicon import Amplicon
from pydna.design import assembly_fragments, primer_design
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer
from pydna.tm import tm_default

# Fragments up to this length are carried by the tails of the flanking primers
# instead of being amplified themselves, which is what makes a pad a pad.
MAX_PAD_LENGTH = 40

PRIMER_PRICE_PER_BASE = 1.8  # DKK


class DesignAssembly:
    """A combinatorial library of overlapping PCR products.

    Parameters
    ----------
    list_of_seqs : Sequence[Sequence[Dseqrecord]]
        One list of interchangeable parts per position in the construct, in the
        order they should appear. ``[[p1, p2], [c1]]`` therefore describes two
        variants. Parts are labelled by their ``name`` in the reported tables.
    list_of_pads : Sequence[Dseqrecord]
        Short sequences, at most ``MAX_PAD_LENGTH`` bp, to splice into every
        variant. A pad rides along in the primer tails, so it adds sequence to
        the construct without adding a PCR to the plan.
    positions_of_pads : Sequence[int]
        Index in ``list_of_seqs`` each pad is inserted before, so ``[0]`` puts a
        pad in front of the first part and ``[len(list_of_seqs)]`` after the last.
    target_tm, limit : float, int
        Target melting temperature and minimum length of the primer footprints.
    overlap : int
        Length of the homology neighbouring amplicons share.
    tm_func : Callable[[str], float]
        Melting temperature function, also used for the reported temperatures.

    Attributes
    ----------
    list_of_assemblies : list[list[Amplicon]]
        The amplicons of each variant, in construct order.
    primers : list[Primer]
        The distinct primers of the whole library, numbered ``P001`` upwards.
    amplicons : list[Amplicon]
        The distinct PCR products of the whole library.
    combinatorial_list_of_names : list[tuple[str, ...]]
        The part names of each variant.
    systematic_names : list[tuple[int, ...]]
        The 1-based part index per position of each variant, e.g. ``(1, 2, 1)``.
    """

    def __init__(
        self,
        list_of_seqs: Sequence[Sequence[Dseqrecord]],
        list_of_pads: Sequence[Dseqrecord] = (),
        positions_of_pads: Sequence[int] = (),
        target_tm: float = 55.0,
        limit: int = 13,
        overlap: int = 35,
        tm_func: Callable[[str], float] = tm_default,
    ):
        if not list_of_seqs or not all(list_of_seqs):
            raise ValueError(
                "list_of_seqs must hold one non-empty list of parts per position"
            )
        for parts in list_of_seqs:
            for part in parts:
                if not isinstance(part, Dseqrecord):
                    raise TypeError(
                        f"parts must be Dseqrecord objects, got {type(part).__name__}; "
                        "wrap Biopython records with Dseqrecord(record) first"
                    )
        if len(list_of_pads) != len(positions_of_pads):
            raise ValueError(
                "list_of_pads and positions_of_pads must be the same length"
            )
        for pad, position in zip(list_of_pads, positions_of_pads):
            if len(pad) > MAX_PAD_LENGTH:
                raise ValueError(
                    f"pad {pad.name!r} is {len(pad)} bp; pads over {MAX_PAD_LENGTH} bp "
                    "do not fit in the primer tails and have to be parts of their own"
                )
            if not 0 <= position <= len(list_of_seqs):
                raise ValueError(
                    f"pad position {position} is outside the construct, which has "
                    f"{len(list_of_seqs)} positions"
                )

        self._tm_func = tm_func
        self.combinatorial_list_of_names = list(
            itertools.product(
                *([part.name for part in parts] for parts in list_of_seqs)
            )
        )
        self.systematic_names = list(
            itertools.product(*(range(1, len(parts) + 1) for parts in list_of_seqs))
        )

        amplicons = [
            [
                primer_design(part, target_tm=target_tm, limit=limit, tm_func=tm_func)
                for part in parts
            ]
            for parts in list_of_seqs
        ]
        # Insert the pads back to front so the earlier positions stay valid.
        pads = sorted(
            zip(positions_of_pads, list_of_pads), key=lambda pad: pad[0], reverse=True
        )

        self.list_of_assemblies = []
        for combination in itertools.product(*amplicons):
            fragments = list(combination)
            for position, pad in pads:
                fragments.insert(position, pad)
            assembly = assembly_fragments(
                fragments, overlap=overlap, maxlink=MAX_PAD_LENGTH
            )
            _describe(assembly)
            self.list_of_assemblies.append(assembly)

        self.primers = _share_primers(self.list_of_assemblies)
        self.amplicons = list(
            {
                str(amplicon.seq): amplicon
                for assembly in self.list_of_assemblies
                for amplicon in assembly
            }.values()
        )

    def show_variants_lib_df(self) -> pd.DataFrame:
        """One row per variant: its parts, its systematic name and its number."""
        variants = pd.DataFrame(self.combinatorial_list_of_names)
        variants["Systematic_name"] = self.systematic_names
        variants["Variant"] = range(len(variants))
        return variants

    def primer_list_to_dataframe(self) -> pd.DataFrame:
        """One row per primer to order."""
        return pd.DataFrame(
            [
                {
                    "id": primer.id,
                    "anneals to": primer.name,
                    "sequence": str(primer.seq),
                    "annealing temperature": self._footprint_tm(primer),
                    "length": len(primer),
                    "price(DKK)": len(primer) * PRIMER_PRICE_PER_BASE,
                    "description": primer.description,
                    "footprint": str(primer.footprint),
                    "len_footprint": len(primer.footprint),
                }
                for primer in self.primers
            ]
        )

    def pcr_list_to_dataframe(self) -> pd.DataFrame:
        """One row per PCR to run."""
        return pd.DataFrame(
            [
                {
                    "pcr_number": f"PCR{number}",
                    "template": amplicon.name,
                    "forward_primer": amplicon.forward_primer.id,
                    "reverse_primer": amplicon.reverse_primer.id,
                    "f_tm": self._footprint_tm(amplicon.forward_primer),
                    "r_tm": self._footprint_tm(amplicon.reverse_primer),
                }
                for number, amplicon in enumerate(self.amplicons, start=1)
            ]
        )

    def _footprint_tm(self, primer: Primer) -> float:
        """Melting temperature of the footprint; the assembly tail does not anneal."""
        return round(self._tm_func(str(primer.footprint)), 2)


def _describe(assembly: list[Amplicon]) -> None:
    """Name the amplicons of one variant, and its primers, after the parts they amplify."""
    for amplicon in assembly:
        amplicon.name = amplicon.template.name
    for position, amplicon in enumerate(assembly):
        upstream = assembly[position - 1].name if position else None
        downstream = (
            assembly[position + 1].name if position + 1 < len(assembly) else None
        )
        for primer, neighbour in (
            (amplicon.forward_primer, upstream),
            (amplicon.reverse_primer, downstream),
        ):
            primer.name = amplicon.name
            primer.description = f"Anneals to {amplicon.name}"
            if neighbour is not None:
                primer.description += f", overlaps to {neighbour}"


def _share_primers(assemblies: list[list[Amplicon]]) -> list[Primer]:
    """Number the distinct primers P001 upwards and point every amplicon at them.

    Variants that share a junction share a primer, which should only be ordered once.
    """
    by_sequence: dict[str, Primer] = {}
    for assembly in assemblies:
        for amplicon in assembly:
            amplicon.forward_primer = by_sequence.setdefault(
                str(amplicon.forward_primer.seq), amplicon.forward_primer
            )
            amplicon.reverse_primer = by_sequence.setdefault(
                str(amplicon.reverse_primer.seq), amplicon.reverse_primer
            )
    primers = list(by_sequence.values())
    for number, primer in enumerate(primers, start=1):
        primer.id = f"P{number:03}"
    return primers


def count_unique_parts(
    predictions_df: pd.DataFrame, max_combinations: int
) -> dict[str, list | str]:
    """Take the top-ranked predictions needed to span ``max_combinations`` constructs.

    Each row of ``predictions_df`` is one predicted construct and each column one
    position in it. Reading from the top, parts are collected until they combine
    into at least ``max_combinations`` constructs, which keeps the library that has
    to be built as small as the target allows.

    Returns the parts collected per column, plus ``sum_of_parts`` (the number of
    combinations they span) and ``prediction_number`` (the number of rows read),
    both as strings.
    """
    parts: dict[str, list] = {column: [] for column in predictions_df.columns}
    combinations = 0
    for rows_read, (_, prediction) in enumerate(predictions_df.iterrows(), start=1):
        for column, collected in parts.items():
            if prediction[column] not in collected:
                collected.append(prediction[column])
        combinations = math.prod(len(collected) for collected in parts.values())
        if combinations >= max_combinations:
            return {
                **parts,
                "sum_of_parts": str(combinations),
                "prediction_number": str(rows_read),
            }
    raise ValueError(
        f"the {len(predictions_df)} predictions span only {combinations} combinations, "
        f"short of the {max_combinations} requested"
    )
