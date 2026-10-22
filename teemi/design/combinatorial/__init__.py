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

"""Combinatorial DNA assembly library planning.

A JSON design spec in, a verified primer and PCR plan out::

    import teemi.design.combinatorial as cb
    from teemi.design.combinatorial.sources import SourceResolver

    result = cb.plan(spec, sources=SourceResolver(base_dir=repo), offline=True,
                     include_designs=True)
    result.raise_for_status()
    cb.primer_dataframe(result)
    cb.write_genbank_files(result.designs, out_dir)

The biology, in the order it happens:

1. **Footprints** come from one ``pydna.design.primer_design`` call per part,
   with NEB as the Tm function and the local model as pydna's cheap first
   estimate.  They are immutable afterwards.
2. **Overhangs** are split, not resized: at a junction with homology ``L`` the
   downstream forward primer carries the upstream part's last ``t_f`` bases and
   the upstream reverse primer the reverse complement of the downstream part's
   first ``t_r``, with ``t_f + t_r = L`` exactly.  Only the split is searched.
3. **Secondary structure** comes from primer3 and only from primer3, on the
   complete tailed oligo, never for a Tm.
4. **Temperatures** come from the local pydna models during the search and from
   one batched NEB request for the authoritative numbers.  The annealing
   temperature is always derived from the **footprint** pair.
5. **Verification** reconstructs the amplicons, assembles them with pydna and
   compares strand- and rotation-invariantly against the intended construct.

One number from one engine is never checked against a limit calibrated for
another: a temperature is a
:class:`~teemi.design.combinatorial.tm.Tm` carrying its engine and having no
arithmetic, and a recorded breach has to name the single source behind both
sides.

The legacy :class:`teemi.design.combinatorial_design.DesignAssembly` is
untouched and still available.

Module map::

    flow.py       FLOW's design document as the input boundary
    spec.py       the backend's own spec: the policy FLOW never sends
    sources.py    part sources -> Dseqrecord; the untrusted-input boundary
    tm.py         Tm and Ta, each tagged with the engine behind it
    design.py     footprints, the split sweep, variants, PCR checks
    simulate.py   assemble it and compare with the intended construct
    report.py     plan(), the result, the three tables
    export.py     the file-writing edge: GenBank and assembly figures
"""

from teemi.design.combinatorial.design import Diagnostic, Violation
from teemi.design.combinatorial.flow import (
    FlowDocument,
    UndecidedPosition,
    UnknownPart,
    from_flow,
    plan_flow,
)
from teemi.design.combinatorial.export import (
    write_assembly_figures,
    write_genbank_files,
)
from teemi.design.combinatorial.report import (
    DesignResult,
    InfeasibleDesignError,
    pcr_dataframe,
    plan,
    primer_dataframe,
    variant_overview_dataframe,
)
from teemi.design.combinatorial.simulate import Verification
from teemi.design.combinatorial.sources import SourceError, SourceResolver
from teemi.design.combinatorial.spec import DesignSpec, json_schema, load_spec

__all__ = [
    "plan",
    "DesignResult",
    "InfeasibleDesignError",
    "variant_overview_dataframe",
    "primer_dataframe",
    "pcr_dataframe",
    "write_genbank_files",
    "write_assembly_figures",
    "DesignSpec",
    "load_spec",
    "json_schema",
    "SourceResolver",
    "SourceError",
    "Diagnostic",
    "Violation",
    "Verification",
    "FlowDocument",
    "from_flow",
    "plan_flow",
    "UndecidedPosition",
    "UnknownPart",
]
