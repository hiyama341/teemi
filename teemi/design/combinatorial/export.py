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

import copy
import hashlib
import importlib
import json
import platform
import re
import shutil
import sys
import tempfile
import zipfile
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from teemi.design.combinatorial import idt
from teemi.design.combinatorial.report import VariantPlan

__all__ = [
    "write_genbank_files",
    "write_assembly_figures",
    "write_bundle",
    "write_bundle_zip",
    "build_summary",
    "bundle_layout",
    "constructs_genbank",
    "construct_maps",
    "environment_record",
    "ALL_CONSTRUCTS",
]


def write_genbank_files(
    designs: Iterable[VariantPlan], out_dir: Union[str, Path]
) -> Tuple[Path, ...]:
    """Write one GenBank file per design and return the paths written.

    A variant with no assembled construct -- verification found no product or
    more than one -- is not written, so the returned paths say what exists
    rather than a count having to be trusted.

    The name is the construct's number and then the variant id, so the files
    sort into library order and a map can be asked for by number; the id is
    kept because it says what the construct is made of.  See :func:`_stems`
    for what happens when the id cannot fit in a filename.
    """
    plans = [design for design in designs if design.construct is not None]
    names = {design.id: f"{design.label}_{design.id}" for design in plans}
    stems = _stems(list(names.values()), ".gb")
    written: List[Path] = []
    for design in plans:
        path = _target(out_dir, stems[names[design.id]], ".gb")
        path.write_text(design.construct.format("genbank"))
        written.append(path)
    return tuple(written)


def write_assembly_figures(
    designs: Iterable[VariantPlan], out_dir: Union[str, Path]
) -> Tuple[Path, ...]:
    """Write pydna's own assembly figure per design, one text file each."""
    plans = [design for design in designs if design.construct is not None]
    names = {design.id: f"{design.label}_{design.id}" for design in plans}
    stems = _stems(list(names.values()), ".txt")
    written: List[Path] = []
    for design in plans:
        path = _target(out_dir, stems[names[design.id]], ".txt")
        path.write_text(
            f"construct {design.label}  {design.id}  {design.length_bp} bp  "
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


# --------------------------------------------------------------------------- #
# The whole build, in one folder
# --------------------------------------------------------------------------- #

#: The one file in ``05_constructs`` that is NOT a map of a single construct.
#: Above :data:`BUNDLE_ORDER` because that table names it.  ``00_`` so it sorts
#: above the numbered maps it is the concatenation of.
ALL_CONSTRUCTS = "00_all_constructs.gb"

#: What a bundle contains, and in what order a person should read it.
#: The order to read the bundle in, which is why the entries are numbered: a
#: folder listing is alphabetical and would otherwise open on ``constructs/``.
#: Each entry is its key, its path, and the lines describing it -- a folder
#: holding two files worth naming gets two.  ``00_LAYOUT.txt`` describes
#: itself by being the map, so it has none.
BUNDLE_ORDER = (
    ("layout", "00_LAYOUT.txt", ()),
    (
        "summary",
        "01_SUMMARY.txt",
        (("", "read this first: what to order, run, look at twice"),),
    ),
    (
        "primers",
        "02_primers",
        (
            ("/primers.csv", "the order sheet, one row per distinct molecule"),
            (
                "/oligo_order_idt.csv",
                "the same oligos as an ORDER: paste into IDT's bulk tube form "
                "(Sequence is 5'-3')",
            ),
            (
                "/oligo_order_idt_plate.xlsx",
                "or upload this; one worksheet per 96-well plate, wells filled "
                "down the columns (Sequence is 5'-3')",
            ),
            ("/primers.fasta", "every oligo as FASTA, for whatever takes no CSV"),
            (
                "/analysis.csv",
                "every column the primer table has, not just the readable ones",
            ),
            ("/uses.csv", "which junctions each oligo serves, one row each"),
        ),
    ),
    (
        "pcrs",
        "03_pcrs",
        (
            ("/pcrs.csv", "one row per reaction"),
            (
                "/worklist_by_ta.csv",
                "the same reactions sorted by the temperature they anneal at",
            ),
        ),
    ),
    (
        "parts",
        "04_parts",
        (("/parts.csv", "what to source, how often each is used, and what is known"),),
    ),
    (
        "constructs",
        "05_constructs",
        (
            ("/constructs.csv", "what each construct is made of, one row each"),
            ("/NNN_*.gb", "one annotated map per verified construct, numbered"),
            (
                f"/{ALL_CONSTRUCTS}",
                "every one of those maps in a single file, to import in one step",
            ),
            ("/assemblies.txt", "every assembly figure, in one file"),
        ),
    ),
    (
        "request",
        "06_request.json",
        (("", "the document as it was handed in, before validation"),),
    ),
    (
        "spec",
        "07_spec.json",
        (("", "every setting that was used, defaults filled in; re-runs this"),),
    ),
    (
        "design",
        "08_design.json",
        (("", "the whole result: junctions, primers, PCRs, diagnostics"),),
    ),
    (
        "run",
        "09_run.json",
        (("", "when this was built, by which versions, against which library"),),
    ),
)

#: Where each thing is, by the key :data:`BUNDLE_ORDER` gives it.
BUNDLE_PATHS = {key: path for key, path, _ in BUNDLE_ORDER}

#: An entry named the way this module names them.  Used to recognise an entry
#: from an older numbering, which no current name would match.
_NUMBERED = re.compile(r"^\d\d_")


def construct_maps(folder: Union[str, Path]) -> Tuple[Path, ...]:
    """The per-construct maps in a written bundle, in filename order.

    Filename order is numbered order up to 999 constructs and not past it: the
    number is zero-padded to three digits (``report``'s ``:03d``), so a library
    of 1000 sorts ``1000_`` ahead of ``100_``.  Said here rather than fixed
    here, because the fix worth making is the padding -- which renames every
    map in every bundle -- and not a sort this function alone would get right.

    A glob for ``*.gb`` is the obvious way to ask this and it is WRONG, because
    :data:`ALL_CONSTRUCTS` sits in the same folder and is every construct at
    once.  Three separate tests asked that way, at four call sites, and all
    three broke the day the combined file arrived -- two counted nine maps for
    eight constructs, and the third handed a multi-record file to
    ``SeqIO.read``, which takes exactly one, at both of its call sites.

    So the question has an owner.  ``folder`` is a bundle root.
    """
    constructs = Path(folder) / BUNDLE_PATHS["constructs"]
    return tuple(
        path
        for path in sorted(constructs.glob("*.gb"))
        if path.name != ALL_CONSTRUCTS
    )


def bundle_layout(
    keys: Optional[Iterable[str]] = None, *, omit: Iterable[str] = ()
) -> str:
    """The map of the bundle, describing the entries named by ``keys``.

    Built from what was actually written rather than from a fixed list,
    because a map naming a file the bundle did not produce is worse for
    someone finding their way around than no map at all -- and three of the
    entries are conditional: there is no spec for a refused design, no
    submitted document unless the caller passed one, and no run record unless
    the caller said when and against what.

    ``omit`` drops INDIVIDUAL LINES from an entry that was itself written,
    which ``keys`` cannot express.  One file needs it: the IDT plate workbook
    lives in ``02_primers`` beside the primer tables, so that entry is written
    whenever there are primers at all, but a result with no primers has
    nothing to order and openpyxl cannot write a workbook with no worksheets
    in it.  Paths are given as they appear here, e.g.
    ``02_primers/oligo_order_idt_plate.xlsx``.

    The gutter is 40 rather than the 30 it was: the plate workbook's path is
    37 characters and a narrower column pushed its description out of
    alignment with every other line.
    """
    wanted = set(BUNDLE_PATHS) if keys is None else set(keys)
    dropped = set(omit)
    lines = [
        f"{path + suffix:<40} {text}"
        for key, path, described in BUNDLE_ORDER
        if key in wanted
        for suffix, text in described
        if path + suffix not in dropped
    ]
    return "\n".join(lines) + "\n" if lines else ""


def spec_json(spec) -> str:
    """The spec as it was actually used, as JSON text.

    Not the document the caller handed in: the *validated* spec, so every
    default the caller never typed is written down -- the Tm the footprints
    were walked to, the dG floors a flag was raised against, which polymerase
    named the annealing temperature.  That is what "which settings were used"
    means, and a document that leaves them out cannot answer it.

    Sources are kept whole, inline sequences included, so the file re-runs the
    design on its own: :func:`load_spec` accepts it back unchanged.  For a
    67-gene library that is about 100 kB, which is a fifth of the result beside
    it and the price of an archive that does not depend on the machine it was
    built on.
    """
    return json.dumps(spec.model_dump(mode="json", exclude_none=True), indent=2)


def constructs_genbank(designs: Iterable[VariantPlan]) -> str:
    """Every verified construct in one GenBank file, in the numbered order.

    One file to import rather than one import per construct, which is the
    difference between a click and a morning when a library has sixty of them.

    IT RE-STAMPS THE LOCUS NAME, and that is the whole reason this is not a
    `cat` of the files :func:`write_genbank_files` writes.  A construct's LOCUS
    is set to ``variant.id[:16]`` when it is assembled, and a variant id is the
    part ids joined in track order -- so every variant that differs only at a
    LATER position has the same first sixteen characters.  The shipped
    combinatorial fixtures collapse eight constructs onto four LOCUS names and
    twelve onto eight.  Separate files hide it because the FILENAMES carry the
    construct number and the full id; one file would show duplicate LOCUS
    lines, and a reader with no way to tell two records apart.

    So the LOCUS becomes the construct number first and as much of the id as
    fits.  The number is unique within a plan by construction, so the records
    are distinguishable however the ids collide.

    The records are COPIED before they are touched.  ``result.designs`` is held
    for the caller's own use -- FLOW draws its on-screen tables from the same
    object -- and a writer that renamed them in place would change what the
    next reader sees.
    """
    plans = [design for design in designs if design.construct is not None]
    out = []
    for design in plans:
        record = copy.copy(design.construct)
        record.annotations = dict(record.annotations)
        record.annotations.setdefault("molecule_type", "DNA")
        # 16 characters is the LOCUS field; the number leads so it cannot clash.
        record.name = f"{design.label}_{design.id}"[:16]
        text = record.format("genbank")
        # pydna ends a record at `//` with NO trailing newline, so joining the
        # formatted records directly gives `//LOCUS` on one line and a file
        # that holds one record as far as any parser is concerned.  The test
        # for this reads the file back with SeqIO rather than counting LOCUS
        # lines, because counting is what missed it.
        out.append(text if text.endswith("\n") else text + "\n")
    return "".join(out)


#: Packages whose version changes what a design comes out as.
_MEASURED = ("teemi", "pydna", "primer3", "pandas", "Bio")


def environment_record() -> Dict[str, object]:
    """The versions that produced this build, read rather than declared.

    A bundle is meant to make a run repeatable.  The spec beside it says what
    was asked for; this says what answered -- and a primer set is a function of
    the Tm engine and the structure predictor as much as of the spec, so
    "which teemi" and "which primer3" are part of the result.
    """
    versions: Dict[str, object] = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
    }
    for name in _MEASURED:
        try:
            module = importlib.import_module(name)
        except Exception:  # pragma: no cover - an optional package absent
            continue
        versions[name] = str(getattr(module, "__version__", "unknown"))
    return versions


def run_json(run: Optional[Mapping[str, object]]) -> str:
    """The run record: what the caller knows, then what this package measures.

    ``run`` is provenance teemi cannot know -- when the design was computed,
    which commit of a part library resolved the ids, which service asked --
    exactly the argument ``part_info`` already makes, so it arrives the same
    way rather than being invented here.

    TEEMI'S OWN KEYS ARE WRITTEN LAST so a caller cannot shadow them, and a
    caller that tries RAISES.  A bundle reporting a teemi version its own code
    did not produce is worse than one reporting none, and the clash is far more
    likely to be a mistake than an intention.

    Nothing here reads the clock.  A timestamp is the caller's: it is the one
    fact that has to agree with the caller's own log lines and filenames, and a
    second reading of the clock here would differ from the first.
    """
    measured = environment_record()
    supplied = dict(run or {})
    clash = sorted(set(supplied) & set(measured))
    if clash:
        raise ValueError(
            f"run may not set {', '.join(clash)}: this package measures those "
            "itself, and a bundle naming a version it was not built with "
            "cannot be checked"
        )
    return json.dumps({**supplied, **measured}, indent=2, default=str)


def write_bundle(
    result,
    out_dir: Union[str, Path],
    *,
    part_info: Optional[Mapping[str, Mapping[str, object]]] = None,
    request: Optional[Mapping[str, object]] = None,
    run: Optional[Mapping[str, object]] = None,
    overwrite: bool = False,
) -> Path:
    """Write the whole build to one folder and return it.

    Everything a person needs to go and do the work, in the order they need it:
    what to order, what to run, what to source, what the constructs should come
    out as, and the summary that says which of it is worth looking at twice.
    The entries are numbered in that order, because a folder listing is
    alphabetical and would otherwise open on ``constructs/``.
    :func:`bundle_layout` is the map, written into the folder from the entries
    that were actually produced.

    ``part_info`` is provenance this package cannot know -- an Addgene number, a
    DOI, which plasmid a part came from.  The design carries sequences and ids;
    where a part can be *obtained* lives in whatever library the caller
    resolved it from, so the caller passes it in, keyed by part id.  Any keys
    are accepted and become columns, so a different library's fields need no
    change here.

    Only verified constructs are written, because an unverified one has no
    single product to draw.  The summary says how many were left out.

    The spec goes in beside the result, because the result does not contain
    it: a reader asking which floors, which polymerase or which flags produced
    these primers had nowhere to look, and nothing in the folder could be
    re-run.

    ``run`` is when this was built and against what: a timestamp, a library
    version, whatever names the caller's own run.  Same argument as
    ``part_info`` -- provenance this package cannot know -- and
    :func:`run_json` adds the versions it can measure.  Omitted, there is no
    ``09_run.json``, because a run record that could only say "some teemi,
    some day" is not one.

    ``request`` is the document as the caller received it -- a FLOW document,
    say -- archived unchanged.  The spec beside it is what teemi made of that
    document, which is not the same thing: validation fills in defaults,
    resolves sources and renames an inline variant.  Keeping both is what lets
    someone see that a setting they never wrote was nonetheless used.

    A non-empty folder is refused unless ``overwrite`` says otherwise, because
    a second run writing into the first leaves both: drop a part and its maps
    stay, and someone reads half of one design and half of another with
    nothing to tell them.

    ``overwrite`` removes the entries in :data:`BUNDLE_PATHS` **and anything
    else numbered the way this function numbers things**, which is how a
    renamed entry gets cleaned up: adding ``06_spec.json`` pushed the result
    to ``07_design.json`` and left the old ``06_design.json`` sitting in the
    folder, which is the same two-runs-in-one-folder problem by another route.
    Everything not matching that convention is left alone.
    """
    # Serialized before the folder is touched, so a document that will not
    # turn into JSON fails here rather than halfway through a bundle.  The run
    # record is built here for the same reason and one more: it REFUSES a
    # caller key that would shadow a measured version, and that refusal has to
    # land before anything is written.
    request_text = None if request is None else json.dumps(request, indent=2) + "\n"
    run_text = None if run is None else run_json(run) + "\n"

    root = Path(out_dir).resolve()
    if root.is_dir() and any(root.iterdir()):
        if not overwrite:
            raise FileExistsError(
                f"{root} is not empty; pass overwrite=True to replace the "
                "bundle in it"
            )
        ours = {root / name for name in BUNDLE_PATHS.values()}
        ours |= {path for path in root.iterdir() if _NUMBERED.match(path.name)}
        for stale in sorted(ours):
            if stale.is_dir():
                shutil.rmtree(stale)
            elif stale.exists():
                stale.unlink()
    for key in ("primers", "pcrs", "parts", "constructs"):
        (root / BUNDLE_PATHS[key]).mkdir(parents=True, exist_ok=True)

    written = {"layout", "summary", "primers", "pcrs", "parts", "constructs", "design"}
    (root / BUNDLE_PATHS["design"]).write_text(
        json.dumps(result.to_dict(), indent=2) + "\n"
    )
    if result.spec is not None:
        (root / BUNDLE_PATHS["spec"]).write_text(spec_json(result.spec) + "\n")
        written.add("spec")
    if request_text is not None:
        (root / BUNDLE_PATHS["request"]).write_text(request_text)
        written.add("request")
    if run_text is not None:
        (root / BUNDLE_PATHS["run"]).write_text(run_text)
        written.add("run")

    primers = root / BUNDLE_PATHS["primers"]
    _frame(result, "primer").to_csv(primers / "primers.csv", index=False)
    _frame(result, "primer", full=True).to_csv(primers / "analysis.csv", index=False)
    _frame(result, "uses").to_csv(primers / "uses.csv", index=False)
    _frame(result, "pcr").to_csv(root / BUNDLE_PATHS["pcrs"] / "pcrs.csv", index=False)
    _frame(result, "worklist").to_csv(
        root / BUNDLE_PATHS["pcrs"] / "worklist_by_ta.csv", index=False
    )
    _parts_frame(result, part_info).to_csv(
        root / BUNDLE_PATHS["parts"] / "parts.csv", index=False
    )

    # The order, in the two shapes a vendor's form takes.  Built from the same
    # `result.primers` the sheet above is, so the three cannot disagree about
    # what is in the order -- which is asserted on SEQUENCES, because an oligo
    # is identified by its sequence and the right count of the wrong molecules
    # is exactly what a count would pass.
    idt.tube_frame(result).to_csv(primers / "oligo_order_idt.csv", index=False)
    (primers / "primers.fasta").write_text(idt.primers_fasta(result))
    plates = idt.plate_frames(result)
    omit = ()
    if plates:
        with pd_excel_writer(primers / "oligo_order_idt_plate.xlsx") as writer:
            for index, plate in enumerate(plates, start=1):
                plate.to_excel(writer, sheet_name=f"Plate{index}", index=False)
    else:
        # No primers, so nothing to order.  The entry is still written -- the
        # primer tables live in it -- so the one missing line is dropped from
        # the map rather than the whole shelf.
        omit = (f"{BUNDLE_PATHS['primers']}/oligo_order_idt_plate.xlsx",)

    (root / BUNDLE_PATHS["layout"]).write_text(bundle_layout(written, omit=omit))

    designs = list(result.designs)
    constructs = root / BUNDLE_PATHS["constructs"]
    _frame(result, "overview").to_csv(constructs / "constructs.csv", index=False)
    write_genbank_files(designs, constructs)
    (constructs / ALL_CONSTRUCTS).write_text(constructs_genbank(designs))
    figures = write_assembly_figures(designs, constructs)
    # A blank line and a rule between figures, or they run into each other.
    joined = ("\n" + "-" * 72 + "\n\n").join(
        path.read_text().rstrip() + "\n" for path in figures
    )
    (constructs / "assemblies.txt").write_text(joined + "\n" if joined else "")
    for path in figures:
        path.unlink()

    (root / BUNDLE_PATHS["summary"]).write_text(build_summary(result, part_info, run=run))
    return root


def write_bundle_zip(
    result,
    zip_path: Union[str, Path],
    *,
    part_info: Optional[Mapping[str, Mapping[str, object]]] = None,
    request: Optional[Mapping[str, object]] = None,
    run: Optional[Mapping[str, object]] = None,
) -> Path:
    """:func:`write_bundle` into a temporary folder, zipped, for handing over.

    What a backend returns: one file, with the folder structure intact inside
    it.  The folder is built on disk first rather than streamed into the
    archive, so a half-written bundle is a half-written temporary directory and
    never a corrupt download.  That directory is new every time, so there is
    no ``overwrite`` here: an existing *archive* is replaced outright.
    """
    target = Path(zip_path).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    stem = target.name[:-4] if target.name.endswith(".zip") else target.name
    with tempfile.TemporaryDirectory() as scratch:
        folder = write_bundle(
            result,
            Path(scratch) / stem,
            part_info=part_info,
            request=request,
            run=run,
        )
        with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as archive:
            for path in sorted(folder.rglob("*")):
                if path.is_file():
                    archive.write(path, Path(stem) / path.relative_to(folder))
    return target


def pd_excel_writer(path: Path):
    """``pd.ExcelWriter`` on openpyxl, named so the one import has one home.

    A workbook is the only binary this module writes, and it is written because
    IDT's plate upload takes Excel and not CSV.
    """
    import pandas as pd

    return pd.ExcelWriter(path, engine="openpyxl")


def _frame(result, kind: str, full: bool = False):
    """A table, imported here so :mod:`.report` can import this module freely."""
    from teemi.design.combinatorial.report import (
        pcr_dataframe,
        pcr_worklist_dataframe,
        primer_dataframe,
        primer_uses_dataframe,
        variant_overview_dataframe,
    )

    if kind == "uses":
        return primer_uses_dataframe(result)
    if kind == "worklist":
        return pcr_worklist_dataframe(result)
    if kind == "overview":
        return variant_overview_dataframe(result)
    build = primer_dataframe if kind == "primer" else pcr_dataframe
    return build(result, full=full)


def _parts_frame(result, part_info: Optional[Mapping[str, Mapping[str, object]]]):
    """One row per part the build needs, with how hard it is working.

    ``times_used`` is what matters for sourcing: a part at a combinatorial
    position is amplified once per distinct neighbour pair, not once per
    variant, so the number of reactions it appears in is lower than the number
    of constructs containing it.  Both are given, because the first is how many
    PCRs to set up and the second is how much of the library depends on it.
    """
    import pandas as pd

    payload = result.to_dict()
    footprints = {f["part_id"]: f for f in payload["footprints"]}
    reactions = Counter(pcr["part_id"] for pcr in payload["pcrs"])
    in_variants = Counter(
        part_id
        for variant in payload["variants"]
        for part_id in variant["fragment_ids"]
    )
    positions: Dict[str, str] = {}
    lengths: Dict[str, int] = {}
    if result.spec is not None:
        for position in result.spec.positions:
            for part in position.parts:
                positions[part.id] = position.name or position.id
        for pad in result.spec.pads:
            positions[pad.id] = "pad"

    rows = []
    for part_id in sorted(set(reactions) | set(in_variants) | set(footprints)):
        row: Dict[str, object] = {
            "part_id": part_id,
            "position": positions.get(part_id, ""),
            "reactions": reactions.get(part_id, 0),
            "in_constructs": in_variants.get(part_id, 0),
        }
        footprint = footprints.get(part_id)
        if footprint:
            # Named for the orientation, like every other column of DNA here:
            # these are the two primers' footprints, in each primer's own
            # 5'-3' direction, not two reads of the same strand.
            row["forward_footprint(5'-3')"] = footprint["forward_footprint"]
            row["reverse_footprint(5'-3')"] = footprint["reverse_footprint"]
        row.update(dict((part_info or {}).get(part_id, {})))
        rows.append(row)
    return pd.DataFrame(rows)


def build_summary(
    result,
    part_info: Optional[Mapping[str, Mapping[str, object]]] = None,
    *,
    run: Optional[Mapping[str, object]] = None,
) -> str:
    """The page to read first: what this build is, and what to look at twice.

    Written rather than tabulated because the useful facts are not one shape --
    a count of reactions, a reason a reaction is flagged, how often a part is
    reused.  Everything here is derived from the result, so it cannot drift from
    the tables beside it.

    ``run``'s date is printed under the heading.  The archive's filename
    carries it too, but a filename is lost the moment someone unzips the folder
    or renames it, and this is the page they actually read.

    It is read from ONE key, ``built_at``, and no other spelling is looked for:
    a caller passing ``timestamp`` or ``date`` gets no line and no complaint,
    because every other key in ``run`` is the caller's to name and this page
    cannot tell a misspelling from a key it was never meant to print.
    """
    if result.spec is None:
        return "no design: " + "; ".join(d.message for d in result.diagnostics) + "\n"

    payload = result.to_dict()
    primers, pcrs = payload["primers"], payload["pcrs"]
    flagged = [p for p in pcrs if p["flagged_by"]]
    causes = Counter(
        cause for pcr in flagged for cause in pcr["flagged_by"].split(",") if cause
    )
    shared = [p for p in primers if p["n_uses"] > 1]
    blocking = [d for d in result.diagnostics if d.severity in ("error", "critical")]
    unverified = [v for v in result.variants if not v.verification.ok]

    lines = [
        result.spec.design.name,
        "=" * len(result.spec.design.name),
    ]
    built_at = (run or {}).get("built_at")
    if built_at:
        lines.append(f"built {built_at}")
    lines += [
        "",
        result.summary(),
        "",
        "WHAT TO ORDER",
        f"  {len(primers)} distinct oligos, "
        f"{sum(p['length_nt'] for p in primers)} bases, "
        f"{result.price:.0f} {result.spec.objective.currency}",
    ]
    if primers:
        plates = idt.plate_count(len(primers))
        where = f"{BUNDLE_PATHS['primers']}/oligo_order_idt.csv"
        lines.append(
            f"  as tubes: {where}  "
            f"({idt.DEFAULT_SCALE}, {idt.DEFAULT_PURIFICATION})"
        )
        if len(primers) < idt.IDT_PLATE_MINIMUM:
            # Said with both numbers rather than as advice, the way a breached
            # bound is reported: the plate file is still written, because a
            # person with 16 oligos and no plate file cannot tell a too-small
            # order from a broken feature.
            lines.append(
                f"  as a plate: {len(primers)} oligos is under IDT's "
                f"{idt.IDT_PLATE_MINIMUM} for a 96-well plate, so this order "
                "goes as tubes"
            )
        else:
            lines.append(
                f"  as a plate: {plates} "
                f"{'plate' if plates == 1 else 'plates'}, filled down the "
                f"columns -- {BUNDLE_PATHS['primers']}/oligo_order_idt_plate.xlsx"
            )
    if shared:
        lines.append(
            f"  {len(shared)} of them are shared, saving "
            f"{sum(p['n_uses'] - 1 for p in shared)} oligos:"
        )
        for primer in sorted(shared, key=lambda p: -p["n_uses"])[:10]:
            lines.append(
                f"     {primer['primer']}  used {primer['n_uses']}x  "
                f"{primer['length_nt']} nt  anneals to {primer['part_id']}"
            )
    longest = max(primers, key=lambda p: p["length_nt"]) if primers else None
    if longest:
        lines.append(
            f"  longest is {longest['primer']} at {longest['length_nt']} nt"
            + (
                "  (past 60 nt, so its structure numbers describe the 3' end only)"
                if longest["length_nt"] > 60
                else ""
            )
        )

    lines += [
        "",
        "WHAT TO RUN",
        f"  {len(pcrs)} reactions for {len(result.variants)} constructs",
        f"  {len(pcrs) - len(flagged)} clean, {len(flagged)} flagged",
    ]
    for cause, count in causes.most_common():
        lines.append(f"     {count:4d}  {cause}")
    if flagged:
        lines.append("  flagged reactions:")
        for pcr in flagged[:20]:
            lines.append(
                f"     {pcr['part_id'][:34]:34s} {pcr['forward_primer']}/"
                f"{pcr['reverse_primer']}  {pcr['flagged_by']}"
            )
        if len(flagged) > 20:
            lines.append(
                f"     ... and {len(flagged) - 20} more, see "
                f"{BUNDLE_PATHS['pcrs']}/pcrs.csv"
            )
        lines.append(
            "  A flag is a prediction, not a refusal: primer3's forecast that a"
        )
        lines.append(
            "  structure may form.  Worth a look before ordering, not a reason to stop."
        )

    folded = [p for p in primers if (p.get("hairpin_folded_at_ta") or 0) > 0.25]
    blocked = [p for p in primers if p.get("hairpin_3prime_paired")]
    if folded or blocked:
        lines += ["", "WORTH A SECOND LOOK"]
        for primer in sorted(folded, key=lambda p: -p["hairpin_folded_at_ta"])[:10]:
            lines.append(
                f"  {primer['primer']} is {primer['hairpin_folded_at_ta']:.0%} folded "
                f"at its annealing temperature"
            )
        for primer in blocked[:10]:
            lines.append(
                f"  {primer['primer']} has a hairpin stem over its 3' end, which a "
                "polymerase cannot extend through"
            )

    parts = _parts_frame(result, part_info)
    lines += [
        "",
        "WHAT TO SOURCE",
        f"  {len(parts)} parts",
    ]
    if not parts.empty:
        busiest = parts.sort_values("reactions", ascending=False).head(5)
        for _, row in busiest.iterrows():
            lines.append(
                f"     {str(row['part_id'])[:30]:30s} {int(row['reactions']):3d} "
                f"reactions, in {int(row['in_constructs'])} constructs"
            )
        missing = [
            part_id
            for part_id in parts["part_id"]
            if not (part_info or {}).get(part_id)
        ]
        if part_info is not None and missing:
            lines.append(
                f"  no provenance given for {len(missing)}: "
                f"{', '.join(str(m) for m in missing[:6])}"
            )

    lines += [
        "",
        "WHAT WAS CHECKED",
        f"  {len(result.variants)} constructs simulated, "
        f"{len(result.variants) - len(unverified)} assemble to the intended sequence",
    ]
    for variant in unverified[:10]:
        lines.append(
            f"     {variant.label}  {variant.id[:34]:34s} "
            f"{variant.verification.verdict}: "
            f"{' '.join(variant.verification.message.split())[-70:]}"
        )
    lines.append(
        f"  {len(list(result.designs))} annotated maps written to "
        f"{BUNDLE_PATHS['constructs']}/"
    )
    if unverified:
        lines.append(
            f"  {len(unverified)} not written: an unverified assembly has no single "
            "product to draw"
        )
    lines += [
        "",
        f"  blocking errors: {len(blocking)}",
    ]
    for diagnostic in blocking[:10]:
        lines.append(f"     {' '.join(diagnostic.message.split())[:96]}")
    if not blocking:
        lines.append("  nothing here cannot be built as specified.")
    return "\n".join(lines) + "\n"
