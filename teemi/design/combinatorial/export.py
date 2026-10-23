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
import json
import shutil
import tempfile
import zipfile
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from teemi.design.combinatorial.report import VariantPlan

__all__ = [
    "write_genbank_files",
    "write_assembly_figures",
    "write_bundle",
    "write_bundle_zip",
    "build_summary",
    "BUNDLE_LAYOUT",
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

#: What a bundle contains, and in what order a person should read it.
#: The order to read the bundle in, which is why the entries are numbered: a
#: folder listing is alphabetical and would otherwise open on ``constructs/``.
BUNDLE_ORDER = (
    ("summary", "01_SUMMARY.txt"),
    ("primers", "02_primers"),
    ("pcrs", "03_pcrs"),
    ("parts", "04_parts"),
    ("constructs", "05_constructs"),
    ("design", "06_design.json"),
    ("layout", "00_LAYOUT.txt"),
)

#: Where each thing is, by the key :data:`BUNDLE_ORDER` gives it.
BUNDLE_PATHS = dict(BUNDLE_ORDER)

BUNDLE_LAYOUT = f"""\
{BUNDLE_PATHS["summary"]:<30} read this first: what to order, run, look at twice
{BUNDLE_PATHS["primers"] + "/primers.csv":<30} the order sheet, one row per distinct molecule
{BUNDLE_PATHS["primers"] + "/analysis.csv":<30} every column the primer table has, not just the readable ones
{BUNDLE_PATHS["pcrs"] + "/pcrs.csv":<30} one row per reaction
{BUNDLE_PATHS["parts"] + "/parts.csv":<30} what to source, how often each is used, and what is known
{BUNDLE_PATHS["constructs"] + "/NNN_*.gb":<30} one annotated map per verified construct, numbered
{BUNDLE_PATHS["constructs"] + "/assemblies.txt":<30} every assembly figure, in one file
{BUNDLE_PATHS["design"]:<30} the whole result, machine-readable; re-runs the design
"""


def write_bundle(
    result,
    out_dir: Union[str, Path],
    *,
    part_info: Optional[Mapping[str, Mapping[str, object]]] = None,
    overwrite: bool = False,
) -> Path:
    """Write the whole build to one folder and return it.

    Everything a person needs to go and do the work, in the order they need it:
    what to order, what to run, what to source, what the constructs should come
    out as, and the summary that says which of it is worth looking at twice.
    The entries are numbered in that order, because a folder listing is
    alphabetical and would otherwise open on ``constructs/``.
    :data:`BUNDLE_LAYOUT` is the map, and is written into the folder too.

    ``part_info`` is provenance this package cannot know -- an Addgene number, a
    DOI, which plasmid a part came from.  The design carries sequences and ids;
    where a part can be *obtained* lives in whatever library the caller
    resolved it from, so the caller passes it in, keyed by part id.  Any keys
    are accepted and become columns, so a different library's fields need no
    change here.

    Only verified constructs are written, because an unverified one has no
    single product to draw.  The summary says how many were left out.

    A non-empty folder is refused unless ``overwrite`` says otherwise, because
    a second run writing into the first leaves both: drop a part and its maps
    stay, and someone reads half of one design and half of another with
    nothing to tell them.  With ``overwrite`` only the entries this function
    writes are removed -- by the names in :data:`BUNDLE_PATHS`, so whatever
    else is in the folder is left alone.
    """
    root = Path(out_dir).resolve()
    if root.is_dir() and any(root.iterdir()):
        if not overwrite:
            raise FileExistsError(
                f"{root} is not empty; pass overwrite=True to replace the "
                "bundle in it"
            )
        for name in BUNDLE_PATHS.values():
            stale = root / name
            if stale.is_dir():
                shutil.rmtree(stale)
            elif stale.exists():
                stale.unlink()
    for key in ("primers", "pcrs", "parts", "constructs"):
        (root / BUNDLE_PATHS[key]).mkdir(parents=True, exist_ok=True)

    (root / BUNDLE_PATHS["layout"]).write_text(BUNDLE_LAYOUT)
    (root / BUNDLE_PATHS["design"]).write_text(
        json.dumps(result.to_dict(), indent=2) + "\n"
    )

    primers = root / BUNDLE_PATHS["primers"]
    _frame(result, "primer").to_csv(primers / "primers.csv", index=False)
    _frame(result, "primer", full=True).to_csv(primers / "analysis.csv", index=False)
    _frame(result, "pcr").to_csv(root / BUNDLE_PATHS["pcrs"] / "pcrs.csv", index=False)
    _parts_frame(result, part_info).to_csv(
        root / BUNDLE_PATHS["parts"] / "parts.csv", index=False
    )

    designs = list(result.designs)
    constructs = root / BUNDLE_PATHS["constructs"]
    write_genbank_files(designs, constructs)
    figures = write_assembly_figures(designs, constructs)
    # A blank line and a rule between figures, or they run into each other.
    joined = ("\n" + "-" * 72 + "\n\n").join(
        path.read_text().rstrip() + "\n" for path in figures
    )
    (constructs / "assemblies.txt").write_text(joined + "\n" if joined else "")
    for path in figures:
        path.unlink()

    (root / BUNDLE_PATHS["summary"]).write_text(build_summary(result, part_info))
    return root


def write_bundle_zip(
    result,
    zip_path: Union[str, Path],
    *,
    part_info: Optional[Mapping[str, Mapping[str, object]]] = None,
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
        folder = write_bundle(result, Path(scratch) / stem, part_info=part_info)
        with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as archive:
            for path in sorted(folder.rglob("*")):
                if path.is_file():
                    archive.write(path, Path(stem) / path.relative_to(folder))
    return target


def _frame(result, kind: str, full: bool = False):
    """A table, imported here so :mod:`.report` can import this module freely."""
    from teemi.design.combinatorial.report import pcr_dataframe, primer_dataframe

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
            row["forward_footprint"] = footprint["forward_footprint"]
            row["reverse_footprint"] = footprint["reverse_footprint"]
        row.update(dict((part_info or {}).get(part_id, {})))
        rows.append(row)
    return pd.DataFrame(rows)


def build_summary(
    result, part_info: Optional[Mapping[str, Mapping[str, object]]] = None
) -> str:
    """The page to read first: what this build is, and what to look at twice.

    Written rather than tabulated because the useful facts are not one shape --
    a count of reactions, a reason a reaction is flagged, how often a part is
    reused.  Everything here is derived from the result, so it cannot drift from
    the tables beside it.
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
        "",
        result.summary(),
        "",
        "WHAT TO ORDER",
        f"  {len(primers)} distinct oligos, "
        f"{sum(p['length_nt'] for p in primers)} bases, "
        f"{result.price:.0f} {result.spec.objective.currency}",
    ]
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
