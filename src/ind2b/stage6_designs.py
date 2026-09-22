"""Stage 6 - read the binders a BindCraft2 campaign actually designed.

This stage reads a campaign's output folder. It does not search for other
people's binders and it does not synthesise anything: if no campaign has been
run, this stage reports that no designs exist, and the report says so.

Layout and column conventions below follow the upstream BindCraft2 output
documentation rather than assumption; the file names and column names are the
ones BindCraft2 writes:

``<project_folder>/``
    ``3_Ranked/!_Ranked.csv``       every accepted sequence, best-first
    ``3_Ranked/<design>.cif``       the accepted predicted complex
    ``2_Refolded/!_Refolded.csv``   every scored candidate, with failed filters
    ``1_Trajectories/!_Trajectories.csv``  one row per design attempt
    ``summary.csv``                 counts and metric summaries
    ``campaign_metadata.json``      resolved settings and input provenance

What this stage is careful about, because each is a documented way to
misread these tables:

*Acceptance is not affinity.* An accepted design passed the configured
computational filters. It does not establish affinity, specificity or
experimental stability, and this stage propagates that wording rather than
letting "accepted" imply a working binder.

*A blank is not a zero.* Missing readings mean the measurement was not taken;
``failed_filters`` distinguishes "not measured" from a poor score. Nothing
here fills a blank with a number.

*Rank is not probability.* ``rank`` is a position in one particular ranking,
not a probability of experimental success.

*Some attempts ran an easier task.* The ``autotuned`` column records settings
that differed from the campaign's own. Where it names ``initial_guess``,
``target_flexibility``, ``validation_model`` or a raised ``design_recycles``,
that attempt ran on the desperation ladder - against an easier problem than
the campaign asked for. This stage flags those designs so a reader does not
compare them with the rest as equals.

*Designed interface residues are in the written structure's numbering*, which
is the same author numbering the stage 5 specification requested hotspots in,
so the two are directly comparable.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any, Iterable, Sequence

import pandas as pd

from . import __version__
from .config import DESIGNS_FILE, InterfaceParams

log = logging.getLogger("ind2b.stage6_designs")

SCHEMA_VERSION = 1
RECORD_KIND = "bindcraft2_designs"

RANKED_CSV = "3_Ranked/!_Ranked.csv"
REFOLDED_CSV = "2_Refolded/!_Refolded.csv"
TRAJECTORIES_CSV = "1_Trajectories/!_Trajectories.csv"
SUMMARY_CSV = "summary.csv"
METADATA_JSON = "campaign_metadata.json"

# Settings whose appearance in `autotuned` means the attempt ran against an
# easier task than the campaign specified (the "desperation ladder").
LADDER_SETTINGS = ("initial_guess", "target_flexibility", "validation_model",
                   "design_recycles")

ACCEPTANCE_MEANING = (
    "Acceptance means the campaign's configured computational checks passed. "
    "It does not establish affinity, specificity or experimental stability. "
    "No design here has been tested in a laboratory."
)


def _read_csv(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        return pd.read_csv(path)
    except Exception as exc:  # noqa: BLE001
        log.warning("could not read %s: %s", path, exc)
        return None


def _read_json(path: Path) -> dict | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except Exception as exc:  # noqa: BLE001
        log.warning("could not read %s: %s", path, exc)
        return None


def discover_campaigns(
    search_roots: Iterable[Path],
    *,
    max_depth: int = 3,
) -> list[Path]:
    """Find BindCraft2 project folders under the given roots.

    A campaign folder is identified by the records BindCraft2 writes, not by
    its name: a folder holding ``campaign_metadata.json`` or a ranked table.
    Stage folders appear only once their first result is written, so a
    campaign that has run but accepted nothing is still found.
    """
    found: list[Path] = []
    for root in search_roots:
        root = Path(root)
        if not root.is_dir():
            continue
        candidates = [root, *(p for p in root.rglob("*")
                              if p.is_dir()
                              and len(p.relative_to(root).parts) <= max_depth)]
        for folder in candidates:
            if ((folder / METADATA_JSON).exists()
                    or (folder / RANKED_CSV).exists()
                    or (folder / TRAJECTORIES_CSV).exists()):
                if folder not in found:
                    found.append(folder)
    return found


def _ladder_flags(autotuned: Any) -> list[str]:
    """Which desperation-ladder settings an attempt used, if any."""
    if autotuned is None or (isinstance(autotuned, float) and pd.isna(autotuned)):
        return []
    text = str(autotuned)
    return [s for s in LADDER_SETTINGS if re.search(rf"\b{s}\b", text)]


def parse_interface_residues(value: Any) -> list[str]:
    """Residue tokens from an ``Interface_*_Residues`` cell.

    Binder chains are separated by ``/`` in the same order as the sequence;
    residues within a chain are comma-separated. A blank means the residues
    were not recorded - returned as an empty list, never as zero contacts.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    out: list[str] = []
    for chunk in str(value).replace("/", ",").split(","):
        tok = chunk.strip()
        if tok:
            out.append(tok)
    return out


def epitope_engagement(
    design_interface_target: Any,
    requested_hotspots: Sequence[str],
) -> dict[str, Any]:
    """Did the design engage the hotspots the specification asked for?

    Both sides are in the written structure's author numbering, which is what
    the stage 5 specification requested, so they compare directly. A design
    that binds elsewhere on the target is a real outcome and is reported as
    such rather than smoothed over.
    """
    got = {t for t in parse_interface_residues(design_interface_target)}
    want = {str(h).strip() for h in requested_hotspots if str(h).strip()}
    if not got:
        return {"comparable": False,
                "reason": "the design's target-side interface residues were not "
                          "recorded in the ranked table"}
    if not want:
        return {"comparable": False,
                "reason": "no requested hotspots were available to compare against"}
    hit = sorted(want & got, key=lambda s: (len(s), s))
    return {
        "comparable": True,
        "n_requested": len(want),
        "n_hotspots_engaged": len(hit),
        "fraction_hotspots_engaged": round(len(hit) / len(want), 3),
        "hotspots_engaged": hit,
        "hotspots_missed": sorted(want - got, key=lambda s: (len(s), s)),
        "n_interface_residues": len(got),
    }


def read_campaign(
    folder: Path,
    *,
    requested_hotspots: Sequence[str] = (),
    max_designs: int = 25,
) -> dict[str, Any]:
    """Read one campaign folder into a record of what it actually produced."""
    folder = Path(folder)
    meta = _read_json(folder / METADATA_JSON)
    ranked = _read_csv(folder / RANKED_CSV)
    refolded = _read_csv(folder / REFOLDED_CSV)
    trajectories = _read_csv(folder / TRAJECTORIES_CSV)
    summary = _read_csv(folder / SUMMARY_CSV)

    record: dict[str, Any] = {
        "project_folder": str(folder),
        "campaign_name": (meta or {}).get("campaign_name"),
        "has_metadata": meta is not None,
        # A stage folder exists only once it has written a result, so absence
        # is informative: it says the campaign never got that far.
        "stages_present": {
            "trajectories": (folder / TRAJECTORIES_CSV).exists(),
            "refolded": (folder / REFOLDED_CSV).exists(),
            "ranked": (folder / RANKED_CSV).exists(),
        },
        "n_attempts": None if trajectories is None else len(trajectories),
        "n_candidates_scored": None if refolded is None else len(refolded),
        "n_accepted": None if ranked is None else len(ranked),
        "acceptance_meaning": ACCEPTANCE_MEANING,
        "designs": [],
    }
    if meta:
        resolved = meta.get("resolved") or {}
        record["settings"] = {
            k: meta.get(k) for k in
            ("target", "modality", "binder_lengths", "number_of_final_designs")
            if k in meta
        }
        record["resolved"] = {k: resolved.get(k) for k in list(resolved)[:12]}
        record["source_revision"] = meta.get("source_revision")

    # Why nothing was accepted is a question the tables can answer.
    if refolded is not None and "failed_filters" in refolded.columns:
        tally: dict[str, int] = {}
        for cell in refolded["failed_filters"].dropna():
            for f in str(cell).split(","):
                f = f.strip()
                if f:
                    tally[f] = tally.get(f, 0) + 1
        record["failed_filter_tally"] = dict(
            sorted(tally.items(), key=lambda kv: -kv[1])[:15])
    if trajectories is not None and "terminated" in trajectories.columns:
        term = trajectories["terminated"].dropna().astype(str)
        record["termination_tally"] = (
            term[term.str.strip() != ""].value_counts().head(10).to_dict())

    if ranked is None or ranked.empty:
        record["designs_note"] = (
            "No accepted designs. Where a campaign ran, inspect the refolded "
            "table's failed filters and the trajectory terminations above; "
            "absence of an accepted design is not a statement about the target."
        )
        return record

    record["ranked_columns"] = list(ranked.columns)
    for _, row in ranked.head(max_designs).iterrows():
        name = row.get("design")
        ladder = _ladder_flags(row.get("autotuned"))
        seq = row.get("Binder_Sequence")
        seq = None if (isinstance(seq, float) and pd.isna(seq)) else seq
        cif = None
        if isinstance(name, str):
            for cand in sorted(folder.glob(f"3_Ranked/{name}*.cif")):
                if "_monomer" not in cand.name:
                    cif = str(cand.relative_to(folder))
                    break
        design = {
            "design": name,
            "hash": row.get("hash"),
            "rank": row.get("rank"),
            "binder_sequence": seq,
            "binder_chains": None if not isinstance(seq, str) else seq.count("/") + 1,
            "binder_length": None if not isinstance(seq, str)
            else len(seq.replace("/", "")),
            "structure_file": cif,
            "interface_target_residues": parse_interface_residues(
                row.get("Interface_Target_Residues")),
            "interface_binder_residues": parse_interface_residues(
                row.get("Interface_Binder_Residues")),
            "failed_filters": (None if pd.isna(row.get("failed_filters"))
                               else str(row.get("failed_filters"))),
            "ran_on_desperation_ladder": bool(ladder),
            "ladder_settings": ladder,
            "autotuned": (None if pd.isna(row.get("autotuned"))
                          else str(row.get("autotuned"))),
        }
        # Every remaining column is a metric; keep them all rather than
        # choosing a favourite, and keep blanks as None rather than zero.
        design["metrics"] = {
            c: (None if pd.isna(row[c]) else
                (float(row[c]) if isinstance(row[c], (int, float)) else str(row[c])))
            for c in ranked.columns
            if c not in ("design", "hash", "rank", "Binder_Sequence",
                         "Interface_Target_Residues", "Interface_Binder_Residues",
                         "failed_filters", "autotuned")
        }
        design["vs_requested_epitope"] = epitope_engagement(
            row.get("Interface_Target_Residues"), requested_hotspots)
        record["designs"].append(design)

    record["n_on_ladder"] = sum(
        1 for d in record["designs"] if d["ran_on_desperation_ladder"])
    return record


def run(
    run_dir: Path,
    *,
    project_folders: Sequence[Path] | None = None,
    specs_manifest: dict | None = None,
    params: InterfaceParams | None = None,
    max_designs: int = 25,
) -> dict[str, Any]:
    """Collect every campaign's designs, matched to the spec that requested them."""
    run_dir = Path(run_dir)
    hotspots_by_target: dict[str, list[str]] = {}
    spec_meta: dict[str, dict] = {}
    for spec in (specs_manifest or {}).get("specs") or []:
        if spec.get("status") != "written":
            continue
        name = spec.get("name")
        hs = [h for h in str(spec.get("hotspots_auth") or "").split(",") if h]
        if name:
            hotspots_by_target[name] = hs
            spec_meta[name] = spec

    roots = [Path(p) for p in (project_folders or [])] or [
        run_dir / "results", run_dir / "campaigns", run_dir / "stage6_designs"]
    folders = discover_campaigns(roots)

    campaigns = []
    for folder in folders:
        # Match the campaign to its specification by name where possible, so
        # the requested hotspots compared against are the ones this campaign
        # was actually given.
        meta = _read_json(folder / METADATA_JSON) or {}
        key = meta.get("target") or folder.name
        hs = hotspots_by_target.get(key) or hotspots_by_target.get(folder.name) or []
        rec = read_campaign(folder, requested_hotspots=hs, max_designs=max_designs)
        rec["matched_spec"] = key if key in hotspots_by_target else None
        rec["requested_hotspots"] = hs
        if rec["matched_spec"] is None:
            rec["match_note"] = (
                "This campaign could not be matched to a stage 5 specification, "
                "so no requested-hotspot comparison is made for it.")
        campaigns.append(rec)

    total = sum(c["n_accepted"] or 0 for c in campaigns)
    return {
        "schema_version": SCHEMA_VERSION,
        "record_kind": RECORD_KIND,
        "package_version": __version__,
        "searched_roots": [str(r) for r in roots],
        "n_campaigns": len(campaigns),
        "n_accepted_total": total,
        "acceptance_meaning": ACCEPTANCE_MEANING,
        "no_campaign_note": (
            "No BindCraft2 campaign output was found under the searched roots, "
            "so this pipeline has produced no binder. Stage 5 wrote "
            "specifications; running them is a separate, explicit step on a GPU "
            "host." if not campaigns else None
        ),
        "campaigns": campaigns,
    }


def write(record: dict[str, Any], run_dir: Path) -> Path:
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    path = run_dir / DESIGNS_FILE
    path.write_text(json.dumps(record, indent=2))
    return path


def read(run_dir: Path) -> dict[str, Any]:
    path = Path(run_dir) / DESIGNS_FILE
    if not path.exists():
        raise FileNotFoundError(f"{path} not found - run stage 6 first")
    return json.loads(path.read_text())
