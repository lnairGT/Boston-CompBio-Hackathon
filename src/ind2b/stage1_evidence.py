"""Stage 1 - assemble the evidence layers for a disease's candidate targets.

Reads ``stage0_disease.json`` and writes ``stage1_evidence.json``: one record
per candidate target holding the raw values stage 2 will score. No scoring,
no ranking and no filtering-out happens here - this stage only gathers, so
that re-weighting later never requires re-fetching.

Fetch ordering
--------------
The association pull and the bulk target detail are two requests for the whole
pool. The remaining layers (baseline expression, known drugs, interaction
partners) are one request *per target*, so they are fetched only for targets
that pass the surface-accessibility check - a de novo mini-binder cannot reach
an intracellular target without a delivery strategy, and those targets will be
dropped by stage 3 regardless. Pass ``all_layers=True`` to fetch them for
every target anyway.

Surface accessibility is recorded as a fact plus the locations that matched,
never as a silent exclusion: the record for an intracellular target is still
written, with ``surface_accessible: false``.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from . import __version__
from .config import STAGE_FILES, SURFACE_TERMS
from .sources import chembl
from .sources import opentargets as ot
from .sources import trials
from .sources import uniprot as up

log = logging.getLogger("ind2b.stage1")

SCHEMA_VERSION = 1
_PAGE_SIZE = 100
_DETAIL_CHUNK = 50


def fetch_pool(disease_id: str, n_targets: int) -> tuple[list[dict], int]:
    """Top ``n_targets`` associated targets, paged, with the total count."""
    rows: list[dict] = []
    index = 0
    total = 0
    while len(rows) < n_targets:
        block = ot.associated_targets(disease_id, size=_PAGE_SIZE, index=index)
        total = block.get("count", 0)
        page = block.get("rows") or []
        if not page:
            break
        rows.extend(page)
        index += 1
        if len(rows) >= total:
            break
    return rows[:n_targets], total


def surface_status(locations: Iterable[dict] | None) -> tuple[bool, list[dict]]:
    """Decide whether a binder could reach this target from outside the cell.

    Matches ``config.SURFACE_TERMS`` against the free-text location and
    ``labelSL`` fields, returning the matches so the call can be audited.
    """
    matched: list[dict] = []
    for loc in locations or []:
        text = " ".join(
            str(loc.get(k) or "") for k in ("location", "labelSL", "termSL")
        ).lower()
        for term in SURFACE_TERMS:
            if term in text:
                matched.append(
                    {
                        "term": term,
                        "location": loc.get("location"),
                        "source": loc.get("source"),
                    }
                )
                break
    return bool(matched), matched


def _tractability(entries: Iterable[dict] | None) -> dict[str, list[str]]:
    """Group tractability buckets by modality, keeping only positive calls."""
    out: dict[str, list[str]] = {}
    for e in entries or []:
        if e.get("value"):
            out.setdefault(e.get("modality") or "?", []).append(e.get("label") or "?")
    return out


def _expression_specificity(rows: Iterable[dict] | None) -> dict[str, Any]:
    """Summarise upstream expression specificity.

    Open Targets ships ``specificity_score`` per row, so this reports the
    maximum across rows rather than recomputing a specificity measure.
    """
    best: dict[str, Any] | None = None
    n = 0
    for r in rows or []:
        n += 1
        score = r.get("specificity_score")
        if score is None:
            continue
        if best is None or score > best["specificity_score"]:
            best = {
                "specificity_score": score,
                "distribution_score": r.get("distribution_score"),
                "tissue": r.get("tissueBiosampleFromSource"),
                "celltype": r.get("celltypeBiosampleFromSource"),
                "datasource": r.get("datasourceId"),
            }
    return {"n_rows": n, "max_specificity": best}


def _known_drugs(block: dict, top_n: int = 5) -> dict[str, Any]:
    """Drug precedent plus the modality breakdown.

    ``drug_types`` is the modality-precedent signal: a target already engaged
    by an approved antibody is far better evidence that a protein binder can
    work on it than one engaged only by small molecules.
    """
    rows = block.get("rows") or []
    stages = [r.get("maxClinicalStage") for r in rows if r.get("maxClinicalStage")]
    drug_types: dict[str, int] = {}
    for r in rows:
        dt = (r.get("drug") or {}).get("drugType") or "unknown"
        drug_types[dt] = drug_types.get(dt, 0) + 1
    return {
        "count": block.get("count", len(rows)),
        "has_approved": any(s == "APPROVAL" for s in stages),
        "drug_types": drug_types,
        "has_biologic_precedent": any(
            k.lower() in ("antibody", "protein", "enzyme", "oligosaccharide")
            for k in drug_types
        ),
        "top": [
            {
                "drug": (r.get("drug") or {}).get("name"),
                "drug_type": (r.get("drug") or {}).get("drugType"),
                "max_clinical_stage": r.get("maxClinicalStage"),
            }
            for r in rows[:top_n]
        ],
    }


def _interactions(rows: Iterable[dict] | None, top_n: int = 25) -> dict[str, Any]:
    """Deduplicate interaction partners, keeping the best-scoring evidence."""
    by_symbol: dict[str, dict] = {}
    for r in rows or []:
        tb = r.get("targetB") or {}
        symbol = tb.get("approvedSymbol")
        if not symbol:
            continue
        score = r.get("score") or 0.0
        prev = by_symbol.get(symbol)
        if prev is None or score > (prev.get("score") or 0.0):
            by_symbol[symbol] = {
                "symbol": symbol,
                "ensembl_id": tb.get("id"),
                "score": score,
                "source": r.get("sourceDatabase"),
                "n_evidence": r.get("count"),
            }
    partners = sorted(by_symbol.values(), key=lambda d: -(d["score"] or 0.0))
    return {"n_partners": len(partners), "partners": partners[:top_n]}


def _uniprot_accessions(protein_ids: Iterable[dict] | None) -> list[str]:
    """Reviewed (Swiss-Prot) accessions first; these drive the PDB search."""
    swiss = [p["id"] for p in protein_ids or [] if p.get("source") == "uniprot_swissprot"]
    if swiss:
        return swiss
    return [p["id"] for p in protein_ids or [] if str(p.get("source", "")).startswith("uniprot")]


def build_records(
    assoc_rows: list[dict],
    *,
    disease_name: str | None = None,
    all_layers: bool = False,
) -> list[dict]:
    """Merge association rows with bulk detail and the per-target layers.

    ``disease_name`` is the human-readable disease label, needed to query the
    trial registry by condition; without it the clinical-precedent layer is
    skipped rather than guessed.
    """
    ids = [r["target"]["id"] for r in assoc_rows]
    details: dict[str, dict] = {}
    for i in range(0, len(ids), _DETAIL_CHUNK):
        chunk = ids[i : i + _DETAIL_CHUNK]
        for d in ot.targets_details(chunk):
            details[d["id"]] = d
        log.info("bulk detail %d/%d", min(i + _DETAIL_CHUNK, len(ids)), len(ids))

    records: list[dict] = []
    for rank, row in enumerate(assoc_rows, start=1):
        tid = row["target"]["id"]
        det = details.get(tid, {})

        # Accessibility is decided by UniProt topology, not by the location
        # keyword hint: the hint calls KRAS accessible because it is annotated
        # at the cell membrane, though it faces the cytoplasm. Both are
        # recorded so the disagreement is visible.
        hint, hint_matches = surface_status(det.get("subcellularLocations"))
        accessions = _uniprot_accessions(det.get("proteinIds"))
        if accessions:
            access = up.accessibility(accessions[0])
        else:
            access = {
                "accessible": False,
                "mode": None,
                "extracellular_spans": [],
                "reason": "no_uniprot_accession",
            }
        accessible = bool(access["accessible"])

        rec: dict[str, Any] = {
            "rank_by_ot_score": rank,
            "ensembl_id": tid,
            "symbol": row["target"].get("approvedSymbol"),
            "name": row["target"].get("approvedName"),
            "biotype": row["target"].get("biotype"),
            "ot_overall_score": row.get("score"),
            "datatype_scores": {
                d["id"]: d["score"] for d in (row.get("datatypeScores") or [])
            },
            "uniprot": accessions,
            "surface_accessible": accessible,
            "accessibility": {
                "mode": access.get("mode"),
                "extracellular_spans": access.get("extracellular_spans") or [],
                "signal_peptide": access.get("signal_peptide"),
                "transmembrane": access.get("transmembrane") or [],
                "sequence_length": access.get("sequence_length"),
                "reason": access.get("reason"),
                "source": "UniProt topological domain / subcellular location",
            },
            "location_keyword_hint": {
                "hint": hint,
                "agrees_with_topology": hint == accessible,
                "matches": hint_matches,
            },
            "n_locations": len(det.get("subcellularLocations") or []),
            "tractability": _tractability(det.get("tractability")),
            "safety_liabilities": {
                "count": len(det.get("safetyLiabilities") or []),
                "events": sorted(
                    {
                        e.get("event")
                        for e in (det.get("safetyLiabilities") or [])
                        if e.get("event")
                    }
                )[:10],
            },
            "genetic_constraint": det.get("geneticConstraint") or [],
            "target_class": [
                c.get("label") for c in (det.get("targetClass") or []) if c.get("label")
            ],
            "layers_fetched": False,
        }

        if accessible or all_layers:
            rec["expression"] = _expression_specificity(ot.target_expression(tid))
            rec["known_drugs"] = _known_drugs(ot.target_known_drugs(tid))
            rec["interactions"] = _interactions(ot.target_interactions(tid))
            rec["chembl"] = (
                chembl.precedent(accessions[0]) if accessions else {"n_mechanisms": 0}
            )
            if disease_name:
                rec["trials"] = trials.precedent(disease_name, rec["symbol"] or "")
            # Literature support comes from the Open Targets `literature` and
            # `genetic_literature` datatype scores already captured in
            # `datatype_scores`; no separate bibliographic source is queried.
            rec["layers_fetched"] = True
            log.info("layers %s (%s)", rec["symbol"], tid)

        records.append(rec)
    return records


def fetch(
    disease_id: str,
    *,
    disease_name: str | None = None,
    n_targets: int = 100,
    all_layers: bool = False,
) -> dict[str, Any]:
    assoc_rows, total = fetch_pool(disease_id, n_targets)
    log.info("pool: %d of %d associated targets", len(assoc_rows), total)
    records = build_records(
        assoc_rows, disease_name=disease_name, all_layers=all_layers
    )
    n_surface = sum(1 for r in records if r["surface_accessible"])
    return {
        "schema_version": SCHEMA_VERSION,
        "package_version": __version__,
        "fetched_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "disease_id": disease_id,
        "disease_name": disease_name,
        "association_total": total,
        "pool_size": len(records),
        "n_surface_accessible": n_surface,
        "all_layers": all_layers,
        "sources": [
            "Open Targets Platform: associations (incl. literature and "
            "genetic_literature datatype scores), tractability, safety, "
            "genetic constraint, target class, baseline expression, "
            "interaction partners",
            "UniProt: topology-based accessibility and extracellular spans",
            "ChEMBL: mechanism and action-type precedent",
            "ClinicalTrials.gov v2: clinical precedent including stopped trials",
        ],
        "targets": records,
    }


def write(record: dict[str, Any], run_dir: Path) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)
    path = run_dir / STAGE_FILES[1]
    path.write_text(json.dumps(record, indent=2))
    return path


def read(run_dir: Path) -> dict[str, Any]:
    path = run_dir / STAGE_FILES[1]
    if not path.exists():
        raise FileNotFoundError(f"{path} not found - run stage 1 first")
    return json.loads(path.read_text())
