"""Stage 3 - the experimental-complex gate.

A target is admitted only when the PDB holds an experimentally determined
complex of that target with a *mechanistically relevant* partner - one of its
known interaction partners, not an antibody fragment and not a crystallisation
chaperone. Targets with no qualifying complex are dropped here, and the reason
is recorded rather than the target silently disappearing.

Three distinctions this stage has to make, all of them learned from real data:

*Peptide fragments are not epitopes.* EGFR appears in the PDB as a 13-residue
phosphopeptide bound to CBL-C. The entry is a genuine experimental complex of
two real proteins, and useless for designing a binder against the receptor.
Hence ``min_target_chain_residues``.

*Antibody complexes are evidence, not mechanism.* An anti-EGFR Fab complex
shows that an epitope is targetable by a biologic, which is informative, but
blocking a defined interaction needs the native partner. Fab chains usually
carry no UniProt accession at all, so they are detected by description and
classified separately rather than counted as partners.

*A native partner does not guarantee a reachable interface.* CBL-C binds the
intracellular tail of EGFR. This stage therefore carries the extracellular
spans forward with each candidate; stage 4 enforces the topological check
against real residue numbers once the interface has been computed, because
entity sequence length alone cannot say where in the chain the interface lies.
"""

from __future__ import annotations

import json
import logging
import re
from dataclasses import asdict
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from . import __version__
from .config import STAGE_FILES, StructureGate
from .sources import opentargets as ot
from .sources import rcsb

log = logging.getLogger("ind2b.stage3")

SCHEMA_VERSION = 1

# Descriptions indicating an immunoglobulin-derived chain rather than a
# physiological binding partner.
_ANTIBODY_RE = re.compile(
    r"\b(fab|f\(ab\)|antibody|antibodies|immunoglobulin|ig[gahmde]\b|nanobody|"
    r"sybody|scfv|vhh|single[- ]chain (?:variable|fv)|heavy chain|light chain|"
    r"fv fragment|variable domain|monobody|affibody|designed ankyrin|darpin)\b",
    re.IGNORECASE,
)

# Descriptions indicating a crystallisation or expression aid.
_CHAPERONE_RE = re.compile(
    r"\b(lysozyme|maltose[- ]binding|glutathione s-transferase|thioredoxin|"
    r"green fluorescent|ubiquitin-like protease|nanodisc|bril)\b",
    re.IGNORECASE,
)


def _accessions(entity: dict) -> list[str]:
    ci = entity.get("rcsb_polymer_entity_container_identifiers") or {}
    return [
        r["database_accession"]
        for r in (ci.get("reference_sequence_identifiers") or [])
        if r.get("database_name") == "UniProt" and r.get("database_accession")
    ]


def _chains(entity: dict) -> list[str]:
    ci = entity.get("rcsb_polymer_entity_container_identifiers") or {}
    return list(ci.get("auth_asym_ids") or [])


def _length(entity: dict) -> int:
    ep = entity.get("entity_poly") or {}
    return int(ep.get("rcsb_sample_sequence_length") or 0)


def _description(entity: dict) -> str:
    return ((entity.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "").strip()


def partner_accession_map(evidence: dict) -> dict[str, dict]:
    """Map each interaction partner's Ensembl id to its symbol and accessions.

    One bulk Open Targets call per 50 partners rather than one call per
    partner; the PDB identifies chains by UniProt accession while Open Targets
    reports partners by Ensembl id, so the two must be joined.
    """
    ids: set[str] = set()
    for t in evidence["targets"]:
        for p in ((t.get("interactions") or {}).get("partners") or []):
            if p.get("ensembl_id"):
                ids.add(p["ensembl_id"])
    id_list = sorted(ids)
    out: dict[str, dict] = {}
    for i in range(0, len(id_list), 50):
        chunk = id_list[i : i + 50]
        for d in ot.targets_details(chunk):
            accs = [
                p["id"]
                for p in (d.get("proteinIds") or [])
                if str(p.get("source", "")).startswith("uniprot")
            ]
            out[d["id"]] = {"symbol": d.get("approvedSymbol"), "accessions": accs}
        log.info("partner mapping %d/%d", min(i + 50, len(id_list)), len(id_list))
    return out


def classify_entry(
    entry: dict,
    target_accession: str,
    partner_index: dict[str, dict],
) -> dict[str, Any]:
    """Split an entry's polymer entities into target, partners and other roles."""
    target_entities, others = [], []
    for pe in entry.get("polymer_entities") or []:
        accs = _accessions(pe)
        rec = {
            "entity_id": pe.get("rcsb_id"),
            "accessions": accs,
            "chains": _chains(pe),
            "length": _length(pe),
            "description": _description(pe),
            "polymer_type": (pe.get("entity_poly") or {}).get("type"),
        }
        if target_accession in accs:
            target_entities.append(rec)
        else:
            others.append(rec)

    for rec in others:
        hit = next(
            (
                partner_index[a]
                for a in rec["accessions"]
                if a in partner_index
            ),
            None,
        )
        if hit:
            rec["role"] = "native_partner"
            rec["partner_symbol"] = hit["symbol"]
        elif _ANTIBODY_RE.search(rec["description"]):
            rec["role"] = "antibody_fragment"
        elif _CHAPERONE_RE.search(rec["description"]):
            rec["role"] = "crystallisation_aid"
        elif not rec["accessions"]:
            rec["role"] = "unmapped_chain"
        else:
            rec["role"] = "other_protein"

    info = entry.get("rcsb_entry_info") or {}
    res = info.get("resolution_combined") or []
    return {
        "entry_id": entry.get("rcsb_id"),
        "title": (entry.get("struct") or {}).get("title"),
        "methods": [m["method"] for m in (entry.get("exptl") or []) if m.get("method")],
        "resolution": min(res) if res else None,
        "target_entities": target_entities,
        "other_entities": others,
    }


def gate_entry(classified: dict, gate: StructureGate) -> tuple[bool, str]:
    """Apply the admission rules to a classified entry."""
    methods = classified["methods"]
    if not any(m in gate.accepted_methods for m in methods):
        return False, f"method_not_accepted: {methods}"

    res = classified["resolution"]
    if res is None:
        if not gate.allow_missing_resolution:
            return False, "no_reported_resolution"
    elif res > gate.max_resolution:
        return False, f"resolution_{res}_above_max_{gate.max_resolution}"

    targets = classified["target_entities"]
    if not targets:
        return False, "target_chain_absent"
    longest = max(t["length"] for t in targets)
    if longest < gate.min_target_chain_residues:
        return False, (
            f"target_chain_too_short: {longest}aa < "
            f"{gate.min_target_chain_residues}aa (peptide fragment, not an epitope)"
        )

    partners = [o for o in classified["other_entities"] if o["role"] == "native_partner"]
    if not partners:
        roles = sorted({o["role"] for o in classified["other_entities"]}) or ["none"]
        return False, f"no_native_partner (present: {', '.join(roles)})"

    return True, "admitted"


def find_complexes(
    ranked: pd.DataFrame,
    evidence: dict,
    *,
    gate: StructureGate | None = None,
    rows_per_target: int = 200,
) -> dict[str, Any]:
    """Run the gate over every target that passed stage 2's accessibility filter."""
    gate = gate or StructureGate()
    partner_index = partner_accession_map(evidence)
    by_symbol = {t["symbol"]: t for t in evidence["targets"]}

    results: list[dict[str, Any]] = []
    for _, row in ranked[ranked["passes_accessibility"]].iterrows():
        symbol, accession = row["symbol"], row["uniprot"]
        ev = by_symbol.get(symbol) or {}
        spans = (ev.get("accessibility") or {}).get("extracellular_spans") or []
        partners_for_target = {
            a: partner_index[p["ensembl_id"]]
            for p in ((ev.get("interactions") or {}).get("partners") or [])
            if p.get("ensembl_id") in partner_index
            for a in partner_index[p["ensembl_id"]]["accessions"]
        }

        if not isinstance(accession, str) or not accession:
            results.append(
                {
                    "symbol": symbol,
                    "uniprot": None,
                    "admitted": [],
                    "n_candidates": 0,
                    "drop_reason": "no_uniprot_accession",
                }
            )
            continue

        found = rcsb.search_complexes(
            accession, max_resolution=gate.max_resolution, rows=rows_per_target
        )
        entries = rcsb.entry_details(found["ids"]) if found["ids"] else []

        admitted, rejected = [], []
        for entry in entries:
            cl = classify_entry(entry, accession, partners_for_target)
            ok, reason = gate_entry(cl, gate)
            if ok:
                target_entity = max(cl["target_entities"], key=lambda t: t["length"])
                partner_entities = [
                    o for o in cl["other_entities"] if o["role"] == "native_partner"
                ]
                admitted.append(
                    {
                        "entry_id": cl["entry_id"],
                        "title": cl["title"],
                        "method": cl["methods"][0] if cl["methods"] else None,
                        "resolution": cl["resolution"],
                        "target_entity": target_entity,
                        "partners": partner_entities,
                        "partner_symbols": sorted(
                            {p.get("partner_symbol") for p in partner_entities if p.get("partner_symbol")}
                        ),
                        "other_roles": sorted(
                            {o["role"] for o in cl["other_entities"] if o["role"] != "native_partner"}
                        ),
                    }
                )
            else:
                rejected.append({"entry_id": cl["entry_id"], "reason": reason})

        # Best-resolved complexes first; a reviewer reads the top of this list.
        admitted.sort(key=lambda a: (a["resolution"] if a["resolution"] is not None else 99.0))
        reason_counts: dict[str, int] = {}
        for r in rejected:
            key = r["reason"].split(":")[0].split(" (")[0]
            reason_counts[key] = reason_counts.get(key, 0) + 1

        results.append(
            {
                "symbol": symbol,
                "uniprot": accession,
                "composite_score": float(row["composite_score"]),
                "extracellular_spans": spans,
                "n_searched": found.get("total", 0),
                "n_candidates": len(entries),
                "n_admitted": len(admitted),
                "admitted": admitted,
                "rejection_reason_counts": reason_counts,
                "rejected_examples": rejected[:5],
                "drop_reason": None if admitted else "no_qualifying_complex",
            }
        )
        log.info(
            "%s: %d searched -> %d admitted", symbol, found.get("total", 0), len(admitted)
        )

    kept = [r for r in results if r["n_admitted"]]
    return {
        "schema_version": SCHEMA_VERSION,
        "package_version": __version__,
        "disease_id": evidence.get("disease_id"),
        "disease_name": evidence.get("disease_name"),
        "gate": asdict(gate),
        "n_targets_in": len(results),
        "n_targets_admitted": len(kept),
        "note": (
            "A native-partner complex does not guarantee a reachable interface; "
            "stage 4 enforces the extracellular-span check on real residue "
            "numbers once interfaces are computed."
        ),
        "targets": results,
    }


def write(record: dict[str, Any], run_dir: Path) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)
    path = run_dir / STAGE_FILES[3]
    path.write_text(json.dumps(record, indent=2))
    return path


def read(run_dir: Path) -> dict[str, Any]:
    path = run_dir / STAGE_FILES[3]
    if not path.exists():
        raise FileNotFoundError(f"{path} not found - run stage 3 first")
    return json.loads(path.read_text())
