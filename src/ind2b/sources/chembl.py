"""ChEMBL client - mechanism and modality precedent for a target.

Used for the *pharmacological* half of precedent: what mechanisms have been
tried against this protein, and with what action type. Modality precedent
(antibody vs small molecule) comes from the Open Targets known-drug layer,
which already carries ``drugType`` per drug.
"""

from __future__ import annotations

import logging
from typing import Any

import requests

from ..config import CHEMBL_BASE
from ..http import request_json

log = logging.getLogger("ind2b.chembl")


def target_ids_for_accession(accession: str, limit: int = 5) -> list[dict]:
    """Map a UniProt accession to ChEMBL target ids (single-protein first)."""
    try:
        data = request_json(
            f"{CHEMBL_BASE}/target.json",
            params={
                "target_components__accession": accession,
                "limit": limit,
            },
        )
    except (requests.HTTPError, RuntimeError) as exc:
        log.warning("chembl target lookup failed for %s: %s", accession, exc)
        return []
    out = []
    for t in data.get("targets", []) or []:
        out.append(
            {
                "target_chembl_id": t.get("target_chembl_id"),
                "target_type": t.get("target_type"),
                "pref_name": t.get("pref_name"),
                "organism": t.get("organism"),
            }
        )
    # Single-protein targets are the informative ones; families dilute.
    out.sort(key=lambda d: 0 if d.get("target_type") == "SINGLE PROTEIN" else 1)
    return out


def mechanisms(target_chembl_id: str, limit: int = 50) -> list[dict]:
    """Known mechanisms of action recorded against a ChEMBL target."""
    try:
        data = request_json(
            f"{CHEMBL_BASE}/mechanism.json",
            params={"target_chembl_id": target_chembl_id, "limit": limit},
        )
    except (requests.HTTPError, RuntimeError) as exc:
        log.warning("chembl mechanism lookup failed for %s: %s", target_chembl_id, exc)
        return []
    return [
        {
            "mechanism_of_action": m.get("mechanism_of_action"),
            "action_type": m.get("action_type"),
            "molecule_chembl_id": m.get("molecule_chembl_id"),
            "max_phase": m.get("max_phase"),
        }
        for m in (data.get("mechanisms", []) or [])
    ]


def precedent(accession: str) -> dict[str, Any]:
    """Summarise ChEMBL precedent for one UniProt accession."""
    targets = target_ids_for_accession(accession)
    if not targets:
        return {
            "chembl_target_id": None,
            "n_mechanisms": 0,
            "action_types": {},
            "mechanisms": [],
            "note": "no ChEMBL target for accession",
        }
    chosen = targets[0]
    mech = mechanisms(chosen["target_chembl_id"])
    action_counts: dict[str, int] = {}
    for m in mech:
        at = m.get("action_type") or "UNKNOWN"
        action_counts[at] = action_counts.get(at, 0) + 1
    distinct = sorted({m["mechanism_of_action"] for m in mech if m.get("mechanism_of_action")})
    return {
        "chembl_target_id": chosen["target_chembl_id"],
        "chembl_target_type": chosen.get("target_type"),
        "n_mechanisms": len(mech),
        "n_distinct_mechanisms": len(distinct),
        "action_types": action_counts,
        "mechanisms": distinct[:8],
    }
