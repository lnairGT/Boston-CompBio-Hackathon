"""ClinicalTrials.gov v2 client - clinical precedent, including failures.

A terminated or withdrawn trial is signal, not noise: it says something was
tried in patients. This client therefore returns the status breakdown rather
than a single count, and records NCT identifiers so any claim about a
programme can be checked against the registry.
"""

from __future__ import annotations

import logging
from typing import Any

import requests

from ..config import CLINICALTRIALS_V2
from ..http import request_json

log = logging.getLogger("ind2b.trials")

# Statuses that indicate a programme stopped early.
STOPPED = {"TERMINATED", "WITHDRAWN", "SUSPENDED"}

_FIELDS = ",".join(
    [
        "NCTId",
        "BriefTitle",
        "OverallStatus",
        "Phase",
        "WhyStopped",
        "InterventionName",
    ]
)


def search(
    condition: str,
    term: str,
    *,
    page_size: int = 20,
) -> dict[str, Any]:
    """Trials matching a condition and a free-text term (usually a gene symbol)."""
    try:
        data = request_json(
            CLINICALTRIALS_V2,
            params={
                "query.cond": condition,
                "query.term": term,
                "fields": _FIELDS,
                "countTotal": "true",
                "pageSize": page_size,
            },
        )
    except (requests.HTTPError, RuntimeError) as exc:
        log.warning("trials lookup failed for %s / %s: %s", condition, term, exc)
        return {"total": 0, "studies": [], "note": f"lookup failed: {type(exc).__name__}"}

    studies = []
    for s in data.get("studies", []) or []:
        proto = s.get("protocolSection") or {}
        ident = proto.get("identificationModule") or {}
        status = proto.get("statusModule") or {}
        design = proto.get("designModule") or {}
        studies.append(
            {
                "nct_id": ident.get("nctId"),
                "title": ident.get("briefTitle"),
                "status": status.get("overallStatus"),
                "why_stopped": status.get("whyStopped"),
                "phases": design.get("phases") or [],
            }
        )
    return {"total": data.get("totalCount", len(studies)), "studies": studies}


def precedent(condition: str, symbol: str) -> dict[str, Any]:
    """Summarise clinical precedent for one target in one indication."""
    res = search(condition, symbol)
    studies = res.get("studies") or []
    status_counts: dict[str, int] = {}
    for s in studies:
        st = s.get("status") or "UNKNOWN"
        status_counts[st] = status_counts.get(st, 0) + 1
    stopped = [
        {"nct_id": s["nct_id"], "status": s["status"], "why_stopped": s.get("why_stopped")}
        for s in studies
        if (s.get("status") or "") in STOPPED
    ]
    phases = sorted({p for s in studies for p in (s.get("phases") or [])})
    return {
        "total_trials": res.get("total", 0),
        "n_sampled": len(studies),
        "status_counts": status_counts,
        "phases_seen": phases,
        "stopped_trials": stopped[:5],
        "example_nct_ids": [s["nct_id"] for s in studies[:5] if s.get("nct_id")],
        "note": res.get("note"),
    }
