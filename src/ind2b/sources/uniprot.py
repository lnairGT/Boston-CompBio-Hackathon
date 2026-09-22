"""UniProt client - topology-based accessibility for a protein binder.

Why this module exists
----------------------
Whether a de novo mini-binder can reach a target is a question about
*topology*, not about location keywords. "Cell membrane" is true of both a
receptor ectodomain and a lipid-anchored intracellular GTPase: KRAS is
annotated ``Cell membrane | Lipid-anchor | Cytoplasmic side`` and is not
reachable from outside the cell, while ALK is annotated ``Cell membrane |
Single-pass type I membrane protein`` with ``Topological domain:
Extracellular 19-1038`` and is. A keyword match on the location string calls
both accessible, which would put undruggable-by-this-modality targets at the
top of the ranking.

So accessibility is decided from UniProt sequence features:

* an **extracellular topological domain** of at least ``min_span`` residues
  -> reachable, and the span bounds the epitope search;
* otherwise a **Secreted** subcellular location -> reachable, with the mature
  chain (after signal-peptide cleavage) as the span;
* otherwise **not reachable**, with the reason recorded.

The returned spans are used again in stages 3-4 to reject candidate epitopes
that fall on an intracellular portion of a structure.
"""

from __future__ import annotations

import logging
from typing import Any

import requests

from ..config import MIN_ECTODOMAIN_SPAN
from ..http import request_json

log = logging.getLogger("ind2b.uniprot")

ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"
ENTRY_FIELDS = (
    "ft_topo_dom,ft_transmem,ft_signal,ft_lipid,ft_carbohyd,"
    "cc_subcellular_location,sequence,protein_name,gene_names"
)


def entry(accession: str) -> dict:
    """Fetch a UniProtKB entry (cached)."""
    return request_json(
        ENTRY_URL.format(acc=accession), params={"fields": ENTRY_FIELDS}
    )


def _features(data: dict, type_name: str) -> list[dict]:
    out = []
    for f in data.get("features", []) or []:
        if f.get("type") != type_name:
            continue
        loc = f.get("location") or {}
        start = (loc.get("start") or {}).get("value")
        end = (loc.get("end") or {}).get("value")
        if start is None or end is None:
            continue
        out.append(
            {"start": start, "end": end, "description": f.get("description") or ""}
        )
    return out


def _subcellular(data: dict) -> list[str]:
    out: list[str] = []
    for c in data.get("comments", []) or []:
        if c.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        for sl in c.get("subcellularLocations", []) or []:
            parts = [
                (sl.get(k) or {}).get("value")
                for k in ("location", "topology", "orientation")
            ]
            joined = " | ".join(p for p in parts if p)
            if joined:
                out.append(joined)
    return out


def accessibility(
    accession: str, *, min_span: int = MIN_ECTODOMAIN_SPAN
) -> dict[str, Any]:
    """Decide whether a binder can reach this protein, and where.

    Returns a record with ``accessible``, a ``mode``
    (``membrane_ectodomain`` / ``secreted``), the reachable
    ``extracellular_spans`` as inclusive 1-based residue ranges, and - when
    not accessible - a ``reason``. Never raises for a missing or obsolete
    accession; that becomes ``accessible: false`` with a reason.
    """
    base: dict[str, Any] = {
        "accession": accession,
        "accessible": False,
        "mode": None,
        "extracellular_spans": [],
        "signal_peptide": None,
        "transmembrane": [],
        "subcellular_locations": [],
        "sequence_length": None,
        "reason": None,
    }

    try:
        data = entry(accession)
    except (requests.HTTPError, RuntimeError) as exc:
        base["reason"] = f"uniprot_lookup_failed: {type(exc).__name__}"
        return base

    seq_len = (data.get("sequence") or {}).get("length")
    topo = _features(data, "Topological domain")
    transmem = _features(data, "Transmembrane")
    signal = _features(data, "Signal")
    lipid = _features(data, "Lipidation")
    subcell = _subcellular(data)

    base.update(
        {
            "sequence_length": seq_len,
            "transmembrane": [(t["start"], t["end"]) for t in transmem],
            "signal_peptide": (signal[0]["start"], signal[0]["end"]) if signal else None,
            "subcellular_locations": subcell,
        }
    )

    ecto = [
        t for t in topo if "extracellular" in (t["description"] or "").lower()
    ]
    ecto_ok = [t for t in ecto if (t["end"] - t["start"] + 1) >= min_span]
    if ecto_ok:
        base.update(
            {
                "accessible": True,
                "mode": "membrane_ectodomain",
                "extracellular_spans": [(t["start"], t["end"]) for t in ecto_ok],
            }
        )
        return base

    secreted = any("secreted" in s.lower() for s in subcell)
    if secreted:
        start = (signal[0]["end"] + 1) if signal else 1
        span = [(start, seq_len)] if seq_len else []
        base.update(
            {"accessible": True, "mode": "secreted", "extracellular_spans": span}
        )
        return base

    # Not reachable - record which negative signal decided it.
    if ecto:
        base["reason"] = (
            f"extracellular_topology_shorter_than_{min_span}aa"
        )
    elif any("cytoplasmic side" in s.lower() for s in subcell) or lipid:
        base["reason"] = "membrane_associated_on_cytoplasmic_side"
    elif transmem:
        base["reason"] = "transmembrane_but_no_annotated_extracellular_domain"
    else:
        base["reason"] = "no_extracellular_topology_or_secretion"
    return base


def sequence(accession: str) -> str | None:
    """Canonical sequence, used by later stages for numbering checks."""
    try:
        return (entry(accession).get("sequence") or {}).get("value")
    except (requests.HTTPError, RuntimeError):
        return None
