"""Resolving *which residues* the cross-reactivity analysis should compare.

The real answer comes from the design lane: ``contracts.BindingRegion.residues``
(or ``hotspot_residues``), in UniProt numbering. ``epitope_from_binding_region``
is the production path and is a pure translation.

``epitope_proxy_from_accessibility`` exists for one purpose: to exercise the
analysis end-to-end on real sequences *before* a real binding region exists. It
picks a contiguous window inside the topologically accessible region of the
target and labels the result ``is_proxy=True`` with ``basis`` saying exactly how
it was chosen. A proxy epitope is a placeholder for plumbing, not a claim about
where a binder would engage, and every downstream record carries the flag.
"""

from __future__ import annotations

from typing import Any

PROXY_WIDTH_DEFAULT = 25


def epitope_from_binding_region(binding_region: dict[str, Any] | None) -> dict[str, Any]:
    """Translate a ``contracts.BindingRegion``-shaped dict into epitope positions.

    Hotspot residues are preferred when present, because they are the positions the
    design was actually driven to engage; otherwise the full residue list is used.
    """
    if not binding_region:
        return {"residues": [], "is_proxy": False, "status": "unspecified",
                "basis": "no binding region supplied by the design lane"}
    hot = sorted(set(binding_region.get("hotspot_residues") or []))
    allr = sorted(set(binding_region.get("residues") or []))
    if hot:
        return {"residues": hot, "is_proxy": False, "status": "hotspots",
                "numbering_scheme": binding_region.get("numbering_scheme", "author"),
                "basis": "BindingRegion.hotspot_residues from the design lane",
                "chain": binding_region.get("chain")}
    if allr:
        return {"residues": allr, "is_proxy": False, "status": "region",
                "numbering_scheme": binding_region.get("numbering_scheme", "author"),
                "basis": "BindingRegion.residues from the design lane (no hotspots declared)",
                "chain": binding_region.get("chain")}
    return {"residues": [], "is_proxy": False, "status": "unspecified",
            "basis": "BindingRegion present but empty"}


def epitope_proxy_from_accessibility(
    protein: dict[str, Any],
    accessibility: dict[str, Any],
    width: int = PROXY_WIDTH_DEFAULT,
) -> dict[str, Any]:
    """A placeholder engaged region for pre-design plumbing tests.

    Takes the largest topologically engageable region reported by
    ``adapters.structures.assess_accessibility`` and returns a `width`-residue
    window at its midpoint. This is NOT a predicted epitope.
    """
    regions = [r for r in (accessibility.get("engageable_regions") or [])
               if r.get("start") and r.get("end") and r["end"] >= r["start"]]
    if not regions:
        return {
            "residues": [], "is_proxy": True, "status": "unavailable",
            "basis": f"no engageable region reported (accessibility verdict "
                     f"'{accessibility.get('verdict')}'); no proxy epitope can be formed",
        }
    biggest = max(regions, key=lambda r: r["end"] - r["start"])
    span = biggest["end"] - biggest["start"] + 1
    w = min(width, span)
    mid = (biggest["start"] + biggest["end"]) // 2
    start = max(biggest["start"], mid - w // 2)
    end = min(biggest["end"], start + w - 1)
    return {
        "residues": list(range(start, end + 1)),
        "is_proxy": True,
        "status": "proxy_window",
        "numbering_scheme": "uniprot_canonical",
        "region_kind": biggest.get("kind"),
        "window": {"start": start, "end": end, "width": end - start + 1},
        "basis": (
            f"{end - start + 1}-residue window at the midpoint of the largest topologically "
            f"engageable region ({biggest.get('kind')} {biggest['start']}-{biggest['end']}) from "
            f"UniProt topology. PLACEHOLDER for plumbing only: no model predicted this epitope and "
            f"no binder was designed against it."
        ),
        "accessibility_verdict": accessibility.get("verdict"),
    }
