"""Reference binders - existing molecules that already bind these targets.

**This is not a pipeline stage and it is not where designed binders come
from.** The pipeline's own binders come from a BindCraft2 campaign and are
read by ``stage6_designs``. This module exists for a narrower purpose: to
answer, for an illustration or a sanity check, *has anyone already made a
binder against this target, and does it engage the surface we selected?*

Run it explicitly (``ind2b reference-binders``) when you want that
comparison. It is deliberately not wired into ``ind2b all``, because mining
the PDB for other people's binders is not a step in designing your own, and a
report that showed them by default would invite them to be mistaken for
pipeline output.

Two kinds of existing binder are collected, both from records the run already
has or can fetch from the PDB:

**designed** - a de novo binder, miniprotein or designed scaffold (DARPin,
designed scFv) in complex with the target. These carry no UniProt accession
because they do not exist in nature, which is exactly why stage 3 set them
aside as non-partners. They are the closest published analogue to what a
BindCraft2 campaign would produce.

**antibody** - a Fab, scFv or nanobody complex. Clinically validated for
several of these targets, and evidence that the epitope is engageable by a
biologic - though an antibody's footprint is much larger than a mini-binder's.

For each, the target-side footprint is computed with the same interface
definition stage 4 uses, mapped into UniProt numbering through the SIFTS
alignment, and compared with the epitope this pipeline selected. Overlap is
**supporting evidence that the chosen surface is bindable, not proof that a
new binder will work**: a designed binder existing against the same patch says
the patch is tractable; it says nothing about affinity or specificity of a
binder not yet made.

A 2-D projection of each example complex is computed here and persisted, so a
report can draw the structure without re-reading coordinate files.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from . import __version__
from .config import REFERENCE_BINDERS_FILE, InterfaceParams
from .sources import rcsb
from .stage3_complexes import classify_entry, partner_accession_map
from .stage4_interface import auth_to_uniprot, interface_residues

log = logging.getLogger("ind2b.reference_binders")

SCHEMA_VERSION = 2

RECORD_KIND = "reference_binders"

# Descriptions of chains that are engineered rather than natural. Matched only
# on chains with no UniProt accession, so a natural protein whose name happens
# to contain "designed" is not miscalled.
DESIGNED_RE = re.compile(
    r"\b(de novo|de-novo|designed|design(?:ed)? (?:protein|binder|mini)|"
    r"mini[- ]?binder|miniprotein|mini[- ]?protein|computationally designed|"
    r"synthetic (?:binder|protein)|hallucinat|ankyrin repeat|darpin|affibody|"
    r"monobody)",
    re.IGNORECASE,
)

# A designed mini-binder is small; a DARPin or scFv is not. Recorded so a
# reader can tell a mini-binder analogue from a large designed scaffold.
MINIBINDER_MAX_RESIDUES = 100


def _binder_kind(entity: dict) -> str | None:
    """Classify a non-target chain as a designed binder, an antibody, or neither."""
    desc = entity.get("description") or ""
    if entity.get("accessions"):
        return None  # a natural protein: partner or bystander, not a made binder
    if DESIGNED_RE.search(desc):
        return "designed"
    if entity.get("role") == "antibody_fragment":
        return "antibody"
    return None


def find_binder_complexes(
    symbol: str,
    accession: str,
    partner_index: dict,
    *,
    max_resolution: float = 4.0,
    rows: int = 200,
    min_target_residues: int = 40,
) -> list[dict[str, Any]]:
    """PDB entries pairing ``accession`` with a designed binder or an antibody."""
    found = rcsb.search_complexes(accession, max_resolution=max_resolution, rows=rows)
    if not found["ids"]:
        return []
    out: list[dict[str, Any]] = []
    for entry in rcsb.entry_details(found["ids"]):
        cl = classify_entry(entry, accession, partner_index)
        targets = [t for t in cl["target_entities"] if t["length"] >= min_target_residues]
        if not targets:
            continue
        target_entity = max(targets, key=lambda t: t["length"])
        groups: dict[str, list[dict]] = {}
        for other in cl["other_entities"]:
            kind = _binder_kind(other)
            if kind:
                groups.setdefault(kind, []).append(other)
        for kind, chains in groups.items():
            binder_len = max(c["length"] for c in chains)
            out.append({
                "symbol": symbol,
                "entry_id": cl["entry_id"],
                "kind": kind,
                "is_minibinder_scale": kind == "designed"
                and binder_len <= MINIBINDER_MAX_RESIDUES,
                "resolution": cl["resolution"],
                "method": cl["methods"][0] if cl["methods"] else None,
                "title": cl["title"],
                "target_chain": (target_entity["chains"] or [None])[0],
                "target_length": target_entity["length"],
                "binder_chains": [c for ch in chains for c in (ch["chains"] or [])],
                "binder_length": binder_len,
                "binder_descriptions": sorted({c["description"] for c in chains}),
                "n_binder_entities": len(chains),
            })
    # Best-resolved first: that is the structure a reader should look at.
    out.sort(key=lambda r: (r["resolution"] if r["resolution"] is not None else 99.0))
    return out


def footprint(
    complex_rec: dict,
    accession: str,
    spans: Sequence[tuple[int, int]],
    struct_dir: Path,
    params: InterfaceParams,
) -> dict[str, Any] | None:
    """Target-side footprint of the binder, in author and UniProt numbering."""
    entry_id = complex_rec["entry_id"]
    target_chain = complex_rec["target_chain"]
    binder_chains = complex_rec["binder_chains"]
    if not target_chain or not binder_chains:
        return None
    try:
        path = rcsb.download_structure(entry_id, struct_dir)
        residues = interface_residues(path, target_chain, binder_chains[:1], params)
    except Exception as exc:  # noqa: BLE001
        log.warning("%s footprint failed: %s", entry_id, exc)
        return None
    if not residues:
        return None

    aligns = rcsb.entry_alignments([entry_id])
    numbering = auth_to_uniprot(aligns.get(entry_id, {}), accession, target_chain)
    for r in residues:
        pos = numbering.get(r["auth_seq_id"])
        r["uniprot_pos"] = pos
        r["in_ectodomain"] = (
            any(lo <= pos <= hi for lo, hi in spans) if pos is not None else None
        )
    mapped = [r for r in residues if r["uniprot_pos"] is not None]
    covered = sorted({p for p in numbering.values() if p is not None})
    return {
        "structure_file": f"structures/{path.name}",
        # Which part of the protein this construct actually contains. A zero
        # overlap means something different when the entry does not include the
        # selected epitope at all: that is silence about the surface, not
        # evidence against it.
        "construct_uniprot_range": [covered[0], covered[-1]] if covered else None,
        "n_construct_residues_mapped": len(covered),
        "_covered_uniprot": covered,
        "binder_chain_used": binder_chains[0],
        "n_footprint_residues": len(residues),
        "n_mapped": len(mapped),
        "buried_area_total": round(sum(r["delta_sasa"] for r in residues), 1),
        "footprint_uniprot": sorted(r["uniprot_pos"] for r in mapped),
        "footprint_auth": [r["auth_seq_id"] for r in residues],
        "top_residues": [
            {k: r.get(k) for k in ("residue", "auth_seq_id", "uniprot_pos",
                                    "delta_sasa", "is_contact", "in_ectodomain")}
            for r in residues[:12]
        ],
    }


def compare_to_epitope(
    footprint_positions: Sequence[int],
    epitope_positions: Sequence[int],
    hotspot_positions: Sequence[int] = (),
) -> dict[str, Any]:
    """Overlap between an existing binder's footprint and the selected epitope."""
    fp, ep = set(footprint_positions), set(epitope_positions)
    hs = set(hotspot_positions)
    if not fp or not ep:
        return {"comparable": False,
                "reason": "footprint or epitope could not be mapped to UniProt numbering"}
    inter = fp & ep
    return {
        "comparable": True,
        "n_shared": len(inter),
        "jaccard": round(len(inter) / len(fp | ep), 3),
        "fraction_of_epitope_covered": round(len(inter) / len(ep), 3),
        "fraction_of_hotspots_covered": (
            round(len(hs & fp) / len(hs), 3) if hs else None
        ),
        "hotspots_covered": sorted(hs & fp),
        "hotspots_missed": sorted(hs - fp),
        "shared_positions": sorted(inter),
    }


def project_complex(
    structure_path: Path,
    target_chain: str,
    binder_chains: Sequence[str],
    *,
    epitope_auth: Sequence[str] = (),
    hotspot_auth: Sequence[str] = (),
) -> dict[str, Any] | None:
    """2-D projection of CA traces for drawing, computed once and persisted.

    The projection is the first two principal components of the complex's CA
    coordinates, which puts the largest spread in the plane of the page. This
    is a schematic of the backbone, not a rendering of the molecular surface,
    and is labelled as such wherever it is drawn.
    """
    from Bio.PDB import MMCIFParser

    try:
        structure = MMCIFParser(QUIET=True).get_structure("s", str(structure_path))
    except Exception as exc:  # noqa: BLE001
        log.warning("projection parse failed for %s: %s", structure_path, exc)
        return None
    model = next(iter(structure))
    keep = [target_chain, *binder_chains]
    chains: dict[str, list[tuple[str, np.ndarray]]] = {}
    for chain in model:
        if chain.id not in keep:
            continue
        pts = []
        for res in chain:
            if res.id[0] != " " or "CA" not in res:
                continue
            key = str(res.id[1]) + (res.id[2].strip() or "")
            pts.append((key, res["CA"].coord.astype(float)))
        if pts:
            chains[chain.id] = pts
    if target_chain not in chains:
        return None

    allc = np.array([c for pts in chains.values() for _, c in pts])
    centre = allc.mean(axis=0)
    # SVD on the centred cloud; rows of Vt are the principal directions.
    _, _, vt = np.linalg.svd(allc - centre, full_matrices=False)

    epi, hot = set(epitope_auth), set(hotspot_auth)
    out_chains = []
    for cid, pts in chains.items():
        # All three components, in the principal-axis frame: the first two are
        # the flat projection used for print, the third is depth. Storing the
        # full frame is what lets a viewer rotate the complex without needing
        # the original coordinate file.
        full = (np.array([c for _, c in pts]) - centre) @ vt.T
        role = "target" if cid == target_chain else "binder"
        out_chains.append({
            "chain": cid,
            "role": role,
            "points": [[round(float(x), 1), round(float(y), 1)] for x, y, _ in full],
            "xyz": [[round(float(x), 1), round(float(y), 1), round(float(z), 1)]
                    for x, y, z in full],
            "residues": [k for k, _ in pts],
            "flags": [
                ("hotspot" if k in hot else "epitope" if k in epi else "")
                for k, _ in pts
            ] if role == "target" else [""] * len(pts),
        })
    xs = [p[0] for c in out_chains for p in c["points"]]
    ys = [p[1] for c in out_chains for p in c["points"]]
    radius = float(np.linalg.norm(
        np.array([p for c in out_chains for p in c["xyz"]]), axis=1).max())
    return {
        "projection": "principal-axis frame of the complex CA coordinates",
        "is_schematic": True,
        "bounds": [round(min(xs), 1), round(min(ys), 1), round(max(xs), 1), round(max(ys), 1)],
        "radius": round(radius, 1),
        "rotatable": True,
        "chains": out_chains,
    }


def run(
    complexes: dict,
    epitopes: pd.DataFrame,
    evidence: dict,
    *,
    run_dir: Path,
    params: InterfaceParams | None = None,
    max_per_kind: int = 1,
) -> dict[str, Any]:
    """Collect existing binders for every target that cleared stage 3."""
    params = params or InterfaceParams()
    struct_dir = run_dir / "structures"
    partner_index = partner_accession_map(evidence)
    ev_by_symbol = {t["symbol"]: t for t in evidence.get("targets") or []}
    accepted = epitopes[epitopes["status"] == "accepted"]
    best_epitope = (
        accepted.sort_values("combined_score", ascending=False)
        .groupby("symbol", as_index=False)
        .first()
        .set_index("symbol")
    )

    results: list[dict[str, Any]] = []
    for t in complexes.get("targets") or []:
        symbol, accession = t["symbol"], t.get("uniprot")
        if not t.get("n_admitted") or not accession:
            continue
        spans = [tuple(s) for s in (t.get("extracellular_spans") or [])]
        ev = ev_by_symbol.get(symbol) or {}
        drugs = (ev.get("known_drugs") or {})
        clinical = [
            d for d in (drugs.get("top") or [])
            if str(d.get("drug_type", "")).lower() in
            ("antibody", "protein", "enzyme", "oligosaccharide")
        ]

        found = find_binder_complexes(symbol, accession, partner_index)
        selected: list[dict[str, Any]] = []
        for kind in ("designed", "antibody"):
            for rec in [r for r in found if r["kind"] == kind][:max_per_kind]:
                fp = footprint(rec, accession, spans, struct_dir, params)
                rec["footprint"] = fp
                if fp and symbol in best_epitope.index:
                    row = best_epitope.loc[symbol]
                    epi_pos = [int(x) for x in str(row.get("interface_uniprot", "")).split(";")
                               if x.isdigit()]
                    hot_pos = [int(x) for x in str(row.get("hotspots_uniprot", "")).split(";")
                               if x.isdigit()]
                    rec["vs_selected_epitope"] = compare_to_epitope(
                        fp["footprint_uniprot"], epi_pos, hot_pos)
                    covered = set(fp.pop("_covered_uniprot", []))
                    present = sorted(set(epi_pos) & covered)
                    rec["vs_selected_epitope"].update({
                        "n_epitope_positions_in_construct": len(present),
                        "epitope_present_in_construct": bool(present),
                        "interpretation": (
                            "the construct in this entry does not contain the "
                            "selected epitope, so its zero overlap says nothing "
                            "about that surface"
                            if not present else
                            "the selected epitope is present in this construct, so "
                            "the overlap is informative about whether this binder "
                            "engages it"
                        ),
                    })
                if fp is not None:
                    fp.pop("_covered_uniprot", None)
                    rec["selected_epitope"] = {
                        "entry_id": row.get("entry_id"),
                        "partner_symbol": row.get("partner_symbol"),
                        "hotspots_uniprot": row.get("hotspots_uniprot"),
                    }
                selected.append(rec)

        results.append({
            "symbol": symbol,
            "uniprot": accession,
            "composite_score": t.get("composite_score"),
            "clinical_binders": clinical,
            "n_clinical_binders": len(clinical),
            "has_approved_biologic": bool(drugs.get("has_approved")
                                          and drugs.get("has_biologic_precedent")),
            "n_designed_complexes": sum(1 for r in found if r["kind"] == "designed"),
            "n_antibody_complexes": sum(1 for r in found if r["kind"] == "antibody"),
            "binders": selected,
        })
        log.info("%s: %d designed, %d antibody complexes", symbol,
                 results[-1]["n_designed_complexes"], results[-1]["n_antibody_complexes"])

    return {
        "schema_version": SCHEMA_VERSION,
        "record_kind": RECORD_KIND,
        "package_version": __version__,
        "disease_id": complexes.get("disease_id"),
        "disease_name": complexes.get("disease_name"),
        "minibinder_max_residues": MINIBINDER_MAX_RESIDUES,
        "max_per_kind": max_per_kind,
        "n_targets": len(results),
        "n_with_designed": sum(1 for r in results if r["n_designed_complexes"]),
        "not_our_designs": (
            "Every binder in this record is an existing molecule from the PDB or "
            "from clinical data. None was designed by this pipeline, which has "
            "produced no binder. Overlap with a selected epitope is evidence that "
            "the surface is bindable, not that a new binder against it would work."
        ),
        "targets": results,
    }


def add_projections(
    record: dict,
    run_dir: Path,
    epitopes: pd.DataFrame,
    *,
    max_examples: int = 6,
) -> dict:
    """Attach 2-D CA-trace projections for the best example complexes."""
    accepted = epitopes[epitopes["status"] == "accepted"]
    best = (accepted.sort_values("combined_score", ascending=False)
            .groupby("symbol", as_index=False).first().set_index("symbol"))
    # Designed binders first, best-resolved: those are the informative examples.
    candidates = [
        (t, b) for t in record["targets"] for b in t["binders"]
        if b.get("footprint") and b["kind"] == "designed"
    ] + [
        (t, b) for t in record["targets"] for b in t["binders"]
        if b.get("footprint") and b["kind"] == "antibody"
    ]
    candidates.sort(key=lambda tb: (tb[1]["kind"] != "designed",
                                    tb[1]["resolution"] or 99.0))
    n = 0
    for t, b in candidates:
        if n >= max_examples:
            break
        path = run_dir / b["footprint"]["structure_file"]
        if not path.exists():
            continue
        epi_auth, hot_auth, source = [], [], "binder footprint in this entry"
        if t["symbol"] in best.index:
            row = best.loc[t["symbol"]]
            if row.get("entry_id") == b["entry_id"]:
                epi_auth = [x for x in str(row.get("interface_auth", "")).split(";") if x]
                hot_auth = [x for x in str(row.get("hotspots_auth", "")).split(";") if x]
                source = "epitope selected by this pipeline (same entry)"
            else:
                # The selected epitope came from a different PDB entry, whose
                # author numbering is unrelated to this one. Translate through
                # UniProt numbering - the only frame the two entries share -
                # rather than dropping the comparison or, worse, reusing author
                # ids across entries as if they matched.
                aligns = rcsb.entry_alignments([b["entry_id"]])
                fwd = auth_to_uniprot(aligns.get(b["entry_id"], {}),
                                      t["uniprot"], b["target_chain"])
                back: dict[int, str] = {}
                for auth_id, up in fwd.items():
                    back.setdefault(up, auth_id)
                epi_up = [int(x) for x in str(row.get("interface_uniprot", "")).split(";")
                          if x.isdigit()]
                hot_up = [int(x) for x in str(row.get("hotspots_uniprot", "")).split(";")
                          if x.isdigit()]
                epi_auth = [back[u] for u in epi_up if u in back]
                hot_auth = [back[u] for u in hot_up if u in back]
                if epi_auth:
                    source = (f'epitope selected from {row.get("entry_id")}, mapped '
                              f'into this entry through UniProt numbering')
        if not epi_auth:
            epi_auth = b["footprint"]["footprint_auth"]
            source = "binder footprint in this entry"
        # Only the chain that actually makes the footprint. An asymmetric unit
        # often holds several copies of the binder (8ZNL has four), and drawing
        # the non-contacting copies over one target would imply contacts that
        # are not there.
        proj = project_complex(path, b["target_chain"],
                               [b["footprint"]["binder_chain_used"]],
                               epitope_auth=epi_auth, hotspot_auth=hot_auth)
        if proj:
            proj["epitope_source"] = source
            proj["n_hotspots_drawn"] = len(hot_auth)
            b["projection"] = proj
            n += 1
    record["n_projections"] = n
    return record


def write(record: dict[str, Any], run_dir: Path) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)
    path = run_dir / REFERENCE_BINDERS_FILE
    path.write_text(json.dumps(record, indent=2))
    return path


def read(run_dir: Path) -> dict[str, Any]:
    path = run_dir / REFERENCE_BINDERS_FILE
    if not path.exists():
        raise FileNotFoundError(f"{path} not found - run `ind2b reference-binders` first")
    return json.loads(path.read_text())
