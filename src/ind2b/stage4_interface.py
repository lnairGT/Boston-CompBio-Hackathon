"""Stage 4 - compute the target-side interface and select hotspots.

For each admitted complex the target-side interface is defined by two
independent criteria, both recorded so the definition can be audited:

* **contact** - a heavy atom within ``contact_cutoff`` A of any partner heavy
  atom;
* **burial** - at least ``min_delta_sasa`` A^2 of solvent-accessible surface
  lost on complex formation, computed with Shrake-Rupley.

A residue satisfying either is called interface. Burial is what ranks them,
because buried area is the quantity that correlates with binding energetics,
while contact counting is sensitive to the cutoff.

The topological check
---------------------
Being a native-partner complex is not sufficient: CBL-C binds the
*intracellular* tail of EGFR, so that interface is unreachable by an
extracellular binder. Each interface residue is therefore mapped from PDB
author numbering into UniProt numbering via the SIFTS alignment and tested
against the extracellular spans from stage 1. A complex whose interface falls
mostly outside those spans is rejected here, with the fraction recorded.

Author numbering is never assumed to equal UniProt numbering - that
assumption silently corrupts the check whenever a construct is numbered from
its own residue 1.
"""

from __future__ import annotations

import copy
import json
import logging
from pathlib import Path
from typing import Any, Iterable, Sequence

import pandas as pd
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley

from . import __version__
from .config import STAGE_FILES, InterfaceParams
from .sources import rcsb

log = logging.getLogger("ind2b.stage4")

SCHEMA_VERSION = 1

# Fraction of interface residues that must lie inside an annotated
# extracellular span for the complex to be usable for an extracellular binder.
MIN_FRACTION_IN_ECTODOMAIN = 0.60

# Minimum buried area for a contact to count as an epitope at all. Splitting
# interfaces per partner chain surfaces incidental contacts - a few tens of
# square angstrom between chains that merely touch in the lattice - and those
# are crystal packing, not binding sites.
MIN_BURIED_AREA = 150.0


def auth_to_uniprot(entry_align: dict, accession: str, chain: str) -> dict[str, int]:
    """Build an author-residue-id -> UniProt-position map for one chain."""
    for entity in entry_align.get("entities") or []:
        align = next(
            (a for a in entity["aligns"] if a["accession"] == accession), None
        )
        if not align or chain not in entity["chains"]:
            continue
        auth_ids = entity["chains"][chain]
        mapping: dict[str, int] = {}
        for region in align["regions"]:
            for k in range(region["length"]):
                entity_pos = region["entity_beg"] + k  # 1-based
                idx = entity_pos - 1
                if 0 <= idx < len(auth_ids):
                    mapping[str(auth_ids[idx])] = region["ref_beg"] + k
        return mapping
    return {}


# Radius (A, CB-CB) within which hotspots are considered one patch. A binder
# engages a contiguous surface, so hotspots must describe a single patch rather
# than the highest-burial residues scattered across the whole interface.
HOTSPOT_PATCH_RADIUS = 12.0


def select_hotspot_patch(
    candidates: list[dict],
    n_wanted: int,
    patch_radius: float = HOTSPOT_PATCH_RADIUS,
) -> list[dict]:
    """Pick a spatially compact set of hotspot residues.

    Seeded on the most-buried contact residue, then grown by proximity to the
    growing patch centroid. Taking the top N by buried area instead can spread
    hotspots across opposite edges of a large interface, which asks the
    designer for a binder that bridges both - a harder problem than the one
    intended. The reference BindCraft2 PD-L1 target definition likewise uses a
    single compact pocket (residues 54, 56, 66, 115) rather than the
    highest-burial residues overall.
    """
    usable = [c for c in candidates if c.get("coord")]
    if not usable:
        return candidates[:n_wanted]

    seed = max(usable, key=lambda c: c["delta_sasa"])
    patch = [seed]
    centroid = list(seed["coord"])

    while len(patch) < n_wanted:
        remaining = [c for c in usable if c not in patch]
        if not remaining:
            break

        def dist(c):
            return sum((a - b) ** 2 for a, b in zip(c["coord"], centroid)) ** 0.5

        nearest = min(remaining, key=dist)
        if dist(nearest) > patch_radius:
            break
        patch.append(nearest)
        n = len(patch)
        centroid = [
            sum(c["coord"][i] for c in patch) / n for i in range(3)
        ]
    return patch


def _load_model(path: Path, keep_chains: Iterable[str], ignore_hetero: bool):
    """Parse a structure and reduce it to the chains of interest."""
    structure = MMCIFParser(QUIET=True).get_structure("s", str(path))
    model = next(iter(structure))  # first model; NMR ensembles use model 0
    keep = set(keep_chains)
    for chain in list(model):
        if chain.id not in keep:
            model.detach_child(chain.id)
    if ignore_hetero:
        for chain in model:
            for res in list(chain):
                if res.id[0] != " ":  # waters, ions, ligands
                    chain.detach_child(res.id)
    return model


def _residue_sasa(model) -> dict[tuple[str, str], float]:
    ShrakeRupley().compute(model, level="R")
    return {
        (chain.id, str(res.id[1]) + (res.id[2].strip() or "")): float(res.sasa)
        for chain in model
        for res in chain
    }


def interface_residues(
    path: Path,
    target_chain: str,
    partner_chains: Sequence[str],
    params: InterfaceParams,
) -> list[dict[str, Any]]:
    """Per-residue interface record for the target side of one complex."""
    complex_model = _load_model(
        path, [target_chain, *partner_chains], params.ignore_hetero
    )
    if target_chain not in [c.id for c in complex_model]:
        return []

    # Burial: SASA alone minus SASA in complex, on identical atom selections.
    alone_model = copy.deepcopy(complex_model)
    for chain in list(alone_model):
        if chain.id != target_chain:
            alone_model.detach_child(chain.id)

    sasa_complex = _residue_sasa(copy.deepcopy(complex_model))
    sasa_alone = _residue_sasa(alone_model)

    # Contacts: target heavy atoms against partner heavy atoms.
    target_atoms, partner_atoms = [], []
    for chain in complex_model:
        for res in chain:
            for atom in res:
                if atom.element == "H":
                    continue
                (target_atoms if chain.id == target_chain else partner_atoms).append(atom)
    if not partner_atoms or not target_atoms:
        return []

    search = NeighborSearch(partner_atoms)
    min_dist: dict[str, float] = {}
    for atom in target_atoms:
        close = search.search(atom.coord, params.contact_cutoff)
        if not close:
            continue
        res = atom.get_parent()
        key = str(res.id[1]) + (res.id[2].strip() or "")
        d = min(float((atom - other)) for other in close)
        if key not in min_dist or d < min_dist[key]:
            min_dist[key] = d

    rows: list[dict[str, Any]] = []
    for chain in complex_model:
        if chain.id != target_chain:
            continue
        for res in chain:
            key = str(res.id[1]) + (res.id[2].strip() or "")
            delta = sasa_alone.get((target_chain, key), 0.0) - sasa_complex.get(
                (target_chain, key), 0.0
            )
            contact = key in min_dist
            if not contact and delta < params.min_delta_sasa:
                continue
            anchor = res["CB"] if "CB" in res else (res["CA"] if "CA" in res else None)
            rows.append(
                {
                    "auth_seq_id": key,
                    "residue": res.get_resname(),
                    "coord": [float(c) for c in anchor.coord] if anchor else None,
                    "delta_sasa": round(delta, 2),
                    "min_distance": round(min_dist.get(key, float("nan")), 2)
                    if contact
                    else None,
                    "is_contact": contact,
                    "is_buried": delta >= params.min_delta_sasa,
                }
            )
    rows.sort(key=lambda r: -r["delta_sasa"])
    return rows


def analyse_target(
    target: dict,
    *,
    struct_dir: Path,
    params: InterfaceParams,
    max_complexes: int = 2,
) -> list[dict[str, Any]]:
    """Compute interfaces for the best complexes of one admitted target."""
    spans = [tuple(s) for s in (target.get("extracellular_spans") or [])]
    accession = target["uniprot"]
    candidates = target["admitted"][:max_complexes]
    if not candidates:
        return []

    aligns = rcsb.entry_alignments([c["entry_id"] for c in candidates])
    out: list[dict[str, Any]] = []

    for cand in candidates:
        entry_id = cand["entry_id"]
        target_chains = cand["target_entity"]["chains"] or []
        partner_chain_map = {
            c: p.get("partner_symbol")
            for p in cand["partners"]
            for c in (p["chains"] or [])
        }
        partner_chains = list(partner_chain_map)
        if not target_chains or not partner_chains:
            out.append({"entry_id": entry_id, "symbol": target["symbol"],
                         "status": "rejected", "reason": "missing chain assignment"})
            continue

        try:
            path = rcsb.download_structure(entry_id, struct_dir)
        except Exception as exc:  # noqa: BLE001
            out.append({"entry_id": entry_id, "symbol": target["symbol"],
                         "status": "rejected", "reason": f"download failed: {type(exc).__name__}"})
            continue

        target_chain = target_chains[0]
        numbering = auth_to_uniprot(aligns.get(entry_id, {}), accession, target_chain)

        # One epitope per partner CHAIN, not per entry. A target contacting two
        # partner chains at once (a receptor dimer plus its ligand, say) yields
        # a composite interface of several thousand square angstrom that no
        # single mini-binder reproduces; the designable unit is the patch
        # engaged by one partner.
        for partner_chain in partner_chains:
            try:
                residues = interface_residues(
                    path, target_chain, [partner_chain], params
                )
            except Exception as exc:  # noqa: BLE001
                out.append({"entry_id": entry_id, "symbol": target["symbol"],
                             "partner_chain": partner_chain, "status": "rejected",
                             "reason": f"parse/SASA failed: {type(exc).__name__}: {exc}"[:120]})
                continue
            if not residues:
                out.append({"entry_id": entry_id, "symbol": target["symbol"],
                             "partner_chain": partner_chain, "status": "rejected",
                             "reason": "no interface residues found"})
                continue
            out.append(
                _epitope_record(
                    target=target,
                    cand=cand,
                    accession=accession,
                    target_chain=target_chain,
                    partner_chain=partner_chain,
                    partner_symbol=partner_chain_map.get(partner_chain),
                    residues=residues,
                    numbering=numbering,
                    spans=spans,
                    path=path,
                    params=params,
                )
            )
    return out


def _epitope_record(
    *,
    target: dict,
    cand: dict,
    accession: str,
    target_chain: str,
    partner_chain: str,
    partner_symbol: str | None,
    residues: list[dict],
    numbering: dict[str, int],
    spans: list[tuple[int, int]],
    path: Path,
    params: InterfaceParams,
) -> dict[str, Any]:
    """Assemble one epitope record and apply the topological check."""
    entry_id = cand["entry_id"]
    for r in residues:
        pos = numbering.get(r["auth_seq_id"])
        r["uniprot_pos"] = pos
        r["in_ectodomain"] = (
            any(lo <= pos <= hi for lo, hi in spans) if pos is not None else None
        )

    mapped = [r for r in residues if r["uniprot_pos"] is not None]
    inside = [r for r in mapped if r["in_ectodomain"]]
    frac = len(inside) / len(mapped) if mapped else 0.0

    record: dict[str, Any] = {
        "symbol": target["symbol"],
        "uniprot": accession,
        "entry_id": entry_id,
        "title": cand.get("title"),
        "method": cand.get("method"),
        "resolution": cand.get("resolution"),
        "partner_symbol": partner_symbol,
        "partner_symbols": cand.get("partner_symbols") or [],
        "partner_chain": partner_chain,
        "target_chain": target_chain,
        "n_interface_residues": len(residues),
        "n_mapped_to_uniprot": len(mapped),
        "fraction_in_ectodomain": round(frac, 3),
        "buried_area_total": round(sum(r["delta_sasa"] for r in residues), 1),
        # Recorded relative to the run directory, not the working directory:
        # an absolute or cwd-relative path here makes the output artifact
        # unusable from anywhere else, including by the next stage.
        "structure_file": f"structures/{path.name}",
        "residues": residues,
    }

    buried_total = record["buried_area_total"]
    if not mapped:
        record.update(status="rejected", reason="no residue mapped to UniProt numbering")
    elif buried_total < MIN_BURIED_AREA:
        record.update(
            status="rejected",
            reason=(
                f"interface too small: {buried_total:.0f} A^2 < {MIN_BURIED_AREA:.0f} "
                f"A^2 (incidental or lattice contact, not a binding site)"
            ),
        )
    elif frac < MIN_FRACTION_IN_ECTODOMAIN:
        record.update(
            status="rejected",
            reason=(
                f"interface outside ectodomain: only {frac:.0%} of mapped "
                f"interface residues fall in {spans} (intracellular or "
                f"non-ectodomain interface)"
            ),
        )
    else:
        contacts_inside = [r for r in inside if r["is_contact"]]
        hotspots = select_hotspot_patch(contacts_inside, params.n_hotspots)
        hotspots.sort(key=lambda h: -h["delta_sasa"])
        coords = [h["coord"] for h in hotspots if h.get("coord")]
        if len(coords) > 1:
            centre = [sum(c[i] for c in coords) / len(coords) for i in range(3)]
            spread = max(
                sum((c[i] - centre[i]) ** 2 for i in range(3)) ** 0.5 for c in coords
            )
        else:
            spread = 0.0
        record.update(
            status="accepted",
            reason=None,
            hotspot_patch_radius=round(spread, 2),
            hotspots=[
                {
                    "auth_seq_id": h["auth_seq_id"],
                    "uniprot_pos": h["uniprot_pos"],
                    "residue": h["residue"],
                    "delta_sasa": h["delta_sasa"],
                }
                for h in hotspots
            ],
        )
    return record


# A designed mini-binder typically buries on the order of 600-1000 A^2 at its
# interface. Epitopes far below that band offer too little surface to engage;
# epitopes far above it are usually composite - a target contacting several
# partner chains, or a large obligate assembly - and no single small binder
# reproduces them. So epitope size is scored as a preference band rather than
# "more is better", which is why buried area alone must not drive the ranking.
IDEAL_BURIED_LOW = 500.0
IDEAL_BURIED_HIGH = 1200.0


def designability(
    buried_area: float,
    resolution: float | None,
    n_hotspot_contacts: int,
    n_hotspots_wanted: int,
) -> dict[str, float]:
    """Score how tractable an epitope looks for de novo binder design.

    This is a heuristic prior on epitope geometry, not a prediction of binding.
    Components are returned individually so a reviewer can disagree with the
    weighting without re-deriving the inputs.
    """
    if buried_area <= 0:
        size = 0.0
    elif buried_area < IDEAL_BURIED_LOW:
        size = buried_area / IDEAL_BURIED_LOW
    elif buried_area <= IDEAL_BURIED_HIGH:
        size = 1.0
    else:
        # Composite interfaces decay rather than being excluded outright.
        size = max(0.25, IDEAL_BURIED_HIGH / buried_area)

    if resolution is None:
        res_score = 0.5
    elif resolution <= 2.5:
        res_score = 1.0
    else:
        res_score = max(0.2, min(1.0, (4.5 - resolution) / 2.0))

    hotspot_score = (
        min(1.0, n_hotspot_contacts / n_hotspots_wanted) if n_hotspots_wanted else 0.0
    )

    score = 0.5 * size + 0.2 * res_score + 0.3 * hotspot_score
    return {
        "size_score": round(size, 4),
        "resolution_score": round(res_score, 4),
        "hotspot_score": round(hotspot_score, 4),
        "designability": round(score, 4),
    }


def run(
    complexes: dict,
    *,
    run_dir: Path,
    params: InterfaceParams | None = None,
    max_complexes: int = 2,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Compute interfaces for every admitted target; returns epitopes, residues, meta."""
    params = params or InterfaceParams()
    struct_dir = run_dir / "structures"
    epitopes: list[dict[str, Any]] = []

    for target in complexes["targets"]:
        if not target.get("n_admitted"):
            continue
        for rec in analyse_target(
            target, struct_dir=struct_dir, params=params, max_complexes=max_complexes
        ):
            epitopes.append(rec)
        log.info("interfaces done: %s", target["symbol"])

    residue_rows = [
        {
            "symbol": e["symbol"],
            "entry_id": e["entry_id"],
            "target_chain": e.get("target_chain"),
            **{k: r.get(k) for k in (
                "auth_seq_id", "uniprot_pos", "residue", "delta_sasa",
                "min_distance", "is_contact", "is_buried", "in_ectodomain",
            )},
        }
        for e in epitopes
        for r in e.get("residues", [])
    ]

    target_scores = {
        t["symbol"]: t.get("composite_score") for t in complexes["targets"]
    }

    epi_rows = []
    for e in epitopes:
        res_in = [r for r in e.get("residues", []) if r.get("in_ectodomain")]
        hotspots = e.get("hotspots", [])
        des = designability(
            e.get("buried_area_total") or 0.0,
            e.get("resolution"),
            len(hotspots),
            params.n_hotspots,
        )
        tscore = target_scores.get(e["symbol"]) or 0.0
        epi_rows.append(
            {
                "symbol": e["symbol"],
                "uniprot": e.get("uniprot"),
                "entry_id": e["entry_id"],
                "status": e.get("status"),
                "reason": e.get("reason"),
                "target_score": round(float(tscore), 4),
                **des,
                "combined_score": round(float(tscore) * des["designability"], 4)
                if e.get("status") == "accepted"
                else 0.0,
                "resolution": e.get("resolution"),
                "method": e.get("method"),
                "partner_symbol": e.get("partner_symbol"),
                "partner_symbols": ";".join(e.get("partner_symbols") or []),
                "target_chain": e.get("target_chain"),
                "partner_chain": e.get("partner_chain"),
                "n_interface_residues": e.get("n_interface_residues"),
                "buried_area_total": e.get("buried_area_total"),
                "fraction_in_ectodomain": e.get("fraction_in_ectodomain"),
                "hotspot_patch_radius": e.get("hotspot_patch_radius"),
                "interface_auth": ";".join(r["auth_seq_id"] for r in res_in),
                "interface_uniprot": ";".join(
                    str(r["uniprot_pos"]) for r in res_in if r.get("uniprot_pos")
                ),
                "hotspots_auth": ";".join(h["auth_seq_id"] for h in e.get("hotspots", [])),
                "hotspots_uniprot": ";".join(
                    str(h["uniprot_pos"]) for h in e.get("hotspots", [])
                ),
                "title": e.get("title"),
                "structure_file": e.get("structure_file"),
            }
        )

    epi_df = pd.DataFrame(epi_rows)
    if not epi_df.empty:
        epi_df = epi_df.sort_values(
            ["status", "combined_score"], ascending=[True, False]
        ).reset_index(drop=True)
        epi_df.insert(
            0,
            "epitope_rank",
            [
                i + 1 if s == "accepted" else None
                for i, s in enumerate(epi_df["status"])
            ],
        )
    res_df = pd.DataFrame(residue_rows)

    meta = {
        "schema_version": SCHEMA_VERSION,
        "package_version": __version__,
        "disease_id": complexes.get("disease_id"),
        "interface_params": {
            "contact_cutoff": params.contact_cutoff,
            "min_delta_sasa": params.min_delta_sasa,
            "sasa_probe_radius": params.sasa_probe_radius,
            "n_hotspots": params.n_hotspots,
            "ignore_hetero": params.ignore_hetero,
            "min_fraction_in_ectodomain": MIN_FRACTION_IN_ECTODOMAIN,
        },
        "max_complexes_per_target": max_complexes,
        "n_epitopes": len(epi_rows),
        "n_accepted": int((epi_df["status"] == "accepted").sum()) if not epi_df.empty else 0,
        "definitions": {
            "interface_residue": "heavy-atom contact within cutoff OR buried SASA >= min_delta_sasa",
            "hotspot": "interface residue in an extracellular span, making a contact, ranked by buried SASA",
            "numbering": "author numbering mapped to UniProt via SIFTS alignment; never assumed equal",
        },
    }
    return epi_df, res_df, meta


def write(
    epi_df: pd.DataFrame, res_df: pd.DataFrame, meta: dict, run_dir: Path
) -> tuple[Path, Path, Path]:
    run_dir.mkdir(parents=True, exist_ok=True)
    epi_path = run_dir / STAGE_FILES[4]
    res_path = run_dir / "stage4_interface_residues.csv"
    meta_path = run_dir / "stage4_interface_meta.json"
    epi_df.to_csv(epi_path, index=False)
    res_df.to_csv(res_path, index=False)
    meta_path.write_text(json.dumps(meta, indent=2))
    return epi_path, res_path, meta_path


def read(run_dir: Path) -> pd.DataFrame:
    path = run_dir / STAGE_FILES[4]
    if not path.exists():
        raise FileNotFoundError(f"{path} not found - run stage 4 first")
    return pd.read_csv(path)
