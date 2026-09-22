"""Stage 5 - emit BindCraft2 target and campaign specifications.

This stage writes files and stops. It never launches a design job: running
BindCraft2 needs a GPU host and real money, and that decision belongs to a
human looking at the epitope shortlist, not to a pipeline.

Schema
------
The layout follows BindCraft2's own convention, in which a *target definition*
and a *campaign* are separate files (verified against the upstream repository
rather than assumed):

``settings/target/<NAME>.json``::

    {"description": ..., "targets": [{"name": ..., "target_path": ...,
      "chains": "A", "hotspots": "54,56,66,115"}]}

``examples/<campaign>.json``::

    {"target": "<NAME>", "modality": "binder", "campaign_name": ...,
     "binder_lengths": [60, 100], "number_of_final_designs": 10,
     "project_folder": "results/..."}

Hotspots are emitted in the **author numbering of the PDB file written
alongside them**, not UniProt numbering, because that is the numbering
BindCraft2 reads from the structure. Both are recorded in the manifest so the
mapping is never lost.

The emitted structure contains the target chain only. The partner is
deliberately stripped: the binder must be designed against the free surface,
and leaving the partner in place would both occlude the epitope and let a
design score well by mimicking the partner's own contacts.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any

import pandas as pd
from Bio.PDB import MMCIFParser, PDBIO, Select

from . import __version__
from .config import STAGE_FILES, BinderSpec

log = logging.getLogger("ind2b.stage5")

SCHEMA_VERSION = 1

_SAFE_NAME = re.compile(r"[^A-Za-z0-9_]+")


class _ChainOnly(Select):
    """Keep one polymer chain, dropping waters, ligands and other chains."""

    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain):  # noqa: D102
        return chain.id == self.chain_id

    def accept_residue(self, residue):  # noqa: D102
        return residue.id[0] == " "

    def accept_atom(self, atom):  # noqa: D102
        # Drop alternate locations beyond the first to keep the file clean.
        return (not atom.is_disordered()) or atom.get_altloc() in (" ", "A")


def spec_name(symbol: str, entry_id: str, chain: str) -> str:
    return _SAFE_NAME.sub("_", f"{symbol}_{entry_id}_{chain}")


def write_target_structure(cif_path: Path, chain_id: str, dest: Path) -> int:
    """Write a PDB containing only ``chain_id``; returns the residue count."""
    structure = MMCIFParser(QUIET=True).get_structure("t", str(cif_path))
    model = next(iter(structure))
    if chain_id not in [c.id for c in model]:
        raise ValueError(f"chain {chain_id} absent from {cif_path.name}")
    n_res = sum(1 for r in model[chain_id] if r.id[0] == " ")
    dest.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(dest), select=_ChainOnly(chain_id))
    return n_res


def build_specs(
    epitopes: pd.DataFrame,
    *,
    run_dir: Path,
    binder: BinderSpec | None = None,
    max_trajectories: int | None = None,
    one_per_target: bool = True,
    top_n: int | None = None,
) -> dict[str, Any]:
    """Emit target/campaign specs for the accepted epitopes."""
    binder = binder or BinderSpec()
    accepted = epitopes[epitopes["status"] == "accepted"].copy()
    if accepted.empty:
        return {"n_specs": 0, "specs": [], "note": "no accepted epitopes"}

    accepted = accepted.sort_values("combined_score", ascending=False)
    if one_per_target:
        accepted = accepted.groupby("symbol", as_index=False).first()
        accepted = accepted.sort_values("combined_score", ascending=False)
    if top_n:
        accepted = accepted.head(top_n)

    out_dir = run_dir / STAGE_FILES[5]
    struct_out = out_dir / "structures"
    target_out = out_dir / "settings" / "target"
    camp_out = out_dir / "campaigns"
    for d in (struct_out, target_out, camp_out):
        d.mkdir(parents=True, exist_ok=True)

    specs: list[dict[str, Any]] = []
    for _, row in accepted.iterrows():
        name = spec_name(row["symbol"], row["entry_id"], row["target_chain"])
        pdb_path = struct_out / f"{name}.pdb"

        # Stage 4 records paths relative to the run directory; fall back to the
        # literal value for records written before that convention.
        recorded = str(row["structure_file"])
        cif = run_dir / recorded
        if not cif.exists():
            cif = Path(recorded)

        if not cif.exists():
            specs.append({"name": name, "status": "skipped",
                           "reason": f"structure missing: {cif}"})
            continue
        try:
            n_res = write_target_structure(cif, row["target_chain"], pdb_path)
        except Exception as exc:  # noqa: BLE001
            specs.append({"name": name, "status": "skipped",
                           "reason": f"{type(exc).__name__}: {exc}"})
            continue

        hotspots_auth = [h for h in str(row["hotspots_auth"]).split(";") if h]
        if not hotspots_auth:
            specs.append({"name": name, "status": "skipped",
                           "reason": "no hotspots selected"})
            continue

        target_def = {
            "description": (
                f"{row['symbol']} ({row['uniprot']}) - interface with "
                f"{row['partner_symbol']} from PDB {row['entry_id']} "
                f"chain {row['target_chain']}, "
                f"{row['buried_area_total']:.0f} A^2 buried. Hotspots are a "
                f"compact patch in the author numbering of this structure."
            ),
            "targets": [
                {
                    "name": name,
                    "target_path": f"structures/{name}.pdb",
                    "chains": str(row["target_chain"]),
                    "hotspots": ",".join(hotspots_auth),
                }
            ],
        }
        campaign = {
            "target": name,
            "modality": "binder",
            "campaign_name": name.lower(),
            "binder_lengths": [binder.binder_length_min, binder.binder_length_max],
            "number_of_final_designs": binder.n_designs,
            "project_folder": f"results/{name.lower()}",
        }
        if max_trajectories:
            campaign["max_trajectories"] = max_trajectories

        (target_out / f"{name}.json").write_text(json.dumps(target_def, indent=2))
        (camp_out / f"{name}.json").write_text(json.dumps(campaign, indent=2))

        specs.append(
            {
                "name": name,
                "status": "written",
                "symbol": row["symbol"],
                "uniprot": row["uniprot"],
                "pdb_entry": row["entry_id"],
                "partner": row["partner_symbol"],
                "target_chain": row["target_chain"],
                "resolution": row["resolution"],
                "buried_area": row["buried_area_total"],
                "n_structure_residues": n_res,
                "hotspots_auth": ",".join(hotspots_auth),
                "hotspots_uniprot": str(row["hotspots_uniprot"]),
                "combined_score": row["combined_score"],
                "target_file": str(target_out / f"{name}.json"),
                "campaign_file": str(camp_out / f"{name}.json"),
                "structure_file": str(pdb_path),
            }
        )
        log.info("spec written: %s", name)

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "package_version": __version__,
        "n_specs": sum(1 for s in specs if s["status"] == "written"),
        "binder_lengths": [binder.binder_length_min, binder.binder_length_max],
        "number_of_final_designs": binder.n_designs,
        "max_trajectories": max_trajectories,
        "numbering": (
            "hotspots are in the author numbering of the emitted PDB; the "
            "UniProt equivalents are in hotspots_uniprot"
        ),
        "partner_stripped": True,
        "not_submitted": (
            "These are specifications only. No design job has been run and none "
            "is launched by this package."
        ),
        "specs": specs,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    _write_readme(out_dir, manifest)
    return manifest


def _write_readme(out_dir: Path, manifest: dict) -> None:
    written = [s for s in manifest["specs"] if s["status"] == "written"]
    lines = [
        "# BindCraft2 specifications",
        "",
        f"{len(written)} target/campaign pairs emitted by `ind2b stage5`.",
        "**Nothing here has been run.** These are inputs for a design campaign.",
        "",
        "## Layout",
        "```",
        "settings/target/<NAME>.json   target definition (chain + hotspots)",
        "campaigns/<NAME>.json         campaign settings",
        "structures/<NAME>.pdb         target chain only, partner stripped",
        "```",
        "",
        "## Running one",
        "BindCraft2 needs Linux, Python 3.12+ and an NVIDIA GPU; there is no CPU",
        "path. Place the target definition where BindCraft2 looks for targets and",
        "the structure at the `target_path` it names, then:",
        "```bash",
        "bindcraft design campaigns/<NAME>.json",
        "```",
        "Calibrate before committing budget: cost per campaign is dominated by the",
        "target-dependent acceptance rate, which is not knowable in advance.",
        "",
        "## Numbering",
        "Hotspots are in the **author numbering of the emitted PDB**. UniProt",
        "equivalents are in `manifest.json` as `hotspots_uniprot`.",
        "",
        "## Specs",
        "",
        "| name | target | partner | PDB | buried A^2 | hotspots (auth) |",
        "|------|--------|---------|-----|-----------:|-----------------|",
    ]
    for s in written:
        lines.append(
            f"| {s['name']} | {s['symbol']} | {s['partner']} | {s['pdb_entry']} "
            f"| {s['buried_area']:.0f} | {s['hotspots_auth']} |"
        )
    (out_dir / "README.md").write_text("\n".join(lines) + "\n")


def read(run_dir: Path) -> dict[str, Any]:
    path = run_dir / STAGE_FILES[5] / "manifest.json"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found - run stage 5 first")
    return json.loads(path.read_text())
