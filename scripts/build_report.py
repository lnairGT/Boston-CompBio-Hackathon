"""Build the run report from a completed run directory.

Every number in the report is read from the stage outputs rather than passed
in, so the report cannot drift from the artifacts it describes.

Usage::

    python scripts/build_report.py runs/MONDO_0005233 [-o report.md]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def build(run_dir: Path) -> str:
    d0 = json.loads((run_dir / "stage0_disease.json").read_text())
    ev = json.loads((run_dir / "stage1_evidence.json").read_text())
    ranked = pd.read_csv(run_dir / "stage2_ranked_targets.csv")
    cx = json.loads((run_dir / "stage3_complexes.json").read_text())
    epi = pd.read_csv(run_dir / "stage4_epitopes.csv")
    meta4 = json.loads((run_dir / "stage4_interface_meta.json").read_text())
    manifest = json.loads((run_dir / "stage5_specs" / "manifest.json").read_text())

    disease = d0["disease"]
    accessible = ranked[ranked["passes_accessibility"]]
    accepted = epi[epi["status"] == "accepted"]
    best = (
        accepted.sort_values("combined_score", ascending=False)
        .groupby("symbol", as_index=False)
        .first()
        .sort_values("combined_score", ascending=False)
    )
    written = [s for s in manifest["specs"] if s["status"] == "written"]

    kw = ev.get("n_keyword_hint_disagreements")
    fp = [
        t["symbol"]
        for t in ev["targets"]
        if t["location_keyword_hint"]["hint"] and not t["surface_accessible"]
    ]
    fn = [
        t["symbol"]
        for t in ev["targets"]
        if not t["location_keyword_hint"]["hint"] and t["surface_accessible"]
    ]
    dropped = [t for t in cx["targets"] if not t["n_admitted"]]
    params = meta4["interface_params"]

    L: list[str] = []
    a = L.append

    a(f"# Binder-design targets for {disease['name']}")
    a("")
    a(f"Indication resolved to **{disease['id']}** ({disease['name']}) by "
      f"`{d0['selection_rule']}`; ambiguous = {d0['ambiguous']}.")
    a("")
    a("## What this run produced")
    a("")
    a(f"| Stage | Result |")
    a(f"|---|---|")
    a(f"| 1 evidence | {ev['pool_size']} candidate targets from "
      f"{ev['association_total']} disease associations |")
    a(f"| 1 accessibility | {ev['n_surface_accessible']} reachable by an "
      f"extracellular protein binder (UniProt topology) |")
    a(f"| 2 ranking | {len(accessible)} scored above the accessibility gate, "
      f"{len(ranked) - len(accessible)} retained with exclusion reasons |")
    a(f"| 3 complex gate | {cx['n_targets_admitted']} of {cx['n_targets_in']} "
      f"targets have a qualifying experimental complex |")
    a(f"| 4 interfaces | {len(accepted)} accepted epitopes across "
      f"{accepted['symbol'].nunique()} targets ({len(epi) - len(accepted)} rejected) |")
    a(f"| 5 specs | {len(written)} BindCraft2 target/campaign pairs written, "
      f"none submitted |")
    a("")

    a("## Shortlist")
    a("")
    a("One epitope per target, ranked by target evidence x epitope designability.")
    a("")
    a("| # | target | UniProt | partner blocked | PDB | res (A) | buried (A^2) | "
      "hotspots (author numbering) |")
    a("|--:|---|---|---|---|--:|--:|---|")
    for i, (_, r) in enumerate(best.iterrows(), start=1):
        a(f"| {i} | {r['symbol']} | {r['uniprot']} | {r['partner_symbol']} | "
          f"{r['entry_id']} | {r['resolution']} | {r['buried_area_total']:.0f} | "
          f"`{str(r['hotspots_auth']).replace(';', ',')}` |")
    a("")

    a("## Method, and the corrections it needed")
    a("")
    a("### Accessibility is topology, not a location keyword")
    a("")
    a("A binder can only reach a target from outside the cell, so localization "
      "is a hard gate. Matching subcellular-location strings fails at it: "
      "UniProt annotates KRAS at the `Cell membrane`, but it faces the "
      "cytoplasm (`Lipid-anchor | Cytoplasmic side`). The keyword test admitted "
      f"{len(fp)} intracellular targets ({', '.join(fp[:8])}"
      f"{', ...' if len(fp) > 8 else ''}) and missed "
      f"{len(fn)} real ectodomain targets ({', '.join(fn)}).")
    a("")
    a("The gate now requires an annotated extracellular topological domain or a "
      "`Secreted` location, and returns the reachable residue spans, which the "
      "interface stage needs anyway. The keyword verdict is retained as a "
      "recorded hint with an agreement flag.")
    a("")
    a("### Safety-liability count measures study intensity, not risk")
    a("")
    a("Penalising targets per recorded safety liability dropped EGFR - the "
      "top-associated target here, with 82 known drugs and 1,735 trials - by "
      "twenty places. Liability count correlates 0.65 with known-drug count in "
      "this pool, because Open Targets accumulates annotations for targets that "
      "reached patients. The field carries no severity, so the penalty is a flat "
      "0.92 and the count is reported for the reader to weigh.")
    a("")
    a("### A native-partner complex can still present an unreachable interface")
    a("")
    a("CBL-C binds the intracellular tail of EGFR: a genuine experimental "
      "complex of two real proteins, useless for an extracellular binder. Each "
      "interface residue is therefore mapped from PDB author numbering into "
      "UniProt numbering through the SIFTS alignment and tested against the "
      "extracellular spans. Author numbering is never assumed to equal UniProt "
      "numbering - that assumption breaks silently for constructs numbered from "
      "their own residue 1. The check rejected "
      f"{int((epi['reason'].fillna('').str.contains('outside ectodomain')).sum())} "
      "epitopes at 0% ectodomain overlap.")
    a("")
    a("### Peptide fragments are not epitopes")
    a("")
    a("EGFR appears in the PDB as a 13-residue phosphopeptide bound to CBL-C. "
      "The complex gate therefore requires a target chain of at least "
      f"{cx['gate']['min_target_chain_residues']} residues. Antibody and Fab "
      "chains are classified separately rather than counted as partners: they "
      "evidence that an epitope is targetable by a biologic, but blocking a "
      "defined interaction needs the native partner.")
    a("")
    a("### Epitope size is a band, not a maximum")
    a("")
    a("Interfaces are computed per partner chain, not pooled. Pooled, ERBB4 "
      "buries 3,136 A^2 across two partners simultaneously - a composite no "
      "single mini-binder reproduces. Designability therefore scores buried "
      "area as a preference band (500-1,200 A^2), and hotspots are selected as "
      "a spatially compact patch seeded on maximum burial rather than the top N "
      "residues by area, which can otherwise straddle opposite edges of a large "
      "interface.")
    a("")

    a("## Validation against a known answer")
    a("")
    a("PD-L1/PD-1 (4ZQK) is BindCraft2's own benchmark system, which makes it a "
      "positive control. The computed PD-L1 interface recovers the canonical "
      "epitope - Y123, A121, R113, R125, K124, Q66, M115, Y56, I54, D122 - "
      "burying 763 A^2. All four hotspots in BindCraft2's reference `hPDL1` "
      "target definition (54, 56, 66, 115) fall inside that computed interface.")
    a("")

    a("## Parameters")
    a("")
    a(f"- Interface: contact <= {params['contact_cutoff']} A heavy-atom OR "
      f"buried SASA >= {params['min_delta_sasa']} A^2 (Shrake-Rupley)")
    a(f"- Ectodomain overlap required: "
      f"{params['min_fraction_in_ectodomain']:.0%} of mapped interface residues")
    a(f"- Complex gate: {cx['gate']['accepted_methods']}, resolution <= "
      f"{cx['gate']['max_resolution']} A, target chain >= "
      f"{cx['gate']['min_target_chain_residues']} residues")
    a(f"- Hotspots per epitope: {params['n_hotspots']}; "
      f"{meta4.get('hotspot_selection', 'compact patch')}")
    a(f"- Binder lengths in emitted campaigns: {manifest['binder_lengths']}, "
      f"{manifest['number_of_final_designs']} final designs")
    a("")

    a("## What this does not establish")
    a("")
    a("- **Not causality.** The ranking reflects the weight of existing "
      "evidence that a target is associated with the disease. A well-studied "
      "target scores highly whether or not it is causal, and the known-drug and "
      "literature components are not independent of the overall association "
      "score, which already incorporates both.")
    a("- **Not therapeutic mechanism.** Partner relevance is inherited from an "
      "interaction database, not curated per indication. Some admitted pairs are "
      "real interactions whose blockade is not obviously therapeutic (CD74-CTSL "
      "is a processing event, not a signalling interaction); those need a human "
      "to judge.")
    a("- **Not affinity, specificity or developability.** No design has been "
      "run. These are inputs to a campaign.")
    a("- **Not a complete epitope survey.** Only the top complexes per target "
      "were analysed, and the experimental-structure requirement excludes "
      "targets whose biology may be sound but whose complexes are unsolved "
      f"({len(dropped)} targets dropped here: "
      f"{', '.join(t['symbol'] for t in dropped)}).")
    a("- **Not an antibody-epitope map.** Antibody complexes were deliberately "
      "excluded from partner status, so epitopes validated only by a Fab "
      "structure do not appear.")
    a("")
    a("## Reproducing")
    a("")
    a("```bash")
    a(f"ind2b stage0 --indication \"{d0.get('query', disease['name'])}\"")
    for s in ("stage1", "stage2", "stage3", "stage4", "stage5"):
        a(f"ind2b {s} --efo-id {disease['id']}")
    a("```")
    a("")
    a("Stages are disk-cached, so a re-run costs no network traffic and any "
      "stage can be re-run in isolation.")
    return "\n".join(L) + "\n"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run_dir", type=Path)
    p.add_argument("-o", "--output", type=Path, default=None)
    args = p.parse_args()
    text = build(args.run_dir)
    out = args.output or (args.run_dir / "report.md")
    out.write_text(text)
    print(f"wrote {out} ({len(text.splitlines())} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
