"""Command-line entry point.

Each pipeline stage is a subcommand that reads the previous stage's output
from the run directory and writes its own. Stages are added to this CLI as
they are implemented; ``ind2b --help`` lists what is currently wired.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from . import __version__
from .config import RUNS_DIR, RunConfig


def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "--run-dir",
        type=Path,
        default=None,
        help="run directory (default: runs/<efo_id or indication slug>)",
    )
    p.add_argument("-v", "--verbose", action="store_true", help="debug logging")


def _resolve_run_dir(args, efo_id: str | None = None, indication: str | None = None) -> Path:
    if args.run_dir:
        return args.run_dir
    cfg = RunConfig(indication=indication or "", efo_id=efo_id)
    return cfg.run_dir(RUNS_DIR)


def cmd_stage0(args) -> int:
    from . import stage0_resolve

    record = stage0_resolve.resolve(
        args.indication, efo_id=args.efo_id, n_candidates=args.n_candidates
    )
    run_dir = _resolve_run_dir(
        args, efo_id=record["disease"]["id"], indication=args.indication
    )
    path = stage0_resolve.write(record, run_dir)

    d = record["disease"]
    print(f"resolved: {d['id']}  {d['name']}")
    print(f"  rule           : {record['selection_rule']}")
    print(f"  ambiguous      : {record['ambiguous']}")
    print(f"  therapeutic    : {', '.join(a['name'] for a in d['therapeutic_areas']) or '-'}")
    print(f"  children       : {len(record['children'])}")
    print(f"  descendants    : {record['descendants']['count']}")
    print(f"  written        : {path}")
    if record["ambiguous"]:
        print("\n  candidates considered:", file=sys.stderr)
        for c in record["candidates_considered"]:
            print(f"    {c['id']:<16} {c['name']}", file=sys.stderr)
        print("  pin the term with --efo-id if the choice is wrong.", file=sys.stderr)
    return 0


def cmd_stage1(args) -> int:
    from . import stage0_resolve, stage1_evidence

    run_dir = _resolve_run_dir(args, efo_id=args.efo_id)
    d0 = stage0_resolve.read(run_dir)
    record = stage1_evidence.fetch(
        d0["disease"]["id"],
        disease_name=d0["disease"]["name"],
        n_targets=args.n_targets,
        all_layers=args.all_layers,
    )
    path = stage1_evidence.write(record, run_dir)

    print(f"disease   : {d0['disease']['id']}  {d0['disease']['name']}")
    print(f"pool      : {record['pool_size']} of {record['association_total']} associated targets")
    print(f"reachable : {record['n_surface_accessible']} (UniProt topology)")
    disagree = [
        t for t in record["targets"]
        if not t["location_keyword_hint"]["agrees_with_topology"]
    ]
    print(f"  keyword-hint disagreements: {len(disagree)}")
    print(f"written   : {path}")
    return 0


def cmd_stage2(args) -> int:
    from . import stage1_evidence, stage2_score

    run_dir = _resolve_run_dir(args, efo_id=args.efo_id)
    evidence = stage1_evidence.read(run_dir)
    df, weights = stage2_score.score_targets(evidence, weights_file=args.weights)
    csv_path, meta_path = stage2_score.write(df, weights, evidence, run_dir)

    passing = df[df["passes_accessibility"]]
    print(f"scored {len(df)} targets; {len(passing)} pass the accessibility gate")
    print(f"\n{'#':>3} {'symbol':<10} {'score':>7}  {'ectodomain':>10}  class")
    for _, row in passing.head(args.top).iterrows():
        cls = (row["target_class"] or "").split(";")[0][:28]
        print(f"{int(row['rank']):>3} {row['symbol']:<10} {row['composite_score']:>7.3f}  "
              f"{int(row['ectodomain_residues'] or 0):>7} aa  {cls}")
    print(f"\nwritten: {csv_path}\n         {meta_path}")
    return 0


def cmd_stage3(args) -> int:
    from . import stage1_evidence, stage2_score, stage3_complexes

    run_dir = _resolve_run_dir(args, efo_id=args.efo_id)
    evidence = stage1_evidence.read(run_dir)
    ranked = stage2_score.read(run_dir)
    record = stage3_complexes.find_complexes(
        ranked, evidence, rows_per_target=args.rows_per_target
    )
    path = stage3_complexes.write(record, run_dir)

    print(f"targets in: {record['n_targets_in']} | admitted: {record['n_targets_admitted']}")
    for t in record["targets"]:
        if t["n_admitted"]:
            best = t["admitted"][0]
            print(f"  {t['symbol']:<9} {t['n_admitted']:>3} complexes  best {best['entry_id']} "
                  f"{best['resolution']}A + {','.join(best['partner_symbols']) or '?'}")
        else:
            print(f"  {t['symbol']:<9} dropped: {t['drop_reason']} {t['rejection_reason_counts']}")
    print(f"\nwritten: {path}")
    return 0


def cmd_stage4(args) -> int:
    from . import stage3_complexes, stage4_interface

    run_dir = _resolve_run_dir(args, efo_id=args.efo_id)
    complexes = stage3_complexes.read(run_dir)
    epi_df, res_df, meta = stage4_interface.run(
        complexes, run_dir=run_dir, max_complexes=args.max_complexes
    )
    epi_path, res_path, meta_path = stage4_interface.write(epi_df, res_df, meta, run_dir)

    accepted = epi_df[epi_df["status"] == "accepted"]
    print(f"epitopes: {len(epi_df)} computed | {len(accepted)} accepted "
          f"| {accepted['symbol'].nunique()} targets")
    print(f"\n{'#':>3} {'symbol':<9} {'PDB':<6} {'partner':<9} {'buried':>7}  hotspots")
    for _, r in accepted.head(args.top).iterrows():
        print(f"{int(r['epitope_rank']):>3} {r['symbol']:<9} {r['entry_id']:<6} "
              f"{str(r['partner_symbol'])[:8]:<9} {r['buried_area_total']:>7.0f}  {r['hotspots_auth']}")
    print(f"\nwritten: {epi_path}\n         {res_path}\n         {meta_path}")
    return 0


def cmd_stage5(args) -> int:
    from . import stage4_interface, stage5_specs
    from .config import BinderSpec

    run_dir = _resolve_run_dir(args, efo_id=args.efo_id)
    epitopes = stage4_interface.read(run_dir)
    binder = BinderSpec(
        binder_length_min=args.binder_min,
        binder_length_max=args.binder_max,
        n_designs=args.n_designs,
    )
    manifest = stage5_specs.build_specs(
        epitopes,
        run_dir=run_dir,
        binder=binder,
        max_trajectories=args.max_trajectories,
        one_per_target=not args.all_epitopes,
        top_n=args.top_n,
    )
    print(f"{manifest['n_specs']} BindCraft2 specification(s) written to "
          f"{run_dir / 'stage5_specs'}")
    for s in manifest["specs"]:
        if s["status"] == "written":
            print(f"  {s['name']:<22} {s['pdb_entry']} vs {s['partner']:<9} "
                  f"hotspots {s['hotspots_auth']}")
        else:
            print(f"  {s['name']:<22} SKIPPED: {s['reason']}")
    print("\nNothing has been submitted. Review the specs, then run BindCraft2 "
          "yourself on a GPU host.")
    return 0


def cmd_stage6(args) -> int:
    from . import stage1_evidence, stage3_complexes, stage4_interface, stage6_binders

    run_dir = _resolve_run_dir(args, efo_id=args.efo_id)
    complexes = stage3_complexes.read(run_dir)
    epitopes = stage4_interface.read(run_dir)
    evidence = stage1_evidence.read(run_dir)
    record = stage6_binders.run(
        complexes, epitopes, evidence, run_dir=run_dir, max_per_kind=args.max_per_kind
    )
    record = stage6_binders.add_projections(
        record, run_dir, epitopes, max_examples=args.max_projections
    )
    path = stage6_binders.write(record, run_dir)

    print(f"{record['n_with_designed']} of {record['n_targets']} targets have a "
          f"designed-binder complex in the PDB")
    for t in record["targets"]:
        for b in t["binders"]:
            cmp_ = b.get("vs_selected_epitope") or {}
            shared = (f'{cmp_["n_shared"]} shared, '
                      f'{cmp_["fraction_of_hotspots_covered"]:.0%} of hotspots'
                      if cmp_.get("comparable") and
                      cmp_.get("fraction_of_hotspots_covered") is not None
                      else "not comparable")
            print(f"  {t['symbol']:<8} {b['kind']:<9} {b['entry_id']} "
                  f"{b['resolution']}A  binder {b['binder_length']}aa  -> {shared}")
    print(f"\nNo binder here was designed by this pipeline. {path}")
    return 0


def cmd_schema_check(args) -> int:
    """Verify the Open Targets field names this package queries still exist."""
    from .sources import opentargets as ot

    expected = {
        "Disease": ["id", "name", "description", "therapeuticAreas", "parents",
                     "children", "descendants", "associatedTargets"],
        "Target": ["id", "approvedSymbol", "approvedName", "biotype", "proteinIds",
                    "functionDescriptions", "subcellularLocations", "tractability",
                    "safetyLiabilities", "geneticConstraint", "interactions",
                    "drugAndClinicalCandidates"],
    }
    failures = 0
    for type_name, fields in expected.items():
        live = set(ot.introspect_type(type_name))
        missing = [f for f in fields if f not in live]
        status = "ok" if not missing else f"MISSING {missing}"
        print(f"{type_name:<10} {len(live):>3} fields  {status}")
        failures += bool(missing)
    if failures:
        print("\nschema drift detected - update src/ind2b/sources/opentargets.py",
              file=sys.stderr)
    return 1 if failures else 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="ind2b",
        description="Indication -> ranked targets -> epitope -> BindCraft2 spec.",
    )
    p.add_argument("--version", action="version", version=f"indication2binder {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    s0 = sub.add_parser("stage0", help="resolve an indication to a disease ontology id")
    s0.add_argument("--indication", default="", help="free-text disease name")
    s0.add_argument("--efo-id", default=None, help="pin the ontology id, skipping text search")
    s0.add_argument("--n-candidates", type=int, default=10, help="ontology hits to consider")
    _add_common(s0)
    s0.set_defaults(func=cmd_stage0)

    s1 = sub.add_parser("stage1", help="gather evidence layers for the candidate targets")
    s1.add_argument("--efo-id", default=None, help="disease id (locates the run directory)")
    s1.add_argument("--n-targets", type=int, default=100, help="size of the candidate pool")
    s1.add_argument(
        "--all-layers",
        action="store_true",
        help="fetch per-target layers for unreachable targets too (slower)",
    )
    _add_common(s1)
    s1.set_defaults(func=cmd_stage1)

    s2 = sub.add_parser("stage2", help="score and rank the candidate targets")
    s2.add_argument("--efo-id", default=None, help="disease id (locates the run directory)")
    s2.add_argument("--weights", default=None, help="JSON file overriding scoring weights")
    s2.add_argument("--top", type=int, default=15, help="rows to print")
    _add_common(s2)
    s2.set_defaults(func=cmd_stage2)

    s3 = sub.add_parser("stage3", help="gate on experimentally determined complexes")
    s3.add_argument("--efo-id", default=None, help="disease id (locates the run directory)")
    s3.add_argument("--rows-per-target", type=int, default=200,
                     help="max PDB entries to inspect per target")
    _add_common(s3)
    s3.set_defaults(func=cmd_stage3)

    s4 = sub.add_parser("stage4", help="compute interfaces and select hotspots")
    s4.add_argument("--efo-id", default=None, help="disease id (locates the run directory)")
    s4.add_argument("--max-complexes", type=int, default=2,
                     help="complexes to analyse per target")
    s4.add_argument("--top", type=int, default=20, help="rows to print")
    _add_common(s4)
    s4.set_defaults(func=cmd_stage4)

    s5 = sub.add_parser("stage5", help="emit BindCraft2 specs (writes files, runs nothing)")
    s5.add_argument("--efo-id", default=None, help="disease id (locates the run directory)")
    s5.add_argument("--binder-min", type=int, default=55, help="minimum binder length")
    s5.add_argument("--binder-max", type=int, default=120, help="maximum binder length")
    s5.add_argument("--n-designs", type=int, default=10, help="final designs per campaign")
    s5.add_argument("--max-trajectories", type=int, default=None,
                     help="cap trajectories per campaign (useful for calibration runs)")
    s5.add_argument("--all-epitopes", action="store_true",
                     help="emit every accepted epitope, not just the best per target")
    s5.add_argument("--top-n", type=int, default=None, help="emit only the top N")
    _add_common(s5)
    s5.set_defaults(func=cmd_stage5)

    s6 = sub.add_parser(
        "stage6", help="find existing binders (designed and antibody) of the targets")
    s6.add_argument("--efo-id", default=None, help="disease id (locates the run directory)")
    s6.add_argument("--max-per-kind", type=int, default=1,
                     help="complexes to analyse per binder kind per target")
    s6.add_argument("--max-projections", type=int, default=6,
                     help="example complexes to project for drawing")
    _add_common(s6)
    s6.set_defaults(func=cmd_stage6)

    sc = sub.add_parser("schema-check", help="verify upstream GraphQL field names")
    _add_common(sc)
    sc.set_defaults(func=cmd_schema_check)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if getattr(args, "verbose", False) else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
