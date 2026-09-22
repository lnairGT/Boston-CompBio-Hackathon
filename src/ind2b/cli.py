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
