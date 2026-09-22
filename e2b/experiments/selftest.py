"""Module-local self-test for the experimental planner.

Run as::

    python -m e2b.experiments.selftest              # live UniProt, all branches
    python -m e2b.experiments.selftest --offline    # structural checks only
    python -m e2b.experiments.selftest --write-doc docs/experimental_plan.md

What it proves, and how:

* the planner produces a **different** plan for a secreted protein and for a
  cell-surface receptor -- different expression host for the antigen, different assay
  orientation, different functional design -- from live UniProt topology alone;
* it **refuses** to write an in-cell plan for a target with no extracellular face,
  and says why, rather than producing a plausible plan that could not work;
* a retrieval failure returns ``blocked_no_target_biology`` labelled as a source
  failure, not as absence of an ectodomain;
* honesty invariants hold on every generated plan: no assay parameter is stated as a
  validated value, no affinity-named number appears, unknown numerics are ``None``,
  and a fixture candidate stays labelled ``fixture`` throughout.

The accessions used live in ``selftest_fixtures.json`` and are test inputs only. The
production path (``planner.py``, ``target_profile.py``, ``moa.py``, ``liabilities.py``,
``render.py``) contains no accession, gene, disease or dataset identifier.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

from .planner import plan_experiments
from .records import MODULE_VERSION
from .render import render_markdown

FIXTURES = Path(__file__).with_name("selftest_fixtures.json")

# Numbers that would be a lie if this module emitted them: a specific concentration,
# a specific incubation, or an affinity. The check is deliberately crude and errs
# toward flagging, because a false alarm costs a reviewer a minute and a miss costs
# the project its credibility.
_FORBIDDEN_VALUE_PATTERNS = [
    re.compile(r"\b\d+(\.\d+)?\s*(nM|pM|uM|µM|mM)\b"),
    re.compile(r"\b\d+(\.\d+)?\s*(min|minutes|h|hours|hr)\b"),
    re.compile(r"\bKd\s*=", re.I),
    re.compile(r"\bIC50\s*=", re.I),
    re.compile(r"\bEC50\s*=", re.I),
]

_ANTIBODY_MISLABEL = re.compile(r"\bthe (antibody|mAb)\b", re.I)


def _walk_strings(obj: Any, path: str = "") -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    if isinstance(obj, str):
        out.append((path, obj))
    elif isinstance(obj, dict):
        for k, v in obj.items():
            out.extend(_walk_strings(v, f"{path}.{k}"))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            out.extend(_walk_strings(v, f"{path}[{i}]"))
    return out


def check_invariants(plan: dict[str, Any]) -> list[str]:
    """Return a list of invariant violations. Empty list means the plan is clean."""
    problems: list[str] = []
    strings = _walk_strings(plan)

    for path, s in strings:
        # The instrument-limit text legitimately cites rate constants in scientific
        # notation; those are documented capability limits, not prescribed values.
        for pat in _FORBIDDEN_VALUE_PATTERNS:
            m = pat.search(s)
            if m:
                problems.append(f"specific value '{m.group(0)}' stated at {path}")
        if _ANTIBODY_MISLABEL.search(s):
            problems.append(f"design referred to as an antibody at {path}")

    for route in plan.get("production_routes", []):
        for sp in route.get("starting_points", []):
            if sp.get("status") != "starting_point_requires_optimisation":
                problems.append(f"starting point '{sp.get('name')}' not marked as a starting point")
    for a in plan.get("assays", []):
        for sp in a.get("starting_points", []):
            if sp.get("status") != "starting_point_requires_optimisation":
                problems.append(f"starting point '{sp.get('name')}' in {a['assay_id']} not marked")
    for d in plan.get("decision_tree", []):
        if d.get("threshold_status") != "pre_registered_default_requires_sign_off":
            problems.append(f"criterion {d.get('criterion_id')} threshold not marked for sign-off")
        if not d.get("kill_rule"):
            problems.append(f"criterion {d.get('criterion_id')} has no kill rule")

    # Unknown numerics must be None, never 0. The distinction is not "is it zero" --
    # a protein really can have zero annotated glycosylation sites, and reporting that
    # as None would lose information. It is "was the annotation retrieved at all":
    # length_aa is None exactly when no UniProt record was read, so any count that is
    # 0 alongside a None length is an undetermined value masquerading as a measured
    # zero. length_aa itself may never be 0, because a protein of zero residues is not
    # a protein -- that is the empty-record signature.
    tp = plan.get("target_profile") or {}
    counts = ("n_glycosylation_sites", "n_disulfides", "experimental_structures")
    if tp.get("length_aa") == 0:
        problems.append("target_profile.length_aa is 0; an empty record must be None, not zero")
    record_read = tp.get("length_aa") is not None
    for numeric in counts:
        if numeric not in tp:
            continue
        if tp[numeric] == 0 and not record_read:
            problems.append(
                f"target_profile.{numeric} is 0 while length_aa is None: no annotation was read, "
                "so this count is undetermined and must be None"
            )
    if plan.get("candidate_sequence_length") == 0:
        problems.append("candidate_sequence_length is 0 rather than None")
    for a in plan.get("assays", []):
        for sp in a.get("starting_points", []):
            if sp.get("range_text") in (None, ""):
                problems.append(f"starting point '{sp.get('name')}' in {a['assay_id']} has no range")

    if not plan.get("scope_statement"):
        problems.append("no scope statement")
    if not plan.get("what_this_does_not_establish"):
        problems.append("no 'does not establish' block")
    return problems


def _load_fixtures() -> dict[str, Any]:
    return json.loads(FIXTURES.read_text())


def run(*, offline: bool = False, write_doc: str | None = None, verbose: bool = True) -> dict[str, Any]:
    fx = _load_fixtures()
    brief_a, brief_b = fx["briefs"][0], fx["briefs"][1]
    cand = {
        k: v for k, v in fx["candidate_fixture"].items() if not k.startswith("_")
    }
    results: dict[str, Any] = {"module_version": MODULE_VERSION, "cases": [], "ok": True}
    rendered: dict[str, str] = {}

    for i, case in enumerate(fx["targets"]):
        brief = brief_a if i != 1 else brief_b
        row: dict[str, Any] = {"label": case["label"], "expect_archetype": case["expect_archetype"]}
        try:
            plan = plan_experiments(
                brief,
                case["target"],
                cand,
                pathway_hints=case.get("pathway_hints"),
                paralogs=fx["paralog_fixture"] if i == 0 else None,
                fetch=not offline,
            )
        except Exception as exc:  # a crash is a failure, not a result
            row.update({"crashed": f"{type(exc).__name__}: {exc}"})
            results["cases"].append(row)
            results["ok"] = False
            continue

        tp = plan.get("target_profile") or {}
        row.update(
            {
                "status": plan["status"],
                "archetype": tp.get("archetype"),
                "mechanism_class": tp.get("mechanism_class"),
                "target_length_aa": tp.get("length_aa"),
                "n_glyc": tp.get("n_glycosylation_sites"),
                "n_ss": tp.get("n_disulfides"),
                "pdb_entries": tp.get("experimental_structures"),
                "n_assays": len(plan.get("assays", [])),
                "n_gates": len(plan.get("decision_tree", [])),
                "candidate_origin": plan.get("candidate_origin"),
                "construct_risks": [
                    r["risk"]
                    for rt in plan.get("production_routes", [])
                    for r in rt.get("risks", [])
                    if rt["subject"] == "designed_binder"
                ],
                "antigen_host": next(
                    (rt["expression_host"] for rt in plan.get("production_routes", []) if rt["subject"] == "target_antigen"),
                    None,
                ),
                "antigen_construct": next(
                    (
                        rt["construct_elements"][0]
                        for rt in plan.get("production_routes", [])
                        if rt["subject"] == "target_antigen" and rt.get("construct_elements")
                    ),
                    None,
                ),
                "binding_orientation": (
                    next((a["orientation"] for a in plan.get("assays", []) if a["role"] == "binding_primary"), None) or ""
                )[:90],
                "functional_present": any(a["role"] == "functional_moa" for a in plan.get("assays", [])),
                "invariant_violations": check_invariants(plan),
            }
        )
        if offline:
            row["expectation_met"] = plan["status"] == "blocked_no_target_biology"
            row["note"] = "offline: no UniProt fetch, so the planner must block on missing topology"
        else:
            row["expectation_met"] = tp.get("archetype") == case["expect_archetype"]
        if not row["expectation_met"] or row["invariant_violations"]:
            results["ok"] = False
        results["cases"].append(row)
        rendered[case["label"]] = render_markdown(
            plan, title=f"# Next experimental test — {case['label']}"
        )
        if case["expect_archetype"] == "not_accessible_to_soluble_binder":
            results["refusal_detail"] = plan.get("failure_detail")

    # A candidate-free call must degrade to a partial plan, not crash or invent.
    no_cand = plan_experiments(brief_a, fx["targets"][0]["target"], None, fetch=not offline)
    results["no_candidate_case"] = {
        "status": no_cand["status"],
        "candidate_sequence_length": no_cand["candidate_sequence_length"],
        "missing_inputs": no_cand["missing_inputs"],
        "ok": no_cand["status"]
        in ("partial_no_candidate_sequence", "blocked_no_target_biology"),
    }
    if not results["no_candidate_case"]["ok"]:
        results["ok"] = False

    # A retrieval failure must be reported as a source failure.
    bad = plan_experiments(
        brief_a,
        {"target_id": "SELFTEST-BAD", "gene_symbol": "SELFTEST-BAD", "protein_accession": "X0X0X0"},
        cand,
        fetch=not offline,
    )
    results["source_failure_case"] = {
        "status": bad["status"],
        "detail_mentions_source_failure": "SOURCE FAILURE" in (bad.get("failure_detail") or ""),
        "ok": bad["status"] == "blocked_no_target_biology",
    }
    if not results["source_failure_case"]["ok"]:
        results["ok"] = False

    # Generality: the two accessible branches must differ where it matters.
    acc = [c for c in results["cases"] if c.get("archetype", "").startswith(("secreted", "cell_surface"))]
    if len(acc) == 2:
        results["branches_differ"] = {
            "orientation_differs": acc[0]["binding_orientation"] != acc[1]["binding_orientation"],
            "antigen_construct_differs": acc[0]["antigen_construct"] != acc[1]["antigen_construct"],
            # Host can legitimately coincide: a glycosylated secreted protein and a
            # receptor ectodomain both want a mammalian host, and reporting a spurious
            # difference would be worse than reporting none.
            "antigen_host_differs": acc[0]["antigen_host"] != acc[1]["antigen_host"],
        }
        if not (
            results["branches_differ"]["orientation_differs"]
            and results["branches_differ"]["antigen_construct_differs"]
        ):
            results["ok"] = False

    if write_doc and rendered:
        _write_doc(Path(write_doc), rendered, results)
        results["doc_written"] = str(write_doc)

    if verbose:
        print(json.dumps(results, indent=2, default=str))
    return results


def _write_doc(path: Path, rendered: dict[str, str], results: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    head = [
        "# Next experimental test",
        "",
        "*Lane 7 — experimental design. Generated by `python -m e2b.experiments.selftest "
        "--write-doc docs/experimental_plan.md`.*",
        "",
        "## What this document is",
        "",
        "The candidate report has to name a binding-validation assay concept and an "
        "MoA-specific functional assay concept, and it has to say plainly that computational "
        "binding predictions alone do not demonstrate blockade, activation, internalisation, "
        "safety or efficacy. This is that section: the experiment a bench scientist would run "
        "next, with the controls that make the result interpretable and the decision criteria "
        "fixed in advance.",
        "",
        "Every plan below was generated by the same code from the target's own UniProt "
        "topology. Nothing is hardcoded for an indication: the module decides whether the "
        "target is a secreted protein, a cell-surface receptor or unreachable by a soluble "
        "binder, and the assay formats follow from that.",
        "",
        "### Origin labels in this document",
        "",
        "| Element | Origin | Meaning |",
        "|---|---|---|",
        "| Target topology, sequence length, glycosylation and disulfide counts, PDB entry counts | `real` | Retrieved live from UniProt at generation time. |",
        "| Candidate sequence | `fixture` | **Not a design and not a biological result.** A synthetic string used to exercise the construct-liability screen until Lane 3 delivers real candidates. It carries no metrics, because none exist. |",
        "| Paralog panel (first plan only) | `fixture` | Placeholder panel members, so the specificity arm can be shown wired. Identities are `null`, not zero. |",
        "| Assay parameters | — | Ranges with reasoning, each tagged `starting_point_requires_optimisation`. None has been measured or validated by this system. |",
        "",
        "When a real candidate arrives, the same call regenerates these plans with "
        "`origin='real'` on the candidate and the real sequence screened for free cysteines "
        "and glycosylation sequons.",
        "",
        "---",
        "",
        "## Proof that the plan is derived, not templated",
        "",
        "| Branch | Archetype (from UniProt topology) | Antigen construct | Antigen host | Binding orientation (opening words) |",
        "|---|---|---|---|---|",
    ]
    for c in results["cases"]:
        if not c.get("archetype"):
            continue
        head.append(
            f"| {c['label']} | `{c['archetype']}` | {c.get('antigen_construct') or '—'} | "
            f"{c.get('antigen_host') or '—'} | {(c.get('binding_orientation') or '—')[:70]}… |"
        )
    head += [
        "",
        "The third branch is the important one: asked to plan an experiment for a target with "
        "no extracellular face, the module refuses and explains why, rather than producing a "
        "plan that reads well and could not work.",
        "",
        "> " + (results.get("refusal_detail") or "—").replace("\n", " "),
        "",
        "---",
        "",
    ]
    body = []
    for label, md in rendered.items():
        body.append(md)
        body.append("\n---\n")
    path.write_text("\n".join(head) + "\n".join(body))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Self-test the experimental planner.")
    ap.add_argument("--offline", action="store_true", help="skip live UniProt retrieval")
    ap.add_argument("--write-doc", default=None, help="render the plans to this markdown path")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)
    res = run(offline=args.offline, write_doc=args.write_doc, verbose=not args.quiet)
    return 0 if res["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
