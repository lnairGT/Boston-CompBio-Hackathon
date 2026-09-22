"""Headless entrypoint: run one investigation from the command line.

    python -m e2b.cli "atopic dermatitis" \
        --desired-effect "suppress type 2 inflammatory signalling" \
        --moa antagonise

Exits non-zero only on an internal failure. An honest scientific dead end -- no Census
coverage, no binder-feasible target, insufficient evidence -- is a SUCCESSFUL run that
reports a limitation, and exits 0. Confusing the two would make the CLI unusable in a
pipeline, because "this indication has no accessible target" is a result, not a crash.
"""

from __future__ import annotations

import argparse
import json
import sys

from .contracts import Budgets
from .agent.loop import AgentLoop, start_run

OUTCOME_NOTES = {
    "completed_candidate_report": "Design candidates were generated and evaluated.",
    "evidence_report_only": "Evidence and a reviewed target proposal, but no design run completed.",
    "needs_clarification": "The indication was ambiguous enough to change which targets get retrieved.",
    "insufficient_evidence": "Retrieved evidence could not support any ranking.",
    "no_binder_feasible_target": "Every candidate was unreachable by the requested modality.",
    "source_unavailable": "A required source failed. This is NOT evidence of absence.",
    "budget_exhausted": "Bounds were reached; this is the best supported partial result.",
    "design_failed": "A design job ran and did not produce an acceptable candidate.",
    "embedding_unavailable": "No suitable Census release/embedding combination was found.",
}


def _make_llm(model: str | None):
    """Prefer the in-kernel bridge when running inside Claude Science, else the API."""
    try:  # pragma: no cover - environment-dependent
        host = __builtins__["host"] if isinstance(__builtins__, dict) else __builtins__.host  # type: ignore[attr-defined]
        from .agent.loop import HostLLMClient
        return HostLLMClient(host, model=model)
    except Exception:
        from .agent.loop import AnthropicClient
        return AnthropicClient(model=model or "claude-sonnet-4-5")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="e2b", description="Indication-agnostic evidence-to-binder investigation.")
    p.add_argument("indication", help="Free text. Any indication; it is resolved at runtime.")
    p.add_argument("--desired-effect", default=None, help="The biological outcome you want.")
    p.add_argument("--moa", default=None, help="Intervention direction, e.g. 'antagonise'.")
    p.add_argument("--modality", default="de novo minibinder")
    p.add_argument("--pool", type=int, default=20, help="Candidate pool size (bounded).")
    p.add_argument("--shortlist", type=int, default=3)
    p.add_argument("--max-decisions", type=int, default=30)
    p.add_argument("--allow-design", action="store_true",
                   help="Permit the ONE paid design batch. Omit for a dry evidence run.")
    p.add_argument("--model", default=None)
    p.add_argument("--runs-dir", default="runs")
    p.add_argument("--json", action="store_true", help="Print the manifest as JSON and nothing else.")
    args = p.parse_args(argv)

    budgets = Budgets(
        max_tool_decisions=args.max_decisions,
        max_candidate_pool=args.pool,
        max_paid_design_batches=1 if args.allow_design else 0,
    )
    loop, store = start_run(
        args.indication, llm=_make_llm(args.model), desired_effect=args.desired_effect,
        moa=args.moa, modality=args.modality, budgets=budgets, root=args.runs_dir)

    result = loop.run(verbose=True)

    if args.json:
        print(json.dumps(result["manifest"], indent=2, default=str))
        return 0

    m = result["manifest"]
    out = sys.stdout
    out.write(f"\nrun {m['run_id']}  ({m['tool_decisions_used']} tool decisions)\n")
    out.write(f"indication  {m['raw_indication']!r} -> {m['resolved_disease']['id']} "
              f"({m['resolved_disease']['label']})\n")
    scope = m.get("candidate_pool_scope") or {}
    out.write(f"pool        {scope.get('retrieved')} of {scope.get('total_associated_in_release')} "
              f"associated targets in {scope.get('release')} (bounded, not exhaustive)\n")
    out.write(f"screened    {m['structures_assessed']} protein(s); {len(m['excluded_targets'])} excluded\n")

    comp = store.comparison
    if comp:
        sym = {t: store.targets[t].gene_symbol for t in store.targets}
        out.write(f"rubric      {m['rubric_version']}\n")
        for i, tid in enumerate(comp.ordered_target_ids[: args.shortlist], 1):
            a = store.assessments[tid]
            ratings = " ".join(f"{c.criterion.split('_')[0]}={c.rating}" for c in a.criteria)
            out.write(f"  {i}. {sym.get(tid, tid):9s} {comp.scores[tid]:.3f}  {ratings}\n")
            out.write(f"     feasibility={a.design_feasibility.rating} (reported separately)\n")
        if comp.ties:
            out.write(f"  TIE: {[[sym.get(t, t) for t in g] for g in comp.ties]} "
                      f"-- not broken automatically\n")
        if comp.pending_screen:
            out.write(f"  {len(comp.pending_screen)} target(s) not accessibility-screened, "
                      f"therefore not eligible\n")
        for tid, why in list(comp.excluded.items())[:4]:
            out.write(f"  excluded {sym.get(tid, tid):9s} {why[:88]}\n")

    out.write(f"\noutcome     {result['outcome']}\n")
    out.write(f"            {OUTCOME_NOTES.get(result['outcome'], '')}\n")
    out.write(f"artifacts   {store.root}/\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
