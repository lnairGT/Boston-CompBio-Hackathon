"""Pure record -> view-dict transforms.

Everything the templates render passes through here, so the honesty rules live in one
testable place. No function in this module performs I/O, and none of them invents a
value: where a record is missing a field, the view dict carries an explicit
``None``/"not measured"/"no data" marker that the template renders as such.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any

from e2b.config import CRITERIA, RUBRIC_VERSION
from e2b.contracts import JobStatus
from e2b.ui.honesty import (
    count,
    freshness,
    metric_caption,
    modality_phrase,
    num,
    origin_badge,
    ts,
)
from e2b.ui.state import STAGE_LABELS, STAGE_ORDER, RunState

MAX_COMPARE_CARDS = 3

RATING_CSS = {
    "strong": "rating-strong",
    "moderate": "rating-moderate",
    "weak": "rating-weak",
    "unknown": "rating-unknown",
}

RATING_TEXT = {
    "strong": "strong",
    "moderate": "moderate",
    "weak": "weak",
    "unknown": "no evidence retrieved",
}

JOB_STATE_VIEW: dict[str, dict[str, str]] = {
    "draft": {"css": "job-idle", "text": "Draft — not submitted."},
    "awaiting_selection": {"css": "job-idle", "text": "Waiting for a target and site to be selected."},
    "ready": {"css": "job-idle", "text": "Ready to submit. Nothing has been spent yet."},
    "queued": {"css": "job-pending", "text": "Queued on the execution route. Not started."},
    "running": {"css": "job-pending", "text": "Running."},
    "succeeded": {"css": "job-ok", "text": "Completed."},
    "failed": {"css": "job-bad", "text": "Failed. No candidates from this run are usable."},
    "timed_out": {
        "css": "job-bad",
        "text": "Timed out. This is an execution failure, not a result — it says nothing "
                "about whether a binder is designable.",
    },
    "cancelled": {"css": "job-bad", "text": "Cancelled before completion."},
}


def _now() -> datetime:
    return datetime.now(timezone.utc)


def provenance_view(prov: Any) -> dict[str, Any] | None:
    if prov is None:
        return None
    return {
        "source": prov.source,
        "release": prov.source_release,
        "endpoint": prov.endpoint,
        "query": prov.query,
        "variables": prov.variables,
        "retrieved_at": ts(prov.retrieved_at),
        "elapsed": num(prov.elapsed_ms, precision=0, unit="ms", missing="timing not recorded"),
        "badge": origin_badge(prov.origin),
        "notes": prov.notes,
    }


# ---------------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------------


def stage_view(state: RunState) -> dict[str, Any]:
    """Current agent stage with a timestamp and an evidence-grounded 'next'.

    There is deliberately no percentage anywhere in this dict. The only quantities are
    counts the records actually carry: tool decisions used against the budget, and how
    many targets/evidence items exist.
    """
    run = state.agent_run
    if run is None:
        return {
            "known": False,
            "label": "No agent run attached",
            "stage_key": None,
            "updated": ts(None),
            "age": "",
            "next": "The spine has not started an agent run for this brief, or has not "
                    "registered it with the UI.",
            "budget_text": None,
            "questions": [],
            "last_action": None,
            "status": None,
            "outcome": None,
            "steps": [],
        }

    key = run.current_stage
    label = STAGE_LABELS.get(key, key)
    b = run.budgets
    last = run.tool_action_log[-1] if run.tool_action_log else None

    steps = []
    reached = True
    for s in STAGE_ORDER:
        if s == key:
            steps.append({"key": s, "label": STAGE_LABELS[s], "state": "current"})
            reached = False
        else:
            steps.append({"key": s, "label": STAGE_LABELS[s], "state": "done" if reached else "todo"})
    if key not in STAGE_ORDER:
        steps.append({"key": key, "label": label, "state": "current"})

    return {
        "known": True,
        "label": label,
        "stage_key": key,
        "status": run.status,
        "outcome": run.final_outcome,
        "outcome_detail": run.outcome_detail,
        "updated": ts(run.updated_at),
        "age": freshness(run)["age"],
        "started": ts(run.started_at),
        "next": _next_explanation(state, run),
        # A real count against a real budget. Rendered as text, never as a progress bar.
        "budget_text": (
            f"{b.tool_decisions_used} of {b.max_tool_decisions} tool decisions used; "
            f"{b.paid_design_batches_used} of {b.max_paid_design_batches} paid design batches used"
        ),
        "questions": list(run.pending_questions),
        "plan": run.plan_summary,
        "last_action": (
            {
                "step": last.step,
                "tool": last.tool_name,
                "reason": last.reason,
                "summary": last.result_summary,
                "error": last.error,
                "at": ts(last.started_at),
                "elapsed": num(last.elapsed_ms, precision=0, unit="ms", missing="timing not recorded"),
            }
            if last
            else None
        ),
        "steps": steps,
    }


def _next_explanation(state: RunState, run: Any) -> str:
    """Describe the next action from what the records actually contain."""
    if run.pending_questions:
        return (
            "Paused for your answer: "
            + run.pending_questions[0]
            + " Nothing downstream runs until this is resolved."
        )
    if run.status == "completed":
        return f"Run finished with outcome '{run.final_outcome}'. " + (run.outcome_detail or "")
    if run.status == "failed":
        return f"Run failed. {run.outcome_detail or 'No detail recorded.'}"
    if run.budgets.exhausted():
        return (
            f"Tool-decision budget is exhausted ({run.budgets.tool_decisions_used}/"
            f"{run.budgets.max_tool_decisions}). The agent will stop and report what it has."
        )
    if run.plan_summary:
        return run.plan_summary

    n_t, n_a, n_e = len(state.targets), len(state.assessments), len(state.evidence)
    key = run.current_stage
    basis = f"Basis so far: {n_t} target(s), {n_a} assessment(s), {n_e} evidence item(s) in hand."
    mapping = {
        "resolve_indication": "Next: resolve the free-text indication to an ontology term "
                              "and present the candidates rather than picking one silently.",
        "gather_targets": "Next: pull disease-associated targets from Open Targets and bound "
                          "the candidate pool.",
        "assess_targets": "Next: rate each target against the rubric criteria, recording what "
                          "was looked for and not found.",
        "cellular_context": "Next: summarise expression in mechanism-relevant populations "
                            "per donor from the Census.",
        "rank_targets": f"Next: apply the deterministic rubric ({RUBRIC_VERSION}) and expose "
                        "ties instead of forcing a winner.",
        "structure_and_site": "Next: check structure availability and propose a binding site "
                              "with an MoA rationale.",
        "design_request": "Next: assemble a DesignRequest; it cannot launch until it is "
                          "approved and complete.",
        "design_run": "Next: poll the design run for real status. No result is shown until "
                      "the executor reports one.",
        "report": "Next: assemble the candidate report with metric definitions and limitations.",
    }
    return f"{mapping.get(key, 'Stage ' + key + ' is in progress.')} {basis}"


# ---------------------------------------------------------------------------------
# Brief
# ---------------------------------------------------------------------------------


def brief_view(state: RunState) -> dict[str, Any]:
    brief = state.brief
    if brief is None:
        return {"exists": False, "resolution_state": "no_brief"}

    cands = [
        {
            "disease_id": c.disease_id,
            "label": c.label,
            "match_type": c.match_type,
            "score": num(c.score, precision=3, missing="no score returned"),
            "description": c.description,
            "selected": c.disease_id == brief.resolved_disease_id,
        }
        for c in brief.resolution_candidates
    ]

    if brief.resolved_disease_id and brief.status == "confirmed":
        res_state = "confirmed"
    elif brief.resolved_disease_id and len(cands) > 1:
        res_state = "ambiguous"
    elif brief.resolved_disease_id:
        res_state = "resolved_unconfirmed"
    elif cands:
        res_state = "ambiguous_unresolved"
    else:
        res_state = "unresolved"

    return {
        "exists": True,
        "brief_id": brief.brief_id,
        "raw_indication": brief.raw_indication,
        "desired_effect": brief.desired_effect,
        "moa": brief.moa,
        "modality_raw": brief.modality,
        "modality": modality_phrase(brief.modality),
        "species": brief.species,
        "species_taxon_id": brief.species_taxon_id,
        "resolved_disease_id": brief.resolved_disease_id,
        "resolved_disease_label": brief.resolved_disease_label,
        "resolution_rationale": brief.resolution_rationale,
        "resolution_candidates": cands,
        "resolution_state": res_state,
        "ambiguous": res_state in ("ambiguous", "ambiguous_unresolved"),
        "census_label": brief.census_disease_label,
        "census_ontology_id": brief.census_disease_ontology_id,
        "census_mapping": provenance_view(brief.census_mapping_provenance),
        "status": brief.status,
        "provenance": provenance_view(brief.provenance),
        "badge": origin_badge(brief.provenance.origin if brief.provenance else None),
    }


# ---------------------------------------------------------------------------------
# Compare
# ---------------------------------------------------------------------------------


def _criterion_view(rating: Any, evidence: dict[str, Any]) -> dict[str, Any]:
    items = [evidence[e] for e in rating.evidence_ids if e in evidence]
    return {
        "criterion": rating.criterion,
        "label": next(
            (c["label"] for c in CRITERIA if c["key"] == rating.criterion), rating.criterion
        ),
        "weight": next((c["weight"] for c in CRITERIA if c["key"] == rating.criterion), None),
        "rating": rating.rating,
        "rating_text": RATING_TEXT.get(rating.rating, rating.rating),
        "css": RATING_CSS.get(rating.rating, "rating-unknown"),
        "rationale": rating.rationale,
        "missingness": rating.missingness,
        "evidence_ids": list(rating.evidence_ids),
        "evidence_count": len(rating.evidence_ids),
        "resolved_evidence_count": len(items),
        "is_gap": rating.rating == "unknown",
    }


def _feasibility_view(f: Any, modality: str | None) -> dict[str, Any]:
    if f.accessible_to_modality is None:
        access_text = "unknown — no topology or location annotation retrieved"
        access_css = "rating-unknown"
    elif f.accessible_to_modality:
        access_text = f"accessible to a {modality_phrase(modality)}"
        access_css = "rating-strong"
    else:
        access_text = f"NOT accessible to a {modality_phrase(modality)}"
        access_css = "rating-bad"
    if f.structure_available is None:
        struct_text = "unknown — structure availability not checked"
    elif f.structure_available:
        struct_text = "structure available: " + (", ".join(f.structure_ids) or "unnamed")
    else:
        struct_text = "no usable structure found"
    return {
        "rating": f.rating,
        "rating_text": RATING_TEXT.get(f.rating, f.rating),
        "css": RATING_CSS.get(f.rating, "rating-unknown"),
        "accessible": f.accessible_to_modality,
        "access_text": access_text,
        "access_css": access_css,
        "localization_basis": list(f.localization_basis),
        "structure_available": f.structure_available,
        "structure_text": struct_text,
        "structure_ids": list(f.structure_ids),
        "blocking_reasons": list(f.blocking_reasons),
        "notes": f.notes,
    }


def compare_view(state: RunState) -> dict[str, Any]:
    comp = state.comparison
    modality = state.brief.modality if state.brief else None

    if comp is not None and comp.ordered_target_ids:
        order = comp.ordered_target_ids
    else:
        order = sorted(state.assessments.keys())

    tie_of: dict[str, list[str]] = {}
    if comp is not None:
        for group in comp.ties:
            for tid in group:
                tie_of[tid] = [g for g in group if g != tid]

    cards: list[dict[str, Any]] = []
    for rank, tid in enumerate(order[:MAX_COMPARE_CARDS], start=1):
        target = state.targets.get(tid)
        a = state.assessments.get(tid)
        if target is None and a is None:
            continue
        crits = [_criterion_view(c, state.evidence) for c in (a.criteria if a else [])]
        score = comp.scores.get(tid) if comp else None
        cards.append(
            {
                "rank": rank,
                "target_id": tid,
                "gene_symbol": target.gene_symbol if target else tid,
                "approved_name": target.approved_name if target else None,
                "protein_accession": target.protein_accession if target else None,
                "biotype": target.biotype if target else None,
                "badge": origin_badge(
                    target.provenance.origin if (target and target.provenance) else None
                ),
                "provenance": provenance_view(target.provenance if target else None),
                "score": num(score, precision=3, missing="not scored"),
                "criteria": crits,
                "gaps": [c for c in crits if c["is_gap"]],
                "liabilities": list(a.liabilities) if a else [],
                "direction_compatibility": a.direction_compatibility if a else "unknown",
                "disease_rationale": a.disease_rationale if a else None,
                "moa_rationale": a.moa_rationale if a else None,
                "review_status": a.review_status if a else None,
                "feasibility": _feasibility_view(a.design_feasibility, modality) if a else None,
                "tied_with": tie_of.get(tid, []),
                "selected": tid == state.selected_target_id,
                "unresolved_criteria": (comp.unresolved_criteria.get(tid, []) if comp else []),
            }
        )

    excluded = []
    if comp is not None:
        for tid, reason in comp.excluded.items():
            t = state.targets.get(tid)
            excluded.append(
                {"target_id": tid, "gene_symbol": t.gene_symbol if t else tid, "reason": reason}
            )
    for tid, a in state.assessments.items():
        if a.excluded and not any(e["target_id"] == tid for e in excluded):
            t = state.targets.get(tid)
            excluded.append(
                {
                    "target_id": tid,
                    "gene_symbol": t.gene_symbol if t else tid,
                    "reason": a.exclusion_reason or "no reason recorded",
                }
            )

    return {
        "exists": bool(cards) or bool(excluded) or comp is not None,
        "cards": cards,
        "excluded": excluded,
        "ties": [g for g in (comp.ties if comp else [])],
        "insufficient_evidence": bool(comp.insufficient_evidence) if comp else False,
        "explanation": comp.explanation if comp else None,
        "rubric_version": comp.rubric_version if comp else RUBRIC_VERSION,
        "pool_scope": comp.candidate_pool_scope if comp else {},
        "n_assessed": len(state.assessments),
        "n_shown": len(cards),
        "score_caveat": (
            "The rubric score is a weighted sum of ordinal evidence ratings under "
            f"{comp.rubric_version if comp else RUBRIC_VERSION}. It is a bookkeeping device "
            "for comparing what was retrieved — not a probability of clinical success, and "
            "not a measured quantity."
        ),
    }


# ---------------------------------------------------------------------------------
# Inspect
# ---------------------------------------------------------------------------------


def _evidence_view(e: Any) -> dict[str, Any]:
    return {
        "evidence_id": e.evidence_id,
        "subject_id": e.subject_id,
        "predicate": e.predicate,
        "object_id": e.object_id,
        "object_value": e.object_value,
        "evidence_type": e.evidence_type,
        "direction": e.direction,
        "source": e.source,
        "source_url": e.source_url,
        "doi": e.doi,
        "locator": e.source_locator,
        "excerpt": e.supporting_excerpt,
        "measurement": e.measurement,
        "disease_scope": e.disease_scope,
        "scope_disease_id": e.scope_disease_id,
        "scope_warning": (
            "Propagated from a related disease term, not observed for the resolved disease "
            "itself." if e.disease_scope == "indirect_propagated" else None
        ),
        "limitations": list(e.limitations),
        "retrieved_at": ts(e.retrieved_at),
        "badge": origin_badge(e.origin),
        "triple": f"{e.subject_id} —{e.predicate}→ "
                  f"{e.object_id if e.object_id is not None else e.object_value}",
    }


def _cell_summary_view(cs: Any) -> dict[str, Any]:
    per_donor = [
        {"donor": d, "value": num(v, precision=3, missing="not measurable in this donor"),
         "cells": count(cs.donor_counts.get(d))}
        for d, v in sorted(cs.per_donor_values.items())
    ] or [
        {"donor": k, "value": num(v, precision=3), "cells": count(cs.donor_counts.get(k))}
        for k, v in sorted(cs.values.items())
    ]
    return {
        "summary_id": cs.summary_id,
        "target_id": cs.target_id,
        "census_release": cs.census_release,
        "dataset_ids": list(cs.dataset_ids),
        "tissue": cs.tissue,
        "condition": cs.condition,
        "population": cs.canonical_population_label or cs.canonical_population_id
                      or "population not labelled",
        "original_labels": list(cs.original_labels),
        "n_donors": count(cs.n_donors, missing="donor count not recorded"),
        "cell_counts": count(cs.cell_counts, missing="cell count not recorded"),
        "metric": cs.expression_metric,
        "method": cs.expression_method,
        "per_donor": per_donor,
        "primary_data_handling": cs.primary_data_handling,
        "limitations": list(cs.limitations),
        "query_filter": cs.query_filter,
        "badge": origin_badge(cs.origin),
        "provenance": provenance_view(cs.provenance),
    }


def _neighborhood_view(n: Any) -> dict[str, Any]:
    return {
        "evidence_id": n.evidence_id,
        "query_population": n.query_population,
        "reference_universe": n.reference_universe,
        "metric": n.distance_metric,
        "preprocessing": n.preprocessing,
        "k": n.k,
        "n_neighbors": len(n.neighbor_ids),
        "distance_range": (
            f"{min(n.distances):.4f} – {max(n.distances):.4f}" if n.distances else "no distances"
        ),
        "filters": n.biological_filters,
        "caveats": list(n.caveats),
    }


def inspect_view(state: RunState) -> dict[str, Any]:
    ev = [_evidence_view(e) for e in state.evidence.values()]
    ev.sort(key=lambda d: (d["subject_id"], d["evidence_type"], d["evidence_id"]))

    by_target: dict[str, list[dict[str, Any]]] = {}
    for item in ev:
        by_target.setdefault(item["subject_id"], []).append(item)

    gene_of = {tid: t.gene_symbol for tid, t in state.targets.items()}
    groups = [
        {
            "target_id": tid,
            "gene_symbol": gene_of.get(tid, tid),
            "records": items,
            "n": len(items),
            "sources": sorted({i["source"] for i in items}),
        }
        for tid, items in sorted(by_target.items(), key=lambda kv: gene_of.get(kv[0], kv[0]))
    ]

    return {
        "exists": bool(ev) or bool(state.cell_summaries) or bool(state.source_failures),
        "groups": groups,
        "n_evidence": len(ev),
        "cell_summaries": [_cell_summary_view(c) for c in state.cell_summaries],
        "neighborhoods": [_neighborhood_view(n) for n in state.neighborhoods],
        "failures": [f.view() for f in state.source_failures],
        "n_failures": len(state.source_failures),
        "graph_edges": [
            {"triple": i["triple"], "direction": i["direction"], "type": i["evidence_type"],
             "evidence_id": i["evidence_id"], "badge": i["badge"]}
            for i in ev
        ],
    }


# ---------------------------------------------------------------------------------
# Design
# ---------------------------------------------------------------------------------


def job_view(run: Any) -> dict[str, Any]:
    """Real job status. No synthetic percentage — stage name and timestamps only."""
    if run is None:
        return {
            "exists": False,
            "state": None,
            "css": "job-idle",
            "text": "No design run has been submitted for this brief.",
        }
    status = run.status.value if isinstance(run.status, JobStatus) else str(run.status)
    meta = JOB_STATE_VIEW.get(status, {"css": "job-idle", "text": f"State '{status}'."})
    fresh = freshness(run)
    return {
        "exists": True,
        "run_id": run.run_id,
        "design_id": run.design_id,
        "state": status,
        "css": meta["css"],
        "text": meta["text"],
        "stage": run.stage or "stage not reported by the executor",
        "route": run.execution_route,
        "handle": run.execution_handle,
        "started": ts(run.started_at, missing="not started"),
        "updated": ts(run.updated_at, missing="no update recorded"),
        "finished": ts(run.finished_at, missing="not finished"),
        "age": fresh["age"],
        "badge": fresh["badge"],
        "replay": fresh["replay"],
        "replay_note": fresh["note"],
        "error": run.error,
        "model": run.model_name,
        "model_version": run.model_version,
        "parameters": run.parameters,
        "seed": run.seed,
        "git_commit": run.git_commit,
        "cost": num(run.cost_usd_estimate, precision=2, unit="USD (estimate)",
                    missing="cost not recorded"),
        "artifact_refs": list(run.artifact_refs),
        "is_terminal": status in ("succeeded", "failed", "timed_out", "cancelled"),
        "is_pending": status in ("queued", "running"),
    }


def design_view(state: RunState) -> dict[str, Any]:
    req = state.design_request
    target = state.targets.get(state.selected_target_id) if state.selected_target_id else None
    modality = state.brief.modality if state.brief else None

    problems: list[str] = []
    if req is None:
        problems.append("No DesignRequest has been prepared by the spine for this run.")
    else:
        try:
            req.assert_launchable()
        except ValueError as exc:
            problems = [p.strip() for p in str(exc).split(":", 1)[-1].split(";")]

    blocked_by_feasibility: list[str] = []
    a = state.assessments.get(state.selected_target_id) if state.selected_target_id else None
    if a is not None:
        f = a.design_feasibility
        if f.accessible_to_modality is False:
            blocked_by_feasibility.append(
                f"{target.gene_symbol if target else state.selected_target_id} is annotated as not "
                f"accessible to a soluble binder."
            )
        blocked_by_feasibility.extend(f.blocking_reasons)

    ticket = None
    for t in sorted(state.launch_tickets.values(), key=lambda x: x.requested_at, reverse=True):
        if t.outcome != "pending":
            ticket = t
            break

    region = req.binding_region if req else None
    can_launch = (
        req is not None
        and not problems
        and not blocked_by_feasibility
        and state.inputs_confirmed
        and state.selected_target_id is not None
    )

    return {
        "exists": req is not None or state.selected_target_id is not None,
        "selected_target_id": state.selected_target_id,
        "gene_symbol": target.gene_symbol if target else None,
        "protein_accession": target.protein_accession if target else None,
        "modality": modality_phrase(modality),
        "design_id": req.design_id if req else None,
        "review_status": req.review_status if req else None,
        "structure_id": req.structure_id if req else None,
        "structure_version": req.structure_version if req else None,
        "sequence_checksum": req.protein_sequence_checksum if req else None,
        "n_designs": req.n_designs if req else None,
        "model_name": req.model_name if req else None,
        "model_config": req.model_config_ if req else {},
        "moa_rationale": req.moa_rationale if req else None,
        "evidence_ids": list(req.evidence_ids) if req else [],
        "input_hash": (req.input_hash or req.compute_input_hash()) if req else None,
        "region": (
            {
                "description": region.description,
                "chain": region.chain,
                "residues": region.residues,
                "hotspots": region.hotspot_residues,
                "numbering": region.numbering_scheme,
                "rationale": region.rationale,
            }
            if region
            else None
        ),
        "problems": [p for p in problems if p],
        "blocked_by_feasibility": blocked_by_feasibility,
        "inputs_confirmed": state.inputs_confirmed,
        "can_launch": can_launch,
        "launch_reason": (
            None if can_launch
            else "Confirm the inputs to enable the launch control."
            if (req is not None and not problems and not blocked_by_feasibility)
            else "The request is not launchable yet; see the blocking items above."
        ),
        "ticket": (
            {"outcome": ticket.outcome, "message": ticket.message,
             "at": ts(ticket.requested_at), "run_id": ticket.run_id}
            if ticket else None
        ),
        "job": job_view(state.design_run),
    }


# ---------------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------------

DEFAULT_LIMITATIONS = [
    "Every number on this page that came from a design or folding model is a computational "
    "score. None of them is a measured binding affinity.",
    "No experimental validation has been performed.",
]


def _candidate_view(c: Any, modality: str | None) -> dict[str, Any]:
    seq = c.sequence
    return {
        "candidate_id": c.candidate_id,
        "run_id": c.run_id,
        "design_id": c.design_id,
        "target_id": c.target_id,
        "sequence": seq,
        "sequence_present": seq is not None,
        "sequence_note": None if seq else "No sequence was returned for this candidate.",
        "length": count(c.sequence_length, missing="no sequence"),
        "structure_ref": c.structure_ref,
        "structure_note": None if c.structure_ref else "No structure file was returned.",
        "metrics": [metric_caption(m) for m in c.metrics],
        "n_metrics": len(c.metrics),
        "evaluation_status": c.evaluation_status,
        "limitations": list(c.limitations),
        "proposed_validation": list(c.proposed_validation),
        "badge": origin_badge(c.origin),
        "modality": modality_phrase(modality),
    }


def results_view(state: RunState) -> dict[str, Any]:
    modality = state.brief.modality if state.brief else None
    job = job_view(state.design_run)
    cands = [_candidate_view(c, modality) for c in state.candidates]

    if cands:
        status = "have_results"
    elif job["exists"] and job.get("is_pending"):
        status = "pending"
    elif job["exists"] and job.get("state") in ("failed", "timed_out", "cancelled"):
        status = job["state"]
    elif job["exists"] and job.get("state") == "succeeded":
        status = "succeeded_no_candidates"
    else:
        status = "no_data"

    metric_names = sorted({m["name"] for c in cands for m in c["metrics"]})
    definitions = []
    for name in metric_names:
        example = next(m for c in cands for m in c["metrics"] if m["name"] == name)
        definitions.append(
            {
                "name": name,
                "method": example["method"],
                "model": example["model"],
                "direction": example["direction"],
                "interpretation": example["interpretation"],
                "disclaimer": example["disclaimer"],
            }
        )

    limitations = list(DEFAULT_LIMITATIONS)
    for c in cands:
        for lim in c["limitations"]:
            if lim not in limitations:
                limitations.append(lim)
    if any(c["badge"]["label"] == "fixture" for c in cands):
        limitations.insert(
            0,
            "One or more candidates below are FIXTURES. They are interface stand-ins, not "
            "design output, and must not be reported as results.",
        )

    next_tests: list[str] = []
    for c in cands:
        for t in c["proposed_validation"]:
            if t not in next_tests:
                next_tests.append(t)

    return {
        "status": status,
        "candidates": cands,
        "n_candidates": len(cands),
        "job": job,
        "metric_definitions": definitions,
        "limitations": limitations,
        "next_tests": next_tests,
        "next_tests_note": (
            None if next_tests
            else "No experimental next step is recorded on these candidates."
        ),
        "modality": modality_phrase(modality),
    }


# ---------------------------------------------------------------------------------
# Whole page
# ---------------------------------------------------------------------------------


def page_view(state: RunState) -> dict[str, Any]:
    origins = sorted(state.origins_present)
    return {
        "run_key": state.run_key,
        "updated": ts(state.updated_at),
        "scenario_label": state.scenario_label,
        "notice": state.notice,
        "origins_present": [origin_badge(o) for o in origins],
        "is_fixture_run": state.is_fixture_run,
        "stage": stage_view(state),
        "brief": brief_view(state),
        "compare": compare_view(state),
        "inspect": inspect_view(state),
        "design": design_view(state),
        "results": results_view(state),
        "rubric_version": RUBRIC_VERSION,
        "now": ts(_now()),
    }


__all__ = [
    "brief_view",
    "compare_view",
    "design_view",
    "inspect_view",
    "job_view",
    "page_view",
    "results_view",
    "stage_view",
]
