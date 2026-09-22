"""Typed application tools exposed to the Claude tool-use loop.

The division the spec insists on: **Claude chooses actions; deterministic services
calculate metrics, enforce budgets and validate scientific preconditions.** So every
handler here is ordinary Python. None of them asks a model for a number. The model's job is
to decide *which* tool to call next and to explain the result using retrieved evidence IDs.

Capability degradation is deliberate. Lanes 2 (Census) and 3 (design) land their adapters
independently, and until they do, their tools return
``{"status": "capability_unavailable", ...}`` rather than raising. That is an honest tool
result the agent can branch on -- it maps onto the ``source_unavailable`` and
``embedding_unavailable`` outcomes -- and it is emphatically **not** a fixture: no
biological value is invented, and the agent is told plainly that the capability is absent.

Tool results are **compact summaries plus artifact references**. Full expression matrices
and embedding vectors never enter the model's context; they stay on disk and are referenced
by ID.
"""

from __future__ import annotations

import importlib
import json
import time
from dataclasses import dataclass, field
from typing import Any, Callable

from ..adapters.open_targets import OpenTargetsAdapter, OpenTargetsError
from ..adapters import structures as struct
from ..config import CENSUS_BOUNDS, CRITERIA, DEFAULT_BUDGETS, RUBRIC_VERSION
from ..contracts import ToolAction, utcnow
from ..pipeline import ranking as rk

TOOLS_VERSION = "tools/1.0.0"


class PreconditionFailed(RuntimeError):
    """A scientific precondition was not satisfied. Not retryable by the agent."""


class BudgetExhausted(RuntimeError):
    pass


@dataclass
class ToolSpec:
    name: str
    description: str
    input_schema: dict[str, Any]
    handler: Callable[..., dict[str, Any]]
    version: str = TOOLS_VERSION
    is_paid: bool = False
    requires_approval: bool = False

    def anthropic_schema(self) -> dict[str, Any]:
        return {"name": self.name, "description": self.description, "input_schema": self.input_schema}


def _optional(module_path: str) -> Any | None:
    """Import a lane's adapter if it has landed, else None."""
    try:
        return importlib.import_module(module_path, package="e2b")
    except Exception:
        return None


class ToolRegistry:
    """Holds the tools, the run state they mutate, and the budget they spend."""

    def __init__(self, *, store: "RunStore", budgets: dict[str, int] | None = None):
        self.store = store
        self.budgets = dict(DEFAULT_BUDGETS if budgets is None else budgets)
        self.ot = OpenTargetsAdapter()
        self.census = _optional("..adapters.cellxgene")
        self.design = _optional("..adapters.protein_design")
        self._specs: dict[str, ToolSpec] = {}
        self._register_all()

    # -- capability reporting -------------------------------------------------------

    def capabilities(self) -> dict[str, Any]:
        return {
            "open_targets": {"available": True, "release": self.ot.release()},
            "structures_uniprot_alphafold": {"available": True},
            "cellxgene_census": {
                "available": self.census is not None,
                "detail": None if self.census else "adapters/cellxgene.py has not landed (lane 2).",
            },
            "protein_design": {
                "available": self.design is not None,
                "detail": None if self.design else "adapters/protein_design.py has not landed (lane 3).",
            },
            "tools_registered": sorted(self._specs),
        }

    def anthropic_tools(self) -> list[dict[str, Any]]:
        return [s.anthropic_schema() for s in self._specs.values()]

    # -- dispatch -------------------------------------------------------------------

    def call(self, name: str, arguments: dict[str, Any], *, reason: str, step: int) -> tuple[dict[str, Any], ToolAction]:
        spec = self._specs.get(name)
        t0 = time.time()
        if spec is None:
            action = ToolAction(step=step, tool_name=name, arguments=arguments, reason=reason,
                                error=f"unknown tool '{name}'")
            return {"status": "error", "error": f"Unknown tool '{name}'."}, action

        if spec.is_paid:
            used = self.store.run.budgets.paid_design_batches_used
            cap = self.store.run.budgets.max_paid_design_batches
            if used >= cap:
                action = ToolAction(step=step, tool_name=name, tool_version=spec.version,
                                    arguments=arguments, reason=reason,
                                    error=f"paid batch budget exhausted ({used}/{cap})")
                return {
                    "status": "budget_exhausted",
                    "error": f"The paid design budget is spent ({used}/{cap} batches). "
                             "A further run requires explicit re-authorisation.",
                }, action

        # Models occasionally attach a parameter that is not in the declared schema. Dropping
        # the extras is better than a TypeError that burns a decision on a self-correcting
        # retry; the discarded keys are recorded on the action so it stays visible.
        declared = set((spec.input_schema.get("properties") or {}).keys())
        accepted = {k: v for k, v in arguments.items() if k in declared}
        discarded = sorted(set(arguments) - declared)

        retries = self.budgets.get("max_transient_retries", 2)
        attempt, last_err = 0, None
        while attempt <= retries:
            try:
                result = spec.handler(**accepted)
                if discarded:
                    result = dict(result)
                    result["_ignored_arguments"] = discarded
                action = ToolAction(
                    step=step, tool_name=name, tool_version=spec.version,
                    arguments=_redact(arguments), reason=reason,
                    result_summary=_summarise(result),
                    result_refs=list(result.get("artifact_refs", []) or []),
                    elapsed_ms=int((time.time() - t0) * 1000),
                )
                return result, action
            except (OpenTargetsError, struct.StructureAdapterError) as exc:
                if getattr(exc, "transient", False) and attempt < retries:
                    attempt += 1
                    time.sleep(1.5 * attempt)
                    last_err = exc
                    continue
                # A source failure is a source failure -- never evidence of absence.
                action = ToolAction(step=step, tool_name=name, tool_version=spec.version,
                                    arguments=_redact(arguments), reason=reason,
                                    error=f"{type(exc).__name__}: {exc}",
                                    elapsed_ms=int((time.time() - t0) * 1000))
                return {
                    "status": "source_unavailable",
                    "error": str(exc),
                    "retries_attempted": attempt,
                    "interpretation": "This is a SOURCE FAILURE, not evidence that no results exist. "
                                      "Do not conclude absence from it.",
                }, action
            except PreconditionFailed as exc:
                action = ToolAction(step=step, tool_name=name, tool_version=spec.version,
                                    arguments=_redact(arguments), reason=reason,
                                    error=f"precondition: {exc}",
                                    elapsed_ms=int((time.time() - t0) * 1000))
                return {"status": "precondition_failed", "error": str(exc)}, action
            except Exception as exc:  # noqa: BLE001 - surfaced to the agent, not swallowed
                action = ToolAction(step=step, tool_name=name, tool_version=spec.version,
                                    arguments=_redact(arguments), reason=reason,
                                    error=f"{type(exc).__name__}: {exc}",
                                    elapsed_ms=int((time.time() - t0) * 1000))
                return {"status": "error", "error": f"{type(exc).__name__}: {exc}"}, action

        return {"status": "source_unavailable", "error": str(last_err)}, ToolAction(
            step=step, tool_name=name, reason=reason, error=str(last_err))

    # -- registration ---------------------------------------------------------------

    def _add(self, spec: ToolSpec) -> None:
        self._specs[spec.name] = spec

    def _register_all(self) -> None:
        self._add(ToolSpec(
            name="resolve_indication",
            description=(
                "Resolve free-text disease input to ranked ontology matches with Open Targets IDs. "
                "Distinguishes exact, synonym and related matches. ALWAYS call this first; never "
                "assume an ontology ID string is valid."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "raw_text": {"type": "string", "description": "The user's indication text, verbatim."},
                    "size": {"type": "integer", "default": 8},
                },
                "required": ["raw_text"],
            },
            handler=self._resolve_indication,
        ))

        self._add(ToolSpec(
            name="discover_targets",
            description=(
                "Retrieve a BOUNDED pool of targets associated with a resolved disease ID, with "
                "per-datatype association score breakdown. Records retrieval order, page and limit. "
                "This is not exhaustive discovery."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "disease_id": {"type": "string"},
                    "limit": {"type": "integer", "default": DEFAULT_BUDGETS["max_candidate_pool"]},
                    "page_index": {"type": "integer", "default": 0},
                    "enable_indirect": {"type": "boolean", "default": True,
                                        "description": "Include ontology-propagated evidence."},
                },
                "required": ["disease_id"],
            },
            handler=self._discover_targets,
        ))

        self._add(ToolSpec(
            name="assess_structure",
            description=(
                "For a UniProt accession: verify protein identity, decide soluble-binder "
                "accessibility from TOPOLOGY (signal peptide / transmembrane / extracellular "
                "topological domain), and list experimental and predicted structures. Use this to "
                "screen modality suitability BEFORE investing in deep evidence retrieval."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "accession": {"type": "string"},
                    "mechanism_constraint": {"type": "string"},
                },
                "required": ["accession"],
            },
            handler=self._assess_structure,
        ))

        self._add(ToolSpec(
            name="get_target_evidence",
            description=(
                "Deep evidence for one target-disease pair: underlying evidence rows, known drugs "
                "with mechanism and clinical stage, tractability, safety liabilities. Preserves "
                "direct versus ontology-propagated scope per item."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "target_id": {"type": "string", "description": "Ensembl gene ID."},
                    "disease_id": {"type": "string"},
                    "evidence_size": {"type": "integer", "default": 50},
                },
                "required": ["target_id", "disease_id"],
            },
            handler=self._get_target_evidence,
        ))

        self._add(ToolSpec(
            name="compare_targets",
            description=(
                f"Run the DETERMINISTIC versioned rubric ({RUBRIC_VERSION}) over the assessments "
                "built so far. Returns an ordered shortlist, explicit ties, exclusions with reasons, "
                "and unresolved criteria. You explain this result; you do not invent weights. "
                f"Criteria: {[c['key'] for c in CRITERIA]}."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "shortlist": {"type": "integer", "default": DEFAULT_BUDGETS["depth_shortlist"]},
                },
            },
            handler=self._compare_targets,
        ))

        self._add(ToolSpec(
            name="summarize_cell_context",
            description=(
                "CELLxGENE Census donor-aware expression summary for genes in selected populations. "
                "Returns per-donor summaries before condition-level ones, with documented "
                "percent-positive definition and coverage gaps."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "disease_label": {"type": "string"},
                    "gene_symbols": {"type": "array", "items": {"type": "string"}},
                    "max_cells": {"type": "integer", "default": CENSUS_BOUNDS["max_cells_interactive"]},
                    "tissue": {"type": "string"},
                },
                "required": ["disease_label", "gene_symbols"],
            },
            handler=self._summarize_cell_context,
        ))

        self._add(ToolSpec(
            name="discover_census_embeddings",
            description="List precomputed cell embeddings available for a pinned Census release and organism.",
            input_schema={
                "type": "object",
                "properties": {
                    "census_release": {"type": "string"},
                    "organism": {"type": "string", "default": "homo_sapiens"},
                },
            },
            handler=self._discover_census_embeddings,
        ))

        self._add(ToolSpec(
            name="explore_cell_neighborhoods",
            description=(
                "Bounded nearest-neighbour search over an already-retrieved embedding slice. States "
                "the ACTUAL search universe. Similarity is context, never proof of identity or causality."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "embedding_slice_id": {"type": "string"},
                    "query_population": {"type": "string"},
                    "k": {"type": "integer", "default": CENSUS_BOUNDS["neighborhood_k"]},
                    "metric": {"type": "string", "default": "cosine"},
                },
                "required": ["embedding_slice_id", "query_population"],
            },
            handler=self._explore_neighborhoods,
        ))

        self._add(ToolSpec(
            name="submit_design",
            description=(
                "PAID. Submit a reviewed design request to the real protein-design route. Requires an "
                "APPROVED DesignRequest with a structure, binding region with hotspot residues, an MoA "
                "rationale and linked evidence IDs. Asynchronous: returns a run_id to poll. Identical "
                "requests are deduplicated by input hash and will not start a second paid job."
            ),
            input_schema={
                "type": "object",
                "properties": {"design_id": {"type": "string"}},
                "required": ["design_id"],
            },
            handler=self._submit_design,
            is_paid=True,
            requires_approval=True,
        ))

        self._add(ToolSpec(
            name="get_run",
            description="Poll a design run: status, stage, error and artifact references. Cheap; not budget-counted.",
            input_schema={"type": "object", "properties": {"run_id": {"type": "string"}}, "required": ["run_id"]},
            handler=self._get_run,
        ))

        self._add(ToolSpec(
            name="evaluate_candidates",
            description=(
                "Model-specific metrics for a completed run, with method and model version per metric. "
                "A succeeded compute job may legitimately have zero candidates passing evaluation."
            ),
            input_schema={"type": "object", "properties": {"run_id": {"type": "string"}}, "required": ["run_id"]},
            handler=self._evaluate_candidates,
        ))

    # -- handlers: implemented ------------------------------------------------------

    def _resolve_indication(self, raw_text: str, size: int = 8) -> dict[str, Any]:
        matches, prov = self.ot.resolve_indication(raw_text, size=size)
        if not matches:
            return {
                "status": "no_match",
                "raw_text": raw_text,
                "interpretation": "No disease in the Open Targets ontology matched this text. "
                                  "Ask the user to rephrase rather than guessing a nearby disease.",
            }
        self.store.record_resolution(raw_text, matches, prov)
        return {
            "status": "ok",
            "release": prov.source_release,
            "needs_clarification": self.ot.needs_clarification(matches),
            "matches": [
                {"disease_id": m.disease_id, "label": m.label, "match_type": m.match_type, "score": m.score}
                for m in matches
            ],
            "guidance": "Ask for clarification ONLY when plausible alternatives would materially "
                        "change which targets get retrieved.",
        }

    def _discover_targets(self, disease_id: str, limit: int = 20, page_index: int = 0,
                          enable_indirect: bool = True) -> dict[str, Any]:
        limit = min(limit, self.budgets.get("max_candidate_pool", 20))
        tgts, meta, prov = self.ot.discover_targets(
            disease_id, limit=limit, page_index=page_index, enable_indirect=enable_indirect)
        self.store.record_pool(tgts, meta, prov)
        return {
            "status": "ok",
            "scope": meta["scope"],
            "targets": [
                {
                    "target_id": t.target_id,
                    "gene_symbol": t.gene_symbol,
                    "protein_accession": t.protein_accession,
                    "overall_score": meta["associations"][t.target_id]["overall_score"],
                    "rank": meta["associations"][t.target_id]["rank"],
                    "datatype_scores": meta["associations"][t.target_id]["datatype_scores"],
                    "subcellular_locations": t.subcellular_locations[:3],
                    "tractability_modalities": sorted(t.tractability),
                }
                for t in tgts
            ],
        }

    def _assess_structure(self, accession: str, mechanism_constraint: str | None = None) -> dict[str, Any]:
        assessment, evs, prov = struct.assess_structure(accession, mechanism_constraint)
        self.store.record_structure(accession, assessment, evs)
        acc = assessment["accessibility"]
        return {
            "status": "ok",
            "accession": accession,
            "protein_name": assessment["protein_name"],
            "length": assessment["length"],
            "accessibility_verdict": acc["verdict"],
            "accessible": acc["accessible"],
            "basis": acc["basis"],
            "engageable_regions": acc["engageable_regions"],
            "n_experimental_structures": assessment["n_experimental_structures"],
            "top_structures": assessment["experimental_structures"][:4],
            "alphafold_available": assessment["alphafold"] is not None,
            "design_ready": assessment["design_ready"],
            "limitations": assessment["limitations"],
        }

    def _get_target_evidence(self, target_id: str, disease_id: str, evidence_size: int = 50) -> dict[str, Any]:
        ann, evs, prov = self.ot.get_target_evidence(target_id, disease_id, evidence_size=evidence_size)
        self.store.record_evidence(target_id, ann, evs, prov)
        direct = sum(1 for e in evs if e.disease_scope == "direct")
        return {
            "status": "ok",
            "target_id": target_id,
            "gene_symbol": ann["gene_symbol"],
            "evidence_retrieved": len(evs),
            "evidence_in_release": ann["evidence_count_in_release"],
            "direct_vs_propagated": {"direct": direct, "indirect_propagated": len(evs) - direct},
            "datatypes": sorted({e.evidence_type for e in evs}),
            "known_drugs": [
                {
                    "name": d["drug_name"], "type": d["drug_type"],
                    "mechanisms": [m["action_type"] for m in d["mechanisms"]],
                    "max_stage": d["max_clinical_stage"], "indications": d["indications"][:3],
                }
                for d in (ann["known_drugs"] or [])[:8]
            ],
            "safety_liabilities": [s["event"] for s in ann["safety_liabilities"]][:8],
            "tractability": ann["tractability"],
            "caveat": "Association and trial presence do not establish that the requested direction "
                      "of modulation is beneficial, nor efficacy for this indication.",
        }

    def _compare_targets(self, shortlist: int = 3) -> dict[str, Any]:
        result = self.store.run_comparison(shortlist=shortlist)
        symbols = {tid: t.gene_symbol for tid, t in self.store.targets.items()}
        return {
            "status": "ok",
            "rubric_version": result.rubric_version,
            "ordered": [f"{symbols.get(t, t)}({t})" for t in result.ordered_target_ids],
            "scores": {symbols.get(t, t): v for t, v in result.scores.items()},
            "ties": [[symbols.get(t, t) for t in g] for g in result.ties],
            "excluded": {symbols.get(t, t): v for t, v in result.excluded.items()},
            "pending_screen": {symbols.get(t, t): v for t, v in result.pending_screen.items()},
            "unresolved_criteria": {symbols.get(t, t): v for t, v in result.unresolved_criteria.items()},
            "insufficient_evidence": result.insufficient_evidence,
            "explanation": result.explanation,
            "instruction": (
                "Explain this ordering using the criteria ratings and evidence IDs. If a tie is "
                "present, say so and identify which missing evidence would break it -- do not invent "
                "a tiebreaker. Targets under `pending_screen` are NOT candidates and NOT ruled out: "
                "they have simply not been accessibility-screened. If a pending target has a high "
                "association score, call assess_structure on it before finalising, or state "
                "explicitly that you chose not to and why."
            ),
        }

    # -- handlers: lane-dependent ---------------------------------------------------

    def _census_unavailable(self, extra: str = "") -> dict[str, Any]:
        return {
            "status": "capability_unavailable",
            "capability": "cellxgene_census",
            "error": "The Census adapter (adapters/cellxgene.py, lane 2) has not landed in this build."
                     + (f" {extra}" if extra else ""),
            "interpretation": "Cellular context is UNAVAILABLE, which is different from 'the target is "
                              "not expressed'. Leave the cellular_context criterion unknown, record the "
                              "gap, and do not substitute bulk expression for single-cell evidence.",
        }

    def _summarize_cell_context(self, disease_label: str, gene_symbols: list[str],
                                max_cells: int = 10000, tissue: str | None = None) -> dict[str, Any]:
        if self.census is None:
            return self._census_unavailable()
        fn = getattr(self.census, "summarize_cell_context", None)
        if fn is None:
            return self._census_unavailable("Adapter present but summarize_cell_context is absent.")
        out = fn(disease_label=disease_label, gene_symbols=gene_symbols,
                 max_cells=min(max_cells, CENSUS_BOUNDS["max_cells_interactive"]), tissue=tissue)
        self.store.record_cell_summaries(out)
        return out

    def _discover_census_embeddings(self, census_release: str | None = None,
                                    organism: str = "homo_sapiens") -> dict[str, Any]:
        if self.census is None:
            return self._census_unavailable()
        fn = getattr(self.census, "discover_census_embeddings", None)
        if fn is None:
            return self._census_unavailable("Adapter present but discover_census_embeddings is absent.")
        return fn(census_release=census_release, organism=organism)

    def _explore_neighborhoods(self, embedding_slice_id: str, query_population: str,
                               k: int = 25, metric: str = "cosine") -> dict[str, Any]:
        if self.census is None:
            return self._census_unavailable()
        fn = getattr(self.census, "explore_cell_neighborhoods", None)
        if fn is None:
            return self._census_unavailable("Adapter present but explore_cell_neighborhoods is absent.")
        return fn(embedding_slice_id=embedding_slice_id, query_population=query_population, k=k, metric=metric)

    def _design_unavailable(self) -> dict[str, Any]:
        return {
            "status": "capability_unavailable",
            "capability": "protein_design",
            "error": "The design adapter (adapters/protein_design.py, lane 3) has not landed in this build.",
            "interpretation": "No design run can be launched. Return an evidence_report_only outcome "
                              "with the reviewed target and binding-region proposal. Do NOT describe any "
                              "sequence as a generated candidate.",
        }

    def _submit_design(self, design_id: str) -> dict[str, Any]:
        request = self.store.design_requests.get(design_id)
        if request is None:
            raise PreconditionFailed(f"No DesignRequest '{design_id}' exists.")
        # The scientific decision gate, enforced deterministically and not by the model.
        try:
            request.assert_launchable()
        except ValueError as exc:
            raise PreconditionFailed(str(exc)) from exc

        request.input_hash = request.compute_input_hash()
        prior = self.store.find_run_by_hash(request.input_hash)
        if prior is not None:
            return {
                "status": "reused_existing",
                "run_id": prior.run_id,
                "job_status": prior.status.value,
                "input_hash": request.input_hash,
                "note": "An identical request has already been submitted. Reusing that run; no second "
                        "paid job was started.",
            }
        if self.design is None:
            return self._design_unavailable()
        fn = getattr(self.design, "submit_design", None)
        if fn is None:
            return self._design_unavailable()

        out = fn(request.model_dump(by_alias=True))
        self.store.record_design_run(request, out)
        self.store.run.budgets.paid_design_batches_used += 1
        return out

    def _get_run(self, run_id: str) -> dict[str, Any]:
        if self.design is None:
            local = self.store.design_runs.get(run_id)
            if local is None:
                return self._design_unavailable()
            return {"status": "ok", "run_id": run_id, "job_status": local.status.value, "stage": local.stage,
                    "error": local.error, "artifact_refs": local.artifact_refs}
        fn = getattr(self.design, "get_run", None)
        return fn(run_id) if fn else self._design_unavailable()

    def _evaluate_candidates(self, run_id: str) -> dict[str, Any]:
        if self.design is None:
            return self._design_unavailable()
        fn = getattr(self.design, "evaluate_candidates", None)
        return fn(run_id) if fn else self._design_unavailable()


# ---------------------------------------------------------------------------------


def _redact(arguments: dict[str, Any]) -> dict[str, Any]:
    out = {}
    for k, v in arguments.items():
        if any(s in k.lower() for s in ("token", "key", "secret", "password")):
            out[k] = "<redacted>"
        elif isinstance(v, (list, tuple)) and len(v) > 12:
            out[k] = f"<{len(v)} items>"
        else:
            out[k] = v
    return out


def _summarise(result: dict[str, Any], limit: int = 400) -> str:
    keys = ("status", "error", "accessibility_verdict", "evidence_retrieved", "rubric_version",
            "run_id", "job_status", "needs_clarification")
    bits = [f"{k}={result[k]}" for k in keys if k in result]
    if "targets" in result:
        bits.append(f"targets={len(result['targets'])}")
    if "matches" in result:
        bits.append(f"matches={len(result['matches'])}")
    s = " ".join(bits) or json.dumps(result, default=str)
    return s[:limit]
