"""Persisted agent-run state.

The agent loop observes *persisted state*, not its own conversation. That distinction is
what makes a run inspectable and resumable: the store is the source of truth, the
transcript is a by-product. A user can reopen a saved run after resolving a blocker and the
loop picks up from ``current_stage`` with the evidence graph intact.

One writer per run directory, and each run gets its own immutable subdirectory. There is
deliberately no shared SQLite database here -- the spec warns against using one on a shared
volume as a job coordinator, and per-run manifests with a single writer are enough.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from ..config import RUBRIC_VERSION
from ..contracts import (
    AgentRun,
    Assessment,
    Budgets,
    CellSummary,
    ComparisonResult,
    DesignRequest,
    DesignRun,
    DiseaseMatch,
    Evidence,
    JobStatus,
    Provenance,
    ResearchBrief,
    Target,
    ToolAction,
    utcnow,
)
from ..pipeline import graph as G
from ..pipeline import ranking as rk


class RunStore:
    def __init__(self, run_id: str, brief: ResearchBrief, root: str | Path = "runs",
                 budgets: Budgets | None = None):
        self.root = Path(root) / run_id
        self.root.mkdir(parents=True, exist_ok=True)
        self.brief = brief
        self.run = AgentRun(run_id=run_id, brief_id=brief.brief_id, budgets=budgets or Budgets())

        self.targets: dict[str, Target] = {}
        self.associations: dict[str, dict[str, Any]] = {}
        self.pool_scope: dict[str, Any] = {}
        self.structures: dict[str, dict[str, Any]] = {}
        self.annotations: dict[str, dict[str, Any]] = {}
        self.evidence: dict[str, list[Evidence]] = {}
        self.cell_summaries: list[CellSummary] = []
        self.assessments: dict[str, Assessment] = {}
        self.comparison: ComparisonResult | None = None
        self.design_requests: dict[str, DesignRequest] = {}
        self.design_runs: dict[str, DesignRun] = {}
        self.candidates: dict[str, list[Any]] = {}
        self.binding_regions: dict[str, dict[str, Any]] = {}
        self.triage: dict[str, dict[str, Any]] = {}
        self.experiment_plans: dict[str, dict[str, Any]] = {}
        self.provenance: list[Provenance] = []

        self.graph = G.EvidenceGraph()
        self.graph.add_node(brief.brief_id, "Brief", f"Brief: {brief.raw_indication}",
                            raw_indication=brief.raw_indication, moa=brief.moa,
                            desired_effect=brief.desired_effect, modality=brief.modality)

    # -- recording ------------------------------------------------------------------

    def record_resolution(self, raw_text: str, matches: list[DiseaseMatch], prov: Provenance) -> None:
        self.brief.resolution_candidates = matches
        self.provenance.append(prov)
        top = matches[0]
        self.brief.resolved_disease_id = top.disease_id
        self.brief.resolved_disease_label = top.label
        self.brief.resolution_rationale = (
            f"Top match '{top.label}' ({top.disease_id}, {top.match_type}); "
            f"{len(matches)} candidate(s) considered."
        )
        self.graph.add_node(top.disease_id, "Indication", top.label, match_type=top.match_type)
        self.graph.add_relation(self.brief.brief_id, "brief_resolves_to", top.disease_id,
                                claim_class="structural", origin="real")
        self._flush("brief.json", self.brief)

    def record_pool(self, targets: list[Target], meta: dict[str, Any], prov: Provenance) -> None:
        self.provenance.append(prov)
        self.pool_scope = meta["scope"]
        self.run.candidate_pool_scope = meta["scope"]
        for t in targets:
            self.targets[t.target_id] = t
            self.associations[t.target_id] = meta["associations"][t.target_id]
            self.graph.add_node(t.target_id, "Target", t.gene_symbol,
                                approved_name=t.approved_name, accession=t.protein_accession,
                                locations=t.subcellular_locations[:4])
        self._flush("pool.json", {"scope": meta["scope"], "associations": meta["associations"]})

    def record_structure(self, accession: str, assessment: dict[str, Any], evs: list[Evidence]) -> None:
        self.structures[accession] = assessment
        self.graph.add_node(accession, "Protein", assessment.get("protein_name") or accession,
                            length=assessment.get("length"),
                            accessibility=assessment["accessibility"]["verdict"])
        for ev in evs:
            self.graph.add_evidence(ev)
            self.run.evidence_ids.append(ev.evidence_id)
        tid = next((t.target_id for t in self.targets.values() if t.protein_accession == accession), None)
        if tid:
            self.graph.add_relation(tid, "encodes_protein", accession,
                                    evidence_ids=[e.evidence_id for e in evs],
                                    claim_class="observed_metadata")
        for s in (assessment.get("experimental_structures") or [])[:4]:
            if s.get("pdb_id"):
                self.graph.add_node(s["pdb_id"], "Structure", f"PDB {s['pdb_id']}",
                                    method=s.get("method"), resolution_A=s.get("resolution_A"))
                self.graph.add_relation(s["pdb_id"], "structure_represents", accession,
                                        claim_class="structural")
        self._flush(f"structure_{accession}.json", assessment)

    def record_evidence(self, target_id: str, annotation: dict[str, Any],
                        evs: list[Evidence], prov: Provenance) -> None:
        self.provenance.append(prov)
        deduped, collapsed = G.dedupe_evidence(evs)
        self.annotations[target_id] = annotation
        self.evidence[target_id] = deduped
        disease_id = self.brief.resolved_disease_id
        for ev in deduped:
            self.graph.add_evidence(ev)
            self.run.evidence_ids.append(ev.evidence_id)
            if disease_id and disease_id in self.graph.nodes:
                self.graph.add_relation(
                    target_id, "associated_with", disease_id, evidence_ids=[ev.evidence_id],
                    claim_class=("asserted_literature" if ev.evidence_type == "literature"
                                 else "observed_metadata"),
                    direction=ev.direction, disease_scope=ev.disease_scope,
                )
        for drug in (annotation.get("known_drugs") or []):
            did = drug.get("drug_id")
            if not did:
                continue
            self.graph.add_node(did, "Drug", drug.get("drug_name") or did,
                                drug_type=drug.get("drug_type"),
                                max_clinical_stage=drug.get("max_clinical_stage"))
            drug_ev = [e.evidence_id for e in deduped if e.evidence_type in ("clinical", "known_drug")]
            if drug_ev:
                self.graph.add_relation(did, "acts_on", target_id, evidence_ids=drug_ev[:5],
                                        claim_class="observed_metadata",
                                        mechanisms=[m["action_type"] for m in drug.get("mechanisms") or []])
        self._flush(f"evidence_{target_id}.json",
                    {"annotation": annotation, "n_evidence": len(deduped), "collapsed_duplicates": collapsed})

    def record_cell_summaries(self, payload: dict[str, Any]) -> None:
        for raw in payload.get("cell_summaries", []) or []:
            try:
                cs = CellSummary(**raw) if isinstance(raw, dict) else raw
            except Exception:
                continue
            self.cell_summaries.append(cs)
            pop = self.graph.population_node_id(
                cs.canonical_population_label or "unlabelled",
                (cs.dataset_ids or ["unknown_dataset"])[0], cs.condition, cs.census_release)
            self.graph.add_node(pop, "Population", cs.canonical_population_label or "unlabelled",
                                condition=cs.condition, tissue=cs.tissue,
                                census_release=cs.census_release, n_donors=cs.n_donors,
                                n_cells=cs.cell_counts)
            ev = Evidence(
                evidence_id=f"census:{cs.summary_id}",
                subject_id=cs.target_id, predicate="expressed_in", object_id=pop,
                object_value=max([v for v in cs.values.values() if v is not None], default=None),
                evidence_type="single_cell_expression", direction="context",
                source=f"cellxgene_census:{cs.census_release}",
                source_locator=f"metric={cs.expression_metric} donors={cs.n_donors}",
                measurement={"metric": cs.expression_metric, "method": cs.expression_method,
                             "values": cs.values, "per_donor": cs.per_donor_values},
                limitations=cs.limitations, origin=cs.origin,
            )
            self.graph.add_evidence(ev)
            self.run.evidence_ids.append(ev.evidence_id)
            if cs.target_id in self.graph.nodes:
                self.graph.add_relation(cs.target_id, "expressed_in", pop, evidence_ids=[ev.evidence_id],
                                        claim_class="computed_measurement")
        self._flush("cell_summaries.json", [json.loads(c.model_dump_json()) for c in self.cell_summaries])

    def record_candidates(self, run_id: str, candidates: list[Any],
                          binding_region: dict[str, Any] | None = None) -> None:
        """Attach real design output to its run. Fixtures are rejected here, not downstream."""
        real = []
        for cand in candidates:
            origin = getattr(cand, "origin", None) or (cand.get("origin") if isinstance(cand, dict) else None)
            if origin == "fixture":
                # A fixture may never reach triage, a ranking, or a report. Dropping it here
                # means no downstream module has to defend against it.
                continue
            real.append(cand)
        self.candidates[run_id] = real
        if binding_region:
            self.binding_regions[run_id] = binding_region
        for cand in real:
            cid = getattr(cand, "candidate_id", None) or (cand.get("candidate_id") if isinstance(cand, dict) else None)
            if not cid:
                continue
            self.graph.add_node(cid, "Candidate", cid, run_id=run_id)
            self.graph.add_relation(cid, "candidate_generated_by", run_id if run_id in self.graph.nodes
                                    else self.graph.add_node(run_id, "DesignRun", run_id),
                                    claim_class="structural")
        self._flush(f"candidates_{run_id}.json",
                    [c.model_dump() if hasattr(c, "model_dump") else c for c in real])

    def record_triage(self, run_id: str, report: dict[str, Any]) -> None:
        self.triage[run_id] = report
        self._flush(f"triage_{run_id}.json", report)

    def record_experiment_plan(self, target_id: str, plan: dict[str, Any]) -> None:
        self.experiment_plans[target_id] = plan
        self._flush(f"experiment_plan_{target_id}.json", plan)

    # -- derived --------------------------------------------------------------------

    def build_assessments(self) -> list[Assessment]:
        self.assessments = {}
        for tid, target in self.targets.items():
            acc = target.protein_accession
            assessment = rk.build_assessment(
                brief=self.brief, target=target, association=self.associations.get(tid, {}),
                ot_annotation=self.annotations.get(tid), evidences=self.evidence.get(tid, []),
                structure_assessment=self.structures.get(acc) if acc else None,
                cell_summaries=self.cell_summaries,
            )
            self.assessments[tid] = assessment
            if assessment.excluded:
                self.run.excluded_targets[tid] = assessment.exclusion_reason or "excluded"
        return list(self.assessments.values())

    def run_comparison(self, shortlist: int = 3) -> ComparisonResult:
        self.build_assessments()
        self.comparison = rk.compare_targets(
            list(self.assessments.values()), self.brief.brief_id,
            shortlist=shortlist, pool_scope=self.pool_scope)
        self._flush("comparison.json", self.comparison)
        self._flush("assessments.json", [json.loads(a.model_dump_json()) for a in self.assessments.values()])
        return self.comparison

    def find_run_by_hash(self, input_hash: str) -> DesignRun | None:
        for r in self.design_runs.values():
            if r.input_hash == input_hash and r.status != JobStatus.FAILED:
                return r
        return None

    def add_design_request(self, request: DesignRequest) -> DesignRequest:
        request.input_hash = request.compute_input_hash()
        self.design_requests[request.design_id] = request
        self.graph.add_node(request.design_id, "DesignRun", f"Design {request.design_id}",
                            model=request.model_name, review_status=request.review_status,
                            input_hash=request.input_hash)
        self.graph.add_relation(request.design_id, "design_run_submitted_for",
                                request.selected_target_id, claim_class="structural")
        # The design is pinned to the evidence that justified it, so a later change of
        # recommendation cannot retroactively re-attribute an existing run.
        for eid in request.evidence_ids:
            if eid in self.graph.evidence:
                self.graph.add_relation(request.design_id, "justified_by", eid,
                                        evidence_ids=[eid], claim_class="observed_metadata")
        self._flush(f"design_request_{request.design_id}.json", request)
        return request

    def record_design_run(self, request: DesignRequest, submit_result: dict[str, Any]) -> DesignRun:
        run = DesignRun(
            run_id=submit_result["run_id"], design_id=request.design_id,
            execution_route=submit_result.get("execution_route", "unknown"),
            execution_handle=submit_result.get("execution_handle"),
            status=JobStatus(submit_result.get("job_status", "queued")),
            stage=submit_result.get("stage"), started_at=utcnow(), updated_at=utcnow(),
            input_hash=request.input_hash, model_name=submit_result.get("model_name"),
            model_version=submit_result.get("model_version"),
            parameters=submit_result.get("parameters", {}), seed=submit_result.get("seed"),
            git_commit=os.environ.get("E2B_GIT_COMMIT"),
        )
        self.design_runs[run.run_id] = run
        self.run.design_run_ids.append(run.run_id)
        self._flush(f"design_run_{run.run_id}.json", run)
        return run

    # -- bookkeeping ----------------------------------------------------------------

    def log_action(self, action: ToolAction) -> None:
        self.run.tool_action_log.append(action)
        self.run.budgets.tool_decisions_used += 1
        self.run.updated_at = utcnow()

    def set_stage(self, stage: str) -> None:
        self.run.current_stage = stage
        self.run.updated_at = utcnow()

    def finish(self, outcome: str, detail: str | None = None) -> None:
        self.run.final_outcome = outcome  # type: ignore[assignment]
        self.run.outcome_detail = detail
        self.run.status = "completed"
        self.run.updated_at = utcnow()
        self.save()

    def save(self) -> dict[str, str]:
        paths = {
            "agent_run": self._flush("agent_run.json", self.run),
            "graph": self._write("graph.json", self.graph.to_json()),
            "manifest": self._flush("manifest.json", self.manifest()),
        }
        return paths

    def manifest(self) -> dict[str, Any]:
        return {
            "run_id": self.run.run_id,
            "brief_id": self.brief.brief_id,
            "raw_indication": self.brief.raw_indication,
            "resolved_disease": {"id": self.brief.resolved_disease_id, "label": self.brief.resolved_disease_label},
            "moa": self.brief.moa,
            "modality": self.brief.modality,
            "rubric_version": RUBRIC_VERSION,
            "candidate_pool_scope": self.pool_scope,
            "targets_retrieved": len(self.targets),
            "structures_assessed": len(self.structures),
            "targets_investigated_in_depth": len(self.evidence),
            "evidence_records": len(self.graph.evidence),
            "cell_summaries": len(self.cell_summaries),
            "excluded_targets": self.run.excluded_targets,
            "shortlist": (self.comparison.ordered_target_ids[:3] if self.comparison else []),
            "ties": (self.comparison.ties if self.comparison else []),
            "design_runs": list(self.design_runs),
            "tool_decisions_used": self.run.budgets.tool_decisions_used,
            "budgets_remaining": self.run.budgets.remaining(),
            "final_outcome": self.run.final_outcome,
            "graph_stats": self.graph.stats(),
            "sources": sorted({p.source for p in self.provenance}),
            "generated_at": utcnow().isoformat(),
        }

    def _flush(self, name: str, obj: Any) -> str:
        payload = obj.model_dump_json(indent=2) if hasattr(obj, "model_dump_json") else json.dumps(
            obj, indent=2, default=str)
        return self._write(name, payload)

    def _write(self, name: str, text: str) -> str:
        # Write-then-rename so a reader never sees a half-written manifest.
        tmp = self.root / f".{name}.tmp"
        tmp.write_text(text)
        final = self.root / name
        os.replace(tmp, final)
        return str(final)
