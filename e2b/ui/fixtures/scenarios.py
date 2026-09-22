"""Scenario builders for the UI walkthrough.

Six scenarios, all ``origin="fixture"``:

``walkthrough``   a complete run through to candidates (job succeeded)
``ambiguous``     indication resolution is uncertain; the choices are shown
``pending``       design job queued/running — pending state
``no_data``       brief resolved, no evidence retrieved anywhere — no-data state
``failed``        design job failed with a real error string
``timed_out``     retrieval timed out AND the design job timed out

A scenario is loaded into the store under its own run key, so a demo can switch
between them without losing the others.
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from typing import Any, Callable

from e2b.contracts import (
    AgentRun,
    Assessment,
    BindingRegion,
    Budgets,
    Candidate,
    CellSummary,
    ComparisonResult,
    CriterionRating,
    DesignFeasibility,
    DesignRequest,
    DesignRun,
    DiseaseMatch,
    Evidence,
    JobStatus,
    Metric,
    NeighborhoodEvidence,
    Provenance,
    ResearchBrief,
    Target,
    ToolAction,
    stable_hash,
)

FIXTURE = "fixture"
_T0 = datetime(2026, 9, 22, 14, 5, 0, tzinfo=timezone.utc)

#: Prefix on every fixture identifier so a fixture id is recognisable anywhere it leaks.
FX = "FIXTURE-"

_NOT_RETRIEVED = (
    "Fixture record. This was written by hand to exercise the interface and was never "
    "retrieved from the named source."
)


def _prov(source: str, endpoint: str | None = None, release: str | None = None,
          at: datetime | None = None, query: str | None = None) -> Provenance:
    return Provenance(
        source=source,
        source_release=release,
        endpoint=endpoint,
        query=query,
        retrieved_at=at or _T0,
        elapsed_ms=None,
        origin=FIXTURE,
        notes=_NOT_RETRIEVED,
    )


# ---------------------------------------------------------------------------------
# Shared building blocks
# ---------------------------------------------------------------------------------


def _brief(
    *,
    brief_id: str,
    resolved: bool = True,
    confirmed: bool = True,
    candidates: bool = True,
) -> ResearchBrief:
    matches = [
        DiseaseMatch(
            disease_id="MONDO_0004980",
            label="atopic eczema",
            match_type="synonym",
            score=0.94,
            description="Fixture stand-in for an ontology match.",
        ),
        DiseaseMatch(
            disease_id="MONDO_0011292",
            label="eczema herpeticum",
            match_type="child",
            score=0.41,
            description="Fixture stand-in for a narrower child term.",
        ),
        DiseaseMatch(
            disease_id="MONDO_0005271",
            label="dermatitis",
            match_type="parent",
            score=0.38,
            description="Fixture stand-in for a broader parent term.",
        ),
    ]
    return ResearchBrief(
        brief_id=brief_id,
        raw_indication="atopic dermatitis",
        resolved_disease_id="MONDO_0004980" if resolved else None,
        resolved_disease_label="atopic eczema" if resolved else None,
        resolution_candidates=matches if candidates else [],
        resolution_rationale=(
            "Fixture: three ontology candidates are shown so the ambiguity UI is exercised. "
            "In a real run this text states why one term was preferred."
        ),
        census_disease_label="atopic dermatitis" if resolved else None,
        census_disease_ontology_id="MONDO_0004980" if resolved else None,
        census_mapping_provenance=_prov("census_disease_map") if resolved else None,
        desired_effect="suppress type 2 inflammatory signalling",
        moa="antagonise",
        modality="de novo minibinder",
        species="homo_sapiens",
        species_taxon_id="NCBITaxon:9606",
        status="confirmed" if (resolved and confirmed) else "pending",
        provenance=_prov("ui_fixture"),
    )


_TARGET_SPECS = [
    ("ENSG00000077238", "IL4R", "interleukin 4 receptor", "P24394", True),
    ("ENSG00000145777", "TSLP", "thymic stromal lymphopoietin", "Q969D9", True),
    ("ENSG00000105974", "CAV1", "caveolin 1", "Q03135", False),
]


def _targets() -> list[Target]:
    out = []
    for gid, sym, name, acc, accessible in _TARGET_SPECS:
        out.append(
            Target(
                target_id=gid,
                gene_symbol=sym,
                stable_gene_id=gid,
                protein_accession=acc,
                proposed_modulation="block",
                approved_name=name,
                biotype="protein_coding",
                subcellular_locations=(
                    ["Cell membrane", "Secreted"] if accessible else ["Cell membrane; cytoplasmic side"]
                ),
                tractability={"AB": ["Predicted Tractable"]} if accessible else {},
                provenance=_prov(
                    "open_targets",
                    endpoint="https://api.platform.opentargets.org/api/v4/graphql",
                    release="fixture (no release)",
                ),
            )
        )
    return out


def _evidence() -> list[Evidence]:
    items: list[Evidence] = []
    spec: list[tuple[str, str, str, str, str, str, str | None]] = [
        ("ENSG00000077238", "genetic_association", "supports", "ot_genetics_portal",
         "associated_with", "MONDO_0004980",
         "Fixture stand-in for a genetic association record."),
        ("ENSG00000077238", "known_drug", "supports", "chembl",
         "targeted_by_approved_drug", "MONDO_0004980",
         "Fixture stand-in for an approved-drug record."),
        ("ENSG00000077238", "affected_pathway", "supports", "reactome",
         "participates_in", "type 2 cytokine signalling",
         "Fixture stand-in for a pathway membership record."),
        ("ENSG00000145777", "genetic_association", "supports", "ot_genetics_portal",
         "associated_with", "MONDO_0004980",
         "Fixture stand-in for a genetic association record."),
        ("ENSG00000145777", "animal_model", "supports", "impc",
         "model_phenotype", "dermatitis-like phenotype",
         "Fixture stand-in for an animal-model record."),
        ("ENSG00000105974", "rna_expression", "context", "expression_atlas",
         "expressed_in", "skin",
         "Fixture stand-in for a bulk expression record."),
        ("ENSG00000105974", "literature", "contradicts", "europepmc",
         "reported_liability_in", "broad tissue expression",
         "Fixture stand-in: a contradicting item, so the trail is not all one-sided."),
    ]
    for i, (subj, etype, direction, source, pred, obj, excerpt) in enumerate(spec, start=1):
        is_id = obj.startswith("MONDO") or obj.startswith("ENSG")
        items.append(
            Evidence(
                evidence_id=f"{FX}EV-{i:03d}",
                subject_id=subj,
                predicate=pred,
                object_id=obj if is_id else None,
                object_value=None if is_id else obj,
                evidence_type=etype,  # type: ignore[arg-type]
                direction=direction,  # type: ignore[arg-type]
                source=source,
                source_url=f"https://platform.opentargets.org/target/{subj}/associations",
                doi=None,
                source_locator=f"datasource:{source}",
                supporting_excerpt=excerpt,
                measurement=None,
                disease_scope="direct" if is_id else "unknown",
                scope_disease_id="MONDO_0004980" if is_id else None,
                limitations=[_NOT_RETRIEVED],
                retrieved_at=_T0,
                origin=FIXTURE,
            )
        )
    return items


def _assessments(brief_id: str) -> list[Assessment]:
    def rating(crit: str, r: str, eids: list[str], why: str, missing: str | None = None):
        return CriterionRating(
            criterion=crit, rating=r, evidence_ids=eids, rationale=why, missingness=missing  # type: ignore[arg-type]
        )

    a1 = Assessment(
        target_id="ENSG00000077238",
        brief_id=brief_id,
        criteria=[
            rating("clinical_causal_support", "strong",
                   [f"{FX}EV-001", f"{FX}EV-002"], "Fixture: association plus approved drug."),
            rating("mechanistic_support", "moderate", [f"{FX}EV-003"],
                   "Fixture: pathway membership only."),
            rating("direction_fit", "moderate", [f"{FX}EV-002"],
                   "Fixture: an approved antagonist exists, so the direction is attested."),
            rating("cellular_context", "unknown", [], None,
                   "Fixture: no Census summary was attached to this scenario."),
        ],
        disease_rationale="Fixture rationale for the disease link.",
        moa_rationale="Fixture rationale for the mechanism link.",
        direction_compatibility="compatible",
        liabilities=["Fixture liability: receptor is broadly expressed."],
        design_feasibility=DesignFeasibility(
            accessible_to_modality=True,
            localization_basis=["Signal peptide 1-25 (fixture)", "Extracellular topological domain (fixture)"],
            structure_available=True,
            structure_ids=["AF-P24394-F1 (fixture reference)"],
            rating="moderate",
            notes="Fixture: accessibility asserted from topology fields, not from a keyword screen.",
        ),
        rubric_version="rubric-2026.09.22-v1",
        review_status="draft",
    )
    a2 = Assessment(
        target_id="ENSG00000145777",
        brief_id=brief_id,
        criteria=[
            rating("clinical_causal_support", "moderate", [f"{FX}EV-004"],
                   "Fixture: association without an approved drug."),
            rating("mechanistic_support", "moderate", [f"{FX}EV-005"],
                   "Fixture: animal-model support."),
            rating("direction_fit", "unknown", [], None,
                   "Fixture: no evidence that antagonism is the beneficial direction."),
            rating("cellular_context", "unknown", [], None,
                   "Fixture: no Census summary was attached to this scenario."),
        ],
        disease_rationale="Fixture rationale.",
        moa_rationale="Fixture rationale.",
        direction_compatibility="unknown",
        liabilities=[],
        design_feasibility=DesignFeasibility(
            accessible_to_modality=True,
            localization_basis=["Secreted (fixture)"],
            structure_available=True,
            structure_ids=["AF-Q969D9-F1 (fixture reference)"],
            rating="moderate",
            notes="Fixture.",
        ),
        rubric_version="rubric-2026.09.22-v1",
        review_status="draft",
    )
    a3 = Assessment(
        target_id="ENSG00000105974",
        brief_id=brief_id,
        criteria=[
            rating("clinical_causal_support", "unknown", [], None,
                   "Fixture: no genetic or clinical item found for this disease."),
            rating("mechanistic_support", "weak", [f"{FX}EV-006"],
                   "Fixture: expression context only."),
            rating("direction_fit", "unknown", [], None, "Fixture: nothing attests a direction."),
            rating("cellular_context", "unknown", [], None, "Fixture: no Census summary."),
        ],
        liabilities=["Fixture liability: broad tissue expression."],
        direction_compatibility="unknown",
        design_feasibility=DesignFeasibility(
            accessible_to_modality=False,
            localization_basis=["Cell membrane; cytoplasmic side (fixture) — faces the cytosol"],
            structure_available=None,
            blocking_reasons=[
                "Fixture: annotated on the cytoplasmic face, so a soluble binder cannot engage it."
            ],
            rating="weak",
        ),
        rubric_version="rubric-2026.09.22-v1",
        review_status="draft",
    )
    return [a1, a2, a3]


def _comparison(brief_id: str, *, insufficient: bool = False) -> ComparisonResult:
    return ComparisonResult(
        brief_id=brief_id,
        rubric_version="rubric-2026.09.22-v1",
        ordered_target_ids=["ENSG00000077238", "ENSG00000145777", "ENSG00000105974"],
        scores={
            "ENSG00000077238": 0.685,
            "ENSG00000145777": 0.355,
            "ENSG00000105974": 0.0625,
        },
        ties=[],
        excluded={},
        unresolved_criteria={
            "ENSG00000077238": ["cellular_context"],
            "ENSG00000145777": ["direction_fit", "cellular_context"],
            "ENSG00000105974": ["clinical_causal_support", "direction_fit", "cellular_context"],
        },
        insufficient_evidence=insufficient,
        explanation=(
            "Fixture ranking. The ordering reflects hand-written ratings, not retrieved "
            "evidence, and must not be read as a target recommendation."
        ),
        candidate_pool_scope={
            "pool": "3 fixture targets",
            "note": "Not a search over Open Targets.",
        },
    )


def _cell_summary() -> CellSummary:
    return CellSummary(
        summary_id=f"{FX}CS-001",
        target_id="ENSG00000077238",
        dataset_ids=["fixture-dataset-a", "fixture-dataset-b"],
        census_release="fixture (no Census release)",
        query_filter={"tissue_general": "skin of body", "disease": "atopic dermatitis"},
        tissue="skin of body",
        condition="atopic dermatitis",
        canonical_population_label="type 2 helper T cell",
        original_labels=["Th2", "CD4+ Th2 cell"],
        donor_counts={"fx-donor-1": 812, "fx-donor-2": 640, "fx-donor-3": 129},
        n_donors=3,
        cell_counts=1581,
        expression_metric="percent_positive_raw_counts_gt0",
        expression_method="Fixture value. No Census query was executed.",
        values={},
        per_donor_values={"fx-donor-1": 0.41, "fx-donor-2": 0.36, "fx-donor-3": None},
        dispersion={},
        primary_data_handling="Fixture: is_primary_data filtering was not applied because no "
                              "query ran.",
        limitations=[_NOT_RETRIEVED, "Donor 3 has no value: too few cells to summarise, not zero."],
        provenance=_prov("cellxgene_census"),
        origin=FIXTURE,
    )


def _neighborhood() -> NeighborhoodEvidence:
    return NeighborhoodEvidence(
        evidence_id=f"{FX}NB-001",
        embedding_slice_id=f"{FX}ES-001",
        query_population="type 2 helper T cell",
        reference_universe="1,581 fixture cells (not the Census)",
        distance_metric="cosine",
        preprocessing="Fixture: no preprocessing was performed.",
        k=5,
        neighbor_ids=[11, 27, 44, 58, 91],
        distances=[0.0812, 0.0975, 0.1103, 0.1250, 0.1388],
        biological_filters={"is_primary_data": True},
        caveats=[_NOT_RETRIEVED],
    )


def _design_request(brief_id: str, *, approved: bool) -> DesignRequest:
    req = DesignRequest(
        design_id=f"{FX}DS-001",
        brief_id=brief_id,
        selected_target_id="ENSG00000077238",
        protein_sequence_ref="fixture://sequence/P24394",
        protein_sequence_checksum=stable_hash("fixture-sequence-P24394"),
        structure_id="AF-P24394-F1 (fixture reference)",
        structure_version="fixture",
        chain_map={"A": "IL4R ectodomain (fixture)"},
        binding_region=BindingRegion(
            description="Fixture binding region on the receptor ectodomain.",
            chain="A",
            residues=list(range(90, 121)),
            hotspot_residues=[97, 101, 104, 118],
            numbering_scheme="author",
            rationale=(
                "Fixture rationale: engaging this face would block the cytokine-receptor "
                "interface and so antagonise the pathway. Not derived from a real analysis."
            ),
        ),
        moa_rationale=(
            "Fixture: antagonise type 2 signalling by occluding the cytokine-binding face."
        ),
        evidence_ids=[f"{FX}EV-001", f"{FX}EV-002", f"{FX}EV-003"],
        model_name="fixture-design-model",
        n_designs=8,
        review_status="approved" if approved else "draft",
    )
    req.input_hash = req.compute_input_hash()
    return req


def _design_run(status: JobStatus, *, origin: str = FIXTURE, **kw: Any) -> DesignRun:
    base: dict[str, Any] = dict(
        run_id=f"{FX}RUN-001",
        design_id=f"{FX}DS-001",
        execution_route="fixture:none (nothing was executed)",
        execution_handle=None,
        status=status,
        started_at=_T0,
        updated_at=_T0 + timedelta(minutes=4),
        model_name="fixture-design-model",
        model_version="fixture",
        parameters={"n_designs": 8},
        seed=20260922,
        cost_usd_estimate=None,
        origin=origin,
    )
    base.update(kw)
    return DesignRun(**base)


def _candidates() -> list[Candidate]:
    seq = (
        "SEEELKKLAEELKKLAEEIKKLAEELKKLGGSPEELLKLAEELAKKAEELAKKLEELAKKLG"
    )
    return [
        Candidate(
            candidate_id=f"{FX}CAND-001",
            run_id=f"{FX}RUN-001",
            design_id=f"{FX}DS-001",
            target_id="ENSG00000077238",
            sequence=seq,
            structure_ref="fixture://structure/cand-001.pdb",
            metrics=[
                Metric(
                    metric_name="predicted_aligned_error_interface",
                    value=8.4,
                    method="Fixture value. No model was run.",
                    model_name="fixture-structure-model",
                    model_version="fixture",
                    higher_is_better=False,
                    interpretation="Lower interface PAE indicates a more confident predicted pose.",
                ),
                Metric(
                    metric_name="plddt_binder",
                    value=82.1,
                    method="Fixture value. No model was run.",
                    model_name="fixture-structure-model",
                    model_version="fixture",
                    higher_is_better=True,
                    interpretation="Per-residue confidence of the predicted binder fold.",
                ),
                Metric(
                    metric_name="rosetta_ddg",
                    value=None,
                    method="Not computed in this fixture.",
                    model_name="fixture-energy-model",
                    higher_is_better=False,
                    interpretation="Left as null to show that an uncomputed metric is not zero.",
                ),
            ],
            evaluation_status="not_evaluated",
            limitations=[
                _NOT_RETRIEVED,
                "This sequence is a hand-written helical placeholder, not a design output.",
            ],
            proposed_validation=[
                "Fixture suggestion: biolayer interferometry against the purified ectodomain "
                "to obtain an actual Kd — the model scores above are not affinities.",
                "Fixture suggestion: a cell-based reporter assay to test whether binding "
                "antagonises signalling rather than merely occurring.",
            ],
            origin=FIXTURE,
        ),
        Candidate(
            candidate_id=f"{FX}CAND-002",
            run_id=f"{FX}RUN-001",
            design_id=f"{FX}DS-001",
            target_id="ENSG00000077238",
            sequence=None,
            structure_ref=None,
            metrics=[],
            evaluation_status="errored",
            limitations=[
                _NOT_RETRIEVED,
                "Deliberately empty: shows a candidate slot with no sequence and no metrics.",
            ],
            proposed_validation=[],
            origin=FIXTURE,
        ),
    ]


def _agent_run(brief_id: str, stage: str, *, status: str = "running",
               outcome: str = "in_progress", detail: str | None = None,
               questions: list[str] | None = None, used: int = 0,
               plan: str | None = None) -> AgentRun:
    return AgentRun(
        run_id=f"{FX}AG-{stage}",
        brief_id=brief_id,
        current_stage=stage,
        status=status,  # type: ignore[arg-type]
        plan_summary=plan,
        pending_questions=questions or [],
        tool_action_log=[
            ToolAction(
                step=1,
                tool_name="resolve_indication",
                arguments={"raw_indication": "atopic dermatitis"},
                reason="Fixture action log entry.",
                result_summary="Fixture: 3 ontology candidates returned.",
                elapsed_ms=None,
                started_at=_T0,
            )
        ],
        budgets=Budgets(tool_decisions_used=used),
        final_outcome=outcome,  # type: ignore[arg-type]
        outcome_detail=detail,
        started_at=_T0,
        updated_at=_T0 + timedelta(minutes=used or 1),
    )


# ---------------------------------------------------------------------------------
# Scenarios
# ---------------------------------------------------------------------------------


def _s_walkthrough() -> dict[str, Any]:
    bid = f"{FX}BRIEF-walkthrough"
    return {
        "scenario_label": "Fixture walkthrough — complete run, job succeeded",
        "brief": _brief(brief_id=bid),
        "agent_run": _agent_run(
            bid, "report", status="completed", outcome="completed_candidate_report",
            detail="Fixture run reached a candidate report.", used=14,
            plan="Fixture: report assembled from two candidates returned by the design run.",
        ),
        "targets": _targets(),
        "assessments": _assessments(bid),
        "comparison": _comparison(bid),
        "evidence": _evidence(),
        "cell_summaries": [_cell_summary()],
        "neighborhoods": [_neighborhood()],
        "design_request": _design_request(bid, approved=True),
        "design_run": _design_run(
            JobStatus.SUCCEEDED, stage="complete",
            finished_at=_T0 + timedelta(minutes=12),
            artifact_refs=["fixture://run/RUN-001/candidates.json"],
        ),
        "candidates": _candidates(),
        "selected_target_id": "ENSG00000077238",
        "confirmed": True,
    }


def _s_cached() -> dict[str, Any]:
    """A replayed run. Everything shows its ORIGINAL time and a cached badge."""
    bid = f"{FX}BRIEF-cached"
    older = _T0 - timedelta(days=1, hours=3)
    run = _design_run(
        JobStatus.SUCCEEDED,
        origin="cached_real",
        run_id=f"{FX}RUN-CACHED",
        started_at=older,
        updated_at=older + timedelta(minutes=9),
        finished_at=older + timedelta(minutes=9),
        stage="complete",
        execution_route="fixture:replay",
    )
    cands = _candidates()
    for c in cands:
        c.origin = "cached_real"
        c.run_id = f"{FX}RUN-CACHED"
    return {
        "scenario_label": "Fixture replay — cached result, shown with its original run time",
        "brief": _brief(brief_id=bid),
        "agent_run": _agent_run(bid, "report", status="completed",
                                outcome="completed_candidate_report",
                                detail="Replayed from cache.", used=14),
        "targets": _targets(),
        "assessments": _assessments(bid),
        "comparison": _comparison(bid),
        "evidence": _evidence(),
        "cell_summaries": [_cell_summary()],
        "design_request": _design_request(bid, approved=True),
        "design_run": run,
        "candidates": cands,
        "selected_target_id": "ENSG00000077238",
        "confirmed": True,
    }


def _s_ambiguous() -> dict[str, Any]:
    bid = f"{FX}BRIEF-ambiguous"
    return {
        "scenario_label": "Fixture — indication resolution is ambiguous, awaiting a choice",
        "brief": _brief(brief_id=bid, resolved=False, confirmed=False),
        "agent_run": _agent_run(
            bid, "resolve_indication", used=1,
            questions=[
                "Three ontology terms match 'atopic dermatitis'. Which one do you mean?"
            ],
        ),
        "targets": [],
        "evidence": [],
    }


def _s_pending() -> dict[str, Any]:
    bid = f"{FX}BRIEF-pending"
    return {
        "scenario_label": "Fixture — design job running (pending state)",
        "brief": _brief(brief_id=bid),
        "agent_run": _agent_run(
            bid, "design_run", used=11,
            plan="Fixture: polling the design run. No result is shown until one is reported.",
        ),
        "targets": _targets(),
        "assessments": _assessments(bid),
        "comparison": _comparison(bid),
        "evidence": _evidence(),
        "cell_summaries": [_cell_summary()],
        "design_request": _design_request(bid, approved=True),
        "design_run": _design_run(
            JobStatus.RUNNING, stage="backbone generation",
            updated_at=_T0 + timedelta(minutes=2), finished_at=None,
        ),
        "candidates": [],
        "selected_target_id": "ENSG00000077238",
        "confirmed": True,
    }


def _s_no_data() -> dict[str, Any]:
    bid = f"{FX}BRIEF-nodata"
    return {
        "scenario_label": "Fixture — brief resolved, nothing retrieved (no-data state)",
        "brief": _brief(brief_id=bid),
        "agent_run": _agent_run(
            bid, "assess_targets", status="completed", outcome="insufficient_evidence",
            detail="Fixture: no evidence item met the admissibility bar.", used=6,
        ),
        "targets": [],
        "assessments": [],
        "comparison": ComparisonResult(
            brief_id=bid,
            rubric_version="rubric-2026.09.22-v1",
            ordered_target_ids=[],
            insufficient_evidence=True,
            explanation="Fixture: nothing to rank. No winner is forced.",
            candidate_pool_scope={"pool": "0 targets"},
        ),
        "evidence": [],
        "failures": [
            {"source": "cellxgene_census", "kind": "empty_result",
             "detail": "Fixture: the query ran and matched no cells.",
             "query": "tissue_general == 'skin of body' AND disease == 'atopic dermatitis'"},
        ],
    }


def _s_failed() -> dict[str, Any]:
    bid = f"{FX}BRIEF-failed"
    return {
        "scenario_label": "Fixture — design job failed",
        "brief": _brief(brief_id=bid),
        "agent_run": _agent_run(
            bid, "design_run", status="failed", outcome="design_failed",
            detail="Fixture: the design executor returned a non-zero exit.", used=12,
        ),
        "targets": _targets(),
        "assessments": _assessments(bid),
        "comparison": _comparison(bid),
        "evidence": _evidence(),
        "design_request": _design_request(bid, approved=True),
        "design_run": _design_run(
            JobStatus.FAILED, stage="backbone generation",
            error="Fixture error string: executor exited 1 during backbone generation.",
            finished_at=_T0 + timedelta(minutes=3),
        ),
        "candidates": [],
        "selected_target_id": "ENSG00000077238",
        "confirmed": True,
    }


def _s_timed_out() -> dict[str, Any]:
    bid = f"{FX}BRIEF-timeout"
    return {
        "scenario_label": "Fixture — retrieval timeout and design timeout",
        "brief": _brief(brief_id=bid),
        "agent_run": _agent_run(
            bid, "design_run", status="failed", outcome="design_failed",
            detail="Fixture: wall-clock budget elapsed while the job was still running.",
            used=19,
        ),
        "targets": _targets(),
        "assessments": _assessments(bid),
        "comparison": _comparison(bid),
        "evidence": _evidence(),
        "design_request": _design_request(bid, approved=True),
        "design_run": _design_run(
            JobStatus.TIMED_OUT, stage="sequence design",
            error="Fixture: polling exceeded poll_wall_clock_s=1800 with the job still running.",
            finished_at=_T0 + timedelta(minutes=30),
        ),
        "candidates": [],
        "selected_target_id": "ENSG00000077238",
        "confirmed": True,
        "failures": [
            {"source": "cellxgene_census", "kind": "timeout",
             "detail": "Fixture: no response within the adapter timeout.",
             "query": "Census expression summary for ENSG00000077238"},
            {"source": "open_targets", "kind": "http_error",
             "detail": "Fixture: HTTP 503 from the GraphQL endpoint.",
             "query": "associatedTargets(efoId: MONDO_0004980)"},
        ],
    }


SCENARIOS: dict[str, Callable[[], dict[str, Any]]] = {
    "walkthrough": _s_walkthrough,
    "cached": _s_cached,
    "ambiguous": _s_ambiguous,
    "pending": _s_pending,
    "no_data": _s_no_data,
    "failed": _s_failed,
    "timed_out": _s_timed_out,
}

SCENARIO_DESCRIPTIONS: dict[str, str] = {
    "walkthrough": "Complete run through to candidates",
    "cached": "Replayed result, original run time shown",
    "ambiguous": "Indication resolution uncertain",
    "pending": "Design job running",
    "no_data": "Nothing retrieved",
    "failed": "Design job failed",
    "timed_out": "Retrieval and design timed out",
}


def available_scenarios() -> list[dict[str, str]]:
    return [
        {"key": k, "description": SCENARIO_DESCRIPTIONS.get(k, k)} for k in SCENARIOS
    ]


def load_scenario(store: Any, name: str, *, run_key: str | None = None) -> str:
    """Load a fixture scenario into ``store`` and return its run key."""
    if name not in SCENARIOS:
        raise KeyError(f"Unknown fixture scenario '{name}'. Known: {sorted(SCENARIOS)}")
    data = SCENARIOS[name]()
    key = run_key or f"fixture-{name}"
    store.delete(key)
    failures = data.pop("failures", [])
    confirmed = data.pop("confirmed", False)
    store.upsert(
        key,
        brief=data.get("brief"),
        agent_run=data.get("agent_run"),
        targets=data.get("targets"),
        assessments=data.get("assessments"),
        comparison=data.get("comparison"),
        evidence=data.get("evidence"),
        cell_summaries=data.get("cell_summaries"),
        neighborhoods=data.get("neighborhoods"),
        design_request=data.get("design_request"),
        design_run=data.get("design_run"),
        candidates=data.get("candidates"),
        selected_target_id=data.get("selected_target_id"),
        scenario_label=data.get("scenario_label"),
        replace_collections=True,
    )
    for f in failures:
        store.note_source_failure(key, **f)
    if confirmed:
        store.set_confirmed(key, True)
    return key
