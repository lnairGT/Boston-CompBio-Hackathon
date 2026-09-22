"""Shared typed records for the evidence-to-binder workflow.

Every record carries ``schema_version`` and, where it asserts something about the
world, a ``provenance`` block and an ``origin`` label.  The three origin labels are
load-bearing and are checked by ``tests/test_contracts.py``:

``real``
    Retrieved or computed now, from the named source, in this run.
``cached_real``
    A previously-retrieved *real* record replayed from cache.  Carries the original
    ``retrieved_at`` so the UI can show its age; it is never presented as live.
``fixture``
    Hand-written stand-in used only to unblock UI wiring.  A fixture is never a
    biological result and may never reach a scientific conclusion or a design run.

Two conventions that the spec calls out explicitly and that are enforced here rather
than left to reviewer discipline:

* Unknown numeric values are ``None``, never ``0.0``.  ``0.0`` means "measured, and
  the measurement was zero".
* Missing evidence can never become a favourable score.  ``Assessment`` therefore
  distinguishes ``unknown`` (no evidence retrieved) from ``weak`` (evidence retrieved,
  and it is weak).
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from enum import Enum
from typing import Annotated, Any, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    NonNegativeFloat,
    field_validator,
    model_validator,
)

SCHEMA_VERSION = "1.0.0"

Origin = Literal["real", "cached_real", "fixture"]
EvidenceDirection = Literal["supports", "contradicts", "context"]
Strength = Literal["strong", "moderate", "weak", "unknown"]
ReviewStatus = Literal["draft", "pending_review", "approved", "rejected"]


def utcnow() -> datetime:
    return datetime.now(timezone.utc)


def stable_hash(payload: Any) -> str:
    """Deterministic hash used for job idempotency and cache keys."""
    blob = json.dumps(payload, sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.sha256(blob.encode()).hexdigest()[:32]


class Record(BaseModel):
    """Base for every shared record."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    schema_version: str = SCHEMA_VERSION


class Provenance(Record):
    """Where a record came from, precisely enough to re-fetch it."""

    source: str = Field(description="Logical source, e.g. 'open_targets', 'cellxgene_census'.")
    source_release: str | None = Field(
        default=None,
        description="Release/version of the source, e.g. OT '26.06' or a Census build date.",
    )
    endpoint: str | None = None
    query: str | None = Field(default=None, description="Query text or tool call, verbatim.")
    variables: dict[str, Any] | None = None
    retrieved_at: datetime = Field(default_factory=utcnow)
    elapsed_ms: int | None = None
    origin: Origin = "real"
    notes: str | None = None


# --------------------------------------------------------------------------------------
# Brief
# --------------------------------------------------------------------------------------


class DiseaseMatch(Record):
    """One candidate resolution of the user's free-text indication."""

    disease_id: str
    label: str
    match_type: Literal["exact", "synonym", "parent", "child", "subtype", "related"]
    score: float | None = None
    description: str | None = None


class ResearchBrief(Record):
    brief_id: str
    raw_indication: str = Field(description="The user's own words, preserved verbatim.")

    resolved_disease_id: str | None = None
    resolved_disease_label: str | None = None
    resolution_candidates: list[DiseaseMatch] = Field(default_factory=list)
    resolution_rationale: str | None = None

    census_disease_label: str | None = Field(
        default=None,
        description="Census disease label this maps to. Resolved through an explicit mapping, "
        "never by assuming Open Targets and Census use identical strings.",
    )
    census_disease_ontology_id: str | None = None
    census_mapping_provenance: Provenance | None = None

    desired_effect: str | None = None
    moa: str | None = None
    modality: str | None = None
    species: str = "homo_sapiens"
    species_taxon_id: str = "NCBITaxon:9606"

    status: Literal["pending", "confirmed"] = "pending"
    provenance: Provenance | None = None

    @model_validator(mode="after")
    def _confirmed_needs_resolution(self) -> "ResearchBrief":
        if self.status == "confirmed" and not self.resolved_disease_id:
            raise ValueError("A confirmed brief must carry a resolved_disease_id.")
        return self


# --------------------------------------------------------------------------------------
# Target and evidence
# --------------------------------------------------------------------------------------


class Target(Record):
    target_id: str = Field(description="Open Targets / Ensembl stable gene ID.")
    gene_symbol: str
    stable_gene_id: str
    protein_accession: str | None = Field(
        default=None, description="UniProt accession. None when no reviewed accession mapped."
    )
    species: str = "homo_sapiens"
    proposed_modulation: Literal["inhibit", "activate", "degrade", "block", "unknown"] = "unknown"

    approved_name: str | None = None
    biotype: str | None = None
    subcellular_locations: list[str] = Field(default_factory=list)
    tractability: dict[str, list[str]] = Field(
        default_factory=dict, description="modality -> list of satisfied tractability buckets."
    )
    provenance: Provenance | None = None

    @field_validator("protein_accession")
    @classmethod
    def _accession_shape(cls, v: str | None) -> str | None:
        if v is not None and not v.strip():
            raise ValueError("protein_accession must be a real accession or None, not empty.")
        return v


class Evidence(Record):
    """One source-backed claim. Every scientific relation in the graph points at one of these."""

    evidence_id: str
    subject_id: str
    predicate: str
    object_id: str | None = None
    object_value: float | str | None = None

    evidence_type: Literal[
        "genetic_association",
        "somatic_mutation",
        "known_drug",
        "clinical",
        "affected_pathway",
        "rna_expression",
        "animal_model",
        "literature",
        "single_cell_expression",
        "structural",
        "annotation",
        "computed",
    ]
    direction: EvidenceDirection

    source: str
    source_url: str | None = None
    doi: str | None = None
    source_locator: str | None = Field(
        default=None, description="Where in the source, e.g. a section, figure or datasource ID."
    )
    supporting_excerpt: str | None = None
    measurement: dict[str, Any] | None = None

    disease_scope: Literal["direct", "indirect_propagated", "unknown"] = "unknown"
    scope_disease_id: str | None = Field(
        default=None, description="The disease this item ACTUALLY refers to, before propagation."
    )

    limitations: list[str] = Field(default_factory=list)
    retrieved_at: datetime = Field(default_factory=utcnow)
    origin: Origin = "real"

    @model_validator(mode="after")
    def _needs_a_referent(self) -> "Evidence":
        if self.object_id is None and self.object_value is None:
            raise ValueError(
                f"Evidence {self.evidence_id} asserts nothing: object_id and object_value are both None."
            )
        if self.origin in ("real", "cached_real") and not (self.source_url or self.doi or self.source_locator):
            raise ValueError(
                f"Evidence {self.evidence_id} claims origin={self.origin} but has no locator "
                "(source_url, doi or source_locator). Unsourced claims are not admissible."
            )
        return self


# --------------------------------------------------------------------------------------
# Cellular context
# --------------------------------------------------------------------------------------


class CellSummary(Record):
    summary_id: str
    target_id: str
    dataset_ids: list[str] = Field(default_factory=list)
    dataset_version: str | None = None
    census_release: str

    query_filter: dict[str, Any] = Field(default_factory=dict)
    tissue: str | None = None
    condition: str | None = None

    canonical_population_id: str | None = None
    canonical_population_label: str | None = None
    original_labels: list[str] = Field(default_factory=list)

    donor_counts: dict[str, int] = Field(
        default_factory=dict, description="Donor ID (namespaced by dataset) -> cell count."
    )
    n_donors: int | None = None
    cell_counts: int | None = None

    expression_metric: str = Field(
        description="e.g. 'percent_positive_raw_counts_gt0' or 'mean_log1p_cp10k'."
    )
    expression_method: str
    values: dict[str, float | None] = Field(
        default_factory=dict, description="Per-donor or per-population values; None where not measurable."
    )
    per_donor_values: dict[str, float | None] = Field(default_factory=dict)
    dispersion: dict[str, float | None] = Field(default_factory=dict)

    primary_data_handling: str = Field(
        description="How is_primary_data was applied, to avoid duplicate cell counting."
    )
    limitations: list[str] = Field(default_factory=list)
    provenance: Provenance | None = None
    origin: Origin = "real"

    @model_validator(mode="after")
    def _rna_caveat_present(self) -> "CellSummary":
        joined = " ".join(self.limitations).lower()
        if "surface protein" not in joined:
            self.limitations.append(
                "RNA expression is not surface protein abundance; undetected transcript is not proof of absence."
            )
        return self


class CellSlice(Record):
    slice_id: str
    census_release: str
    organism: str
    dataset_ids: list[str] = Field(default_factory=list)
    filters: dict[str, Any] = Field(default_factory=dict)
    soma_joinids_ref: str | None = Field(default=None, description="Artifact reference, not inline IDs.")
    eligible_cell_count: int | None = None
    sampled_cell_count: int | None = None
    sampling_method: str | None = None
    sampling_seed: int | None = None
    gene_ids: list[str] = Field(default_factory=list)
    queried_at: datetime = Field(default_factory=utcnow)
    provenance: Provenance | None = None


class EmbeddingSlice(Record):
    embedding_slice_id: str
    cell_slice_id: str
    embedding_name: str
    embedding_metadata: dict[str, Any] = Field(default_factory=dict)
    embedding_uri: str | None = None
    census_release: str
    model_provenance: str | None = None
    dimension: int | None = None
    aligned_cell_ids_ref: str | None = None
    vectors_ref: str | None = None
    valid_row_mask_ref: str | None = None
    requested_cells: int | None = None
    retrieved_cells: int | None = None
    valid_embedding_cells: int | None = None
    coverage_by_condition: dict[str, int] = Field(default_factory=dict)
    coverage_by_donor: dict[str, int] = Field(default_factory=dict)
    coverage_by_population: dict[str, int] = Field(default_factory=dict)
    checksum: str | None = None


class NeighborhoodEvidence(Record):
    evidence_id: str
    embedding_slice_id: str
    query_population: str
    reference_universe: str = Field(
        description="The ACTUAL search universe, e.g. '9,842 retrieved cells'. Never 'the Census'."
    )
    distance_metric: str
    preprocessing: str
    k: int
    neighbor_ids: list[int] = Field(default_factory=list)
    distances: list[float] = Field(default_factory=list)
    biological_filters: dict[str, Any] = Field(default_factory=dict)
    caveats: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def _similarity_is_not_causality(self) -> "NeighborhoodEvidence":
        joined = " ".join(self.caveats).lower()
        if "not proof" not in joined:
            self.caveats.append(
                "Embedding similarity is computed context. It is not proof of shared identity, "
                "target causality, target specificity or therapeutic response."
            )
        if len(self.neighbor_ids) != len(self.distances):
            raise ValueError("neighbor_ids and distances must be the same length.")
        return self


# --------------------------------------------------------------------------------------
# Assessment and ranking
# --------------------------------------------------------------------------------------


class CriterionRating(Record):
    criterion: str
    rating: Strength
    evidence_ids: list[str] = Field(default_factory=list)
    missingness: str | None = Field(
        default=None, description="What was looked for and not found, when rating is 'unknown'."
    )
    rationale: str | None = None

    @model_validator(mode="after")
    def _rating_needs_backing(self) -> "CriterionRating":
        if self.rating in ("strong", "moderate") and not self.evidence_ids:
            raise ValueError(
                f"Criterion '{self.criterion}' is rated {self.rating} with no evidence_ids. "
                "Missing evidence cannot become a favourable score."
            )
        if self.rating == "unknown" and not self.missingness:
            raise ValueError(
                f"Criterion '{self.criterion}' is 'unknown' but does not say what was missing."
            )
        return self


class DesignFeasibility(Record):
    """Kept deliberately separate from biological support, per spec section 9."""

    accessible_to_modality: bool | None = None
    localization_basis: list[str] = Field(default_factory=list)
    structure_available: bool | None = None
    structure_ids: list[str] = Field(default_factory=list)
    blocking_reasons: list[str] = Field(default_factory=list)
    rating: Strength = "unknown"
    notes: str | None = None


class Assessment(Record):
    target_id: str
    brief_id: str
    criteria: list[CriterionRating] = Field(default_factory=list)
    disease_rationale: str | None = None
    moa_rationale: str | None = None
    direction_compatibility: Literal["compatible", "incompatible", "unknown"] = "unknown"
    liabilities: list[str] = Field(default_factory=list)
    design_feasibility: DesignFeasibility = Field(default_factory=DesignFeasibility)
    rubric_version: str
    review_status: ReviewStatus = "draft"
    excluded: bool = False
    exclusion_reason: str | None = None

    @model_validator(mode="after")
    def _exclusion_needs_reason(self) -> "Assessment":
        if self.excluded and not self.exclusion_reason:
            raise ValueError("An excluded target must record why it was excluded.")
        return self


class ComparisonResult(Record):
    brief_id: str
    rubric_version: str
    ordered_target_ids: list[str] = Field(default_factory=list)
    scores: dict[str, float] = Field(default_factory=dict)
    ties: list[list[str]] = Field(default_factory=list)
    excluded: dict[str, str] = Field(
        default_factory=dict, description="Target -> reason. Screened and ruled out."
    )
    pending_screen: dict[str, str] = Field(
        default_factory=dict,
        description="Target -> what is still missing. NOT yet accessibility-screened, so not "
        "eligible for the shortlist. Distinct from `excluded`, which means screened and "
        "rejected: an unscreened target must never be presented as a finalist on the strength "
        "of having no adverse evidence yet.",
    )
    unresolved_criteria: dict[str, list[str]] = Field(default_factory=dict)
    insufficient_evidence: bool = False
    explanation: str | None = None
    candidate_pool_scope: dict[str, Any] = Field(default_factory=dict)


# --------------------------------------------------------------------------------------
# Design
# --------------------------------------------------------------------------------------


class JobStatus(str, Enum):
    DRAFT = "draft"
    AWAITING_SELECTION = "awaiting_selection"
    READY = "ready"
    QUEUED = "queued"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    TIMED_OUT = "timed_out"
    CANCELLED = "cancelled"


TERMINAL_STATES = {JobStatus.SUCCEEDED, JobStatus.FAILED, JobStatus.TIMED_OUT, JobStatus.CANCELLED}


class BindingRegion(Record):
    description: str
    chain: str
    residues: list[int] = Field(default_factory=list)
    hotspot_residues: list[int] = Field(default_factory=list)
    numbering_scheme: str = "author"
    rationale: str = Field(description="Why engagement here could produce the intended MoA.")


class DesignRequest(Record):
    design_id: str
    brief_id: str
    selected_target_id: str

    protein_sequence_ref: str | None = None
    protein_sequence_checksum: str | None = None
    structure_id: str | None = None
    structure_version: str | None = None
    chain_map: dict[str, str] = Field(default_factory=dict)
    binding_region: BindingRegion | None = None

    moa_rationale: str | None = None
    evidence_ids: list[str] = Field(default_factory=list)

    model_name: str | None = None
    model_config_: dict[str, Any] = Field(default_factory=dict, alias="model_config")
    n_designs: int = 8

    review_status: ReviewStatus = "draft"
    reviewed_by: str | None = None
    reviewed_at: datetime | None = None
    input_hash: str | None = None

    model_config = ConfigDict(extra="forbid", validate_assignment=True, populate_by_name=True)

    def compute_input_hash(self) -> str:
        return stable_hash(
            {
                "target": self.selected_target_id,
                "structure": (self.structure_id, self.structure_version),
                "region": self.binding_region.model_dump() if self.binding_region else None,
                "model": self.model_name,
                "config": self.model_config_,
                "n": self.n_designs,
            }
        )

    def assert_launchable(self) -> None:
        """The scientific decision gate. Only an approved, complete request may spend money."""
        problems: list[str] = []
        if self.review_status != "approved":
            problems.append(f"review_status is '{self.review_status}', not 'approved'")
        if not self.structure_id:
            problems.append("no structure_id")
        if self.binding_region is None:
            problems.append("no binding_region")
        elif not self.binding_region.hotspot_residues:
            problems.append("binding_region has no hotspot residues")
        elif not self.binding_region.rationale.strip():
            problems.append("binding_region has no MoA rationale")
        if not self.evidence_ids:
            problems.append("no evidence_ids linking the request to the evidence that justified it")
        if not self.protein_sequence_checksum:
            problems.append("no protein_sequence_checksum")
        if problems:
            raise ValueError("DesignRequest is not launchable: " + "; ".join(problems))


class DesignRun(Record):
    run_id: str
    design_id: str
    execution_route: str = Field(description="e.g. 'modal:gpu:a100' or 'external:<service>'.")
    execution_handle: str | None = None
    status: JobStatus = JobStatus.DRAFT
    stage: str | None = None
    started_at: datetime | None = None
    updated_at: datetime | None = None
    finished_at: datetime | None = None
    error: str | None = None
    input_hash: str | None = None
    model_name: str | None = None
    model_version: str | None = None
    parameters: dict[str, Any] = Field(default_factory=dict)
    seed: int | None = None
    git_commit: str | None = None
    artifact_refs: list[str] = Field(default_factory=list)
    cost_usd_estimate: float | None = None
    origin: Origin = "real"

    @model_validator(mode="after")
    def _terminal_consistency(self) -> "DesignRun":
        if self.status in (JobStatus.FAILED, JobStatus.TIMED_OUT) and not self.error:
            raise ValueError(f"Run {self.run_id} is {self.status.value} but records no error.")
        if self.status == JobStatus.SUCCEEDED and self.error:
            raise ValueError(f"Run {self.run_id} succeeded but carries an error string.")
        return self


class Metric(Record):
    metric_name: str
    value: float | None = Field(default=None, description="None when not computed. Never 0.0 as a stand-in.")
    method: str
    model_name: str
    model_version: str | None = None
    higher_is_better: bool | None = None
    interpretation: str | None = None

    @model_validator(mode="after")
    def _no_fabricated_affinity(self) -> "Metric":
        banned = ("kd", "ic50", "ec50", "affinity_nm", "binding_affinity")
        if self.metric_name.strip().lower().replace(" ", "_") in banned:
            raise ValueError(
                f"'{self.metric_name}' is a measured binding quantity. A computational model score "
                "may not be recorded under that name. Use the model's own score name."
            )
        return self


class Candidate(Record):
    candidate_id: str
    run_id: str
    design_id: str | None = None
    target_id: str | None = Field(
        default=None,
        description="The target ACTUALLY submitted. Stays pinned even if later evidence changes the "
        "recommendation.",
    )
    sequence: str | None = None
    sequence_length: int | None = None
    structure_ref: str | None = None
    metrics: list[Metric] = Field(default_factory=list)
    evaluation_status: Literal["pending", "passed", "failed_filters", "not_evaluated", "errored"] = "pending"
    limitations: list[str] = Field(default_factory=list)
    proposed_validation: list[str] = Field(default_factory=list)
    origin: Origin = "real"

    @field_validator("sequence")
    @classmethod
    def _valid_aa(cls, v: str | None) -> str | None:
        if v is None:
            return v
        allowed = set("ACDEFGHIKLMNPQRSTVWY")
        bad = sorted(set(v.upper()) - allowed)
        if bad:
            raise ValueError(f"Sequence contains non-standard residues: {bad}")
        return v.upper()

    @model_validator(mode="after")
    def _length_matches(self) -> "Candidate":
        if self.sequence is not None:
            n = len(self.sequence)
            if self.sequence_length is None:
                self.sequence_length = n
            elif self.sequence_length != n:
                raise ValueError("sequence_length does not match the sequence.")
        return self


# --------------------------------------------------------------------------------------
# Agent run
# --------------------------------------------------------------------------------------


AgentOutcome = Literal[
    "completed_candidate_report",
    "evidence_report_only",
    "needs_clarification",
    "insufficient_evidence",
    "no_binder_feasible_target",
    "source_unavailable",
    "budget_exhausted",
    "design_failed",
    "embedding_unavailable",
    "in_progress",
]


class ToolAction(Record):
    step: int
    tool_name: str
    tool_version: str | None = None
    arguments: dict[str, Any] = Field(default_factory=dict)
    reason: str = Field(description="Brief action rationale. Not private chain-of-thought.")
    result_refs: list[str] = Field(default_factory=list)
    result_summary: str | None = None
    elapsed_ms: int | None = None
    error: str | None = None
    started_at: datetime = Field(default_factory=utcnow)


class Budgets(Record):
    max_tool_decisions: int = 30
    tool_decisions_used: int = 0
    max_transient_retries: int = 2
    max_paid_design_batches: int = 1
    paid_design_batches_used: int = 0
    max_candidate_pool: int = 20
    poll_wall_clock_s: int = 1800

    def remaining(self) -> dict[str, int]:
        return {
            "tool_decisions": self.max_tool_decisions - self.tool_decisions_used,
            "paid_design_batches": self.max_paid_design_batches - self.paid_design_batches_used,
        }

    def exhausted(self) -> bool:
        return self.tool_decisions_used >= self.max_tool_decisions


class AgentRun(Record):
    run_id: str
    brief_id: str
    current_stage: str = "resolve_indication"
    status: Literal["running", "paused", "completed", "failed"] = "running"
    plan_summary: str | None = None
    pending_questions: list[str] = Field(default_factory=list)
    tool_action_log: list[ToolAction] = Field(default_factory=list)
    evidence_ids: list[str] = Field(default_factory=list)
    excluded_targets: dict[str, str] = Field(default_factory=dict)
    candidate_pool_scope: dict[str, Any] = Field(default_factory=dict)
    budgets: Budgets = Field(default_factory=Budgets)
    design_run_ids: list[str] = Field(default_factory=list)
    final_outcome: AgentOutcome = "in_progress"
    outcome_detail: str | None = None
    started_at: datetime = Field(default_factory=utcnow)
    updated_at: datetime = Field(default_factory=utcnow)

    @model_validator(mode="after")
    def _completed_needs_outcome(self) -> "AgentRun":
        if self.status == "completed" and self.final_outcome == "in_progress":
            raise ValueError("A completed AgentRun must record a real final_outcome.")
        return self


__all__ = [n for n in dir() if n[0].isupper() or n in ("SCHEMA_VERSION", "stable_hash", "utcnow")]
