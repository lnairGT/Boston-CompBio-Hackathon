"""Module-local typed records for the experimental plan.

These are plain dataclasses on purpose. Lane 7 returns JSON-able dicts at its
boundary so the spine can adapt without Lane 7 blocking on a shared record type,
and nothing here duplicates a concept ``e2b/contracts.py`` already owns.

Two honesty rules are enforced by the shapes themselves rather than left to
reviewer discipline:

``StartingPoint``
    Any assay parameter a bench scientist needs -- concentration span, incubation,
    surface density, temperature -- can only be expressed as a *range with
    reasoning*, tagged ``starting_point_requires_optimisation``. There is no field
    in which to record a specific validated concentration or time, because this
    module has validated none.

unknown numerics are ``None``
    Every numeric field defaults to ``None``, meaning "not determined". ``0`` would
    mean "determined, and it is zero".
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any, Literal

MODULE_VERSION = "experiments/1.0.0"

Origin = Literal["real", "cached_real", "fixture"]

PlanStatus = Literal[
    # A binding-validation route and an MoA-specific functional route were both derived.
    "plan_complete",
    # Binding route derived; the MoA readout was not derivable from the given inputs.
    "partial_no_functional_readout",
    # No candidate sequence supplied, so the construct-liability screen did not run.
    "partial_no_candidate_sequence",
    # The target is not reachable by a soluble binder, so an in-cell blockade
    # experiment cannot honestly be designed for this modality.
    "blocked_target_not_accessible",
    # Target biology could not be retrieved, so no plan is derivable.
    "blocked_no_target_biology",
    # The requested direction of modulation is not achievable by this modality.
    "blocked_direction_not_addressable",
]


def _utcnow_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass
class StartingPoint:
    """An assay parameter as a range plus the reason for that range."""

    name: str
    range_text: str
    reasoning: str
    status: str = "starting_point_requires_optimisation"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class Control:
    """One control, with what it rules out and what its failure would mean."""

    name: str
    kind: Literal["negative", "positive", "buffer_or_system", "specificity", "reagent_activity"]
    purpose: str
    interpretation_if_it_fails: str
    available: bool | None = None
    availability_note: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ConstructRisk:
    """A specific, sequence-derived liability in a construct."""

    risk: str
    evidence: str
    positions: list[int] = field(default_factory=list)
    consequence: str = ""
    mitigation: str = ""
    severity: Literal["blocking", "high", "moderate", "watch"] = "moderate"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ProductionRoute:
    """How the design and the antigen reagent get made, and the QC that gates belief."""

    subject: Literal["designed_binder", "target_antigen"]
    expression_host: str
    host_rationale: str
    construct_elements: list[str] = field(default_factory=list)
    construct_rationale: list[str] = field(default_factory=list)
    purification_steps: list[str] = field(default_factory=list)
    release_qc: list[str] = field(default_factory=list)
    qc_gates: list[str] = field(default_factory=list)
    risks: list[ConstructRisk] = field(default_factory=list)
    starting_points: list[StartingPoint] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class AssayConcept:
    """One assay concept: format, orientation, what it resolves, and its controls."""

    assay_id: str
    role: Literal[
        "binding_primary",
        "binding_orthogonal",
        "functional_moa",
        "specificity",
        "biophysical_qc",
    ]
    title: str
    format: str
    principle: str
    rationale: str
    orientation: str | None = None
    orientation_rationale: str | None = None
    readout: str | None = None
    readout_reports_mechanism_because: str | None = None
    system: str | None = None
    system_requirement: str | None = None
    resolvable_range: str | None = None
    resolvable_range_reasoning: str | None = None
    starting_points: list[StartingPoint] = field(default_factory=list)
    controls: list[Control] = field(default_factory=list)
    reports: list[str] = field(default_factory=list)
    does_not_report: list[str] = field(default_factory=list)
    derived_from: list[str] = field(default_factory=list)
    limitations: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class DecisionCriterion:
    """A pre-registered criterion: fixed before the experiment, not after the data."""

    criterion_id: str
    stage: str
    question: str
    pass_rule: str
    kill_rule: str
    rule_basis: str
    on_pass: str = ""
    on_kill: str = ""
    on_ambiguous: str | None = None
    threshold_status: str = "pre_registered_default_requires_sign_off"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class TargetAssayProfile:
    """The assay-relevant reading of the target's biology. Every field is sourced."""

    gene_symbol: str | None = None
    protein_accession: str | None = None
    protein_name: str | None = None
    archetype: str = "unknown"
    archetype_basis: list[str] = field(default_factory=list)
    mechanism_class: str = "unknown"
    mechanism_class_label: str | None = None
    mechanism_class_basis: list[str] = field(default_factory=list)
    engageable_regions: list[dict[str, Any]] = field(default_factory=list)
    length_aa: int | None = None
    n_glycosylation_sites: int | None = None
    n_disulfides: int | None = None
    experimental_structures: int | None = None
    origin: Origin = "real"
    source_notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ExperimentalPlan:
    """The Lane 7 deliverable: a wet-lab plan derived from brief + target + candidate."""

    plan_id: str
    status: str
    module_version: str = MODULE_VERSION
    generated_at: str = field(default_factory=_utcnow_iso)

    brief_id: str | None = None
    raw_indication: str | None = None
    desired_effect: str | None = None
    direction: str | None = None
    modality: str | None = None
    species: str | None = None

    target_id: str | None = None
    candidate_id: str | None = None
    run_id: str | None = None
    candidate_origin: str | None = None
    candidate_sequence_length: int | None = None

    target_profile: TargetAssayProfile | None = None
    production_routes: list[ProductionRoute] = field(default_factory=list)
    assays: list[AssayConcept] = field(default_factory=list)
    decision_tree: list[DecisionCriterion] = field(default_factory=list)

    what_this_does_not_establish: list[str] = field(default_factory=list)
    scope_statement: str = ""
    failure_detail: str | None = None
    missing_inputs: list[str] = field(default_factory=list)
    provisional_claims: list[str] = field(default_factory=list)
    inputs_seen: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def assay(self, role: str) -> AssayConcept | None:
        for a in self.assays:
            if a.role == role:
                return a
        return None
