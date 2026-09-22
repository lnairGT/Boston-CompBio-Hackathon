"""Lane 7 -- experimental design.

Takes a ``ResearchBrief``, a ``Target`` and a ``Candidate`` and generates the wet-lab
plan that would actually test the design: production route and QC, a binding-validation
concept with its orientation reasoning and controls, an orthogonal confirmation, an
MoA-specific functional concept, a specificity arm, and a pre-registered go/no-go tree.

Entry points::

    from e2b.experiments import plan_experiments, render_markdown

    plan = plan_experiments(brief, target, candidate, protein=protein,
                            accessibility=accessibility, paralogs=paralogs,
                            pathway_hints=pathways)     # -> plain dict
    md = render_markdown(plan)                           # -> markdown string

``plan_experiments`` accepts pydantic records from ``e2b.contracts`` or plain dicts with
the same field names, and always returns a JSON-able dict. It never raises on missing
inputs: it returns a ``status`` of ``plan_complete``, ``partial_*`` or ``blocked_*`` with
``failure_detail`` and ``missing_inputs`` saying exactly what was absent.
"""

from .liabilities import screen_binder_sequence, sequence_properties
from .moa import MECHANISM_CLASSES, classify_mechanism, normalise_direction
from .planner import plan_experiments
from .records import (
    AssayConcept,
    Control,
    ConstructRisk,
    DecisionCriterion,
    ExperimentalPlan,
    ProductionRoute,
    StartingPoint,
    TargetAssayProfile,
    MODULE_VERSION,
)
from .render import render_markdown
from .target_profile import TargetBiologyUnavailable, build_target_profile

__all__ = [
    "plan_experiments",
    "render_markdown",
    "build_target_profile",
    "classify_mechanism",
    "normalise_direction",
    "screen_binder_sequence",
    "sequence_properties",
    "MECHANISM_CLASSES",
    "MODULE_VERSION",
    "ExperimentalPlan",
    "TargetAssayProfile",
    "ProductionRoute",
    "AssayConcept",
    "Control",
    "ConstructRisk",
    "DecisionCriterion",
    "StartingPoint",
    "TargetBiologyUnavailable",
]
