"""Runtime budgets and the versioned ranking rubric.

The rubric lives here, in data, for two reasons the spec is explicit about: the
comparison must be *deterministic and versioned*, and the language model must
**explain** a ranking rather than silently invent numeric weights. Nothing in this
module is indication-specific, and there is no target, disease or dataset allowlist
anywhere in the production path.
"""

from __future__ import annotations

from typing import Any

# Bumped from v1 when the datatype mapping was corrected against the live 26.06 vocabulary
# and the text-mined-literature cap was added. Rankings recorded under v1 are not
# comparable to v2, which is exactly why this string travels with every Assessment.
RUBRIC_VERSION = "rubric-2026.09.22-v2"

# Ordinal map for the four strength levels. `unknown` scores 0.0 and can never lift a
# target above one with retrieved evidence -- missing evidence is not a favourable score.
STRENGTH_POINTS: dict[str, float] = {
    "strong": 1.0,
    "moderate": 0.6,
    "weak": 0.25,
    "unknown": 0.0,
}

# Criteria are scored independently and combined by a fixed weighted sum.
# `design_feasibility` is deliberately NOT in this list: section 9 requires design
# feasibility to be shown separately from biological support, and it acts as a gate
# rather than as a contribution to the biology score.
# The datatype IDs below were enumerated from the LIVE release (26.06) rather than taken
# from documentation examples: this release emits `clinical` and `genetic_literature`, and
# does NOT emit `known_drug`. `known_drug`, `somatic_mutation` and `affected_pathway` are
# retained because they appear for other indications (e.g. oncology), and an ID that is
# absent for a given disease simply contributes nothing.
CRITERIA: list[dict[str, Any]] = [
    {
        "key": "clinical_causal_support",
        "label": "Clinical and causal support",
        "weight": 0.30,
        "description": "Human genetic association and clinical/drug-precedence evidence that the "
        "target is causally involved in this indication.",
        "ot_datatypes": [
            "genetic_association",
            "genetic_literature",
            "clinical",
            "known_drug",
            "somatic_mutation",
        ],
    },
    {
        "key": "mechanistic_support",
        "label": "Mechanistic support",
        "weight": 0.25,
        "description": "Pathway, animal-model and literature evidence that the target acts through "
        "the requested mechanism.",
        "ot_datatypes": ["affected_pathway", "animal_model", "literature"],
    },
    {
        "key": "direction_fit",
        "label": "Fit to intended direction of modulation",
        "weight": 0.25,
        "description": "Evidence that modulating the target in the requested direction is the "
        "beneficial direction. Association alone never establishes this.",
        "ot_datatypes": [],
    },
    {
        "key": "cellular_context",
        "label": "Cellular context",
        "weight": 0.20,
        "description": "Expression in mechanism-relevant populations in the disease context, from "
        "CELLxGENE Census, summarised per donor.",
        "ot_datatypes": ["rna_expression"],
    },
]

# Liabilities subtract, but are capped so that a single annotation cannot dominate.
LIABILITY_PENALTY_PER_ITEM = 0.05
LIABILITY_PENALTY_CAP = 0.20

# "Literature volume alone should not determine the winner" (spec section 17). Text-mined
# co-occurrence scales with how much a gene has been written about, which tracks research
# attention rather than causal support. So when the ONLY datatype contributing to a
# criterion is text-mined literature, the rating is capped here regardless of its score.
TEXT_MINED_DATATYPES: frozenset[str] = frozenset({"literature"})
LITERATURE_ONLY_RATING_CAP = "moderate"

# Baseline expression in healthy tissue is not disease context. When the Census has no
# diseased slice for an indication -- atopic dermatitis has none in skin -- a reference-only
# summary still shows the target is present in the mechanism-relevant populations, which is
# real and useful, but it cannot show disease-associated change. Cap rather than discard.
REFERENCE_CONDITION_LABELS: frozenset[str] = frozenset(
    {"normal", "healthy", "reference", "control", "na", "unspecified", ""}
)
REFERENCE_ONLY_RATING_CAP = "moderate"

# Thresholds mapping a raw Open Targets datatype score to an ordinal rating. These are
# the rubric's only numeric cut-points and they are versioned with it.
OT_SCORE_TO_STRENGTH: list[tuple[float, str]] = [
    (0.60, "strong"),
    (0.30, "moderate"),
    (0.05, "weak"),
]

# ---------------------------------------------------------------------------------
# Modality accessibility
# ---------------------------------------------------------------------------------
# A binder reaches extracellular and cell-surface epitopes. These substrings are matched
# against UniProt subcellular-location strings to screen accessibility BEFORE expensive
# design, per the section 17 selection policy. This is a modality screen, not a target
# allowlist: it encodes a property of the binder modality, not a preference for any gene.

ACCESSIBLE_LOCATION_HINTS: tuple[str, ...] = (
    "secreted",
    "cell membrane",
    "plasma membrane",
    "cell surface",
    "extracellular",
    "membrane raft",
    "apical cell membrane",
    "basolateral cell membrane",
)

# Locations that are membrane-associated but face the cytosol, so a soluble binder
# cannot engage them. Checked after the accessible hints and they override.
INACCESSIBLE_QUALIFIERS: tuple[str, ...] = (
    "peripheral membrane protein",
    "cytoplasmic side",
    "mitochondri",
    "nucle",
    "endoplasmic reticulum lumen",
    "golgi",
    "lysosom",
    "peroxisom",
)

MODALITY_TRACTABILITY_BUCKET: dict[str, str] = {
    "de novo minibinder": "AB",
    "minibinder": "AB",
    "antibody": "AB",
    "biologic": "AB",
    "small molecule": "SM",
    "protac": "PR",
    "oligonucleotide": "OC",
}


# ---------------------------------------------------------------------------------
# Bounds
# ---------------------------------------------------------------------------------

DEFAULT_BUDGETS: dict[str, int] = {
    "max_candidate_pool": 20,
    "depth_shortlist": 3,
    "max_tool_decisions": 30,
    "max_transient_retries": 2,
    "max_paid_design_batches": 1,
    "poll_wall_clock_s": 1800,
}

CENSUS_BOUNDS: dict[str, int] = {
    "max_cells_interactive": 10_000,
    "default_sampling_seed": 20260922,
    "max_genes_per_query": 50,
    "neighborhood_k": 25,
}

OPEN_TARGETS_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
OPEN_TARGETS_TIMEOUT_S = 60

# Optional demo presets. These are INPUT CONVENIENCES ONLY -- they prefill the brief form
# and are never consulted by the pipeline, the agent loop or the ranking code. Removing
# this dict must not change any production behaviour; tests/test_no_allowlist.py asserts
# that no module other than the UI imports it.
DEMO_PRESETS: list[dict[str, str]] = [
    {
        "raw_indication": "atopic dermatitis",
        "desired_effect": "suppress type 2 inflammatory signalling",
        "moa": "antagonise",
        "modality": "de novo minibinder",
    },
    {
        "raw_indication": "non-small cell lung carcinoma",
        "desired_effect": "restore anti-tumour immune recognition",
        "moa": "antagonise",
        "modality": "de novo minibinder",
    },
]


def rating_from_score(score: float | None) -> str:
    """Map a raw association datatype score onto the rubric's ordinal scale."""
    if score is None:
        return "unknown"
    for cut, label in OT_SCORE_TO_STRENGTH:
        if score >= cut:
            return label
    return "weak" if score > 0 else "unknown"


def is_accessible_to_binder(locations: list[str]) -> tuple[bool | None, list[str]]:
    """REMOVED -- this screen was scientifically wrong and must not be revived.

    It matched keywords against free-text subcellular-location strings. On live data it
    passes PDE4A (P27815), PDE4B (Q07343) and PDE4D (Q08499) as binder-accessible, because
    each is annotated "plasma membrane" or "apical cell membrane" while docking on the
    CYTOPLASMIC face. PDE4D additionally has ~122 experimental structures, so a
    structure-coverage heuristic makes it look like an excellent binder target. It is a
    cytosolic enzyme and no soluble binder can reach it.

    Use ``e2b.adapters.structures.assess_accessibility``, which decides on UniProt
    TOPOLOGY -- signal peptide without transmembrane segment (secreted), or an explicit
    ``Extracellular`` topological domain (ectodomain). Location strings are retained there
    as reported context and never decide the verdict.

    Kept as a raising stub rather than deleted so that reviving it is a loud failure
    instead of a silent one.
    """
    raise NotImplementedError(
        "is_accessible_to_binder() was removed: keyword matching on subcellular-location "
        "strings wrongly passes cytosolic PDE4A/PDE4B/PDE4D as binder-accessible. "
        "Use e2b.adapters.structures.assess_accessibility (UniProt topology) instead."
    )
