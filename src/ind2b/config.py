"""Paths, endpoints and run configuration.

Every tunable that affects a reported number lives here or in a weights JSON
file, never inline in the stage code, so a reviewer can see what was assumed.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any

# --- endpoints -------------------------------------------------------------
# All public; no authentication required except OpenAlex (see literature.py).
OPEN_TARGETS_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA_GRAPHQL = "https://data.rcsb.org/graphql"
RCSB_FILES = "https://files.rcsb.org/download"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
CLINICALTRIALS_V2 = "https://clinicaltrials.gov/api/v2/studies"
OPENALEX_WORKS = "https://api.openalex.org/works"

# A plain descriptive agent string. Do not impersonate a browser.
USER_AGENT = "indication2binder/0.1 (python-requests)"

# --- filesystem ------------------------------------------------------------
CACHE_DIR = Path(os.environ.get("IND2B_CACHE", ".ind2b_cache"))
RUNS_DIR = Path(os.environ.get("IND2B_RUNS", "runs"))

STAGE_FILES = {
    0: "stage0_disease.json",
    1: "stage1_evidence.json",
    2: "stage2_ranked_targets.csv",
    3: "stage3_complexes.json",
    4: "stage4_epitopes.csv",
    5: "stage5_specs",  # a directory
}


@dataclass(frozen=True)
class StructureGate:
    """Stage 3 admission rules for an experimental complex.

    ``max_resolution`` is in angstrom and applies to diffraction and cryo-EM
    entries; entries with no reported resolution (e.g. solution NMR) are
    admitted only when ``allow_missing_resolution`` is true.
    """

    max_resolution: float = 4.0
    allow_missing_resolution: bool = False
    min_target_chain_residues: int = 40
    accepted_methods: tuple[str, ...] = (
        "X-RAY DIFFRACTION",
        "ELECTRON MICROSCOPY",
        "ELECTRON CRYSTALLOGRAPHY",
        "SOLUTION NMR",
        "NEUTRON DIFFRACTION",
    )


@dataclass(frozen=True)
class InterfaceParams:
    """Stage 4 interface-definition parameters.

    A residue is called *interface* when it either makes a heavy-atom contact
    with the partner chain within ``contact_cutoff`` angstrom, or buries at
    least ``min_delta_sasa`` A^2 of solvent-accessible surface on complex
    formation. Both criteria are recorded separately in the output so the
    definition used can be audited.
    """

    contact_cutoff: float = 5.0
    min_delta_sasa: float = 10.0
    sasa_probe_radius: float = 1.40
    sasa_n_points: int = 200
    n_hotspots: int = 6
    ignore_hetero: bool = True


@dataclass(frozen=True)
class BinderSpec:
    """Stage 5 BindCraft2 specification defaults."""

    binder_length_min: int = 55
    binder_length_max: int = 120
    n_designs: int = 10


@dataclass
class RunConfig:
    indication: str
    efo_id: str | None = None
    n_targets: int = 100
    include_descendants: bool = True
    weights_file: str | None = None
    gate: StructureGate = field(default_factory=StructureGate)
    interface: InterfaceParams = field(default_factory=InterfaceParams)
    binder: BinderSpec = field(default_factory=BinderSpec)

    def run_dir(self, root: Path | None = None) -> Path:
        root = root or RUNS_DIR
        slug = (self.efo_id or self.indication).replace(" ", "_").replace("/", "-")
        return root / slug

    def to_json(self) -> dict[str, Any]:
        return asdict(self)


# --- scoring ---------------------------------------------------------------
# Weights are the judgement calls in this pipeline. They are data, not code:
# override with --weights to re-rank without touching the source.
DEFAULT_WEIGHTS: dict[str, float] = {
    "ot_overall": 0.22,
    "genetic_evidence": 0.20,
    "known_drug": 0.14,
    "literature": 0.10,
    "expression_specificity": 0.10,
    "pathway_and_model": 0.06,
    "clinical_precedent": 0.08,
    "antibody_tractability": 0.10,
}

# Multiplicative penalties, applied after the weighted sum.
DEFAULT_PENALTIES: dict[str, float] = {
    "safety_liability": 0.85,   # per Open Targets safety liability entry
    "safety_floor": 0.60,       # penalty cannot push below this multiplier
    "broad_expression": 0.90,   # low tissue specificity
}

# Subcellular locations implying a de novo protein binder can physically
# reach the target from outside the cell.
SURFACE_TERMS = (
    "cell membrane",
    "plasma membrane",
    "cell surface",
    "secreted",
    "extracellular",
    "membrane raft",
    "apical plasma membrane",
    "basolateral plasma membrane",
)


def load_weights(path: str | Path | None) -> dict[str, float]:
    """Load scoring weights, falling back to ``DEFAULT_WEIGHTS``."""
    if path is None:
        return dict(DEFAULT_WEIGHTS)
    data = json.loads(Path(path).read_text())
    unknown = set(data) - set(DEFAULT_WEIGHTS)
    if unknown:
        raise ValueError(f"unknown weight keys: {sorted(unknown)}")
    merged = dict(DEFAULT_WEIGHTS)
    merged.update(data)
    return merged
