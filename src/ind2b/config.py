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

# Antibody-tractability buckets graded by how much they actually evidence
# that a *biologic* can engage the target. The clinical buckets are real
# precedent; the location buckets only restate that the protein is on the
# surface, which the topology filter already establishes, so they score low.
AB_BUCKET_SCORES: dict[str, float] = {
    "Approved Drug": 1.00,
    "Advanced Clinical": 0.85,
    "Phase 1 Clinical": 0.70,
    "UniProt loc high conf": 0.35,
    "GO CC high conf": 0.35,
    "UniProt SigP or TMHMM": 0.30,
    "Human Protein Atlas loc": 0.25,
    "UniProt loc med conf": 0.20,
    "GO CC med conf": 0.20,
}

# Saturation points for count-like evidence: the count at which the
# normalised score reaches ~1.0. Counts are compressed logarithmically
# because the difference between 1 and 10 trials matters far more than the
# difference between 500 and 1000.
SATURATION = {
    "known_drugs": 40,
    "trials": 300,
    "chembl_mechanisms": 6,
}

# Multiplicative penalties, applied after the weighted sum.
#
# The safety penalty is deliberately FLAT rather than scaled by the number of
# recorded liabilities. In this dataset the liability count correlates ~0.65
# with the known-drug count: Open Targets accumulates liability annotations
# for targets that have been taken into patients, so a count-scaled penalty
# measures how well studied a target is, not how dangerous it is. Scaling it
# dropped EGFR - the top-associated NSCLC target, with 82 drugs and 1,735
# trials - twenty places. Severity is not available from this field, so the
# honest treatment is a mild flag plus the count reported in the table for
# the reader to judge.
DEFAULT_PENALTIES: dict[str, float] = {
    "safety_liability_any": 0.92,  # applied once if any liability is recorded
    "broad_expression": 0.90,      # low tissue specificity
}

# Minimum length (residues) of an annotated extracellular topological domain
# for a target to count as reachable by a protein binder. Guards against
# single-pass proteins whose "extracellular" annotation is a few residues of
# linker rather than a foldable, targetable domain.
MIN_ECTODOMAIN_SPAN = 30

# Location keywords used only as a secondary, recorded hint. Accessibility is
# decided from UniProt topology (see sources/uniprot.py) because a location
# string cannot distinguish which side of the membrane a protein sits on:
# KRAS is annotated "Cell membrane" yet faces the cytoplasm.
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
