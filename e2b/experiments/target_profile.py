"""Read the target's biology into the two facts that decide the assay formats.

Which assays are even possible is decided by *topology*, not by the indication:

* a **secreted ligand** is a soluble analyte -- it can be titrated in solution, and
  blocking it is a ligand-trap experiment run on cells that express its receptor;
* a **cell-surface receptor** has to be produced as an ectodomain construct, which
  raises glycosylation and oligomeric-state questions, and blocking it is run on
  cells that express the receptor itself;
* a target with **no extracellular face** cannot be engaged by a soluble binder at
  all, so no in-cell blockade experiment can honestly be designed for this modality.

Lane 1 already decides accessibility from UniProt topology in
``e2b/adapters/structures.py``. This module reuses that verdict when it is handed in
or when the adapter is importable, and never re-derives accessibility from location
keyword strings -- the PDE4 case in that adapter's docstring is exactly why.
"""

from __future__ import annotations

from typing import Any

from .coerce import as_dict, get
from .moa import classify_mechanism
from .records import TargetAssayProfile

# UniProt topology verdict -> assay archetype.
_VERDICT_TO_ARCHETYPE = {
    "accessible_secreted": "secreted_ligand",
    "accessible_ectodomain": "cell_surface_receptor",
    "not_accessible": "not_accessible_to_soluble_binder",
}


class TargetBiologyUnavailable(RuntimeError):
    """Raised when target biology could not be retrieved.

    This is a *source failure*, not evidence that the target has no extracellular
    face. The caller must report it as such.
    """

    def __init__(self, message: str, *, transient: bool = False):
        super().__init__(message)
        self.transient = transient


def fetch_target_biology(accession: str) -> tuple[dict[str, Any], dict[str, Any], str]:
    """Fetch UniProt protein + accessibility via the Lane 1 adapter (read-only use).

    Returns ``(protein, accessibility, origin)``. Raises ``TargetBiologyUnavailable``
    on any retrieval failure, so a timeout is never silently turned into "no
    extracellular domain".
    """
    try:
        from ..adapters import structures as _structures
    except Exception as exc:  # pragma: no cover - import-time environment problem
        raise TargetBiologyUnavailable(
            f"structure adapter not importable: {exc}", transient=False
        ) from exc
    try:
        protein, _prov = _structures.fetch_protein(accession)
    except Exception as exc:
        transient = bool(getattr(exc, "transient", False))
        raise TargetBiologyUnavailable(
            f"UniProt retrieval for {accession} failed: {exc}", transient=transient
        ) from exc
    if is_empty_protein(protein):
        # UniProt answers a malformed or unknown accession with HTTP 200 and an empty
        # body rather than a 404, so the adapter returns a record with length 0 and no
        # features. An empty record has no topology, and a topology screen run on it
        # would return "no extracellular face" -- which would be a fabricated negative.
        raise TargetBiologyUnavailable(
            f"UniProt returned an empty record for '{accession}' (no sequence, no features). "
            "The accession is unknown or malformed. This is a retrieval/identity failure and "
            "must not be read as the target having no extracellular domain.",
            transient=False,
        )
    accessibility = _structures.assess_accessibility(protein)
    return protein, accessibility, "real"


def is_empty_protein(protein: dict[str, Any] | None) -> bool:
    """True when a protein record carries no sequence, so no topology can be read."""
    if not protein:
        return True
    return not (protein.get("sequence") or protein.get("length"))


def build_target_profile(
    target: Any,
    *,
    protein: dict[str, Any] | None = None,
    accessibility: dict[str, Any] | None = None,
    pathway_hints: list[str] | None = None,
    desired_effect: str | None = None,
    fetch: bool = False,
    origin: str = "real",
) -> tuple[TargetAssayProfile, dict[str, Any]]:
    """Assemble the assay-relevant profile of the target, and the mechanism-class record.

    ``protein``/``accessibility`` are the dicts produced by
    ``e2b.adapters.structures``. When they are absent and ``fetch`` is true, they are
    retrieved live from the target's accession. When they are absent and ``fetch`` is
    false, the archetype stays ``unknown`` -- which the planner reports as a blocking
    missing input rather than guessing.
    """
    t = as_dict(target)
    accession = get(t, "protein_accession", "accession", "uniprot")
    notes: list[str] = []

    if protein is None and fetch and accession:
        protein, accessibility, origin = fetch_target_biology(accession)
        notes.append(f"UniProt entry retrieved live for {accession}")

    if protein is not None and is_empty_protein(protein):
        raise TargetBiologyUnavailable(
            "The supplied UniProt record is empty (no sequence, no features), so no topology "
            "can be read from it. Treating this as 'no extracellular domain' would be a "
            "fabricated negative.",
            transient=False,
        )

    protein = protein or {}
    if accessibility is None and protein:
        try:
            from ..adapters import structures as _structures

            accessibility = _structures.assess_accessibility(protein)
            notes.append("accessibility verdict computed from supplied UniProt topology")
        except Exception:
            accessibility = None
    accessibility = accessibility or {}

    verdict = accessibility.get("verdict")
    archetype = _VERDICT_TO_ARCHETYPE.get(verdict or "", "unknown")
    basis = list(accessibility.get("basis") or [])
    if archetype == "cell_surface_receptor":
        regions = accessibility.get("engageable_regions") or []
        if any(r.get("kind") == "inferred_ectodomain" for r in regions):
            archetype = "cell_surface_receptor_inferred"
            basis.append(
                "ectodomain boundaries were inferred from signal peptide and first "
                "transmembrane segment, not read from an explicit 'Extracellular' "
                "topological domain: treat the construct boundary as provisional"
            )
    if archetype == "unknown":
        basis.append(
            "no UniProt topology verdict was available, so the assay archetype is "
            "undetermined; this is missing information, not a finding"
        )

    mech = classify_mechanism(
        protein_name=protein.get("protein_name") or get(t, "approved_name"),
        gene_symbol=get(t, "gene_symbol"),
        keywords=protein.get("keywords") or [],
        pathway_hints=pathway_hints or [],
        desired_effect=desired_effect,
        extra_text=[str(x) for x in (get(t, "subcellular_locations", default=[]) or [])],
    )

    pdbs = protein.get("pdb_entries")
    return TargetAssayProfile(
        gene_symbol=get(t, "gene_symbol"),
        protein_accession=accession,
        protein_name=protein.get("protein_name") or get(t, "approved_name"),
        archetype=archetype,
        archetype_basis=basis,
        mechanism_class=mech["key"],
        mechanism_class_label=mech.get("label"),
        mechanism_class_basis=list(mech.get("basis") or []) + (
            [mech["missing"]] if mech.get("missing") else []
        ),
        engageable_regions=list(accessibility.get("engageable_regions") or []),
        length_aa=protein.get("length"),
        n_glycosylation_sites=(
            len(protein.get("glycosylation_sites")) if protein.get("glycosylation_sites") is not None else None
        ),
        n_disulfides=(
            len(protein.get("disulfide_bonds")) if protein.get("disulfide_bonds") is not None else None
        ),
        experimental_structures=(len(pdbs) if pdbs is not None else None),
        origin=origin if origin in ("real", "cached_real", "fixture") else "real",
        source_notes=notes + ([accessibility["caveat"]] if accessibility.get("caveat") else []),
    ), mech


def antigen_reagent_plan(profile: TargetAssayProfile, protein: dict[str, Any] | None) -> dict[str, Any]:
    """What the target-side reagent has to be, given the archetype.

    Returns the host, the construct shape and the specific reasons -- this is where a
    secreted cytokine and a receptor ectodomain diverge.
    """
    protein = protein or {}
    n_glyc = profile.n_glycosylation_sites
    n_ss = profile.n_disulfides
    regions = profile.engageable_regions

    if profile.archetype == "secreted_ligand":
        span = next((r for r in regions if r.get("kind") == "mature_secreted"), None)
        return {
            "construct": "mature secreted protein (signal peptide removed)",
            "span": span,
            "host": (
                "mammalian (HEK293-derived) secretory expression"
                if (n_glyc or 0) > 0
                else "E. coli periplasmic or refolded expression, with mammalian expression as the fallback"
            ),
            "host_rationale": (
                f"{n_glyc} N-glycosylation site(s) annotated: express in a mammalian host so the "
                "glycan occupancy resembles the native ligand, because glycans near an epitope "
                "change binder access"
                if (n_glyc or 0) > 0
                else "no annotated glycosylation, so a bacterial route is viable; a mammalian "
                "control preparation is still worth having to rule out a folding artefact"
            ),
            "reasons": [
                f"{n_ss if n_ss is not None else 'unknown number of'} annotated disulfide bond(s): "
                "oxidative folding must be verified, since a reduced or scrambled ligand is not the "
                "analyte you think you are titrating",
                "the ligand is the analyte in the binding assay and the stimulus in the functional "
                "assay, so the same lot should be used for both and its bioactivity established "
                "against a reference preparation",
            ],
            "activity_requirement": (
                "the recombinant ligand lot must produce a dose-dependent pathway response in the "
                "functional assay before any binding number from it is believed"
            ),
        }

    if profile.archetype.startswith("cell_surface_receptor"):
        span = next((r for r in regions if r.get("kind") in ("ectodomain", "inferred_ectodomain")), None)
        return {
            "construct": "ectodomain construct, boundaries taken from the topology annotation",
            "span": span,
            "host": "mammalian (HEK293-derived) secretory expression",
            "host_rationale": (
                "receptor ectodomains are usually glycosylated and disulfide-bonded "
                f"({n_glyc if n_glyc is not None else 'unknown'} annotated N-glycosylation site(s), "
                f"{n_ss if n_ss is not None else 'unknown'} disulfide bond(s)); a bacterial "
                "preparation risks a non-native fold that a binder may or may not recognise"
            ),
            "reasons": [
                "produce BOTH a monomeric tagged ectodomain and, if an Fc or other dimerising "
                "fusion is used for capture, treat that reagent as bivalent: a bivalent analyte "
                "gives an avidity-inflated apparent affinity, which is not the monovalent affinity",
                "confirm the ectodomain is correctly folded by binding a reference ligand or a "
                "conformation-dependent antibody, not only by SDS-PAGE purity",
            ],
            "activity_requirement": (
                "the ectodomain lot must bind its native ligand at the expected order of magnitude "
                "before it is used to judge a designed binder"
            ),
        }

    return {
        "construct": None,
        "span": None,
        "host": None,
        "host_rationale": None,
        "reasons": [
            "the target has no annotated extracellular face, so there is no antigen reagent that "
            "a soluble binder could engage in its native context"
        ],
        "activity_requirement": None,
    }
