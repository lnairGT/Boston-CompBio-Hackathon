"""Deterministic, versioned target comparison.

The division of labour the spec insists on: **this module computes, the language model
explains.** Nothing here consults an LLM, and the agent loop is not permitted to invent
weights — it reads ``RUBRIC_VERSION`` and the criteria out of ``e2b.config`` and describes
what the arithmetic did.

Three properties that tests/test_ranking.py pins down:

* Determinism. Identical inputs give an identical order, independent of dict insertion
  order, because the final sort key is ``(-score, target_id)``.
* Missing evidence is never favourable. ``unknown`` contributes 0.0 and additionally
  records *what was looked for and not found*, so an unevidenced target cannot outrank an
  evidenced one on the strength of silence.
* Biological support and design feasibility never mix. Feasibility is computed, reported
  and gates the shortlist, but it contributes nothing to the biology score — so a target
  is never made to look biologically better because it happens to be druggable.
"""

from __future__ import annotations

from typing import Any

from ..config import (
    CRITERIA,
    LIABILITY_PENALTY_CAP,
    LIABILITY_PENALTY_PER_ITEM,
    LITERATURE_ONLY_RATING_CAP,
    REFERENCE_CONDITION_LABELS,
    REFERENCE_ONLY_RATING_CAP,
    RUBRIC_VERSION,
    STRENGTH_POINTS,
    TEXT_MINED_DATATYPES,
    rating_from_score,
)
from ..contracts import (
    Assessment,
    CellSummary,
    ComparisonResult,
    CriterionRating,
    DesignFeasibility,
    Evidence,
    ResearchBrief,
    Target,
)

# Below this, an "extracellular" region is too small to present a useful epitope for a
# designed binder. HRH1 for instance is annotated with an extracellular topological
# domain of 29 residues (a 7TM receptor N-terminus): real, but not a design surface.
MIN_EPITOPE_RESIDUES = 40

# Action types that count as evidence FOR each requested direction of modulation.
DIRECTION_ACTIONS: dict[str, set[str]] = {
    "antagonise": {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE MODULATOR", "INVERSE AGONIST", "DEGRADER"},
    "inhibit": {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE MODULATOR", "DEGRADER"},
    "block": {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE MODULATOR"},
    "agonise": {"AGONIST", "POSITIVE MODULATOR", "ACTIVATOR", "PARTIAL AGONIST"},
    "activate": {"AGONIST", "POSITIVE MODULATOR", "ACTIVATOR"},
}

# Enumerated from the live release: Open Targets 26.06 emits SCREAMING_SNAKE stage values
# (APPROVAL, PHASE_3, ...), not the "Phase III" title-case strings older examples show.
# Both spellings are accepted so the mapping survives a release that changes presentation.
CLINICAL_STAGE_RANK = {
    "APPROVAL": 4, "PHASE_4": 4, "PHASE_3": 3, "PHASE_2": 2, "PHASE_1": 1,
    "PRECLINICAL": 0, "UNKNOWN": 0,
    "Approved": 4, "Phase IV": 4, "Phase III": 3, "Phase II": 2, "Phase I": 1,
    "Preclinical": 0, "Early Phase I": 1,
}


def _normalise_direction(moa: str | None) -> str | None:
    if not moa:
        return None
    m = moa.strip().lower()
    for key in DIRECTION_ACTIONS:
        if key in m:
            return key
    if "antagon" in m or "inhib" in m or "block" in m or "suppress" in m or "neutralis" in m or "neutraliz" in m:
        return "antagonise"
    if "agon" in m or "activat" in m or "restor" in m or "enhanc" in m:
        return "agonise"
    return None


# ---------------------------------------------------------------------------------
# Feasibility
# ---------------------------------------------------------------------------------


def build_feasibility(structure_assessment: dict[str, Any] | None, modality: str | None) -> DesignFeasibility:
    """Design feasibility for the requested modality. Deliberately independent of biology."""
    if structure_assessment is None:
        return DesignFeasibility(
            rating="unknown",
            notes="No structure or accessibility assessment was retrieved for this target.",
        )

    access = structure_assessment.get("accessibility") or {}
    accessible = access.get("accessible")
    regions = access.get("engageable_regions") or []
    experimental = structure_assessment.get("experimental_structures") or []
    af = structure_assessment.get("alphafold")

    blocking: list[str] = []
    largest = 0
    for r in regions:
        s, e = r.get("start"), r.get("end")
        if isinstance(s, int) and isinstance(e, int):
            largest = max(largest, e - s + 1)

    if accessible is None:
        rating = "unknown"
        blocking.append("No subcellular topology annotation retrieved; accessibility unknown.")
    elif not accessible:
        rating = "weak"
        blocking.append(
            "Not accessible to a soluble binder: "
            + "; ".join(access.get("basis") or ["no extracellular topology"])
        )
    else:
        if largest and largest < MIN_EPITOPE_RESIDUES:
            rating = "weak"
            blocking.append(
                f"Largest engageable region is only {largest} residues "
                f"(minimum useful epitope surface taken as {MIN_EPITOPE_RESIDUES}); "
                "extracellular but not a practical design surface."
            )
        elif experimental:
            rating = "strong"
        elif af:
            rating = "moderate"
            blocking.append("No experimental structure; a predicted model would be the design input.")
        else:
            rating = "weak"
            blocking.append("Accessible, but no experimental or predicted structure available.")

    if not experimental and not af:
        blocking.append("No structure of any kind cross-referenced for this accession.")

    n_glyc = access.get("n_glycosylation_sites") or 0
    if n_glyc:
        blocking.append(f"{n_glyc} annotated glycosylation site(s) may occlude an epitope.")

    return DesignFeasibility(
        accessible_to_modality=accessible,
        localization_basis=list(access.get("basis") or []),
        structure_available=bool(experimental or af),
        structure_ids=[s["pdb_id"] for s in experimental[:6] if s.get("pdb_id")]
        + ([af["model_id"]] if af and af.get("model_id") else []),
        blocking_reasons=blocking,
        rating=rating,
        notes=(
            f"Modality '{modality or 'unspecified'}'. Largest engageable region "
            f"{largest or 'unknown'} residues. "
            + (access.get("caveat") or "")
        ),
    )


# ---------------------------------------------------------------------------------
# Assessment
# ---------------------------------------------------------------------------------


def build_assessment(
    *,
    brief: ResearchBrief,
    target: Target,
    association: dict[str, Any],
    ot_annotation: dict[str, Any] | None,
    evidences: list[Evidence],
    structure_assessment: dict[str, Any] | None,
    cell_summaries: list[CellSummary] | None = None,
    rubric_version: str = RUBRIC_VERSION,
) -> Assessment:
    """Turn retrieved evidence into a rated Assessment. Pure function of its inputs."""
    datatype_scores: dict[str, float] = association.get("datatype_scores") or {}
    by_type: dict[str, list[str]] = {}
    for ev in evidences:
        by_type.setdefault(ev.evidence_type, []).append(ev.evidence_id)

    ratings: list[CriterionRating] = []

    for crit in CRITERIA:
        key = crit["key"]

        if key == "direction_fit":
            ratings.append(_rate_direction(brief, ot_annotation, evidences))
            continue

        if key == "cellular_context":
            ratings.append(_rate_cellular(target, cell_summaries, datatype_scores, by_type))
            continue

        # Criteria backed by Open Targets datatype scores.
        relevant = [datatype_scores.get(dt) for dt in crit["ot_datatypes"] if datatype_scores.get(dt) is not None]
        ev_ids = [eid for dt in crit["ot_datatypes"] for eid in by_type.get(dt, [])]
        if not relevant:
            ratings.append(
                CriterionRating(
                    criterion=key,
                    rating="unknown",
                    evidence_ids=[],
                    missingness=(
                        f"No Open Targets evidence of type(s) {crit['ot_datatypes']} was returned for this "
                        f"target-disease pair in the queried release."
                    ),
                    rationale=None,
                )
            )
            continue

        best = max(relevant)
        rating = rating_from_score(best)
        caveats: list[str] = []

        contributing = {
            dt for dt in crit["ot_datatypes"]
            if datatype_scores.get(dt) is not None and datatype_scores[dt] == best
        }
        supporting = {dt for dt in crit["ot_datatypes"] if datatype_scores.get(dt)}
        if supporting and supporting <= TEXT_MINED_DATATYPES and rating == "strong":
            # Text-mined co-occurrence tracks research attention, not causal support.
            rating = LITERATURE_ONLY_RATING_CAP
            caveats.append(
                f"Capped at '{LITERATURE_ONLY_RATING_CAP}': the only contributing datatype is "
                f"text-mined literature, whose score scales with publication volume."
            )

        if rating in ("strong", "moderate") and not ev_ids:
            # An aggregate datatype score exists but no individual evidence rows were
            # retrieved for it. Do not award a favourable rating we cannot point at.
            rating = "weak"
            caveats.append("Downgraded: aggregate datatype score present but no citable evidence rows retrieved.")

        ratings.append(
            CriterionRating(
                criterion=key,
                rating=rating,
                evidence_ids=sorted(set(ev_ids))[:25],
                missingness=None if ev_ids else "Aggregate datatype score present but no individual evidence rows retrieved.",
                rationale=(
                    f"Highest contributing Open Targets datatype score {best:.3f} from "
                    f"{sorted(contributing) or crit['ot_datatypes']}; "
                    f"{len(set(ev_ids))} evidence item(s) retrieved."
                    + ("  " + " ".join(caveats) if caveats else "")
                ),
            )
        )

    liabilities = [
        f"{s.get('event')} ({s.get('datasource')})"
        for s in ((ot_annotation or {}).get("safety_liabilities") or [])
        if s.get("event")
    ]

    feasibility = build_feasibility(structure_assessment, brief.modality)
    direction_rating = next(r for r in ratings if r.criterion == "direction_fit")

    excluded = False
    reason: str | None = None
    if feasibility.accessible_to_modality is False:
        excluded = True
        reason = (
            f"Modality mismatch: {brief.modality or 'binder'} cannot engage this target. "
            + "; ".join(feasibility.blocking_reasons[:2])
        )
    elif feasibility.rating == "weak" and feasibility.accessible_to_modality and any(
        "practical design surface" in b for b in feasibility.blocking_reasons
    ):
        excluded = True
        reason = "Extracellular region too small to present a useful epitope: " + "; ".join(
            b for b in feasibility.blocking_reasons if "practical design surface" in b
        )

    return Assessment(
        target_id=target.target_id,
        brief_id=brief.brief_id,
        criteria=ratings,
        disease_rationale=(
            f"Open Targets overall association {association.get('overall_score')} "
            f"(rank {association.get('rank')} of the retrieved pool) with "
            f"{brief.resolved_disease_label or brief.resolved_disease_id}."
        ),
        moa_rationale=direction_rating.rationale,
        direction_compatibility=(
            "compatible" if direction_rating.rating in ("strong", "moderate")
            else "unknown" if direction_rating.rating == "unknown"
            else "incompatible" if direction_rating.missingness and "opposite" in direction_rating.missingness
            else "unknown"
        ),
        liabilities=liabilities,
        design_feasibility=feasibility,
        rubric_version=rubric_version,
        review_status="draft",
        excluded=excluded,
        exclusion_reason=reason,
    )


def _rate_direction(
    brief: ResearchBrief, ot_annotation: dict[str, Any] | None, evidences: list[Evidence]
) -> CriterionRating:
    """Direction of modulation must be verified separately. Association never establishes it."""
    requested = _normalise_direction(brief.moa) or _normalise_direction(brief.desired_effect)
    if requested is None:
        return CriterionRating(
            criterion="direction_fit",
            rating="unknown",
            missingness="The brief does not state an intervention direction, so fit cannot be assessed.",
        )

    wanted = DIRECTION_ACTIONS.get(requested, set())
    opposite: set[str] = set()
    for k, v in DIRECTION_ACTIONS.items():
        if k != requested and not (v & wanted):
            opposite |= v

    # A drug's precedent counts far more when it was pursued FOR THIS INDICATION. Section 17:
    # "A drug's success for one indication, target combination or modality does not
    # automatically transfer to this request." So matches are partitioned, not pooled.
    this_disease = {brief.resolved_disease_id, (brief.resolved_disease_label or "").lower()} - {None, ""}

    on_indication: list[tuple[str, str, str | None]] = []
    off_indication: list[tuple[str, str, str | None, str]] = []
    conflicting: list[tuple[str, str]] = []
    for drug in ((ot_annotation or {}).get("known_drugs") or []):
        stage = drug.get("max_clinical_stage")
        labels = {str(x).lower() for x in (drug.get("indications") or []) if x}
        ids = {x for x in (drug.get("indication_ids") or []) if x}
        matches_indication = bool((labels | ids) & this_disease)
        for m in drug.get("mechanisms") or []:
            action = (m.get("action_type") or "").upper()
            name = drug.get("drug_name") or "?"
            if action in wanted:
                if matches_indication:
                    on_indication.append((name, action, stage))
                else:
                    off_indication.append((name, action, stage, ", ".join(sorted(labels)[:2]) or "unstated"))
                break
            if action in opposite:
                conflicting.append((name, action))
                break

    modulation = {
        (ev.measurement or {}).get("target_modulation")
        for ev in evidences
        if (ev.measurement or {}).get("target_modulation")
    }

    ev_ids = [
        ev.evidence_id for ev in evidences if ev.evidence_type in ("known_drug", "clinical")
    ][:20]

    if on_indication or off_indication:
        best_on = max((CLINICAL_STAGE_RANK.get(s or "", 0) for _, _, s in on_indication), default=-1)
        if best_on >= 3:
            rating = "strong"
        elif best_on >= 1:
            rating = "moderate"
        elif on_indication:
            rating = "weak"
        else:
            # Direction precedent exists, but only in other indications.
            rating = "weak"

        if rating in ("strong", "moderate") and not ev_ids:
            rating = "weak"

        on_txt = ", ".join(f"{n} ({a}, {s})" for n, a, s in on_indication[:4]) or "none"
        off_txt = ", ".join(f"{n} ({a}, {s}, for {i})" for n, a, s, i in off_indication[:3]) or "none"
        return CriterionRating(
            criterion="direction_fit",
            rating=rating,
            evidence_ids=ev_ids,
            missingness=(
                None
                if on_indication
                else (
                    f"Direction precedent exists only in OTHER indications ({off_txt}); no agent with a "
                    f"'{requested}'-compatible mechanism has been pursued for "
                    f"{brief.resolved_disease_label or brief.resolved_disease_id}. Precedent does not transfer."
                )
            ),
            rationale=(
                f"Requested direction '{requested}'. Agents matching that direction FOR THIS INDICATION: "
                f"{on_txt}. Matching agents for other indications: {off_txt}. "
                f"A trial entry is not evidence of efficacy or approval, and a result in one indication "
                f"does not transfer to this request."
            ),
        )

    if conflicting:
        return CriterionRating(
            criterion="direction_fit",
            rating="unknown",
            evidence_ids=ev_ids,
            missingness=(
                f"Only agents with the opposite direction are annotated "
                f"({', '.join(f'{n} ({a})' for n, a in conflicting[:3])}); no evidence that "
                f"'{requested}' is the beneficial direction."
            ),
            rationale=None,
        )

    return CriterionRating(
        criterion="direction_fit",
        rating="unknown",
        evidence_ids=[],
        missingness=(
            f"No annotated agent with a '{requested}'-compatible mechanism, and no target_modulation "
            f"direction in the retrieved evidence"
            + (f" (observed modulation values: {sorted(modulation)})" if modulation else "")
            + ". Direction of benefit is unresolved."
        ),
    )


def _rate_cellular(
    target: Target,
    cell_summaries: list[CellSummary] | None,
    datatype_scores: dict[str, float],
    by_type: dict[str, list[str]],
) -> CriterionRating:
    """Census-derived cellular context. Falls back to 'unknown', never to the bulk RNA score."""
    mine = [c for c in (cell_summaries or []) if c.target_id == target.target_id]
    if not mine:
        return CriterionRating(
            criterion="cellular_context",
            rating="unknown",
            evidence_ids=[],
            missingness=(
                "No CELLxGENE Census cell summary was retrieved for this target in the selected slice. "
                "Note the Open Targets rna_expression datatype score is bulk-tissue aggregate evidence "
                "and is deliberately NOT substituted for single-cell disease context here."
            ),
        )

    donor_counts = [c.n_donors for c in mine if c.n_donors is not None]
    measured = [v for c in mine for v in c.values.values() if v is not None]
    if not measured:
        return CriterionRating(
            criterion="cellular_context",
            rating="unknown",
            evidence_ids=[c.summary_id for c in mine],
            missingness="A cell summary exists but contains no measurable expression values.",
        )

    peak = max(measured)
    n_donors = max(donor_counts) if donor_counts else 0
    # Single-donor observations cannot support a 'strong' rating however high the value.
    if n_donors < 2:
        rating = "weak"
    elif peak >= 0.25:
        rating = "strong"
    elif peak >= 0.10:
        rating = "moderate"
    elif peak > 0:
        rating = "weak"
    else:
        rating = "unknown"

    # Baseline expression in healthy tissue is NOT disease context. It says the target is
    # present in the mechanism-relevant populations; it says nothing about whether the
    # disease changes it. Many indications simply have no diseased tissue in the Census --
    # atopic dermatitis has none in skin -- and in that case the honest ceiling is
    # 'moderate', with the gap named. Silently rating baseline data as disease context
    # would let an absent contrast masquerade as a positive finding.
    conditions = {(c.condition or "").strip().lower() for c in mine}
    disease_matched = any(cond not in REFERENCE_CONDITION_LABELS for cond in conditions)
    reference_only_note = None
    if not disease_matched and rating == "strong":
        rating = REFERENCE_ONLY_RATING_CAP
        reference_only_note = (
            f"Capped at '{REFERENCE_ONLY_RATING_CAP}': every retrieved cell summary is from "
            f"reference/normal tissue ({sorted(conditions)}), so this measures baseline presence in "
            f"the relevant populations, not disease-associated change. No diseased-tissue slice was "
            f"available for this indication."
        )

    return CriterionRating(
        criterion="cellular_context",
        rating=rating,
        evidence_ids=[c.summary_id for c in mine],
        missingness=None if rating != "unknown" else "Expression measured as zero in all sampled populations.",
        rationale=(
            f"Peak {mine[0].expression_metric} of {peak:.3f} across {len(mine)} population summary(ies), "
            f"{n_donors} donor(s), condition(s) {sorted(conditions)}. "
            + ("Single donor, so rated no higher than weak. " if n_donors < 2 else "")
            + (reference_only_note + " " if reference_only_note else "")
            + "RNA expression is not surface protein abundance, and undetected expression is not "
              "proof of absence."
        ),
    )


# ---------------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------------


def score_assessment(assessment: Assessment) -> tuple[float, dict[str, float]]:
    """Weighted sum over the versioned criteria. Liability penalty capped."""
    parts: dict[str, float] = {}
    total = 0.0
    weights = {c["key"]: c["weight"] for c in CRITERIA}
    for rating in sorted(assessment.criteria, key=lambda r: r.criterion):
        w = weights.get(rating.criterion)
        if w is None:
            continue
        contribution = STRENGTH_POINTS[rating.rating] * w
        parts[rating.criterion] = round(contribution, 6)
        total += contribution
    penalty = min(len(assessment.liabilities) * LIABILITY_PENALTY_PER_ITEM, LIABILITY_PENALTY_CAP)
    parts["liability_penalty"] = -round(penalty, 6)
    return round(max(total - penalty, 0.0), 6), parts


def compare_targets(
    assessments: list[Assessment],
    brief_id: str,
    *,
    rubric_version: str = RUBRIC_VERSION,
    shortlist: int = 3,
    pool_scope: dict[str, Any] | None = None,
    tie_epsilon: float = 1e-9,
) -> ComparisonResult:
    """Deterministic comparison. Excluded targets keep their biology visible but do not rank."""
    excluded = {a.target_id: (a.exclusion_reason or "excluded") for a in assessments if a.excluded}

    # A target whose accessibility has never been screened is NOT a candidate. Without this
    # gate an unscreened target rides into the shortlist on the strength of having no adverse
    # evidence yet -- which is precisely missing evidence becoming a favourable position.
    pending: dict[str, str] = {}
    included: list[Assessment] = []
    for a in assessments:
        if a.excluded:
            continue
        if a.design_feasibility.accessible_to_modality is None:
            pending[a.target_id] = (
                "Accessibility to the requested modality has not been screened "
                "(assess_structure was not called for this target), so it is not eligible "
                "for the shortlist."
            )
            continue
        included.append(a)

    scored: list[tuple[float, str, Assessment]] = []
    for a in included:
        s, _ = score_assessment(a)
        scored.append((s, a.target_id, a))
    # Stable, reproducible: descending score, then ascending target_id.
    scored.sort(key=lambda t: (-t[0], t[1]))

    ordered = [tid for _, tid, _ in scored]
    scores = {tid: s for s, tid, _ in scored}

    ties: list[list[str]] = []
    group: list[str] = []
    for i, (s, tid, _) in enumerate(scored):
        if i and abs(s - scored[i - 1][0]) <= tie_epsilon:
            if not group:
                group = [scored[i - 1][1]]
            group.append(tid)
        elif group:
            ties.append(group)
            group = []
    if group:
        ties.append(group)

    unresolved: dict[str, list[str]] = {}
    for a in assessments:
        gaps = [r.criterion for r in a.criteria if r.rating == "unknown"]
        if gaps:
            unresolved[a.target_id] = sorted(gaps)

    top = scored[:shortlist]
    insufficient = (not included) or all(s == 0.0 for s, _, _ in top)

    if not included:
        explanation = (
            f"No screened target in the retrieved pool survived the modality screen. "
            f"{len(excluded)} target(s) were excluded: "
            + "; ".join(f"{k} ({v[:90]})" for k, v in list(excluded.items())[:5])
            + (f". A further {len(pending)} target(s) have not been accessibility-screened yet "
               f"and are not eligible; screen them before concluding no target is feasible."
               if pending else "")
        )
    elif insufficient:
        explanation = (
            "Every shortlisted target scored 0.0 under "
            f"{rubric_version}: retrieved evidence was insufficient to support any ranking."
        )
    else:
        lead = top[0]
        parts = score_assessment(lead[2])[1]
        drivers = ", ".join(
            f"{k}={v:+.3f}" for k, v in sorted(parts.items(), key=lambda kv: -abs(kv[1])) if abs(v) > 1e-9
        )
        explanation = (
            f"Under {rubric_version}, {lead[1]} leads with {lead[0]:.3f} ({drivers}). "
            f"{len(included)} of {len(assessments)} assessed targets entered the ranking; "
            f"{len(excluded)} were excluded on modality or accessibility grounds with recorded reasons"
            + (f"; {len(pending)} have not been accessibility-screened yet and are therefore not "
               f"eligible for the shortlist" if pending else "")
            + ". Design feasibility is reported separately and contributes nothing to this score."
        )
        if ties:
            explanation += f" Unbroken tie(s) present: {ties}."

    return ComparisonResult(
        brief_id=brief_id,
        rubric_version=rubric_version,
        ordered_target_ids=ordered,
        scores=scores,
        ties=ties,
        excluded=excluded,
        pending_screen=pending,
        unresolved_criteria=unresolved,
        insufficient_evidence=insufficient,
        explanation=explanation,
        candidate_pool_scope=pool_scope or {},
    )
