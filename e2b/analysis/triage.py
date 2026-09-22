"""Deterministic, versioned candidate triage.

Given ``contracts.Candidate``-shaped records from the design lane, this module
produces, per candidate: computed sequence-developability indicators, an
interpretation of the design tool's own interface metrics, the paralogue
cross-reactivity risk for the engaged epitope, and a verdict of
``advance`` / ``conditional`` / ``deprioritise`` with the specific flags that
drove it.

Design rules that are enforced in code, not left to discipline:

* **Deterministic.** No randomness, no model calls, no wall-clock in the decision
  path. Same inputs -> same verdict, same ordering, same ``TRIAGE_VERSION``.
* **Missing makes it worse, never better.** Every indicator that could not be
  computed becomes a flag of severity ``unknown``, and an ``unknown`` can only
  hold a candidate back. There is no code path in which absent data improves a
  verdict.
* **No composite score.** Candidates are ordered by an explicit lexicographic key
  over categorical fields (stated in ``RANKING_KEY_DESCRIPTION``), not by a number
  blended from incomparable scales.
* **No measurement language.** Nothing emitted here is a Kd, IC50, EC50 or
  affinity, and no output asserts that a candidate binds.
"""

from __future__ import annotations

from typing import Any

from .developability import as_contract_metric_dicts, compute_developability
from .interface import claim_language, interpret_interface_metrics
from .paralogs import paralog_cross_reactivity_risk

TRIAGE_VERSION = "triage-2026.09.22-v1.1.0"

SEVERITY_ORDER = ("blocker", "major", "unknown", "minor", "info")

# All decision thresholds in one versioned block. Each is a convention; every one
# is defined with its rationale and its limitation in docs/analysis_methods.md.
TRIAGE_THRESHOLDS: dict[str, Any] = {
    "length_min_expected": 40,
    "length_max_expected": 150,
    "length_min_hard": 25,
    "length_max_hard": 220,
    "gravy_minor": 0.0,
    "gravy_major": 0.5,
    "hydrophobic_patch_fraction_minor": 0.10,
    "hydrophobic_patch_fraction_major": 0.25,
    "pi_neutral_window": (6.5, 8.0),
    "n_glyc_minor": 1,
    "n_glyc_major": 3,
    "deamidation_ng_minor": 1,
    "deamidation_ng_major": 3,
    "oxidation_density_minor": 0.06,
    "low_complexity_fraction_minor": 0.20,
    "epitope_identity_high_risk_pct": 70.0,
    "epitope_identity_moderate_risk_pct": 40.0,
    "max_majors_before_deprioritise": 3,
}

RANKING_KEY_DESCRIPTION = (
    "Candidates are ordered lexicographically by: (1) verdict rank "
    "advance < conditional < deprioritise; (2) interface confidence band "
    "high < moderate < low < very_low < unknown; (3) number of blocker flags, "
    "ascending; (4) number of major flags, ascending; (5) number of unknown flags, "
    "ascending; (6) epitope cross-reactivity risk band low < moderate < high < unknown; "
    "(7) candidate_id, ascending, as a deterministic tie-break. No numeric score is "
    "formed across metrics of different scales."
)

_BAND_RANK = {"high": 0, "moderate": 1, "low": 2, "very_low": 3, None: 4}
_RISK_RANK = {"low": 0, "moderate": 1, "high": 2, None: 3}
_VERDICT_RANK = {"advance": 0, "conditional": 1, "deprioritise": 2}


def _flag(code: str, severity: str, message: str, indicator: str | None = None, value: Any = None) -> dict[str, Any]:
    return {"code": code, "severity": severity, "message": message, "indicator": indicator, "value": value}


# ---------------------------------------------------------------------------------
# Flagging
# ---------------------------------------------------------------------------------


def flag_developability(dev: dict[str, Any]) -> list[dict[str, Any]]:
    """Turn computed sequence indicators into severity-tagged flags."""
    T = TRIAGE_THRESHOLDS
    flags: list[dict[str, Any]] = []
    if not dev.get("computed"):
        flags.append(
            _flag(
                "sequence_indicators_unavailable", "unknown",
                "No sequence was available, so no developability indicator could be computed. "
                + "; ".join(dev.get("notes") or []),
                indicator="all_sequence_indicators",
            )
        )
        return flags

    ind = dev["indicators"]

    def val(k: str):
        return ind.get(k, {}).get("value")

    n = val("length")
    if n is not None:
        if n < T["length_min_hard"] or n > T["length_max_hard"]:
            flags.append(_flag("length_out_of_range", "major",
                               f"Length {n} aa is outside the plausible de novo minibinder range "
                               f"{T['length_min_hard']}-{T['length_max_hard']} aa.", "length", n))
        elif n < T["length_min_expected"] or n > T["length_max_expected"]:
            flags.append(_flag("length_atypical", "minor",
                               f"Length {n} aa is outside the typical de novo minibinder range "
                               f"{T['length_min_expected']}-{T['length_max_expected']} aa.", "length", n))

    up = val("unpaired_cysteine")
    cys_n = val("cysteine_count")
    if up == 1:
        flags.append(_flag("unpaired_cysteine", "major",
                           f"At least one free thiol inferred ({cys_n} cysteines, odd count). A free "
                           "cysteine in a secreted binder invites disulfide scrambling, covalent dimer "
                           "and lot-to-lot heterogeneity.", "unpaired_cysteine", cys_n))
    elif up is None and (cys_n or 0) > 0:
        flags.append(_flag("cysteine_pairing_unknown", "unknown",
                           f"{cys_n} cysteines with an even count; pairing cannot be determined from "
                           "sequence alone and was not supplied. Treated as unresolved, not as paired.",
                           "unpaired_cysteine", cys_n))

    g = val("gravy")
    if g is not None:
        if g >= T["gravy_major"]:
            flags.append(_flag("high_hydrophobicity", "major",
                               f"GRAVY {g} indicates a strongly hydrophobic sequence (solubility and "
                               "aggregation risk).", "gravy", g))
        elif g >= T["gravy_minor"]:
            flags.append(_flag("net_hydrophobic", "minor",
                               f"GRAVY {g} is net-hydrophobic for a soluble binder.", "gravy", g))

    ap = val("hydrophobic_patch_fraction")
    if ap is not None:
        if ap >= T["hydrophobic_patch_fraction_major"]:
            flags.append(_flag("large_hydrophobic_patches", "major",
                               f"{ap:.0%} of residues sit in uncharged high-hydropathy windows "
                               "(aggregation-prone-region proxy).", "hydrophobic_patch_fraction", ap))
        elif ap >= T["hydrophobic_patch_fraction_minor"]:
            flags.append(_flag("hydrophobic_patches", "minor",
                               f"{ap:.0%} of residues sit in uncharged high-hydropathy windows.",
                               "hydrophobic_patch_fraction", ap))

    pi = val("theoretical_pi")
    lo, hi = T["pi_neutral_window"]
    if pi is not None and lo <= pi <= hi:
        flags.append(_flag("pi_near_physiological", "minor",
                           f"Theoretical pI {pi} lies in {lo}-{hi}; net charge near zero at formulation "
                           "and serum pH tends to reduce solubility.", "theoretical_pi", pi))

    gl = val("n_glycosylation_sequons")
    if gl is not None:
        if gl >= T["n_glyc_major"]:
            flags.append(_flag("many_n_glyc_sequons", "major",
                               f"{gl} N-X-S/T sequons; heterogeneous glycosylation in a eukaryotic host, "
                               "and sequons at the interface can block engagement.",
                               "n_glycosylation_sequons", gl))
        elif gl >= T["n_glyc_minor"]:
            flags.append(_flag("n_glyc_sequon", "minor",
                               f"{gl} N-X-S/T sequon(s) present (presence, not occupancy).",
                               "n_glycosylation_sequons", gl))

    deam_detail = ind.get("deamidation_motifs", {}).get("detail", {}) or {}
    ng = len(deam_detail.get("ng", []) or [])
    if ng >= T["deamidation_ng_major"]:
        flags.append(_flag("many_deamidation_motifs", "major",
                           f"{ng} NG motifs (fastest deamidation context); chemical instability risk.",
                           "deamidation_motifs", ng))
    elif ng >= T["deamidation_ng_minor"]:
        flags.append(_flag("deamidation_motif", "minor", f"{ng} NG motif(s) present.",
                           "deamidation_motifs", ng))

    ox = val("oxidation_prone_sites")
    if ox is not None and n:
        dens = ox / n
        if dens >= T["oxidation_density_minor"]:
            flags.append(_flag("oxidation_prone_density", "minor",
                               f"{ox} Met/Trp in {n} aa ({dens:.1%}); oxidation-prone side chains. "
                               "Solvent exposure not assessed.", "oxidation_prone_sites", ox))

    lc = val("fraction_low_complexity")
    if lc is not None and lc >= T["low_complexity_fraction_minor"]:
        flags.append(_flag("low_complexity", "minor",
                           f"{lc:.0%} of residues lie in low-complexity windows.",
                           "fraction_low_complexity", lc))
    return flags


def flag_interface(itf: dict[str, Any]) -> list[dict[str, Any]]:
    flags: list[dict[str, Any]] = []
    band = itf.get("confidence_band")
    if band is None:
        flags.append(_flag("interface_confidence_unknown", "unknown",
                           "No recognised interface confidence metric carried a value. Interface support "
                           "is unknown; this is not evidence of a poor interface, and it is not evidence "
                           "of a good one.", "confidence_band", None))
    elif band == "very_low":
        flags.append(_flag("interface_confidence_very_low", "blocker",
                           "Interface confidence is in the lowest band on the weakest recognised metric.",
                           "confidence_band", band))
    elif band == "low":
        flags.append(_flag("interface_confidence_low", "major",
                           "Interface confidence is low on the weakest recognised metric.",
                           "confidence_band", band))
    elif band == "moderate":
        flags.append(_flag("interface_confidence_moderate", "minor",
                           "Interface confidence is moderate on the weakest recognised metric.",
                           "confidence_band", band))
    for gf in itf.get("geometry_flags") or []:
        flags.append(_flag("interface_" + gf["flag"],
                           "major" if gf["severity"] == "major" else "minor",
                           gf["detail"], "geometry", gf.get("value")))
    if itf.get("n_confidence_missing"):
        flags.append(_flag("interface_metric_value_missing", "unknown",
                           f"{itf['n_confidence_missing']} recognised confidence metric(s) were present "
                           "but had no value.", "confidence_metrics", None))
    if itf.get("uninterpreted"):
        names = [u["metric_name"] for u in itf["uninterpreted"]]
        flags.append(_flag("uninterpreted_metrics", "info",
                           f"Metrics reported verbatim but not scored by this triage version: {names}.",
                           "uninterpreted", len(names)))
    return flags


def flag_cross_reactivity(par: dict[str, Any] | None) -> list[dict[str, Any]]:
    flags: list[dict[str, Any]] = []
    if par is None:
        flags.append(_flag("paralog_analysis_not_run", "unknown",
                           "Paralogue cross-reactivity analysis was not run for this candidate.",
                           "epitope_risk_band", None))
        return flags
    status = par.get("status")
    if status == "source_failure":
        flags.append(_flag("paralog_retrieval_failed", "unknown",
                           "Family retrieval failed against the source. A retrieval failure is a source "
                           "failure, not evidence that the target has no paralogues.",
                           "epitope_risk_band", None))
        return flags
    if par.get("epitope_status") != "specified":
        flags.append(_flag("epitope_unspecified", "unknown",
                           "No engaged epitope was supplied, so epitope-level identity could not be "
                           "computed. Cross-reactivity risk is unknown.", "epitope_risk_band", None))
        return flags
    if status == "no_family_members_found":
        flags.append(_flag("no_paralogues_found", "info",
                           "No reviewed human family member was returned by any route that succeeded. "
                           "This is an absence reported by the routes queried, not a proteome-wide "
                           "specificity guarantee.", "epitope_risk_band", None))
        return flags
    band = (par.get("summary") or {}).get("risk_band")
    pct = (par.get("summary") or {}).get("max_epitope_identity_pct")
    who = (par.get("summary") or {}).get("max_epitope_identity_member")
    if band == "high":
        flags.append(_flag("epitope_conserved_across_paralogues", "major",
                           f"Highest epitope-region identity is {pct}% ({who}). A conserved epitope is a "
                           "specificity risk indicator; it is not measured cross-reactivity.",
                           "epitope_risk_band", pct))
    elif band == "moderate":
        flags.append(_flag("epitope_partially_conserved", "minor",
                           f"Highest epitope-region identity is {pct}% ({who}).", "epitope_risk_band", pct))
    elif band == "low":
        flags.append(_flag("epitope_divergent", "info",
                           f"Highest epitope-region identity is {pct}% ({who}); the epitope diverges "
                           "across the family members retrieved, which is a specificity opportunity.",
                           "epitope_risk_band", pct))
    else:
        flags.append(_flag("epitope_identity_uncomputable", "unknown",
                           "No epitope position could be aligned to any family member.",
                           "epitope_risk_band", None))
    return flags


# ---------------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------------


def decide_verdict(flags: list[dict[str, Any]], assessable: bool) -> dict[str, Any]:
    """Deterministic verdict from the flag set. Missing data can only hold a
    candidate back, never push it forward."""
    by_sev = {s: [f for f in flags if f["severity"] == s] for s in SEVERITY_ORDER}
    blockers, majors, unknowns = by_sev["blocker"], by_sev["major"], by_sev["unknown"]

    if not assessable:
        verdict, driver = "deprioritise", "insufficient_information"
        reason = ("Neither a sequence nor any interpretable interface metric was available, so there is "
                  "nothing to advance on. This reflects absence of information, not adverse information.")
    elif blockers:
        verdict, driver = "deprioritise", "blocker"
        reason = "Blocking flag(s): " + "; ".join(f["code"] for f in blockers)
    elif len(majors) >= TRIAGE_THRESHOLDS["max_majors_before_deprioritise"]:
        verdict, driver = "deprioritise", "accumulated_major_liabilities"
        reason = (f"{len(majors)} major liabilities: " + "; ".join(f["code"] for f in majors))
    elif majors or unknowns:
        verdict, driver = "conditional", "major_or_unknown"
        parts = []
        if majors:
            parts.append("major: " + ", ".join(f["code"] for f in majors))
        if unknowns:
            parts.append("unresolved: " + ", ".join(f["code"] for f in unknowns))
        reason = ("Advanceable only if these are resolved or engineered out -- " + "; ".join(parts))
    else:
        verdict, driver = "advance", "clean"
        reason = ("No blocking, major or unresolved flags; interface confidence is in an interpretable "
                  "band and every indicator required by this triage version was computed.")
    return {
        "verdict": verdict,
        "driver": driver,
        "reason": reason,
        "n_blocker": len(blockers),
        "n_major": len(majors),
        "n_unknown": len(unknowns),
        "n_minor": len(by_sev["minor"]),
        "driving_flags": [f["code"] for f in blockers + majors + unknowns],
    }


# ---------------------------------------------------------------------------------
# Per-candidate and batch entry points
# ---------------------------------------------------------------------------------


def triage_candidate(
    candidate: dict[str, Any],
    paralog_report: dict[str, Any] | None = None,
    known_disulfide_pairs: list[tuple[int, int]] | None = None,
) -> dict[str, Any]:
    """Triage one ``contracts.Candidate``-shaped dict.

    `candidate` needs: ``candidate_id``, optionally ``sequence``, ``metrics``
    (a list of ``contracts.Metric``-shaped dicts), ``run_id``, ``target_id``,
    ``origin``. `paralog_report` is the output of
    ``paralogs.paralog_cross_reactivity_risk`` for the *target* -- it is computed
    once per target, not once per candidate, and passed in here.

    Returns a JSON-able dict. Nothing in it is a measurement.
    """
    cid = candidate.get("candidate_id") or "<unnamed>"
    seq = candidate.get("sequence")
    origin = candidate.get("origin", "real")

    dev = compute_developability(seq, known_disulfide_pairs=known_disulfide_pairs)
    itf = interpret_interface_metrics(candidate.get("metrics"))
    flags = flag_developability(dev) + flag_interface(itf) + flag_cross_reactivity(paralog_report)
    assessable = bool(dev.get("computed")) or itf.get("confidence_band") is not None
    verdict = decide_verdict(flags, assessable)

    limitations = [
        "Every indicator here is computed from sequence or interpreted from a model score. None is a "
        "measurement.",
        "No binding, blockade, potency, cellular selectivity, immunogenicity or safety has been "
        "demonstrated for this candidate.",
        "Developability thresholds are conventions chosen for this triage version, not validated "
        "decision boundaries.",
    ]
    if origin == "fixture":
        limitations.insert(0, "origin=fixture: this record is a hand-written stand-in, not a design "
                              "result, and its triage output is not a biological result.")

    return {
        "triage_version": TRIAGE_VERSION,
        "candidate_id": cid,
        "run_id": candidate.get("run_id"),
        "design_id": candidate.get("design_id"),
        "target_id": candidate.get("target_id"),
        "origin": origin,
        "is_measurement": False,
        "verdict": verdict["verdict"],
        "verdict_detail": verdict,
        "claim": claim_language(itf.get("confidence_band")),
        "developability": dev,
        "interface": itf,
        "cross_reactivity": _slim_paralog(paralog_report),
        "flags": sorted(flags, key=lambda f: (SEVERITY_ORDER.index(f["severity"]), f["code"])),
        "contract_metric_dicts": as_contract_metric_dicts(dev),
        "limitations": limitations,
    }


def _slim_paralog(par: dict[str, Any] | None) -> dict[str, Any] | None:
    """Per-candidate copy of the target-level paralogue analysis, without repeating
    every per-position table for every candidate."""
    if par is None:
        return None
    return {
        "method_version": par.get("method_version"),
        "status": par.get("status"),
        "epitope_status": par.get("epitope_status"),
        "epitope_residues": par.get("epitope_residues"),
        "n_members_compared": (par.get("summary") or {}).get("n_members_compared"),
        "summary": par.get("summary"),
        "is_measurement": False,
        "limitations": par.get("limitations"),
    }


def _sort_key(row: dict[str, Any]) -> tuple:
    vd = row["verdict_detail"]
    band = row["interface"].get("confidence_band")
    risk = ((row.get("cross_reactivity") or {}).get("summary") or {}).get("risk_band")
    return (
        _VERDICT_RANK[row["verdict"]],
        _BAND_RANK.get(band, 4),
        vd["n_blocker"],
        vd["n_major"],
        vd["n_unknown"],
        _RISK_RANK.get(risk, 3),
        str(row["candidate_id"]),
    )


def triage_candidates(
    candidates: list[dict[str, Any]],
    target_accession: str | None = None,
    epitope_residues: list[int] | None = None,
    ensembl_gene_id: str | None = None,
    run_paralog_analysis: bool = True,
    max_paralog_members: int = 12,
    paralog_report: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Triage a batch of candidates against one target. Main entry point for the spine.

    The paralogue analysis is run ONCE for the target (it depends on the target and
    the engaged epitope, not on the individual binder) and attached to every
    candidate. Pass a precomputed `paralog_report` to skip the network entirely.

    Returns a JSON-able dict with a ranked list and a batch summary.
    """
    par = paralog_report
    par_error = None
    if par is None and run_paralog_analysis and target_accession:
        try:
            par = paralog_cross_reactivity_risk(
                target_accession,
                epitope_residues,
                ensembl_gene_id=ensembl_gene_id,
                max_members=max_paralog_members,
            )
        except Exception as exc:  # network / parsing; must not take the triage down
            par_error = f"{type(exc).__name__}: {exc}"
            par = None

    rows = [triage_candidate(c, paralog_report=par) for c in candidates]
    rows.sort(key=_sort_key)
    for i, r in enumerate(rows, start=1):
        r["rank"] = i

    counts = {v: sum(1 for r in rows if r["verdict"] == v) for v in ("advance", "conditional", "deprioritise")}
    origins = sorted({r["origin"] for r in rows})
    return {
        "triage_version": TRIAGE_VERSION,
        "thresholds": TRIAGE_THRESHOLDS,
        "ranking_key": RANKING_KEY_DESCRIPTION,
        "target_accession": target_accession,
        "epitope_residues": sorted(set(epitope_residues)) if epitope_residues else [],
        "n_candidates": len(rows),
        "verdict_counts": counts,
        "origins_present": origins,
        "contains_only_fixtures": origins == ["fixture"] if origins else False,
        "paralog_analysis": par,
        "paralog_analysis_error": par_error,
        "candidates": rows,
        "headline_caveat": (
            "Triage output ranks candidates by computed liabilities and by the design tool's own "
            "confidence scores. It does not demonstrate binding, blockade, potency, selectivity in "
            "cells, immunogenicity or safety, and it contains no measured affinity of any kind."
        ),
    }


def render_markdown(report: dict[str, Any], max_rows: int | None = None) -> str:
    """Compact human-readable triage table for the UI lane and the demo."""
    L: list[str] = []
    L.append(f"# Candidate triage ({report['triage_version']})")
    L.append("")
    if report.get("contains_only_fixtures"):
        L.append("> **Every record below is origin=`fixture`.** These are hand-written stand-ins used to "
                 "exercise the triage code. They are not design results and carry no biological meaning.")
        L.append("")
    counts = report["verdict_counts"]
    L.append(f"{report['n_candidates']} candidate(s): "
             f"{counts['advance']} advance, {counts['conditional']} conditional, "
             f"{counts['deprioritise']} deprioritise.")
    L.append("")
    par = report.get("paralog_analysis")
    if par:
        s = par.get("summary") or {}
        L.append(f"**Paralogue cross-reactivity (target {par.get('target_accession')}):** "
                 f"status `{par.get('status')}`, {s.get('n_members_compared')} family member(s) compared. "
                 f"{s.get('risk_basis')}")
        L.append("")
    L.append("| rank | candidate | origin | verdict | interface band | epitope risk | blockers | major | unknown |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    rows = report["candidates"][: max_rows or len(report["candidates"])]
    for r in rows:
        vd = r["verdict_detail"]
        band = r["interface"].get("confidence_band") or "unknown"
        risk = ((r.get("cross_reactivity") or {}).get("summary") or {}).get("risk_band") or "unknown"
        L.append(f"| {r['rank']} | {r['candidate_id']} | {r['origin']} | **{r['verdict']}** | {band} | "
                 f"{risk} | {vd['n_blocker']} | {vd['n_major']} | {vd['n_unknown']} |")
    L.append("")
    for r in rows:
        L.append(f"### {r['rank']}. {r['candidate_id']} — {r['verdict']}")
        L.append(f"{r['verdict_detail']['reason']}")
        L.append("")
        L.append(f"{r['claim']}")
        L.append("")
        for f in r["flags"]:
            L.append(f"- `{f['severity']}` **{f['code']}** — {f['message']}")
        L.append("")
    L.append("---")
    L.append(report["headline_caveat"])
    return "\n".join(L)
