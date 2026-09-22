"""Stage 2 - score and rank the candidate targets.

Design rules
------------
*Nothing is hidden.* Every target from the stage 1 pool appears in the output
table, including the ones that fail the accessibility filter; they carry
``passes_accessibility`` false and the reason. A filtered-out target is a
recorded decision, not a missing row.

*Every rank is traceable.* The table carries the raw value and the weighted
contribution of each component, so any score can be reconstructed by hand
from its columns. Re-weighting needs no re-fetching: pass a different
weights file and re-run.

*Accessibility is a gate, not a term.* A binder that cannot reach its target
is not a worse binder, it is not a binder. So accessibility multiplies the
composite by zero rather than subtracting from it - which is why the column
is boolean and the reason is preserved.

A caution on the known-drug and literature components: Open Targets
association scores already incorporate known-drug and literature evidence, so
those components partly re-express information also present in
``ot_overall``. This is deliberate (they are the modality- and
maturity-relevant parts of that evidence) but it means the composite is not a
sum of independent signals, and a heavily drugged target ranks partly on its
own pharmacology. The per-component columns are there so you can see it.
"""

from __future__ import annotations

import json
import logging
import math
from pathlib import Path
from typing import Any

import pandas as pd

from . import __version__
from .config import (
    AB_BUCKET_SCORES,
    DEFAULT_PENALTIES,
    SATURATION,
    STAGE_FILES,
    load_weights,
)

log = logging.getLogger("ind2b.stage2")

SCHEMA_VERSION = 1


def _log_norm(count: float | None, saturation: float) -> float:
    """Compress a count to 0-1, reaching ~1.0 at ``saturation``."""
    if not count or count <= 0:
        return 0.0
    return min(1.0, math.log1p(count) / math.log1p(saturation))


def _mean(values: list[float]) -> float:
    vals = [v for v in values if v is not None]
    return sum(vals) / len(vals) if vals else 0.0


def components(target: dict) -> dict[str, float]:
    """Raw 0-1 component scores for one stage 1 target record."""
    ds = target.get("datatype_scores") or {}
    kd = target.get("known_drugs") or {}
    ch = target.get("chembl") or {}
    tr = target.get("trials") or {}
    expr = (target.get("expression") or {}).get("max_specificity") or {}

    # Genetic support: for an oncology indication somatic evidence carries
    # most of the weight, so germline and somatic are averaged over whatever
    # is present rather than requiring both.
    genetic = _mean(
        [ds.get("genetic_association"), ds.get("somatic_mutation")]
    )

    ab_buckets = (target.get("tractability") or {}).get("AB") or []
    ab_score = max((AB_BUCKET_SCORES.get(b, 0.0) for b in ab_buckets), default=0.0)
    if kd.get("has_biologic_precedent"):
        ab_score = max(ab_score, 0.9)

    known_drug = _log_norm(kd.get("count"), SATURATION["known_drugs"])
    if kd.get("has_approved"):
        known_drug = max(known_drug, 0.8)

    return {
        "ot_overall": float(target.get("ot_overall_score") or 0.0),
        "genetic_evidence": genetic,
        "known_drug": known_drug,
        "literature": float(ds.get("literature") or 0.0),
        "expression_specificity": float(expr.get("specificity_score") or 0.0),
        "pathway_and_model": _mean(
            [ds.get("affected_pathway"), ds.get("animal_model")]
        ),
        "clinical_precedent": max(
            _log_norm(tr.get("total_trials"), SATURATION["trials"]),
            _log_norm(ch.get("n_distinct_mechanisms"), SATURATION["chembl_mechanisms"]),
        ),
        "antibody_tractability": ab_score,
    }


def penalty(target: dict) -> tuple[float, list[str]]:
    """Multiplicative penalty for safety liabilities and broad expression.

    The safety term is flat, not per-liability: see the note on
    ``DEFAULT_PENALTIES`` in ``config`` for why compounding it measures study
    intensity rather than risk.
    """
    reasons: list[str] = []
    mult = 1.0

    n_safety = (target.get("safety_liabilities") or {}).get("count") or 0
    if n_safety:
        safety_mult = DEFAULT_PENALTIES["safety_liability_any"]
        mult *= safety_mult
        reasons.append(
            f"{n_safety} safety liabilit{'y' if n_safety == 1 else 'ies'} "
            f"recorded (flat x{safety_mult:.2f})"
        )

    expr = (target.get("expression") or {}).get("max_specificity") or {}
    spec = expr.get("specificity_score")
    if spec is not None and spec < 0.10:
        mult *= DEFAULT_PENALTIES["broad_expression"]
        reasons.append(
            f"low tissue specificity {spec:.2f} (x{DEFAULT_PENALTIES['broad_expression']:.2f})"
        )

    return mult, reasons


def score_targets(
    evidence: dict, *, weights_file: str | Path | None = None
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Score every stage 1 target; returns the table and the weights used."""
    weights = load_weights(weights_file)
    total_w = sum(weights.values())
    if abs(total_w - 1.0) > 1e-6:
        log.warning("weights sum to %.4f, normalising", total_w)

    rows: list[dict[str, Any]] = []
    for t in evidence["targets"]:
        access = t.get("accessibility") or {}
        passes = bool(t.get("surface_accessible"))
        comp = components(t)
        pen_mult, pen_reasons = penalty(t)

        weighted = {k: comp[k] * weights[k] / total_w for k in weights}
        composite = sum(weighted.values()) * pen_mult
        if not passes:
            composite = 0.0

        spans = access.get("extracellular_spans") or []
        row: dict[str, Any] = {
            "symbol": t.get("symbol"),
            "ensembl_id": t.get("ensembl_id"),
            "uniprot": (t.get("uniprot") or [None])[0],
            "name": t.get("name"),
            "composite_score": round(composite, 5),
            "passes_accessibility": passes,
            "accessibility_mode": access.get("mode"),
            "accessibility_reason": access.get("reason"),
            "extracellular_spans": ";".join(f"{s}-{e}" for s, e in spans),
            "ectodomain_residues": sum(e - s + 1 for s, e in spans),
            "penalty_multiplier": round(pen_mult, 4),
            "penalty_reasons": "; ".join(pen_reasons),
            "target_class": "; ".join(t.get("target_class") or []),
            "n_safety_liabilities": (t.get("safety_liabilities") or {}).get("count"),
            "n_known_drugs": (t.get("known_drugs") or {}).get("count"),
            "has_approved_drug": (t.get("known_drugs") or {}).get("has_approved"),
            "has_biologic_precedent": (t.get("known_drugs") or {}).get(
                "has_biologic_precedent"
            ),
            "drug_types": json.dumps((t.get("known_drugs") or {}).get("drug_types") or {}),
            "n_chembl_mechanisms": (t.get("chembl") or {}).get("n_distinct_mechanisms"),
            "chembl_mechanisms": "; ".join((t.get("chembl") or {}).get("mechanisms") or []),
            "n_trials": (t.get("trials") or {}).get("total_trials"),
            "n_stopped_trials": len((t.get("trials") or {}).get("stopped_trials") or []),
            "example_nct_ids": ";".join((t.get("trials") or {}).get("example_nct_ids") or []),
            "n_interaction_partners": (t.get("interactions") or {}).get("n_partners"),
            "top_partners": ";".join(
                p["symbol"] for p in ((t.get("interactions") or {}).get("partners") or [])[:8]
            ),
            "ot_rank": t.get("rank_by_ot_score"),
            "layers_fetched": t.get("layers_fetched"),
        }
        row.update({f"raw_{k}": round(v, 5) for k, v in comp.items()})
        row.update({f"wtd_{k}": round(v, 5) for k, v in weighted.items()})
        rows.append(row)

    df = pd.DataFrame(rows)
    df = df.sort_values(
        ["passes_accessibility", "composite_score"], ascending=[False, False]
    ).reset_index(drop=True)
    df.insert(0, "rank", [i + 1 if p else None for i, p in enumerate(df["passes_accessibility"])])
    return df, weights


def write(
    df: pd.DataFrame, weights: dict[str, float], evidence: dict, run_dir: Path
) -> tuple[Path, Path]:
    run_dir.mkdir(parents=True, exist_ok=True)
    csv_path = run_dir / STAGE_FILES[2]
    df.to_csv(csv_path, index=False)

    meta_path = run_dir / "stage2_scoring_meta.json"
    n_pass = int(df["passes_accessibility"].sum())
    meta_path.write_text(
        json.dumps(
            {
                "schema_version": SCHEMA_VERSION,
                "package_version": __version__,
                "disease_id": evidence.get("disease_id"),
                "disease_name": evidence.get("disease_name"),
                "pool_size": len(df),
                "n_passing_accessibility": n_pass,
                "weights": weights,
                "penalties": DEFAULT_PENALTIES,
                "ab_bucket_scores": AB_BUCKET_SCORES,
                "saturation": SATURATION,
                "caveats": [
                    "Open Targets association scores already include known-drug "
                    "and literature evidence, so those components are not "
                    "independent of ot_overall.",
                    "Ranking reflects weight of existing evidence, not causality.",
                    "Accessibility is a hard gate: failing targets score 0.",
                ],
            },
            indent=2,
        )
    )
    return csv_path, meta_path


def read(run_dir: Path) -> pd.DataFrame:
    path = run_dir / STAGE_FILES[2]
    if not path.exists():
        raise FileNotFoundError(f"{path} not found - run stage 2 first")
    return pd.read_csv(path)
