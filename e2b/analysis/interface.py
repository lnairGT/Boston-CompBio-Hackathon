"""Interpretation of the interface metrics the design lane supplies.

This module does **not** recompute anything structural. It takes the metrics that
came out of the co-folding / design tools, recognises the ones it has a stated
interpretation for, bands each one on its own scale, and refuses to touch the rest.

Two rules that shape the whole file:

1. **No composite score across incomparable metrics.** ipTM (0-1, unitless,
   AlphaFold-Multimer-family) and interface PAE (angstroms, same family) and a
   self-consistency RMSD (angstroms, from a different pipeline entirely) are not
   commensurable, and averaging them would invent a quantity with no meaning. The
   only aggregation performed here is a **conservative categorical floor**: the
   overall interface support is the *worst* band among the confidence metrics
   present. That is aggregation by "the weakest link governs", which is an
   explicit and defensible decision rule, not an arithmetic blend.

2. **An unrecognised metric is reported verbatim and never scored.** If the design
   lane emits a metric this registry does not know, it appears under
   ``uninterpreted`` with its value, and it pushes the confidence coverage down
   rather than being silently ignored.

Nothing here is a measurement of binding. The strongest statement any band
supports is "the model predicts this interface with <band> confidence".
"""

from __future__ import annotations

from typing import Any

INTERFACE_VERSION = "interface/1.1.0"

BANDS = ("high", "moderate", "low", "very_low")
_BAND_RANK = {"high": 3, "moderate": 2, "low": 1, "very_low": 0}


def _band_from_thresholds(value: float, thresholds: list[tuple[float, str]], higher_is_better: bool) -> str:
    """thresholds: descending cut list [(cut, band), ...] applied on the oriented value."""
    v = value if higher_is_better else -value
    for cut, band in thresholds:
        if v >= cut:
            return band
    return thresholds[-1][1]


# name (normalised) -> spec. `role`:
#   "confidence" -> participates in the conservative floor
#   "geometry"   -> produces a hard geometry flag, never a confidence band
#   "context"    -> reported, never scored
METRIC_REGISTRY: dict[str, dict[str, Any]] = {
    "iptm": {
        "role": "confidence",
        "units": "unitless 0-1",
        "higher_is_better": True,
        "thresholds": [(0.80, "high"), (0.60, "moderate"), (0.50, "low"), (float("-inf"), "very_low")],
        "source_family": "AlphaFold-Multimer / AF3 / Boltz / Chai co-folding confidence",
        "caveat": "ipTM is a predicted-accuracy score for the predicted complex. It correlates with, "
                  "but does not measure, binding; high ipTM on a designed binder is routinely "
                  "obtained for sequences that do not bind experimentally.",
    },
    "ptm": {
        "role": "context",
        "units": "unitless 0-1",
        "higher_is_better": True,
        "source_family": "AlphaFold-family global predicted TM-score",
        "caveat": "Whole-complex score dominated by the larger chain; not interface-specific.",
    },
    "pae_interface": {
        "role": "confidence",
        "units": "angstrom",
        "higher_is_better": False,
        "thresholds": [(-5.0, "high"), (-10.0, "moderate"), (-15.0, "low"), (float("-inf"), "very_low")],
        "source_family": "AlphaFold-family predicted aligned error, restricted to inter-chain pairs",
        "caveat": "Expected positional error between chains, in angstroms. Low PAE means the model is "
                  "confident about the relative placement, not that the placement is correct.",
    },
    "plddt_binder": {
        "role": "confidence",
        "units": "unitless 0-100",
        "higher_is_better": True,
        "thresholds": [(90.0, "high"), (70.0, "moderate"), (50.0, "low"), (float("-inf"), "very_low")],
        "source_family": "AlphaFold-family per-residue confidence, averaged over the binder chain",
        "caveat": "Per-residue local confidence. A confidently-folded binder can still be docked wrongly.",
    },
    "plddt_complex": {
        "role": "context",
        "units": "unitless 0-100",
        "higher_is_better": True,
        "source_family": "AlphaFold-family mean pLDDT over the complex",
        "caveat": "Dominated by the target chain when the target is much larger than the binder.",
    },
    "self_consistency_rmsd": {
        "role": "confidence",
        "units": "angstrom",
        "higher_is_better": False,
        "thresholds": [(-1.0, "high"), (-2.0, "moderate"), (-4.0, "low"), (float("-inf"), "very_low")],
        "source_family": "design self-consistency: refolded designed sequence vs designed backbone "
                         "(RFdiffusion/ProteinMPNN convention, in-silico success usually scRMSD < 2 A)",
        "caveat": "Measures agreement between two models, not agreement with reality.",
    },
    "n_clashes": {
        "role": "geometry",
        "units": "count",
        "higher_is_better": False,
        "source_family": "steric clash count at the modelled interface",
        "caveat": "Clash definition is tool-specific; a nonzero count from one tool is not comparable "
                  "to a nonzero count from another.",
    },
    "n_interface_contacts": {
        "role": "context",
        "units": "count",
        "higher_is_better": True,
        "source_family": "inter-chain residue contacts under the tool's distance cutoff",
        "caveat": "Contact counts scale with the cutoff used; a small count may mean a small interface "
                  "or a strict cutoff.",
    },
    "interface_buried_sasa": {
        "role": "context",
        "units": "angstrom^2",
        "higher_is_better": True,
        "source_family": "buried solvent-accessible surface area at the interface",
        "caveat": "Computed on a model, so it inherits the model's error.",
    },
}

# Accepted aliases -> canonical registry key. Design tools disagree on names.
ALIASES = {
    "iptm": "iptm", "ipTM": "iptm", "iptm_score": "iptm", "interface_ptm": "iptm",
    "ptm": "ptm", "pTM": "ptm",
    "pae_interface": "pae_interface", "interface_pae": "pae_interface", "ipae": "pae_interface",
    "pae_int": "pae_interface", "mean_interface_pae": "pae_interface",
    "plddt_binder": "plddt_binder", "binder_plddt": "plddt_binder", "plddt_design": "plddt_binder",
    "plddt": "plddt_complex", "plddt_complex": "plddt_complex", "mean_plddt": "plddt_complex",
    "self_consistency_rmsd": "self_consistency_rmsd", "scrmsd": "self_consistency_rmsd",
    "rmsd_self_consistency": "self_consistency_rmsd", "design_rmsd": "self_consistency_rmsd",
    "n_clashes": "n_clashes", "clashes": "n_clashes", "clash_count": "n_clashes",
    "n_interface_contacts": "n_interface_contacts", "interface_contacts": "n_interface_contacts",
    "contact_count": "n_interface_contacts",
    "interface_buried_sasa": "interface_buried_sasa", "buried_sasa": "interface_buried_sasa",
    "dsasa": "interface_buried_sasa",
}

CLASH_MAJOR_THRESHOLD = 5


def normalise_name(name: str) -> str | None:
    key = (name or "").strip()
    if key in ALIASES:
        return ALIASES[key]
    key2 = key.lower().replace(" ", "_").replace("-", "_")
    return ALIASES.get(key2)


def interpret_interface_metrics(metrics: list[dict[str, Any]] | None) -> dict[str, Any]:
    """Band each supplied interface metric on its own scale.

    `metrics` is a list of ``contracts.Metric``-shaped dicts (``metric_name``,
    ``value``, ``method``, ``model_name``, ``model_version`` ...), i.e. exactly what
    ``Candidate.metrics`` carries. Values of None are treated as *not computed*
    and lower the coverage; they are never treated as zero.

    Returns a JSON-able dict with per-metric bands, the conservative overall band,
    geometry flags, and the metrics it could not interpret.
    """
    out: dict[str, Any] = {
        "method_version": INTERFACE_VERSION,
        "is_measurement": False,
        "per_metric": [],
        "uninterpreted": [],
        "geometry_flags": [],
        "confidence_band": None,
        "confidence_basis": [],
        "n_confidence_metrics": 0,
        "n_confidence_missing": 0,
        "aggregation_rule": (
            "conservative categorical floor: overall band = worst band among the recognised "
            "confidence metrics. No arithmetic combination of different scales is performed."
        ),
        "models_reported": [],
        "notes": [],
    }
    if not metrics:
        out["notes"].append(
            "no interface metrics supplied; interface support is UNKNOWN (not weak, not strong)"
        )
        return out

    models: set[str] = set()
    bands: list[tuple[str, str]] = []
    for m in metrics:
        raw_name = m.get("metric_name") or m.get("name") or ""
        value = m.get("value")
        model = m.get("model_name")
        if model:
            models.add(str(model) + (f":{m['model_version']}" if m.get("model_version") else ""))
        canon = normalise_name(raw_name)
        if canon is None:
            out["uninterpreted"].append(
                {
                    "metric_name": raw_name,
                    "value": value,
                    "model_name": model,
                    "reason": "no stated interpretation for this metric name in this triage version; "
                              "reported verbatim, not scored",
                }
            )
            continue
        spec = METRIC_REGISTRY[canon]
        row: dict[str, Any] = {
            "metric_name": raw_name,
            "canonical_name": canon,
            "value": value,
            "units": spec["units"],
            "role": spec["role"],
            "higher_is_better": spec["higher_is_better"],
            "model_name": model,
            "model_version": m.get("model_version"),
            "source_family": spec["source_family"],
            "caveat": spec["caveat"],
            "band": None,
        }
        if value is None:
            row["band"] = None
            row["note"] = "value is None (not computed); treated as missing, never as zero"
            if spec["role"] == "confidence":
                out["n_confidence_missing"] += 1
        elif spec["role"] == "confidence":
            row["band"] = _band_from_thresholds(float(value), spec["thresholds"], spec["higher_is_better"])
            row["thresholds"] = [
                {"at_or_better_than": (c if spec["higher_is_better"] else -c), "band": b}
                for c, b in spec["thresholds"]
                if c not in (float("-inf"),)
            ]
            bands.append((canon, row["band"]))
            out["n_confidence_metrics"] += 1
        elif spec["role"] == "geometry":
            if canon == "n_clashes" and float(value) > 0:
                sev = "major" if float(value) >= CLASH_MAJOR_THRESHOLD else "minor"
                out["geometry_flags"].append(
                    {
                        "flag": "steric_clashes_at_interface",
                        "severity": sev,
                        "value": value,
                        "detail": f"{value} clashes reported by {model or 'the design tool'}; "
                                  f">= {CLASH_MAJOR_THRESHOLD} treated as major",
                    }
                )
        out["per_metric"].append(row)

    if bands:
        worst = min(bands, key=lambda kb: _BAND_RANK[kb[1]])
        out["confidence_band"] = worst[1]
        out["confidence_basis"] = [{"metric": k, "band": b} for k, b in bands]
        out["notes"].append(
            f"overall band '{worst[1]}' is set by '{worst[0]}', the weakest recognised confidence metric"
        )
    else:
        out["notes"].append(
            "no recognised confidence metric had a value; interface support is UNKNOWN"
        )
    if out["n_confidence_missing"]:
        out["notes"].append(
            f"{out['n_confidence_missing']} recognised confidence metric(s) had value None"
        )
    if out["uninterpreted"]:
        out["notes"].append(
            f"{len(out['uninterpreted'])} metric(s) had no stated interpretation and were not scored"
        )
    out["models_reported"] = sorted(models)
    return out


def claim_language(band: str | None) -> str:
    """The strongest honest sentence for a given band. Used by the report layer so
    that the wording cannot drift into 'binds'."""
    if band is None:
        return (
            "No interface confidence was supplied, so nothing can be said about a predicted "
            "interface for this candidate."
        )
    phrase = {
        "high": "with high model confidence",
        "moderate": "with moderate model confidence",
        "low": "with low model confidence",
        "very_low": "with very low model confidence",
    }[band]
    return (
        f"A structure-prediction model predicts an interface between this designed sequence and the "
        f"target {phrase}. This is a model prediction, not a measured interaction: no binding, "
        f"blockade, potency or selectivity has been demonstrated."
    )
