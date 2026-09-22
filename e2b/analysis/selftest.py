"""Self-test for the Lane 6 triage.

Two parts, deliberately separated:

1. ``run_offline_selftest()`` -- pure-function checks on fixture sequences. No
   network. Asserts the properties the triage promises: determinism, that missing
   data never improves a verdict, that no banned affinity name is emitted, that an
   unrecognised metric is reported but not scored. Fixture output is never a
   biological result and is labelled ``origin="fixture"`` throughout.

2. ``run_online_selftest(...)`` -- the real-data path. Takes a *runtime-discovered*
   target: either a UniProt accession supplied by the caller, or an indication
   string that is resolved through the Open Targets adapter at runtime. No
   indication, disease id, gene or accession is hardcoded anywhere in this file or
   in the package: the CLI requires one to be supplied.

CLI:
    python -m e2b.analysis.selftest --offline
    python -m e2b.analysis.selftest --accession <UNIPROT_ACC> [--epitope 100-140]
    python -m e2b.analysis.selftest --indication "<free text>"
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any

from .developability import compute_developability
from .fixtures import fixture_candidates
from .interface import interpret_interface_metrics
from .paralogs import paralog_cross_reactivity_risk
from .triage import TRIAGE_VERSION, triage_candidate, triage_candidates

BANNED_METRIC_NAMES = ("kd", "ic50", "ec50", "affinity_nm", "binding_affinity")

# Assertive binding claims. A match inside an explicit negation ("no measured
# affinity", "has not demonstrated binding") is *the wording we want* and is not a
# violation, so the check strips negated occurrences before failing.
BANNED_CLAIM_PATTERNS = (
    r"binds\s+(?:to\s+)?(?:the\s+)?target",
    r"demonstrat(?:es|ed)\s+binding",
    r"measured\s+(?:kd|affinity|potency)",
    r"confirmed\s+binder",
)
_NEGATORS = ("no ", "not ", "never ", "cannot ", "without ", "nothing ")


def _assertive_claim_hits(text: str) -> list[str]:
    import re as _re

    hits: list[str] = []
    for pat in BANNED_CLAIM_PATTERNS:
        for m in _re.finditer(pat, text):
            prefix = text[max(0, m.start() - 24) : m.start()]
            if any(neg in prefix for neg in _NEGATORS):
                continue
            hits.append(m.group(0))
    return hits


# ---------------------------------------------------------------------------------
# Offline
# ---------------------------------------------------------------------------------


def run_offline_selftest() -> dict[str, Any]:
    checks: list[dict[str, Any]] = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        checks.append({"check": name, "passed": bool(ok), "detail": detail})

    cands = fixture_candidates()
    rep = triage_candidates(cands, run_paralog_analysis=False)

    # 1. determinism
    rep2 = triage_candidates(fixture_candidates(), run_paralog_analysis=False)
    same = [c["verdict"] for c in rep["candidates"]] == [c["verdict"] for c in rep2["candidates"]]
    same &= [c["candidate_id"] for c in rep["candidates"]] == [c["candidate_id"] for c in rep2["candidates"]]
    check("deterministic_verdicts_and_order", same,
          f"order={[c['candidate_id'] for c in rep['candidates']]}")

    # 2. fixture batch is labelled as such
    check("fixture_batch_flagged", rep["contains_only_fixtures"] is True,
          f"origins={rep['origins_present']}")

    # 3. no banned metric names anywhere in the output
    blob = json.dumps(rep).lower()
    offenders = [b for b in BANNED_METRIC_NAMES if f'"metric_name": "{b}"' in blob or f'"seq_{b}"' in blob]
    check("no_banned_affinity_metric_names", not offenders, f"offenders={offenders}")

    # 4. no assertive claim that a candidate binds
    offenders2 = _assertive_claim_hits(blob)
    check("no_binding_claims", not offenders2, f"offenders={offenders2}")

    # 4b. the report does carry the explicit "not a measurement" caveat
    check("caveat_present",
          "does not demonstrate binding" in rep["headline_caveat"]
          and all(c["is_measurement"] is False for c in rep["candidates"]),
          "headline_caveat and per-candidate is_measurement=False")

    # 5. missing interface metrics -> not advance
    byid = {c["candidate_id"]: c for c in rep["candidates"]}
    nm = byid["FIX-004-no-metrics"]
    check("missing_metrics_never_advance", nm["verdict"] != "advance",
          f"verdict={nm['verdict']} driver={nm['verdict_detail']['driver']}")

    # 6. removing information can only make the verdict equal or worse
    rank = {"advance": 0, "conditional": 1, "deprioritise": 2}
    clean = dict(cands[0])
    stripped = dict(clean)
    stripped["metrics"] = []
    a = triage_candidate(clean)
    b = triage_candidate(stripped)
    check("information_removal_is_monotone", rank[b["verdict"]] >= rank[a["verdict"]],
          f"with_metrics={a['verdict']} without_metrics={b['verdict']}")

    # 7. free cysteine detected on an odd count, and unresolved on an even count
    fc = compute_developability(cands[1]["sequence"])
    odd_ok = fc["indicators"]["unpaired_cysteine"]["value"] == 1
    ev = compute_developability("CAAAAC")
    even_unknown = ev["indicators"]["unpaired_cysteine"]["value"] is None
    check("free_cysteine_logic", odd_ok and even_unknown,
          f"odd_flags_free={odd_ok} even_is_unknown={even_unknown}")

    # 8. unrecognised metric reported but not scored
    vl = byid["FIX-005-verylow-nosequence"]
    unint = [u["metric_name"] for u in vl["interface"]["uninterpreted"]]
    check("unknown_metric_reported_not_scored", unint == ["some_new_tool_score"], f"uninterpreted={unint}")

    # 9. a None-valued metric is never treated as zero
    z = interpret_interface_metrics([{"metric_name": "iptm", "value": None, "model_name": "m"}])
    check("none_value_is_not_zero", z["confidence_band"] is None and z["n_confidence_missing"] == 1,
          f"band={z['confidence_band']}")

    # 10. no sequence -> every indicator None, and candidate is not advanced
    ns = compute_developability(None)
    check("no_sequence_yields_none_not_zero",
          ns["computed"] is False and ns["indicators"] == {} and vl["verdict"] != "advance",
          f"verdict={vl['verdict']}")

    # 11. very low confidence blocks
    check("very_low_confidence_blocks", vl["verdict"] == "deprioritise",
          f"verdict={vl['verdict']} driver={vl['verdict_detail']['driver']}")

    # 12. pI/GRAVY sanity on sequences with a known sign
    acidic = compute_developability("DDDDEEEEDDDDEEEE")["indicators"]
    basic = compute_developability("KKKKRRRRKKKKRRRR")["indicators"]
    check("pi_direction_sane",
          acidic["theoretical_pi"]["value"] < 5.0 < basic["theoretical_pi"]["value"],
          f"acidic_pI={acidic['theoretical_pi']['value']} basic_pI={basic['theoretical_pi']['value']}")
    check("net_charge_direction_sane",
          acidic["net_charge_ph7_4"]["value"] < 0 < basic["net_charge_ph7_4"]["value"],
          f"acidic_q={acidic['net_charge_ph7_4']['value']} basic_q={basic['net_charge_ph7_4']['value']}")

    passed = sum(1 for c in checks if c["passed"])
    return {
        "mode": "offline",
        "triage_version": TRIAGE_VERSION,
        "origin_of_inputs": "fixture",
        "n_checks": len(checks),
        "n_passed": passed,
        "all_passed": passed == len(checks),
        "checks": checks,
        "verdicts": {c["candidate_id"]: c["verdict"] for c in rep["candidates"]},
        "note": "Fixture sequences only. No biological conclusion follows from any of this.",
    }


# ---------------------------------------------------------------------------------
# Online
# ---------------------------------------------------------------------------------


def _resolve_accession_from_indication(indication: str) -> dict[str, Any]:
    """Discover a target at runtime through the Open Targets adapter. Nothing about
    the indication, the disease id or the gene is known to this module in advance."""
    from ..adapters.open_targets import OpenTargetsAdapter

    from ..adapters.structures import assess_accessibility, fetch_protein

    ot = OpenTargetsAdapter()
    matches, _ = ot.resolve_indication(indication)
    if not matches:
        return {"ok": False, "error": f"no disease match for {indication!r}"}
    disease_id = matches[0].disease_id
    targets, assoc, _ = ot.discover_targets(disease_id, limit=8)
    scan: list[dict[str, Any]] = []
    for t in targets:
        if not t.protein_accession:
            scan.append({"gene": t.gene_symbol, "skipped": "no UniProt accession"})
            continue
        try:
            prot, _p = fetch_protein(t.protein_accession)
            verdict = assess_accessibility(prot)
        except Exception as exc:
            scan.append({"gene": t.gene_symbol, "accession": t.protein_accession,
                         "skipped": f"topology lookup failed: {exc}"})
            continue
        scan.append({"gene": t.gene_symbol, "accession": t.protein_accession,
                     "accessibility": verdict["verdict"]})
        if verdict["accessible"]:
            return {
                "ok": True,
                "accession": t.protein_accession,
                "gene_symbol": t.gene_symbol,
                "ensembl_gene_id": t.stable_gene_id,
                "disease_id": disease_id,
                "disease_label": matches[0].label,
                "selection_rule": (
                    "highest-ranked Open Targets associated target that carries a UniProt accession "
                    "AND is topologically accessible to a soluble binder"
                ),
                "scan": scan,
                "pool_scope": assoc.get("scope"),
            }
    return {"ok": False, "error": "no binder-accessible target in the retrieved pool", "scan": scan}


def run_online_selftest(
    accession: str | None = None,
    indication: str | None = None,
    ensembl_gene_id: str | None = None,
    epitope: list[int] | None = None,
    max_members: int = 8,
) -> dict[str, Any]:
    """Exercise the real-data path: real UniProt sequences, real family retrieval,
    real epitope-restricted identity. Candidate *sequences* are still fixtures --
    the design lane has not produced real candidates yet, and that is stated in
    the output rather than papered over."""
    out: dict[str, Any] = {
        "mode": "online",
        "triage_version": TRIAGE_VERSION,
        "target_discovery": None,
        "candidate_sequences_origin": "fixture",
        "candidate_sequences_note": (
            "No real design candidates existed when this self-test was run. The target protein, its "
            "family and every epitope identity below are REAL UniProt/Ensembl data; the binder "
            "sequences are fixtures and their developability numbers are therefore not results."
        ),
        "paralog_analysis": None,
        "triage": None,
        "errors": [],
    }
    if accession is None:
        if not indication:
            out["errors"].append("supply either an accession or an indication; nothing is hardcoded")
            return out
        disc = _resolve_accession_from_indication(indication)
        out["target_discovery"] = disc
        if not disc.get("ok"):
            out["errors"].append(disc.get("error", "target discovery failed"))
            return out
        accession = disc["accession"]
        ensembl_gene_id = ensembl_gene_id or disc.get("ensembl_gene_id")
    else:
        out["target_discovery"] = {"ok": True, "accession": accession,
                                   "selection_rule": "accession supplied by the caller"}

    if epitope is None:
        from ..adapters.structures import assess_accessibility, fetch_protein
        from .epitope import epitope_proxy_from_accessibility

        try:
            prot, _ = fetch_protein(accession)
            acc_verdict = assess_accessibility(prot)
            prox = epitope_proxy_from_accessibility(prot, acc_verdict)
            out["epitope_resolution"] = prox
            out["accessibility"] = {k: acc_verdict[k] for k in ("verdict", "accessible", "basis")}
            epitope = prox["residues"] or None
        except Exception as exc:
            out["errors"].append(f"epitope proxy unavailable: {type(exc).__name__}: {exc}")
    else:
        out["epitope_resolution"] = {"residues": epitope, "is_proxy": False,
                                     "basis": "epitope positions supplied by the caller"}

    par = paralog_cross_reactivity_risk(
        accession, epitope, ensembl_gene_id=ensembl_gene_id, max_members=max_members
    )
    out["paralog_analysis"] = par
    rep = triage_candidates(
        fixture_candidates(),
        target_accession=accession,
        epitope_residues=epitope,
        paralog_report=par,
    )
    out["triage"] = rep
    return out


def _parse_epitope(spec: str | None) -> list[int] | None:
    if not spec:
        return None
    positions: list[int] = []
    for part in spec.split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-", 1)
            positions.extend(range(int(a), int(b) + 1))
        elif part:
            positions.append(int(part))
    return sorted(set(positions))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Lane 6 triage self-test")
    ap.add_argument("--offline", action="store_true", help="fixture-only checks, no network")
    ap.add_argument("--accession", default=None, help="UniProt accession of the target")
    ap.add_argument("--indication", default=None, help="free-text indication, resolved at runtime")
    ap.add_argument("--ensembl-gene", default=None)
    ap.add_argument("--epitope", default=None, help="1-based target positions, e.g. '110-125,140'")
    ap.add_argument("--max-members", type=int, default=8)
    ap.add_argument("--json-out", default=None)
    args = ap.parse_args(argv)

    if args.offline or (not args.accession and not args.indication):
        res = run_offline_selftest()
    else:
        res = run_online_selftest(
            accession=args.accession,
            indication=args.indication,
            ensembl_gene_id=args.ensembl_gene,
            epitope=_parse_epitope(args.epitope),
            max_members=args.max_members,
        )
    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump(res, fh, indent=2, default=str)
    print(json.dumps({k: v for k, v in res.items() if k not in ("triage", "paralog_analysis", "checks")},
                     indent=2, default=str))
    if res.get("mode") == "offline":
        for c in res["checks"]:
            print(("PASS " if c["passed"] else "FAIL ") + c["check"] + "  " + c["detail"])
        return 0 if res["all_passed"] else 1
    return 0 if not res.get("errors") else 1


if __name__ == "__main__":
    sys.exit(main())
