"""Sequence-derived developability indicators for a designed minibinder.

Every value here is **computed from the amino-acid sequence alone**. None of it is
a measurement. A computed pI is not a measured pI; a hydrophobic-patch fraction is
not an observed aggregation rate; a sequon count is not observed glycosylation.
They are liability *flags* that tell a team what to look at first, and they are
reported with the method that produced them so a reviewer can disagree with the
method rather than with a bare number.

Implemented here rather than imported so that each formula is visible and
versioned with the triage. No third-party dependency is used.
"""

from __future__ import annotations

import math
import re
from dataclasses import asdict, dataclass, field
from typing import Any

AA = "ACDEFGHIKLMNPQRSTVWY"

# Kyte & Doolittle (1982) J Mol Biol 157:105-132 hydropathy index.
KD_HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# EMBOSS (iep) pKa set. Stated explicitly because pI values differ by ~0.2-0.5 pH
# units between the EMBOSS, Bjellqvist/ExPASy and Sillero sets for the same sequence.
PKA_EMBOSS = {
    "n_term": 8.6, "c_term": 3.6,
    "D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1, "H": 6.5, "K": 10.8, "R": 12.5,
}
_POSITIVE = ("K", "R", "H")
_NEGATIVE = ("D", "E", "C", "Y")

METHOD_VERSION = "developability/1.1.0"

# Thresholds used by the flagging layer. Versioned with the module; every one of
# them is a convention, not a law of nature, and docs/analysis_methods.md says so.
APR_WINDOW = 5
APR_MEAN_KD = 1.5
LOWCOMPLEXITY_WINDOW = 12   # SEG default window
LOWCOMPLEXITY_ENTROPY_BITS = 2.2  # SEG default trigger complexity


@dataclass
class Indicator:
    """One computed indicator. `value` is None when it could not be computed."""

    name: str
    value: float | int | None
    units: str | None
    method: str
    higher_is: str  # "better" | "worse" | "context"
    detail: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _charge_at_ph(seq: str, ph: float) -> float:
    counts = {a: seq.count(a) for a in set(seq)}
    pos = 1.0 / (1.0 + 10 ** (ph - PKA_EMBOSS["n_term"]))
    for a in _POSITIVE:
        if counts.get(a):
            pos += counts[a] / (1.0 + 10 ** (ph - PKA_EMBOSS[a]))
    neg = 1.0 / (1.0 + 10 ** (PKA_EMBOSS["c_term"] - ph))
    for a in _NEGATIVE:
        if counts.get(a):
            neg += counts[a] / (1.0 + 10 ** (PKA_EMBOSS[a] - ph))
    return pos - neg


def theoretical_pi(seq: str, tol: float = 1e-4) -> float | None:
    """pH at which the Henderson-Hasselbalch net charge is zero, by bisection."""
    if not seq:
        return None
    lo, hi = 0.0, 14.0
    for _ in range(200):
        mid = (lo + hi) / 2.0
        c = _charge_at_ph(seq, mid)
        if abs(c) < tol or (hi - lo) < tol:
            return round(mid, 3)
        if c > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2.0, 3)


def gravy(seq: str) -> float | None:
    """Grand average of hydropathy: mean Kyte-Doolittle index over all residues."""
    vals = [KD_HYDROPATHY[a] for a in seq if a in KD_HYDROPATHY]
    if not vals:
        return None
    return round(sum(vals) / len(vals), 4)


def hydrophobic_patches(seq: str, window: int = APR_WINDOW, mean_kd: float = APR_MEAN_KD) -> dict[str, Any]:
    """Contiguous windows with high mean hydropathy and no formal charge.

    Heuristic proxy for aggregation-prone regions: a window of `window` residues
    whose mean Kyte-Doolittle hydropathy is >= `mean_kd` and which contains no
    D/E/K/R. This is a coarse stand-in for a trained APR predictor (TANGO,
    AGGRESCAN, CamSol) and will both over- and under-call relative to them.
    """
    n = len(seq)
    if n < window:
        return {"n_patches": 0, "patches": [], "fraction_residues_in_patch": 0.0 if n else None}
    covered = set()
    raw: list[tuple[int, int]] = []
    for i in range(0, n - window + 1):
        w = seq[i : i + window]
        if any(c in w for c in "DEKR"):
            continue
        vals = [KD_HYDROPATHY.get(a, 0.0) for a in w]
        if sum(vals) / window >= mean_kd:
            raw.append((i + 1, i + window))
            covered.update(range(i + 1, i + window + 1))
    merged: list[list[int]] = []
    for s, e in raw:
        if merged and s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return {
        "n_patches": len(merged),
        "patches": [{"start": s, "end": e, "sequence": seq[s - 1 : e]} for s, e in merged],
        "fraction_residues_in_patch": round(len(covered) / n, 4),
    }


def cysteines(seq: str, known_disulfide_pairs: list[tuple[int, int]] | None = None) -> dict[str, Any]:
    """Cysteine census and an explicitly-qualified unpaired-cysteine call.

    A free thiol in a secreted protein is a real manufacturability liability
    (disulfide scrambling, covalent dimer, heterogeneity). Sequence alone cannot
    prove pairing: an ODD cysteine count guarantees at least one unpaired
    cysteine, while an EVEN count proves nothing. So `unpaired_inferred` is
    True on odd counts and None on even counts unless pairing is supplied.
    """
    pos = [i + 1 for i, a in enumerate(seq) if a == "C"]
    n = len(pos)
    result: dict[str, Any] = {
        "count": n,
        "positions": pos,
        "parity": "even" if n % 2 == 0 else "odd",
        "unpaired_inferred": None,
        "n_unpaired": None,
        "basis": None,
    }
    if known_disulfide_pairs:
        paired = {p for pair in known_disulfide_pairs for p in pair}
        free = [p for p in pos if p not in paired]
        result["unpaired_inferred"] = len(free) > 0
        result["n_unpaired"] = len(free)
        result["free_positions"] = free
        result["basis"] = "explicit disulfide pairing supplied by caller"
    elif n == 0:
        result["unpaired_inferred"] = False
        result["n_unpaired"] = 0
        result["basis"] = "no cysteines present"
    elif n % 2 == 1:
        result["unpaired_inferred"] = True
        result["n_unpaired"] = 1  # at least
        result["basis"] = "odd cysteine count guarantees >=1 free thiol"
    else:
        result["basis"] = (
            "even cysteine count; pairing cannot be determined from sequence alone "
            "(treated as unknown, not as paired)"
        )
    return result


_SEQUON = re.compile(r"(?=(N[^P][ST]))")


def n_glycosylation_sequons(seq: str) -> dict[str, Any]:
    """N-X-S/T with X != P (Asn-linked sequon). Presence of a sequon is not
    evidence of occupancy; occupancy depends on structure and expression host."""
    sites = [{"position": m.start() + 1, "motif": m.group(1)} for m in _SEQUON.finditer(seq)]
    return {"count": len(sites), "sites": sites}


def oxidation_sites(seq: str) -> dict[str, Any]:
    """Methionine and tryptophan positions (oxidation-prone side chains).

    Solvent exposure dominates real oxidation risk and is NOT assessed here;
    these are candidate positions, not predicted oxidation."""
    met = [i + 1 for i, a in enumerate(seq) if a == "M"]
    trp = [i + 1 for i, a in enumerate(seq) if a == "W"]
    return {
        "met_positions": met,
        "trp_positions": trp,
        "count": len(met) + len(trp),
        "exposure_assessed": False,
    }


def deamidation_motifs(seq: str) -> dict[str, Any]:
    """NG (high-rate) and NS/NT (lower-rate) asparagine deamidation motifs, plus
    DG/DP aspartate isomerisation motifs reported separately as context."""
    def _find(pat: str) -> list[int]:
        return [m.start() + 1 for m in re.finditer(f"(?={pat})", seq)]

    ng = _find("NG")
    ns = _find("N[ST]")
    dg = _find("D[GP]")
    return {
        "ng_positions": ng,
        "ns_nt_positions": ns,
        "high_risk_count": len(ng),
        "moderate_risk_count": len(ns),
        "isomerisation_dg_dp_positions": dg,
        "count": len(ng) + len(ns),
    }


def _shannon_bits(window: str) -> float:
    n = len(window)
    if n == 0:
        return 0.0
    h = 0.0
    for a in set(window):
        p = window.count(a) / n
        h -= p * math.log2(p)
    return h


def low_complexity(seq: str, window: int = LOWCOMPLEXITY_WINDOW, bits: float = LOWCOMPLEXITY_ENTROPY_BITS) -> dict[str, Any]:
    """Low-complexity regions by Shannon entropy over a sliding window (SEG-like),
    plus exact tandem repeats of period 1-5 with >= 3 copies."""
    n = len(seq)
    covered = set()
    regions: list[list[int]] = []
    if n >= window:
        for i in range(0, n - window + 1):
            if _shannon_bits(seq[i : i + window]) < bits:
                s, e = i + 1, i + window
                covered.update(range(s, e + 1))
                if regions and s <= regions[-1][1] + 1:
                    regions[-1][1] = max(regions[-1][1], e)
                else:
                    regions.append([s, e])
    repeats = []
    for period in range(1, 6):
        i = 0
        while i + period * 3 <= n:
            unit = seq[i : i + period]
            copies = 1
            j = i + period
            while seq[j : j + period] == unit:
                copies += 1
                j += period
            if copies >= 3:
                repeats.append({"start": i + 1, "end": j, "unit": unit, "copies": copies})
                i = j
            else:
                i += 1
    return {
        "n_low_complexity_regions": len(regions),
        "low_complexity_regions": [{"start": s, "end": e, "sequence": seq[s - 1 : e]} for s, e in regions],
        "fraction_low_complexity": round(len(covered) / n, 4) if n else None,
        "tandem_repeats": repeats,
    }


def compute_developability(sequence: str | None, known_disulfide_pairs: list[tuple[int, int]] | None = None) -> dict[str, Any]:
    """All sequence-derived indicators for one candidate sequence.

    Returns a JSON-able dict. If `sequence` is None or empty every indicator is
    None and `computed` is False -- the triage layer then treats every one of
    them as missing, which makes the verdict more conservative, never less.
    """
    out: dict[str, Any] = {
        "method_version": METHOD_VERSION,
        "computed": False,
        "is_measurement": False,
        "sequence_length": None,
        "indicators": {},
        "raw": {},
        "notes": [],
    }
    if not sequence:
        out["notes"].append("no sequence supplied; all sequence indicators are None")
        return out
    seq = sequence.strip().upper()
    bad = sorted(set(seq) - set(AA))
    if bad:
        out["notes"].append(f"non-standard residues present: {bad}; indicators not computed")
        return out

    n = len(seq)
    pi = theoretical_pi(seq)
    gv = gravy(seq)
    q74 = round(_charge_at_ph(seq, 7.4), 3)
    apr = hydrophobic_patches(seq)
    cys = cysteines(seq, known_disulfide_pairs)
    gly = n_glycosylation_sequons(seq)
    ox = oxidation_sites(seq)
    deam = deamidation_motifs(seq)
    lc = low_complexity(seq)

    ind = {
        "length": Indicator("length", n, "residues", "len(sequence)", "context"),
        "theoretical_pi": Indicator(
            "theoretical_pi", pi, "pH units",
            "Henderson-Hasselbalch net charge zero-crossing, bisection, EMBOSS pKa set", "context",
            {"pka_set": "EMBOSS"},
        ),
        "gravy": Indicator(
            "gravy", gv, "Kyte-Doolittle units",
            "mean Kyte-Doolittle (1982) hydropathy over all residues", "worse",
        ),
        "net_charge_ph7_4": Indicator(
            "net_charge_ph7_4", q74, "elementary charges",
            "Henderson-Hasselbalch sum over ionisable groups at pH 7.4, EMBOSS pKa set", "context",
            {"pka_set": "EMBOSS"},
        ),
        "hydrophobic_patch_fraction": Indicator(
            "hydrophobic_patch_fraction", apr["fraction_residues_in_patch"], "fraction of residues",
            f"fraction of residues inside any {APR_WINDOW}-residue window with mean KD hydropathy "
            f">= {APR_MEAN_KD} and no D/E/K/R", "worse",
            {"n_patches": apr["n_patches"]},
        ),
        "cysteine_count": Indicator("cysteine_count", cys["count"], "residues", "count of C", "context"),
        "unpaired_cysteine": Indicator(
            "unpaired_cysteine",
            None if cys["unpaired_inferred"] is None else int(cys["unpaired_inferred"]),
            "boolean (1 = at least one free thiol inferred)",
            cys["basis"] or "", "worse",
            {"parity": cys["parity"], "positions": cys["positions"]},
        ),
        "n_glycosylation_sequons": Indicator(
            "n_glycosylation_sequons", gly["count"], "count",
            "regex N[^P][ST]; sequon presence, not occupancy", "worse",
            {"sites": gly["sites"]},
        ),
        "oxidation_prone_sites": Indicator(
            "oxidation_prone_sites", ox["count"], "count",
            "count of Met and Trp residues; solvent exposure not assessed", "worse",
            {"met": ox["met_positions"], "trp": ox["trp_positions"]},
        ),
        "deamidation_motifs": Indicator(
            "deamidation_motifs", deam["count"], "count",
            "NG (high rate) plus NS/NT (lower rate) motif occurrences", "worse",
            {"ng": deam["ng_positions"], "ns_nt": deam["ns_nt_positions"],
             "isomerisation_dg_dp": deam["isomerisation_dg_dp_positions"]},
        ),
        "fraction_low_complexity": Indicator(
            "fraction_low_complexity", lc["fraction_low_complexity"], "fraction of residues",
            f"fraction inside any {LOWCOMPLEXITY_WINDOW}-residue window with Shannon entropy "
            f"< {LOWCOMPLEXITY_ENTROPY_BITS} bits (SEG-like)", "worse",
            {"n_regions": lc["n_low_complexity_regions"], "tandem_repeats": lc["tandem_repeats"]},
        ),
    }

    out["computed"] = True
    out["sequence_length"] = n
    out["indicators"] = {k: v.to_dict() for k, v in ind.items()}
    out["raw"] = {
        "hydrophobic_patches": apr,
        "cysteines": cys,
        "n_glycosylation": gly,
        "oxidation": ox,
        "deamidation": deam,
        "low_complexity": lc,
    }
    return out


def as_contract_metric_dicts(dev: dict[str, Any]) -> list[dict[str, Any]]:
    """Render computed indicators in the shape of ``contracts.Metric``.

    ``model_name`` is this module, not a predictor, so nothing here can be read as
    a model score or a measurement. None of these names collide with the banned
    affinity names that ``contracts.Metric`` rejects.
    """
    rows: list[dict[str, Any]] = []
    for key, ind in (dev.get("indicators") or {}).items():
        hib = {"better": True, "worse": False}.get(ind["higher_is"])
        rows.append(
            {
                "metric_name": f"seq_{key}",
                "value": None if ind["value"] is None else float(ind["value"]),
                "method": ind["method"],
                "model_name": f"e2b.analysis.{METHOD_VERSION}",
                "model_version": METHOD_VERSION.split("/")[-1],
                "higher_is_better": hib,
                "interpretation": "sequence-derived developability indicator, computed, not measured",
            }
        )
    return rows
