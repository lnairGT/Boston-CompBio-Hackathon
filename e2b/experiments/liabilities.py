"""Sequence-derived construct liabilities for a designed binder.

Everything here is computed from the candidate's own sequence. Nothing is measured,
and the outputs are labelled as computed sequence properties so they cannot be read
as biophysical data. Where the sequence cannot settle a question -- whether two
cysteines actually form a disulfide, for instance -- the screen says what experiment
settles it instead of guessing.
"""

from __future__ import annotations

from typing import Any

from .records import ConstructRisk

# Kyte-Doolittle hydropathy, used only to compute a mean hydropathy index. This is a
# computed sequence property, not a measured solubility.
_KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def sequence_properties(seq: str) -> dict[str, Any]:
    """Computed (not measured) properties of a binder sequence."""
    s = (seq or "").upper()
    if not s:
        return {
            "length": None,
            "n_cys": None,
            "net_charge_ph7_4": None,
            "mean_hydropathy": None,
            "basis": "no sequence supplied",
        }
    pos = s.count("K") + s.count("R") + 0.1 * s.count("H")
    neg = s.count("D") + s.count("E")
    return {
        "length": len(s),
        "n_cys": s.count("C"),
        "net_charge_ph7_4": round(pos - neg, 1),
        "mean_hydropathy": round(sum(_KD.get(a, 0.0) for a in s) / len(s), 2),
        "basis": (
            "computed from the sequence: net charge from K/R/(H x 0.1) minus D/E, hydropathy "
            "as the mean Kyte-Doolittle index. Both are computed properties and neither is a "
            "measured solubility, stability or aggregation propensity."
        ),
    }


def _sequons(s: str) -> list[int]:
    """1-based positions of N in N-X-S/T sequons with X != P."""
    out = []
    for i in range(len(s) - 2):
        if s[i] == "N" and s[i + 1] != "P" and s[i + 2] in ("S", "T"):
            out.append(i + 1)
    return out


def screen_binder_sequence(seq: str | None, *, expression_host: str) -> list[ConstructRisk]:
    """Liabilities that a specific sequence creates in a specific host."""
    risks: list[ConstructRisk] = []
    if not seq:
        return risks
    s = seq.upper()
    host_is_eukaryotic = any(k in expression_host.lower() for k in ("hek", "cho", "mammalian", "insect", "pichia", "yeast"))

    cys = [i + 1 for i, a in enumerate(s) if a == "C"]
    if len(cys) == 1:
        risks.append(
            ConstructRisk(
                risk="single free cysteine",
                evidence=f"one cysteine at position {cys[0]} of {len(s)}",
                positions=cys,
                consequence=(
                    "an unpaired thiol drives intermolecular disulfide formation on air oxidation, "
                    "so the purified protein can present as a covalent dimer. A dimer is bivalent, "
                    "and a bivalent analyte in a binding assay gives an avidity-inflated apparent "
                    "affinity that is not the monovalent affinity of the design."
                ),
                mitigation=(
                    "if the design model does not need the thiol, substitute to serine or alanine and "
                    "re-score the design; if it is kept, purify under reducing conditions, cap the "
                    "thiol, or confirm the monomer fraction by non-reduced SDS-PAGE and SEC on every lot"
                ),
                severity="high",
            )
        )
    elif len(cys) > 1 and len(cys) % 2 == 1:
        risks.append(
            ConstructRisk(
                risk="odd cysteine count: at least one thiol must be unpaired",
                evidence=f"{len(cys)} cysteines at positions {cys}",
                positions=cys,
                consequence=(
                    "at least one free thiol is unavoidable, giving mixed redox species and "
                    "lot-to-lot variability in the binding readout"
                ),
                mitigation=(
                    "determine the free-thiol count experimentally (Ellman assay, or reduced versus "
                    "non-reduced mass) rather than assuming the pairing, and fix the redox state in "
                    "the purification before any affinity measurement"
                ),
                severity="high",
            )
        )
    elif len(cys) >= 2:
        risks.append(
            ConstructRisk(
                risk="cysteines present, pairing unknown from sequence",
                evidence=f"{len(cys)} cysteines at positions {cys}",
                positions=cys,
                consequence=(
                    "an even count does not prove correct pairing; scrambled isomers co-purify and "
                    "are not separable by SDS-PAGE"
                ),
                mitigation=(
                    "confirm the disulfide state by non-reduced versus reduced mass and, if the design "
                    "specifies particular pairs, by peptide mapping. Cytoplasmic E. coli expression "
                    "keeps cysteines reduced, so periplasmic export or an oxidising strain is needed "
                    "if the design depends on a disulfide."
                ),
                severity="moderate",
            )
        )

    seq_pos = _sequons(s)
    if seq_pos:
        risks.append(
            ConstructRisk(
                risk="N-linked glycosylation sequon(s) in the designed sequence",
                evidence=f"N-X-S/T sequon(s) with X != P at N{', N'.join(str(p) for p in seq_pos)}",
                positions=seq_pos,
                consequence=(
                    "occupied in a eukaryotic host, these add heterogeneous glycan mass that can mask "
                    "the designed interface and shifts the observed mass away from the predicted one"
                    if host_is_eukaryotic
                    else "not glycosylated in a prokaryotic host, so they are silent for this material, "
                    "but they become live the moment the construct moves to a mammalian host, an in "
                    "vivo experiment or a fusion format"
                ),
                mitigation=(
                    "if the host is eukaryotic and the sequon is not needed for the design, substitute "
                    "the sequon (N to Q, or S/T to A) and re-score; otherwise record the host dependence "
                    "so the mass and any later host change are interpreted correctly"
                ),
                severity="high" if host_is_eukaryotic else "watch",
            )
        )

    deam = [i + 1 for i in range(len(s) - 1) if s[i] == "N" and s[i + 1] == "G"]
    isom = [i + 1 for i in range(len(s) - 1) if s[i] == "D" and s[i + 1] in ("G", "P")]
    if deam or isom:
        risks.append(
            ConstructRisk(
                risk="chemical-degradation motifs",
                evidence=f"NG at {deam or 'none'}; DG/DP at {isom or 'none'}",
                positions=sorted(deam + isom),
                consequence=(
                    "asparagine deamidation and aspartate isomerisation change charge and can change "
                    "binding over storage, so an affinity measured on aged material may not reproduce"
                ),
                mitigation=(
                    "run a forced-degradation arm (elevated temperature and pH) and re-measure binding "
                    "on stressed material; if binding drops, the motif is engineerable"
                ),
                severity="watch",
            )
        )

    props = sequence_properties(s)
    if props["mean_hydropathy"] is not None and props["mean_hydropathy"] > 0.5:
        risks.append(
            ConstructRisk(
                risk="high computed mean hydropathy",
                evidence=f"mean Kyte-Doolittle index {props['mean_hydropathy']} over {props['length']} residues",
                consequence=(
                    "hydrophobic designs are the ones that aggregate and that stick non-specifically to "
                    "assay surfaces, which shows up as a binding signal that is not binding"
                ),
                mitigation=(
                    "gate on SEC monodispersity before any binding assay and include a non-specific "
                    "binding control surface in the binding experiment"
                ),
                severity="moderate",
            )
        )
    if props["net_charge_ph7_4"] is not None and abs(props["net_charge_ph7_4"]) >= 6:
        risks.append(
            ConstructRisk(
                risk="high computed net charge",
                evidence=f"computed net charge {props['net_charge_ph7_4']} at pH 7.4",
                consequence=(
                    "strongly charged binders bind charged surfaces and matrices electrostatically, "
                    "producing ionic-strength-dependent artefacts in surface-based assays"
                ),
                mitigation=(
                    "titrate ionic strength in the binding buffer and confirm the signal survives a "
                    "physiological salt concentration; include a charge-matched non-binding control"
                ),
                severity="moderate",
            )
        )
    return risks
