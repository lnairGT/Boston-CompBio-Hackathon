"""Derive the mechanism class, and therefore the functional readout, from inputs.

This is the part that makes the plan MoA-specific rather than generic, and it is the
part most at risk of becoming an indication allowlist. It is not one: the table below
is keyed on **protein-class signalling biology** -- what kind of pathway the target
drives and therefore what a blockade of it would change in a cell. The trigger terms
are protein-family and pathway vocabulary (``interleukin``, ``jak``, ``integrin``,
``protease``), never disease names, and the class is chosen from the target's own
UniProt/Open Targets annotation plus any pathway evidence the spine passes in.

If nothing matches, the class is ``unknown`` and the caller gets an explicit
missingness statement. It does not fall back to a plausible-sounding assay: a
functional assay aimed at the wrong pathway is worse than an admitted gap.
"""

from __future__ import annotations

from typing import Any

from .coerce import text_blob

# --------------------------------------------------------------------------------------
# Mechanism classes
# --------------------------------------------------------------------------------------
# Each entry:
#   label                 human-readable class name
#   triggers              protein-family / pathway vocabulary that selects the class
#   proximal_readout      the mechanism-proximal readout, the thing that moves first
#   readout_because       why that readout reports the mechanism rather than correlating
#   distal_readout        a slower, more physiological confirmation
#   system_requirement    what the cell system must be shown to express
#   stimulus              what drives the pathway in the assay
#   agonist_note          what changes if the requested direction is activation
# --------------------------------------------------------------------------------------

MECHANISM_CLASSES: list[dict[str, Any]] = [
    {
        "key": "cytokine_jak_stat",
        "label": "Cytokine / interleukin receptor signalling (JAK-STAT)",
        "triggers": [
            "interleukin", "cytokine", "jak", "stat", "il-4", "il-13", "il4", "il13",
            "interferon", "type 2 inflammation", "type 2 inflammatory", "th2",
            "common gamma chain", "receptor subunit alpha",
        ],
        "proximal_readout": (
            "receptor-proximal STAT phosphorylation of the STAT isoform named by the "
            "target's own annotated JAK-STAT pathway, by phospho-flow or a quantitative "
            "immunoassay on lysates"
        ),
        "readout_because": (
            "STAT phosphorylation is the first receptor-dependent step after ligand "
            "engagement, so a binder that blocks engagement must reduce it; a downstream "
            "transcript could move for reasons unrelated to this receptor"
        ),
        "distal_readout": (
            "a STAT-response-element reporter line and, separately, an endogenous "
            "pathway-dependent transcript or secreted protein measured in the same cells"
        ),
        "system_requirement": (
            "verified expression of every chain of the receptor complex the target "
            "signals through, confirmed in-house by flow cytometry or qPCR before use"
        ),
        "stimulus": "recombinant target ligand titrated to its own EC50-EC80 in that line",
        "agonist_note": (
            "for an activating direction the stimulus is omitted and the binder itself must "
            "raise the readout above unstimulated baseline"
        ),
    },
    {
        "key": "chemokine_chemotaxis",
        "label": "Chemokine-driven chemotaxis",
        "triggers": ["chemokine", "ccl", "cxcl", "ccr", "cxcr", "chemotaxis", "chemoattractant"],
        "proximal_readout": (
            "receptor-proximal calcium flux or beta-arrestin recruitment in a receptor-"
            "expressing reporter line"
        ),
        "readout_because": (
            "both are immediate consequences of receptor occupancy, so blockade of the "
            "chemokine-receptor interaction must suppress them"
        ),
        "distal_readout": "transwell or real-time migration of primary cells toward the chemokine",
        "system_requirement": "verified expression of the cognate chemokine receptor",
        "stimulus": "recombinant chemokine at its own EC50-EC80",
        "agonist_note": "an activating direction requires demonstrating flux without chemokine",
    },
    {
        "key": "rtk_growth_factor",
        "label": "Receptor tyrosine kinase / growth-factor signalling",
        "triggers": [
            "receptor tyrosine kinase", "growth factor", "egf", "erbb", "her2", "fgf",
            "vegf", "pdgf", "hgf", "insulin receptor", "tyrosine-protein kinase receptor",
            "autophosphorylation",
        ],
        "proximal_readout": "receptor autophosphorylation, then phospho-ERK and phospho-AKT",
        "readout_because": (
            "receptor autophosphorylation is the direct physical consequence of ligand-driven "
            "receptor dimerisation, which is what a blocking binder is meant to prevent"
        ),
        "distal_readout": "ligand-dependent proliferation or colony formation over several days",
        "system_requirement": "verified receptor expression and demonstrated ligand dependence",
        "stimulus": "recombinant cognate growth factor at its own EC50-EC80",
        "agonist_note": "an activating direction must show phosphorylation without ligand",
    },
    {
        "key": "tnfsf_nfkb",
        "label": "TNF-superfamily / NF-kB and death-receptor signalling",
        "triggers": ["tumor necrosis factor", "tnf", "nf-kappa", "nfkb", "trail", "fas ligand", "cd40", "ox40", "baff"],
        "proximal_readout": "IkB degradation or an NF-kB transcriptional reporter",
        "readout_because": "NF-kB activation is the canonical immediate transcriptional consequence of receptor ligation",
        "distal_readout": "caspase activation and viability for death receptors, or cytokine secretion for costimulatory receptors",
        "system_requirement": "verified receptor expression and a demonstrated ligand-dependent response window",
        "stimulus": "recombinant cognate TNF-superfamily ligand at its own EC50-EC80",
        "agonist_note": "an agonist arm must be run with and without receptor clustering, since valency changes the answer",
    },
    {
        "key": "immune_checkpoint",
        "label": "Immune checkpoint / T-cell costimulation",
        "triggers": ["checkpoint", "pd-1", "pdcd1", "pd-l1", "ctla", "lag-3", "tigit", "tim-3", "t-cell activation", "immune recognition"],
        "proximal_readout": "IL-2 or IFN-gamma secretion and activation-marker upregulation in an antigen-specific T-cell co-culture",
        "readout_because": "checkpoint blockade acts by relieving inhibition of TCR signalling, so the mechanism only appears where a TCR signal exists to be relieved",
        "distal_readout": "T-cell proliferation and antigen-specific killing of the partner cell",
        "system_requirement": "a co-culture in which the antigen-presenting partner expresses the checkpoint ligand and the T cells the receptor, both verified",
        "stimulus": "antigen-specific or anti-CD3 TCR stimulation held sub-maximal so a relief of inhibition is measurable",
        "agonist_note": "an agonist direction on an inhibitory receptor would suppress, not raise, the readout, and must be pre-registered as such",
    },
    {
        "key": "protease",
        "label": "Protease activity",
        "triggers": ["protease", "peptidase", "proteinase", "metalloprotease", "cathepsin", "caspase", "convertase"],
        "proximal_readout": "cleavage of a defined substrate, by a fluorogenic peptide assay and by gel or mass-spectrometric detection of the real substrate",
        "readout_because": "loss of catalytic turnover is the mechanism itself, not a proxy for it",
        "distal_readout": "processing of the endogenous substrate in a cell system",
        "system_requirement": "a cell system in which the endogenous substrate is expressed and processed",
        "stimulus": "substrate titration around its own Km",
        "agonist_note": "activation of a protease by an exosite binder must be shown as a change in kcat/Km, not only as endpoint signal",
    },
    {
        "key": "gpcr",
        "label": "G-protein-coupled receptor signalling",
        "triggers": ["g-protein coupled", "gpcr", "seven-transmembrane", "7tm", "adenylate cyclase", "beta-arrestin"],
        "proximal_readout": "cAMP or calcium second-messenger change, plus beta-arrestin recruitment to separate G-protein from arrestin bias",
        "readout_because": "second-messenger change is the immediate transduction step after receptor occupancy",
        "distal_readout": "a pathway-dependent transcriptional reporter",
        "system_requirement": "verified receptor expression and a validated agonist response window",
        "stimulus": "reference agonist at its own EC50-EC80",
        "agonist_note": "for an agonist direction, report both potency and efficacy relative to the reference agonist",
    },
    {
        "key": "adhesion_integrin",
        "label": "Adhesion / integrin-mediated interaction",
        "triggers": ["integrin", "adhesion", "selectin", "cadherin", "icam", "vcam"],
        "proximal_readout": "cell adhesion to the immobilised counter-receptor under defined shear",
        "readout_because": "adhesion is the direct function of the interaction the binder is meant to block",
        "distal_readout": "transmigration or spreading over hours",
        "system_requirement": "verified expression of the adhesion receptor and its counter-receptor on the partner surface",
        "stimulus": "counter-receptor coated surface, with activation state controlled",
        "agonist_note": "activating integrin binders exist and must be distinguished from blockade by a conformation-reporter arm",
    },
    {
        "key": "enzyme_other",
        "label": "Other enzymatic activity",
        "triggers": ["hydrolase", "transferase", "kinase", "phosphatase", "oxidoreductase", "lyase", "ec 1.", "ec 2.", "ec 3."],
        "proximal_readout": "a direct activity assay on purified enzyme with defined substrate and product detection",
        "readout_because": "the mechanism is catalysis, so the readout must be turnover",
        "distal_readout": "the pathway metabolite or substrate pool in cells",
        "system_requirement": "a cell system with measurable basal pathway flux",
        "stimulus": "substrate titration around Km",
        "agonist_note": "activation must be reported as a kinetic parameter change",
    },
]

DIRECTION_SYNONYMS: dict[str, list[str]] = {
    "antagonise": ["antagonis", "antagoniz", "block", "inhibit", "suppress", "neutralis", "neutraliz", "reduce", "dampen"],
    "agonise": ["agonis", "agoniz", "activat", "stimulat", "restore", "enhance", "potentiat", "increase"],
    "degrade": ["degrad", "deplete", "clear"],
}


def normalise_direction(*sources: Any) -> tuple[str, list[str]]:
    """Map free-text direction / desired-effect language onto a canonical direction."""
    blob = text_blob(*sources)
    hits: list[str] = []
    scores: dict[str, int] = {}
    for canon, stems in DIRECTION_SYNONYMS.items():
        matched = [s for s in stems if s in blob]
        if matched:
            scores[canon] = len(matched)
            hits.extend(f"'{m}' -> {canon}" for m in matched)
    if not scores:
        return "unknown", ["no direction vocabulary found in the brief"]
    # "restore anti-tumour immune recognition" by blocking a checkpoint is an
    # antagonist at the molecular level; molecular blockade language wins when both
    # appear, because the assay format follows the molecular action.
    if "antagonise" in scores and "agonise" in scores:
        hits.append(
            "both blockade and activation language present; taking the molecular action as "
            "blockade because the assay format follows what the binder does to the target, "
            "not the organism-level goal"
        )
        return "antagonise", hits
    best = max(scores, key=lambda k: scores[k])
    return best, hits


def classify_mechanism(
    *,
    protein_name: str | None = None,
    gene_symbol: str | None = None,
    keywords: list[str] | None = None,
    pathway_hints: list[str] | None = None,
    desired_effect: str | None = None,
    extra_text: list[str] | None = None,
) -> dict[str, Any]:
    """Pick a mechanism class from target annotation and brief language.

    Target-derived evidence (protein name, UniProt keywords, pathway evidence) is
    weighted above brief language, because the brief says what the user wants and the
    annotation says what the protein does. Returns ``key == "unknown"`` with an
    explicit ``missing`` string when nothing matches.
    """
    target_blob = text_blob(protein_name, gene_symbol, keywords, pathway_hints, extra_text)
    brief_blob = text_blob(desired_effect)

    scored: list[tuple[float, dict[str, Any], list[str]]] = []
    for spec in MECHANISM_CLASSES:
        basis: list[str] = []
        score = 0.0
        for t in spec["triggers"]:
            if t in target_blob:
                score += 2.0
                basis.append(f"target annotation contains '{t}'")
            elif t in brief_blob:
                score += 0.75
                basis.append(f"brief desired_effect contains '{t}'")
        if score > 0:
            scored.append((score, spec, basis))

    if not scored:
        return {
            "key": "unknown",
            "label": "Undetermined",
            "basis": [],
            "missing": (
                "No protein-family or pathway vocabulary in the target's annotation or the "
                "brief matched a known signalling mechanism class. A mechanism-specific "
                "functional assay cannot be specified without knowing which pathway the "
                "target drives. Supply pathway evidence (for example Reactome pathway names "
                "from the target's Open Targets record) and re-run."
            ),
            "alternatives": [],
        }

    scored.sort(key=lambda x: (-x[0], x[1]["key"]))
    top_score, spec, basis = scored[0]
    alts = [{"key": s["key"], "label": s["label"], "score": sc} for sc, s, _ in scored[1:4]]
    ambiguous = len(scored) > 1 and scored[1][0] >= top_score
    out = dict(spec)
    out["basis"] = basis
    out["score"] = top_score
    out["alternatives"] = alts
    out["ambiguous"] = ambiguous
    if ambiguous:
        out["basis"] = basis + [
            f"tie with '{scored[1][1]['key']}' at the same trigger weight; the functional assay "
            "below is written for the reported class and the alternative should be reviewed"
        ]
    return out
