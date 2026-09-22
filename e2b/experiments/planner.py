"""Generate an experimental plan from a brief, a target and a candidate.

The plan is *derived*, not templated. Two decisions drive every branch:

``archetype``   whether the target is a secreted ligand, a cell-surface receptor, or
                unreachable by a soluble binder. This sets the antigen reagent, the
                assay orientation and what the functional assay stimulates.
``mechanism``   which signalling class the target drives. This sets the functional
                readout and why that readout reports the mechanism.

Both come from the target's own annotation, so the same code produces a ligand-trap
plan for a secreted cytokine and a receptor-blockade plan for a cell-surface receptor,
and refuses to produce an in-cell plan for a cytosolic target.

What this module will not do: state a concentration, an incubation or a protocol step
as a validated value, report a model score as an affinity, or invent a control that
does not exist. Where a required input is missing, the plan says so and returns a
partial result.
"""

from __future__ import annotations

import hashlib
import json
from typing import Any

from .coerce import as_dict, get
from .liabilities import screen_binder_sequence, sequence_properties
from .moa import normalise_direction
from .records import (
    AssayConcept,
    Control,
    ConstructRisk,
    DecisionCriterion,
    ExperimentalPlan,
    ProductionRoute,
    StartingPoint,
)
from .target_profile import (
    TargetBiologyUnavailable,
    antigen_reagent_plan,
    build_target_profile,
)

SCOPE_STATEMENT = (
    "These are assay CONCEPTS with rationale and controls, not validated protocols. Every "
    "concentration, incubation and surface density below is a starting point that the "
    "responsible scientist must optimise, and no number here has been measured by this "
    "system. Computational binding predictions alone do not demonstrate blockade, "
    "activation, internalisation, safety or efficacy; nothing in this plan, and no result "
    "obtained from it, establishes clinical efficacy or safety."
)

DOES_NOT_ESTABLISH = [
    "Computational binding predictions alone do not demonstrate blockade, activation, "
    "internalisation, safety or efficacy.",
    "A computational model score is not an affinity. No Kd, IC50 or EC50 appears in this plan "
    "as a predicted value, and none may be reported until it is measured.",
    "A positive binding result in one format is not binding; it is one format's signal, which "
    "is why an orthogonal, different-principle confirmation is part of the plan.",
    "In vitro potency is not in vivo efficacy. This plan contains no pharmacokinetic, "
    "immunogenicity, biodistribution or toxicology experiment and says nothing about them.",
    "A designed minibinder is not an antibody; it has no Fc, no effector function and a "
    "different clearance profile, so antibody precedent does not transfer.",
]


def _plan_id(brief: dict, target: dict, candidate: dict) -> str:
    blob = json.dumps(
        {
            "brief": get(brief, "brief_id"),
            "target": get(target, "target_id", "stable_gene_id", "gene_symbol"),
            "candidate": get(candidate, "candidate_id"),
        },
        sort_keys=True,
        default=str,
    )
    return "plan-" + hashlib.sha256(blob.encode()).hexdigest()[:12]


# --------------------------------------------------------------------------------------
# Production
# --------------------------------------------------------------------------------------


def _binder_production_route(
    candidate: dict[str, Any], profile, *, needs_disulfide: bool
) -> ProductionRoute:
    seq = get(candidate, "sequence")
    props = sequence_properties(seq or "")
    n_cys = props.get("n_cys")

    if needs_disulfide:
        host = "E. coli periplasmic secretion (or an oxidising cytoplasmic strain)"
        host_rationale = (
            "the design carries cysteines that are expected to pair, and the reducing E. coli "
            "cytoplasm will not form that bond; periplasmic export provides the oxidising "
            "compartment. A mammalian secretory host is the fallback if folding yield is poor."
        )
    else:
        host = "E. coli cytoplasmic expression"
        host_rationale = (
            "a small, single-domain, non-glycosylated de novo design expresses at high yield in "
            "E. coli, and a bacterial host keeps the protein free of the heterogeneous glycans "
            "that would complicate mass confirmation. The host is chosen from the sequence, not "
            "by default: cysteine pairing or a sequon that must be occupied would move it."
        )

    route = ProductionRoute(
        subject="designed_binder",
        expression_host=host,
        host_rationale=host_rationale,
        construct_elements=[
            "the designed sequence, unmodified in the interface region",
            "a C-terminal (or otherwise interface-distal) purification tag, His6 being the default",
            "a specific protease cleavage site between the design and the tag (TEV or SUMO), so the "
            "tag can be removed",
            "optionally, an interface-distal single biotin-acceptor site or SpyTag for site-specific, "
            "oriented immobilisation in the binding assay",
        ],
        construct_rationale=[
            "the tag is placed away from the designed interface because a tag in or adjacent to the "
            "interface can abolish or, worse, create binding",
            "the tag is cleavable because an uncleaved tag is itself a binding epitope for the capture "
            "reagent, and because the tagged and untagged protein must give the same answer",
            "a site-specific immobilisation handle is preferred over amine coupling: amine coupling "
            "reacts with lysines anywhere, including in the designed interface, and produces a "
            "heterogeneous surface whose affinity is not the molecule's affinity",
        ],
        purification_steps=[
            "immobilised metal affinity chromatography (IMAC) on clarified lysate",
            "tag cleavage, then a subtractive IMAC pass to remove the tag, the protease and uncleaved material",
            "preparative size-exclusion chromatography as the final polish, in the assay buffer",
            "if the sequence screen flagged free thiols: a defined redox step before the final column, "
            "so the redox state is fixed rather than drifting",
        ],
        release_qc=[
            "analytical SEC: retention consistent with a monomer, single symmetric peak. Oligomeric "
            "state is the first gate because a multimer binds multivalently and reports an avidity, "
            "not an affinity",
            "SDS-PAGE, reduced and non-reduced: the non-reduced lane is what reveals a disulfide-linked "
            "dimer that the reduced lane hides",
            "mass confirmation by ESI-MS against the calculated mass of the intended construct, which "
            "is the only check that the protein in the tube is the sequence that was designed",
            "concentration by A280 using the extinction coefficient computed from the sequence, because "
            "every affinity number is only as accurate as the analyte concentration",
            "endotoxin by LAL before any cell-based assay: bacterial endotoxin carried through from an "
            "E. coli prep activates innate signalling and will be read as a functional effect",
        ],
        qc_gates=[
            "no binding assay result is believed before analytical SEC, SDS-PAGE and mass confirmation pass",
            "no cell-based assay result is believed before the endotoxin result passes",
            "a lot that fails QC is re-made, not carried forward with a caveat",
        ],
        starting_points=[
            StartingPoint(
                name="protein needed for the full plan",
                range_text="roughly 1-5 mg of the design and of each control protein, more if a "
                "calorimetric orthogonal method is chosen",
                reasoning="surface-based binding and cell assays consume little; solution calorimetry "
                "is the material-hungry step, so the quantity is set by which orthogonal method is chosen",
            ),
            StartingPoint(
                name="storage and handling",
                range_text="a neutral-pH buffer at physiological ionic strength, single-use aliquots, "
                "with a freeze-thaw and a short-term stability arm run once",
                reasoning="designed miniproteins are usually robust, but binding measured on "
                "repeatedly thawed material is a common source of irreproducible affinities",
            ),
        ],
        risks=screen_binder_sequence(seq, expression_host=host),
        notes=[
            f"computed sequence properties (not measurements): {props['basis']}",
            f"length {props['length']}, cysteines {n_cys}, computed net charge {props['net_charge_ph7_4']} "
            f"at pH 7.4, mean Kyte-Doolittle hydropathy {props['mean_hydropathy']}"
            if props["length"]
            else "no candidate sequence supplied, so no sequence-derived liability was screened",
        ],
    )
    return route


def _antigen_production_route(profile, protein: dict[str, Any] | None) -> ProductionRoute | None:
    reagent = antigen_reagent_plan(profile, protein)
    if not reagent.get("construct"):
        return None
    span = reagent.get("span") or {}
    span_text = (
        f"residues {span.get('start')}-{span.get('end')} ({span.get('kind')})"
        if span.get("start")
        else "boundaries not resolved from the topology annotation"
    )
    risks: list[ConstructRisk] = []
    if (profile.n_glycosylation_sites or 0) > 0:
        risks.append(
            ConstructRisk(
                risk="glycosylation on the target reagent",
                evidence=f"{profile.n_glycosylation_sites} annotated N-glycosylation site(s) on the target",
                consequence=(
                    "glycans add heterogeneous mass and can shield the designed epitope, so a design "
                    "that binds a bacterially produced, unglycosylated antigen may not bind the native "
                    "glycoprotein"
                ),
                mitigation=(
                    "produce the antigen in a mammalian host and, if the epitope sits near a sequon, "
                    "compare binding on glycosylated and enzymatically deglycosylated antigen"
                ),
                severity="high",
            )
        )
    if profile.n_disulfides:
        risks.append(
            ConstructRisk(
                risk="disulfide-dependent fold on the target reagent",
                evidence=f"{profile.n_disulfides} annotated disulfide bond(s)",
                consequence="a misfolded antigen lot gives both false negatives and false positives",
                mitigation="confirm the lot binds its native partner before using it to judge a design",
                severity="high",
            )
        )
    return ProductionRoute(
        subject="target_antigen",
        expression_host=reagent["host"],
        host_rationale=reagent["host_rationale"],
        construct_elements=[
            f"{reagent['construct']}: {span_text}",
            "an interface-distal tag for oriented capture (His, biotin-acceptor or Fc), chosen so the "
            "capture chemistry does not touch the region the design engages",
            "for a receptor: a monomeric tagged form AND, if an Fc fusion is used, an explicit note "
            "that the Fc form is bivalent",
        ],
        construct_rationale=reagent.get("reasons", []),
        purification_steps=[
            "affinity capture on the tag, then preparative SEC in the assay buffer",
            "for a glycoprotein, confirm homogeneity by analytical SEC rather than assuming it",
        ],
        release_qc=[
            "analytical SEC for oligomeric state, SDS-PAGE, and mass confirmation",
            "an activity check: "
            + (reagent.get("activity_requirement") or "the lot must be shown to be functional"),
        ],
        qc_gates=[
            "a binding result against an antigen lot that has not been shown to be active is not "
            "evidence about the design; it is evidence about the lot"
        ],
        risks=risks,
        notes=[profile.source_notes[-1]] if profile.source_notes else [],
    )


# --------------------------------------------------------------------------------------
# Binding validation
# --------------------------------------------------------------------------------------


def _binding_assay(profile, candidate: dict[str, Any]) -> AssayConcept:
    cand_len = get(candidate, "sequence_length") or len(get(candidate, "sequence") or "") or None
    target_len = profile.length_aa

    common_controls = [
        Control(
            name="sequence-matched non-binding negative",
            kind="negative",
            purpose=(
                "a scrambled or interface-ablated variant of the same sequence, at the same length and "
                "amino-acid composition, expressed and purified by the same route with the same tag. It "
                "is the only control that separates designed binding from the behaviour of a protein of "
                "this composition on this surface"
            ),
            interpretation_if_it_fails=(
                "if the negative gives the same signal, the signal is not designed binding: it is "
                "non-specific surface or matrix interaction, and the primary result is void"
            ),
            available=True,
            availability_note="must be made alongside the design, not sourced later",
        ),
        Control(
            name="reference binder positive",
            kind="positive",
            purpose=(
                "a molecule already known to bind this target -- its native binding partner is usually "
                "the best available choice because its affinity is documented -- run in the same "
                "orientation to show the assay can detect binding at all"
            ),
            interpretation_if_it_fails=(
                "if the reference does not bind, the immobilised partner is inactive or the surface is "
                "wrong, and a negative result on the design is uninterpretable rather than a failure of "
                "the design"
            ),
            available=None,
            availability_note=(
                "availability is target-specific and is not asserted here. Confirm whether a validated "
                "binder or a recombinant native partner exists for this target before the run; if none "
                "exists, the antigen-activity control below carries this burden and the limitation is "
                "recorded"
            ),
        ),
        Control(
            name="buffer, blank and reference-channel subtraction",
            kind="buffer_or_system",
            purpose=(
                "buffer-only injections on the active surface and identical injections on a reference "
                "surface, double-referenced, to remove bulk refractive-index shift and drift"
            ),
            interpretation_if_it_fails=(
                "uncorrected bulk shifts and drift are routinely mistaken for fast, weak binding"
            ),
            available=True,
        ),
        Control(
            name="surface stability and regeneration control",
            kind="buffer_or_system",
            purpose=(
                "one concentration repeated at the start, middle and end of the series, with the "
                "regeneration condition held constant"
            ),
            interpretation_if_it_fails=(
                "a drifting repeat injection means the surface is losing activity, so the fitted "
                "kinetics describe a changing surface and not the interaction"
            ),
            available=True,
        ),
        Control(
            name="analyte concentration accuracy",
            kind="reagent_activity",
            purpose=(
                "analyte concentration determined by A280 with the sequence-computed extinction "
                "coefficient (amino-acid analysis if aromatic content is too low for A280)"
            ),
            interpretation_if_it_fails=(
                "an affinity is only as accurate as the concentration axis; a systematically wrong "
                "concentration produces a systematically wrong Kd with a convincing-looking fit"
            ),
            available=True,
        ),
    ]

    conc = StartingPoint(
        name="analyte concentration series",
        range_text=(
            "a starting bracket spanning at least two and preferably three orders of magnitude, in "
            "half-log steps, with the top concentration high enough to approach saturation; re-bracket "
            "empirically after the first pass so the series straddles the observed midpoint"
        ),
        reasoning=(
            "a titration resolves an affinity only inside the range it covers and only if the top "
            "concentration approaches saturation. The design's model score does NOT predict an "
            "affinity and must not be used to set this bracket, so the first pass is deliberately "
            "wide and is narrowed by data, not by expectation"
        ),
    )
    resolvable = StartingPoint(
        name="what the format can resolve",
        range_text=(
            "surface plasmon resonance resolves association rates up to roughly 1e7 M-1 s-1 and "
            "dissociation rates down to roughly 1e-5 s-1; outside that window a kinetic fit is not "
            "trustworthy and a steady-state (equilibrium) affinity should be reported instead, which "
            "requires reaching saturation"
        ),
        reasoning=(
            "instrument limits, not the analysis software, set the range in which a kinetic constant "
            "means anything. Very tight binders hit the dissociation-rate floor and very weak ones "
            "dissociate within the injection"
        ),
    )

    if profile.archetype == "secreted_ligand":
        orientation = (
            "Immobilise the DESIGNED BINDER through its interface-distal site-specific handle "
            "(streptavidin capture of a single biotin, or anti-tag capture); flow the recombinant "
            "target ligand as the analyte."
        )
        orientation_rationale = (
            "Three reasons, in order of weight. (1) Both partners are monomeric here"
            + (f" (design {cand_len} aa, target {target_len} aa)" if cand_len and target_len else "")
            + ", so a 1:1 interaction model is physically justified in this orientation and no avidity "
            "term is introduced. (2) Surface-plasmon response scales with analyte mass, and the "
            "secreted target is the heavier partner, so flowing the target gives the larger response "
            "per binding event and a better signal-to-noise at low occupancy. (3) The target ligand is "
            "the reagent whose native fold matters most and is the harder one to replace; keeping it in "
            "solution avoids immobilisation chemistry on the molecule whose activity the whole "
            "experiment depends on. Capture rather than amine coupling is used for the design so that "
            "lysines in the designed interface are not chemically modified."
        )
        density = StartingPoint(
            name="immobilisation density",
            range_text="the lowest capture level that still gives a usable maximum response, deliberately "
            "well below a saturated surface",
            reasoning="high surface density causes mass-transport limitation and analyte rebinding, both "
            "of which flatten the apparent dissociation and make a binder look tighter than it is",
        )
    else:
        orientation = (
            "Capture the TARGET ECTODOMAIN on the surface through its tag (anti-His, streptavidin-biotin "
            "or anti-Fc capture, never random amine coupling); flow the monomeric designed binder as the "
            "analyte."
        )
        orientation_rationale = (
            "Three reasons. (1) Receptor ectodomain reagents are frequently dimeric, either natively or "
            "because they are Fc-fused; a bivalent analyte binds a surface bivalently and returns an "
            "avidity-inflated apparent affinity that is not the monovalent affinity of the design. "
            "Keeping the potentially bivalent partner on the surface and the strictly monomeric design "
            "in solution preserves a 1:1 interaction. (2) Capture chemistry, rather than amine coupling, "
            "keeps the glycosylated, disulfide-bonded ectodomain in its folded state and avoids "
            "modifying lysines in the epitope the design was built against. (3) The captured surface can "
            "be rebuilt between cycles, so a fresh, demonstrably active antigen presents on every cycle "
            "instead of one surface degrading across the titration. The cost of this orientation is a "
            "small analyte mass and therefore a low maximum response, which is why the reciprocal "
            "orientation is run as a cross-check below."
        )
        density = StartingPoint(
            name="capture density",
            range_text="the lowest antigen capture level that yields a detectable maximum response for a "
            "small analyte; expect to trade signal against artefact and to report which way that was traded",
            reasoning="a small analyte forces higher ligand density for signal, which is exactly the "
            "condition that creates rebinding artefacts, so the trade-off is explicit and reported rather "
            "than hidden",
        )

    return AssayConcept(
        assay_id="bind-primary",
        role="binding_primary",
        title="Primary binding validation by surface plasmon resonance",
        format="Surface plasmon resonance (SPR), single-cycle or multi-cycle kinetics with a "
        "reference-subtracted, double-referenced series. Biolayer interferometry (BLI) is an "
        "acceptable substitute where SPR is unavailable.",
        principle="label-free detection of mass accumulating at an optical sensor surface",
        rationale=(
            "A surface-based kinetic method is the right primary because it reports association and "
            "dissociation separately, and because a designed binder's dissociation rate is the property "
            "that matters most and the one an endpoint assay cannot see."
        ),
        orientation=orientation,
        orientation_rationale=orientation_rationale,
        readout="association and dissociation traces, fitted to a 1:1 interaction model, reported with "
        "residuals and with the fitted maximum response compared against the theoretical maximum",
        readout_reports_mechanism_because=(
            "this assay reports binding only. It does not report blockade, and it must not be described "
            "as showing inhibition"
        ),
        system="purified recombinant proteins in a defined buffer, at controlled temperature",
        system_requirement="both partners QC-passed, monomeric where claimed, and the antigen lot shown "
        "to be active",
        resolvable_range=resolvable.range_text,
        resolvable_range_reasoning=resolvable.reasoning,
        starting_points=[conc, density, resolvable],
        controls=common_controls
        + [
            Control(
                name="reciprocal-orientation cross-check",
                kind="specificity",
                purpose="repeat the measurement with the partners swapped; the two affinities should agree "
                "within a few fold",
                interpretation_if_it_fails="a large disagreement between orientations means one "
                "orientation carries an artefact (avidity, immobilisation damage or mass-transport), and "
                "neither number should be quoted until it is resolved",
                available=True,
            )
        ],
        reports=[
            "whether a concentration-dependent, saturable interaction exists",
            "association and dissociation rate constants and the equilibrium dissociation constant, "
            "within the instrument's resolvable window",
        ],
        does_not_report=[
            "blockade, agonism or any functional consequence",
            "binding in a cellular context, where the target is glycosylated, clustered and surrounded "
            "by other surface proteins",
            "stoichiometry, which needs a solution method",
        ],
        derived_from=[
            f"target archetype: {profile.archetype}",
            f"archetype basis: {'; '.join(profile.archetype_basis[:2]) if profile.archetype_basis else 'none'}",
        ],
        limitations=[
            "a single-format binding signal is not binding; the orthogonal confirmation below is part of "
            "the claim, not an optional extra",
            "immobilisation, however gentle, is a perturbation",
        ],
    )


def _orthogonal_assay(profile) -> AssayConcept:
    return AssayConcept(
        assay_id="bind-orthogonal",
        role="binding_orthogonal",
        title="Orthogonal, immobilisation-free confirmation",
        format="Two solution methods of different physical principle: (a) mass photometry or SEC-MALS "
        "for complex formation and stoichiometry; (b) isothermal titration calorimetry for a "
        "surface-free affinity where material allows, or a fluorescence-based solution competition "
        "assay where it does not.",
        principle="single-molecule mass measurement / light scattering / calorimetry -- none of which "
        "involves attaching either partner to a surface",
        rationale=(
            "Biolayer interferometry is NOT an orthogonal confirmation of surface plasmon resonance: "
            "both immobilise a partner and both detect mass at an optical interface, so they share the "
            "artefacts that most often produce a false positive. A genuine cross-check has to remove "
            "immobilisation from the experiment. Mass photometry or SEC-MALS additionally answers a "
            "question the surface method cannot: whether the complex is 1:1 or something else."
        ),
        orientation="none: both partners free in solution",
        orientation_rationale="removing the surface removes avidity, mass-transport and "
        "coupling-chemistry artefacts as explanations for the primary signal",
        readout="complex mass and stoichiometry; binding enthalpy and affinity, or competition-shifted "
        "signal",
        readout_reports_mechanism_because="it still reports binding only, but binding that survives the "
        "removal of the surface",
        system="purified proteins in solution",
        system_requirement="both lots QC-passed; calorimetry additionally needs matched buffers, since "
        "buffer mismatch produces heat that looks like binding",
        resolvable_range="calorimetry is practical for roughly nanomolar to high-micromolar affinities "
        "and needs substantially more material; mass photometry reports stoichiometry and "
        "complex formation rather than a precise affinity",
        resolvable_range_reasoning="each method has a window, and choosing one commits to what can be "
        "concluded from it",
        starting_points=[
            StartingPoint(
                name="titration design",
                range_text="a concentration series centred on the affinity estimated from the primary "
                "assay, with the stoichiometry arm run near and above the estimated affinity",
                reasoning="the primary assay's estimate sets the useful range; running blind wastes "
                "material on a method that consumes the most of it",
            )
        ],
        controls=[
            Control(
                name="same sequence-matched negative",
                kind="negative",
                purpose="the negative from the primary assay, run in the orthogonal format",
                interpretation_if_it_fails="signal shared with the negative in both formats means a "
                "non-specific interaction, not a design that works",
                available=True,
            ),
            Control(
                name="buffer-mismatch and dilution-heat control",
                kind="buffer_or_system",
                purpose="titrate buffer into protein and protein into buffer",
                interpretation_if_it_fails="dilution and mismatch heats are routinely misfitted as a "
                "weak binding isotherm",
                available=True,
            ),
        ],
        reports=["whether the complex forms without a surface", "the stoichiometry of the complex"],
        does_not_report=["any functional consequence", "behaviour on a cell"],
        derived_from=[f"target archetype: {profile.archetype}"],
        limitations=[
            "agreement between two formats raises confidence in binding; it still says nothing about "
            "mechanism"
        ],
    )


# --------------------------------------------------------------------------------------
# Functional MoA
# --------------------------------------------------------------------------------------


def _functional_assay(profile, mech: dict[str, Any], direction: str, brief: dict[str, Any]) -> AssayConcept | None:
    if mech.get("key") == "unknown":
        return None

    desired = get(brief, "desired_effect") or "the effect named in the brief"
    is_antagonist = direction == "antagonise"

    if profile.archetype == "secreted_ligand":
        design = (
            "The target is a secreted ligand, so this is a ligand-trap experiment. Responder cells "
            "expressing the target's RECEPTOR (not the target) are stimulated with the recombinant "
            "target ligand held at a fixed, sub-maximal concentration near its own half-maximal "
            "response, and the designed binder is titrated against it. Run two addition orders -- "
            "binder pre-incubated with the ligand before addition, and binder added to cells "
            "simultaneously with the ligand -- because a sequestering trap and a receptor-competitive "
            "blocker behave differently between the two, and the difference is mechanistic information."
        )
        system = (
            "a cell line or primary cell population with verified expression of every chain of the "
            "receptor complex the target signals through"
        )
        extra_control = Control(
            name="ligand dose-response shift (competition analysis)",
            kind="specificity",
            purpose="repeat the full ligand dose-response at two or three fixed binder concentrations; a "
            "competitive trap shifts the ligand curve rightward without lowering the maximum response",
            interpretation_if_it_fails="a suppressed maximum rather than a rightward shift indicates "
            "something other than simple competition -- non-competitive interference, cytotoxicity or an "
            "assay artefact -- and must be resolved before the result is called blockade",
            available=True,
        )
    elif profile.archetype.startswith("cell_surface_receptor"):
        design = (
            "The target is a cell-surface receptor, so this is a receptor-blockade experiment. Cells "
            "expressing the TARGET receptor are stimulated with its cognate ligand at a fixed, "
            "sub-maximal concentration, and the designed binder is titrated against it. Include a "
            "binder-alone arm without ligand: a binder against a receptor can be an unintended agonist, "
            "especially if it clusters or stabilises an active conformation, and that has to be measured "
            "rather than assumed absent."
        )
        system = (
            "a cell line or primary population with verified surface expression of the target receptor "
            "and a demonstrated ligand-dependent response window"
        )
        extra_control = Control(
            name="binder-alone (unintended agonism) arm",
            kind="specificity",
            purpose="the full binder titration with no stimulating ligand present",
            interpretation_if_it_fails="a signal in the binder-alone arm means the molecule is an "
            "agonist or partial agonist, which for a blocking programme is a kill or a redesign, not a "
            "footnote",
            available=True,
        )
    else:
        return None

    if not is_antagonist:
        design += " " + (mech.get("agonist_note") or "")

    return AssayConcept(
        assay_id="func-moa",
        role="functional_moa",
        title=f"Mechanism-specific functional assay: {mech['label']}",
        format="cell-based, dose-response, with a mechanism-proximal readout and an independent "
        "downstream confirmation",
        principle=mech["label"],
        rationale=(
            f"The brief asks to {direction} the target in order to {desired}. This assay is written "
            f"against the mechanism the target actually drives -- {mech['label']} -- because a generic "
            "viability or reporter readout would move for dozens of reasons unrelated to this target, "
            "and could not distinguish blockade from toxicity."
        ),
        orientation=design,
        orientation_rationale=(
            "stimulus and cell system follow the archetype: a secreted target is the stimulus and the "
            "receptor-bearing cell is the responder, whereas a receptor target is on the responder and "
            "its own ligand is the stimulus"
        ),
        readout=mech["proximal_readout"],
        readout_reports_mechanism_because=mech["readout_because"],
        system=system,
        system_requirement=(
            mech["system_requirement"]
            + ". Expression must be confirmed in the actual cells used, in-house, before the experiment: "
            "a published expression claim about a line is not a property of your vial. Where cellular "
            "evidence identified the disease-relevant population, a primary-cell confirmation in that "
            "population is the stronger version of this experiment."
        ),
        resolvable_range=(
            "a dose-response spanning several orders of magnitude around the affinity measured in the "
            "binding assay, in half-log steps, with enough points above and below the midpoint to fit a "
            "four-parameter curve and enough replicates to state a confidence interval on the midpoint"
        ),
        resolvable_range_reasoning=(
            "the functional midpoint is expected near or above the measured affinity; a functional "
            "potency far tighter than the measured affinity is a red flag for an artefact rather than a "
            "good result"
        ),
        starting_points=[
            StartingPoint(
                name="stimulus level",
                range_text=mech["stimulus"]
                + ", established by running the stimulus dose-response in these cells first",
                reasoning="a saturating stimulus makes any competitive blocker look inactive; a "
                "sub-maximal stimulus is what gives the assay a window to inhibit",
            ),
            StartingPoint(
                name="binder titration",
                range_text="8-12 points, half-log spacing, bracketing the measured affinity by at least "
                "two logs in each direction, in at least triplicate on independent days",
                reasoning="a midpoint quoted from a single plate is not a potency; day-to-day variation "
                "in cell state is usually larger than the fit error",
            ),
            StartingPoint(
                name="incubation and readout timing",
                range_text="for a phosphorylation readout, minutes after stimulation; for a "
                "transcriptional or secreted readout, hours; the exact time must be found from a "
                "time-course in these cells",
                reasoning="proximal signalling is transient, so a readout taken at the wrong time "
                "reports the decay of the signal rather than its magnitude",
            ),
        ],
        controls=[
            Control(
                name="sequence-matched non-binding negative, same lots",
                kind="negative",
                purpose="the scrambled or interface-ablated variant carried through the same production "
                "and the same endotoxin release test",
                interpretation_if_it_fails="an effect from the negative points at endotoxin, buffer or "
                "protein load rather than at the designed interaction",
                available=True,
            ),
            Control(
                name="pathway-positive control inhibitor",
                kind="positive",
                purpose="a molecule known to block this pathway at or below the receptor, to prove the "
                "assay can register inhibition on this plate",
                interpretation_if_it_fails="without a working positive, a flat binder curve cannot be "
                "distinguished from an assay that was never able to show inhibition",
                available=None,
                availability_note="a pathway-level inhibitor usually exists for a characterised "
                "signalling class, but existence for this specific pathway is not asserted here and must "
                "be confirmed",
            ),
            Control(
                name="target-negative cell system",
                kind="specificity",
                purpose="the same experiment in cells that do not express the target (or in which it is "
                "knocked out), which is the cleanest on-target demonstration available in vitro",
                interpretation_if_it_fails="an effect that persists without the target is off-target, "
                "whatever the binding data show",
                available=True,
            ),
            Control(
                name="excess soluble target rescue",
                kind="specificity",
                purpose="pre-absorb the binder with excess recombinant target; the functional effect "
                "should disappear",
                interpretation_if_it_fails="a non-rescuable effect is not mediated by target engagement",
                available=True,
            ),
            Control(
                name="viability and endotoxin confound control",
                kind="buffer_or_system",
                purpose="viability measured in parallel at every dose, with the endotoxin release result "
                "on record for each lot",
                interpretation_if_it_fails="cytotoxicity reads as inhibition of every pathway at once, "
                "and endotoxin reads as activation of innate signalling",
                available=True,
            ),
            extra_control,
        ],
        reports=[
            "whether target engagement produces the pathway-level consequence the brief asks for",
            "a functional potency, with a confidence interval, in this cell system",
        ],
        does_not_report=[
            "efficacy in disease, in an animal or in a human",
            "durability, exposure, immunogenicity or safety",
            "an effect in the disease-relevant primary cell type, unless that is the system used",
        ],
        derived_from=[
            f"mechanism class: {mech['key']} ({mech['label']})",
            "class basis: " + ("; ".join(mech.get("basis") or []) or "none"),
            f"direction: {direction}",
            f"desired effect from brief: {desired}",
        ],
        limitations=[
            "a downstream confirmation is required alongside the proximal readout: "
            + mech["distal_readout"],
            "a cell line is not the disease; the assay establishes mechanism, not benefit",
        ]
        + (["mechanism class was ambiguous between two classes; review before running"] if mech.get("ambiguous") else []),
    )


def _specificity_assay(profile, paralogs: list[dict[str, Any]] | None) -> AssayConcept:
    paralogs = paralogs or []
    have = bool(paralogs)
    names = ", ".join(
        str(p.get("gene_symbol") or p.get("accession") or p.get("target_id")) for p in paralogs[:8]
    )
    return AssayConcept(
        assay_id="spec-panel",
        role="specificity",
        title="Specificity and selectivity arm",
        format="binding counter-screen against a paralog panel, plus a polyspecificity screen, plus the "
        "cellular on-target controls from the functional assay",
        principle="a binder is only useful if what it does NOT bind is also known",
        rationale=(
            "Selectivity has to be measured against the proteins most likely to be cross-bound -- the "
            "target's closest paralogs, which share the fold and often the epitope surface -- not "
            "against an arbitrary panel of unrelated proteins, which is an easy test to pass."
        ),
        orientation="same format and orientation as the primary binding assay, so the numbers are "
        "comparable; each paralog lot carries its own activity control",
        orientation_rationale="a selectivity ratio computed across two different assay formats is not a "
        "selectivity ratio",
        readout="affinity, or an explicit lower bound where no binding is detected, for every panel "
        "member; selectivity reported as a ratio with the bound stated",
        readout_reports_mechanism_because="cross-reactivity to a paralog that carries a different "
        "function is the most likely route to an unintended effect",
        system="purified paralog ectodomains or mature proteins, produced by the route the archetype "
        "requires",
        system_requirement=(
            f"panel supplied by the paralog analysis: {names}"
            if have
            else "NO paralog panel was supplied. The panel must come from a sequence and family "
            "analysis of this target before this arm can be specified; it is not guessed here"
        ),
        resolvable_range="where no binding is seen, report the highest concentration tested as a "
        "lower bound on the dissociation constant; 'no binding' without a bound is not a result",
        resolvable_range_reasoning="an undetected interaction is a limit of detection, not a zero",
        starting_points=[
            StartingPoint(
                name="panel composition",
                range_text="the closest paralogs by sequence identity over the engaged region, plus any "
                "family member with a known opposing function",
                reasoning="identity over the region the binder engages predicts cross-reactivity far "
                "better than whole-protein identity",
            ),
            StartingPoint(
                name="polyspecificity screen",
                range_text="a small panel of unrelated proteins and a cell-surface polyspecificity "
                "screen, run at a concentration well above the target affinity",
                reasoning="designed miniproteins can be sticky; a polyspecific molecule produces "
                "reproducible, meaningless signals in later assays",
            ),
        ],
        controls=[
            Control(
                name="paralog lot activity control",
                kind="reagent_activity",
                purpose="each paralog preparation shown to bind its own native partner",
                interpretation_if_it_fails="apparent selectivity against an inactive paralog lot is not "
                "selectivity; it is a failed reagent",
                available=True,
            )
        ],
        reports=["selectivity over the closest paralogs, with bounds", "general stickiness"],
        does_not_report=["proteome-wide selectivity, which this panel does not sample"],
        derived_from=(
            ["paralog panel supplied by the downstream analysis lane"]
            if have
            else ["no paralog panel available at plan time"]
        ),
        limitations=(
            []
            if have
            else [
                "PROVISIONAL: without a paralog panel this arm names a method but not its members, and "
                "any selectivity statement from it would be unsupported"
            ]
        ),
    )


# --------------------------------------------------------------------------------------
# Decision tree
# --------------------------------------------------------------------------------------


def _decision_tree(profile, mech: dict[str, Any], direction: str, has_functional: bool) -> list[DecisionCriterion]:
    tree = [
        DecisionCriterion(
            criterion_id="G0",
            stage="Production and QC",
            question="Is the material in the tube the designed molecule, monomeric, and clean enough to "
            "believe an assay run on it?",
            pass_rule="a single symmetric analytical-SEC peak at monomer retention with the great "
            "majority of the mass in it; observed mass matching the calculated mass of the intended "
            "construct within the instrument's stated accuracy; no disulfide-linked species in the "
            "non-reduced gel; endotoxin below the assay's threshold if cell work follows",
            kill_rule="mass does not match the intended sequence after re-sequencing the construct, or "
            "the material is irreducibly multimeric or aggregated across independent preparations",
            rule_basis="oligomeric state and identity are prerequisites, not results: a multimer reports "
            "avidity and a mass mismatch means the molecule tested is not the molecule designed",
            on_pass="proceed to G1",
            on_kill="stop. Report as a production failure of this design, not as a binding failure, and "
            "do not report any affinity from this material",
            on_ambiguous="a partially aggregated prep is re-purified or re-made; it is not carried "
            "forward with a caveat",
        ),
        DecisionCriterion(
            criterion_id="G1",
            stage="Primary binding",
            question="Is there a concentration-dependent, saturable interaction that the "
            "sequence-matched negative does not show?",
            pass_rule="a dose-dependent signal that saturates, fits a 1:1 model with acceptable "
            "residuals, agrees within a few fold between the two orientations, and exceeds the "
            "sequence-matched negative by a margin fixed before the run",
            kill_rule="no signal at the top concentration in both orientations WHILE the antigen-activity "
            "control and, where available, the reference-binder positive both work",
            rule_basis="the kill is conditional on the controls working, because a dead surface produces "
            "the same trace as a non-binder and the two must not be confused",
            on_pass="proceed to G2",
            on_kill="kill the candidate for binding. This is a real, reportable negative result about the "
            "design",
            on_ambiguous="if the positive or activity control failed, the result is INCONCLUSIVE: fix the "
            "reagent and repeat. An inconclusive run is never recorded as a pass",
        ),
        DecisionCriterion(
            criterion_id="G2",
            stage="Orthogonal confirmation",
            question="Does the complex form with no surface involved, at the expected stoichiometry?",
            pass_rule="complex detected by the solution method at a concentration consistent with the "
            "affinity from G1, with stoichiometry consistent with the design intent",
            kill_rule="no complex in solution despite a clean surface-based signal, once buffer-mismatch "
            "and dilution-heat controls are accounted for",
            rule_basis="a signal that exists only when a partner is immobilised is most often an "
            "immobilisation artefact, and publishing it as binding is the error this gate exists to catch",
            on_pass="proceed to G3",
            on_kill="kill the binding claim. Record that the surface signal did not reproduce without a "
            "surface",
            on_ambiguous="a stoichiometry other than the intended one is not automatically a kill, but it "
            "changes the mechanism being tested and the functional design must be revisited first",
        ),
    ]

    if has_functional:
        tree.append(
            DecisionCriterion(
                criterion_id="G3",
                stage="Mechanism-specific function",
                question=f"Does the design {direction} the pathway through the mechanism-proximal readout "
                f"({mech['label']}), in a dose-dependent way?",
                pass_rule="a monotonic dose-response with a fitted midpoint and confidence interval, a "
                "clear effect size at the top dose, reproduced on independent days, with the pathway "
                "positive control working, no effect from the sequence-matched negative, no viability "
                "loss at the active doses, and the effect abolished in target-negative cells and by "
                "pre-absorption with excess target",
                kill_rule="binding confirmed at G1 and G2, but no dose-dependent functional effect at "
                "concentrations two or more logs above the measured affinity, while the pathway positive "
                "control does produce the expected inhibition",
                rule_basis="a design can bind a real epitope that does not participate in signalling. "
                "Binding without function at saturating occupancy is the signature of a non-functional "
                "epitope, and the honest conclusion is that this epitope is the wrong one",
                on_pass="proceed to G4",
                on_kill="kill this design for this mechanism. The evidence supports redesign against a "
                "different epitope, not a re-interpretation of the functional data",
                on_ambiguous="a partial effect that plateaus below the positive control is reported as "
                "partial, with the plateau value, and is not rounded up to blockade",
            )
        )
        tree.append(
            DecisionCriterion(
                criterion_id="G4",
                stage="Specificity",
                question="Is the effect selective enough over the closest paralogs to be worth pursuing?",
                pass_rule="a selectivity margin over each panel paralog that was fixed before the run, "
                "with bounds reported where no binding was detected, and no polyspecific behaviour",
                kill_rule="comparable potency on a paralog whose inhibition is undesirable, or a "
                "functional effect that persists in target-negative cells",
                rule_basis="selectivity requirements depend on paralog biology and must be set with the "
                "target's family in view, before the data exist",
                on_pass="candidate is supported for the next stage of characterisation, which this plan "
                "does not cover",
                on_kill="kill or redesign toward a target-distinguishing epitope",
                on_ambiguous="if no paralog panel exists, G4 cannot be evaluated and the candidate "
                "carries an explicit, unresolved selectivity risk. It is not passed by default",
            )
        )
    else:
        tree.append(
            DecisionCriterion(
                criterion_id="G3",
                stage="Mechanism-specific function",
                question="NOT SPECIFIABLE: the mechanism class could not be derived, so no functional "
                "gate can be pre-registered",
                pass_rule="none; this gate is undefined until the target's signalling mechanism is known",
                kill_rule="none",
                rule_basis="a functional gate written against an unknown pathway would be a guess, and a "
                "guessed criterion is worse than a stated gap because it looks pre-registered",
                on_pass="",
                on_kill="",
                on_ambiguous="obtain pathway evidence for this target and regenerate the plan",
            )
        )
    return tree


# --------------------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------------------


def plan_experiments(
    brief: Any,
    target: Any,
    candidate: Any | None = None,
    *,
    protein: dict[str, Any] | None = None,
    accessibility: dict[str, Any] | None = None,
    paralogs: list[dict[str, Any]] | None = None,
    pathway_hints: list[str] | None = None,
    cell_context: list[Any] | None = None,
    fetch: bool = False,
    plan_id: str | None = None,
) -> dict[str, Any]:
    """Build an experimental plan. Returns a plain JSON-able dict.

    Parameters
    ----------
    brief, target, candidate
        ``ResearchBrief``, ``Target`` and ``Candidate`` records from ``e2b.contracts``,
        or plain dicts with the same field names. ``candidate`` may be ``None``, in
        which case the plan is returned with ``status='partial_no_candidate_sequence'``
        and the construct-liability screen is skipped rather than invented.
    protein, accessibility
        Output of ``e2b.adapters.structures.fetch_protein`` / ``assess_accessibility``.
        Pass them in when the spine already has them.
    fetch
        When true and ``protein`` is absent, retrieve UniProt live from the target's
        accession. A retrieval failure returns ``status='blocked_no_target_biology'``
        with the transport error, and is never reported as absence of an ectodomain.
    paralogs
        Paralog panel from the downstream analysis lane, as dicts with at least
        ``gene_symbol``. Absent means the specificity arm is returned as provisional.
    pathway_hints
        Pathway names or descriptions from the target's evidence, used to derive the
        mechanism class.
    """
    b, t, c = as_dict(brief), as_dict(target), as_dict(candidate)
    pid = plan_id or _plan_id(b, t, c)
    direction, direction_basis = normalise_direction(
        get(b, "moa"), get(b, "direction"), get(b, "desired_effect"), get(t, "proposed_modulation")
    )
    missing: list[str] = []
    provisional: list[str] = []

    inputs_seen = {
        "brief_id": get(b, "brief_id"),
        "target_id": get(t, "target_id"),
        "target_accession": get(t, "protein_accession"),
        "candidate_id": get(c, "candidate_id"),
        "candidate_origin": get(c, "origin"),
        "candidate_sequence_supplied": bool(get(c, "sequence")),
        "protein_record_supplied": bool(protein),
        "paralog_panel_supplied": bool(paralogs),
        "pathway_hints_supplied": len(pathway_hints or []),
        "cell_context_records": len(cell_context or []),
        "direction_basis": direction_basis,
    }

    base = dict(
        plan_id=pid,
        brief_id=get(b, "brief_id"),
        raw_indication=get(b, "raw_indication"),
        desired_effect=get(b, "desired_effect"),
        direction=direction,
        modality=get(b, "modality"),
        species=get(b, "species"),
        target_id=get(t, "target_id"),
        candidate_id=get(c, "candidate_id"),
        run_id=get(c, "run_id"),
        candidate_origin=get(c, "origin"),
        candidate_sequence_length=(len(get(c, "sequence")) if get(c, "sequence") else None),
        scope_statement=SCOPE_STATEMENT,
        what_this_does_not_establish=list(DOES_NOT_ESTABLISH),
        inputs_seen=inputs_seen,
    )

    # ---- target biology -------------------------------------------------------------
    try:
        profile, mech = build_target_profile(
            target,
            protein=protein,
            accessibility=accessibility,
            pathway_hints=pathway_hints,
            desired_effect=get(b, "desired_effect"),
            fetch=fetch,
        )
    except TargetBiologyUnavailable as exc:
        return ExperimentalPlan(
            status="blocked_no_target_biology",
            failure_detail=(
                f"Target biology could not be retrieved: {exc}. This is a SOURCE FAILURE, not evidence "
                "that the target lacks an extracellular domain, and no plan is derivable from it. "
                f"Transient: {getattr(exc, 'transient', None)}."
            ),
            missing_inputs=["UniProt topology for the target accession"],
            **base,
        ).to_dict()

    if profile.archetype == "unknown":
        return ExperimentalPlan(
            status="blocked_no_target_biology",
            target_profile=profile,
            failure_detail=(
                "No UniProt topology verdict was available for this target, so the assay archetype "
                "(secreted ligand versus cell-surface receptor) is undetermined. The two need different "
                "antigen reagents, different assay orientations and different functional designs, so a "
                "plan written without it would be a guess. Pass `protein`/`accessibility` from the "
                "structure adapter, or call with fetch=True and a valid accession."
            ),
            missing_inputs=["UniProt topology / accessibility verdict for the target"],
            **base,
        ).to_dict()

    if profile.archetype == "not_accessible_to_soluble_binder":
        return ExperimentalPlan(
            status="blocked_target_not_accessible",
            target_profile=profile,
            failure_detail=(
                "UniProt topology gives this target no extracellular face: no signal peptide, no "
                "transmembrane segment and no extracellular topological domain. A soluble binder cannot "
                "reach it in a cell, so no functional blockade experiment can honestly be designed for "
                "this modality. A biochemical binding assay against recombinant protein would still "
                "produce numbers, and those numbers would NOT constitute evidence of mechanism, because "
                "the binder never reaches the protein in its native compartment. The correct action is "
                "to change target or change modality, not to run the assay and caveat it. Basis: "
                + "; ".join(profile.archetype_basis)
            ),
            missing_inputs=[],
            **base,
        ).to_dict()

    if direction == "degrade":
        return ExperimentalPlan(
            status="blocked_direction_not_addressable",
            target_profile=profile,
            failure_detail=(
                "The brief asks for degradation. A binder on its own does not degrade its target: it "
                "would need a degradation-competent format (a bispecific to a degradation receptor, or a "
                "conjugate), and the assay that proves degradation is a measurement of target protein "
                "loss, which is a different experiment from the ones below. Confirm the intended format "
                "before a plan is written."
            ),
            missing_inputs=["the degradation-competent format intended for this modality"],
            **base,
        ).to_dict()

    # ---- routes and assays ----------------------------------------------------------
    seq = get(c, "sequence")
    if not seq:
        missing.append(
            "candidate sequence (the construct-liability screen, free-cysteine and sequon checks did "
            "not run)"
        )
    needs_disulfide = bool(seq) and (seq.upper().count("C") >= 2)
    routes = [_binder_production_route(c, profile, needs_disulfide=needs_disulfide)]
    antigen = _antigen_production_route(profile, protein)
    if antigen:
        routes.append(antigen)

    assays = [_binding_assay(profile, c), _orthogonal_assay(profile)]
    functional = _functional_assay(profile, mech, direction, b)
    if functional:
        assays.append(functional)
    else:
        missing.append(
            "mechanism class for the target: "
            + (mech.get("missing") or "no signalling-class vocabulary matched")
        )
    spec = _specificity_assay(profile, paralogs)
    assays.append(spec)
    if not paralogs:
        provisional.append(
            "Selectivity: no paralog panel was available, so the specificity arm names a method but not "
            "its members and no selectivity claim can be made from this plan as it stands."
        )

    if cell_context:
        provisional.append(
            f"{len(cell_context)} cellular-evidence record(s) were supplied and are referenced only to "
            "propose which primary population the functional assay should be confirmed in. Transcript "
            "detection is not surface protein abundance, so the cell system still has to be verified "
            "in-house."
        )

    if not seq:
        status = "partial_no_candidate_sequence"
    elif functional is None:
        status = "partial_no_functional_readout"
    else:
        status = "plan_complete"

    plan = ExperimentalPlan(
        status=status,
        target_profile=profile,
        production_routes=routes,
        assays=assays,
        decision_tree=_decision_tree(profile, mech, direction, functional is not None),
        missing_inputs=missing,
        provisional_claims=provisional,
        failure_detail=(
            None
            if status == "plan_complete"
            else "Plan returned as partial. Missing: " + "; ".join(missing)
        ),
        **base,
    )
    return plan.to_dict()
