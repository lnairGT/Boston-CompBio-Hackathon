# Lane 7 — experimental design: capability report

*Boston Computational Biology Hackathon, 22 September 2026. Module `experiments/1.0.0`.
Owned paths: `e2b/experiments/`, `docs/experimental_plan.md`, this file.*

## What is proven

`e2b.experiments.plan_experiments()` takes a `ResearchBrief`, a `Target` and a
`Candidate` (pydantic records from `e2b/contracts.py`, or plain dicts with the same
field names) and returns a JSON-able experimental plan: production route and release
QC, a binding-validation concept with orientation reasoning and controls, an
orthogonal different-principle confirmation, an MoA-specific functional concept, a
specificity arm, and a pre-registered go/no-go tree with explicit kill rules.

The plan is **derived from the target's biology**, not templated. Two facts decide
every branch, and both are read at runtime from the target's own annotation:

| Decision | Source | Consequence in the plan |
|---|---|---|
| **Archetype** — secreted ligand / cell-surface receptor / unreachable | UniProt topology, via Lane 1's `adapters/structures.py` (signal peptide, transmembrane segments, `Extracellular` topological domains) | antigen construct and host, assay orientation, what the functional assay stimulates and which cell is the responder |
| **Mechanism class** — which signalling class the target drives | UniProt protein name and keywords, plus any pathway evidence the spine passes in, plus the brief's `desired_effect` | the mechanism-proximal readout, why that readout reports the mechanism, the stimulus, and the downstream confirmation |

### Self-test results (live, 22 Sep 2026, `python -m e2b.experiments.selftest`)

Three real UniProt records, chosen for topology class only, exercised end to end.
All checks passed (`ok: true`).

| Branch | Archetype derived | Status | Antigen construct | Binding orientation |
|---|---|---|---|---|
| secreted protein (P35225, 146 aa, 4 glyco sites, 2 disulfides, 13 PDB entries) | `secreted_ligand` | `plan_complete` | mature secreted protein, residues 25–146 | immobilise the **binder**, flow the target |
| cell-surface receptor (P24394, 825 aa, 6 glyco sites, 2 disulfides, 10 PDB entries) | `cell_surface_receptor` | `plan_complete` | ectodomain, residues 26–232 | capture the **target ectodomain**, flow the monomeric binder |
| cytosolic enzyme (Q07343, 736 aa, 42 PDB entries) | `not_accessible_to_soluble_binder` | `blocked_target_not_accessible` | none — refused | none — refused |

Both accessible branches produced a complete plan with 4 assay concepts and 5
pre-registered gates; the mechanism class was derived as cytokine/JAK-STAT for both
from their own annotation. The two orientations and the two antigen constructs differ,
which is the check that the module is branching rather than templating. The antigen
*host* is mammalian for both, correctly: a glycosylated secreted protein and a
receptor ectodomain both want a mammalian secretory host, and reporting a spurious
difference there would be worse than reporting none.

Mechanism-class generality was additionally checked on shaped records: a secreted
protease routed to `protease` (substrate-cleavage readout), a secreted growth factor
to `rtk_growth_factor` (receptor autophosphorylation, then pERK/pAKT), and an
uncharacterised secreted protein to `unknown`, which returns
`partial_no_functional_readout` and states what evidence is missing rather than
inventing a plausible assay.

### Failure states, proven rather than asserted

| Input condition | Returned status | Behaviour |
|---|---|---|
| target with no extracellular face | `blocked_target_not_accessible` | refuses to design an in-cell blockade experiment; says a biochemical binding number from such a target would not be evidence of mechanism |
| unknown / malformed accession (UniProt answers HTTP 200 with an empty body) | `blocked_no_target_biology` | labelled **source failure**; an empty record is explicitly not read as "no ectodomain" |
| no candidate supplied | `partial_no_candidate_sequence` | plan returned; `candidate_sequence_length` is `null`, not `0`; the construct-liability screen is recorded as not run |
| mechanism class not derivable | `partial_no_functional_readout` | binding routes returned, functional gate `G3` explicitly marked not specifiable |
| brief asks for degradation | `blocked_direction_not_addressable` | a binder alone does not degrade; asks for the intended format |
| no paralog panel | `plan_complete` with a provisional flag | the specificity arm names the method but not its members, and says no selectivity claim can be made |

### Honesty invariants, checked programmatically

`selftest.check_invariants(plan)` runs on every generated plan and returned zero
violations for all cases. It fails a plan that:

- states a specific concentration or incubation time (regex over every string in the
  plan for `N nM|pM|uM|mM` and `N min|h`), or writes `Kd =`, `IC50 =`, `EC50 =`;
- carries an assay parameter not tagged `starting_point_requires_optimisation`;
- carries a decision criterion without a kill rule, or with a threshold not tagged
  `pre_registered_default_requires_sign_off`;
- refers to the design as an antibody;
- omits the scope statement or the "does not establish" block.

On the unknown-numerics rule the check is narrower than "no zeros anywhere", and the
distinction matters: a protein really can have zero annotated glycosylation sites, and
turning that into `null` would lose information. What is checked is whether a zero is
an *undetermined* value in disguise. `length_aa` is `None` exactly when no UniProt
record was read, so the check fires when `n_glycosylation_sites`, `n_disulfides` or
`experimental_structures` is `0` while `length_aa` is `None`; when `length_aa` is `0`
(the empty-record signature); and when `candidate_sequence_length` is `0`. Legitimate
measured zeros alongside a real record pass, by design. Verified with negative
controls: each of those four conditions fires, and a real plan carrying three genuine
zero counts does not.

Not covered by this check: a count that is spuriously zero *while a record was
successfully read* would pass. That case would require the upstream annotation itself
to be wrong, which this module cannot detect.

Structurally, the `StartingPoint` record has **no value field** — only a range and the
reasoning for it — so this module cannot emit a validated-looking number even by
accident.

## Entry points

```python
from e2b.experiments import plan_experiments, render_markdown

plan = plan_experiments(
    brief,                  # ResearchBrief record or dict
    target,                 # Target record or dict (protein_accession is what matters)
    candidate,              # Candidate record or dict; None is allowed -> partial plan
    protein=protein,        # dict from adapters.structures.fetch_protein (preferred)
    accessibility=acc,      # dict from adapters.structures.assess_accessibility
    paralogs=[{...}],       # Lane 6 paralog panel; None -> specificity arm provisional
    pathway_hints=[...],    # pathway names from the target's evidence
    cell_context=[...],     # CellSummary records; referenced for primary-cell confirmation
    fetch=False,            # True -> retrieve UniProt live from target.protein_accession
)                           # -> plain JSON-able dict, never raises on missing inputs

md = render_markdown(plan, title="# Next experimental test")   # -> markdown string
```

Support functions, if the spine wants them directly:
`build_target_profile(target, protein=..., ...) -> (TargetAssayProfile, mechanism_dict)`,
`classify_mechanism(...) -> dict`, `normalise_direction(*texts) -> (direction, basis)`,
`screen_binder_sequence(seq, expression_host=...) -> list[ConstructRisk]`,
`sequence_properties(seq) -> dict`.

Plan dict shape: `status`, `failure_detail`, `missing_inputs`, `provisional_claims`,
`target_profile`, `production_routes[]`, `assays[]` (roles `binding_primary`,
`binding_orthogonal`, `functional_moa`, `specificity`), `decision_tree[]`,
`scope_statement`, `what_this_does_not_establish[]`, `inputs_seen`.

The module reads Lane 1 adapters only (`adapters/structures.py`) and only when
`fetch=True` or when it needs `assess_accessibility` on a supplied record. It edits no
file outside `e2b/experiments/` and `docs/`.

## Not proven / out of scope

- **No real candidate has been planned yet.** The self-test candidate is a synthetic
  string labelled `origin='fixture'`, used only to exercise the free-cysteine and
  glycosylation-sequon screens. It carries no metrics because none exist. When Lane 3
  delivers, the same call regenerates the plans with the real sequence.
- **No paralog panel exists yet**, so the specificity arm is method-only and flagged
  provisional in every generated plan.
- **No pathway evidence from Open Targets has been wired in.** The mechanism class was
  derived from UniProt annotation plus brief language; `pathway_hints` is accepted and
  used but has not been exercised with real Open Targets pathway records.
- **Nothing here was executed on Modal**, and nothing here needs a GPU.
- The plan names assay concepts. It is not a protocol, no parameter in it has been
  measured, and no result obtained from it would establish clinical efficacy or safety.

## Notes for the lead

1. `adapters/structures.fetch_protein()` returns an **empty record** (length 0, no
   features) for an unknown or malformed accession, because UniProt answers HTTP 200
   with an empty body instead of 404. Any accessibility screen run on that record will
   report "not accessible", which is a fabricated negative. Lane 7 guards against it
   (`target_profile.is_empty_protein`), but the guard belongs in the adapter so every
   lane gets it.
2. If the UI wants a short panel rather than the full document, render only
   `assays[role='functional_moa']` and the `decision_tree` — the kill rules are the
   part a reviewer looks for.
