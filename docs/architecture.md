# System design

Evidence-to-binder · Boston Computational Biology Hackathon · 22 September 2026

Written for the whole team. No part of it assumes you worked on the lane it describes.

---

## 1. What the system is

A scientist types an indication and the mechanism they want. The system decides **which
protein to target, why, and whether a binder design job is worth running** — then runs it
and reports what actually came back.

```
"atopic dermatitis"                                    a reviewed recommendation
"suppress type 2 inflammatory signalling"    ───▶      with its evidence trail,
"antagonise"                                           or an honest limitation
"de novo minibinder"
```

Two properties define it, and everything in the design follows from them.

**It is indication-agnostic.** There is no disease list, no target shortlist and no dataset
hardcoded in the production path. Changing the indication string changes everything
downstream. Verified: atopic dermatitis resolves to MONDO_0004980 (3,207 associated
targets) and non-small cell lung carcinoma to MONDO_0005233 (12,475) with **zero overlap**
in their top-20 pools and no code change.

**It is built to be honest about its own limits.** A target-selection tool that overstates
its evidence is worse than no tool, because it spends other people's months. So the
honesty rules are not review guidelines — they are validators, caps and gates in code.
Section 7 lists them.

---

## 2. The shape of the whole thing

```
                        ┌──────────────────────────────────┐
     scientist ────────▶│  UI  (e2b/ui)                    │◀──── reads records,
                        │  brief · compare · inspect ·     │      never computes
                        │  design · results                │
                        └───────────────┬──────────────────┘
                                        │ ui_bridge.publish()
                        ┌───────────────▼──────────────────┐
                        │  AGENT LOOP  (e2b/agent)         │
                        │  observe → choose → execute →    │
                        │  validate → update → repeat      │
                        │                                  │
                        │  the model chooses ACTIONS       │
                        │  Python computes every NUMBER    │
                        └───────────────┬──────────────────┘
                                        │ 13 typed tools
        ┌───────────────┬───────────────┼───────────────┬──────────────┐
        ▼               ▼               ▼               ▼              ▼
  ┌──────────┐  ┌─────────────┐  ┌────────────┐  ┌───────────┐  ┌────────────┐
  │ adapters │  │  pipeline   │  │  adapters  │  │ analysis  │  │experiments │
  │          │  │             │  │            │  │           │  │            │
  │open_     │  │ ranking     │  │cellxgene   │  │ triage    │  │ wet-lab    │
  │targets   │  │ graph       │  │protein_    │  │ paralogs  │  │ plan       │
  │structures│  │ ui_bridge   │  │design      │  │           │  │            │
  └────┬─────┘  └──────┬──────┘  └─────┬──────┘  └─────┬─────┘  └──────┬─────┘
       │               │               │                │               │
       ▼               ▼               ▼                ▼               ▼
  Open Targets    deterministic    CELLxGENE      UniProt/Ensembl   UniProt
  UniProt         rubric +         Census         paralogues        topology
  AlphaFold       evidence graph   Modal GPU
  PDBe SIFTS
                                        │
                        ┌───────────────▼──────────────────┐
                        │  CONTRACTS  (e2b/contracts.py)   │
                        │  21 typed records — the common   │
                        │  language every layer speaks     │
                        └──────────────────────────────────┘
```

**The one-sentence version of the architecture:** every layer communicates in typed
records, the agent decides *what to do next*, deterministic Python decides *what
everything is worth*, and the UI renders the result without recomputing any of it.

---

## 3. The journey, end to end

What actually happens on a run, and which module owns each step.

| # | Step | Owner | In | Out |
|---|---|---|---|---|
| 1 | Resolve the indication to an ontology ID | `adapters/open_targets` | free text | ranked `DiseaseMatch` list |
| 2 | Retrieve a **bounded** target pool | `adapters/open_targets` | disease ID | `Target` records + association scores |
| 3 | Screen which targets a binder can physically reach | `adapters/structures` | UniProt accession | topology verdict + engageable regions |
| 4 | Pull deep evidence for survivors | `adapters/open_targets` | target + disease | `Evidence` records, drugs, liabilities |
| 5 | Pull disease-context expression | `adapters/cellxgene` | disease + genes | `CellSummary` per donor |
| 6 | Score and compare, deterministically | `pipeline/ranking` | all the above | `Assessment` + `ComparisonResult` |
| 7 | Human reviews target and binding site | UI | comparison | approved `DesignRequest` |
| 8 | Run the design job | `adapters/protein_design` | design request | `DesignRun` → `Candidate` + `Metric` |
| 9 | Triage the candidates | `analysis/` | candidates + epitope | verdicts, cross-reactivity risk |
| 10 | Derive the wet-lab plan | `experiments/` | target + candidate | assays, controls, kill rules |

Steps 1–6 always run. Step 7 is a **human gate** — nothing paid happens without it. Steps
8–10 run only if that gate is passed, and the run ends honestly at step 6 if it is not.

---

## 4. Module by module

### `e2b/contracts.py` — the common language
**21 typed records**, Pydantic-validated. Everything that crosses a module boundary is one
of these. The important ones:

| Record | Carries | Produced by | Consumed by |
|---|---|---|---|
| `ResearchBrief` | the question: indication, desired effect, MoA, modality | UI / CLI | everything |
| `DiseaseMatch` | one candidate ontology resolution + match type | open_targets | brief, UI |
| `Target` | gene, protein accession, locations, tractability | open_targets | ranking, structures |
| `Evidence` | one sourced claim: subject, predicate, object, locator | open_targets, structures, cellxgene | graph, ranking |
| `CellSummary` | per-donor expression in one population | cellxgene | ranking, UI |
| `CriterionRating` | one criterion's rating + the evidence IDs behind it | ranking | Assessment |
| `DesignFeasibility` | can a binder reach this, and is there a structure | ranking | Assessment, UI |
| `Assessment` | everything known about one target, scored | ranking | comparison, UI |
| `ComparisonResult` | the ordering, ties, exclusions, pending screens | ranking | UI, agent |
| `BindingRegion` | where on the protein, in which numbering | structures | design, analysis |
| `DesignRequest` | the reviewed, approved ask | UI + human | design |
| `DesignRun` / `Candidate` / `Metric` | what the job did and produced | design | analysis, UI |
| `AgentRun` / `ToolAction` / `Budgets` | what the agent did and what it spent | agent | UI, audit |

Every record carries `schema_version` and an **origin label**:

- `real` — retrieved or computed now, this run
- `cached_real` — a real record replayed, with its original timestamp shown
- `fixture` — a hand-written stand-in for wiring a UI, **never a biological result**

### `e2b/config.py` — the rubric, as data
The scoring rubric lives here as data with a version string, not as code, so the model can
read and explain it but cannot invent it. Four criteria, weights summing to 1.0:

| Criterion | Weight | Evidence it reads |
|---|---|---|
| `clinical_causal_support` | 0.30 | genetic association, genetic literature, clinical, known drug, somatic mutation |
| `mechanistic_support` | 0.25 | affected pathway, animal model, literature |
| `direction_fit` | 0.25 | annotated drug mechanisms, partitioned by indication |
| `cellular_context` | 0.20 | single-cell expression (never bulk, as a substitute) |

Ratings are ordinal: `strong` 1.0, `moderate` 0.6, `weak` 0.25, **`unknown` 0.0**.
Liabilities subtract 0.05 each, capped at 0.20. Current version: **`rubric-2026.09.22-v2`**.

**Design feasibility is deliberately not a criterion.** It is reported separately and acts
as a gate, so a target never looks *biologically* better because it happens to be
druggable.

### `e2b/adapters/` — the outside world
Each adapter owns one source, returns typed records, and attaches a `Provenance` block.

- **`open_targets.py`** — GraphQL. Resolution, bounded discovery, deep evidence, drugs,
  tractability, liabilities. `verify_schema()` re-checks at runtime that the fields it
  depends on still exist, because the schema moves between releases.
- **`structures.py`** — UniProt, AlphaFold, PDBe SIFTS. Protein identity, **topology-based
  binder accessibility**, structure coverage, and PDB↔UniProt residue-numbering
  reconciliation.
- **`ontology.py`** — EBI OLS4 / MONDO. **Disease identity across sources.** Different
  sources name the same disease differently, and this is the module that stops that from
  losing data. See §8.
- **`cellxgene.py`** (lane 2) — Census cells, donor-aware expression, precomputed embeddings.
- **`protein_design.py`** (lane 3) — the GPU design route on Modal.

### `e2b/pipeline/` — deterministic computation
- **`ranking.py`** — turns evidence into `Assessment`s and an ordering. **No LLM is
  consulted anywhere in this file.** Identical inputs give an identical order; the sort key
  is `(-score, target_id)` so dict ordering cannot influence it.
- **`graph.py`** — the evidence graph: nodes, edges, and the `Evidence` behind every edge.
  Refuses to add a scientific relation with no evidence. Each edge carries a `claim_class`
  (`asserted_literature`, `observed_metadata`, `computed_measurement`, `model_prediction`)
  so a model prediction never reads as a measurement.
- **`ui_bridge.py`** — pushes spine records into the UI store, including failures.

### `e2b/agent/` — the investigation loop
- **`tools.py`** — the 13 typed tools, budget enforcement, retries, precondition checks.
- **`state.py`** — `RunStore`: the persisted run. Every run writes an immutable directory.
- **`loop.py`** — the loop itself, plus the injectable LLM clients.

### `e2b/ui/` — the review interface
Server-rendered FastAPI + Jinja2. No build step, no npm, no client framework; readable with
JavaScript off. Five sections on one page: **Brief · Compare · Inspect · Design · Results**.
It is a *reader* — it renders records and shows their origin, and computes nothing.
17 honesty assertions run against the rendered HTML (`python -m e2b.ui.selfcheck`).

### `e2b/analysis/` and `e2b/experiments/` — downstream
- **`analysis/`** — deterministic candidate triage: sequence developability, interpretation
  of interface metrics, and paralogue cross-reactivity over the *engaged epitope* rather
  than whole-sequence identity.
- **`experiments/`** — derives a wet-lab plan from the target's own biology: production
  route and QC, a binding assay with orientation reasoning, an orthogonal confirmation, a
  mechanism-specific functional assay, a specificity arm, and pre-registered go/no-go rules.
  Every assay states what it **does not** report.

---

## 5. The agent loop

```
   ┌──────────────────────────────────────────────────────┐
   │  observe        read PERSISTED STATE, not chat history│
   │  choose         pick the single most useful next tool │
   │  execute        run one typed tool                    │
   │  validate       check the result, record failures     │
   │  update         write records, extend the graph       │
   └───────────────────────┬──────────────────────────────┘
                           │ continue · revise · ask · stop
```

**The division of labour is the whole point.**

| The model does | The model never does |
|---|---|
| chooses which tool to call next | computes a score |
| explains a result using evidence IDs | invents a rubric weight |
| identifies which missing evidence matters | decides a design job may launch |
| asks a clarifying question when it changes the answer | declares an outcome unchecked |

**Bounds.** 20 candidate targets, 3 investigated in depth, 30 tool decisions, 2 retries on
transient failures, **1 paid design batch**. Budget exhaustion returns the best supported
partial result.

**Stages.** Nine, and the UI names them rather than naming tools, so renaming a tool does
not change what a scientist sees: `resolve_indication → gather_targets → assess_targets →
cellular_context → rank_targets → structure_and_site → design_request → design_run →
report`.

**Outcomes.** Every run ends in exactly one, and `_reconcile()` checks the model's claim
against what actually happened:

| Outcome | Meaning |
|---|---|
| `completed_candidate_report` | candidates were generated and evaluated |
| `evidence_report_only` | a reviewed target proposal, no design run completed |
| `needs_clarification` | the indication was ambiguous enough to change the answer |
| `insufficient_evidence` | retrieved evidence could not support any ranking |
| `no_binder_feasible_target` | every candidate was unreachable by this modality |
| `source_unavailable` | a source failed — **not** evidence of absence |
| `budget_exhausted` | bounds reached; best supported partial result |
| `design_failed` | a job ran and produced nothing acceptable |
| `embedding_unavailable` | no suitable Census release/embedding combination |

**The LLM is injected, not imported**, so the same loop runs three ways:
`HostLLMClient` (inside Claude Science, no API key in the app), `AnthropicClient` (from a
Modal container, key from a Secret), or any object with `.complete(messages, tools, system)`
for tests.

---

## 6. The 13 tools

| Tool | Required in | Returns | Paid |
|---|---|---|---|
| `resolve_indication` | `raw_text` | ranked matches, `needs_clarification` flag | |
| `discover_targets` | `disease_id` | bounded pool + association breakdown | |
| `assess_structure` | `accession` | accessibility verdict, engageable regions, structures | |
| `get_target_evidence` | `target_id`, `disease_id` | evidence rows, drugs, liabilities | |
| `compare_targets` | — | ordering, ties, exclusions, pending screens | |
| `summarize_cell_context` | `disease_label`, `gene_symbols` | donor-aware expression | |
| `discover_census_embeddings` | — | available precomputed embeddings | |
| `explore_cell_neighborhoods` | `embedding_slice_id`, `query_population` | neighbours + search universe | |
| `submit_design` | `design_id` | run ID | **yes** |
| `get_run` | `run_id` | status, stage, artifacts | |
| `evaluate_candidates` | `run_id` | per-candidate metrics | |
| `triage_candidates` | `run_id` | developability + specificity verdicts | |
| `plan_experiments` | `target_id` | assays, controls, decision criteria | |

Tool results are **compact summaries plus artifact references**. Expression matrices and
embedding vectors stay on disk; they never enter the model's context.

A tool whose lane has not landed returns `{"status": "capability_unavailable"}` — a result
the agent branches on, not a crash, and emphatically not a plausible fake value.

---

## 7. How honesty is enforced

Not by review discipline — by code. Each of these is a specific mechanism you can go read.

| Rule | Where it lives |
|---|---|
| Evidence that asserts nothing is rejected | `contracts.Evidence` validator |
| A `real` record with no locator is rejected | `contracts.Evidence` validator |
| A favourable rating with no evidence IDs is rejected | `contracts.CriterionRating` validator |
| `unknown` must say what was looked for and not found | `contracts.CriterionRating` validator |
| A metric named `Kd`, `IC50` or `affinity` raises | `contracts.Metric` validator |
| Unknown numbers are `None`, never `0.0` | convention, checked in tests |
| Missing evidence scores 0.0 and cannot outrank evidence | `STRENGTH_POINTS["unknown"] = 0.0` |
| Text-mined literature alone caps at `moderate` | `LITERATURE_ONLY_RATING_CAP` |
| Healthy-tissue-only expression caps at `moderate` | `REFERENCE_ONLY_RATING_CAP` |
| A single donor caps at `weak` | `ranking._rate_cellular` |
| Drug precedent in another indication does not transfer | `ranking._rate_direction` partition |
| Unscreened ≠ not-excluded | `ComparisonResult.pending_screen` |
| A design job needs an approved, complete request | `DesignRequest.assert_launchable()` |
| Identical requests do not start a second paid job | `compute_input_hash()` dedupe |
| Fixtures are dropped before they reach triage | `RunStore.record_candidates()` |
| The model cannot claim candidates that do not exist | `AgentLoop._reconcile()` |
| A failed source is rendered as a failure, not a blank | `ui_bridge` → `SourceFailure` |
| Ties are shown as ties | `ComparisonResult.ties` |

---

## 8. Three design decisions worth understanding

**Disease identity is resolved on ontology terms, never on labels.** Open Targets resolves
"atopic dermatitis" to **MONDO_0004980, whose preferred label is "atopic eczema"** —
"atopic dermatitis" is only a synonym. A source search on the literal input string finds
nothing and concludes the disease is absent. That happened here, on real data: a string
search of CELLxGENE reported no atopic dermatitis anywhere, while the ontology-grounded
mapping finds `MONDO:0004980` immediately — one dataset, 280,518 cells of skin, carrying
both the disease and `normal` in the same study.

`adapters/ontology.py` therefore maps a disease onto any source's vocabulary by identifier
and records **how** each match was reached:

| Match | Meaning | Same disease? |
|---|---|---|
| `exact` | same term | yes |
| `cross_reference` | another ontology's ID for it — `EFO:0000274` ↔ `MONDO:0004980` | yes |
| `synonym_label` | matched a recorded synonym; weaker than an ID match | yes, flagged |
| `descendant` | a **sub-type** — narrows the question | no |
| `ancestor` | a **parent** — broadens it, includes other patients | no |

Descendants and ancestors are returned, never silently folded in. For this query the
Census vocabulary of 331 disease terms yields exactly one same-disease term and two
broader parents (`inflammatory disease`, `immune system disorder`) which are *not* the
queried population. This also explains an earlier puzzle cleanly: `EFO:0000274` returning
null from Open Targets was an identifier migration between releases, not a different
disease — the ontology records the two as the same thing.



**Accessibility is decided on topology, not location strings.** A soluble binder can only
reach a secreted protein or an extracellular domain. The obvious implementation — match
keywords against subcellular-location annotations — is wrong: PDE4A, PDE4B and PDE4D are
all annotated "plasma membrane" because they dock on the *cytoplasmic* face, so a keyword
screen calls them accessible and would aim a binder at a cytosolic enzyme. The screen
therefore requires a UniProt **signal peptide with no transmembrane segment** (secreted) or
an explicit **`Extracellular` topological domain** (ectodomain). A region smaller than
`MIN_EPITOPE_RESIDUES = 40` is excluded as too small to present a useful epitope — this
caught HRH1, whose extracellular N-terminus is 29 residues.

**Residue numbering is reconciled explicitly.** Design tools address residues in PDB author
numbering; sequence analysis addresses them in UniProt canonical numbering. These differ —
for 1IAR the offset is **+25**, so hotspot "59" is really residue 84. `structures.py`
resolves the chain *from the accession* (1IAR chain A is IL-4, **not** IL-4R; 4I77's IL-13
chain is Z) and records which numbering basis the offset applies to, because 3L5X and 4I77
are the same protein with different offsets.

---

## 9. What a run leaves behind

```
runs/<run_id>/
  brief.json                 the question, and how the indication resolved
  pool.json                  what was retrieved, and out of how many
  structure_<acc>.json       topology, accessibility, structures
  evidence_<target>.json     evidence rows, drugs, liabilities, duplicates collapsed
  cell_summaries.json        donor-aware expression
  assessments.json           per-target ratings with evidence IDs
  comparison.json            ordering, ties, exclusions, pending screens
  design_request_*.json      what was approved, by whom, with what hash
  design_run_*.json          model, version, parameters, seed, git commit
  candidates_*.json          real candidates only
  triage_*.json              developability and specificity verdicts
  experiment_plan_*.json     the wet-lab plan
  graph.json                 the whole evidence graph
  agent_run.json             every tool call, argument, reason and error
  manifest.json              the summary
```

**These files are the record; the chat transcript is a by-product.** Written
write-then-rename so a reader never sees a half-written manifest.

---

## 10. Lane ownership

| # | Lane | Owns |
|---|---|---|
| 1 | Integration / spine | `contracts.py`, `config.py`, `pipeline/`, `agent/`, `adapters/open_targets.py`, `adapters/structures.py` |
| 2 | Cellular evidence | `adapters/cellxgene.py`, `modal_app/census_worker.py` |
| 3 | Binder design | `adapters/protein_design.py`, `modal_app/design_worker.py` |
| 4 | Evaluation / QA | `tests/`, `e2b/eval/`, `docs/verification.md` |
| 5 | Experience / UI | `e2b/ui/`, `modal_app/web.py`, `docs/demo_script.md` |
| 6 | Downstream analysis | `e2b/analysis/`, `docs/analysis_methods.md` |
| 7 | Experimental design | `e2b/experiments/`, `docs/experimental_plan.md` |

Each lane creates and edits only its own paths, reads anything, and never edits the shared
files. Lanes hand off through typed records, never by calling into each other's modules —
which is why a lane that has not landed cannot block one that has.

---

## 11. Extending it

**Add a data source** — write an adapter that returns `Evidence` with a `Provenance` block
and real `source_locator`s. Do not invent a new record type for a concept `contracts.py`
already covers.

**Add a tool** — add a `ToolSpec` in `agent/tools.py` with a JSON schema and a handler.
The handler is ordinary Python; it must not consult an LLM. Map it to a stage in
`TOOL_TO_STAGE`.

**Change the rubric** — edit `CRITERIA` in `config.py` and **bump `RUBRIC_VERSION`**.
Rankings recorded under the old version are not comparable, which is why the version
travels on every `Assessment` and `ComparisonResult`.

**Add a honesty rule** — prefer a validator in `contracts.py` or a cap in `config.py` over
a comment. If it can be enforced structurally, enforce it structurally.

---

## 12. Status

Verified working: Open Targets (release 26.06), UniProt, AlphaFold, PDBe SIFTS,
deterministic ranking, evidence graph, the agent loop end to end, the UI on real records,
downstream triage, experimental planning.

Not yet verified at the time of writing: Census retrieval, design execution through the
spine, Modal deployment of the integrated app. **No design candidate exists yet and none is
claimed.**

See `docs/capabilities.md` for the evidence behind each of those, including the six bugs
that verification caught, and `docs/decision_log.md` for who decided what and why.
