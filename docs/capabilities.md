# Capabilities, execution routes and blockers

Boston Computational Biology Hackathon · 22 September 2026 · recorded by the integration lane.

This file records what was **actually verified in this environment**, with the evidence that
verified it. Anything unverified says so. Platform documentation establishes what a platform
can do in general; it does not establish access in our workspace, so nothing below is
asserted on the strength of documentation alone.

## Verified working

| Capability | Route | Evidence |
|---|---|---|
| Open Targets GraphQL | `POST https://api.platform.opentargets.org/api/v4/graphql` | Live. API 26.6.3, data version **26.06**. Every field the adapter uses was introspected against the live schema; `verify_schema()` re-checks them at runtime and reports any that vanish. |
| Disease resolution | `search(entityNames:["disease"])` | `"atopic dermatitis"` → `MONDO_0004980` ("atopic eczema"), 5.5× score margin. `"non-small cell lung carcinoma"` → `MONDO_0005233`. |
| Bounded target discovery | `disease.associatedTargets(page:)` | AD: 3,207 associated, pool of 20 retrieved. NSCLC: 12,475 associated, pool of 20. |
| Target evidence + drugs | `disease.evidences(ensemblIds:)`, `target.drugAndClinicalCandidates` | 844 evidence items for IL13/AD; retrieved 50. Indication-matched agents recovered with mechanism and stage. |
| UniProt identity + topology | `https://rest.uniprot.org/uniprotkb/{acc}.json` | Signal peptide, transmembrane and topological-domain features retrieved for all 20 pool members. |
| AlphaFold coverage | `https://alphafold.ebi.ac.uk/api/prediction/{acc}` | Reachable; queried as fallback when no experimental structure exists. |
| Deterministic ranking | `e2b.pipeline.ranking` | Ran over the real 20-target pool: 7 ranked, 13 excluded with reasons, three-way tie surfaced rather than broken. |
| Evidence graph | `e2b.pipeline.graph` | 563 nodes, 517 edges, 345 source-backed evidence records, 6 preserved as ontology-propagated. |
| Agent tool-use loop | `e2b.agent.loop` via `host.llm` | Ran end to end on the real brief: 10 tool decisions, 0 errors, correct `evidence_report_only` outcome. |

## Network access

Granted on request in this sandbox (both were initially blocked):

- `census.cziscience.com` — Census release manifest.
- `cellxgene-census-public-us-west-2.s3.us-west-2.amazonaws.com` — the Census SOMA store.

Already allowlisted: `api.platform.opentargets.org`, `rest.uniprot.org`, `alphafold.ebi.ac.uk`,
`api.cellxgene.cziscience.com`, PyPI.

## Capability gaps found, and what was done about them

**No de novo backbone generator in the tool catalog.** Searched it directly: there is no
RFdiffusion and no BindCraft skill. What exists is ProteinMPNN, SolubleMPNN, LigandMPNN
(inverse folding — they need a backbone as input), and Boltz-2, Chai-1, AlphaFold2, ESMFold2
(structure prediction / validation). The team's decision was to install the generator on Modal
and run it there, which is the route lane 3 is building. Until that lands, no design claim of
any kind is made: `adapters/protein_design.py` being absent makes the design tools return
`capability_unavailable`, and the agent's outcome reconciler downgrades any claimed
`completed_candidate_report` to `evidence_report_only`.

**No GitHub write access.** The repo `lnairGT/Boston-CompBio-Hackathon` is public and was
clonable, but it is **empty — zero commits**, and no GitHub credential is configured, so
nothing can be pushed. Per spec section 5, an empty repo means defaulting to a small Python
app with a server-rendered interface, which is what was built. **Blocker: a GitHub credential
with repo write scope is needed before any of this reaches the remote.**

## Bugs found by verifying rather than assuming

These are recorded because each would have produced confidently wrong science:

1. **`EFO_0000274` returns `null`** for atopic dermatitis in release 26.06. The live ID is
   `MONDO_0004980`. Hardcoding an ontology ID from memory retrieves nothing, silently.
2. **`ClinicalDiseaseListItem` has no `id`/`name`** — it exposes `diseaseFromSource` and a
   nested `disease`. The wrong shape returned HTTP 400 and aborted *all* evidence retrieval,
   which presented as every target scoring `weak` in a six-way tie.
3. **Association datatype vocabulary is release-specific.** 26.06 emits `clinical` and
   `genetic_literature`; it does **not** emit `known_drug`. A rubric keyed on the documented
   older names rated every criterion `unknown`.
4. **Clinical stages are `APPROVAL` / `PHASE_3`**, not `"Approved"` / `"Phase III"`. Dupilumab's
   approval initially contributed zero to direction fit.
5. **Keyword-based accessibility screening is wrong.** PDE4A/B/D are annotated "plasma
   membrane" / "apical cell membrane" because they dock on the cytoplasmic face. A keyword
   screen calls them binder-accessible. The screen was rewritten to decide on UniProt topology
   (signal peptide, transmembrane, `Extracellular` topological domain), which correctly rejects
   all four PDE4s, FLG, NR3C1 and the JAK family.
6. **Unscreened targets were entering the shortlist.** With `accessible_to_modality = None`
   they were neither excluded nor penalised, so CYP24A1 (mitochondrial) reached shortlist
   position 3 on the absence of adverse evidence. Now gated into `ComparisonResult.pending_screen`,
   which is explicitly *not* eligible for the shortlist and explicitly *not* a rejection.

## Required secret names (no values)

| Name | Used by | Purpose |
|---|---|---|
| `Modal` | lanes 2, 3, deployment | Modal workspace authentication. |
| `ANTHROPIC_API_KEY` | `e2b.agent.loop.AnthropicClient` | Only needed when the loop runs **outside** Claude Science. Inside it, `HostLLMClient` is used and no key exists in the app. Must be supplied as a Modal Secret, never browser-exposed. |

## Not verified

- Census cells / expression / precomputed embeddings — lane 2 in flight at time of writing.
- Any protein-design execution — lane 3 in flight. **No candidate sequence, structure or
  metric exists yet, and none is claimed.**
- Modal web serving of the integrated app.
- Cross-release Census behaviour; only the pinned release lane 2 selects is in scope.
