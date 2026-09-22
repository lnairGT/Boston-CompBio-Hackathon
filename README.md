# indication2binder

A staged pipeline that goes from a **disease indication** to a **BindCraft2 binder-design
specification**:

1. rank candidate protein targets for the indication from public evidence,
2. keep only targets whose disease-relevant interface is **experimentally determined**,
3. extract the interface residues and select hotspots,
4. emit a BindCraft2 target specification.

The pipeline never submits a design job. Stage 5 writes specification files and stops.

## Status

| Stage | Module | Does | Status |
|------:|--------|------|--------|
| 0 | `stage0_resolve` | indication string -> disease ontology id + subtree | **implemented** |
| 1 | `stage1_evidence` | disease id -> per-target evidence payloads | **implemented** |
| 2 | `stage2_score` | evidence -> ranked target table | **implemented** |
| 3 | `stage3_complexes` | ranked targets -> qualifying experimental complexes | planned |
| 4 | `stage4_interface` | complexes -> interface residues + buried SASA | planned |
| 5 | `stage5_specs` | epitopes -> BindCraft2 target spec JSON | planned |

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

Python 3.11+. Stages 0-2 need only network access. Stages 3-4 additionally use
`biopython` and `gemmi` for structure parsing, and run on a laptop.

## Quickstart

```bash
# check that the upstream GraphQL schema still matches what this package queries
ind2b schema-check

# stage 0: resolve a free-text indication to an ontology id
ind2b stage0 --indication "non-small cell lung carcinoma"

# or pin the term explicitly (the reproducible path once a term is agreed)
ind2b stage0 --efo-id MONDO_0005233

# stage 1: gather the evidence layers for the candidate targets
ind2b stage1 --efo-id MONDO_0005233 --n-targets 100

# stage 2: score and rank (re-run with --weights to re-rank without re-fetching)
ind2b stage2 --efo-id MONDO_0005233 --top 20
```

Output:

```
resolved: MONDO_0005233  non-small cell lung carcinoma
  rule           : exact_name_match
  ambiguous      : False
  therapeutic    : respiratory or thoracic disease, cancer or benign tumor
  children       : 4
  descendants    : 27
  written        : runs/MONDO_0005233/stage0_disease.json
```

## Run layout

Each stage reads the previous stage's file from the run directory and writes its own,
so any stage can be re-run in isolation.

```
runs/<disease_id>/
  stage0_disease.json        resolved term, rejected candidates, selection rule
  stage1_evidence.json       raw evidence payloads per target
  stage2_ranked_targets.csv  ranked table, one column per evidence source
  stage3_complexes.json      qualifying PDB entries per target
  stage4_epitopes.csv        interface residues per epitope
  stage5_specs/              BindCraft2 target specifications
.ind2b_cache/                disk cache of every HTTP response
```

## Design notes

**Standalone clients.** Every data source is reached over its public HTTP API
(`src/ind2b/sources/`). Nothing depends on a hosted agent runtime, so the pipeline runs
anywhere Python and network access are available.

**Everything is cached.** `src/ind2b/http.py` disk-caches every response keyed by
method, URL, params and body, throttles per host, and retries with backoff on 429/5xx.
Re-running a stage costs no network traffic, and a populated `.ind2b_cache/` replays a
whole run offline.

**Schema drift is checked, not assumed.** Open Targets renames GraphQL fields between
releases. `ind2b schema-check` introspects the live schema and fails loudly when a field
this package queries has disappeared. Run it first when something breaks.

**Association scope.** Open Targets association scores for a parent disease already
aggregate evidence from descendant terms (the "indirect" association). Pulling each
descendant separately would count the same evidence repeatedly, so stage 0 records the
subtree for provenance but the default scope stays at the parent term. The scope used is
written to `stage0_disease.json` as `association_scope`.

**Term selection is auditable.** Stage 0 records which ontology hits it rejected and the
rule that picked the winner. When the choice rests on search rank alone it sets
`ambiguous: true` and tells you to pin the term with `--efo-id`.

**Surface accessibility is a filter, not a tiebreak.** A de novo mini-binder is a protein:
absent a delivery strategy it can only reach secreted proteins and receptor ectodomains.
Subcellular localization therefore gates the target list rather than nudging its ranking.

**The experimental-interface gate.** Stage 3 admits a target only when its
disease-relevant interaction has an experimentally determined complex structure. This
cuts the candidate list substantially and biases toward well-studied biology. That is a
deliberate trade: epitope choice largely decides whether a binder campaign succeeds, and
a measured interface is the only input to that choice which is not itself a prediction.

## Data sources

| Source | Used for |
|--------|----------|
| Open Targets Platform | disease ontology, target associations (including the `literature` and `genetic_literature` evidence scores), tractability, safety, expression specificity, interaction partners |
| UniProt | topology-based accessibility, extracellular spans, sequence |
| RCSB PDB | experimental complex search, structure files |
| ChEMBL | mechanism and action-type precedent |
| ClinicalTrials.gov | clinical precedent, including stopped trials |

Literature support is taken from the Open Targets evidence breakdown rather than a
separate bibliographic database, so no additional API key is needed.

Please cite the underlying databases in any work that uses this pipeline.

## What this pipeline does not establish

- **Not causality.** Stage 2 aggregates existing evidence. A high rank means a target is
  well supported in the literature and databases, not that it is causal in the disease.
- **Not affinity.** Stage 4 identifies where a binder should bind. It does not predict
  whether any designed binder will bind, or how tightly.
- **Not a validated design.** Stage 5 emits a specification. Designs produced from it are
  in-silico candidates requiring experimental validation.
- Known-drug evidence inflates association scores, so highly drugged targets rank partly
  on their own pharmacology. Per-source columns in the stage 2 output let you see this.

## BindCraft2

Stage 5 targets [BindCraft2](https://github.com/PacesaLab/BindCraft2). Two things to know
before running it:

- It requires Linux, Python 3.12+, and an NVIDIA GPU. There is no CPU install path, so it
  cannot run on a macOS laptop - use a GPU host.
- It is released under the *BindCraft2 Source-Available License (Hosting-Restricted)*, not
  MIT. Internal use (including commercial), modification and redistribution are permitted,
  but offering it to third parties as a hosted or managed service requires a separate
  commercial licence from the Pacesa Lab. That restriction applies to BindCraft2 itself,
  not to this pipeline.

## Licence

This pipeline is MIT-licensed (see `LICENSE`). The data sources and BindCraft2 carry their
own terms.
