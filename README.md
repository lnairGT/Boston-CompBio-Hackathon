# Evidence-to-binder

An **indication-agnostic agentic research workflow** for target selection and binder design.
You type an indication and the mechanism you want; the agent resolves the disease, retrieves a
bounded candidate pool from Open Targets, screens which targets a binder could physically
reach, pulls disease-context expression from the CELLxGENE Census, compares candidates under a
versioned rubric, and — once a human has reviewed the target and binding site — submits a real
protein-design job and reports actual model output.

There is no disease allowlist, no target shortlist and no dataset hardcoded anywhere in the
production path. Swapping `"atopic dermatitis"` for `"non-small cell lung carcinoma"` changes
nothing but the input, and returns a completely different real target pool.

## What it will and will not tell you

It is built to be **honest about its own limits**, because a target-selection tool that
overstates its evidence is worse than no tool:

- "Highest confidence" means *best supported among the targets actually evaluated, under the
  stated rubric*. It is not a probability of therapeutic success.
- Biological evidence and model confidence are reported separately, and never added together.
- A computational interface score is never reported as an affinity. `contracts.Metric` raises
  on a metric named `Kd`, `IC50` or `affinity`, so it cannot happen by accident.
- Missing evidence is `unknown`, not zero, and cannot improve a target's ranking.
- When the evidence cannot separate two targets, you get a tie and an explanation of which
  missing evidence would break it — not an invented winner.
- A target that is well-associated with the disease but physically unreachable by the chosen
  modality is excluded *with its biology still visible*, and the reason recorded.

## Quick start

```bash
git clone https://github.com/lnairGT/Boston-CompBio-Hackathon.git
cd Boston-CompBio-Hackathon
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Run the web interface:

```bash
uvicorn e2b.modal_app.web:app --reload --port 8000   # then open http://localhost:8000
```

Run one investigation headless:

```bash
export ANTHROPIC_API_KEY=...        # only needed outside Claude Science
python -m e2b.cli "atopic dermatitis" \
    --desired-effect "suppress type 2 inflammatory signalling" \
    --moa antagonise \
    --modality "de novo minibinder"
```

Every run writes an immutable directory under `runs/<run_id>/` containing `brief.json`,
`pool.json`, per-target `evidence_*.json`, `comparison.json`, `assessments.json`,
`agent_run.json`, `graph.json` and `manifest.json`. Those files *are* the record — the chat
transcript is a by-product.

## How the agent works

```
observe persisted state
  -> choose the next evidence-gathering or design action
  -> execute one allowed typed tool
  -> validate the result
  -> update the evidence graph and the remaining uncertainty
  -> continue, revise, ask a targeted question, or stop
```

The division of labour is deliberate and is the reason the numbers can be trusted:
**the model chooses actions and explains results; deterministic Python computes every score,
enforces every budget, and validates every scientific precondition.** The ranking rubric lives
in `config.py` as data with a version string on it. The model cannot invent a weight, cannot
score a target, and cannot authorise a paid design job.

Typed tools available to it: `resolve_indication`, `discover_targets`, `assess_structure`,
`get_target_evidence`, `compare_targets`, `summarize_cell_context`,
`discover_census_embeddings`, `explore_cell_neighborhoods`, `submit_design`, `get_run`,
`evaluate_candidates`.

Runs are bounded: 20 candidate targets, 3 investigated in depth, 30 tool decisions, 2 retries
on transient read failures, **one paid design batch**. Budget exhaustion returns the best
supported partial result rather than stopping silently.

Every run ends in one of: `completed_candidate_report`, `evidence_report_only`,
`needs_clarification`, `insufficient_evidence`, `no_binder_feasible_target`,
`source_unavailable`, `budget_exhausted`, `design_failed`, `embedding_unavailable`. The outcome
is **reconciled against what actually happened** — the agent cannot claim a candidate report
when no design run produced artifacts.

## Layout

```
e2b/contracts.py        typed records, provenance, origin labels, the design gate
e2b/config.py           budgets + the versioned rubric (data, not code)
e2b/adapters/           open_targets, structures (UniProt/AlphaFold), cellxgene, protein_design
e2b/pipeline/           ranking (deterministic), graph (source-backed evidence)
e2b/agent/              tools (typed registry), state (RunStore), loop (tool-use)
e2b/ui/                 server-rendered interface
e2b/modal_app/          CPU services, GPU worker, serving
tests/                  contract and critical-path checks
docs/                   capabilities, decision log, workstreams, verification, demo script
```

`docs/workstreams.md` defines which lane owns which paths. `docs/capabilities.md` records what
was actually verified in our environment, and the bugs that verification caught.

## Where the LLM comes from

Injected, not imported, so the same loop runs in three places:

| Client | Use |
|---|---|
| `HostLLMClient(host)` | Inside Claude Science; no API key exists in the app. |
| `AnthropicClient()` | From a Modal container; key supplied by a Modal **Secret**. |
| any object with `.complete(messages, tools, system)` | Tests, offline replay. |

## Credentials

Secret **names** only; no values are committed, logged, or exposed to a browser.

| Name | Needed for |
|---|---|
| `Modal` | Census processing and the GPU design worker |
| `ANTHROPIC_API_KEY` | the agent loop when running outside Claude Science |

## Data sources

| Source | Access | Pinned |
|---|---|---|
| Open Targets Platform | GraphQL, `api/v4` | release recorded per query |
| UniProtKB | REST | live |
| AlphaFold DB | REST | model version recorded |
| CELLxGENE Census | `cellxgene-census` Python API, `open_soma` | release pinned per run, never `latest` |

Association scores measure aggregated association evidence — not causality, not probability of
success, and not compatibility with a requested mechanism. Direction of modulation is verified
separately. RNA expression is not surface protein abundance, and undetected expression is not
proof of absence. These caveats are attached to the records themselves, not just written here.
