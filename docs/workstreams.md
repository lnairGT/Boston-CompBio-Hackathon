# Workstreams, file ownership and handoffs

Boston Computational Biology Hackathon · 22 September 2026

One integration owner controls the shared schema and the deployment entrypoint. Everything
else is owned by exactly one workstream. **The rule that keeps seven parallel lanes from
colliding: you may create and edit files inside your own paths, and you may only *read*
files owned by another lane.** If you need a change to a shared file, ask the integration
owner rather than editing it.

## Shared, integration-owned — read-only to every lane

| Path | Contents |
|---|---|
| `e2b/contracts.py` | Every typed record. `schema_version`, provenance, origin labels. |
| `e2b/config.py` | Budgets, the versioned ranking rubric, modality accessibility screen. |
| `e2b/pipeline/` | Evidence assembly, deterministic ranking, run coordination. |
| `e2b/agent/` | The Claude tool-use loop, tool registry, `AgentRun` persistence. |

Changes here are proposed to the integration owner, not pushed directly. The rubric is
versioned (`RUBRIC_VERSION`); changing a weight means bumping it, because rankings recorded
under an old rubric must stay interpretable.

## Lanes

| # | Lane | Owns (create/edit) | Reads | First handoff |
|---|---|---|---|---|
| 1 | **Integration / spine** | `contracts.py`, `config.py`, `pipeline/`, `agent/`, `adapters/open_targets.py`, `adapters/structures.py` | everything | Typed contracts + Open Targets and structure adapters |
| 2 | **Cellular evidence** | `adapters/cellxgene.py`, `modal_app/census_worker.py` | contracts, config | Versioned cell summary + aligned embedding slice + coverage report |
| 3 | **Binder design** | `adapters/protein_design.py`, `modal_app/design_worker.py` | contracts, config | Validated design request + real job status |
| 4 | **Evaluation / QA** | `tests/`, `e2b/eval/`, `docs/verification.md` | everything (read-only) | Critical-path check results, pass/fail with evidence |
| 5 | **Experience / UI** | `e2b/ui/`, `modal_app/web.py`, `docs/demo_script.md` | contracts only | Complete flow on labelled fixtures, then real records |
| 6 | **Downstream analysis** | `e2b/analysis/`, `docs/analysis_methods.md` | contracts, candidates | Candidate triage report with metric definitions |
| 7 | **Experimental design** | `e2b/experiments/`, `docs/experimental_plan.md` | contracts, candidates | Binding-validation and MoA functional assay concepts |

## Handoff contracts between lanes

Lanes talk through **typed records, not function calls into each other's modules**. The
records are defined once in `e2b/contracts.py`:

```
Lane 2 → spine        CellSlice, EmbeddingSlice, CellSummary, NeighborhoodEvidence
spine → Lane 3        DesignRequest  (review_status must be "approved")
Lane 3 → spine        DesignRun, Candidate, Metric
spine → Lane 5        ResearchBrief, ComparisonResult, Assessment, AgentRun
Lane 3 → Lane 6, 7    Candidate  (never a raw model output dict)
Lane 6 → Lane 7       triage ranking, so assays are designed for candidates worth testing
Lane 4 → everyone     failing checks, as issues against the owning lane
```

A lane that needs a record shape that does not exist yet builds against the contract and
returns plain JSON-able dicts at its boundary. It does **not** block waiting for another
lane, and it does **not** invent a new parallel record type for the same concept.

## Branch and PR discipline

- Branch per lane: `lane/<n>-<short-name>`, e.g. `lane/4-evaluation`.
- Small PRs with module-specific files only. A PR touching another lane's files gets
  bounced, not merged.
- Never force-push a shared branch. Never commit secrets, weights, datasets or generated
  structures — commit small manifests and reference run artifacts by run ID.
- Schema changes are coordinated before merge, because every lane validates at its boundary.

## Origin labels are not decoration

Every record carries `origin`: `real`, `cached_real` or `fixture`. Lane 5 legitimately
starts on fixtures to unblock the UI, and that is expected — but a `fixture` record may
never reach a scientific conclusion, a ranking, or a design run, and the UI must show its
label. Lane 4 tests this boundary specifically. A cached result shows its original run time
and is never presented as a newly completed model run.
