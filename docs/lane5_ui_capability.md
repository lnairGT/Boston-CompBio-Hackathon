# Lane 5 — experience / UI · capability report

Boston Computational Biology Hackathon · 22 September 2026
Scope: `e2b/ui/`, `e2b/modal_app/web.py`, `docs/demo_script.md`. No file outside these
paths was created or modified.

## What is proven

A single-page, server-rendered review surface over the evidence-to-binder workflow,
built on FastAPI + Jinja2 + one stylesheet. No build step, no npm, no client framework,
no credential in any template. **17 of 17 honesty checks pass against the rendered
HTML**, not against the source:

```
$ python -m e2b.ui.selfcheck
pass  every scenario renders
pass  fixture records are visibly labelled
pass  no synthetic progress indicator anywhere
pass  unknown numerics render as text, never as 0
pass  cached results show their original run time and are not called new
pass  pending, no-data, failed and timed-out states are all reachable
pass  a retrieval timeout is presented as a source failure, not absence
pass  a minibinder is never called an antibody
pass  model scores are never labelled as measured affinity
pass  design feasibility is shown separately from biological support
pass  insufficient evidence is shown rather than a forced winner
pass  launch control is disabled until inputs are confirmed
pass  a second click cannot start a second run
pass  presets prefill but never restrict the form
pass  an unresolved brief is shown as unresolved, not guessed
pass  real records render a real badge and no fixture banner
pass  no secret-shaped material reaches a template

17 passed, 0 failed, 17 total
```

### The five sections

| Section | What it shows | Honesty behaviour |
|---|---|---|
| **Brief** | free-text indication and desired effect, resolved disease identity, every ontology candidate with match type and score, Census disease mapping, current agent stage and next action | unresolved renders as *unresolved* — no term is ever guessed; the ambiguity table is a radio choice, not a silent default; demo presets prefill the form and cannot restrict it (the field is a free text input, there is no `<select>`) |
| **Compare** | up to three target cards: four rubric criteria with ratings and evidence counts, gaps, liabilities, direction compatibility, rubric score, excluded targets with reasons | design feasibility sits in its own block below a rule and is never one of the scored criteria; ties render as ties; `insufficient_evidence` renders a banner instead of a winner; a criterion rated `unknown` renders *no evidence retrieved* plus what was looked for |
| **Inspect** | evidence trail grouped by target — predicate triple, excerpt, source link, locator, disease scope, limitations; per-donor cellular context; neighbourhood context; graph edge table; source-status table | provenance is behind a disclosure, not in the main path; propagated evidence is flagged; a timeout renders in a *source status* table with the text "This is a retrieval failure, not evidence that the data is absent" |
| **Design** | selected target, structure, binding region with hotspots and MoA rationale, linked evidence ids, request review status, real job status | the launch button is `disabled` until a human ticks the confirm box **and** `DesignRequest.assert_launchable()` passes; the blocking reasons are listed verbatim; one outstanding one-shot token per run, so a refresh, a second tab or a double click cannot start a second paid run |
| **Results** | candidate sequence and length, structure reference, model scores, metric definitions, limitations, proposed next experimental test | metrics are captioned as computational scores with the producing model named; a null metric renders *not computed*, never 0; a sequence-less candidate renders a no-data state |

### Job states

All four required states are reachable and each has a fixture URL:
`pending` (`/fixtures/pending`), no-data (`/fixtures/no_data`), `failed`
(`/fixtures/failed`), `timed_out` (`/fixtures/timed_out`). A failed or timed-out run
renders "An execution failure says nothing about whether the target is designable."

### Cached results

`/fixtures/cached` replays a run whose original completion is 27 hours old. The panel
shows a `cached` badge, the original `finished_at`, the age, and the sentence "Replay of
a run that originally completed at <time>. Not a new model run."

### No fabricated progress

The stage strip shows the named stage the run reports, a timestamp, an age, and the
real budget counters (`14 of 30 tool decisions used`). There is no `<progress>` element,
no percentage and no CSS width animation anywhere in the rendered output; the self-check
asserts this by regex over every page.

## Running it

```bash
pip install fastapi jinja2 'pydantic>=2' uvicorn python-multipart
uvicorn e2b.ui.app:app --reload --port 8000     # http://127.0.0.1:8000
python -m e2b.ui.selfcheck                      # 17 checks
```

## Modal serving entrypoint

`e2b/modal_app/web.py` defines `modal.App("evidence-to-binder-ui")` with a
`@modal.asgi_app()` function. Verified locally with **modal 1.5.5**: the module imports,
the image builds as an object, and the function registers (`registered_functions ==
['ui']`). **It has not been deployed and no Modal URL exists.** Nothing in the demo
should claim otherwise.

It needs: the `e2b` package on the image (added with `add_local_python_source`, which
carries `ui/templates` and `ui/static`); `fastapi`, `jinja2`, `pydantic>=2`, `requests`;
and **no secret** — the UI holds no credential and renders none. Pin `max_containers=1`
because the run store is in-memory per container.

## Integration surface for the spine

```python
from e2b.ui import create_app, get_store
store = get_store()
store.upsert(run_key, brief=..., agent_run=..., targets=[...], assessments=[...],
             comparison=..., evidence=[...], cell_summaries=[...], neighborhoods=[...],
             design_request=..., design_run=..., candidates=[...],
             replace_collections=True)
store.attach_resolver(fn)     # dict -> ResearchBrief, called on brief submission
store.attach_launcher(fn)     # DesignRequest -> DesignRun, called by the launch control
store.attach_refresher(key, fn)   # key -> dict of records, called on every render
store.note_source_failure(run_key, source=..., kind=..., detail=..., query=...)
```

All arguments are `e2b.contracts` records. The UI defines no parallel record type; its
only local dataclasses are `SourceFailure` (a source that did not deliver — the
contracts do not model this) and `LaunchTicket`/`RunState` (UI session state).
`GET /run/{key}/state.json` returns the whole view model for Lane 4.

## Limits and what is not proven

* **Every record in the shipped fixtures is `origin="fixture"`.** The real-record path
  is exercised by the self-check (`real records render a real badge and no fixture
  banner`) using synthetic `origin="real"` records, but no live Open Targets, Census or
  design record has been rendered yet — the spine has not handed any over.
* **Not deployed to Modal.** Object graph verified; deployment not attempted.
* **The store is in-memory and per-process.** A restart loses run state. This is
  deliberate — the spine owns the durable record — but it means the UI and the spine
  must run in the same process, or the spine must register a refresher.
* No accessibility audit, no cross-browser testing, no responsive testing below ~320 px.
* The evidence trail renders every evidence item for a target with no pagination; a run
  with thousands of items will produce a very long page.
* `POST /brief` has no CSRF token and no authentication. Acceptable for a laptop demo,
  not for an exposed deployment.
