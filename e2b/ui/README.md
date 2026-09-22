# Lane 5 — review UI

A server-rendered single page over the evidence-to-binder workflow. FastAPI + Jinja2 +
one plain stylesheet. No build step, no npm, no client framework, no credential in a
template. The page is fully readable with JavaScript disabled; the ~40 lines of vanilla
script only disable a submit button and poll one status fragment.

## Run it

```bash
pip install fastapi jinja2 'pydantic>=2' uvicorn python-multipart
uvicorn e2b.ui.app:app --reload --port 8000
# http://127.0.0.1:8000
```

Fixture scenarios are preloaded at boot. `E2B_UI_PRELOAD_FIXTURES=0` serves a clean
instance.

## Check it

```bash
python -m e2b.ui.selfcheck     # 17 assertions over the RENDERED HTML
```

## Wire the spine to it

```python
from e2b.ui import get_store

store = get_store()                    # same process as the web app

# 1. Push records as they are produced. Absent arguments leave existing values alone.
store.upsert(
    run_key,
    brief=brief,                       # ResearchBrief
    agent_run=agent_run,               # AgentRun
    targets=[t1, t2, t3],              # Target
    assessments=[a1, a2, a3],          # Assessment
    comparison=comparison,             # ComparisonResult
    evidence=evidence_items,           # Evidence
    cell_summaries=[cs],               # CellSummary
    neighborhoods=[nb],                # NeighborhoodEvidence
    design_request=req,                # DesignRequest
    design_run=run,                    # DesignRun
    candidates=cands,                  # Candidate
    replace_collections=True,          # replace rather than append
)

# 2. Or let the UI pull on every render.
store.attach_refresher(run_key, lambda key: {"agent_run": spine.agent_run(key)})

# 3. Free-text brief submission calls this. Return a ResearchBrief or raise.
store.attach_resolver(lambda payload: spine.resolve_brief(**payload))

# 4. The launch control calls this, and only after assert_launchable() passes.
store.attach_launcher(lambda design_request: spine.launch(design_request))

# 5. A source that did not deliver. Never rendered as evidence of absence.
store.note_source_failure(run_key, source="cellxgene_census", kind="timeout",
                          detail="no response in 60 s", query="<the query issued>")
```

`kind` is one of `timeout`, `http_error`, `unreachable`, `not_configured`,
`empty_result`, `unsupported`. Only `empty_result` is rendered as telling you anything
about the world, and only about the query as issued.

## Routes

| Route | Purpose |
|---|---|
| `GET /` | brief form, fixture list, runs in this process |
| `POST /brief` | free text in; calls the resolver, or records an unresolved brief |
| `GET /run/{key}` | the page — five sections, progressive detail |
| `GET /run/{key}/state.json` | the same view model as JSON (Lane 4 reads this) |
| `POST /run/{key}/disease` | choose among the resolver's candidates |
| `POST /run/{key}/select` | select a target for design |
| `POST /run/{key}/confirm` | confirm inputs; enables the launch control |
| `POST /run/{key}/launch` | one-shot token; a repeat returns the original outcome |
| `GET /run/{key}/fragment/job` | job-status fragment for polling |
| `GET /fixtures/{name}` | load a labelled fixture scenario |
| `GET /healthz` | liveness plus which hooks are attached |

## What this layer refuses to do

* Guess an ontology term. An unresolved brief renders as unresolved.
* Render a value the record does not carry. `None` is *not measured*, never `0`.
* Present a replay as a new run. `cached_real` renders its original completion time.
* Draw a progress bar. The executor reports a stage, so the page shows a stage.
* Call a designed minibinder an antibody, or a model score an affinity.
* Turn a timeout into "no evidence found".

## Files

```
e2b/ui/__init__.py        public surface: create_app, get_store, UiStore
e2b/ui/app.py             FastAPI routes
e2b/ui/state.py           RunState, UiStore, launch tokens
e2b/ui/viewmodels.py      record -> view dict, pure
e2b/ui/honesty.py         badges, numerics, vocabulary guards, SourceFailure
e2b/ui/selfcheck.py       17 checks over the rendered HTML
e2b/ui/fixtures/          seven labelled scenarios
e2b/ui/templates/         base, run page, five section partials, macros
e2b/ui/static/            app.css, app.js
e2b/modal_app/web.py      Modal ASGI entrypoint (import-checked, NOT deployed)
docs/demo_script.md       timed three-minute script
```
