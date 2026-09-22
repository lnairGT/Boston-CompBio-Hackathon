"""FastAPI application: one server-rendered page with progressive detail.

Run locally::

    uvicorn e2b.ui.app:app --reload --port 8000

No build step, no npm, no client-side framework, no API key ever reaches a template.
The only JavaScript is ~40 lines of vanilla script that disables a submit button and
polls one HTML fragment; the page is fully readable with JavaScript disabled.
"""

from __future__ import annotations

import os
import uuid
from pathlib import Path
from typing import Any

from fastapi import APIRouter, FastAPI, Form, Request
from fastapi.responses import HTMLResponse, JSONResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from e2b.config import DEMO_PRESETS
from e2b.contracts import ResearchBrief, utcnow
from e2b.ui.fixtures import available_scenarios, load_scenario
from e2b.ui.state import UiStore, get_store
from e2b.ui.viewmodels import job_view, page_view

HERE = Path(__file__).resolve().parent
TEMPLATES = Jinja2Templates(directory=str(HERE / "templates"))

#: Presets only ever PREFILL the form. The form accepts any free text; nothing in the
#: request path consults this list, and deleting it changes no behaviour.
PRESETS = DEMO_PRESETS


def _store(request: Request) -> UiStore:
    return request.app.state.store


def _page(request: Request, run_key: str, status_code: int = 200) -> HTMLResponse:
    store = _store(request)
    store.refresh(run_key)
    state = store.get(run_key)
    if state is None:
        return TEMPLATES.TemplateResponse(
            request=request,
            name="missing.html",
            context={"run_key": run_key, "runs": store.list_runs(),
                     "scenarios": available_scenarios()},
            status_code=404,
        )
    ctx = {
        "v": page_view(state),
        "presets": PRESETS,
        "scenarios": available_scenarios(),
        "runs": store.list_runs(),
        "has_launcher": store.has_launcher,
        "has_resolver": store.has_resolver,
        "launch_token": (
            state.pending_launch_token() if state.inputs_confirmed else None
        ),
    }
    return TEMPLATES.TemplateResponse(
        request=request, name="run.html", context=ctx, status_code=status_code
    )


router = APIRouter()


@router.get("/", response_class=HTMLResponse)
def index(request: Request) -> Any:
    store = _store(request)
    return TEMPLATES.TemplateResponse(
        request=request,
        name="index.html",
        context={
            "presets": PRESETS,
            "scenarios": available_scenarios(),
            "runs": store.list_runs(),
            "has_resolver": store.has_resolver,
            "has_launcher": store.has_launcher,
        },
    )


@router.get("/healthz")
def healthz(request: Request) -> JSONResponse:
    store = _store(request)
    return JSONResponse(
        {
            "ok": True,
            "runs": len(store.list_runs()),
            "resolver_attached": store.has_resolver,
            "launcher_attached": store.has_launcher,
        }
    )


@router.post("/brief")
def submit_brief(
    request: Request,
    raw_indication: str = Form(...),
    desired_effect: str = Form(""),
    moa: str = Form(""),
    modality: str = Form(""),
    species_taxon_id: str = Form("NCBITaxon:9606"),
) -> Any:
    """Accept free text and hand it to the spine's resolver.

    The UI never guesses an ontology term. With no resolver attached the run is created
    with an explicitly unresolved brief and the page says so.
    """
    store = _store(request)
    run_key = f"run-{uuid.uuid4().hex[:8]}"
    payload = {
        "raw_indication": raw_indication.strip(),
        "desired_effect": desired_effect.strip() or None,
        "moa": moa.strip() or None,
        "modality": modality.strip() or None,
        "species_taxon_id": species_taxon_id.strip() or "NCBITaxon:9606",
        "run_key": run_key,
    }
    brief, error = store.resolve_brief(payload)
    if brief is None:
        brief = ResearchBrief(
            brief_id=f"brief-{run_key}",
            raw_indication=payload["raw_indication"],
            desired_effect=payload["desired_effect"],
            moa=payload["moa"],
            modality=payload["modality"],
            species_taxon_id=payload["species_taxon_id"],
            status="pending",
        )
        store.upsert(run_key, brief=brief, notice=error)
        store.note_source_failure(
            run_key,
            source="indication_resolver",
            kind="not_configured" if not store.has_resolver else "unreachable",
            detail=error,
            query=payload["raw_indication"],
        )
    else:
        store.upsert(run_key, brief=brief)
    return RedirectResponse(url=f"/run/{run_key}", status_code=303)


@router.get("/run/{run_key}", response_class=HTMLResponse)
def run_page(request: Request, run_key: str) -> Any:
    return _page(request, run_key)


@router.get("/run/{run_key}/state.json")
def run_state_json(request: Request, run_key: str) -> JSONResponse:
    """Machine-readable view model. Used by Lane 4's checks."""
    store = _store(request)
    store.refresh(run_key)
    state = store.get(run_key)
    if state is None:
        return JSONResponse({"error": "unknown run"}, status_code=404)
    return JSONResponse(_jsonable(page_view(state)))


def _jsonable(obj: Any) -> Any:
    from e2b.ui.honesty import Num

    if isinstance(obj, Num):
        return {"text": obj.text, "known": obj.known, "raw": obj.raw}
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if hasattr(obj, "isoformat"):
        return obj.isoformat()
    return obj


@router.post("/run/{run_key}/disease")
def choose_disease(request: Request, run_key: str, disease_id: str = Form(...)) -> Any:
    """Record the user's choice among the resolver's candidates.

    The UI only ever selects from candidates the resolver returned; it never mints an
    ontology identifier.
    """
    store = _store(request)
    state = store.get(run_key)
    if state is None or state.brief is None:
        return RedirectResponse(url=f"/run/{run_key}", status_code=303)
    match = next(
        (c for c in state.brief.resolution_candidates if c.disease_id == disease_id), None
    )
    if match is None:
        store.upsert(
            run_key,
            notice=f"'{disease_id}' is not one of the resolver's candidates; nothing changed.",
        )
        return RedirectResponse(url=f"/run/{run_key}", status_code=303)
    brief = state.brief
    brief.resolved_disease_id = match.disease_id
    brief.resolved_disease_label = match.label
    brief.status = "confirmed"
    brief.resolution_rationale = (
        (brief.resolution_rationale or "")
        + f" Operator selected '{match.label}' ({match.disease_id}) at {utcnow().isoformat()}."
    ).strip()
    store.upsert(run_key, brief=brief, notice=f"Disease identity set to {match.label}.")
    return RedirectResponse(url=f"/run/{run_key}#brief", status_code=303)


@router.post("/run/{run_key}/select")
def select_target(request: Request, run_key: str, target_id: str = Form(...)) -> Any:
    store = _store(request)
    store.set_selection(run_key, target_id)
    return RedirectResponse(url=f"/run/{run_key}#design", status_code=303)


@router.post("/run/{run_key}/confirm")
def confirm_inputs(request: Request, run_key: str, confirm: str = Form("")) -> Any:
    store = _store(request)
    store.set_confirmed(run_key, confirm == "yes")
    return RedirectResponse(url=f"/run/{run_key}#design", status_code=303)


@router.post("/run/{run_key}/launch")
def launch(request: Request, run_key: str, launch_token: str = Form(...)) -> Any:
    """Consume a one-shot launch token.

    A duplicate submission of the same token returns the original outcome without
    calling the launcher again, so a double click cannot start two paid runs.
    """
    store = _store(request)
    store.launch(run_key, launch_token)
    return RedirectResponse(url=f"/run/{run_key}#design", status_code=303)


@router.get("/run/{run_key}/fragment/job", response_class=HTMLResponse)
def job_fragment(request: Request, run_key: str) -> Any:
    """Small HTML fragment for polling real job state. No progress percentage."""
    store = _store(request)
    store.refresh(run_key)
    state = store.get(run_key)
    if state is None:
        return HTMLResponse("<p class='muted'>Run not found.</p>", status_code=404)
    return TEMPLATES.TemplateResponse(
        request=request,
        name="partials/job_status.html",
        context={"job": job_view(state.design_run), "run_key": run_key},
    )


@router.post("/fixtures/{name}")
@router.get("/fixtures/{name}")
def load_fixture(request: Request, name: str) -> Any:
    store = _store(request)
    try:
        key = load_scenario(store, name)
    except KeyError as exc:
        return JSONResponse({"error": str(exc)}, status_code=404)
    return RedirectResponse(url=f"/run/{key}", status_code=303)


def create_app(store: UiStore | None = None) -> FastAPI:
    """Build the app. Pass a store to share state with an in-process spine."""
    app = FastAPI(
        title="CELLxGENE evidence-to-binder",
        description="Server-rendered review surface over the evidence-to-binder workflow.",
        docs_url="/api-docs",
    )
    app.state.store = store or get_store()
    app.mount("/static", StaticFiles(directory=str(HERE / "static")), name="static")
    app.include_router(router)

    if os.environ.get("E2B_UI_PRELOAD_FIXTURES", "1") != "0":
        for s in available_scenarios():
            try:
                load_scenario(app.state.store, s["key"])
            except Exception:  # pragma: no cover - fixtures must never break boot
                pass
    return app


app = create_app()

__all__ = ["app", "create_app", "router"]
