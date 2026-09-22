"""Bridge between the agent's persisted ``RunStore`` and the review UI's store.

The UI is deliberately a *reader*. It never computes a score, never re-derives a rating,
and never decides whether a design may launch -- it renders records the spine produced and
shows their origin label. This module is the one place those two stores meet.

It is also where a source failure becomes visible. When a tool returns
``source_unavailable`` or ``capability_unavailable``, that is pushed to the UI as a
``SourceFailure``, not silently dropped -- otherwise a page with no cellular-context panel
looks identical whether the target is unexpressed or the Census was simply unreachable.
"""

from __future__ import annotations

from typing import Any, Callable

from ..agent.state import RunStore


def publish(store: RunStore, ui_store: Any, run_key: str | None = None) -> str:
    """Push everything the spine currently holds into the UI store. Idempotent."""
    key = run_key or store.run.run_id
    evidence: list[Any] = []
    for items in store.evidence.values():
        evidence.extend(items)

    ui_store.upsert(
        key,
        brief=store.brief,
        agent_run=store.run,
        targets=list(store.targets.values()),
        assessments=list(store.assessments.values()) or None,
        comparison=store.comparison,
        evidence=evidence or None,
        cell_summaries=store.cell_summaries or None,
        design_request=next(iter(store.design_requests.values()), None),
        design_run=next(iter(store.design_runs.values()), None),
        replace_collections=True,
    )

    # Surface every failed retrieval, so "nothing shown" is never ambiguous.
    #
    # Two distinct cases, and missing either one lets an empty panel read as a negative
    # finding. A raised error lands on `action.error`; a capability that is simply absent
    # returns a well-formed result whose *status* says so, with no error at all. The
    # second case is the one that matters most here -- "no cellular context shown" must
    # never be mistaken for "the target is not expressed".
    for action in store.run.tool_action_log:
        summary = action.result_summary or ""
        unavailable = "capability_unavailable" in summary or "source_unavailable" in summary
        if not action.error and not unavailable:
            continue
        detail = action.error or summary
        kind = (
            "unavailable" if "capability_unavailable" in detail or "capability_unavailable" in summary
            else "timeout" if "timeout" in detail.lower()
            else "unreachable"
        )
        ui_store.note_source_failure(
            key, source=action.tool_name, kind=kind, detail=detail[:400],
            query=str(action.arguments)[:300],
        )
    return key


def attach(store: RunStore, ui_store: Any, *, launcher: Callable[..., Any] | None = None,
           resolver: Callable[..., Any] | None = None) -> Callable[[], None]:
    """Wire the spine into a UI store and return a callback for ``AgentLoop(publish=...)``."""
    if resolver is not None:
        ui_store.attach_resolver(resolver)
    if launcher is not None:
        ui_store.attach_launcher(launcher)

    def _publish() -> None:
        try:
            publish(store, ui_store)
        except Exception:  # a display failure must not abort an investigation
            pass

    _publish()
    return _publish
