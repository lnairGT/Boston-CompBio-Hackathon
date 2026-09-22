"""UI-side run state and the store the spine writes into.

The UI owns no scientific logic. It owns a **snapshot** of whatever contracts records
exist for a run, plus UI-only facts the contracts do not model: which section the user
has confirmed, which launch tokens have been consumed, and which sources failed.

The spine has two ways to get real records in front of the user:

1. **Push** — call :meth:`UiStore.upsert` whenever a record is produced::

       store.upsert("run-1", brief=brief, agent_run=agent_run, targets=[t1, t2])

2. **Pull** — register a refresh callable that the UI invokes on every render::

       store.attach_refresher("run-1", lambda key: {"agent_run": spine.agent_run(key)})

Both accept plain contracts objects. Nothing here mutates a record.
"""

from __future__ import annotations

import threading
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Callable, Iterable

from e2b.contracts import (
    AgentRun,
    Assessment,
    Candidate,
    CellSummary,
    ComparisonResult,
    DesignRequest,
    DesignRun,
    Evidence,
    JobStatus,
    NeighborhoodEvidence,
    ResearchBrief,
    Target,
)
from e2b.ui.honesty import SourceFailure

#: Ordered agent stages the UI knows how to name. An unrecognised stage is shown
#: verbatim rather than being forced into this list.
STAGE_ORDER: list[str] = [
    "resolve_indication",
    "gather_targets",
    "assess_targets",
    "cellular_context",
    "rank_targets",
    "structure_and_site",
    "design_request",
    "design_run",
    "report",
]

STAGE_LABELS: dict[str, str] = {
    "resolve_indication": "Resolving the indication",
    "gather_targets": "Gathering associated targets",
    "assess_targets": "Assessing evidence per target",
    "cellular_context": "Retrieving cellular context",
    "rank_targets": "Comparing targets against the rubric",
    "structure_and_site": "Checking structure and binding site",
    "design_request": "Preparing the design request",
    "design_run": "Running binder design",
    "report": "Assembling the report",
}


@dataclass
class LaunchTicket:
    """Records that a launch was requested, so a second click cannot re-spend money."""

    token: str
    requested_at: datetime
    design_id: str | None = None
    run_id: str | None = None
    outcome: str = "pending"  # pending | launched | rejected | no_executor
    message: str | None = None


@dataclass
class RunState:
    """Everything the page needs for one run key."""

    run_key: str
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    updated_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))

    brief: ResearchBrief | None = None
    agent_run: AgentRun | None = None
    targets: dict[str, Target] = field(default_factory=dict)
    assessments: dict[str, Assessment] = field(default_factory=dict)
    comparison: ComparisonResult | None = None
    evidence: dict[str, Evidence] = field(default_factory=dict)
    cell_summaries: list[CellSummary] = field(default_factory=list)
    neighborhoods: list[NeighborhoodEvidence] = field(default_factory=list)
    design_request: DesignRequest | None = None
    design_run: DesignRun | None = None
    candidates: list[Candidate] = field(default_factory=list)

    # UI-only
    selected_target_id: str | None = None
    inputs_confirmed: bool = False
    launch_tickets: dict[str, LaunchTicket] = field(default_factory=dict)
    source_failures: list[SourceFailure] = field(default_factory=list)
    scenario_label: str | None = None
    notice: str | None = None

    # ------------------------------------------------------------------
    def touch(self) -> None:
        self.updated_at = datetime.now(timezone.utc)

    @property
    def origins_present(self) -> set[str]:
        seen: set[str] = set()
        for obj in (
            [self.brief, self.design_run]
            + list(self.targets.values())
            + list(self.evidence.values())
            + list(self.cell_summaries)
            + list(self.candidates)
        ):
            o = getattr(obj, "origin", None)
            if o:
                seen.add(o)
            prov = getattr(obj, "provenance", None)
            if prov is not None and getattr(prov, "origin", None):
                seen.add(prov.origin)
        return seen

    @property
    def is_fixture_run(self) -> bool:
        """True when nothing in this run came from a live source."""
        origins = self.origins_present
        return bool(origins) and origins <= {"fixture"}

    def pending_launch_token(self) -> str:
        """Return the single outstanding launch token, minting one only if none exists.

        Re-rendering the page must NOT mint a fresh token: if it did, two browser tabs
        or a refresh would each hold a usable token and a double submission could start
        two paid runs. One unconsumed token exists at a time, and consuming it (or
        changing the target/confirmation) retires it.
        """
        for t in self.launch_tickets.values():
            if t.outcome == "pending":
                return t.token
        token = uuid.uuid4().hex
        self.launch_tickets[token] = LaunchTicket(
            token=token, requested_at=datetime.now(timezone.utc)
        )
        return token

    def retire_pending_tokens(self) -> None:
        """Invalidate outstanding tokens after an input change."""
        for t in list(self.launch_tickets.values()):
            if t.outcome == "pending":
                t.outcome = "rejected"
                t.message = "Inputs changed after this token was issued; re-confirm to launch."


class UiStore:
    """Thread-safe in-process store of :class:`RunState`.

    Deliberately in-process and unpersisted: the UI is a view, the spine owns the
    durable record. Restarting the web process loses no scientific state.
    """

    def __init__(self) -> None:
        self._lock = threading.RLock()
        self._runs: dict[str, RunState] = {}
        self._refreshers: dict[str, Callable[[str], dict[str, Any]]] = {}
        self._launcher: Callable[[DesignRequest], DesignRun] | None = None
        self._resolver: Callable[[dict[str, Any]], ResearchBrief] | None = None

    # -- registration -------------------------------------------------------
    def attach_refresher(self, run_key: str, fn: Callable[[str], dict[str, Any]]) -> None:
        """Register a callable the UI invokes before each render of ``run_key``.

        It must return a dict whose keys are accepted by :meth:`upsert`.
        """
        with self._lock:
            self._refreshers[run_key] = fn

    def attach_launcher(self, fn: Callable[[DesignRequest], DesignRun]) -> None:
        """Register the spine's design launcher.

        With no launcher attached the Design section shows an explicit
        *no execution route attached* state; it never pretends to launch.
        """
        self._launcher = fn

    def attach_resolver(self, fn: Callable[[dict[str, Any]], ResearchBrief]) -> None:
        """Register the spine's indication resolver, called on brief submission."""
        self._resolver = fn

    @property
    def has_launcher(self) -> bool:
        return self._launcher is not None

    @property
    def has_resolver(self) -> bool:
        return self._resolver is not None

    # -- access -------------------------------------------------------------
    def get(self, run_key: str, *, create: bool = False) -> RunState | None:
        with self._lock:
            state = self._runs.get(run_key)
            if state is None and create:
                state = RunState(run_key=run_key)
                self._runs[run_key] = state
            return state

    def list_runs(self) -> list[RunState]:
        with self._lock:
            return sorted(self._runs.values(), key=lambda s: s.updated_at, reverse=True)

    def delete(self, run_key: str) -> None:
        with self._lock:
            self._runs.pop(run_key, None)
            self._refreshers.pop(run_key, None)

    def refresh(self, run_key: str) -> RunState | None:
        """Pull fresh records from a registered refresher, if any."""
        fn = self._refreshers.get(run_key)
        if fn is None:
            return self.get(run_key)
        try:
            payload = fn(run_key) or {}
        except Exception as exc:  # a broken refresher must not blank the page
            state = self.get(run_key, create=True)
            assert state is not None
            self.note_source_failure(
                run_key, source="spine_refresher", kind="unreachable", detail=str(exc)
            )
            return state
        return self.upsert(run_key, **payload)

    # -- writes -------------------------------------------------------------
    def upsert(
        self,
        run_key: str,
        *,
        brief: ResearchBrief | None = None,
        agent_run: AgentRun | None = None,
        targets: Iterable[Target] | None = None,
        assessments: Iterable[Assessment] | None = None,
        comparison: ComparisonResult | None = None,
        evidence: Iterable[Evidence] | None = None,
        cell_summaries: Iterable[CellSummary] | None = None,
        neighborhoods: Iterable[NeighborhoodEvidence] | None = None,
        design_request: DesignRequest | None = None,
        design_run: DesignRun | None = None,
        candidates: Iterable[Candidate] | None = None,
        selected_target_id: str | None = None,
        scenario_label: str | None = None,
        notice: str | None = None,
        replace_collections: bool = False,
    ) -> RunState:
        """Merge records into a run. Absent arguments leave existing values alone."""
        with self._lock:
            state = self.get(run_key, create=True)
            assert state is not None
            if brief is not None:
                state.brief = brief
            if agent_run is not None:
                state.agent_run = agent_run
            if comparison is not None:
                state.comparison = comparison
            if design_request is not None:
                state.design_request = design_request
            if design_run is not None:
                state.design_run = design_run
            if selected_target_id is not None:
                state.selected_target_id = selected_target_id
            if scenario_label is not None:
                state.scenario_label = scenario_label
            if notice is not None:
                state.notice = notice

            if targets is not None:
                if replace_collections:
                    state.targets = {}
                state.targets.update({t.target_id: t for t in targets})
            if assessments is not None:
                if replace_collections:
                    state.assessments = {}
                state.assessments.update({a.target_id: a for a in assessments})
            if evidence is not None:
                if replace_collections:
                    state.evidence = {}
                state.evidence.update({e.evidence_id: e for e in evidence})
            if cell_summaries is not None:
                state.cell_summaries = (
                    list(cell_summaries) if replace_collections
                    else state.cell_summaries + list(cell_summaries)
                )
            if neighborhoods is not None:
                state.neighborhoods = (
                    list(neighborhoods) if replace_collections
                    else state.neighborhoods + list(neighborhoods)
                )
            if candidates is not None:
                state.candidates = (
                    list(candidates) if replace_collections
                    else state.candidates + list(candidates)
                )
            state.touch()
            return state

    def note_source_failure(
        self,
        run_key: str,
        *,
        source: str,
        kind: str,
        detail: str | None = None,
        query: str | None = None,
    ) -> SourceFailure:
        """Record that a source did not deliver. Never turns into evidence of absence."""
        with self._lock:
            state = self.get(run_key, create=True)
            assert state is not None
            failure = SourceFailure(source=source, kind=kind, detail=detail, query=query)
            state.source_failures.append(failure)
            state.touch()
            return failure

    def set_selection(self, run_key: str, target_id: str | None) -> RunState:
        with self._lock:
            state = self.get(run_key, create=True)
            assert state is not None
            state.selected_target_id = target_id
            state.inputs_confirmed = False  # changing the target invalidates confirmation
            state.retire_pending_tokens()
            state.touch()
            return state

    def set_confirmed(self, run_key: str, confirmed: bool) -> RunState:
        with self._lock:
            state = self.get(run_key, create=True)
            assert state is not None
            if not confirmed:
                state.retire_pending_tokens()
            state.inputs_confirmed = confirmed
            state.touch()
            return state

    # -- launch -------------------------------------------------------------
    def launch(self, run_key: str, token: str) -> LaunchTicket:
        """Consume a launch token exactly once.

        A repeated submission of the same token (double click, refresh, back button)
        returns the ORIGINAL ticket and does not call the launcher again.
        """
        with self._lock:
            state = self.get(run_key, create=True)
            assert state is not None
            ticket = state.launch_tickets.get(token)
            if ticket is None:
                ticket = LaunchTicket(
                    token=token,
                    requested_at=datetime.now(timezone.utc),
                    outcome="rejected",
                    message="Unknown or expired launch token. Re-confirm the inputs to launch.",
                )
                state.launch_tickets[token] = ticket
                return ticket
            if ticket.outcome != "pending":
                return ticket  # idempotent replay

            request = state.design_request
            if request is None:
                ticket.outcome = "rejected"
                ticket.message = "No DesignRequest exists for this run."
                return ticket
            try:
                request.assert_launchable()
            except ValueError as exc:
                ticket.outcome = "rejected"
                ticket.message = str(exc)
                return ticket
            if self._launcher is None:
                ticket.outcome = "no_executor"
                ticket.message = (
                    "No design execution route is attached to this UI process. Nothing was "
                    "submitted and no compute was spent."
                )
                return ticket
            ticket.design_id = request.design_id
            try:
                run = self._launcher(request)
            except Exception as exc:
                ticket.outcome = "rejected"
                ticket.message = f"Launcher raised: {exc}"
                return ticket
            state.design_run = run
            ticket.run_id = run.run_id
            ticket.outcome = "launched"
            ticket.message = f"Submitted as run {run.run_id} via {run.execution_route}."
            state.touch()
            return ticket

    def resolve_brief(self, payload: dict[str, Any]) -> tuple[ResearchBrief | None, str | None]:
        """Ask the registered resolver to turn free text into a ResearchBrief.

        Returns ``(brief, error)``. With no resolver attached this returns
        ``(None, <explanation>)`` — the UI then shows an explicit unresolved state.
        """
        if self._resolver is None:
            return None, (
                "No indication resolver is attached to this UI process, so the free-text "
                "indication has not been resolved to an ontology term. Nothing was guessed."
            )
        try:
            return self._resolver(payload), None
        except Exception as exc:
            return None, f"The resolver failed: {exc}"


_STORE: UiStore | None = None
_STORE_LOCK = threading.Lock()


def get_store() -> UiStore:
    """Process-wide store. The spine and the web app share this instance."""
    global _STORE
    with _STORE_LOCK:
        if _STORE is None:
            _STORE = UiStore()
        return _STORE


def reset_store() -> UiStore:
    """Drop all run state. For tests."""
    global _STORE
    with _STORE_LOCK:
        _STORE = UiStore()
        return _STORE


TERMINAL_JOB_STATES = {
    JobStatus.SUCCEEDED,
    JobStatus.FAILED,
    JobStatus.TIMED_OUT,
    JobStatus.CANCELLED,
}

__all__ = [
    "LaunchTicket",
    "RunState",
    "STAGE_LABELS",
    "STAGE_ORDER",
    "UiStore",
    "get_store",
    "reset_store",
]
