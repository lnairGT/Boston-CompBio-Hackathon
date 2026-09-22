"""Executable checks of the UI honesty rules.

Run::

    python -m e2b.ui.selfcheck

Every check is an assertion about the RENDERED HTML or the view model, not about the
source, so it catches a template that quietly drops a label. Lane 4 can import
:func:`run_checks` and fold the results into its own report.
"""

from __future__ import annotations

import re
import sys
from datetime import datetime, timedelta, timezone
from typing import Any, Callable

from fastapi.testclient import TestClient

from e2b.contracts import (
    Assessment,
    Candidate,
    CriterionRating,
    DesignFeasibility,
    DesignRun,
    Evidence,
    JobStatus,
    Metric,
    Provenance,
    ResearchBrief,
    Target,
)
from e2b.ui.app import create_app
from e2b.ui.fixtures import SCENARIOS
from e2b.ui.state import reset_store

CHECKS: list[tuple[str, Callable[[TestClient, Any], None]]] = []


def check(name: str) -> Callable[[Callable], Callable]:
    def deco(fn: Callable) -> Callable:
        CHECKS.append((name, fn))
        return fn

    return deco


def _pages(client: TestClient) -> dict[str, str]:
    out = {}
    for key in SCENARIOS:
        r = client.get(f"/fixtures/{key}", follow_redirects=True)
        assert r.status_code == 200, f"{key} rendered {r.status_code}"
        out[key] = r.text
    return out


# ---------------------------------------------------------------------------------


@check("every scenario renders")
def _c_render(client: TestClient, store: Any) -> None:
    pages = _pages(client)
    assert len(pages) == len(SCENARIOS)
    for key, html in pages.items():
        assert "<section class=\"card\" id=\"brief\"" in html, f"{key}: no Brief section"
        for sec in ("compare", "inspect", "design", "results"):
            assert f'id="{sec}"' in html, f"{key}: no {sec} section"


@check("fixture records are visibly labelled")
def _c_fixture_badge(client: TestClient, store: Any) -> None:
    for key in SCENARIOS:
        html = client.get(f"/run/fixture-{key}").text
        assert "origin-fixture" in html, f"{key}: no fixture badge in rendered page"
        assert "Fixture run" in html or "fixture" in html.lower()


@check("no synthetic progress indicator anywhere")
def _c_no_progress(client: TestClient, store: Any) -> None:
    bad = (re.compile(r"<progress"), re.compile(r"width:\s*\d+%"),
           re.compile(r"\d+%\s*(complete|done|progress)", re.I))
    for key in SCENARIOS:
        html = client.get(f"/run/fixture-{key}").text
        for pat in bad:
            assert not pat.search(html), f"{key}: progress indicator matched {pat.pattern}"


@check("unknown numerics render as text, never as 0")
def _c_unknown_not_zero(client: TestClient, store: Any) -> None:
    html = client.get("/run/fixture-walkthrough").text
    assert "not computed" in html, "a null metric did not render as 'not computed'"
    assert "not measurable in this donor" in html or "not measured" in html
    j = client.get("/run/fixture-walkthrough/state.json").json()
    for c in j["results"]["candidates"]:
        for m in c["metrics"]:
            if m["value"]["raw"] is None:
                assert m["value"]["known"] is False
                assert m["value"]["text"] != "0"


@check("cached results show their original run time and are not called new")
def _c_cached(client: TestClient, store: Any) -> None:
    html = client.get("/run/fixture-cached").text
    assert "origin-cached" in html, "no cached badge"
    assert "Not a new model run" in html, "cached result does not disclaim being new"
    j = client.get("/run/fixture-cached/state.json").json()
    job = j["design"]["job"]
    assert job["replay"] is True
    assert job["finished"] and "not finished" not in job["finished"]
    assert "originally completed at" in (job["replay_note"] or "")


@check("pending, no-data, failed and timed-out states are all reachable")
def _c_states(client: TestClient, store: Any) -> None:
    want = {
        "pending": "pending",
        "no_data": "no_data",
        "failed": "failed",
        "timed_out": "timed_out",
    }
    for scenario, expect in want.items():
        j = client.get(f"/run/fixture-{scenario}/state.json").json()
        got = j["results"]["status"]
        assert got == expect, f"{scenario}: results status was '{got}', expected '{expect}'"
    html = client.get("/run/fixture-no_data").text
    assert "nodata" in html


@check("a retrieval timeout is presented as a source failure, not absence")
def _c_timeout_language(client: TestClient, store: Any) -> None:
    j = client.get("/run/fixture-timed_out/state.json").json()
    kinds = {f["kind"] for f in j["inspect"]["failures"]}
    assert "timeout" in kinds
    for f in j["inspect"]["failures"]:
        if f["kind"] == "timeout":
            assert f["is_absence_evidence"] is False
            assert "not evidence" in f["meaning"]
    html = client.get("/run/fixture-timed_out").text
    assert "says nothing about whether the target is designable" in html


@check("a minibinder is never called an antibody")
def _c_no_antibody(client: TestClient, store: Any) -> None:
    for key in SCENARIOS:
        html = client.get(f"/run/fixture-{key}").text
        for m in re.finditer(r"antibod\w*", html, re.I):
            window = html[max(0, m.start() - 40): m.start()]
            assert "not an " in window.lower() or "is not a" in window.lower(), (
                f"{key}: 'antibody' used without a negation: ...{window[-40:]}"
                f"{html[m.start():m.end()]}..."
            )


@check("model scores are never labelled as measured affinity")
def _c_no_affinity(client: TestClient, store: Any) -> None:
    j = client.get("/run/fixture-walkthrough/state.json").json()
    for d in j["results"]["metric_definitions"]:
        assert "not a measured" in d["disclaimer"].lower() or \
               "not an experimental" in d["disclaimer"].lower()
    html = client.get("/run/fixture-walkthrough").text
    assert "None of them is a measured binding affinity" in html


@check("design feasibility is shown separately from biological support")
def _c_feasibility_separate(client: TestClient, store: Any) -> None:
    html = client.get("/run/fixture-walkthrough").text
    assert "Design feasibility" in html and "separate from biological support" in html
    j = client.get("/run/fixture-walkthrough/state.json").json()
    for card in j["compare"]["cards"]:
        keys = {c["criterion"] for c in card["criteria"]}
        assert "design_feasibility" not in keys, "feasibility leaked into the biology criteria"
        assert card["feasibility"] is not None


@check("insufficient evidence is shown rather than a forced winner")
def _c_insufficient(client: TestClient, store: Any) -> None:
    j = client.get("/run/fixture-no_data/state.json").json()
    assert j["compare"]["insufficient_evidence"] is True
    assert j["compare"]["cards"] == []
    html = client.get("/run/fixture-no_data").text
    assert "Insufficient evidence to choose" in html


@check("launch control is disabled until inputs are confirmed")
def _c_launch_gate(client: TestClient, store: Any) -> None:
    key = "fixture-walkthrough"
    client.post(f"/run/{key}/confirm", data={"confirm": ""}, follow_redirects=True)
    html = client.get(f"/run/{key}").text
    block = html[html.index('action="/run/' + key + '/launch"'):]
    assert "disabled" in block[:400], "launch button was enabled without confirmation"
    client.post(f"/run/{key}/confirm", data={"confirm": "yes"}, follow_redirects=True)
    html = client.get(f"/run/{key}").text
    block = html[html.index('action="/run/' + key + '/launch"'):]
    assert "disabled" not in block[:400], "launch button stayed disabled after confirmation"


@check("a second click cannot start a second run")
def _c_double_click(client: TestClient, store: Any) -> None:
    calls: list[str] = []

    def launcher(req: Any) -> DesignRun:
        calls.append(req.design_id)
        return DesignRun(
            run_id=f"selfcheck-run-{len(calls)}",
            design_id=req.design_id,
            execution_route="selfcheck:none",
            status=JobStatus.QUEUED,
            started_at=datetime.now(timezone.utc),
            updated_at=datetime.now(timezone.utc),
            origin="fixture",
        )

    store.attach_launcher(launcher)
    key = "fixture-pending"
    client.post(f"/run/{key}/confirm", data={"confirm": "yes"}, follow_redirects=True)
    html = client.get(f"/run/{key}").text
    tok = re.search(r'name="launch_token" value="([0-9a-f]{32})"', html).group(1)
    html2 = client.get(f"/run/{key}").text
    tok2 = re.search(r'name="launch_token" value="([0-9a-f]{32})"', html2).group(1)
    assert tok == tok2, "re-rendering the page minted a second usable launch token"
    for _ in range(3):
        client.post(f"/run/{key}/launch", data={"launch_token": tok}, follow_redirects=True)
    assert len(calls) == 1, f"launcher called {len(calls)} times for one token"
    store._launcher = None


@check("presets prefill but never restrict the form")
def _c_presets(client: TestClient, store: Any) -> None:
    html = client.get("/").text
    assert "prefill only" in html
    assert "<select" not in html, "the indication field must not be a closed list"
    r = client.post(
        "/brief",
        data={"raw_indication": "an indication that is in no preset list",
              "desired_effect": "an arbitrary effect", "moa": "agonise",
              "modality": "de novo minibinder"},
        follow_redirects=True,
    )
    assert r.status_code == 200
    assert "an indication that is in no preset list" in r.text


@check("an unresolved brief is shown as unresolved, not guessed")
def _c_unresolved(client: TestClient, store: Any) -> None:
    r = client.post(
        "/brief",
        data={"raw_indication": "a condition with no resolver attached",
              "desired_effect": "", "moa": "", "modality": "de novo minibinder"},
        follow_redirects=True,
    )
    assert "has <strong>not</strong> been resolved" in r.text or "not</strong> been resolved" in r.text
    assert "No term was guessed" in r.text


@check("real records render a real badge and no fixture banner")
def _c_real_path(client: TestClient, store: Any) -> None:
    key = "selfcheck-real"
    prov = Provenance(
        source="open_targets",
        source_release="26.06",
        endpoint="https://api.platform.opentargets.org/api/v4/graphql",
        origin="real",
    )
    brief = ResearchBrief(
        brief_id="selfcheck-brief",
        raw_indication="a real-labelled indication",
        resolved_disease_id="MONDO_0000001",
        resolved_disease_label="a real-labelled disease",
        modality="de novo minibinder",
        status="confirmed",
        provenance=prov,
    )
    target = Target(
        target_id="ENSG00000000001",
        gene_symbol="SELFCHECK1",
        stable_gene_id="ENSG00000000001",
        protein_accession="P00000",
        provenance=prov,
    )
    ev = Evidence(
        evidence_id="EV-real-1",
        subject_id="ENSG00000000001",
        predicate="associated_with",
        object_id="MONDO_0000001",
        evidence_type="genetic_association",
        direction="supports",
        source="open_targets",
        source_url="https://platform.opentargets.org/",
        origin="real",
    )
    assess = Assessment(
        target_id="ENSG00000000001",
        brief_id="selfcheck-brief",
        criteria=[
            CriterionRating(criterion="clinical_causal_support", rating="moderate",
                            evidence_ids=["EV-real-1"], rationale="real-labelled"),
            CriterionRating(criterion="cellular_context", rating="unknown",
                            missingness="no Census summary retrieved"),
        ],
        design_feasibility=DesignFeasibility(accessible_to_modality=None),
        rubric_version="rubric-2026.09.22-v1",
    )
    store.upsert(key, brief=brief, targets=[target], evidence=[ev], assessments=[assess],
                 replace_collections=True)
    html = client.get(f"/run/{key}").text
    assert "origin-real" in html, "a real record did not render a real badge"
    assert "Fixture run" not in html, "a real run was banner-labelled as a fixture run"
    j = client.get(f"/run/{key}/state.json").json()
    assert j["is_fixture_run"] is False
    store.delete(key)


@check("no secret-shaped material reaches a template")
def _c_no_secrets(client: TestClient, store: Any) -> None:
    patterns = [
        re.compile(r"(api[_-]?key|client[_-]?secret|secret[_-]?key|password"
                   r"|bearer\s+[A-Za-z0-9._-]{12,})", re.I),
        re.compile(r"sk-[A-Za-z0-9]{16,}"),
        re.compile(r"AKIA[0-9A-Z]{12,}"),
    ]
    pages = ["/"] + [f"/run/fixture-{k}" for k in SCENARIOS]
    for path in pages:
        html = client.get(path).text
        for pat in patterns:
            m = pat.search(html)
            assert m is None, f"{path}: possible secret-shaped string {m.group(0)!r}"


# ---------------------------------------------------------------------------------


def run_checks(verbose: bool = True) -> tuple[int, int, list[str]]:
    store = reset_store()
    app = create_app(store)
    client = TestClient(app)
    passed, failed, messages = 0, 0, []
    for name, fn in CHECKS:
        try:
            fn(client, store)
        except AssertionError as exc:
            failed += 1
            messages.append(f"FAIL  {name}: {exc}")
            if verbose:
                print(f"FAIL  {name}\n      {exc}")
        except Exception as exc:  # noqa: BLE001
            failed += 1
            messages.append(f"ERROR {name}: {type(exc).__name__}: {exc}")
            if verbose:
                print(f"ERROR {name}\n      {type(exc).__name__}: {exc}")
        else:
            passed += 1
            if verbose:
                print(f"pass  {name}")
    if verbose:
        print(f"\n{passed} passed, {failed} failed, {len(CHECKS)} total")
    return passed, failed, messages


if __name__ == "__main__":
    _, failed, _ = run_checks()
    sys.exit(1 if failed else 0)
