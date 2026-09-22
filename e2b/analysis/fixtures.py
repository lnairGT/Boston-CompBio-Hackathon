"""Fixture candidates. FOR UNIT TESTS ONLY.

Every record returned by this module carries ``origin="fixture"``. These sequences
were written by hand to exercise specific code paths in the triage (a free
cysteine, a hydrophobic patch, a glycosylation sequon, missing metrics, an
unrecognised metric name). They are **not** design outputs, they were not produced
by any model, they have never been expressed, and no triage verdict computed on
them is a biological result.

The triage report marks a batch whose records are all fixtures
(``contains_only_fixtures``) and the markdown renderer prints a banner, so a
fixture batch cannot be mistaken for a design run in the UI.
"""

from __future__ import annotations

from typing import Any

FIXTURE_NOTE = (
    "hand-written stand-in for unit tests; not a design output, not a biological result"
)


def _metric(name: str, value: float | None, model: str = "fixture-model", version: str = "n/a") -> dict[str, Any]:
    return {
        "metric_name": name,
        "value": value,
        "method": "fixture value, not computed and not measured",
        "model_name": model,
        "model_version": version,
        "higher_is_better": None,
        "interpretation": FIXTURE_NOTE,
    }


def fixture_candidates() -> list[dict[str, Any]]:
    """Five ``contracts.Candidate``-shaped dicts covering the triage decision paths."""
    return [
        {
            # clean-ish, well-behaved designed helical bundle shape, no Cys, no sequon
            "candidate_id": "FIX-001-clean",
            "run_id": "fixture-run",
            "design_id": "fixture-design",
            "target_id": None,
            "origin": "fixture",
            "sequence": (
                "SEEELKKLAEELKKRAEELAKKNPSDEEVLKLAEEAAKLAEKNGDPELLKKAEELIKRAQELMKKGN"
            ),
            "metrics": [
                _metric("iptm", 0.86),
                _metric("pae_interface", 3.9),
                _metric("plddt_binder", 91.2),
                _metric("self_consistency_rmsd", 0.9),
            ],
        },
        {
            # odd cysteine count -> free thiol; plus a sequon
            "candidate_id": "FIX-002-free-cys",
            "run_id": "fixture-run",
            "origin": "fixture",
            "sequence": (
                "SEEELKKLCEELKKRAEELAKKNGSDEEVLKLAEEAAKLAEKNGDPELLKKAEELIKRAQELMKKGN"
            ),
            "metrics": [
                _metric("iptm", 0.72),
                _metric("pae_interface", 7.4),
                _metric("plddt_binder", 84.0),
            ],
        },
        {
            # strongly hydrophobic with large uncharged patches
            "candidate_id": "FIX-003-hydrophobic",
            "run_id": "fixture-run",
            "origin": "fixture",
            "sequence": (
                "MAVLLVLAIGAVLLVGAAILVLGVAALVIGGVLLAVGAILVLGAVLLVGAAILVLGVAALVIGG"
            ),
            "metrics": [
                _metric("iptm", 0.55),
                _metric("pae_interface", 12.1),
                _metric("n_clashes", 7.0),
            ],
        },
        {
            # no interface metrics at all -> unknown, must be conservative
            "candidate_id": "FIX-004-no-metrics",
            "run_id": "fixture-run",
            "origin": "fixture",
            "sequence": (
                "SPEELAKRLAEELKKQAEELAKKNPSDEEVLKLAEEAAKLAEKQGDPELLKRAEELIKRAQELAKKQ"
            ),
            "metrics": [],
        },
        {
            # very low confidence + an unrecognised metric name + no sequence
            "candidate_id": "FIX-005-verylow-nosequence",
            "run_id": "fixture-run",
            "origin": "fixture",
            "sequence": None,
            "metrics": [
                _metric("iptm", 0.31),
                _metric("pae_interface", None),
                _metric("some_new_tool_score", 0.44),
            ],
        },
    ]
