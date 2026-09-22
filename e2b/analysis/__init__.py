"""Lane 6 -- downstream analysis: deterministic candidate triage.

Public surface the spine should call:

    from e2b.analysis import (
        TRIAGE_VERSION,
        triage_candidates,          # batch entry point
        triage_candidate,           # one candidate
        compute_developability,     # sequence indicators only
        interpret_interface_metrics,
        paralog_cross_reactivity_risk,
        render_markdown,
    )

Everything returns plain JSON-able dicts. Nothing in this package is a
measurement, and nothing here asserts that a candidate binds.
"""

from .developability import (
    METHOD_VERSION as DEVELOPABILITY_VERSION,
    as_contract_metric_dicts,
    compute_developability,
)
from .interface import (
    INTERFACE_VERSION,
    METRIC_REGISTRY,
    claim_language,
    interpret_interface_metrics,
)
from .paralogs import PARALOG_VERSION, paralog_cross_reactivity_risk, retrieve_family
from .triage import (
    RANKING_KEY_DESCRIPTION,
    TRIAGE_THRESHOLDS,
    TRIAGE_VERSION,
    render_markdown,
    triage_candidate,
    triage_candidates,
)

__all__ = [
    "TRIAGE_VERSION",
    "TRIAGE_THRESHOLDS",
    "RANKING_KEY_DESCRIPTION",
    "DEVELOPABILITY_VERSION",
    "INTERFACE_VERSION",
    "PARALOG_VERSION",
    "METRIC_REGISTRY",
    "triage_candidates",
    "triage_candidate",
    "compute_developability",
    "as_contract_metric_dicts",
    "interpret_interface_metrics",
    "claim_language",
    "paralog_cross_reactivity_risk",
    "retrieve_family",
    "render_markdown",
]
