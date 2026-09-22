"""Presentation-layer honesty helpers.

Every rule in this module exists because breaking it would make the UI assert
something the underlying record does not support. They are small, pure and unit
testable so Lane 4 can check them directly.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any

# ---------------------------------------------------------------------------------
# Origin labelling
# ---------------------------------------------------------------------------------

ORIGIN_BADGES: dict[str, dict[str, str]] = {
    "real": {
        "label": "real",
        "css": "origin-real",
        "tooltip": "Retrieved or computed in this run, from the named source.",
    },
    "cached_real": {
        "label": "cached",
        "css": "origin-cached",
        "tooltip": "A previously retrieved real record, replayed. Its original time is shown.",
    },
    "fixture": {
        "label": "fixture",
        "css": "origin-fixture",
        "tooltip": (
            "Hand-written stand-in used to exercise the interface. NOT a biological "
            "result and never admissible as evidence."
        ),
    },
}

UNKNOWN_ORIGIN = {
    "label": "origin unrecorded",
    "css": "origin-unknown",
    "tooltip": "This record carries no origin label. Treat it as unverified.",
}


def origin_badge(origin: str | None) -> dict[str, str]:
    """Return the badge descriptor for an ``origin`` value.

    An unrecognised or missing origin is surfaced as ``origin unrecorded`` rather than
    silently defaulting to ``real``.
    """
    if origin is None:
        return dict(UNKNOWN_ORIGIN)
    return dict(ORIGIN_BADGES.get(origin, UNKNOWN_ORIGIN))


def is_fixture(obj: Any) -> bool:
    return getattr(obj, "origin", None) == "fixture"


# ---------------------------------------------------------------------------------
# Numbers: unknown is never zero
# ---------------------------------------------------------------------------------

NOT_MEASURED = "not measured"


@dataclass(frozen=True)
class Num:
    """A number that knows whether it is actually known."""

    text: str
    known: bool
    raw: float | int | None = None

    def __str__(self) -> str:  # pragma: no cover - convenience only
        return self.text


def num(
    value: float | int | None,
    *,
    precision: int = 2,
    unit: str | None = None,
    missing: str = NOT_MEASURED,
    percent: bool = False,
) -> Num:
    """Format a numeric for display.

    ``None`` renders as *not measured*. ``0`` renders as ``0`` — a measured zero is a
    real result and is never conflated with a missing one.
    """
    if value is None:
        return Num(text=missing, known=False, raw=None)
    if percent:
        body = f"{value * 100:.{max(precision - 2, 0)}f}%"
    else:
        body = f"{value:.{precision}f}" if isinstance(value, float) else f"{value:,}"
    if unit:
        body = f"{body} {unit}"
    return Num(text=body, known=True, raw=value)


def count(value: int | None, *, missing: str = "not counted") -> Num:
    return num(value, missing=missing)


# ---------------------------------------------------------------------------------
# Times
# ---------------------------------------------------------------------------------


def _aware(dt: datetime) -> datetime:
    return dt if dt.tzinfo else dt.replace(tzinfo=timezone.utc)


def ts(value: datetime | None, *, missing: str = "no timestamp") -> str:
    if value is None:
        return missing
    return _aware(value).strftime("%Y-%m-%d %H:%M:%S UTC")


def age_text(value: datetime | None, *, now: datetime | None = None) -> str:
    """Human age of a timestamp, e.g. ``3 min ago``. Empty string when unknown."""
    if value is None:
        return ""
    now = now or datetime.now(timezone.utc)
    delta = (now - _aware(value)).total_seconds()
    if delta < 0:
        return "in the future"
    for limit, div, word in ((90, 1, "s"), (5400, 60, "min"), (172800, 3600, "h")):
        if delta < limit:
            return f"{int(delta // div)} {word} ago"
    return f"{int(delta // 86400)} d ago"


def freshness(obj: Any, *, now: datetime | None = None) -> dict[str, Any]:
    """Freshness line for any record that has an origin plus a time.

    For ``cached_real`` the caller MUST have populated the original retrieval or run
    time on the record; this helper renders that original time and says plainly that
    the record is a replay. It never presents a replay as a fresh completion.
    """
    origin = getattr(obj, "origin", None)
    prov = getattr(obj, "provenance", None)
    when = (
        getattr(obj, "finished_at", None)
        or (getattr(prov, "retrieved_at", None) if prov is not None else None)
        or getattr(obj, "retrieved_at", None)
        or getattr(obj, "updated_at", None)
    )
    out: dict[str, Any] = {
        "badge": origin_badge(origin),
        "when": ts(when),
        "age": age_text(when, now=now),
        "replay": origin == "cached_real",
    }
    if origin == "cached_real":
        out["note"] = (
            f"Replay of a run that originally completed at {ts(when)}"
            f"{' (' + out['age'] + ')' if out['age'] else ''}. Not a new model run."
        )
    elif origin == "fixture":
        out["note"] = "Fixture record. Interface stand-in only — not a biological result."
    else:
        out["note"] = None
    return out


# ---------------------------------------------------------------------------------
# Vocabulary guards
# ---------------------------------------------------------------------------------

#: Substrings that would misdescribe a de novo minibinder as an immunoglobulin.
ANTIBODY_WORDS = ("antibody", "antibodies", "mab", "monoclonal", "igg", "nanobody")

#: Measured binding quantities. A computational score may never be shown under these.
MEASURED_AFFINITY_WORDS = ("kd", "ki", "ic50", "ec50", "affinity", "koff", "kon", "nanomolar")

_WORD = re.compile(r"[a-z0-9]+")


def _words(text: str) -> set[str]:
    return set(_WORD.findall(text.lower()))


def modality_phrase(modality: str | None) -> str:
    """Render a modality for display without ever upgrading it to an antibody."""
    if not modality:
        return "modality not specified"
    low = modality.strip().lower()
    if low in ("de novo minibinder", "minibinder"):
        return "de novo minibinder (small designed protein — not an antibody)"
    return modality.strip()


def check_no_antibody_claim(text: str, modality: str | None) -> list[str]:
    """Return warnings if ``text`` calls a non-antibody modality an antibody."""
    if not text or not modality:
        return []
    if "antibody" in modality.lower() or "nanobody" in modality.lower():
        return []
    hits = sorted(_words(text) & set(ANTIBODY_WORDS))
    if hits:
        return [
            f"Text uses antibody vocabulary ({', '.join(hits)}) for modality "
            f"'{modality}'. A designed minibinder is not an antibody."
        ]
    return []


def metric_caption(metric: Any) -> dict[str, Any]:
    """View dict for a :class:`e2b.contracts.Metric`.

    Adds an explicit disclaimer whenever the metric name is close to a measured
    binding quantity, and always names the producing model. ``value=None`` renders as
    *not computed*, never 0.
    """
    name = getattr(metric, "metric_name", "unnamed metric")
    value = getattr(metric, "value", None)
    model = getattr(metric, "model_name", None) or "unnamed model"
    version = getattr(metric, "model_version", None)
    hib = getattr(metric, "higher_is_better", None)
    looks_measured = bool(_words(name) & set(MEASURED_AFFINITY_WORDS))
    return {
        "name": name,
        "value": num(value, missing="not computed"),
        "method": getattr(metric, "method", None),
        "model": f"{model}{' ' + version if version else ''}",
        "direction": (
            "higher is better" if hib is True else "lower is better" if hib is False else
            "direction of merit unrecorded"
        ),
        "interpretation": getattr(metric, "interpretation", None),
        "is_model_score": True,
        "disclaimer": (
            "Computational model score. This is NOT a measured binding affinity and has "
            "no experimental Kd behind it."
            if looks_measured
            else "Computational model score from the named model. Not an experimental measurement."
        ),
    }


# ---------------------------------------------------------------------------------
# Source failure is not evidence of absence
# ---------------------------------------------------------------------------------

FAILURE_KINDS: dict[str, str] = {
    "timeout": "The source did not answer in time. This is a retrieval failure, not evidence "
    "that the data is absent.",
    "http_error": "The source returned an error. Nothing can be concluded about the data.",
    "unreachable": "The source could not be reached. Nothing can be concluded about the data.",
    "not_configured": "No adapter or credential is attached for this source in this deployment.",
    "empty_result": "The source answered and returned no matching records. This IS a real "
    "negative for the query as issued, and the query is shown.",
    "unsupported": "The source cannot answer this query shape.",
}


@dataclass
class SourceFailure:
    """UI-local note that a source did not deliver. Not a contracts record.

    ``kind='empty_result'`` is the only kind that carries information about the world,
    and even then only about the query as issued.
    """

    source: str
    kind: str
    detail: str | None = None
    query: str | None = None
    at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))

    def view(self) -> dict[str, Any]:
        return {
            "source": self.source,
            "kind": self.kind,
            "meaning": FAILURE_KINDS.get(self.kind, "Unrecognised failure kind."),
            "detail": self.detail,
            "query": self.query,
            "at": ts(self.at),
            "is_absence_evidence": self.kind == "empty_result",
        }


__all__ = [
    "ORIGIN_BADGES",
    "Num",
    "SourceFailure",
    "age_text",
    "check_no_antibody_claim",
    "count",
    "freshness",
    "is_fixture",
    "metric_caption",
    "modality_phrase",
    "num",
    "origin_badge",
    "ts",
]
