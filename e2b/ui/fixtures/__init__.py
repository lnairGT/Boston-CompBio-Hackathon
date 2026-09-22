"""Labelled UI fixtures.

Every record built here carries ``origin="fixture"``. Fixtures exist for exactly one
reason: to make all five sections and all four job states walkable before the spine is
wired up. They are **not** biological results, they are never consulted by the
pipeline, the agent loop or the ranking code, and the UI badges them on every card.

No fixture carries a DOI or a literature citation — a fabricated citation would be a
manufactured source. Where a fixture needs a locator it uses a datasource identifier
and a generic platform URL, and every fixture evidence item states in its
``limitations`` that it was not retrieved.
"""

from e2b.ui.fixtures.scenarios import (  # noqa: F401
    SCENARIOS,
    available_scenarios,
    load_scenario,
)

__all__ = ["SCENARIOS", "available_scenarios", "load_scenario"]
