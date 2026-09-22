"""Lane 5 — experience / UI.

A server-rendered single page over the evidence-to-binder workflow. Five sections in
fixed order (Brief, Compare, Inspect, Design, Results) with progressive detail.

The UI layer holds **no scientific logic**. It reads :mod:`e2b.contracts` records that
the spine puts into :class:`e2b.ui.state.UiStore` and renders them; where a record is
absent it renders an explicit *no data* state rather than inventing one.

Public surface for the spine::

    from e2b.ui import create_app, get_store
    store = get_store()
    store.upsert(run_key, brief=brief, agent_run=agent_run, ...)
    app = create_app()                      # uvicorn e2b.ui.app:app

Honesty rules enforced in this layer (Lane 4 tests them):

* every record with an ``origin`` renders a visible badge (real / cached / fixture);
* ``None`` numerics render as "not measured", never as ``0``;
* a ``cached_real`` record renders its ORIGINAL run time, never "just completed";
* no synthetic percentage progress — stage name plus timestamp only;
* a minibinder is never called an antibody, and a model score is never labelled an
  affinity or a Kd.
"""

from e2b.ui.app import create_app  # noqa: F401
from e2b.ui.state import UiStore, get_store  # noqa: F401

__all__ = ["create_app", "UiStore", "get_store"]
