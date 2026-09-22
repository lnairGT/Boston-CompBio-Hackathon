"""Accept either a pydantic record from ``e2b.contracts`` or a plain dict.

Lane 7 must not block on the exact shape of another lane's objects, and it must not
crash the demo because a field arrived as a dict instead of a model. Everything that
enters this module goes through ``as_dict`` first, and every field read goes through
``get`` so a missing field is ``None`` rather than an exception.
"""

from __future__ import annotations

from typing import Any


def as_dict(obj: Any) -> dict[str, Any]:
    """Coerce a pydantic v2 model, a dataclass-like object or a dict into a dict."""
    if obj is None:
        return {}
    if isinstance(obj, dict):
        return dict(obj)
    dump = getattr(obj, "model_dump", None)
    if callable(dump):
        try:
            return dump(mode="python")
        except TypeError:
            return dump()
    d = getattr(obj, "__dict__", None)
    if isinstance(d, dict):
        return dict(d)
    raise TypeError(f"Cannot read {type(obj).__name__} as a record.")


def get(d: dict[str, Any], *names: str, default: Any = None) -> Any:
    """First present, non-None value among ``names``."""
    for n in names:
        if n in d and d[n] is not None:
            return d[n]
    return default


def text_blob(*parts: Any) -> str:
    """Lower-cased concatenation of free-text inputs, for keyword derivation."""
    out: list[str] = []
    for p in parts:
        if p is None:
            continue
        if isinstance(p, str):
            out.append(p)
        elif isinstance(p, (list, tuple, set)):
            out.extend(str(x) for x in p if x is not None)
        elif isinstance(p, dict):
            out.extend(str(x) for x in p.values() if x is not None)
        else:
            out.append(str(p))
    return " ".join(out).lower()
