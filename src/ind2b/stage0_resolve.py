"""Stage 0 - resolve a free-text indication to a disease ontology identifier.

Why this is a stage of its own: every downstream number is conditioned on the
disease id, so guessing an EFO id silently changes the whole result. This
stage resolves the id from the ontology, records the alternatives it rejected
and the rule it used to choose, and flags the ambiguous cases for a human.

A note on the disease subtree
-----------------------------
Open Targets association scores for a parent disease already aggregate
evidence from its descendant terms (the platform's "indirect" association).
Pulling associations for each descendant as well would therefore count the
same evidence repeatedly. So this stage *records* the subtree for provenance
and for optional subtype-specific views, but the default association scope
stays at the parent term. ``association_scope`` in the output states which
was used.
"""

from __future__ import annotations

import json
import logging
import re
from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__
from .config import STAGE_FILES
from .sources import opentargets as ot

log = logging.getLogger("ind2b.stage0")

SCHEMA_VERSION = 1

_TOKEN_RE = re.compile(r"[a-z0-9]+")


def _tokens(text: str) -> set[str]:
    return set(_TOKEN_RE.findall((text or "").lower()))


def select_disease(
    query: str,
    hits: list[dict],
    *,
    synonyms_of: Callable[[list[str]], dict[str, dict[str, list[str]]]] | None = None,
) -> tuple[dict, str, bool]:
    """Choose one ontology hit for ``query``.

    Returns ``(hit, rule, ambiguous)``. The rule is recorded in the stage
    output so a reader can see *why* a term was chosen, and ``ambiguous`` is
    true when the choice rested on search rank alone - the case where a human
    should confirm or pass an explicit id.

    Precedence: exact name match, then exact *synonym* match, then a hit whose
    name contains every query token, then top search rank.

    The synonym step exists because a disease's everyday name is often a MONDO
    synonym rather than its preferred label, and token matching on labels alone
    then picks the wrong term. "atopic dermatitis" is the case that motivated
    it: MONDO's label for that disease is "atopic eczema", which shares no
    token with "dermatitis", while the rare Mendelian term "immunodeficiency
    11b with atopic dermatitis" contains the query verbatim - so token matching
    silently selected a different disease. "atopic dermatitis" is a
    ``hasExactSynonym`` of the intended term, which settles it on ontology
    evidence rather than on string overlap.

    Only ``hasExactSynonym`` is treated as identity. Broad and narrow synonyms
    name a more general or more specific disease and must not collapse into it.

    ``synonyms_of`` is injected so this function stays pure and offline-testable;
    when omitted the synonym step is skipped and the previous precedence holds.
    """
    if not hits:
        raise ValueError(f"no disease hits for query {query!r}")

    q_norm = (query or "").strip().lower()
    for hit in hits:
        if (hit.get("name") or "").strip().lower() == q_norm:
            return hit, "exact_name_match", False

    if synonyms_of is not None and q_norm:
        ids = [h["id"] for h in hits if h.get("id")]
        try:
            syn = synonyms_of(ids)
        except Exception:  # a synonym lookup failure must not break resolution
            log.warning("synonym lookup failed; falling back to name matching", exc_info=True)
            syn = {}
        for hit in hits:
            exact = syn.get(hit.get("id"), {}).get("hasExactSynonym") or []
            if q_norm in {t.strip().lower() for t in exact}:
                return hit, "exact_synonym_match", False

    q_tokens = _tokens(query)
    if q_tokens:
        for hit in hits:
            if q_tokens <= _tokens(hit.get("name")):
                return hit, "all_query_tokens_in_name", False

    return hits[0], "top_search_rank", len(hits) > 1


def resolve(
    indication: str,
    *,
    efo_id: str | None = None,
    n_candidates: int = 10,
) -> dict[str, Any]:
    """Resolve an indication string (or verify an explicit id) against the ontology.

    Passing ``efo_id`` skips text search and verifies that the id exists,
    which is the reproducible path once a term has been agreed.
    """
    candidates: list[dict] = []
    if efo_id:
        rule, ambiguous = "explicit_efo_id", False
        chosen_id = efo_id
    else:
        if not indication:
            raise ValueError("provide either an indication string or an efo_id")
        candidates = ot.search_disease(indication, size=n_candidates)
        hit, rule, ambiguous = select_disease(
            indication, candidates, synonyms_of=ot.diseases_synonyms
        )
        chosen_id = hit["id"]

    info = ot.disease_info(chosen_id)  # raises if the id is unknown

    descendants = info.get("descendants") or []
    record: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "package_version": __version__,
        "resolved_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": "Open Targets Platform GraphQL (EFO/MONDO ontology)",
        "query": indication,
        "selection_rule": rule,
        "ambiguous": ambiguous,
        "disease": {
            "id": info["id"],
            "name": info.get("name"),
            "description": info.get("description"),
            "therapeutic_areas": info.get("therapeuticAreas") or [],
        },
        "parents": info.get("parents") or [],
        "children": info.get("children") or [],
        "descendants": {"count": len(descendants), "ids": descendants},
        "association_scope": "indirect_parent",
        "candidates_considered": [
            {"id": c.get("id"), "name": c.get("name")} for c in candidates
        ],
    }
    if ambiguous:
        log.warning(
            "indication %r resolved by search rank alone to %s (%s); "
            "pass --efo-id to pin it",
            indication, record["disease"]["id"], record["disease"]["name"],
        )
    return record


def write(record: dict[str, Any], run_dir: Path) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)
    path = run_dir / STAGE_FILES[0]
    path.write_text(json.dumps(record, indent=2))
    return path


def read(run_dir: Path) -> dict[str, Any]:
    path = run_dir / STAGE_FILES[0]
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found - run stage 0 first (`ind2b stage0 --indication ...`)"
        )
    return json.loads(path.read_text())
