"""CELLxGENE Discover client - is there single-cell data for this disease, and in what?

The evidence layers in stage 1 answer "is this target associated with the disease".
This module answers a different question that none of them cover: **in this disease's
own tissue, which cell populations are present, and was the diseased state actually
sampled** - or only healthy tissue.

Two things this module is careful about, because both are easy to get wrong and both
change what a result means.

*Discover is not the Census.* CELLxGENE publishes a dataset index (Discover, queried
here) and a periodically-built expression corpus (Census). They are not the same
corpus and they do not carry the same terms: ``MONDO:0004980`` (atopic eczema) is
present in Discover and absent from every Census build checked in September 2026,
which the Census lags by roughly six months. A coverage answer is therefore only
meaningful with the corpus named, and every record here says ``corpus="discover"``.

*Match on the ontology term, never on the label.* Disease labels differ between
sources for the same disease - Open Targets resolves the input "atopic dermatitis" to
``MONDO_0004980``, whose CELLxGENE label is "atopic eczema". A label search finds
nothing and reports absence. Stage 0 already produces the ontology id; this module
takes that id and matches ``ontology_term_id`` directly, so the label never enters the
comparison.

What this module does **not** establish: presence of a dataset is not evidence of
differential expression, RNA is not surface protein abundance, and a disease absent
from Discover is absent *from Discover*, not from biology.
"""

from __future__ import annotations

import logging
from collections import Counter
from typing import Any

import requests

from ..config import CELLXGENE_CURATION
from ..http import request_json

log = logging.getLogger("ind2b.cellxgene")

CORPUS = "discover"

# A dataset with both the disease term and a normal/healthy term is worth far more than
# two datasets with one each: the disease-vs-reference contrast is then within-study, so
# batch and protocol are shared rather than confounded with the comparison.
# PATO:0000461 is the CELLxGENE convention for "normal".
NORMAL_TERM_IDS = frozenset({"PATO:0000461"})


def _datasets() -> list[dict] | None:
    """The full public dataset index. ``None`` distinguishes failure from an empty index."""
    try:
        data = request_json(f"{CELLXGENE_CURATION}/datasets")
    except (requests.HTTPError, RuntimeError) as exc:
        log.warning("cellxgene dataset index unavailable: %s", exc)
        return None
    return data if isinstance(data, list) else None


def _norm(term_id: str) -> str:
    """MONDO_0004980 and MONDO:0004980 are the same term written two ways."""
    return (term_id or "").replace("_", ":").strip()


def coverage(disease_term_id: str) -> dict[str, Any]:
    """Which Discover datasets carry this disease term?

    Returns a record whose ``status`` is one of:

    ``covered``
        At least one dataset carries the term.
    ``not_in_corpus``
        The index was retrieved and the term is not in it. A real, reportable finding.
    ``source_unavailable``
        The index could not be retrieved. **Not** a finding about coverage - the
        question was not answered, and a caller must not treat it as absence.
    """
    wanted = _norm(disease_term_id)
    index = _datasets()
    if index is None:
        return {
            "status": "source_unavailable",
            "corpus": CORPUS,
            "disease_term_id": wanted,
            "note": "Dataset index could not be retrieved; coverage is unknown, not absent.",
        }

    hits = []
    for ds in index:
        terms = {_norm((d or {}).get("ontology_term_id")) for d in (ds.get("disease") or [])}
        if wanted not in terms:
            continue
        # Does this dataset ALSO carry a normal state? If so the contrast is within-study.
        has_normal = bool(terms & NORMAL_TERM_IDS)
        hits.append(
            {
                "dataset_id": ds.get("dataset_id"),
                "dataset_version_id": ds.get("dataset_version_id"),
                "title": ds.get("title"),
                "collection_id": ds.get("collection_id"),
                "cell_count": ds.get("cell_count"),
                "tissues": sorted({(t or {}).get("label") for t in (ds.get("tissue") or []) if t}),
                "assays": sorted({(a or {}).get("label") for a in (ds.get("assay") or []) if a}),
                "disease_terms": sorted(t for t in terms if t),
                "paired_normal_in_same_dataset": has_normal,
                "n_cell_types": len(ds.get("cell_type") or []),
            }
        )

    if not hits:
        return {
            "status": "not_in_corpus",
            "corpus": CORPUS,
            "disease_term_id": wanted,
            "n_datasets_scanned": len(index),
            "note": (
                f"{wanted} is not present in CELLxGENE {CORPUS}. This is absence from this "
                "corpus at this time, not evidence that the disease lacks single-cell data "
                "elsewhere, and not evidence about any target's expression."
            ),
        }

    hits.sort(key=lambda h: -(h["cell_count"] or 0))
    tissues = Counter(t for h in hits for t in h["tissues"])
    return {
        "status": "covered",
        "corpus": CORPUS,
        "disease_term_id": wanted,
        "n_datasets": len(hits),
        "total_cells": sum(h["cell_count"] or 0 for h in hits),
        "tissues": dict(tissues.most_common()),
        "n_with_paired_normal": sum(1 for h in hits if h["paired_normal_in_same_dataset"]),
        "datasets": hits,
        "caveats": [
            "Dataset presence is coverage, not differential expression.",
            "RNA abundance is not surface protein abundance.",
            f"Corpus is CELLxGENE {CORPUS}; the Census is a separate corpus and lags it.",
        ],
    }


def best_contrast_dataset(disease_term_id: str) -> dict[str, Any] | None:
    """The dataset best suited to a disease-vs-reference comparison, or ``None``.

    Preference order: a dataset carrying **both** the disease term and a normal term
    (within-study contrast, shared batch and protocol) over a larger dataset carrying
    only the disease term. Cell count breaks ties within each group. Returns ``None``
    when the disease is not covered or the index was unreachable - the caller must
    distinguish those two by calling :func:`coverage`.
    """
    cov = coverage(disease_term_id)
    if cov["status"] != "covered":
        return None
    ranked = sorted(
        cov["datasets"],
        key=lambda h: (not h["paired_normal_in_same_dataset"], -(h["cell_count"] or 0)),
    )
    best = dict(ranked[0])
    best["contrast"] = (
        "within_study_disease_vs_normal"
        if best["paired_normal_in_same_dataset"]
        else "disease_only_reference_must_come_from_another_study"
    )
    if not best["paired_normal_in_same_dataset"]:
        best["caveat"] = (
            "No normal state in this dataset, so any disease-vs-reference comparison is "
            "between studies and batch is not separable from biology."
        )
    return best