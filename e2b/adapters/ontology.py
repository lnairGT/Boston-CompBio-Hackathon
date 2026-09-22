"""Disease identity across sources, grounded in the MONDO ontology.

Different sources name the same disease differently, and matching on strings loses real
data. The case that motivated this module: Open Targets resolves "atopic dermatitis" to
**MONDO_0004980, whose label is "atopic eczema"**. A Census search for the literal string
"atopic dermatitis" finds nothing and concludes the disease is absent — when it is present
under the ontology's preferred label. "Atopic dermatitis" is merely a *synonym*.

The same lookup shows `EFO:0000274` is a **cross-reference** of MONDO_0004980, not a
different disease. Open Targets migrated its identifiers between releases; the ontology
records that the two refer to the same thing.

So disease identity is resolved on **ontology terms**, never on labels, and every match
records *how* it was reached:

``exact``
    Same term. Safe.
``cross_reference``
    A different ontology's identifier for the same disease (EFO ↔ MONDO). Safe.
``synonym_label``
    Matched a label against a recorded synonym. Usable, but weaker than an ID match, and
    recorded as such because synonym lists are curated and incomplete.
``descendant``
    The vocabulary term is **narrower** than what was asked for — a sub-type. Using it
    answers a *different, more specific* question.
``ancestor``
    The vocabulary term is **broader** — a parent disease. Using it answers a *more
    general* question and will include unrelated patients.

The last two are returned, never silently accepted. Substituting a parent or child disease
without saying so is how a cohort quietly stops matching the brief.

Source: EBI Ontology Lookup Service (OLS4), which serves MONDO directly.
"""

from __future__ import annotations

import time
from typing import Any, Iterable, Literal
from urllib.parse import quote

import requests

ADAPTER_VERSION = "ontology/1.0.0"

OLS4 = "https://www.ebi.ac.uk/ols4/api"
OBO_IRI = "http://purl.obolibrary.org/obo/{term}"
TIMEOUT_S = 45

MatchType = Literal["exact", "cross_reference", "synonym_label", "descendant", "ancestor", "none"]

#: Matches that mean "the same disease". Anything else changes the question being asked.
SAME_DISEASE: frozenset[str] = frozenset({"exact", "cross_reference", "synonym_label"})


class OntologyError(RuntimeError):
    def __init__(self, message: str, *, transient: bool = False):
        super().__init__(message)
        self.transient = transient


def _curie_to_iri(curie: str) -> str:
    """``MONDO:0004980`` / ``MONDO_0004980`` -> the OBO PURL, double-encoded for OLS4."""
    term = curie.strip().replace(":", "_")
    return quote(quote(OBO_IRI.format(term=term), safe=""), safe="")


def _normalise(curie: str) -> str:
    """Canonical CURIE form: ``MONDO:0004980``."""
    c = curie.strip().replace("_", ":")
    parts = c.split(":")
    return f"{parts[0].upper()}:{parts[-1]}" if len(parts) >= 2 else c


def _get(url: str, params: dict | None = None) -> requests.Response:
    try:
        return requests.get(url, params=params or {}, timeout=TIMEOUT_S,
                            headers={"Accept": "application/json"})
    except requests.RequestException as exc:
        raise OntologyError(f"OLS4 unreachable: {exc}", transient=True) from exc


class DiseaseOntology:
    """Read-only MONDO lookups with a small in-process cache."""

    def __init__(self, ontology: str = "mondo"):
        self.ontology = ontology
        self._terms: dict[str, dict[str, Any]] = {}
        self._desc: dict[str, set[str]] = {}
        self._anc: dict[str, set[str]] = {}

    # -- single term ----------------------------------------------------------------

    def term(self, curie: str) -> dict[str, Any]:
        key = _normalise(curie)
        if key in self._terms:
            return self._terms[key]
        r = _get(f"{OLS4}/ontologies/{self.ontology}/terms/{_curie_to_iri(key)}")
        if r.status_code == 404:
            raise OntologyError(f"{key} is not a term in {self.ontology.upper()}.")
        if r.status_code != 200:
            raise OntologyError(f"OLS4 HTTP {r.status_code} for {key}", transient=r.status_code >= 500)
        d = r.json()
        xrefs = []
        for x in d.get("obo_xref") or []:
            db, xid = x.get("database"), x.get("id")
            if db and xid:
                xrefs.append(_normalise(f"{db}:{xid}"))
        out = {
            "curie": key,
            "label": d.get("label"),
            "definition": (d.get("description") or [None])[0],
            "synonyms": sorted({s for s in (d.get("synonyms") or []) if s}),
            "cross_references": sorted(set(xrefs)),
            "is_obsolete": bool(d.get("is_obsolete")),
            "source": "EBI OLS4 / MONDO",
            "retrieved_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        }
        self._terms[key] = out
        return out

    def _related(self, curie: str, relation: str, cache: dict, cap: int = 400) -> set[str]:
        key = _normalise(curie)
        if key in cache:
            return cache[key]
        found: set[str] = set()
        page, size = 0, 100
        while len(found) < cap:
            r = _get(f"{OLS4}/ontologies/{self.ontology}/terms/{_curie_to_iri(key)}/{relation}",
                     {"size": size, "page": page})
            if r.status_code != 200:
                break
            body = r.json()
            terms = (body.get("_embedded") or {}).get("terms") or []
            for t in terms:
                if t.get("obo_id"):
                    found.add(_normalise(t["obo_id"]))
            pg = body.get("page") or {}
            page += 1
            if page >= (pg.get("totalPages") or 1) or not terms:
                break
        cache[key] = found
        return found

    def descendants(self, curie: str) -> set[str]:
        """Narrower terms: sub-types of this disease."""
        return self._related(curie, "hierarchicalDescendants", self._desc)

    def ancestors(self, curie: str) -> set[str]:
        """Broader terms: diseases this one is a kind of."""
        return self._related(curie, "hierarchicalAncestors", self._anc)

    # -- mapping --------------------------------------------------------------------

    def equivalent_ids(self, curie: str) -> set[str]:
        """Every identifier that denotes THIS disease: itself plus its cross-references."""
        t = self.term(curie)
        return {t["curie"], *t["cross_references"]}

    def match(self, query_curie: str, candidate_curie: str | None,
              candidate_label: str | None = None) -> dict[str, Any]:
        """Classify how a vocabulary entry relates to the disease actually asked about."""
        q = _normalise(query_curie)
        qterm = self.term(q)
        cand = _normalise(candidate_curie) if candidate_curie else None

        match_type: MatchType = "none"
        basis = ""

        if cand:
            if cand == q:
                match_type, basis = "exact", f"Same ontology term {q}."
            elif cand in self.equivalent_ids(q):
                match_type, basis = "cross_reference", (
                    f"{cand} is a recorded cross-reference of {q} ({qterm['label']}); "
                    "they denote the same disease.")
            elif cand in self.descendants(q):
                match_type, basis = "descendant", (
                    f"{cand} is a SUB-TYPE of {q} ({qterm['label']}). Using it narrows the "
                    "question to a more specific disease.")
            elif cand in self.ancestors(q):
                match_type, basis = "ancestor", (
                    f"{cand} is a PARENT of {q} ({qterm['label']}). Using it broadens the "
                    "question and will include patients who do not have the queried disease.")

        if match_type == "none" and candidate_label:
            syns = {s.strip().lower() for s in qterm["synonyms"]}
            syns.add((qterm["label"] or "").strip().lower())
            if candidate_label.strip().lower() in syns:
                match_type, basis = "synonym_label", (
                    f"'{candidate_label}' is a recorded synonym of {q} ({qterm['label']}). "
                    "Matched on label, which is weaker than an identifier match.")

        return {
            "query": q,
            "query_label": qterm["label"],
            "candidate": cand,
            "candidate_label": candidate_label,
            "match_type": match_type,
            "same_disease": match_type in SAME_DISEASE,
            "basis": basis or "No ontology relationship found between the two terms.",
        }

    def map_to_vocabulary(
        self,
        query_curie: str,
        vocabulary: Iterable[tuple[str | None, str | None]],
        *,
        include_descendants: bool = True,
        include_ancestors: bool = False,
    ) -> dict[str, Any]:
        """Find every entry of a source's disease vocabulary that relates to this disease.

        ``vocabulary`` is ``(ontology_term_id, label)`` pairs — for CELLxGENE, the
        ``disease_ontology_term_id`` / ``disease`` fields, which are already MONDO.

        ``include_descendants`` defaults to True because a sub-type cohort is usually
        still informative, but every such entry is returned tagged ``descendant`` so the
        caller must decide deliberately. ``include_ancestors`` defaults to False: a parent
        disease is a different, broader population.
        """
        exact, broader, narrower, unrelated = [], [], [], 0
        for term_id, label in vocabulary:
            m = self.match(query_curie, term_id, label)
            if m["same_disease"]:
                exact.append(m)
            elif m["match_type"] == "descendant":
                narrower.append(m)
            elif m["match_type"] == "ancestor":
                broader.append(m)
            else:
                unrelated += 1

        usable = list(exact)
        if include_descendants:
            usable += narrower
        if include_ancestors:
            usable += broader

        return {
            "query": _normalise(query_curie),
            "query_label": self.term(query_curie)["label"],
            "same_disease": exact,
            "narrower_subtypes": narrower,
            "broader_parents": broader,
            "unrelated_count": unrelated,
            "usable": usable,
            "usable_term_ids": sorted({m["candidate"] for m in usable if m["candidate"]}),
            "coverage": (
                "no_coverage" if not usable
                else "exact_only" if not narrower and not broader
                else "includes_related_terms"
            ),
            "caveat": (
                None if not (narrower or broader) else
                "This mapping includes terms that are not the queried disease: "
                + ", ".join(
                    f"{m['candidate']} ({m['candidate_label']}, {m['match_type']})"
                    for m in (narrower + broader)[:6]
                )
                + ". Report which terms contributed to any cohort built from this mapping."
            ),
            "source": "EBI OLS4 / MONDO",
            "adapter_version": ADAPTER_VERSION,
        }


_DEFAULT = DiseaseOntology()


def resolve_disease_synonyms(curie: str) -> dict[str, Any]:
    """Convenience: every label and identifier that denotes this disease."""
    t = _DEFAULT.term(curie)
    return {
        "curie": t["curie"], "label": t["label"],
        "synonyms": t["synonyms"], "cross_references": t["cross_references"],
        "search_strings": sorted({t["label"], *t["synonyms"]} - {None}),
    }


def map_to_vocabulary(query_curie: str, vocabulary: Iterable[tuple[str | None, str | None]],
                      **kw: Any) -> dict[str, Any]:
    return _DEFAULT.map_to_vocabulary(query_curie, vocabulary, **kw)
