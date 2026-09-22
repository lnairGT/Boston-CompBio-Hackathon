"""The evidence graph: connected, source-backed records.

The spec requires the graph as *connected source-backed records*; a graph database and an
interactive network visualisation are explicitly optional. So this is an in-memory
adjacency structure that serialises to JSON, and the demo renders it as a readable
evidence trail rather than a force-directed picture.

Two invariants the rest of the system leans on:

* **Every scientific relation references evidence.** ``add_relation`` refuses an edge whose
  ``evidence_ids`` is empty unless the edge is explicitly declared structural (the
  brief-to-disease link, a run-to-candidate link -- bookkeeping, not biology).
* **Provenance classes stay distinguishable.** Each edge carries ``claim_class``, one of
  ``asserted_literature``, ``observed_metadata``, ``computed_measurement`` or
  ``model_prediction``. Collapsing these is how a model score ends up being read as a
  measurement, so they are kept apart at the data layer, not just in the UI copy.

Population nodes deliberately retain dataset and condition context: a "CD4-positive T cell"
in lesional atopic-dermatitis skin is not the same node as one in a healthy reference, and
merging them would produce a disease-independent entity that supports claims it should not.
"""

from __future__ import annotations

import json
from collections import defaultdict
from typing import Any, Iterable, Literal

from ..contracts import Evidence

NodeKind = Literal[
    "Indication", "Mechanism", "Target", "Gene", "Protein", "Population", "Dataset",
    "Paper", "Structure", "DesignRun", "Candidate", "Drug", "ClinicalEvidence",
    "AssociationEvidence", "CellSlice", "EmbeddingSlice", "Brief",
]

ClaimClass = Literal[
    "asserted_literature", "observed_metadata", "computed_measurement", "model_prediction",
    "structural",
]

# Relations that are bookkeeping rather than scientific claims, and so may carry no evidence.
STRUCTURAL_PREDICATES = frozenset({
    "brief_resolves_to", "run_targets_protein", "candidate_generated_by", "structure_represents",
    "slice_drawn_from", "embedding_of", "assessed_under_rubric", "design_run_submitted_for",
})


class GraphError(ValueError):
    pass


class EvidenceGraph:
    def __init__(self) -> None:
        self.nodes: dict[str, dict[str, Any]] = {}
        self.edges: list[dict[str, Any]] = []
        self.evidence: dict[str, Evidence] = {}
        self._adj: dict[str, list[int]] = defaultdict(list)

    # -- nodes ----------------------------------------------------------------------

    def add_node(self, node_id: str, kind: NodeKind, label: str, **attrs: Any) -> str:
        existing = self.nodes.get(node_id)
        if existing:
            if existing["kind"] != kind:
                raise GraphError(
                    f"Node '{node_id}' already exists as {existing['kind']}, cannot redeclare as {kind}."
                )
            existing["attrs"].update({k: v for k, v in attrs.items() if v is not None})
            return node_id
        self.nodes[node_id] = {"id": node_id, "kind": kind, "label": label, "attrs": dict(attrs)}
        return node_id

    def population_node_id(
        self, label: str, dataset_id: str, condition: str | None, census_release: str
    ) -> str:
        """Population identity is (label, dataset, condition, release) -- never the label alone."""
        cond = (condition or "unspecified").replace(" ", "_")
        return f"pop:{census_release}:{dataset_id}:{cond}:{label.replace(' ', '_')}"

    # -- evidence and edges ---------------------------------------------------------

    def add_evidence(self, ev: Evidence) -> str:
        self.evidence[ev.evidence_id] = ev
        self.add_node(
            ev.evidence_id,
            "AssociationEvidence" if ev.evidence_type != "literature" else "Paper",
            f"{ev.source}:{ev.evidence_type}",
            source_url=ev.source_url,
            doi=ev.doi,
            locator=ev.source_locator,
            disease_scope=ev.disease_scope,
            origin=ev.origin,
        )
        return ev.evidence_id

    def add_relation(
        self,
        subject: str,
        predicate: str,
        obj: str,
        *,
        evidence_ids: Iterable[str] = (),
        claim_class: ClaimClass = "observed_metadata",
        direction: str | None = None,
        disease_scope: str | None = None,
        **attrs: Any,
    ) -> dict[str, Any]:
        ev_ids = list(evidence_ids)
        if not ev_ids and predicate not in STRUCTURAL_PREDICATES:
            raise GraphError(
                f"Relation ({subject})-[{predicate}]->({obj}) carries no evidence_ids. Every "
                f"scientific relation must reference evidence; declare it in STRUCTURAL_PREDICATES "
                f"if it is bookkeeping."
            )
        unknown = [e for e in ev_ids if e not in self.evidence]
        if unknown:
            raise GraphError(f"Relation references unregistered evidence: {unknown}")
        for endpoint in (subject, obj):
            if endpoint not in self.nodes:
                raise GraphError(f"Relation endpoint '{endpoint}' is not a registered node.")

        edge = {
            "subject": subject,
            "predicate": predicate,
            "object": obj,
            "evidence_ids": ev_ids,
            "claim_class": claim_class,
            "direction": direction,
            "disease_scope": disease_scope,
            "attrs": dict(attrs),
        }
        self.edges.append(edge)
        idx = len(self.edges) - 1
        self._adj[subject].append(idx)
        self._adj[obj].append(idx)
        return edge

    def add_similarity(
        self,
        subject: str,
        obj: str,
        *,
        embedding_model: str,
        distance: float,
        metric: str,
        search_universe: str,
        embedding_slice_id: str,
        evidence_ids: Iterable[str],
    ) -> dict[str, Any]:
        """A ``similar_to`` edge is COMPUTED CONTEXT and never a causal edge."""
        return self.add_relation(
            subject,
            "similar_to",
            obj,
            evidence_ids=evidence_ids,
            claim_class="computed_measurement",
            embedding_model=embedding_model,
            distance=distance,
            metric=metric,
            search_universe=search_universe,
            embedding_slice_id=embedding_slice_id,
            causal=False,
            caveat="Computed embedding similarity. Not evidence of shared identity or causality.",
        )

    # -- reading --------------------------------------------------------------------

    def trail(self, node_id: str, max_hops: int = 2) -> list[dict[str, Any]]:
        """A readable evidence trail around one node, for the Inspect view."""
        seen_edges: set[int] = set()
        frontier = {node_id}
        out: list[dict[str, Any]] = []
        for hop in range(max_hops):
            nxt: set[str] = set()
            for n in sorted(frontier):
                for idx in self._adj.get(n, []):
                    if idx in seen_edges:
                        continue
                    seen_edges.add(idx)
                    e = self.edges[idx]
                    out.append(
                        {
                            "hop": hop + 1,
                            "subject": e["subject"],
                            "subject_label": self.nodes[e["subject"]]["label"],
                            "predicate": e["predicate"],
                            "object": e["object"],
                            "object_label": self.nodes[e["object"]]["label"],
                            "claim_class": e["claim_class"],
                            "disease_scope": e["disease_scope"],
                            "evidence": [
                                {
                                    "evidence_id": eid,
                                    "source": self.evidence[eid].source,
                                    "url": self.evidence[eid].source_url,
                                    "locator": self.evidence[eid].source_locator,
                                    "excerpt": self.evidence[eid].supporting_excerpt,
                                    "origin": self.evidence[eid].origin,
                                }
                                for eid in e["evidence_ids"]
                            ],
                            "attrs": e["attrs"],
                        }
                    )
                    nxt.add(e["subject"])
                    nxt.add(e["object"])
            frontier = nxt - frontier
        return out

    def contradictions(self, subject: str | None = None) -> list[dict[str, Any]]:
        """Contradicting evidence is preserved, not dropped. This surfaces it."""
        out = []
        for e in self.edges:
            for eid in e["evidence_ids"]:
                ev = self.evidence[eid]
                if ev.direction == "contradicts" and (subject is None or e["subject"] == subject):
                    out.append(
                        {
                            "subject": e["subject"],
                            "predicate": e["predicate"],
                            "object": e["object"],
                            "evidence_id": eid,
                            "source": ev.source,
                            "url": ev.source_url,
                            "excerpt": ev.supporting_excerpt,
                        }
                    )
        return out

    def stats(self) -> dict[str, Any]:
        by_kind: dict[str, int] = defaultdict(int)
        for n in self.nodes.values():
            by_kind[n["kind"]] += 1
        by_class: dict[str, int] = defaultdict(int)
        for e in self.edges:
            by_class[e["claim_class"]] += 1
        by_scope: dict[str, int] = defaultdict(int)
        for ev in self.evidence.values():
            by_scope[ev.disease_scope] += 1
        return {
            "nodes": len(self.nodes),
            "edges": len(self.edges),
            "evidence_records": len(self.evidence),
            "nodes_by_kind": dict(sorted(by_kind.items())),
            "edges_by_claim_class": dict(sorted(by_class.items())),
            "evidence_by_disease_scope": dict(sorted(by_scope.items())),
            "contradictions": len(self.contradictions()),
        }

    def to_json(self) -> str:
        return json.dumps(
            {
                "nodes": list(self.nodes.values()),
                "edges": self.edges,
                "evidence": {k: json.loads(v.model_dump_json()) for k, v in self.evidence.items()},
                "stats": self.stats(),
            },
            indent=2,
            default=str,
        )


def dedupe_evidence(evidences: list[Evidence]) -> tuple[list[Evidence], dict[str, int]]:
    """Collapse evidence that shares an underlying identifier.

    "Repeated appearances of the same study are not independent replication." Two Open
    Targets rows citing the same PMID for the same target-disease pair are one piece of
    evidence, and literature volume must not be inflated by counting them twice.
    """
    seen: dict[tuple, Evidence] = {}
    collapsed: dict[str, int] = defaultdict(int)
    for ev in evidences:
        pmid = None
        if ev.source_locator and "pmid=" in ev.source_locator:
            pmid = ev.source_locator.split("pmid=")[1].split()[0]
        key = (
            ev.subject_id,
            ev.object_id,
            ev.evidence_type,
            pmid or ev.doi or ev.source_url or ev.evidence_id,
        )
        if key in seen:
            collapsed[str(key[3])] += 1
            continue
        seen[key] = ev
    return list(seen.values()), dict(collapsed)
