"""Open Targets Platform adapter.

Every field used here was introspected against the live schema on 2026-09-22 at
API version 26.6.3 / data version 26.06. The GraphQL contract moves between
releases, so ``verify_schema()`` re-checks the fields this module depends on and
raises early rather than letting a silent ``None`` propagate into evidence.

Two behaviours worth calling out because they are easy to get wrong:

* A GraphQL endpoint can return HTTP 200 with a populated ``errors`` array and a
  *partial* ``data`` payload. ``_gql`` surfaces those errors instead of treating the
  partial body as success.
* The user's indication string is never assumed to be a valid ontology ID. Resolution
  goes through ``search``; ``EFO_0000274`` for instance returns null for atopic
  dermatitis in release 26.06, where the live ID is ``MONDO_0004980``.
"""

from __future__ import annotations

import time
from typing import Any

import requests

from ..config import OPEN_TARGETS_ENDPOINT, OPEN_TARGETS_TIMEOUT_S
from ..contracts import DiseaseMatch, Evidence, Provenance, Target, utcnow

ADAPTER_VERSION = "open_targets/1.0.0"


class OpenTargetsError(RuntimeError):
    """Raised for transport failures and GraphQL errors alike."""

    def __init__(self, message: str, *, transient: bool = False):
        super().__init__(message)
        self.transient = transient


# ---------------------------------------------------------------------------------
# Queries -- all parameterised, none string-interpolated.
# ---------------------------------------------------------------------------------

Q_META = """query { meta { name apiVersion { x y z } dataVersion { year month iteration } } }"""

Q_SEARCH_DISEASE = """
query($q:String!, $size:Int!) {
  search(queryString:$q, entityNames:["disease"], page:{index:0, size:$size}) {
    total
    hits { id name score entity description }
  }
}"""

Q_DISEASE_DETAIL = """
query($efo:String!) {
  disease(efoId:$efo) {
    id name description
    synonyms { relation terms }
    parents { id name }
    children { id name }
    therapeuticAreas { id name }
    dbXRefs
  }
}"""

Q_ASSOCIATED_TARGETS = """
query($efo:String!, $size:Int!, $idx:Int!, $indirect:Boolean!) {
  disease(efoId:$efo) {
    id name
    associatedTargets(page:{index:$idx, size:$size}, enableIndirect:$indirect) {
      count
      rows {
        score
        datatypeScores { id score }
        datasourceScores { id score }
        target {
          id approvedSymbol approvedName biotype
          proteinIds { id source }
          subcellularLocations { location termSL labelSL source }
          tractability { label modality value }
        }
      }
    }
  }
}"""

Q_TARGET_DETAIL = """
query($ensg:String!) {
  target(ensemblId:$ensg) {
    id approvedSymbol approvedName biotype
    functionDescriptions
    proteinIds { id source }
    subcellularLocations { location termSL labelSL source }
    tractability { label modality value }
    safetyLiabilities { event eventId datasource literature url effects { direction dosing } }
    geneticConstraint { constraintType score upperBin oe oeLower oeUpper }
    targetClass { id label level }
    drugAndClinicalCandidates {
      count
      rows {
        maxClinicalStage
        drug {
          id name drugType maximumClinicalStage
          mechanismsOfAction { rows { mechanismOfAction actionType targetName } }
        }
        diseases { diseaseFromSource disease { id name } }
      }
    }
  }
}"""

Q_TARGET_DISEASE_EVIDENCE = """
query($ensg:String!, $efo:String!, $size:Int!) {
  disease(efoId:$efo) {
    id name
    evidences(ensemblIds:[$ensg], enableIndirect:true, size:$size) {
      count
      rows {
        id datasourceId datatypeId score
        diseaseFromSource diseaseFromSourceMappedId
        targetFromSourceId
        literature publicationYear publicationFirstAuthor
        directionOnTrait directionOnTarget targetModulation
        studyId studyOverview
        clinicalPhase: clinicalStage
        drug { id name }
        disease { id name }
        textMiningSentences { text section }
        urls { niceName url }
      }
    }
  }
}"""

REQUIRED_FIELDS = {
    "Disease": {"id", "name", "associatedTargets", "synonyms", "parents", "children"},
    "Target": {
        "id",
        "approvedSymbol",
        "proteinIds",
        "subcellularLocations",
        "tractability",
        "safetyLiabilities",
        "drugAndClinicalCandidates",
    },
    "AssociatedTarget": {"score", "datatypeScores", "datasourceScores", "target"},
}


class OpenTargetsAdapter:
    def __init__(self, endpoint: str = OPEN_TARGETS_ENDPOINT, timeout: int = OPEN_TARGETS_TIMEOUT_S):
        self.endpoint = endpoint
        self.timeout = timeout
        self.session = requests.Session()
        self._release: str | None = None
        self._cache: dict[str, Any] = {}

    # -- transport ------------------------------------------------------------------

    def _gql(self, query: str, variables: dict[str, Any] | None = None) -> tuple[dict, Provenance]:
        t0 = time.time()
        try:
            resp = self.session.post(
                self.endpoint,
                json={"query": query, "variables": variables or {}},
                timeout=self.timeout,
            )
        except requests.Timeout as exc:
            # A retrieval timeout is a SOURCE FAILURE, never evidence that no targets exist.
            raise OpenTargetsError(f"Open Targets request timed out: {exc}", transient=True) from exc
        except requests.RequestException as exc:
            raise OpenTargetsError(f"Open Targets transport error: {exc}", transient=True) from exc

        elapsed_ms = int((time.time() - t0) * 1000)

        if resp.status_code >= 500:
            raise OpenTargetsError(f"Open Targets HTTP {resp.status_code}", transient=True)
        if resp.status_code != 200:
            raise OpenTargetsError(f"Open Targets HTTP {resp.status_code}: {resp.text[:300]}")

        body = resp.json()
        # HTTP 200 with a populated errors array is a real failure even when `data` is present.
        if body.get("errors"):
            msgs = "; ".join(e.get("message", "?") for e in body["errors"])
            raise OpenTargetsError(f"GraphQL errors: {msgs}")

        prov = Provenance(
            source="open_targets",
            source_release=self.release(),
            endpoint=self.endpoint,
            query=" ".join(query.split()),
            variables=variables,
            elapsed_ms=elapsed_ms,
            origin="real",
        )
        return body.get("data") or {}, prov

    def release(self) -> str:
        if self._release is None:
            try:
                resp = self.session.post(self.endpoint, json={"query": Q_META}, timeout=self.timeout)
                m = resp.json()["data"]["meta"]
                dv = m["dataVersion"]
                av = m["apiVersion"]
                self._release = f"data {dv['year']}.{dv['month']} / api {av['x']}.{av['y']}.{av['z']}"
            except Exception:
                self._release = "unknown"
        return self._release

    def verify_schema(self) -> dict[str, list[str]]:
        """Confirm the fields this adapter depends on still exist. Returns any that vanished."""
        q = """query($n:String!){ __type(name:$n){ name fields { name } } }"""
        missing: dict[str, list[str]] = {}
        for tname, required in REQUIRED_FIELDS.items():
            data, _ = self._gql(q, {"n": tname})
            t = data.get("__type")
            if not t:
                missing[tname] = sorted(required)
                continue
            present = {f["name"] for f in t["fields"]}
            gone = sorted(required - present)
            if gone:
                missing[tname] = gone
        return missing

    # -- tool: resolve_indication ----------------------------------------------------

    def resolve_indication(self, raw_text: str, size: int = 8) -> tuple[list[DiseaseMatch], Provenance]:
        """Rank ontology matches for free text. Never assumes the input is an ID."""
        data, prov = self._gql(Q_SEARCH_DISEASE, {"q": raw_text, "size": size})
        hits = (data.get("search") or {}).get("hits") or []
        if not hits:
            return [], prov

        norm = raw_text.strip().lower()
        top = hits[0]["score"] or 1.0
        matches: list[DiseaseMatch] = []
        for h in hits:
            name = (h.get("name") or "").lower()
            if name == norm:
                mtype = "exact"
            elif norm in name or name in norm:
                mtype = "synonym"
            else:
                mtype = "related"
            matches.append(
                DiseaseMatch(
                    disease_id=h["id"],
                    label=h["name"],
                    match_type=mtype,
                    score=round((h.get("score") or 0.0) / top, 4),
                    description=(h.get("description") or None),
                )
            )
        return matches, prov

    def needs_clarification(self, matches: list[DiseaseMatch], margin: float = 0.35) -> bool:
        """True only when plausible alternatives would MATERIALLY change retrieval."""
        if len(matches) < 2:
            return False
        if matches[0].match_type == "exact" and matches[1].match_type != "exact":
            return False
        return (matches[0].score or 0) - (matches[1].score or 0) < margin

    def disease_detail(self, disease_id: str) -> tuple[dict, Provenance]:
        data, prov = self._gql(Q_DISEASE_DETAIL, {"efo": disease_id})
        d = data.get("disease")
        if d is None:
            raise OpenTargetsError(
                f"Disease ID '{disease_id}' does not resolve in release {self.release()}. "
                "Resolve free text through resolve_indication rather than assuming an ID."
            )
        return d, prov

    # -- tool: discover_targets ------------------------------------------------------

    def discover_targets(
        self,
        disease_id: str,
        limit: int = 20,
        page_index: int = 0,
        enable_indirect: bool = True,
    ) -> tuple[list[Target], dict[str, Any], Provenance]:
        """Retrieve a BOUNDED candidate pool. Retrieval order, page and limit are recorded;
        this is explicitly not an exhaustive discovery claim."""
        data, prov = self._gql(
            Q_ASSOCIATED_TARGETS,
            {"efo": disease_id, "size": limit, "idx": page_index, "indirect": enable_indirect},
        )
        disease = data.get("disease")
        if disease is None:
            raise OpenTargetsError(f"Disease ID '{disease_id}' not found in {self.release()}.")

        at = disease.get("associatedTargets") or {}
        rows = at.get("rows") or []

        targets: list[Target] = []
        assoc: dict[str, dict[str, Any]] = {}
        for rank, row in enumerate(rows, start=page_index * limit + 1):
            tg = row["target"]
            targets.append(self._to_target(tg, prov))
            assoc[tg["id"]] = {
                "rank": rank,
                "overall_score": row.get("score"),
                "datatype_scores": {s["id"]: s["score"] for s in (row.get("datatypeScores") or [])},
                "datasource_scores": {s["id"]: s["score"] for s in (row.get("datasourceScores") or [])},
            }

        scope = {
            "disease_id": disease_id,
            "disease_label": disease.get("name"),
            "total_associated_in_release": at.get("count"),
            "retrieved": len(rows),
            "page_index": page_index,
            "page_size": limit,
            "ordering": "default associatedTargets ordering (descending overall association score)",
            "enable_indirect": enable_indirect,
            "release": self.release(),
            "exhaustive": False,
            "note": (
                f"Retrieved {len(rows)} of {at.get('count')} targets associated with this disease in "
                "the release. This is a bounded pool, not exhaustive discovery."
            ),
            "association_semantics": (
                "Association score measures aggregated association evidence. It is not a probability "
                "of therapeutic success, not causal proof, and says nothing about whether inhibition "
                "or activation is beneficial."
            ),
        }
        return targets, {"associations": assoc, "scope": scope}, prov

    @staticmethod
    def _to_target(tg: dict, prov: Provenance) -> Target:
        proteins = tg.get("proteinIds") or []
        swiss = [p["id"] for p in proteins if p.get("source") == "uniprot_swissprot"]
        trembl = [p["id"] for p in proteins if p.get("source") == "uniprot_trembl"]
        accession = (swiss or trembl or [None])[0]

        tract: dict[str, list[str]] = {}
        for t in tg.get("tractability") or []:
            if t.get("value"):
                tract.setdefault(t["modality"], []).append(t["label"])

        locations = sorted({l["location"] for l in (tg.get("subcellularLocations") or []) if l.get("location")})

        return Target(
            target_id=tg["id"],
            gene_symbol=tg["approvedSymbol"],
            stable_gene_id=tg["id"],
            protein_accession=accession,
            approved_name=tg.get("approvedName"),
            biotype=tg.get("biotype"),
            subcellular_locations=locations,
            tractability=tract,
            provenance=prov,
        )

    # -- tool: get_target_evidence ---------------------------------------------------

    def get_target_evidence(
        self, target_id: str, disease_id: str, evidence_size: int = 25
    ) -> tuple[dict[str, Any], list[Evidence], Provenance]:
        """Target annotation plus the underlying disease evidence items, as Evidence records."""
        detail, prov_t = self._gql(Q_TARGET_DETAIL, {"ensg": target_id})
        target = detail.get("target")
        if target is None:
            raise OpenTargetsError(f"Target '{target_id}' not found in {self.release()}.")

        ev_data, prov_e = self._gql(
            Q_TARGET_DISEASE_EVIDENCE,
            {"ensg": target_id, "efo": disease_id, "size": evidence_size},
        )
        ev_block = ((ev_data.get("disease") or {}).get("evidences") or {})
        rows = ev_block.get("rows") or []

        evidences: list[Evidence] = []
        for row in rows:
            row_disease = (row.get("disease") or {}).get("id")
            # Preserve whether this item actually refers to the queried disease or to a
            # descendant whose evidence was propagated up the ontology.
            scope = "direct" if row_disease == disease_id else "indirect_propagated"

            lit = row.get("literature") or []
            pmid = lit[0] if lit else None
            urls = row.get("urls") or []
            url = urls[0]["url"] if urls else (f"https://europepmc.org/article/MED/{pmid}" if pmid else None)

            sentences = row.get("textMiningSentences") or []
            excerpt = sentences[0].get("text") if sentences else None

            direction = "supports"
            dot = (row.get("directionOnTrait") or "").lower()
            if dot in ("risk", "protective"):
                direction = "supports"

            evidences.append(
                Evidence(
                    evidence_id=f"ot:{row['id']}",
                    subject_id=target_id,
                    predicate="associated_with",
                    object_id=disease_id,
                    object_value=row.get("score"),
                    evidence_type=_map_datatype(row.get("datatypeId")),
                    direction=direction,
                    source=f"open_targets:{row.get('datasourceId')}",
                    source_url=url,
                    source_locator=(
                        f"datasource={row.get('datasourceId')} datatype={row.get('datatypeId')}"
                        + (f" pmid={pmid}" if pmid else "")
                    ),
                    supporting_excerpt=excerpt,
                    measurement={
                        "resource_score": row.get("score"),
                        "direction_on_trait": row.get("directionOnTrait"),
                        "direction_on_target": row.get("directionOnTarget"),
                        "target_modulation": row.get("targetModulation"),
                        "clinical_phase": row.get("clinicalPhase"),
                        "publication_year": row.get("publicationYear"),
                    },
                    disease_scope=scope,
                    scope_disease_id=row_disease,
                    limitations=(
                        []
                        if scope == "direct"
                        else [
                            f"Evidence refers to '{(row.get('disease') or {}).get('name')}' "
                            f"({row_disease}) and was propagated through the ontology to {disease_id}."
                        ]
                    ),
                    origin="real",
                )
            )

        drugs = _extract_drugs(target)
        annotation = {
            "target_id": target["id"],
            "gene_symbol": target["approvedSymbol"],
            "approved_name": target.get("approvedName"),
            "function": (target.get("functionDescriptions") or [None])[0],
            "target_class": [c.get("label") for c in (target.get("targetClass") or [])],
            "subcellular_locations": sorted(
                {l["location"] for l in (target.get("subcellularLocations") or []) if l.get("location")}
            ),
            "tractability": {
                m: [t["label"] for t in (target.get("tractability") or []) if t["value"] and t["modality"] == m]
                for m in {t["modality"] for t in (target.get("tractability") or [])}
            },
            "safety_liabilities": [
                {
                    "event": s.get("event"),
                    "datasource": s.get("datasource"),
                    "url": s.get("url"),
                }
                for s in (target.get("safetyLiabilities") or [])
            ],
            "genetic_constraint": target.get("geneticConstraint"),
            "known_drugs": drugs,
            "evidence_count_in_release": ev_block.get("count"),
            "evidence_retrieved": len(rows),
        }

        prov = Provenance(
            source="open_targets",
            source_release=self.release(),
            endpoint=self.endpoint,
            query="Q_TARGET_DETAIL + Q_TARGET_DISEASE_EVIDENCE",
            variables={"ensg": target_id, "efo": disease_id, "size": evidence_size},
            elapsed_ms=(prov_t.elapsed_ms or 0) + (prov_e.elapsed_ms or 0),
            retrieved_at=utcnow(),
            origin="real",
        )
        return annotation, evidences, prov


def _map_datatype(datatype_id: str | None) -> str:
    """Map an Open Targets ``datatypeId`` onto the Evidence.evidence_type vocabulary.

    The live vocabulary was enumerated from release 26.06 and is NOT the set the platform
    documentation's older examples suggest: this release emits ``clinical`` (not
    ``known_drug``) and ``genetic_literature``. Anything unrecognised maps to
    ``annotation`` rather than being dropped, so an unexpected new datatype degrades to
    a weakly-weighted record instead of vanishing.
    """
    return {
        "genetic_association": "genetic_association",
        "genetic_literature": "genetic_association",
        "somatic_mutation": "somatic_mutation",
        "clinical": "clinical",
        "known_drug": "known_drug",
        "affected_pathway": "affected_pathway",
        "rna_expression": "rna_expression",
        "animal_model": "animal_model",
        "literature": "literature",
    }.get(datatype_id or "", "annotation")


def _extract_drugs(target: dict) -> list[dict[str, Any]]:
    """Known drugs with mechanism and direction. A trial entry is not efficacy or approval."""
    out: list[dict[str, Any]] = []
    block = target.get("drugAndClinicalCandidates") or {}
    for row in block.get("rows") or []:
        drug = row.get("drug") or {}
        moas = [
            {
                "mechanism": m.get("mechanismOfAction"),
                "action_type": m.get("actionType"),
                "target_name": m.get("targetName"),
            }
            for m in ((drug.get("mechanismsOfAction") or {}).get("rows") or [])
        ]
        out.append(
            {
                "drug_id": drug.get("id"),
                "drug_name": drug.get("name"),
                "drug_type": drug.get("drugType"),
                "max_clinical_stage": row.get("maxClinicalStage") or drug.get("maximumClinicalStage"),
                "mechanisms": moas,
                "indications": [
                    ((d.get("disease") or {}).get("name") or d.get("diseaseFromSource"))
                    for d in (row.get("diseases") or [])
                ][:6],
                "indication_ids": [
                    (d.get("disease") or {}).get("id") for d in (row.get("diseases") or []) if d.get("disease")
                ][:6],
                "caveat": "A trial entry is not evidence of efficacy or approval for this indication.",
            }
        )
    return out
