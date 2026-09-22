"""Open Targets Platform GraphQL client.

Field names below are checked against the live schema by ``introspect_type``;
the pipeline runs that as a self-test rather than trusting a hard-coded
schema, because Open Targets renames fields between releases
(``knownDrugs`` -> ``drugAndClinicalCandidates``, for example).
"""

from __future__ import annotations

import logging
from typing import Any, Iterable

from ..config import OPEN_TARGETS_GRAPHQL
from ..http import request_json

log = logging.getLogger("ind2b.opentargets")


class OpenTargetsError(RuntimeError):
    pass


def graphql(query: str, variables: dict[str, Any] | None = None, *, use_cache: bool = True) -> dict:
    payload = {"query": query, "variables": variables or {}}
    data = request_json(
        OPEN_TARGETS_GRAPHQL, method="POST", json_body=payload, use_cache=use_cache
    )
    if data.get("errors"):
        raise OpenTargetsError(str(data["errors"])[:500])
    if data.get("data") is None:
        raise OpenTargetsError(f"null data for query: {query[:120]}")
    return data["data"]


def introspect_type(type_name: str) -> list[str]:
    """Return the field names of a GraphQL type (schema self-check)."""
    q = "query($n:String!){ __type(name:$n){ name fields{ name } } }"
    data = graphql(q, {"n": type_name})
    t = data.get("__type")
    if not t:
        raise OpenTargetsError(f"unknown type {type_name}")
    return [f["name"] for f in (t.get("fields") or [])]


# --- stage 0 ---------------------------------------------------------------

def search_disease(query_string: str, size: int = 10) -> list[dict]:
    """Free-text disease search; returns ranked ontology hits."""
    q = """
    query($qs:String!,$size:Int!){
      search(queryString:$qs, entityNames:["disease"], page:{size:$size,index:0}){
        total
        hits{ id name description entity }
      }
    }
    """
    data = graphql(q, {"qs": query_string, "size": size})
    return data["search"]["hits"]


def diseases_synonyms(efo_ids: list[str]) -> dict[str, dict[str, list[str]]]:
    """Fetch ontology synonyms for several disease ids in one request.

    Returns ``{efo_id: {relation: [terms]}}``. The relation is kept rather than
    flattened because the relations do not mean the same thing:
    ``hasExactSynonym`` asserts the terms name the *same* disease, while
    ``hasBroadSynonym`` and ``hasNarrowSynonym`` name a more general or more
    specific one. Only the exact relation can be treated as identity.

    Batched deliberately: term selection compares a query against every
    candidate hit, and one request for all of them keeps that to a single
    cached call instead of one per candidate.
    """
    if not efo_ids:
        return {}
    q = """
    query($ids:[String!]!){
      diseases(efoIds:$ids){ id name synonyms { relation terms } }
    }
    """
    data = graphql(q, {"ids": list(efo_ids)})
    out: dict[str, dict[str, list[str]]] = {}
    for dis in data.get("diseases") or []:
        rel: dict[str, list[str]] = {}
        for block in dis.get("synonyms") or []:
            relation = block.get("relation")
            if relation:
                rel[relation] = list(block.get("terms") or [])
        out[dis["id"]] = rel
    return out


def disease_info(efo_id: str) -> dict:
    """Disease record with the descendant subtree used to widen the pull."""
    q = """
    query($efo:String!){
      disease(efoId:$efo){
        id name description
        therapeuticAreas{ id name }
        parents{ id name }
        children{ id name }
        descendants
      }
    }
    """
    data = graphql(q, {"efo": efo_id})
    if data.get("disease") is None:
        raise OpenTargetsError(f"no disease for id {efo_id}")
    return data["disease"]


# --- stage 1 ---------------------------------------------------------------

def associated_targets(efo_id: str, size: int = 100, index: int = 0) -> dict:
    """Disease-target associations with the per-datatype evidence breakdown."""
    q = """
    query($efo:String!,$size:Int!,$index:Int!){
      disease(efoId:$efo){
        id name
        associatedTargets(page:{size:$size,index:$index}){
          count
          rows{
            score
            datatypeScores{ id score }
            target{ id approvedSymbol approvedName biotype }
          }
        }
      }
    }
    """
    data = graphql(q, {"efo": efo_id, "size": size, "index": index})
    return data["disease"]["associatedTargets"]


def targets_details(ensembl_ids: Iterable[str]) -> list[dict]:
    """Bulk target detail: tractability, safety, localization, constraint.

    Uses the bulk ``targets(ensemblIds:)`` field so a 100-target pull is one
    request rather than 100.
    """
    ids = list(ensembl_ids)
    if not ids:
        return []
    q = """
    query($ids:[String!]!){
      targets(ensemblIds:$ids){
        id approvedSymbol approvedName biotype
        proteinIds{ id source }
        functionDescriptions
        subcellularLocations{ location source termSL labelSL }
        tractability{ label modality value }
        safetyLiabilities{ event eventId datasource }
        geneticConstraint{ constraintType score upperBin }
        targetClass{ id label level }
      }
    }
    """
    data = graphql(q, {"ids": ids})
    return data["targets"]


def target_expression(ensembl_id: str) -> list[dict]:
    """Baseline expression rows per tissue/cell type.

    The schema field is ``baselineExpression`` (there is no ``expressions``
    field). Each row carries an upstream ``specificity_score`` and
    ``distribution_score``, which stage 2 uses directly rather than
    recomputing a specificity measure from per-tissue values.
    """
    q = """
    query($id:String!){
      target(ensemblId:$id){
        id approvedSymbol
        baselineExpression{
          count
          rows{
            datasourceId datatypeId unit
            median q1 q3 min max
            specificity_score distribution_score
            tissueBiosampleFromSource celltypeBiosampleFromSource
          }
        }
      }
    }
    """
    data = graphql(q, {"id": ensembl_id})
    t = data.get("target") or {}
    block = t.get("baselineExpression") or {}
    return block.get("rows") or []


def target_known_drugs(ensembl_id: str, size: int = 50) -> dict:
    """Drugs and clinical candidates against a target.

    The upstream field is unpaginated, so rows are sliced client-side.
    """
    q = """
    query($id:String!){
      target(ensemblId:$id){
        id approvedSymbol
        drugAndClinicalCandidates{
          count
          rows{
            id
            maxClinicalStage
            drug{ id name drugType maximumClinicalStage }
          }
        }
      }
    }
    """
    data = graphql(q, {"id": ensembl_id})
    t = data.get("target") or {}
    block = t.get("drugAndClinicalCandidates") or {"count": 0, "rows": []}
    rows = (block.get("rows") or [])[:size]
    return {"count": block.get("count", len(rows)), "rows": rows}


def target_interactions(ensembl_id: str, size: int = 50) -> list[dict]:
    """Protein-protein interaction partners (IntAct/Reactome/Signor via OT).

    Stage 3 uses these to decide which partner interface is the
    mechanistically relevant one to block.
    """
    q = """
    query($id:String!,$size:Int!){
      target(ensemblId:$id){
        id approvedSymbol
        interactions(page:{size:$size,index:0}){
          count
          rows{
            intA intB
            intABiologicalRole intBBiologicalRole
            sourceDatabase
            score
            count
            targetB{ id approvedSymbol }
          }
        }
      }
    }
    """
    data = graphql(q, {"id": ensembl_id, "size": size})
    t = data.get("target") or {}
    block = t.get("interactions") or {}
    return block.get("rows") or []
