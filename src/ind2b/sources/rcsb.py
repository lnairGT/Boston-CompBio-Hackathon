"""RCSB PDB client - experimental complex search and structure download.

The structural half of the gate is pushed into the search query (method,
resolution, at least two protein entities) so that only plausible complexes
are described in detail. EGFR alone has ~400 PDB entries, most of them
kinase-domain-plus-inhibitor structures that cannot inform a binder against
the ectodomain; filtering server-side keeps the detail queries small.

The mechanistic half - is the partner the one we want to block - cannot be
expressed as a search filter and runs client-side in ``stage3_complexes``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Iterable, Sequence

from ..config import RCSB_DATA_GRAPHQL, RCSB_FILES, RCSB_SEARCH
from ..http import download, request_json

log = logging.getLogger("ind2b.rcsb")

_DETAIL_CHUNK = 40

_ENTRY_QUERY = """
query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    exptl{ method }
    struct{ title }
    rcsb_entry_info{
      resolution_combined
      deposited_polymer_entity_instance_count
      polymer_entity_count_protein
    }
    polymer_entities{
      rcsb_id
      rcsb_polymer_entity{ pdbx_description }
      entity_poly{ rcsb_sample_sequence_length type }
      rcsb_polymer_entity_container_identifiers{
        auth_asym_ids asym_ids
        reference_sequence_identifiers{ database_accession database_name }
      }
    }
  }
}
"""


def search_complexes(
    accession: str,
    *,
    max_resolution: float | None = 4.0,
    min_protein_entities: int = 2,
    rows: int = 100,
) -> dict[str, Any]:
    """Entry ids for experimental complexes containing ``accession``.

    Returns ``{"total": int, "ids": [...]}``. Sorted by resolution so that,
    when the row cap bites, the best-resolved complexes are the ones kept.
    """
    nodes: list[dict[str, Any]] = [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers"
                ".reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": accession,
            },
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers"
                ".reference_sequence_identifiers.database_name",
                "operator": "exact_match",
                "value": "UniProt",
            },
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                "operator": "greater_or_equal",
                "value": min_protein_entities,
            },
        },
    ]
    if max_resolution is not None:
        nodes.append(
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "less_or_equal",
                    "value": max_resolution,
                },
            }
        )

    query = {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": rows},
            "results_verbosity": "compact",
            "sort": [
                {
                    "sort_by": "rcsb_entry_info.resolution_combined",
                    "direction": "asc",
                }
            ],
        },
    }
    try:
        data = request_json(RCSB_SEARCH, method="POST", json_body=query)
    except Exception as exc:  # noqa: BLE001 - reported, not raised
        log.warning("rcsb search failed for %s: %s", accession, exc)
        return {"total": 0, "ids": [], "note": f"search failed: {type(exc).__name__}"}
    if not data:
        return {"total": 0, "ids": []}
    return {
        "total": data.get("total_count", 0),
        "ids": list(data.get("result_set") or []),
    }


def entry_details(entry_ids: Sequence[str]) -> list[dict]:
    """Batched entry metadata: method, resolution, entities, chains, accessions."""
    out: list[dict] = []
    ids = [i for i in entry_ids if i]
    for i in range(0, len(ids), _DETAIL_CHUNK):
        chunk = ids[i : i + _DETAIL_CHUNK]
        data = request_json(
            RCSB_DATA_GRAPHQL,
            method="POST",
            json_body={"query": _ENTRY_QUERY, "variables": {"ids": chunk}},
        )
        if data.get("errors"):
            log.warning("rcsb data errors: %s", str(data["errors"])[:200])
        for e in (data.get("data") or {}).get("entries") or []:
            if e:
                out.append(e)
    return out


_ALIGN_QUERY = """
query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    polymer_entities{
      rcsb_id
      rcsb_polymer_entity_align{
        reference_database_accession reference_database_name
        aligned_regions{ entity_beg_seq_id ref_beg_seq_id length }
      }
      polymer_entity_instances{
        rcsb_polymer_entity_instance_container_identifiers{
          auth_asym_id asym_id auth_to_entity_poly_seq_mapping
        }
      }
    }
  }
}
"""


def entry_alignments(entry_ids: Sequence[str]) -> dict[str, dict]:
    """Per-chain author-numbering to UniProt-numbering mappings.

    PDB author numbering is not guaranteed to equal UniProt numbering -
    constructs, tags and non-canonical isoforms all break the assumption - so
    any claim that an interface lies inside an extracellular span has to be
    made in a consistent coordinate system. This returns, per entry, the
    SIFTS-derived alignment regions and the author-id list per chain needed to
    convert between the two.
    """
    out: dict[str, dict] = {}
    ids = [i for i in entry_ids if i]
    for i in range(0, len(ids), _DETAIL_CHUNK):
        chunk = ids[i : i + _DETAIL_CHUNK]
        data = request_json(
            RCSB_DATA_GRAPHQL,
            method="POST",
            json_body={"query": _ALIGN_QUERY, "variables": {"ids": chunk}},
        )
        if data.get("errors"):
            log.warning("rcsb align errors: %s", str(data["errors"])[:200])
        for e in (data.get("data") or {}).get("entries") or []:
            if not e:
                continue
            entities = []
            for pe in e.get("polymer_entities") or []:
                aligns = []
                for a in pe.get("rcsb_polymer_entity_align") or []:
                    if a.get("reference_database_name") != "UniProt":
                        continue
                    aligns.append(
                        {
                            "accession": a.get("reference_database_accession"),
                            "regions": [
                                {
                                    "entity_beg": r["entity_beg_seq_id"],
                                    "ref_beg": r["ref_beg_seq_id"],
                                    "length": r["length"],
                                }
                                for r in (a.get("aligned_regions") or [])
                                if r.get("entity_beg_seq_id") and r.get("ref_beg_seq_id")
                            ],
                        }
                    )
                chains = {}
                for inst in pe.get("polymer_entity_instances") or []:
                    ci = inst.get("rcsb_polymer_entity_instance_container_identifiers") or {}
                    if ci.get("auth_asym_id"):
                        chains[ci["auth_asym_id"]] = list(
                            ci.get("auth_to_entity_poly_seq_mapping") or []
                        )
                entities.append(
                    {"entity_id": pe.get("rcsb_id"), "aligns": aligns, "chains": chains}
                )
            out[e["rcsb_id"]] = {"entities": entities}
    return out


def download_structure(entry_id: str, dest_dir: Path, fmt: str = "cif") -> Path:
    """Fetch one structure file (mmCIF by default)."""
    name = f"{entry_id.upper()}.{fmt}"
    return download(f"{RCSB_FILES}/{name}", dest_dir / name)
