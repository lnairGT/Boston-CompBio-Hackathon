"""Protein identity, binder accessibility and structure assessment.

Why this module exists as a separate screen: Open Targets ``subcellularLocations``
is free text drawn from several sources, and matching keywords against it gives
false positives that matter. PDE4A, PDE4B and PDE4D are all annotated "plasma
membrane" or "apical cell membrane" because they dock on the *cytoplasmic* face --
a keyword screen calls them accessible and would aim a soluble binder at a cytosolic
phosphodiesterase.

So accessibility is decided on **topology**, from UniProt features, not on location
strings:

* secreted   -- a signal peptide and no transmembrane segment
* ectodomain -- at least one topological domain annotated "Extracellular"

Everything else is inaccessible to a soluble binder, however it is described. The
location strings are retained as a weak corroborating prior and are reported, but
they never decide the verdict.
"""

from __future__ import annotations

import hashlib
import time
from typing import Any

import requests

from ..contracts import Evidence, Provenance

ADAPTER_VERSION = "structures/1.0.0"

UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{acc}.json"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"

UNIPROT_FIELDS = ",".join(
    [
        "accession",
        "id",
        "protein_name",
        "gene_names",
        "organism_name",
        "sequence",
        "ft_signal",
        "ft_topo_dom",
        "ft_transmem",
        "ft_chain",
        "ft_carbohyd",
        "ft_disulfid",
        "ft_domain",
        "cc_subcellular_location",
        "keyword",
        "xref_pdb",
    ]
)


class StructureAdapterError(RuntimeError):
    def __init__(self, message: str, *, transient: bool = False):
        super().__init__(message)
        self.transient = transient


def _get(url: str, **kw) -> requests.Response:
    try:
        return requests.get(url, timeout=45, **kw)
    except requests.Timeout as exc:
        raise StructureAdapterError(f"timeout: {url}", transient=True) from exc
    except requests.RequestException as exc:
        raise StructureAdapterError(f"transport error: {url}: {exc}", transient=True) from exc


# ---------------------------------------------------------------------------------
# UniProt identity and topology
# ---------------------------------------------------------------------------------


def fetch_protein(accession: str) -> tuple[dict[str, Any], Provenance]:
    """Canonical sequence, topology features and PDB cross-references for one accession."""
    t0 = time.time()
    r = _get(UNIPROT_ENTRY.format(acc=accession), params={"fields": UNIPROT_FIELDS})
    if r.status_code == 404:
        raise StructureAdapterError(f"UniProt accession '{accession}' not found.")
    if r.status_code != 200:
        raise StructureAdapterError(f"UniProt HTTP {r.status_code} for {accession}", transient=r.status_code >= 500)

    # UniProt answers an unknown or malformed accession with HTTP 200 and an EMPTY body,
    # not a 404. Returning the resulting empty record would be actively dangerous: a
    # protein with no sequence and no features reads as "no signal peptide, no
    # transmembrane, no extracellular domain", i.e. `not_accessible`, and the target gets
    # excluded for a modality mismatch it may not have. Absence of a record is a source
    # failure and must raise, never degrade into a negative finding.
    try:
        d = r.json()
    except ValueError as exc:
        raise StructureAdapterError(
            f"UniProt returned HTTP 200 with a non-JSON body for '{accession}'.", transient=True
        ) from exc
    if not d or not isinstance(d, dict) or not d.get("primaryAccession"):
        raise StructureAdapterError(
            f"UniProt returned HTTP 200 with an empty record for '{accession}'. The accession "
            "is unknown or malformed. This is a retrieval failure, NOT evidence that the "
            "protein lacks topology features."
        )

    seq = (d.get("sequence") or {}).get("value", "")
    feats = d.get("features") or []
    if not seq:
        raise StructureAdapterError(
            f"UniProt record for '{accession}' carries no sequence, so topology and "
            "accessibility cannot be assessed. Treat as unknown, not inaccessible."
        )

    def _span(f: dict) -> tuple[int | None, int | None]:
        loc = f.get("location") or {}
        return (loc.get("start") or {}).get("value"), (loc.get("end") or {}).get("value")

    signal = [dict(zip(("start", "end"), _span(f))) for f in feats if f["type"] == "Signal"]
    transmem = [
        {"start": _span(f)[0], "end": _span(f)[1], "description": f.get("description")}
        for f in feats
        if f["type"] == "Transmembrane"
    ]
    topo = [
        {"description": f.get("description"), "start": _span(f)[0], "end": _span(f)[1]}
        for f in feats
        if f["type"] == "Topological domain"
    ]
    glyco = [{"start": _span(f)[0], "description": f.get("description")} for f in feats if f["type"] == "Glycosylation"]
    disulf = [{"start": _span(f)[0], "end": _span(f)[1]} for f in feats if f["type"] == "Disulfide bond"]
    chains = [{"start": _span(f)[0], "end": _span(f)[1], "description": f.get("description")} for f in feats if f["type"] == "Chain"]

    keywords = [k.get("name") for k in (d.get("keywords") or [])]

    pdbs: list[dict[str, Any]] = []
    for x in d.get("uniProtKBCrossReferences") or []:
        if x.get("database") != "PDB":
            continue
        props = {p["key"]: p["value"] for p in x.get("properties", [])}
        res = props.get("Resolution", "")
        try:
            resolution = float(res.split()[0]) if res and res[0].isdigit() else None
        except (ValueError, IndexError):
            resolution = None
        pdbs.append(
            {
                "pdb_id": x.get("id"),
                "method": props.get("Method"),
                "resolution_A": resolution,
                "chains": props.get("Chains"),
            }
        )

    subcell = []
    for c in d.get("comments") or []:
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            for loc in c.get("subcellularLocations") or []:
                v = (loc.get("location") or {}).get("value")
                if v:
                    subcell.append(v)

    info = {
        "accession": d.get("primaryAccession", accession),
        "entry_name": d.get("uniProtkbId"),
        "protein_name": (((d.get("proteinDescription") or {}).get("recommendedName") or {}).get("fullName") or {}).get("value"),
        "organism": (d.get("organism") or {}).get("scientificName"),
        "taxon_id": (d.get("organism") or {}).get("taxonId"),
        "sequence": seq,
        "length": len(seq),
        "sequence_checksum": "sha256:" + hashlib.sha256(seq.encode()).hexdigest()[:32],
        "signal_peptide": signal,
        "transmembrane": transmem,
        "topological_domains": topo,
        "glycosylation_sites": glyco,
        "disulfide_bonds": disulf,
        "chains": chains,
        "keywords": keywords,
        "subcellular_locations": subcell,
        "pdb_entries": sorted(pdbs, key=lambda p: (p["resolution_A"] is None, p["resolution_A"] or 99)),
    }
    prov = Provenance(
        source="uniprot",
        source_release="UniProtKB REST (live)",
        endpoint=UNIPROT_ENTRY.format(acc=accession),
        query=UNIPROT_FIELDS,
        elapsed_ms=int((time.time() - t0) * 1000),
        origin="real",
    )
    return info, prov


# ---------------------------------------------------------------------------------
# Accessibility verdict
# ---------------------------------------------------------------------------------


def assess_accessibility(protein: dict[str, Any]) -> dict[str, Any]:
    """Decide soluble-binder accessibility from topology. Returns a verdict plus the
    residue ranges a binder could actually engage."""
    sig = protein.get("signal_peptide") or []
    tms = protein.get("transmembrane") or []
    topo = protein.get("topological_domains") or []
    kws = protein.get("keywords") or []
    length = protein.get("length") or 0

    extracellular = [t for t in topo if "xtracellular" in (t.get("description") or "")]

    verdict: str
    basis: list[str] = []
    engageable: list[dict[str, Any]] = []

    if extracellular:
        verdict = "accessible_ectodomain"
        for t in extracellular:
            basis.append(f"UniProt topological domain 'Extracellular' at {t['start']}-{t['end']}")
            engageable.append({"kind": "ectodomain", "start": t["start"], "end": t["end"]})
        basis.append(f"{len(tms)} transmembrane segment(s) annotated")
    elif sig and not tms:
        verdict = "accessible_secreted"
        mature_start = (sig[0]["end"] or 0) + 1
        basis.append(f"Signal peptide 1-{sig[0]['end']} and no transmembrane segment")
        if "Secreted" in kws:
            basis.append("UniProt keyword 'Secreted'")
        engageable.append({"kind": "mature_secreted", "start": mature_start, "end": length})
    elif sig and tms:
        verdict = "accessible_ectodomain"
        basis.append("Signal peptide with transmembrane segment, but no explicit extracellular "
                     "topological domain annotation; ectodomain boundaries are inferred and unverified")
        first_tm = min(t["start"] for t in tms if t["start"])
        engageable.append({"kind": "inferred_ectodomain", "start": (sig[0]["end"] or 0) + 1, "end": first_tm - 1})
    else:
        verdict = "not_accessible"
        membraneish = [k for k in kws if k in ("Membrane", "Cell membrane")]
        if membraneish:
            basis.append(
                f"UniProt keyword(s) {membraneish} are present, but there is no signal peptide, no "
                "transmembrane segment and no extracellular topological domain: this is a "
                "membrane-ASSOCIATED protein on the cytoplasmic face, not a surface-exposed one."
            )
        else:
            basis.append("No signal peptide, no transmembrane segment and no extracellular topological domain")

    return {
        "verdict": verdict,
        "accessible": verdict != "not_accessible",
        "basis": basis,
        "engageable_regions": engageable,
        "n_glycosylation_sites": len(protein.get("glycosylation_sites") or []),
        "n_disulfides": len(protein.get("disulfide_bonds") or []),
        "caveat": (
            "Accessibility is inferred from UniProt topology annotation, not measured. "
            "Glycan shielding, oligomeric state and epitope occlusion in vivo are not assessed here."
        ),
    }


# ---------------------------------------------------------------------------------
# Structures
# ---------------------------------------------------------------------------------


def alphafold_coverage(accession: str) -> dict[str, Any] | None:
    r = _get(ALPHAFOLD_API.format(acc=accession))
    if r.status_code != 200:
        return None
    entries = r.json()
    if not entries:
        return None
    e = entries[0]
    return {
        "source": "alphafold",
        "model_id": e.get("entryId"),
        "version": e.get("latestVersion"),
        "model_url": e.get("pdbUrl"),
        "cif_url": e.get("cifUrl"),
        "uniprot_start": e.get("uniprotStart"),
        "uniprot_end": e.get("uniprotEnd"),
        "model_created": e.get("modelCreatedDate"),
        "limitation": "Predicted model. Per-residue pLDDT should be checked before using a region as a design target.",
    }


def assess_structure(accession: str, mechanism_constraint: str | None = None) -> tuple[dict[str, Any], list[Evidence], Provenance]:
    """Combined identity + accessibility + structure assessment for one protein."""
    protein, prov = fetch_protein(accession)
    access = assess_accessibility(protein)

    experimental = protein["pdb_entries"]
    af = alphafold_coverage(accession)

    limitations: list[str] = [access["caveat"]]
    if not experimental:
        limitations.append("No experimental structure cross-referenced in UniProt for this accession.")
    if access["n_glycosylation_sites"]:
        limitations.append(
            f"{access['n_glycosylation_sites']} annotated glycosylation site(s); glycans are absent from "
            "both experimental constructs and predicted models and may occlude an epitope."
        )
    if protein["length"] > 1500:
        limitations.append(
            f"Full-length chain is {protein['length']} residues; a design construct must be a defined "
            "domain, not the full chain."
        )

    evidence = [
        Evidence(
            evidence_id=f"uniprot:{accession}:topology",
            subject_id=accession,
            predicate="has_accessibility",
            object_value=access["verdict"],
            evidence_type="annotation",
            direction="context",
            source="uniprot",
            source_url=f"https://www.uniprot.org/uniprotkb/{accession}",
            source_locator="UniProt features: Signal, Transmembrane, Topological domain",
            supporting_excerpt="; ".join(access["basis"])[:600],
            limitations=[access["caveat"]],
            origin="real",
        )
    ]

    assessment = {
        "accession": accession,
        "protein_name": protein["protein_name"],
        "organism": protein["organism"],
        "taxon_id": protein["taxon_id"],
        "length": protein["length"],
        "sequence_checksum": protein["sequence_checksum"],
        "accessibility": access,
        "experimental_structures": experimental[:8],
        "n_experimental_structures": len(experimental),
        "alphafold": af,
        "mechanism_constraint": mechanism_constraint,
        "limitations": limitations,
        "design_ready": bool(access["accessible"] and (experimental or af)),
    }
    return assessment, evidence, prov


def verify_identity(gene_symbol: str, ensembl_id: str, accession: str, expected_taxon: int = 9606) -> dict[str, Any]:
    """Confirm gene -> protein identity rather than joining on symbol alone."""
    protein, _ = fetch_protein(accession)
    problems: list[str] = []
    if protein["taxon_id"] != expected_taxon:
        problems.append(f"Accession {accession} is taxon {protein['taxon_id']}, expected {expected_taxon}.")
    r = _get(UNIPROT_ENTRY.format(acc=accession), params={"fields": "gene_names,xref_ensembl"})
    xrefs, genes = [], []
    if r.status_code == 200:
        d = r.json()
        genes = [(g.get("geneName") or {}).get("value") for g in (d.get("genes") or [])]
        for x in d.get("uniProtKBCrossReferences") or []:
            if x.get("database") == "Ensembl":
                for p in x.get("properties", []):
                    if p.get("key") == "GeneId":
                        xrefs.append(p["value"].split(".")[0])
    if genes and gene_symbol not in genes:
        problems.append(f"Symbol '{gene_symbol}' is not among UniProt gene names {genes}.")
    if xrefs and ensembl_id not in xrefs:
        problems.append(f"Ensembl gene {ensembl_id} is not cross-referenced by {accession} (has {sorted(set(xrefs))[:4]}).")
    return {
        "verified": not problems,
        "problems": problems,
        "uniprot_gene_names": genes,
        "ensembl_xrefs": sorted(set(xrefs))[:8],
        "taxon_id": protein["taxon_id"],
        "sequence_checksum": protein["sequence_checksum"],
        "length": protein["length"],
    }
