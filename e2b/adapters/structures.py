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

from ..contracts import Evidence, Provenance, utcnow

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


PDBE_SIFTS = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb}"


def pdb_uniprot_mapping(pdb_id: str) -> dict[str, Any]:
    """Map PDB author residue numbering onto UniProt canonical numbering, via SIFTS.

    This is not bookkeeping, it is a correctness requirement. Design tools address
    residues in **PDB author numbering**; sequence-level analysis (paralogue epitope
    cross-reactivity, conservation, mutagenesis design) addresses them in **UniProt
    canonical numbering**. For 1IAR the IL4R chain carries an offset of +25, so a hotspot
    written as residue 59 in the structure is residue 84 in the sequence. Passing one set
    of numbers into the other silently analyses a different part of the protein.

    It also answers a question that is easy to get wrong in a complex: *which chain is
    actually my target?* In 1IAR chain A is IL-4 (P05112) and chain B is IL-4Rα (P24394).
    Designing against "chain A" of that entry engages the cytokine, not the receptor.

    Returns ``{pdb_id, chains: {chain: {accession, unp_start, unp_end, author_start,
    offset, identifier}}, accessions: {accession: [chains]}}``.
    """
    r = _get(PDBE_SIFTS.format(pdb=pdb_id.lower()))
    if r.status_code != 200:
        raise StructureAdapterError(
            f"SIFTS mapping unavailable for {pdb_id} (HTTP {r.status_code}). Residue "
            "numbering cannot be reconciled; do not assume author numbering equals UniProt.",
            transient=r.status_code >= 500,
        )
    body = (r.json() or {}).get(pdb_id.lower()) or {}
    blocks = body.get("UniProt") or {}
    chains: dict[str, Any] = {}
    accessions: dict[str, list[str]] = {}
    for acc, blk in blocks.items():
        for m in blk.get("mappings") or []:
            chain = m.get("chain_id") or m.get("struct_asym_id")
            start, end = (m.get("start") or {}), (m.get("end") or {})
            unp_start = m.get("unp_start")
            # SIFTS gives author numbering for most entries, but for some (3L5X, 4I77,
            # 8K4Q among our candidates) `author_residue_number` is null and only the
            # SEQRES index is present. Those are DIFFERENT numbering bases and an offset
            # derived from one is wrong for the other, so the basis is recorded rather
            # than quietly substituted.
            author_start = start.get("author_residue_number")
            if author_start is not None:
                basis, ref_start, ref_end = "author", author_start, end.get("author_residue_number")
            else:
                basis, ref_start, ref_end = "seqres", start.get("residue_number"), end.get("residue_number")
            offset = (unp_start - ref_start) if (unp_start is not None and ref_start is not None) else None
            chains[chain] = {
                "accession": acc,
                "identifier": blk.get("identifier"),
                "unp_start": unp_start,
                "unp_end": m.get("unp_end"),
                "numbering_basis": basis,
                "ref_start": ref_start,
                "ref_end": ref_end,
                "author_start": author_start,
                "offset_to_uniprot": offset,
            }
            accessions.setdefault(acc, []).append(chain)

    # One accession can map to SEVERAL chains in an entry -- 8K4Q maps P24394 to both A and
    # C. `accessions` has always exposed that, but a caller has to go looking, and
    # chain_for_accession() returns a single chain. A signal a caller can miss is not a
    # signal, so multiplicity is surfaced as warnings and alternative_chains here.
    #
    # The offsets happening to agree between copies is a property of a given entry, not a
    # guarantee, so the two cases are reported differently: agreeing offsets are an
    # ambiguity about WHICH COPY to design against, disagreeing offsets additionally make
    # the numbering wrong if the wrong chain is used.
    multiplicity: dict[str, dict[str, Any]] = {}
    mapping_warnings: list[str] = []
    for acc, chs in accessions.items():
        if len(chs) < 2:
            continue
        offsets = {c: chains[c].get("offset_to_uniprot") for c in chs}
        distinct = {o for o in offsets.values() if o is not None}
        multiplicity[acc] = {"chains": sorted(chs), "offsets": offsets,
                             "offsets_agree": len(distinct) <= 1}
        if len(distinct) <= 1:
            mapping_warnings.append(
                f"{acc} maps to {len(chs)} chains in {pdb_id.upper()}: {sorted(chs)}. They share "
                f"offset {next(iter(distinct), None)}, so numbering is unaffected, but WHICH COPY "
                "to design against is ambiguous and should be chosen deliberately."
            )
        else:
            mapping_warnings.append(
                f"{acc} maps to {len(chs)} chains in {pdb_id.upper()} with DIFFERENT offsets "
                f"{offsets}. Residue numbering depends on which chain is used; picking the wrong "
                "one shifts every position."
            )
    for acc, chs in accessions.items():
        for c in chs:
            chains[c]["alternative_chains"] = [x for x in chs if x != c]
    if not chains:
        raise StructureAdapterError(f"SIFTS returned no UniProt mapping for {pdb_id}.")
    return {"pdb_id": pdb_id.upper(), "chains": chains, "accessions": accessions,
            "multiple_chains_per_accession": multiplicity,
            "warnings": mapping_warnings or None,
            "source": "PDBe SIFTS", "retrieved_at": utcnow().isoformat()}


def chain_for_accession(pdb_id: str, accession: str) -> tuple[str | None, dict[str, Any]]:
    """Which chain of this entry is the requested protein? Never guess 'A'.

    When the accession maps to several chains the FIRST is returned, but the mapping
    carries ``warnings`` and each chain carries ``alternative_chains``, so the caller
    cannot silently design against an arbitrary copy without the ambiguity being visible.
    """
    mapping = pdb_uniprot_mapping(pdb_id)
    chains = mapping["accessions"].get(accession) or []
    return (chains[0] if chains else None), mapping


def to_uniprot_numbering(pdb_id: str, chain: str, author_residues: "list[int | str]") -> dict[str, Any]:
    """Convert author-numbered residues (``59`` or ``"B59"``) to UniProt canonical numbering."""
    mapping = pdb_uniprot_mapping(pdb_id)
    info = mapping["chains"].get(chain)
    if info is None:
        raise StructureAdapterError(
            f"Chain '{chain}' has no UniProt mapping in {pdb_id}. Available: {sorted(mapping['chains'])}."
        )
    offset = info["offset_to_uniprot"]
    if offset is None:
        raise StructureAdapterError(f"SIFTS gave no usable offset for {pdb_id} chain {chain}.")
    converted, out_of_range = [], []
    for res in author_residues:
        n = int(str(res).lstrip(chain)) if isinstance(res, str) else int(res)
        u = n + offset
        (converted if (info["unp_start"] <= u <= info["unp_end"]) else out_of_range).append(u)

    caveats = []
    if out_of_range:
        caveats.append("Residues outside the mapped span are not covered by this structure and must "
                       "not be treated as verified positions.")
    if info["numbering_basis"] != "author":
        caveats.append(f"SIFTS exposes no author numbering for this entry; the offset is relative to "
                       f"{info['numbering_basis']} numbering. Confirm which basis the input residues use "
                       f"before relying on the converted positions.")
    return {
        "pdb_id": mapping["pdb_id"], "chain": chain, "accession": info["accession"],
        "numbering_basis": info["numbering_basis"],
        "offset_to_uniprot": offset,
        "uniprot_residues": converted,
        "out_of_mapped_range": out_of_range,
        "mapped_span_uniprot": [info["unp_start"], info["unp_end"]],
        "caveats": caveats or None,
        "source": "PDBe SIFTS",
    }


RCSB_PDB = "https://files.rcsb.org/download/{pdb}.pdb"

# Heavy-atom contact distance. 5.0 A is the usual convention for a residue-level protein
# interface; it is a convention, not a measurement, and is reported with the result.
CONTACT_CUTOFF_A = 5.0


def _parse_pdb_atoms(text: str) -> dict[str, list[tuple[int, str, str, float, float, float]]]:
    """Minimal PDB ATOM parser -> {chain: [(resseq, icode, resname, x, y, z), ...]}.

    Heavy atoms only, first altloc only. Deliberately dependency-free: this runs in the
    same process as the agent and should not drag in a structure-parsing stack.
    """
    out: dict[str, list] = {}
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        element = line[76:78].strip()
        if element == "H":
            continue
        altloc = line[16]
        if altloc not in (" ", "A"):
            continue
        chain = line[21]
        try:
            resseq = int(line[22:26])
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        out.setdefault(chain, []).append((resseq, line[26].strip(), line[17:20].strip(), x, y, z))
    return out


def interface_residues(
    pdb_id: str,
    target_chain: str,
    partner_chains: "list[str] | None" = None,
    *,
    cutoff_a: float = CONTACT_CUTOFF_A,
    top_n: int | None = None,
) -> dict[str, Any]:
    """Derive the natural binding interface on ``target_chain`` from the deposited complex.

    This is how hotspot residues should be chosen: **measured from coordinates**, not
    recalled. For a target whose natural partner is in the entry, the residues its partner
    actually contacts are the epitope a competitive binder must cover.

    Returns residues in BOTH numbering systems -- author numbering for the design tools,
    UniProt canonical for sequence-level analysis -- because those differ and the design
    lane rejects a request that does not say which it is using.

    Contact count is a proxy for interface involvement, not a measured binding energy. A
    residue that contacts the partner is not thereby a hotspot in the thermodynamic sense;
    calling these "hotspots" follows the design-tool convention of "residues to engage".
    """
    r = _get(RCSB_PDB.format(pdb=pdb_id.upper()))
    if r.status_code != 200:
        raise StructureAdapterError(
            f"Could not fetch coordinates for {pdb_id} (HTTP {r.status_code}).",
            transient=r.status_code >= 500,
        )
    chains = _parse_pdb_atoms(r.text)
    if target_chain not in chains:
        raise StructureAdapterError(
            f"Chain '{target_chain}' has no ATOM records in {pdb_id}. Present: {sorted(chains)}."
        )
    partners = partner_chains or [c for c in chains if c != target_chain]
    partners = [c for c in partners if c in chains]
    if not partners:
        raise StructureAdapterError(
            f"{pdb_id} has no partner chain besides '{target_chain}', so no natural interface "
            "can be derived from it. Choose an entry containing the complex."
        )

    import math

    tgt = chains[target_chain]
    cut2 = cutoff_a * cutoff_a
    counts: dict[tuple[int, str], dict[str, Any]] = {}
    # Bucket partner atoms into a coarse grid so this stays linear-ish rather than N*M.
    grid: dict[tuple[int, int, int], list[tuple[float, float, float, str]]] = {}
    cs = cutoff_a
    for pc in partners:
        for (_, _, _, x, y, z) in chains[pc]:
            key = (int(x // cs), int(y // cs), int(z // cs))
            grid.setdefault(key, []).append((x, y, z, pc))

    for (resseq, icode, resname, x, y, z) in tgt:
        gx, gy, gz = int(x // cs), int(y // cs), int(z // cs)
        near = 0
        hit_partners: set[str] = set()
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for (px, py, pz, pc) in grid.get((gx + dx, gy + dy, gz + dz), ()):
                        if (px - x) ** 2 + (py - y) ** 2 + (pz - z) ** 2 <= cut2:
                            near += 1
                            hit_partners.add(pc)
        if near:
            rec = counts.setdefault((resseq, icode), {
                "author_residue": resseq, "insertion_code": icode or None,
                "residue_name": resname, "atom_contacts": 0, "partner_chains": set(),
            })
            rec["atom_contacts"] += near
            rec["partner_chains"] |= hit_partners

    # An interface is only meaningful RELATIVE TO A PARTNER. Identify every partner before
    # reporting residues, because three failure modes look identical in the coordinates:
    #   - the partner is the natural ligand            -> this is the epitope we want
    #   - the partner is another copy of the target    -> homotypic/crystal contact
    #   - the partner has no UniProt mapping           -> usually an antibody or nanobody,
    #                                                     so the "interface" is that
    #                                                     binder's footprint, not the
    #                                                     receptor-binding face
    partner_identity: dict[str, str | None] = {}
    warnings: list[str] = []
    try:
        mp = pdb_uniprot_mapping(pdb_id)
        self_acc = (mp["chains"].get(target_chain) or {}).get("accession")
        for pc in partners:
            info = mp["chains"].get(pc)
            partner_identity[pc] = (info or {}).get("identifier")
            if info is None:
                warnings.append(
                    f"Chain {pc} has no UniProt mapping. Unmapped polymer chains in a complex are "
                    "commonly antibody or nanobody fragments; if so, residues contacting it are that "
                    "binder's epitope, NOT the natural ligand-binding site."
                )
            elif self_acc and info.get("accession") == self_acc:
                warnings.append(
                    f"Chain {pc} is another copy of the target itself ({info.get('identifier')}). "
                    "Contacts with it are homotypic or crystal-packing, not a ligand interface."
                )
    except StructureAdapterError as exc:
        warnings.append(f"Partner identity could not be established: {exc}")

    biologically_relevant = bool(partner_identity) and not warnings

    ranked = sorted(counts.values(), key=lambda d: (-d["atom_contacts"], d["author_residue"]))
    if top_n:
        ranked = ranked[:top_n]
    for rec in ranked:
        rec["partner_chains"] = sorted(rec["partner_chains"])

    author_nums = [rec["author_residue"] for rec in ranked]
    numbering: dict[str, Any] = {}
    try:
        numbering = to_uniprot_numbering(pdb_id, target_chain, author_nums)
        for rec, unp in zip(ranked, numbering["uniprot_residues"] + [None] * len(ranked)):
            rec["uniprot_residue"] = unp
    except StructureAdapterError as exc:
        numbering = {"error": str(exc)}

    return {
        "pdb_id": pdb_id.upper(),
        "target_chain": target_chain,
        "partner_chains": partners,
        "partner_identity": partner_identity,
        "partner_is_biologically_relevant": biologically_relevant,
        "warnings": warnings or None,
        "cutoff_a": cutoff_a,
        "n_interface_residues": len(ranked),
        "residues": ranked,
        "hotspot_residues_author": [f"{target_chain}{rec['author_residue']}" for rec in ranked],
        "hotspot_residues_uniprot": [rec.get("uniprot_residue") for rec in ranked],
        "numbering": numbering,
        "method": (
            f"Heavy-atom contacts within {cutoff_a} A between chain {target_chain} and "
            f"{partners} in deposited entry {pdb_id.upper()}, ranked by contact count. "
            "Contact count indicates interface involvement; it is NOT a measured binding "
            "energy and these are not thermodynamic hotspots."
        ),
        "source": "RCSB PDB coordinates + PDBe SIFTS numbering",
    }


def transfer_epitope(
    from_pdb: str, from_chain: str, author_residues: "list[int | str]",
    to_pdb: str, to_chain: str,
) -> dict[str, Any]:
    """Re-express an epitope defined on one entry in another entry's author numbering.

    Why this exists: the entry that **defines** the right epitope and the entry that gives
    the best **coordinates to design against** are often not the same one.

    For IL-4Rα, 1IAR is the only entry containing the natural ligand (IL-4), so it is the
    only entry from which the cytokine-binding face can be derived -- but its receptor
    chain has two unresolved breaks. 8K4Q covers the same span with none, yet its partner
    chains are another IL-4Rα copy and unmapped chains, so an interface derived from it
    is not the ligand site. The correct move is to take the epitope from 1IAR and the
    coordinates from 8K4Q, with UniProt canonical numbering as the bridge.

    Both entries must map to the same UniProt accession; this refuses otherwise rather
    than transferring residues between different proteins.
    """
    src = to_uniprot_numbering(from_pdb, from_chain, author_residues)
    dst_map = pdb_uniprot_mapping(to_pdb)
    dst = dst_map["chains"].get(to_chain)
    if dst is None:
        raise StructureAdapterError(
            f"Chain '{to_chain}' has no UniProt mapping in {to_pdb}. "
            f"Available: {sorted(dst_map['chains'])}."
        )
    if dst["accession"] != src["accession"]:
        raise StructureAdapterError(
            f"Refusing to transfer: {from_pdb}/{from_chain} is {src['accession']} but "
            f"{to_pdb}/{to_chain} is {dst['accession']}. These are different proteins."
        )
    offset = dst["offset_to_uniprot"]
    if offset is None:
        raise StructureAdapterError(f"No usable offset for {to_pdb} chain {to_chain}.")

    transferred, out_of_span = [], []
    for unp in src["uniprot_residues"]:
        if dst["unp_start"] <= unp <= dst["unp_end"]:
            transferred.append({"uniprot_residue": unp, "author_residue": unp - offset})
        else:
            out_of_span.append(unp)

    return {
        "accession": src["accession"],
        "from": {"pdb_id": src["pdb_id"], "chain": from_chain, "numbering_basis": src["numbering_basis"]},
        "to": {"pdb_id": dst_map["pdb_id"], "chain": to_chain,
               "numbering_basis": dst["numbering_basis"], "offset_to_uniprot": offset,
               "mapped_span_uniprot": [dst["unp_start"], dst["unp_end"]]},
        "uniprot_residues": [t["uniprot_residue"] for t in transferred],
        "hotspot_residues_author": [f"{to_chain}{t['author_residue']}" for t in transferred],
        "not_covered_by_target_entry": out_of_span,
        "caveat": (
            f"{len(out_of_span)} residue(s) of the source epitope are outside the span "
            f"{to_pdb} chain {to_chain} resolves and were dropped: {out_of_span}. The "
            "transferred epitope is therefore partial."
        ) if out_of_span else None,
        "source": "PDBe SIFTS numbering bridge",
    }


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
