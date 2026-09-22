"""Paralogue cross-reactivity RISK for the selected target.

The question this answers: *if the design engages this epitope, which other human
proteins present a similar surface there?* An epitope conserved across close family
members is a specificity risk; an epitope that diverges is a specificity opportunity
and is the thing worth telling a designer about.

What this is NOT: a measurement of cross-reactivity. Nothing here is binding data.
Sequence identity over an epitope region is a *risk indicator* computed from
sequences; a divergent epitope can still be cross-bound (shape and chemistry, not
identity, determine recognition) and a conserved epitope may not be, if the
paralogue is not co-expressed or is not accessible. Both directions of error are
real and are stated in the output.

Family retrieval uses up to three independent, named routes and reports each
route's status separately, so that "no paralogues found" is distinguishable from
"retrieval failed":

* ``ensembl_compara``  gene-level paralogues from Ensembl Compara (the canonical
  definition of a paralogue), mapped to reviewed human UniProt accessions.
* ``uniref50``         reviewed human members of the target's UniRef50 cluster
  (>= ~50% sequence identity to the cluster representative).
* ``pfam_family``      reviewed human proteins sharing a Pfam domain with the target.

Sequences are pulled through the existing structures adapter
(``e2b.adapters.structures.fetch_protein``) so that every sequence used here is the
same real UniProt record the rest of the pipeline uses, with its own provenance.
"""

from __future__ import annotations

import time
from typing import Any

import requests

from ..adapters.structures import StructureAdapterError, fetch_protein
from .align import align_global, region_identity

PARALOG_VERSION = "paralogs/1.1.0"

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIREF_SEARCH = "https://rest.uniprot.org/uniref/search"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{acc}.json"
ENSEMBL_HOMOLOGY = "https://rest.ensembl.org/homology/id/human/{gene}"

HTTP_TIMEOUT_S = 45
ENSEMBL_RETRIES = 3
ENSEMBL_BACKOFF_S = 2.0
DEFAULT_MAX_MEMBERS = 12

# Epitope-identity risk bands. Conventions, stated in docs/analysis_methods.md.
RISK_BANDS = [(70.0, "high"), (40.0, "moderate"), (0.0, "low")]


def _risk_band(pct: float | None) -> str | None:
    if pct is None:
        return None
    for cut, band in RISK_BANDS:
        if pct >= cut:
            return band
    return "low"


def _get(url: str, **kw) -> requests.Response:
    return requests.get(url, timeout=HTTP_TIMEOUT_S, **kw)


# ---------------------------------------------------------------------------------
# Family retrieval routes. Each returns (accessions, status_dict).
# ---------------------------------------------------------------------------------


def _route_ensembl_compara(ensembl_gene_id: str | None) -> tuple[list[str], dict[str, Any]]:
    st: dict[str, Any] = {"route": "ensembl_compara", "status": "skipped", "n_found": 0, "detail": None}
    if not ensembl_gene_id:
        st["detail"] = "no Ensembl gene id supplied for the target"
        return [], st
    t0 = time.time()
    r = None
    attempts = 0
    # Ensembl REST returns 5xx under load; a 5xx is a transient source failure, so it
    # is retried a bounded number of times and then REPORTED as a failure -- never
    # silently converted into "this gene has no paralogues".
    for attempt in range(ENSEMBL_RETRIES):
        attempts = attempt + 1
        try:
            r = _get(
                ENSEMBL_HOMOLOGY.format(gene=ensembl_gene_id),
                params={"type": "paralogues", "format": "condensed"},
                headers={"Accept": "application/json"},
            )
        except requests.RequestException as exc:
            st.update(status="source_failure", detail=f"transport error: {exc}", attempts=attempts)
            if attempt == ENSEMBL_RETRIES - 1:
                return [], st
            time.sleep(ENSEMBL_BACKOFF_S * (attempt + 1))
            continue
        if r.status_code < 500:
            break
        if attempt < ENSEMBL_RETRIES - 1:
            time.sleep(ENSEMBL_BACKOFF_S * (attempt + 1))
    st["attempts"] = attempts
    if r is None:
        st.update(status="source_failure", detail="no response from Ensembl")
        return [], st
    st["elapsed_ms"] = int((time.time() - t0) * 1000)
    st["endpoint"] = ENSEMBL_HOMOLOGY.format(gene=ensembl_gene_id)
    if r.status_code != 200:
        st.update(status="source_failure", detail=f"Ensembl HTTP {r.status_code}")
        return [], st
    try:
        homologies = (r.json().get("data") or [{}])[0].get("homologies") or []
    except (ValueError, IndexError, AttributeError) as exc:
        st.update(status="source_failure", detail=f"unparseable Ensembl payload: {exc}")
        return [], st
    genes = [h.get("id") for h in homologies if h.get("id")]
    if not genes:
        st.update(status="ok", n_found=0, detail="Ensembl Compara reports no paralogues for this gene")
        return [], st
    accs, mapst = _ensembl_genes_to_accessions(genes)
    st.update(status="ok", n_found=len(accs), detail=f"{len(genes)} paralogous genes, {len(accs)} mapped to "
                                                     f"reviewed human accessions", mapping=mapst)
    return accs, st


def _ensembl_genes_to_accessions(gene_ids: list[str], batch: int = 25) -> tuple[list[str], dict[str, Any]]:
    accs: list[str] = []
    errors: list[str] = []
    for i in range(0, len(gene_ids), batch):
        chunk = gene_ids[i : i + batch]
        q = " OR ".join(f"xref:{g}" for g in chunk)
        try:
            r = _get(
                UNIPROT_SEARCH,
                params={
                    "query": f"({q}) AND (reviewed:true) AND (organism_id:9606)",
                    "fields": "accession",
                    "size": 100,
                },
            )
            if r.status_code != 200:
                errors.append(f"HTTP {r.status_code}")
                continue
            accs.extend(x["primaryAccession"] for x in r.json().get("results", []))
        except requests.RequestException as exc:
            errors.append(str(exc))
    return sorted(set(accs)), {"errors": errors}


def _route_uniref50(accession: str) -> tuple[list[str], dict[str, Any]]:
    st: dict[str, Any] = {"route": "uniref50", "status": "skipped", "n_found": 0, "detail": None}
    t0 = time.time()
    try:
        # Resolve the cluster that CONTAINS the accession first: keying the cluster id
        # off the accession only works when the target happens to be the cluster
        # representative, which silently loses family members otherwise.
        rc = _get(
            UNIREF_SEARCH,
            params={"query": f"uniprot_id:{accession} AND identity:0.5", "fields": "id", "size": 5},
        )
        if rc.status_code != 200:
            st.update(status="source_failure", detail=f"UniRef HTTP {rc.status_code} resolving cluster")
            return [], st
        clusters = [x["id"] for x in rc.json().get("results", [])]
        if not clusters:
            st.update(status="ok", n_found=0, detail="accession is in no UniRef50 cluster")
            return [], st
        q = " OR ".join(f"uniref_cluster_50:{c}" for c in clusters)
        r = _get(
            UNIPROT_SEARCH,
            params={
                "query": f"({q}) AND (organism_id:9606) AND (reviewed:true)",
                "fields": "accession",
                "size": 100,
            },
        )
    except requests.RequestException as exc:
        st.update(status="source_failure", detail=f"transport error: {exc}")
        return [], st
    st["elapsed_ms"] = int((time.time() - t0) * 1000)
    if r.status_code != 200:
        st.update(status="source_failure", detail=f"UniProt HTTP {r.status_code}")
        return [], st
    accs = sorted({x["primaryAccession"] for x in r.json().get("results", [])} - {accession})
    st.update(
        status="ok",
        n_found=len(accs),
        clusters=clusters,
        detail=f"reviewed human members of UniRef50 cluster(s) {clusters} containing the target "
               "(>= ~50% identity to the cluster representative)",
    )
    return accs, st


def _route_pfam_family(accession: str) -> tuple[list[str], dict[str, Any]]:
    st: dict[str, Any] = {"route": "pfam_family", "status": "skipped", "n_found": 0, "detail": None}
    t0 = time.time()
    try:
        r = _get(UNIPROT_ENTRY.format(acc=accession), params={"fields": "xref_pfam"})
        if r.status_code != 200:
            st.update(status="source_failure", detail=f"UniProt HTTP {r.status_code} fetching Pfam xrefs")
            return [], st
        pfams = [x["id"] for x in r.json().get("uniProtKBCrossReferences", []) if x.get("database") == "Pfam"]
        if not pfams:
            st.update(status="ok", n_found=0, detail="target has no Pfam cross-reference")
            return [], st
        q = " OR ".join(f"xref:pfam-{p}" for p in pfams)
        r2 = _get(
            UNIPROT_SEARCH,
            params={
                "query": f"({q}) AND (organism_id:9606) AND (reviewed:true)",
                "fields": "accession",
                "size": 100,
            },
        )
        if r2.status_code != 200:
            st.update(status="source_failure", detail=f"UniProt HTTP {r2.status_code} on family search")
            return [], st
        accs = sorted({x["primaryAccession"] for x in r2.json().get("results", [])} - {accession})
    except requests.RequestException as exc:
        st.update(status="source_failure", detail=f"transport error: {exc}")
        return [], st
    st["elapsed_ms"] = int((time.time() - t0) * 1000)
    st.update(status="ok", n_found=len(accs), detail=f"reviewed human proteins sharing Pfam {pfams}",
              pfam_ids=pfams)
    return accs, st


def retrieve_family(
    accession: str,
    ensembl_gene_id: str | None = None,
    routes: tuple[str, ...] = ("ensembl_compara", "uniref50", "pfam_family"),
) -> dict[str, Any]:
    """Union of family members across the named routes, with per-route status."""
    per_route: list[dict[str, Any]] = []
    found: dict[str, list[str]] = {}
    for route in routes:
        if route == "ensembl_compara":
            accs, st = _route_ensembl_compara(ensembl_gene_id)
        elif route == "uniref50":
            accs, st = _route_uniref50(accession)
        elif route == "pfam_family":
            accs, st = _route_pfam_family(accession)
        else:
            accs, st = [], {"route": route, "status": "skipped", "detail": "unknown route name", "n_found": 0}
        per_route.append(st)
        for a in accs:
            found.setdefault(a, []).append(route)
    found.pop(accession, None)
    any_ok = any(s["status"] == "ok" for s in per_route)
    all_failed = per_route and all(s["status"] == "source_failure" for s in per_route)
    return {
        "target_accession": accession,
        "routes": per_route,
        "members": [{"accession": a, "routes": rs} for a, rs in sorted(found.items())],
        "n_members": len(found),
        "retrieval_ok": bool(any_ok),
        "all_routes_failed": bool(all_failed),
    }


# ---------------------------------------------------------------------------------
# Epitope-level comparison
# ---------------------------------------------------------------------------------


def paralog_cross_reactivity_risk(
    target_accession: str,
    epitope_residues: list[int] | None,
    ensembl_gene_id: str | None = None,
    max_members: int = DEFAULT_MAX_MEMBERS,
    routes: tuple[str, ...] = ("ensembl_compara", "uniref50", "pfam_family"),
    target_sequence: str | None = None,
) -> dict[str, Any]:
    """Retrieve the target's close family, align each member, and report identity
    over the engaged epitope region specifically.

    `epitope_residues` are 1-based positions in the target's canonical UniProt
    sequence -- i.e. ``contracts.BindingRegion.residues`` (or ``hotspot_residues``)
    already mapped to UniProt numbering by the caller. If they are None or empty,
    the function still returns full-length identity but sets
    ``epitope_status="unspecified"`` and every epitope figure to None; the triage
    layer then treats cross-reactivity risk as UNKNOWN, which is conservative.

    All returned numbers are computed from real UniProt sequences.
    """
    out: dict[str, Any] = {
        "method_version": PARALOG_VERSION,
        "is_measurement": False,
        "origin": "real",
        "target_accession": target_accession,
        "epitope_residues": sorted(set(int(x) for x in epitope_residues)) if epitope_residues else [],
        "epitope_status": "specified" if epitope_residues else "unspecified",
        "status": "ok",
        "failures": [],
        "family_retrieval": None,
        "comparisons": [],
        "summary": {
            "n_members_compared": 0,
            "max_epitope_identity_pct": None,
            "max_epitope_identity_member": None,
            "risk_band": None,
            "risk_basis": None,
        },
        "limitations": [
            "Sequence identity over an epitope is a RISK INDICATOR, not a measured cross-reactivity.",
            "Recognition depends on shape, chemistry and conformation; a divergent epitope can still "
            "be cross-bound and a conserved one may not be.",
            "Identity is read from a global pairwise alignment (BLOSUM62, affine gaps), which is a "
            "conservative estimate for distantly related pairs.",
            "Only reviewed human proteins retrieved by the named routes are considered; this is not a "
            "proteome-wide off-target screen and says nothing about non-human species used in tox.",
            "Co-expression of the paralogue in the relevant tissue is not assessed here.",
        ],
    }

    # Target sequence
    tgt_prov = None
    if target_sequence is None:
        try:
            info, tgt_prov = fetch_protein(target_accession)
            target_sequence = info["sequence"]
            out["target_length"] = info["length"]
            out["target_gene"] = info.get("protein_name")
        except StructureAdapterError as exc:
            out["status"] = "source_failure"
            out["failures"].append(
                {"stage": "target_sequence", "error": str(exc), "transient": getattr(exc, "transient", False),
                 "note": "retrieval failure, NOT evidence that the target has no family"}
            )
            return out
    else:
        out["target_length"] = len(target_sequence)
    out["target_provenance"] = _prov_dict(tgt_prov)

    if epitope_residues:
        bad = [p for p in out["epitope_residues"] if p < 1 or p > len(target_sequence)]
        if bad:
            out["failures"].append(
                {"stage": "epitope_mapping", "error": f"epitope positions outside the target sequence: {bad}",
                 "note": "positions must be 1-based indices into the canonical UniProt sequence"}
            )
            out["epitope_status"] = "invalid"

    fam = retrieve_family(target_accession, ensembl_gene_id=ensembl_gene_id, routes=routes)
    out["family_retrieval"] = fam
    if fam["all_routes_failed"]:
        out["status"] = "source_failure"
        out["failures"].append(
            {"stage": "family_retrieval", "error": "every family-retrieval route failed",
             "note": "a retrieval failure is a source failure, not evidence that no paralogues exist"}
        )
        return out
    if fam["n_members"] == 0:
        out["status"] = "no_family_members_found"
        out["summary"]["risk_basis"] = (
            "no reviewed human family member was returned by any route that succeeded; "
            "this is an absence reported by the sources queried, not a proteome-wide guarantee"
        )
        return out

    members, prefilter = _prioritise_members(target_sequence, fam["members"], max_members)
    out["member_prefilter"] = prefilter

    for mem in members:
        acc = mem["accession"]
        try:
            info, prov = fetch_protein(acc)
        except StructureAdapterError as exc:
            out["failures"].append({"stage": "member_sequence", "accession": acc, "error": str(exc)})
            continue
        seq = info["sequence"]
        aln = align_global(target_sequence, seq)
        row: dict[str, Any] = {
            "accession": acc,
            "protein_name": info.get("protein_name"),
            "entry_name": info.get("entry_name"),
            "length": info.get("length"),
            "retrieval_routes": mem["routes"],
            "full_length_percent_identity": round(aln.percent_identity, 2) if aln.percent_identity is not None else None,
            "alignment_score_blosum62": aln.score,
            "epitope": None,
            "provenance": _prov_dict(prov),
            "origin": "real",
        }
        if out["epitope_status"] == "specified":
            reg = region_identity(aln, out["epitope_residues"])
            reg["percent_identity_region"] = (
                round(reg["percent_identity_region"], 2) if reg["percent_identity_region"] is not None else None
            )
            reg["risk_band"] = _risk_band(reg["percent_identity_region"])
            row["epitope"] = reg
        out["comparisons"].append(row)

    out["summary"]["n_members_compared"] = len(out["comparisons"])
    if out["epitope_status"] == "specified":
        scored = [(c["epitope"]["percent_identity_region"], c) for c in out["comparisons"]
                  if c["epitope"] and c["epitope"]["percent_identity_region"] is not None]
        if scored:
            best = max(scored, key=lambda t: t[0])
            out["summary"]["max_epitope_identity_pct"] = best[0]
            out["summary"]["max_epitope_identity_member"] = best[1]["accession"]
            out["summary"]["risk_band"] = _risk_band(best[0])
            n_al = best[1]["epitope"]["n_aligned"]
            n_req = len(out["epitope_residues"])
            out["summary"]["risk_basis"] = (
                f"highest epitope-region identity among {len(scored)} compared family member(s) is "
                f"{best[0]:.1f}% ({best[1]['accession']}, {n_al}/{n_req} epitope positions aligned)"
            )
            out["summary"]["epitope_coverage_fraction"] = round(n_al / n_req, 3) if n_req else None
            if n_req and n_al / n_req < 0.6:
                out["summary"]["coverage_caveat"] = (
                    f"only {n_al} of {n_req} epitope positions aligned to this member, so the "
                    "percentage rests on a minority of the epitope and should be read as weak"
                )
        else:
            out["summary"]["risk_basis"] = (
                "no family member had any epitope position alignable; epitope-level risk is UNKNOWN"
            )
    else:
        out["summary"]["risk_basis"] = (
            "no epitope region was supplied, so epitope-level identity could not be computed; "
            "cross-reactivity risk is UNKNOWN and must be treated conservatively"
        )
        maxfl = [c["full_length_percent_identity"] for c in out["comparisons"]
                 if c["full_length_percent_identity"] is not None]
        if maxfl:
            out["summary"]["max_full_length_identity_pct"] = max(maxfl)
    return out


KMER_K = 3


def _kmers(seq: str, k: int = KMER_K) -> set[str]:
    return {seq[i : i + k] for i in range(len(seq) - k + 1)}


def fetch_sequences_bulk(accessions: list[str], batch: int = 100) -> tuple[dict[str, str], list[str]]:
    """One FASTA stream call per batch, used only to PRE-RANK family members.

    The sequences finally compared are re-fetched individually through the
    structures adapter so that each carries its own provenance and features.
    """
    seqs: dict[str, str] = {}
    errors: list[str] = []
    for i in range(0, len(accessions), batch):
        chunk = accessions[i : i + batch]
        q = " OR ".join(f"accession:{a}" for a in chunk)
        try:
            r = _get(
                "https://rest.uniprot.org/uniprotkb/stream",
                params={"query": f"({q})", "format": "fasta"},
            )
            if r.status_code != 200:
                errors.append(f"HTTP {r.status_code}")
                continue
            acc = None
            buf: list[str] = []
            for line in r.text.splitlines():
                if line.startswith(">"):
                    if acc:
                        seqs[acc] = "".join(buf)
                    parts = line[1:].split("|")
                    acc = parts[1] if len(parts) > 2 else parts[0]
                    buf = []
                else:
                    buf.append(line.strip())
            if acc:
                seqs[acc] = "".join(buf)
        except requests.RequestException as exc:
            errors.append(str(exc))
    return seqs, errors


def _prioritise_members(
    target_sequence: str, members: list[dict[str, Any]], max_members: int
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Choose which family members to align, deterministically and by closeness.

    A family route can return a hundred proteins that merely share a domain (every
    human tyrosine kinase shares PF07714 with EGFR). Aligning the alphabetically
    first N of those would compare the target to an arbitrary subset. Members are
    therefore pre-ranked by 3-mer containment similarity against the target -- a
    cheap, deterministic, order-independent proxy for sequence relatedness -- and
    the top `max_members` are fully aligned. The pre-ranking is a FILTER, not a
    reported quantity; every identity figure in the output comes from the full
    Needleman-Wunsch alignment.

    Ties and members whose sequence could not be pre-fetched fall back to: number
    of retrieval routes (descending), Ensembl Compara membership first, then
    accession, so the selection is reproducible.
    """
    info: dict[str, Any] = {
        "method": f"{KMER_K}-mer containment similarity against the target sequence, "
                  "descending; ties broken by number of retrieval routes, Compara membership, accession",
        "purpose": "select which family members to align when the family is larger than max_members",
        "is_reported_quantity": False,
        "n_total_members": len(members),
        "n_selected": min(max_members, len(members)),
        "n_prefetched": 0,
        "errors": [],
    }
    if len(members) <= max_members:
        info["method"] = "no pre-filter applied: the whole family fits within max_members"
        return list(members), info

    accs = [m["accession"] for m in members]
    seqs, errors = fetch_sequences_bulk(accs)
    info["n_prefetched"] = len(seqs)
    info["errors"] = errors
    tk = _kmers(target_sequence)

    def sim(acc: str) -> float:
        s = seqs.get(acc)
        if not s:
            return -1.0
        mk = _kmers(s)
        if not mk or not tk:
            return -1.0
        return len(tk & mk) / len(tk | mk)

    scored = sorted(
        members,
        key=lambda m: (
            -round(sim(m["accession"]), 6),
            -len(m["routes"]),
            0 if "ensembl_compara" in m["routes"] else 1,
            m["accession"],
        ),
    )
    selected = scored[:max_members]
    info["selected"] = [
        {"accession": m["accession"], "prefilter_similarity": round(sim(m["accession"]), 4), "routes": m["routes"]}
        for m in selected
    ]
    info["n_not_compared"] = len(members) - len(selected)
    info["note"] = (
        f"{len(members)} family members were retrieved and {len(selected)} were aligned. The "
        f"{len(members) - len(selected)} not aligned are reported in family_retrieval and are NOT "
        "evidence of low cross-reactivity risk."
    )
    return selected, info


def _prov_dict(prov: Any) -> dict[str, Any] | None:
    if prov is None:
        return None
    try:
        d = prov.model_dump(mode="json")
    except AttributeError:
        return None
    return d
