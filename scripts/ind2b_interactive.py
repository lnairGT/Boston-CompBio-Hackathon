"""Render one indication2binder run directory as a single interactive HTML file.

Self-contained by construction: CSS, JavaScript and every chart are inlined,
charts are hand-emitted SVG rather than raster images, and no asset is fetched
at view time. The file opens from a laptop with no network.

Interaction is deliberately plain DOM work - sortable and filterable tables,
expandable detail rows, hover tooltips, copy buttons, scroll-spy navigation -
so there is no framework to break and nothing to load from a CDN.

Honesty rules, inherited from the static renderer and enforced here too:

* a layer that is absent is reported **absent**, never illustrated with a
  placeholder or a plausible-looking value;
* a null renders as an em dash, never as zero;
* a target or epitope that was *excluded* stays in the table with its reason,
  because a filtered-out row is a recorded decision rather than a missing one;
* the report states that no binder has been designed unless a design run
  record exists and reached ``succeeded``. Emitting a specification is not
  designing a binder, and the report says so.

Usage::

    python scripts/ind2b_interactive.py --run-dir runs/MONDO_0005233 \
        --out nsclc_report.html

Only ``pandas`` and the standard library are required.
"""

from __future__ import annotations

import argparse
import ast
import html
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import pandas as pd

SCHEMA_VERSION = 1

# Layer -> candidate filenames, shallowest match wins. Kept explicit so the
# script runs standalone, outside the skill and outside the repo.
LAYERS: dict[str, list[str]] = {
    "disease": ["stage0_disease.json", "disease.json"],
    "evidence": ["stage1_evidence.json"],
    "ranked": ["stage2_ranked_targets.csv", "ranked_targets.csv"],
    "scoring_meta": ["stage2_scoring_meta.json", "scoring_meta.json"],
    "complexes": ["stage3_complexes.json", "complexes.json"],
    "epitopes": ["stage4_epitopes.csv", "epitopes.csv"],
    "residues": ["stage4_interface_residues.csv", "interface_residues.csv"],
    "interface_meta": ["stage4_interface_meta.json"],
    "specs_manifest": ["stage5_specs/manifest.json"],
    "binders": ["stage6_binders.json"],
    "run_manifest": ["run_manifest.json", "provenance.json"],
    "design_run": ["design_run.json"],
    "candidates": ["candidates.json", "candidates.csv"],
}

PALETTE = [
    "#1f6feb", "#8250df", "#1a7f37", "#bf8700",
    "#cf222e", "#0969da", "#6e7781", "#953800",
]


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------

def load_run(run_dir: Path) -> dict[str, Any]:
    """Read every layer present; absent layers are simply missing from the map."""
    out: dict[str, Any] = {"_run_dir": run_dir, "_files": []}
    for p in sorted(run_dir.rglob("*")):
        if p.is_file():
            out["_files"].append(str(p.relative_to(run_dir)))

    for layer, names in LAYERS.items():
        for name in names:
            path = run_dir / name
            if not path.exists():
                continue
            try:
                if path.suffix == ".json":
                    out[layer] = json.loads(path.read_text())
                else:
                    out[layer] = pd.read_csv(path)
            except Exception as exc:  # noqa: BLE001
                out[layer] = {"_error": f"{type(exc).__name__}: {exc}"}
            out.setdefault("_paths", {})[layer] = str(path.relative_to(run_dir))
            break

    specs_dir = run_dir / "stage5_specs"
    if specs_dir.is_dir():
        out["_specs_dir"] = specs_dir
    return out


def missing_layers(run: dict) -> list[str]:
    return sorted(k for k in LAYERS if k not in run)


# --------------------------------------------------------------------------
# formatting helpers
# --------------------------------------------------------------------------

def esc(v: Any) -> str:
    if v is None:
        return "&mdash;"
    if isinstance(v, float) and math.isnan(v):
        return "&mdash;"
    return html.escape(str(v), quote=True)


def num(v: Any, nd: int = 2) -> str:
    """Format a number; a missing value is an em dash, never a zero."""
    if v is None:
        return "&mdash;"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return esc(v)
    if math.isnan(f):
        return "&mdash;"
    return f"{f:,.{nd}f}".rstrip("0").rstrip(".") if nd else f"{f:,.0f}"


def pct(v: Any) -> str:
    if v is None:
        return "&mdash;"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return esc(v)
    if math.isnan(f):
        return "&mdash;"
    return f"{f * 100:.0f}%"


def trunc(value: Any, n: int) -> str:
    """Escape and shorten a text cell.

    Written as its own helper because ``str(value or "")`` is wrong here:
    ``float('nan')`` is truthy, so a missing pandas cell survives that idiom
    and renders as the literal string "nan" - a missing value disguised as
    data, which is exactly what this report must not do.
    """
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "&mdash;"
    text = str(value)
    if not text or text == "nan":
        return "&mdash;"
    return esc(text[:n] + ("\u2026" if len(text) > n else ""))


def slug(value: Any) -> str:
    """Filesystem/DOM-safe fragment for element ids."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "na"
    return "".join(ch if ch.isalnum() else "-" for ch in str(value)) or "na"


def parse_spans(value: Any) -> list[tuple[int, int]]:
    """Parse extracellular spans from a list or its CSV string form."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    if isinstance(value, str):
        text = value.strip()
        if not text or text in ("[]", "nan"):
            return []
        try:
            value = ast.literal_eval(text)
        except (ValueError, SyntaxError):
            return []
    spans = []
    for item in value or []:
        try:
            lo, hi = int(item[0]), int(item[1])
        except (TypeError, ValueError, IndexError):
            continue
        spans.append((lo, hi))
    return spans


def split_ids(value: Any) -> list[str]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    return [p for p in str(value).split(";") if p and p != "nan"]


# --------------------------------------------------------------------------
# SVG components
# --------------------------------------------------------------------------

def svg_funnel(steps: list[tuple[str, int, str]], width: int = 860) -> str:
    """Horizontal funnel: (label, value, tooltip) per stage."""
    if not steps:
        return ""
    n = len(steps)
    gap, pad = 10, 16
    bar_w = (width - 2 * pad - gap * (n - 1)) / n
    vmax = max(v for _, v, _ in steps) or 1
    h, top, bot = 168, 34, 44
    body = h - top - bot
    parts = [f'<svg viewBox="0 0 {width} {h}" class="chart" role="img" '
             f'aria-label="stage funnel">']
    for i, (label, value, tip) in enumerate(steps):
        x = pad + i * (bar_w + gap)
        bh = max(4.0, body * (value / vmax))
        y = top + (body - bh)
        colour = PALETTE[i % len(PALETTE)]
        parts.append(
            f'<g class="fnl" data-tip="{esc(tip)}">'
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{bh:.1f}" '
            f'rx="4" fill="{colour}" fill-opacity="0.85"/>'
            f'<text x="{x + bar_w / 2:.1f}" y="{y - 8:.1f}" class="svg-val">{value:,}</text>'
            f'<text x="{x + bar_w / 2:.1f}" y="{h - 24:.1f}" class="svg-lbl">{esc(label)}</text>'
            f"</g>"
        )
        if i < n - 1:
            ax = x + bar_w + gap / 2
            parts.append(
                f'<path d="M{ax - 3:.1f},{top + body / 2:.1f} l6,0 m-3,-3 l3,3 l-3,3" '
                f'stroke="#8c959f" stroke-width="1.4" fill="none"/>'
            )
    parts.append("</svg>")
    return "".join(parts)


def svg_score_bars(components: list[tuple[str, float, float]], width: int = 520) -> str:
    """Per-component raw value (light) and weighted contribution (solid)."""
    if not components:
        return ""
    row_h, pad_l, pad_t = 22, 150, 6
    h = pad_t + row_h * len(components) + 8
    inner = width - pad_l - 70
    parts = [f'<svg viewBox="0 0 {width} {h}" class="chart">']
    for i, (name, raw, wtd) in enumerate(components):
        y = pad_t + i * row_h
        colour = PALETTE[i % len(PALETTE)]
        raw_w = max(0.0, min(1.0, raw)) * inner
        wtd_w = max(0.0, min(1.0, wtd / 0.25)) * inner if wtd else 0.0
        parts.append(
            f'<text x="{pad_l - 8}" y="{y + 12}" class="svg-lbl" text-anchor="end">'
            f"{esc(name)}</text>"
            f'<rect x="{pad_l}" y="{y + 4}" width="{raw_w:.1f}" height="6" rx="3" '
            f'fill="{colour}" fill-opacity="0.28"/>'
            f'<rect x="{pad_l}" y="{y + 11}" width="{wtd_w:.1f}" height="6" rx="3" '
            f'fill="{colour}"/>'
            f'<text x="{pad_l + inner + 6}" y="{y + 14}" class="svg-num">'
            f"{raw:.2f}</text>"
        )
    parts.append("</svg>")
    return "".join(parts)


def svg_sequence_track(
    spans: Sequence[tuple[int, int]],
    interface: Sequence[int],
    hotspots: Sequence[int],
    width: int = 820,
) -> str:
    """Positions along the chain: ectodomain spans, interface, hotspot patch."""
    marks = [p for p in list(interface) + list(hotspots) if p]
    if not spans and not marks:
        return ""
    lo = 1
    hi = max([s[1] for s in spans] + marks + [2])
    pad, h = 14, 86
    inner = width - 2 * pad

    def sx(pos: int) -> float:
        return pad + inner * (pos - lo) / max(1, hi - lo)

    parts = [f'<svg viewBox="0 0 {width} {h}" class="chart" role="img" '
             f'aria-label="interface positions along the chain">']
    parts.append(f'<rect x="{pad}" y="30" width="{inner}" height="12" rx="6" '
                 f'fill="#eaeef2"/>')
    for s_lo, s_hi in spans:
        parts.append(
            f'<rect x="{sx(s_lo):.1f}" y="30" width="{max(1.0, sx(s_hi) - sx(s_lo)):.1f}" '
            f'height="12" rx="6" fill="#1a7f37" fill-opacity="0.28" '
            f'data-tip="extracellular span {s_lo}-{s_hi} (UniProt numbering)"/>'
        )
    hot = set(hotspots)
    for pos in sorted(set(interface) - hot):
        parts.append(f'<line x1="{sx(pos):.1f}" y1="26" x2="{sx(pos):.1f}" y2="46" '
                     f'stroke="#0969da" stroke-width="1.5" stroke-opacity="0.55" '
                     f'data-tip="interface residue {pos}"/>')
    for pos in sorted(hot):
        parts.append(f'<line x1="{sx(pos):.1f}" y1="20" x2="{sx(pos):.1f}" y2="52" '
                     f'stroke="#cf222e" stroke-width="2.4" '
                     f'data-tip="hotspot {pos}"/>')
    for pos in (lo, hi):
        parts.append(f'<text x="{sx(pos):.1f}" y="70" class="svg-lbl" '
                     f'text-anchor="middle">{pos}</text>')
    parts.append(
        '<g class="legend">'
        '<rect x="14" y="4" width="10" height="8" rx="2" fill="#1a7f37" fill-opacity="0.28"/>'
        '<text x="29" y="12" class="svg-lbl">ectodomain</text>'
        '<line x1="112" y1="4" x2="112" y2="12" stroke="#0969da" stroke-width="1.5"/>'
        '<text x="120" y="12" class="svg-lbl">interface</text>'
        '<line x1="186" y1="4" x2="186" y2="12" stroke="#cf222e" stroke-width="2.4"/>'
        '<text x="194" y="12" class="svg-lbl">hotspot patch</text>'
        "</g>"
    )
    parts.append("</svg>")
    return "".join(parts)


def svg_structure(projection: dict, width: int = 820, height: int = 380) -> str:
    """Draw a persisted 2-D CA-trace projection of a target-binder complex.

    This is a backbone schematic, not a molecular surface, and is labelled as
    such: the projection was computed upstream (first two principal components
    of the complex's CA coordinates) and stored, so the report draws from a
    record rather than re-reading coordinates.
    """
    chains = projection.get("chains") or []
    if not chains:
        return ""
    x0, y0, x1, y1 = projection.get("bounds") or [0, 0, 1, 1]
    span_x, span_y = max(1e-6, x1 - x0), max(1e-6, y1 - y0)
    pad = 26
    scale = min((width - 2 * pad) / span_x, (height - 2 * pad) / span_y)
    off_x = pad + ((width - 2 * pad) - span_x * scale) / 2
    off_y = pad + ((height - 2 * pad) - span_y * scale) / 2

    def px(p):
        return (off_x + (p[0] - x0) * scale, off_y + (y1 - p[1]) * scale)

    parts = [f'<svg viewBox="0 0 {width} {height}" class="chart struct" role="img" '
             f'aria-label="backbone projection of the target-binder complex">']
    # Binder behind, target in front, so highlighted epitope dots stay visible.
    for chain in sorted(chains, key=lambda c: c["role"] != "binder"):
        pts = [px(p) for p in chain["points"]]
        if len(pts) < 2:
            continue
        path = "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts)
        if chain["role"] == "binder":
            parts.append(
                f'<path d="{path}" fill="none" stroke="#8250df" stroke-width="2.6" '
                f'stroke-opacity="0.85" stroke-linejoin="round" '
                f'data-tip="binder chain {esc(chain["chain"])} '
                f'({len(chain["points"])} residues modelled)"/>'
            )
        else:
            parts.append(
                f'<path d="{path}" fill="none" stroke="#8c959f" stroke-width="1.8" '
                f'stroke-opacity="0.75" stroke-linejoin="round" '
                f'data-tip="target chain {esc(chain["chain"])} '
                f'({len(chain["points"])} residues modelled)"/>'
            )
            flags = chain.get("flags") or []
            resids = chain.get("residues") or []
            for i, (x, y) in enumerate(pts):
                flag = flags[i] if i < len(flags) else ""
                if not flag:
                    continue
                rid = resids[i] if i < len(resids) else "?"
                if flag == "hotspot":
                    parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.2" '
                                 f'fill="#cf222e" data-tip="hotspot {esc(rid)}"/>')
                else:
                    parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.4" '
                                 f'fill="#0969da" fill-opacity="0.75" '
                                 f'data-tip="epitope residue {esc(rid)}"/>')
    parts.append(
        '<g class="legend">'
        '<line x1="14" y1="12" x2="34" y2="12" stroke="#8250df" stroke-width="2.6"/>'
        '<text x="39" y="16" class="svg-lbl">binder backbone</text>'
        '<line x1="146" y1="12" x2="166" y2="12" stroke="#8c959f" stroke-width="1.8"/>'
        '<text x="171" y="16" class="svg-lbl">target backbone</text>'
        '<circle cx="288" cy="12" r="2.4" fill="#0969da" fill-opacity="0.75"/>'
        '<text x="296" y="16" class="svg-lbl">epitope</text>'
        '<circle cx="356" cy="12" r="4.2" fill="#cf222e"/>'
        '<text x="366" y="16" class="svg-lbl">hotspot</text>'
        "</g>"
    )
    parts.append("</svg>")
    return "".join(parts)


def svg_residue_bars(rows: list[dict], width: int = 820, top: int = 16) -> str:
    """Buried area per interface residue, hotspots highlighted."""
    rows = [r for r in rows if r.get("delta_sasa") is not None][:top]
    if not rows:
        return ""
    bar_w, gap, h = 34, 8, 150
    base, plot = 112, 92
    vmax = max(float(r["delta_sasa"]) for r in rows) or 1.0
    w = max(width, len(rows) * (bar_w + gap) + 30)
    parts = [f'<svg viewBox="0 0 {w} {h}" class="chart">']
    for i, r in enumerate(rows):
        v = float(r["delta_sasa"])
        bh = max(2.0, plot * v / vmax)
        x = 18 + i * (bar_w + gap)
        is_hot = bool(r.get("is_hotspot"))
        colour = "#cf222e" if is_hot else "#0969da"
        dist = r.get("min_distance")
        tip = (f"{r.get('residue','?')}{r.get('auth_seq_id','')} "
               f"(UniProt {r.get('uniprot_pos') or '—'}) "
               f"buried {v:.1f} A^2"
               + (f", closest contact {float(dist):.2f} A" if dist and not (
                   isinstance(dist, float) and math.isnan(dist)) else ", no contact within cutoff")
               + (" — hotspot" if is_hot else ""))
        parts.append(
            f'<g data-tip="{esc(tip)}">'
            f'<rect x="{x}" y="{base - bh:.1f}" width="{bar_w}" height="{bh:.1f}" rx="3" '
            f'fill="{colour}" fill-opacity="{0.95 if is_hot else 0.6}"/>'
            f'<text x="{x + bar_w / 2}" y="{base + 13}" class="svg-lbl" '
            f'text-anchor="middle">{esc(str(r.get("residue", ""))[:3])}</text>'
            f'<text x="{x + bar_w / 2}" y="{base + 26}" class="svg-num" '
            f'text-anchor="middle">{esc(r.get("auth_seq_id"))}</text>'
            f"</g>"
        )
    parts.append(f'<text x="18" y="16" class="svg-lbl">buried area per residue '
                 f'(A^2), max {vmax:.0f}</text>')
    parts.append("</svg>")
    return "".join(parts)


# --------------------------------------------------------------------------
# HTML primitives
# --------------------------------------------------------------------------

def card(title: str, body: str, *, note: str | None = None, sub: str | None = None) -> str:
    head = f'<div class="card-head"><h3>{esc(title)}</h3>'
    if sub:
        head += f'<span class="sub">{esc(sub)}</span>'
    head += "</div>"
    note_html = f'<p class="note">{note}</p>' if note else ""
    return f'<section class="card">{head}{note_html}{body}</section>'


def absent_card(title: str, what: str) -> str:
    return card(
        title,
        f'<p class="absent"><strong>Not present in this run.</strong> {esc(what)} '
        f"No value is shown for it, and nothing here stands in for it.</p>",
    )


def kv(pairs: Iterable[tuple[str, str]]) -> str:
    items = "".join(
        f'<div class="kv-row"><dt>{esc(k)}</dt><dd>{v}</dd></div>' for k, v in pairs
    )
    return f'<dl class="kv">{items}</dl>'


def stat_row(stats: list[tuple[str, str, str | None]]) -> str:
    cells = "".join(
        f'<div class="stat"><span class="stat-v">{v}</span>'
        f'<span class="stat-k">{esc(k)}</span>'
        + (f'<span class="stat-n">{esc(n)}</span>' if n else "")
        + "</div>"
        for k, v, n in stats
    )
    return f'<div class="stats">{cells}</div>'


def table(
    tid: str,
    headers: Sequence[tuple[str, str]],
    rows: Sequence[Sequence[str]],
    *,
    details: Sequence[str] | None = None,
    filter_label: str = "Filter",
    initial_sort: int | None = None,
) -> str:
    """Sortable, filterable table; `details` adds one expandable row each."""
    if not rows:
        return '<p class="absent">No rows.</p>'
    ths = "".join(
        f'<th data-type="{t}" onclick="sortTable(\'{tid}\',{i})">'
        f'{esc(h)}<span class="arrow"></span></th>'
        for i, (h, t) in enumerate(headers)
    )
    body = []
    for r, (cells) in enumerate(rows):
        detail = details[r] if details and r < len(details) else None
        rid = f"{tid}-r{r}"
        cls = ' class="has-detail" onclick="toggleDetail(\'%s\')"' % rid if detail else ""
        tds = "".join(f"<td>{c}</td>" for c in cells)
        first = f'<tr data-row{cls}>{tds}</tr>'
        body.append(first)
        if detail:
            body.append(
                f'<tr class="detail" id="{rid}"><td colspan="{len(headers)}">'
                f'<div class="detail-in">{detail}</div></td></tr>'
            )
    return (
        f'<div class="tbl-tools">'
        f'<input class="filter" type="search" placeholder="{esc(filter_label)}" '
        f'oninput="filterTable(\'{tid}\',this.value)" aria-label="{esc(filter_label)}">'
        f'<span class="count" id="{tid}-count"></span></div>'
        f'<div class="tbl-wrap"><table id="{tid}" class="grid" '
        f'data-sort="{initial_sort if initial_sort is not None else ""}">'
        f"<thead><tr>{ths}</tr></thead><tbody>{''.join(body)}</tbody></table></div>"
    )


def copy_btn(text: str, label: str = "copy") -> str:
    return (f'<button class="copy" onclick="copyText(this,\'{esc(text)}\')" '
            f'type="button">{esc(label)}</button>')


def pill(text: str, kind: str = "") -> str:
    return f'<span class="pill {kind}">{esc(text)}</span>'


# --------------------------------------------------------------------------
# sections
# --------------------------------------------------------------------------

def _designed_stat(run: dict) -> tuple[str, str, str]:
    """The 'binders designed' headline, read from the run rather than assumed.

    This is the only place the report states design state, so it must not be
    hardcoded: a run that later gains a design record has to be reported
    honestly, and a design job that ran without succeeding yields no candidates.
    """
    dr = run.get("design_run") or {}
    cands = run.get("candidates")
    n_cand: int | None = None
    if isinstance(cands, pd.DataFrame):
        n_cand = len(cands)
    elif isinstance(cands, list):
        n_cand = len(cands)

    state = dr.get("status") or dr.get("state")
    if not dr and n_cand is None:
        return ("binders designed", "0",
                "specifications only; no design run in this directory")
    if state and state != "succeeded":
        return ("binders designed", "0",
                f"design run {state} - a run that has not succeeded yields no "
                f"candidates, and none are inferred from a partial run")
    if n_cand is None:
        return ("binders designed", "&mdash;",
                f"design run {state or 'present'} but no candidate record found, "
                f"so the count is unknown rather than zero")
    return ("binder candidates", f"{n_cand:,}",
            "a succeeded design job is not the same as a candidate passing "
            "evaluation")


def sec_overview(run: dict) -> str:
    d0 = run.get("disease") or {}
    ev = run.get("evidence") or {}
    cx = run.get("complexes") or {}
    epi = run.get("epitopes")
    sm = run.get("specs_manifest") or {}
    disease = (d0.get("disease") or {})

    n_epi_ok = n_epi = n_epi_targets = None
    if isinstance(epi, pd.DataFrame) and not epi.empty:
        n_epi = len(epi)
        ok = epi[epi["status"] == "accepted"]
        n_epi_ok, n_epi_targets = len(ok), ok["symbol"].nunique()

    steps = []
    if ev.get("pool_size"):
        steps.append(("candidates", int(ev["pool_size"]),
                      f"{ev['pool_size']} of {ev.get('association_total', '?')} "
                      f"disease-associated targets pooled"))
    if ev.get("n_surface_accessible") is not None:
        steps.append(("reachable", int(ev["n_surface_accessible"]),
                      "annotated extracellular domain or secreted, by UniProt topology"))
    if cx.get("n_targets_admitted") is not None:
        steps.append(("has complex", int(cx["n_targets_admitted"]),
                      "experimental complex with a native interaction partner"))
    if n_epi_ok is not None:
        steps.append(("epitopes", int(n_epi_ok),
                      f"{n_epi_ok} accepted of {n_epi} computed, "
                      f"across {n_epi_targets} targets"))
    if sm.get("n_specs") is not None:
        steps.append(("specs", int(sm["n_specs"]),
                      "BindCraft2 target/campaign pairs written - none submitted"))

    stats = [
        ("indication", esc(disease.get("name")), esc(disease.get("id"))),
        ("resolution rule", esc(d0.get("selection_rule")),
         "ambiguous" if d0.get("ambiguous") else "unambiguous"),
        _designed_stat(run),
    ]
    if isinstance(epi, pd.DataFrame) and n_epi_ok:
        ok = epi[epi["status"] == "accepted"]
        stats.insert(2, ("median epitope",
                         f'{ok["buried_area_total"].median():,.0f} &#8491;&sup2;',
                         "buried area"))

    body = stat_row(stats) + svg_funnel(steps)
    body += ('<p class="note">Each bar is a gate, not a score cut. Targets removed at '
             'any stage stay in the tables below with the reason they were removed.</p>')
    return card("Run overview", body,
                sub=disease.get("name") or "indication2binder run")


def sec_provenance(run: dict) -> str:
    rows = []
    paths = run.get("_paths", {})
    manifest = run.get("run_manifest") or {}
    origins = manifest.get("origins") or {}
    default_origin = manifest.get("origin")
    for layer in LAYERS:
        if layer not in run:
            rows.append([esc(layer), pill("absent", "bad"), "&mdash;"])
            continue
        origin = origins.get(layer, default_origin)
        kind = {"real": "ok", "cached_real": "warn", "fixture": "bad"}.get(origin or "", "warn")
        rows.append([
            esc(layer),
            pill(origin or "unlabelled", kind),
            f'<code>{esc(paths.get(layer, ""))}</code>',
        ])
    unlabelled = [
        layer for layer in LAYERS
        if layer in run and not origins.get(layer, default_origin)
    ]
    note = None
    if unlabelled:
        note = ('<span class="warn-text">Provenance unlabelled for '
                f'{len(unlabelled)} layer(s).</span> Write a top-level '
                '<code>run_manifest.json</code> with <code>{"origin": "real"}</code> '
                'or per-layer <code>origins</code> to label them. Until then these '
                'are shown as unlabelled rather than assumed to be real results.')
    return card(
        "Data provenance",
        table("tbl-prov", [("layer", "text"), ("origin", "text"), ("file", "text")], rows,
              filter_label="Filter layers"),
        note=note,
    )


def sec_disease(run: dict) -> str:
    d0 = run.get("disease")
    if not d0:
        return absent_card("Indication resolution",
                           "The stage 0 disease record was not found.")
    disease = d0.get("disease") or {}
    desc = d0.get("descendants") or {}
    pairs = [
        ("resolved id", f'<code>{esc(disease.get("id"))}</code>'),
        ("name", esc(disease.get("name"))),
        ("query", f'<code>{esc(d0.get("query"))}</code>'),
        ("rule", esc(d0.get("selection_rule"))),
        ("ambiguous", pill("yes", "warn") if d0.get("ambiguous") else pill("no", "ok")),
        ("association scope", esc(d0.get("association_scope"))),
        ("children / descendants",
         f'{esc(len(d0.get("children") or []))} / {esc(desc.get("count"))}'),
    ]
    cands = d0.get("candidates_considered") or []
    rows = [[f'<code>{esc(c.get("id"))}</code>', esc(c.get("name"))] for c in cands]
    body = kv(pairs)
    if rows:
        body += ('<h4>Candidates considered</h4>'
                 '<p class="note">Recorded so the rejected alternatives are visible; '
                 'every downstream number is conditioned on the chosen id.</p>'
                 + table("tbl-cand", [("id", "text"), ("name", "text")], rows,
                         filter_label="Filter candidates"))
    return card("Indication resolution", body, sub=disease.get("id"))


def _evidence_index(run: dict) -> dict[str, dict]:
    ev = run.get("evidence") or {}
    return {t["symbol"]: t for t in (ev.get("targets") or []) if t.get("symbol")}


def _target_detail(row: pd.Series, ev: dict | None, weights: dict) -> str:
    comps = []
    for name in weights:
        raw = row.get(f"raw_{name}")
        wtd = row.get(f"wtd_{name}")
        if raw is None or (isinstance(raw, float) and math.isnan(raw)):
            continue
        comps.append((name.replace("_", " "), float(raw), float(wtd or 0.0)))
    body = ""
    if comps:
        body += ('<h5>Score components</h5>'
                 '<p class="note">Light bar is the raw component value (0-1); solid bar '
                 'is its weighted contribution. These reconstruct the composite exactly.'
                 '</p>' + svg_score_bars(comps))
    pen = row.get("penalty_reasons")
    pairs = [
        ("composite", num(row.get("composite_score"), 4)),
        ("penalty multiplier", num(row.get("penalty_multiplier"), 3)),
        ("penalty reasons", esc(pen) if isinstance(pen, str) and pen else "none"),
        ("accessibility", esc(row.get("accessibility_mode")) if row.get("passes_accessibility")
         else f'<span class="warn-text">excluded</span> &mdash; {esc(row.get("accessibility_reason"))}'),
        ("ectodomain residues", num(row.get("ectodomain_residues"), 0)),
        ("target class", esc(row.get("target_class"))),
        ("known drugs", num(row.get("n_known_drugs"), 0)),
        ("biologic precedent",
         pill("yes", "ok") if row.get("has_biologic_precedent") else pill("no")),
        ("trials (stopped)",
         f'{num(row.get("n_trials"), 0)} ({num(row.get("n_stopped_trials"), 0)})'),
        ("safety liabilities", num(row.get("n_safety_liabilities"), 0)),
        ("interaction partners", esc(row.get("top_partners"))),
    ]
    if ev:
        ds = ev.get("datatype_scores") or {}
        if ds:
            pairs.append(("datatype scores", ", ".join(
                f"{k}&nbsp;{v:.2f}" for k, v in sorted(ds.items(), key=lambda x: -x[1])[:6])))
        tract = (ev.get("tractability") or {}).get("AB") or []
        if tract:
            pairs.append(("antibody tractability", ", ".join(esc(t) for t in tract)))
        hint = ev.get("location_keyword_hint") or {}
        if hint and hint.get("agrees_with_topology") is False:
            pairs.append((
                "keyword-hint disagreement",
                '<span class="warn-text">location keywords and topology disagree</span> '
                f'(keyword said {"reachable" if hint.get("hint") else "not reachable"}; '
                f'topology decided)'))
        trials = ev.get("trials") or {}
        stopped = trials.get("stopped_trials") or []
        if stopped:
            pairs.append(("stopped trials", "; ".join(
                f'{esc(s.get("nct_id"))} ({esc(s.get("status"))})' for s in stopped[:3])))
    body += kv(pairs)
    return body


def sec_targets(run: dict) -> str:
    ranked = run.get("ranked")
    if not isinstance(ranked, pd.DataFrame) or ranked.empty:
        return absent_card("Target ranking", "The stage 2 ranked-target table was not found.")
    weights = ((run.get("scoring_meta") or {}).get("weights")) or {}
    evidx = _evidence_index(run)

    headers = [("#", "num"), ("target", "text"), ("score", "num"), ("gate", "text"),
               ("ectodomain", "num"), ("drugs", "num"), ("trials", "num"),
               ("safety", "num"), ("class", "text")]
    rows, details = [], []
    for _, r in ranked.iterrows():
        passes = bool(r.get("passes_accessibility"))
        gate = (pill("reachable", "ok") if passes
                else pill("excluded", "bad"))
        rows.append([
            num(r.get("rank"), 0),
            f'<strong>{esc(r.get("symbol"))}</strong> '
            f'<code class="dim">{esc(r.get("uniprot"))}</code>',
            f'<strong>{num(r.get("composite_score"), 3)}</strong>' if passes
            else num(r.get("composite_score"), 3),
            gate,
            num(r.get("ectodomain_residues"), 0),
            num(r.get("n_known_drugs"), 0),
            num(r.get("n_trials"), 0),
            num(r.get("n_safety_liabilities"), 0),
            trunc(r.get("target_class"), 26),
        ])
        details.append(_target_detail(r, evidx.get(r.get("symbol")), weights))

    n_pass = int(ranked["passes_accessibility"].sum())
    note = (f'All {len(ranked)} pooled targets are listed. {n_pass} pass the '
            f'accessibility gate; the other {len(ranked) - n_pass} are retained with '
            f'their exclusion reason rather than dropped. Click any row for its score '
            f'breakdown and evidence.')
    return card("Target ranking", table(
        "tbl-targets", headers, rows, details=details,
        filter_label="Filter by symbol, class, gate..."), note=note,
        sub=f"{len(ranked)} targets")


def sec_rubric(run: dict) -> str:
    sm = run.get("scoring_meta")
    if not sm:
        return absent_card("Scoring rubric", "The stage 2 scoring metadata was not found.")
    weights = sm.get("weights") or {}
    rows = [[esc(k.replace("_", " ")), num(v, 3)] for k, v in weights.items()]
    body = ""
    if rows:
        body += ("<h4>Weights</h4>"
                 + table("tbl-w", [("component", "text"), ("weight", "num")], rows,
                         filter_label="Filter components"))
    pen = sm.get("penalties") or {}
    if pen:
        body += "<h4>Penalties</h4>" + kv([(k.replace("_", " "), num(v, 3))
                                             for k, v in pen.items()])
    sat = sm.get("saturation") or {}
    if sat:
        body += "<h4>Count saturation</h4>" + kv(
            [(k.replace("_", " "), num(v, 0)) for k, v in sat.items()])
    caveats = sm.get("caveats") or []
    if caveats:
        body += ("<h4>Caveats recorded with the ranking</h4><ul class=\"bullets\">"
                 + "".join(f"<li>{esc(c)}</li>" for c in caveats) + "</ul>")
    return card("Scoring rubric", body)


def sec_complexes(run: dict) -> str:
    cx = run.get("complexes")
    if not cx:
        return absent_card("Experimental-complex gate",
                           "The stage 3 complex record was not found.")
    gate = cx.get("gate") or {}
    headers = [("target", "text"), ("verdict", "text"), ("searched", "num"),
               ("admitted", "num"), ("best complex", "text"), ("partner", "text")]
    rows, details = [], []
    for t in cx.get("targets") or []:
        admitted = t.get("admitted") or []
        best = admitted[0] if admitted else None
        rows.append([
            f'<strong>{esc(t.get("symbol"))}</strong>',
            pill("admitted", "ok") if admitted else pill("dropped", "bad"),
            num(t.get("n_searched"), 0),
            num(t.get("n_admitted"), 0),
            (f'<code>{esc(best["entry_id"])}</code> '
             f'<span class="dim">{num(best.get("resolution"), 2)}&#8491;</span>')
            if best else "&mdash;",
            esc(", ".join(best.get("partner_symbols") or [])) if best else "&mdash;",
        ])

        d = ""
        if not admitted:
            counts = t.get("rejection_reason_counts") or {}
            d += (f'<p class="absent"><strong>Dropped:</strong> '
                  f'{esc(t.get("drop_reason") or "no qualifying complex")}. '
                  + (f'Rejection reasons: {esc(json.dumps(counts))}.' if counts
                     else 'No candidate entries were returned by the structure search.')
                  + " This target is not carried forward; nothing is substituted for it.</p>")
        else:
            crows = [[
                f'<code>{esc(a.get("entry_id"))}</code>',
                num(a.get("resolution"), 2),
                esc(a.get("method")),
                esc(", ".join(a.get("partner_symbols") or [])),
                num((a.get("target_entity") or {}).get("length"), 0),
                esc(", ".join(a.get("other_roles") or []) or "—"),
                trunc(a.get("title"), 70),
            ] for a in admitted]
            d += table(
                f'tbl-cx-{slug(t.get("symbol"))}',
                [("PDB", "text"), ("res &#8491;", "num"), ("method", "text"),
                 ("partner", "text"), ("target aa", "num"),
                 ("other chains", "text"), ("title", "text")],
                crows, filter_label="Filter complexes")
        rej = t.get("rejected_examples") or []
        if rej:
            d += ('<h5>Rejected examples</h5><ul class="bullets">' + "".join(
                f'<li><code>{esc(x.get("entry_id"))}</code> &mdash; {esc(x.get("reason"))}</li>'
                for x in rej) + "</ul>")
        details.append(d)

    note = (f'Gate: method in {esc(", ".join(gate.get("accepted_methods") or []))}, '
            f'resolution &le; {num(gate.get("max_resolution"), 1)}&#8491;, target chain '
            f'&ge; {num(gate.get("min_target_chain_residues"), 0)} residues, and at least '
            f'one chain that is a known interaction partner. Antibody and Fab chains are '
            f'classified separately: they show an epitope is targetable by a biologic, but '
            f'blocking a defined interaction needs the native partner.')
    if cx.get("note"):
        note += f' {esc(cx["note"])}'
    return card("Experimental-complex gate",
                table("tbl-cx", headers, rows, details=details,
                      filter_label="Filter by target or PDB id") ,
                note=note,
                sub=f'{cx.get("n_targets_admitted")} of {cx.get("n_targets_in")} admitted')


def _residue_index(run: dict) -> dict[tuple, list[dict]]:
    res = run.get("residues")
    if not isinstance(res, pd.DataFrame) or res.empty:
        return {}
    has_partner = "partner_chain" in res.columns
    out: dict[tuple, list[dict]] = {}
    for _, r in res.iterrows():
        key = (r.get("symbol"), r.get("entry_id"), r.get("target_chain"),
               r.get("partner_chain") if has_partner else None)
        out.setdefault(key, []).append(r.to_dict())
    for rows in out.values():
        rows.sort(key=lambda x: -(x.get("delta_sasa") or 0))
    return out


def sec_epitopes(run: dict) -> str:
    epi = run.get("epitopes")
    if not isinstance(epi, pd.DataFrame) or epi.empty:
        return absent_card("Epitopes", "The stage 4 epitope table was not found.")
    meta = run.get("interface_meta") or {}
    params = meta.get("interface_params") or {}
    ridx = _residue_index(run)
    spans_by_symbol = {
        t.get("symbol"): parse_spans(t.get("extracellular_spans"))
        for t in ((run.get("complexes") or {}).get("targets") or [])
    }

    headers = [("#", "num"), ("target", "text"), ("partner", "text"), ("PDB", "text"),
               ("res &#8491;", "num"), ("buried &#8491;&sup2;", "num"),
               ("residues", "num"), ("in ecto", "num"), ("score", "num"),
               ("status", "text")]
    rows, details = [], []
    ordered = epi.sort_values(
        ["status", "combined_score"], ascending=[True, False]
    ) if "combined_score" in epi.columns else epi
    for _, r in ordered.iterrows():
        accepted = r.get("status") == "accepted"
        rows.append([
            num(r.get("epitope_rank"), 0) if accepted else "&mdash;",
            f'<strong>{esc(r.get("symbol"))}</strong>',
            esc(r.get("partner_symbol")),
            f'<code>{esc(r.get("entry_id"))}</code>'
            f'<span class="dim">/{esc(r.get("target_chain"))}</span>',
            num(r.get("resolution"), 2),
            num(r.get("buried_area_total"), 0),
            num(r.get("n_interface_residues"), 0),
            pct(r.get("fraction_in_ectodomain")),
            num(r.get("combined_score"), 3),
            pill("accepted", "ok") if accepted else pill("rejected", "bad"),
        ])

        d = ""
        if not accepted:
            d += (f'<p class="absent"><strong>Rejected.</strong> '
                  f'{esc(r.get("reason"))}</p>')
        iface_up = [int(x) for x in split_ids(r.get("interface_uniprot")) if x.isdigit()]
        hot_up = [int(x) for x in split_ids(r.get("hotspots_uniprot")) if x.isdigit()]
        spans = spans_by_symbol.get(r.get("symbol")) or []
        track = svg_sequence_track(spans, iface_up, hot_up)
        if track:
            d += ("<h5>Position along the chain</h5>"
                  '<p class="note">UniProt numbering. The full chain length is not '
                  'recorded in the run outputs, so the axis ends at the last annotated '
                  'or contacting residue.</p>' + track)

        key = (r.get("symbol"), r.get("entry_id"), r.get("target_chain"),
               r.get("partner_chain"))
        rres = ridx.get(key) or ridx.get(
            (r.get("symbol"), r.get("entry_id"), r.get("target_chain"), None)) or []
        if rres:
            hot_auth = set(split_ids(r.get("hotspots_auth")))
            for x in rres:
                x["is_hotspot"] = str(x.get("auth_seq_id")) in hot_auth
            d += "<h5>Buried area per residue</h5>" + svg_residue_bars(rres)
            trows = [[
                f'{esc(x.get("residue"))}<strong>{esc(x.get("auth_seq_id"))}</strong>',
                num(x.get("uniprot_pos"), 0),
                num(x.get("delta_sasa"), 1),
                num(x.get("min_distance"), 2),
                pill("contact", "ok") if x.get("is_contact") else pill("burial only"),
                pill("hotspot", "bad") if x.get("is_hotspot") else "",
                pill("ecto", "ok") if x.get("in_ectodomain") else pill("outside", "warn"),
            ] for x in rres]
            d += table(
                f'tbl-res-{slug(r.get("symbol"))}-{slug(r.get("entry_id"))}-'
                f'{slug(r.get("partner_chain"))}',
                [("residue", "text"), ("UniProt", "num"), ("buried &#8491;&sup2;", "num"),
                 ("min dist &#8491;", "num"), ("criterion", "text"),
                 ("hotspot", "text"), ("topology", "text")],
                trows, filter_label="Filter residues")
        elif accepted:
            d += ('<p class="absent">Per-residue detail is not present for this '
                  'epitope in the run outputs.</p>')

        hs_auth = ",".join(split_ids(r.get("hotspots_auth")))
        d += kv([
            ("hotspots (author numbering)",
             f'<code>{esc(hs_auth)}</code> {copy_btn(hs_auth)}' if hs_auth else "&mdash;"),
            ("hotspots (UniProt)",
             f'<code>{esc(",".join(str(h) for h in hot_up))}</code>' if hot_up else "&mdash;"),
            ("hotspot patch radius",
             f'{num(r.get("hotspot_patch_radius"), 1)}&#8491;'),
            ("designability",
             f'{num(r.get("designability"), 3)} '
             f'<span class="dim">(size {num(r.get("size_score"), 2)}, '
             f'resolution {num(r.get("resolution_score"), 2)}, '
             f'hotspots {num(r.get("hotspot_score"), 2)})</span>'),
            ("method", esc(r.get("method"))),
            ("structure", f'<code>{esc(r.get("structure_file"))}</code>'),
            ("entry title", esc(r.get("title"))),
        ])
        details.append(d)

    n_ok = int((epi["status"] == "accepted").sum())
    note = (f'Interface residue = heavy-atom contact within '
            f'{num(params.get("contact_cutoff"), 1)}&#8491; <em>or</em> buried SASA '
            f'&ge; {num(params.get("min_delta_sasa"), 0)}&#8491;&sup2; (Shrake-Rupley). '
            f'Interfaces are computed per partner chain, not pooled, because a target '
            f'touching several partners at once gives a composite no single mini-binder '
            f'reproduces. An epitope is accepted only if &ge; '
            f'{pct(params.get("min_fraction_in_ectodomain"))} of its mapped residues fall '
            f'inside an annotated extracellular span &mdash; author numbering is mapped to '
            f'UniProt numbering through the SIFTS alignment, never assumed equal. '
            f'{n_ok} of {len(epi)} accepted; rejected rows are kept with their reason.')
    return card("Epitopes and hotspots",
                table("tbl-epi", headers, rows, details=details,
                      filter_label="Filter by target, partner, PDB id..."),
                note=note, sub=f"{n_ok} accepted")


def sec_specs(run: dict) -> str:
    sm = run.get("specs_manifest")
    if not sm:
        return absent_card("BindCraft2 specifications",
                           "No stage 5 specification manifest was found.")
    specs_dir = run.get("_specs_dir")
    headers = [("name", "text"), ("target", "text"), ("partner", "text"), ("PDB", "text"),
               ("buried &#8491;&sup2;", "num"), ("chain aa", "num"),
               ("hotspots", "text"), ("status", "text")]
    rows, details = [], []
    for s in sm.get("specs") or []:
        written = s.get("status") == "written"
        hs = s.get("hotspots_auth") or ""
        rows.append([
            f'<code>{esc(s.get("name"))}</code>',
            esc(s.get("symbol")),
            esc(s.get("partner")),
            f'<code>{esc(s.get("pdb_entry"))}</code>',
            num(s.get("buried_area"), 0),
            num(s.get("n_structure_residues"), 0),
            f'<code>{esc(hs)}</code> {copy_btn(hs)}' if hs else "&mdash;",
            pill("written", "ok") if written else pill("skipped", "bad"),
        ])
        d = ""
        if not written:
            d = f'<p class="absent"><strong>Skipped.</strong> {esc(s.get("reason"))}</p>'
        else:
            d = kv([
                ("UniProt", f'<code>{esc(s.get("uniprot"))}</code>'),
                ("target chain", esc(s.get("target_chain"))),
                ("resolution", f'{num(s.get("resolution"), 2)}&#8491;'),
                ("hotspots (UniProt)", f'<code>{esc(s.get("hotspots_uniprot"))}</code>'),
                ("combined score", num(s.get("combined_score"), 3)),
            ])
            if specs_dir:
                for label, rel in (("target definition", f'settings/target/{s["name"]}.json'),
                                    ("campaign", f'campaigns/{s["name"]}.json')):
                    p = Path(specs_dir) / rel
                    if p.exists():
                        d += (f'<h5>{esc(label)} <code class="dim">{esc(rel)}</code></h5>'
                              f'<pre class="code">{esc(p.read_text())}</pre>')
        details.append(d)

    note = (f'{esc(sm.get("numbering", ""))}. Binder lengths '
            f'{esc(sm.get("binder_lengths"))}, '
            f'{esc(sm.get("number_of_final_designs"))} final designs per campaign'
            + (f', max {esc(sm.get("max_trajectories"))} trajectories'
               if sm.get("max_trajectories") else "") + ". "
            "The emitted structure contains the target chain only: the partner is "
            "stripped so a design cannot score well by reproducing the partner's own "
            "contacts.")
    return card("BindCraft2 specifications",
                table("tbl-specs", headers, rows, details=details,
                      filter_label="Filter specifications"),
                note=note, sub=f'{sm.get("n_specs")} written')


def sec_binders(run: dict) -> str:
    """Existing binders of these targets - none produced by this pipeline."""
    rec = run.get("binders")
    if not rec:
        return absent_card(
            "Existing binders",
            "No stage 6 binder record was found, so no existing binder of these "
            "targets is shown. Run `ind2b stage6` to collect them.")

    lead = (
        '<p class="verdict warn"><strong>None of these binders came from this '
        'pipeline.</strong> Every molecule below already exists &mdash; a designed '
        'binder deposited in the PDB, or a clinical antibody. They are shown '
        'because an existing binder against the same surface is evidence that '
        'the surface is <em>bindable</em>; that is not evidence that a new '
        'binder against it will have useful affinity or specificity.</p>'
    )

    # Mini-binder-scale designs are the ones comparable to a BindCraft2 output.
    minis = [
        (t, b) for t in rec["targets"] for b in t["binders"]
        if b.get("is_minibinder_scale") and (b.get("vs_selected_epitope") or {}).get("comparable")
    ]
    covered = [
        (t, b) for t, b in minis
        if (b["vs_selected_epitope"].get("fraction_of_hotspots_covered") or 0) >= 0.99
    ]
    if minis:
        lead += (
            f'<p class="note">Of the {len(minis)} published de novo binder(s) at '
            f'mini-binder scale (&le; {esc(rec.get("minibinder_max_residues"))} residues) '
            f'whose footprint could be mapped, <strong>{len(covered)}</strong> engage every '
            f'hotspot this pipeline selected for that target. Larger designed scaffolds '
            f'(DARPins and designed scFvs) in the table below mostly do not &mdash; they '
            f'were raised against other surfaces. Read a zero there carefully: where '
            f'the deposited construct does not even contain the selected epitope, '
            f'the zero is silence about that surface rather than evidence against '
            f'it, and the row says which case applies.</p>'
            '<p class="note">Calibrate that agreement honestly: these are the '
            'canonical ligand-binding interfaces of well-studied receptors, which '
            'is where both a human designer and this pipeline would look first. '
            'Convergence on them is a sanity check that the epitope selection is '
            'not eccentric &mdash; it is not evidence that the pipeline finds '
            'non-obvious sites, and it does not extend to the targets here that '
            'have no published designed binder.</p>'
        )

    # --- example structures
    examples = [
        (t, b) for t in rec["targets"] for b in t["binders"] if b.get("projection")
    ]
    examples.sort(key=lambda tb: (not tb[1].get("is_minibinder_scale"),
                                   tb[1].get("resolution") or 99.0))
    gallery = ""
    for t, b in examples:
        cmp_ = b.get("vs_selected_epitope") or {}
        fp = b.get("footprint") or {}
        sel = b.get("selected_epitope") or {}
        cov = cmp_.get("fraction_of_hotspots_covered")
        headline = (
            f'{t["symbol"]} &mdash; {esc(b["binder_descriptions"][0])}'
            if b.get("binder_descriptions") else f'{t["symbol"]} binder'
        )
        badge = (pill(f'{b["binder_length"]} aa', "ok" if b.get("is_minibinder_scale") else "")
                 + " " + pill(b["kind"], "ok" if b["kind"] == "designed" else ""))
        note = None
        if cmp_.get("comparable"):
            if cov is not None and cov >= 0.99:
                note = (f'<span class="ok-text">This binder covers every hotspot '
                        f'selected for {esc(t["symbol"])}</span> '
                        f'({esc(",".join(str(h) for h in cmp_.get("hotspots_covered") or []))}), '
                        f'sharing {cmp_["n_shared"]} interface positions '
                        f'(Jaccard {num(cmp_["jaccard"], 2)}) with the epitope taken from '
                        f'<code>{esc(sel.get("entry_id"))}</code> vs '
                        f'{esc(sel.get("partner_symbol"))}.')
            elif cov is not None:
                note = (f'Covers {pct(cov)} of the selected hotspots '
                        f'({cmp_["n_shared"]} shared positions, Jaccard '
                        f'{num(cmp_["jaccard"], 2)}) &mdash; this binder engages a '
                        f'largely different surface from the one selected here.')
        body = (
            f'<div class="card-head"><h4 style="margin:0">{headline}</h4>'
            f'<span class="sub">{badge}</span></div>'
            + (f'<p class="note">{note}</p>' if note else "")
            + svg_structure(b["projection"])
            + f'<p class="caption">Backbone schematic of <code>{esc(b["entry_id"])}</code> '
              f'at {num(b.get("resolution"), 2)}&#8491;, projected onto the complex\'s two '
              f'principal axes. Lines trace C&alpha; positions only &mdash; this is not a '
              f'molecular surface, and apparent gaps are unmodelled residues. Highlighted '
              f'residues are the {esc(b["projection"].get("epitope_source", "epitope"))}. '
              f'Footprint buries {num(fp.get("buried_area_total"), 0)}&#8491;&sup2; over '
              f'{esc(fp.get("n_footprint_residues"))} residues.</p>'
        )
        gallery += f'<div class="example">{body}</div>'

    if gallery:
        gallery = ('<h4>Example complexes</h4>' + gallery)
    else:
        gallery = ('<p class="absent">No example complex could be projected, so no '
                   'structure is drawn here.</p>')

    # --- full table
    headers = [("target", "text"), ("kind", "text"), ("binder", "text"),
               ("PDB", "text"), ("res &#8491;", "num"), ("binder aa", "num"),
               ("buried &#8491;&sup2;", "num"), ("shared", "num"),
               ("hotspots covered", "num")]
    rows, details = [], []
    for t in rec["targets"]:
        for b in t["binders"]:
            cmp_ = b.get("vs_selected_epitope") or {}
            fp = b.get("footprint") or {}
            cov = cmp_.get("fraction_of_hotspots_covered")
            rows.append([
                f'<strong>{esc(t["symbol"])}</strong>',
                pill(b["kind"], "ok" if b["kind"] == "designed" else ""),
                trunc((b.get("binder_descriptions") or [None])[0], 34),
                f'<code>{esc(b["entry_id"])}</code>',
                num(b.get("resolution"), 2),
                num(b.get("binder_length"), 0),
                num(fp.get("buried_area_total"), 0),
                num(cmp_.get("n_shared"), 0) if cmp_.get("comparable") else "&mdash;",
                pct(cov) if cov is not None else "&mdash;",
            ])
            d = ""
            if not fp:
                d += ('<p class="absent">The target-side footprint could not be '
                      'computed for this complex, so no comparison with the selected '
                      'epitope is shown.</p>')
            elif not cmp_.get("comparable"):
                d += (f'<p class="absent">Not comparable: '
                      f'{esc(cmp_.get("reason", "footprint could not be mapped"))}.</p>')
            else:
                d += kv([
                    ("shared positions",
                     f'<code>{esc(",".join(str(x) for x in cmp_.get("shared_positions") or []))}</code>'),
                    ("hotspots covered",
                     f'<code>{esc(",".join(str(x) for x in cmp_.get("hotspots_covered") or []) or "none")}</code>'),
                    ("hotspots missed",
                     f'<code>{esc(",".join(str(x) for x in cmp_.get("hotspots_missed") or []) or "none")}</code>'),
                    ("Jaccard overlap", num(cmp_.get("jaccard"), 3)),
                    ("epitope covered", pct(cmp_.get("fraction_of_epitope_covered"))),
                    ("selected epitope in this construct",
                     (pill("yes", "ok") if cmp_.get("epitope_present_in_construct")
                      else pill("no", "warn"))
                     + f' <span class="dim">({esc(cmp_.get("n_epitope_positions_in_construct"))} '
                       f'of its positions modelled here)</span>'),
                    ("how to read the overlap", esc(cmp_.get("interpretation"))),
                ])
            if fp:
                top = fp.get("top_residues") or []
                if top:
                    d += ("<h5>Most-buried footprint residues</h5>" + table(
                        f'tbl-bfp-{slug(t["symbol"])}-{slug(b["entry_id"])}-{slug(b["kind"])}',
                        [("residue", "text"), ("UniProt", "num"),
                         ("buried &#8491;&sup2;", "num"), ("contact", "text")],
                        [[f'{esc(x.get("residue"))}<strong>{esc(x.get("auth_seq_id"))}</strong>',
                          num(x.get("uniprot_pos"), 0), num(x.get("delta_sasa"), 1),
                          pill("yes", "ok") if x.get("is_contact") else pill("no")]
                         for x in top],
                        filter_label="Filter residues"))
            d += kv([
                ("entry title", esc(b.get("title"))),
                ("method", esc(b.get("method"))),
                ("binder chains", esc(", ".join(b.get("binder_chains") or []))),
                ("chain used for footprint", esc(fp.get("binder_chain_used"))),
                ("construct UniProt range",
                 (f'{esc((fp.get("construct_uniprot_range") or [None, None])[0])}'
                  f'&ndash;{esc((fp.get("construct_uniprot_range") or [None, None])[1])}'
                  f' <span class="dim">({esc(fp.get("n_construct_residues_mapped"))} '
                  f'residues mapped)</span>')
                 if fp.get("construct_uniprot_range") else "&mdash;"),
                ("mini-binder scale",
                 pill("yes", "ok") if b.get("is_minibinder_scale") else pill("no")),
            ])
            details.append(d)

    clinical_rows = []
    for t in rec["targets"]:
        for d in t.get("clinical_binders") or []:
            clinical_rows.append([
                f'<strong>{esc(t["symbol"])}</strong>',
                esc(d.get("drug")),
                esc(d.get("drug_type")),
                pill(str(d.get("max_clinical_stage")),
                     "ok" if d.get("max_clinical_stage") == "APPROVAL" else "warn"),
            ])

    body = lead + gallery
    if rows:
        body += ("<h4>All analysed binder complexes</h4>"
                 '<p class="note">One complex per binder kind per target, best resolved '
                 'first. Click a row for the shared positions and the most-buried '
                 'footprint residues.</p>'
                 + table("tbl-binders", headers, rows, details=details,
                         filter_label="Filter by target, PDB id, kind..."))
    if clinical_rows:
        body += ("<h4>Clinical binders recorded for these targets</h4>"
                 '<p class="note">Biologic drugs from the run\'s own drug-precedent '
                 'layer. These are molecules in development or approved, not designs, '
                 'and their epitopes are not computed here unless a structure appears '
                 'in the table above.</p>'
                 + table("tbl-clin",
                         [("target", "text"), ("drug", "text"), ("modality", "text"),
                          ("stage", "text")],
                         clinical_rows, filter_label="Filter drugs"))
    return card("Existing binders for these targets", body,
                note=esc(rec.get("not_our_designs")),
                sub=f'{rec.get("n_with_designed")} of {rec.get("n_targets")} '
                    f'targets have a designed binder in the PDB')


def sec_methods(run: dict) -> str:
    meta = run.get("interface_meta") or {}
    params = meta.get("interface_params") or {}
    defs = meta.get("definitions") or {}
    ev = run.get("evidence") or {}
    body = ""
    if params:
        body += "<h4>Interface parameters</h4>" + kv(
            [(k.replace("_", " "), num(v, 2) if isinstance(v, (int, float)) else esc(v))
             for k, v in params.items()])
    if meta.get("hotspot_selection"):
        body += "<h4>Hotspot selection</h4>" + f'<p>{esc(meta["hotspot_selection"])}</p>'
    if defs:
        body += "<h4>Definitions</h4>" + kv(
            [(k.replace("_", " "), esc(v)) for k, v in defs.items()])
    if meta.get("max_complexes_per_target"):
        body += ('<p class="note">Only the top '
                 f'{esc(meta["max_complexes_per_target"])} complexes per target were '
                 'analysed, so this is not a complete epitope survey of the PDB.</p>')
    srcs = ev.get("sources") or []
    if srcs:
        body += ("<h4>Data sources</h4><ul class=\"bullets\">"
                 + "".join(f"<li>{esc(s)}</li>" for s in srcs) + "</ul>")
    if not body:
        return absent_card("Methods", "No stage metadata was found to describe methods.")
    return card("Methods and parameters", body)


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------

CSS = """
:root{--bg:#f6f8fa;--surface:#fff;--border:#d0d7de;--text:#1f2328;
--muted:#656d76;--accent:#0969da;--ok:#1a7f37;--warn:#9a6700;--bad:#cf222e;
--radius:8px;--shadow:0 1px 3px rgba(31,35,40,.08),0 0 1px rgba(31,35,40,.06)}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--text);
font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Inter,Roboto,
"Helvetica Neue",Arial,sans-serif;-webkit-font-smoothing:antialiased}
code,pre,.mono{font-family:ui-monospace,SFMono-Regular,"SF Mono",Menlo,
Consolas,"Liberation Mono",monospace;font-size:.87em}
header.top{position:sticky;top:0;z-index:20;background:rgba(255,255,255,.92);
backdrop-filter:blur(8px);border-bottom:1px solid var(--border)}
.top-in{max-width:1180px;margin:0 auto;padding:14px 24px 0}
.top h1{margin:0;font-size:19px;font-weight:600;letter-spacing:-.01em}
.top .meta{color:var(--muted);font-size:13px;margin:3px 0 0}
nav{display:flex;gap:2px;overflow-x:auto;margin-top:12px}
nav a{padding:8px 12px;font-size:13.5px;color:var(--muted);text-decoration:none;
border-bottom:2px solid transparent;white-space:nowrap}
nav a:hover{color:var(--text)}
nav a.on{color:var(--accent);border-bottom-color:var(--accent);font-weight:550}
main{max-width:1180px;margin:0 auto;padding:22px 24px 80px}
.card{background:var(--surface);border:1px solid var(--border);
border-radius:var(--radius);box-shadow:var(--shadow);padding:18px 20px;margin:0 0 18px}
.card-head{display:flex;align-items:baseline;gap:10px;margin-bottom:6px}
.card h3{margin:0;font-size:16px;font-weight:600}
.card-head .sub{color:var(--muted);font-size:13px}
h4{margin:18px 0 6px;font-size:14px;font-weight:600}
h5{margin:14px 0 4px;font-size:13px;font-weight:600;color:var(--muted);
text-transform:uppercase;letter-spacing:.04em}
p{margin:6px 0}
.note{color:var(--muted);font-size:13.5px;margin:4px 0 12px}
.absent{background:#fff8c5;border:1px solid #d4a72c55;border-radius:6px;
padding:9px 11px;font-size:13.5px;margin:8px 0}
.verdict{border-radius:6px;padding:10px 12px;font-size:14px;margin:4px 0 10px}
.verdict.bad{background:#ffebe9;border:1px solid #ff818266}
.verdict.warn{background:#fff8c5;border:1px solid #d4a72c55}
.verdict.ok{background:#dafbe1;border:1px solid #1a7f3733}
.warn-text{color:var(--bad);font-weight:550}
.banner{background:#ffebe9;border:1px solid #ff818266;border-radius:var(--radius);
padding:12px 14px;margin:0 0 18px;font-size:14px}
.stats{display:flex;flex-wrap:wrap;gap:10px;margin:4px 0 14px}
.stat{flex:1 1 150px;background:var(--bg);border:1px solid var(--border);
border-radius:6px;padding:10px 12px}
.stat-v{display:block;font-size:20px;font-weight:600;letter-spacing:-.02em}
.stat-k{display:block;font-size:12px;color:var(--muted);text-transform:uppercase;
letter-spacing:.04em;margin-top:2px}
.stat-n{display:block;font-size:12px;color:var(--muted);margin-top:3px}
.chart{display:block;width:100%;height:auto;margin:6px 0 4px;overflow:visible}
.svg-lbl{font-size:10.5px;fill:var(--muted)}
.svg-val{font-size:12px;font-weight:600;fill:var(--text);text-anchor:middle}
.svg-num{font-size:10px;fill:var(--muted)}
.fnl rect,.chart g[data-tip] rect{transition:fill-opacity .12s}
.fnl:hover rect,.chart g[data-tip]:hover rect{fill-opacity:1}
.chart [data-tip]{cursor:default}
.tbl-tools{display:flex;align-items:center;gap:10px;margin:8px 0 6px}
.filter{flex:0 1 320px;padding:6px 10px;border:1px solid var(--border);
border-radius:6px;font-size:13.5px;background:var(--surface);color:var(--text)}
.filter:focus{outline:2px solid #0969da33;border-color:var(--accent)}
.count{color:var(--muted);font-size:12.5px}
.tbl-wrap{overflow-x:auto;border:1px solid var(--border);border-radius:6px}
table.grid{border-collapse:collapse;width:100%;font-size:13.5px}
table.grid th{background:var(--bg);text-align:left;padding:8px 10px;
font-weight:600;font-size:12px;text-transform:uppercase;letter-spacing:.03em;
color:var(--muted);cursor:pointer;user-select:none;white-space:nowrap;
border-bottom:1px solid var(--border);position:sticky;top:0}
table.grid th:hover{color:var(--text)}
table.grid th .arrow{margin-left:4px;opacity:.5;font-size:10px}
table.grid td{padding:7px 10px;border-bottom:1px solid #eaeef2;vertical-align:top}
table.grid tbody tr:last-child td{border-bottom:none}
table.grid tbody tr.has-detail{cursor:pointer}
table.grid tbody tr.has-detail:hover{background:#f6f8fa}
table.grid tbody tr.open{background:#ddf4ff55}
tr.detail{display:none}
tr.detail.open{display:table-row}
tr.detail td{background:#fbfcfd;border-bottom:1px solid var(--border)}
.detail-in{padding:6px 2px 10px}
.dim{color:var(--muted)}
.ok-text{color:var(--ok);font-weight:550}
.caption{color:var(--muted);font-size:12.5px;margin:2px 0 10px;max-width:70ch}
.example{border:1px solid var(--border);border-radius:6px;padding:12px 14px;
margin:10px 0;background:#fbfcfd}
.struct{background:#fff;border:1px solid #eaeef2;border-radius:6px}
.pill{display:inline-block;padding:1px 7px;border-radius:99px;font-size:11.5px;
font-weight:550;background:#eaeef2;color:#424a53;white-space:nowrap}
.pill.ok{background:#dafbe1;color:#0f5323}
.pill.warn{background:#fff8c5;color:#7d4e00}
.pill.bad{background:#ffebe9;color:#a40e26}
.kv{margin:6px 0 4px;display:grid;gap:0}
.kv-row{display:grid;grid-template-columns:210px 1fr;gap:12px;padding:5px 0;
border-bottom:1px solid #eaeef2;font-size:13.5px}
.kv-row:last-child{border-bottom:none}
.kv dt{color:var(--muted)}
.kv dd{margin:0;overflow-wrap:anywhere}
.bullets{margin:4px 0 8px;padding-left:20px;font-size:13.5px;color:#424a53}
.bullets li{margin:3px 0}
pre.code{background:var(--bg);border:1px solid var(--border);border-radius:6px;
padding:10px 12px;overflow-x:auto;margin:4px 0 8px;line-height:1.45}
button.copy{margin-left:6px;padding:1px 7px;font-size:11.5px;border-radius:5px;
border:1px solid var(--border);background:var(--surface);color:var(--muted);
cursor:pointer}
button.copy:hover{color:var(--text);border-color:#8c959f}
#tip{position:fixed;z-index:60;pointer-events:none;display:none;max-width:320px;
background:#1f2328;color:#fff;font-size:12.5px;line-height:1.4;padding:6px 9px;
border-radius:6px;box-shadow:0 4px 12px rgba(0,0,0,.18)}
footer{max-width:1180px;margin:0 auto;padding:0 24px 60px;color:var(--muted);
font-size:12.5px}
@media (max-width:720px){.kv-row{grid-template-columns:1fr}
.top-in,main{padding-left:14px;padding-right:14px}}
@media print{header.top{position:static}nav{display:none}
.card{break-inside:avoid;box-shadow:none}tr.detail{display:table-row}}
"""

JS = """
function tipInit(){
  var tip=document.getElementById('tip');
  document.addEventListener('mouseover',function(e){
    var el=e.target.closest('[data-tip]');
    if(!el){tip.style.display='none';return;}
    tip.textContent=el.getAttribute('data-tip');
    tip.style.display='block';
  });
  document.addEventListener('mousemove',function(e){
    if(tip.style.display!=='block')return;
    var x=e.clientX+14, y=e.clientY+14;
    var r=tip.getBoundingClientRect();
    if(x+r.width>window.innerWidth-8)x=e.clientX-r.width-14;
    if(y+r.height>window.innerHeight-8)y=e.clientY-r.height-14;
    tip.style.left=x+'px'; tip.style.top=y+'px';
  });
}
function cellKey(td,type){
  var t=(td.innerText||'').trim();
  if(type==='num'){
    var m=t.replace(/,/g,'').match(/-?\\d+(\\.\\d+)?/);
    return m?parseFloat(m[0]):Number.NEGATIVE_INFINITY;
  }
  return t.toLowerCase();
}
function sortTable(id,col){
  var tb=document.getElementById(id), body=tb.tBodies[0];
  var type=tb.tHead.rows[0].cells[col].getAttribute('data-type')||'text';
  var prev=tb.getAttribute('data-sorted-col'), dir=1;
  if(prev===String(col)) dir=(tb.getAttribute('data-sorted-dir')==='1')?-1:1;
  var groups=[], rows=Array.prototype.slice.call(body.rows);
  for(var i=0;i<rows.length;i++){
    if(rows[i].classList.contains('detail')) continue;
    var g={main:rows[i],detail:null};
    if(i+1<rows.length&&rows[i+1].classList.contains('detail')) g.detail=rows[i+1];
    groups.push(g);
  }
  groups.sort(function(a,b){
    var x=cellKey(a.main.cells[col],type), y=cellKey(b.main.cells[col],type);
    if(x<y)return -dir; if(x>y)return dir; return 0;
  });
  groups.forEach(function(g){ body.appendChild(g.main);
    if(g.detail) body.appendChild(g.detail); });
  tb.setAttribute('data-sorted-col',col);
  tb.setAttribute('data-sorted-dir',dir===1?'1':'-1');
  var ths=tb.tHead.rows[0].cells;
  for(var j=0;j<ths.length;j++){
    var a=ths[j].querySelector('.arrow');
    if(a) a.textContent=(j===col)?(dir===1?'\\u25b2':'\\u25bc'):'';
  }
}
function filterTable(id,q){
  var tb=document.getElementById(id), body=tb.tBodies[0];
  q=(q||'').toLowerCase().trim();
  var shown=0,total=0;
  var rows=Array.prototype.slice.call(body.rows);
  for(var i=0;i<rows.length;i++){
    var r=rows[i];
    if(r.classList.contains('detail')) continue;
    total++;
    var hit=!q||(r.innerText||'').toLowerCase().indexOf(q)>=0;
    r.style.display=hit?'':'none';
    if(hit)shown++;
    var d=rows[i+1];
    if(d&&d.classList.contains('detail')&&!hit){
      d.classList.remove('open'); r.classList.remove('open');
    }
  }
  var c=document.getElementById(id+'-count');
  if(c) c.textContent=q?(shown+' of '+total+' rows'):(total+' rows');
}
function toggleDetail(rid){
  var d=document.getElementById(rid);
  if(!d)return;
  var main=d.previousElementSibling;
  var open=d.classList.toggle('open');
  if(main) main.classList.toggle('open',open);
}
function copyText(btn,txt){
  var done=function(){ var o=btn.textContent; btn.textContent='copied';
    setTimeout(function(){btn.textContent=o;},1100); };
  if(navigator.clipboard&&navigator.clipboard.writeText){
    navigator.clipboard.writeText(txt).then(done,function(){});
  }else{
    var ta=document.createElement('textarea'); ta.value=txt;
    document.body.appendChild(ta); ta.select();
    try{document.execCommand('copy');done();}catch(e){}
    document.body.removeChild(ta);
  }
}
function navInit(){
  var links=Array.prototype.slice.call(document.querySelectorAll('nav a'));
  var secs=links.map(function(a){return document.querySelector(a.getAttribute('href'));});
  function onScroll(){
    var best=0, y=window.scrollY+140;
    for(var i=0;i<secs.length;i++){ if(secs[i]&&secs[i].offsetTop<=y) best=i; }
    links.forEach(function(a,i){ a.classList.toggle('on',i===best); });
  }
  window.addEventListener('scroll',onScroll,{passive:true});
  onScroll();
}
window.addEventListener('DOMContentLoaded',function(){
  tipInit(); navInit();
  document.querySelectorAll('table.grid').forEach(function(t){
    filterTable(t.id,'');
    var s=t.getAttribute('data-sort');
    if(s!==null&&s!=='') sortTable(t.id,parseInt(s,10));
  });
});
"""


def build(run_dir: str | Path, out_path: str | Path, title: str | None = None) -> dict:
    """Render ``run_dir`` to a single self-contained interactive HTML file."""
    run_dir = Path(run_dir)
    if not run_dir.is_dir():
        raise NotADirectoryError(f"{run_dir} is not a directory")
    run = load_run(run_dir)

    disease = ((run.get("disease") or {}).get("disease") or {})
    heading = title or (
        f'indication2binder &mdash; {disease.get("name") or run_dir.name}'
    )
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    sections = [
        ("overview", "Overview", sec_overview(run)),
        ("disease", "Indication", sec_disease(run)),
        ("targets", "Targets", sec_targets(run)),
        ("rubric", "Rubric", sec_rubric(run)),
        ("complexes", "Complex gate", sec_complexes(run)),
        ("epitopes", "Epitopes", sec_epitopes(run)),
        ("specs", "Specifications", sec_specs(run)),
        ("binders", "Existing binders", sec_binders(run)),
        ("methods", "Methods", sec_methods(run)),
        ("provenance", "Provenance", sec_provenance(run)),
    ]

    missing = missing_layers(run)
    banner = ""
    manifest = run.get("run_manifest") or {}
    origins = manifest.get("origins") or {}
    fixtures = [k for k, v in origins.items() if v == "fixture"]
    if fixtures:
        banner += (
            f'<div class="banner"><strong>Fixture data present.</strong> '
            f'Layer(s) {esc(", ".join(fixtures))} are labelled <code>fixture</code>: '
            f'they are placeholders for wiring, not biological results, and nothing '
            f'derived from them is a finding.</div>'
        )

    nav = "".join(
        f'<a href="#{sid}">{esc(label)}</a>' for sid, label, _ in sections
    )
    body = "".join(
        f'<div id="{sid}">{html_block}</div>' for sid, _, html_block in sections
    )

    doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{esc(disease.get("name") or run_dir.name)} &mdash; indication2binder run report</title>
<style>{CSS}</style></head>
<body>
<header class="top"><div class="top-in">
<h1>{heading}</h1>
<p class="meta">Run directory <code>{esc(run_dir)}</code> &middot; rendered
{esc(generated)} &middot; reads persisted run outputs only, no source re-queried</p>
<nav>{nav}</nav></div></header>
<main>{banner}{body}
<footer>
<p>Generated by <code>ind2b_interactive.py</code> (schema v{SCHEMA_VERSION}) from
{len(run.get("_files", []))} files in <code>{esc(run_dir)}</code>. Absent layers are
reported absent rather than illustrated: {esc(", ".join(missing)) or "none"}.</p>
<p>This report describes a computational shortlist. It is not a claim that any
target is causal, that blocking a named partner is therapeutic, or that a binder
exists for any epitope shown.</p>
</footer>
</main>
<div id="tip" role="tooltip"></div>
<script>{JS}</script>
</body></html>
"""
    out_path = Path(out_path)
    out_path.write_text(doc, encoding="utf-8")

    epi = run.get("epitopes")
    return {
        "out": str(out_path),
        "bytes": out_path.stat().st_size,
        "layers": sorted(k for k in run if not k.startswith("_")),
        "missing": missing,
        "n_sections": len(sections),
        "n_epitopes_accepted": (
            int((epi["status"] == "accepted").sum())
            if isinstance(epi, pd.DataFrame) and "status" in epi.columns else None
        ),
        "n_specs": (run.get("specs_manifest") or {}).get("n_specs"),
        "designs_present": bool(run.get("design_run") or run.get("candidates") is not None),
    }


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--run-dir", required=True, type=Path)
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--title", default=None)
    args = p.parse_args(argv)
    out = args.out or Path(f"{args.run_dir.name}_report.html")
    result = build(args.run_dir, out, title=args.title)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
