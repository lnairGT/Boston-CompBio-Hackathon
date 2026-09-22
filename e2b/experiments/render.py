"""Render a plan dict as markdown for the demo's "next experimental test" panel.

Written for a bench scientist reading it once, quickly: what gets made, what gets
measured, what the controls rule out, and -- visibly, not buried -- what result kills
the candidate.
"""

from __future__ import annotations

from typing import Any

_STATUS_LINE = {
    "plan_complete": "Plan complete: binding-validation and mechanism-specific functional routes both derived.",
    "partial_no_functional_readout": "PARTIAL: binding route derived; the functional readout could not be derived from the inputs.",
    "partial_no_candidate_sequence": "PARTIAL: no candidate sequence was supplied, so the construct-liability screen did not run.",
    "blocked_target_not_accessible": "BLOCKED: the target has no extracellular face, so a soluble binder cannot engage it in a cell.",
    "blocked_no_target_biology": "BLOCKED: target biology was unavailable, so no plan is derivable.",
    "blocked_direction_not_addressable": "BLOCKED: the requested direction of modulation is not achievable with this modality.",
}


def _bullets(items: list[Any], indent: str = "") -> str:
    return "\n".join(f"{indent}- {i}" for i in items if i) or f"{indent}- (none recorded)"


def _starting_points(sps: list[dict[str, Any]]) -> str:
    if not sps:
        return ""
    rows = [
        f"- **{s['name']}** — {s['range_text']}  \n  *Why:* {s['reasoning']}  \n  *Status:* `{s['status']}`"
        for s in sps
    ]
    return "\n".join(rows)


def _controls(cs: list[dict[str, Any]]) -> str:
    if not cs:
        return "- (none recorded)"
    out = []
    for c in cs:
        avail = c.get("available")
        tag = {True: "available", False: "not available", None: "availability unconfirmed"}[avail]
        line = f"- **{c['name']}** ({c['kind']}, {tag}) — {c['purpose']}  \n  *If it fails:* {c['interpretation_if_it_fails']}"
        if c.get("availability_note"):
            line += f"  \n  *Note:* {c['availability_note']}"
        out.append(line)
    return "\n".join(out)


def _assay_section(a: dict[str, Any]) -> str:
    L: list[str] = [f"### {a['title']}", "", f"**Format.** {a['format']}", ""]
    L += [f"**Why this format.** {a['rationale']}", ""]
    if a.get("orientation"):
        L += [f"**Orientation / design.** {a['orientation']}", ""]
    if a.get("orientation_rationale"):
        L += [f"**Why this orientation.** {a['orientation_rationale']}", ""]
    if a.get("readout"):
        L += [f"**Readout.** {a['readout']}", ""]
    if a.get("readout_reports_mechanism_because"):
        L += [f"**Why that readout reports the mechanism.** {a['readout_reports_mechanism_because']}", ""]
    if a.get("system"):
        L += [f"**System.** {a['system']}", ""]
    if a.get("system_requirement"):
        L += [f"**System requirement.** {a['system_requirement']}", ""]
    if a.get("resolvable_range"):
        L += [
            f"**What it can actually resolve.** {a['resolvable_range']}",
            "",
            f"*Reasoning:* {a['resolvable_range_reasoning']}",
            "",
        ]
    if a.get("starting_points"):
        L += ["**Starting points (to be optimised — none of these is a validated value):**", "", _starting_points(a["starting_points"]), ""]
    L += ["**Controls:**", "", _controls(a["controls"]), ""]
    L += ["**Reports:**", _bullets(a.get("reports", [])), ""]
    L += ["**Does not report:**", _bullets(a.get("does_not_report", [])), ""]
    if a.get("limitations"):
        L += ["**Limitations:**", _bullets(a["limitations"]), ""]
    if a.get("derived_from"):
        L += ["*Derived from:* " + "; ".join(a["derived_from"]), ""]
    return "\n".join(L)


def _route_section(r: dict[str, Any]) -> str:
    subject = "Designed binder" if r["subject"] == "designed_binder" else "Target antigen reagent"
    L = [f"### {subject}", "", f"**Expression host.** {r['expression_host']}", "", f"*Why:* {r['host_rationale']}", ""]
    L += ["**Construct.**", _bullets(r.get("construct_elements", [])), ""]
    if r.get("construct_rationale"):
        L += ["**Construct reasoning.**", _bullets(r["construct_rationale"]), ""]
    L += ["**Purification.**", _bullets(r.get("purification_steps", [])), ""]
    L += ["**Release QC — nothing downstream is believed until these pass.**", _bullets(r.get("release_qc", [])), ""]
    if r.get("qc_gates"):
        L += ["**QC gates.**", _bullets(r["qc_gates"]), ""]
    if r.get("risks"):
        L += ["**Sequence-derived risks.**", ""]
        for k in r["risks"]:
            pos = f" (positions {k['positions']})" if k.get("positions") else ""
            L += [
                f"- **{k['risk']}** [{k['severity']}] — evidence: {k['evidence']}{pos}  \n"
                f"  *Consequence:* {k['consequence']}  \n  *Mitigation:* {k['mitigation']}"
            ]
        L += [""]
    if r.get("starting_points"):
        L += ["**Starting points.**", "", _starting_points(r["starting_points"]), ""]
    if r.get("notes"):
        L += ["**Notes.**", _bullets(r["notes"]), ""]
    return "\n".join(L)


def render_markdown(plan: dict[str, Any], *, title: str | None = None) -> str:
    """Render a plan dict (output of ``plan_experiments``) as markdown."""
    p = plan
    tp = p.get("target_profile") or {}
    L: list[str] = []
    L.append(title or "# Next experimental test")
    L.append("")
    L.append(f"*Plan `{p['plan_id']}` · generated {p['generated_at']} · module `{p['module_version']}`*")
    L.append("")
    L.append(f"**Status.** {_STATUS_LINE.get(p['status'], p['status'])}")
    L.append("")

    # Provenance block: what this plan was built from, and what each input's origin was.
    L.append("| Input | Value |")
    L.append("|---|---|")
    L.append(f"| Indication (verbatim from brief) | {p.get('raw_indication') or '—'} |")
    L.append(f"| Desired effect | {p.get('desired_effect') or '—'} |")
    L.append(f"| Direction (derived) | {p.get('direction') or '—'} |")
    L.append(f"| Modality | {p.get('modality') or '—'} |")
    L.append(f"| Target | {tp.get('gene_symbol') or p.get('target_id') or '—'} ({tp.get('protein_accession') or 'no accession'}) |")
    L.append(f"| Target biology origin | `{tp.get('origin') or '—'}` |")
    L.append(f"| Candidate | {p.get('candidate_id') or 'none supplied'} (origin `{p.get('candidate_origin') or '—'}`) |")
    L.append(f"| Candidate length | {p.get('candidate_sequence_length') if p.get('candidate_sequence_length') is not None else 'not determined'} |")
    L.append(f"| Assay archetype | {tp.get('archetype') or '—'} |")
    L.append(f"| Mechanism class | {tp.get('mechanism_class_label') or tp.get('mechanism_class') or '—'} |")
    L.append("")

    if p.get("failure_detail"):
        L += ["> **Why this plan is not complete.** " + p["failure_detail"], ""]
    if p.get("missing_inputs"):
        L += ["**Missing inputs.**", _bullets(p["missing_inputs"]), ""]
    if p.get("provisional_claims"):
        L += ["**Provisional — do not report as findings.**", _bullets(p["provisional_claims"]), ""]

    L += ["## How the plan was derived", ""]
    L += [
        "The assay formats below are not a template. They follow two facts read from the "
        "target's own annotation:",
        "",
        f"- **Archetype: `{tp.get('archetype')}`** — " + "; ".join(tp.get("archetype_basis") or ["no basis recorded"]),
        f"- **Mechanism class: `{tp.get('mechanism_class')}`** — " + "; ".join(tp.get("mechanism_class_basis") or ["no basis recorded"]),
        "",
    ]
    if tp.get("engageable_regions"):
        L += ["Regions a soluble binder could engage (from UniProt topology): " + str(tp["engageable_regions"]), ""]

    if p.get("production_routes"):
        L += ["## 1. Protein production", ""]
        for r in p["production_routes"]:
            L.append(_route_section(r))

    order = [
        ("binding_primary", "## 2. Binding validation"),
        ("binding_orthogonal", "## 3. Orthogonal binding confirmation"),
        ("functional_moa", "## 4. Mechanism-specific functional assay"),
        ("specificity", "## 5. Specificity and selectivity"),
    ]
    for role, heading in order:
        a = next((x for x in p.get("assays", []) if x["role"] == role), None)
        if a:
            L += [heading, "", _assay_section(a)]

    if p.get("decision_tree"):
        L += ["## 6. Pre-registered go / no-go criteria", ""]
        L += [
            "These are fixed before the experiment runs. The point of writing them down now is that "
            "the outcome cannot be reinterpreted afterwards. Every threshold below is a "
            "**pre-registered default that the responsible scientist must review and sign off**, not a "
            "validated cut-off.",
            "",
        ]
        for d in p["decision_tree"]:
            L += [
                f"### {d['criterion_id']} — {d['stage']}",
                "",
                f"**Question.** {d['question']}",
                "",
                f"**Pass.** {d['pass_rule']}",
                "",
                f"**KILL.** {d['kill_rule']}",
                "",
                f"*Basis:* {d['rule_basis']}",
                "",
                f"*On pass:* {d['on_pass'] or '—'}  \n*On kill:* {d['on_kill'] or '—'}"
                + (f"  \n*If ambiguous:* {d['on_ambiguous']}" if d.get("on_ambiguous") else ""),
                "",
                f"*Threshold status:* `{d['threshold_status']}`",
                "",
            ]

    L += ["## 7. Scope and what this does not establish", "", p.get("scope_statement", ""), ""]
    L += [_bullets(p.get("what_this_does_not_establish", [])), ""]
    return "\n".join(L)
