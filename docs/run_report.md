# Run report

`scripts/ind2b_interactive.py` renders one run directory as a single
interactive HTML file.

```bash
python scripts/ind2b_interactive.py --run-dir runs/MONDO_0005233 \
    --out nsclc_report.html
```

Only `pandas` and the standard library are needed. CSS, JavaScript and every
chart are inlined and no asset is fetched at view time, so the file opens from
a laptop with no network. Charts are hand-emitted SVG rather than raster
images, so they stay sharp at any zoom and can carry their own tooltips.

## Sections

| Section | Shows |
|---|---|
| Overview | the funnel as a gate sequence: candidates -> reachable -> has complex -> epitopes -> specs |
| Design status | whether anything was actually designed, and what the run does not establish |
| Indication | resolved ontology id, the rule that chose it, and the candidates rejected |
| Targets | every pooled target including excluded ones, expanding to score-component bars and evidence |
| Rubric | weights, penalties, count saturation, and the caveats recorded with the ranking |
| Complex gate | admitted complexes per target; for dropped targets, the rejection reasons |
| Epitopes | accepted and rejected, expanding to a chain position track, per-residue burial bars and the residue table |
| Specifications | the emitted target and campaign JSON inline, with copy-to-clipboard hotspots |
| Methods | interface cutoffs, hotspot-selection rule, definitions, data sources |
| Provenance | every layer with its origin label and source file |

## What it will not do

The renderer reads only what the run persisted. It never re-queries a source
and never fills a gap with a plausible value:

- an absent layer is reported **absent**, not illustrated with a placeholder;
- a null renders as an em dash, never as zero. `str(x or "")` is avoided
  deliberately: `float('nan')` is truthy, so that idiom turns a missing cell
  into the literal string `nan` - a missing value disguised as data;
- targets and epitopes that were **excluded** stay in the tables with their
  reason, because a filtered-out row is a recorded decision, not a missing one;
- it states that no binder has been designed unless a `design_run.json`
  exists and reached `succeeded`. Emitting a specification is not designing
  a binder.

## Provenance labels

Write a top-level `run_manifest.json` to label where the data came from:

```json
{"origin": "real", "origins": {"epitopes": "real"}}
```

Layers left unlabelled are shown as unlabelled rather than assumed to be real
results. Use `fixture` for placeholder data and the report raises a banner
saying nothing derived from it is a finding.

A matching skill, `ind2b-run-report`, wraps this renderer so an agent can
build the report in one call.
