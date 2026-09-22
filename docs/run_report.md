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

## Example output

`docs/examples/nsclc_run_report.html` is the report this renderer produced for the
non-small cell lung carcinoma run, committed so it can be reviewed without running
the pipeline. GitHub will not execute the page in the file view — download it or use
the Raw link. See `docs/examples/README.md` for how it was generated.

## Sections

| Section | Shows |
|---|---|
| Overview | the funnel as a gate sequence: candidates -> reachable -> has complex -> epitopes -> specs, plus the design-state headline read from the run's own design record |
| Indication | resolved ontology id, the rule that chose it, and the candidates rejected |
| Targets | every pooled target including excluded ones, expanding to score-component bars and evidence |
| Rubric | weights, penalties, count saturation, and the caveats recorded with the ranking |
| Complex gate | admitted complexes per target; for dropped targets, the rejection reasons |
| Epitopes | accepted and rejected, expanding to a chain position track, per-residue burial bars and the residue table |
| Specifications | the emitted target and campaign JSON inline, with copy-to-clipboard hotspots |
| Methods | interface cutoffs, hotspot-selection rule, definitions, data sources |
| Existing binders | published designed binders and antibodies of these targets, with backbone figures and epitope-overlap comparison |
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


## Existing binders (stage 6)

`ind2b stage6` collects binders that **already exist** against the shortlisted
targets and compares their footprint with the epitope this pipeline selected:

```bash
ind2b stage6 --efo-id MONDO_0005233
```

Two kinds are collected, both from the PDB search stage 3 already performed:

- **designed** — a de novo binder, miniprotein or designed scaffold (DARPin,
  designed scFv). These carry no UniProt accession because they do not exist in
  nature, which is why stage 3 set them aside as non-partners. They are the
  closest published analogue to a BindCraft2 output.
- **antibody** — Fab, scFv or nanobody complexes, clinically validated for
  several of these targets.

For each, the target-side footprint is computed with stage 4's interface
definition, mapped into UniProt numbering through the SIFTS alignment, and
compared with the selected epitope (shared positions, Jaccard, fraction of
selected hotspots covered). A 2-D projection of the complex's CA trace is
computed here and persisted, so the report draws the structure from a record
rather than re-reading coordinates.

### Reading the overlap correctly

**Nothing in this stage is a design produced by the pipeline.** Overlap with an
existing binder is evidence that a surface is *bindable*; it says nothing about
the affinity or specificity of a binder not yet made.

A zero overlap has two very different causes, and the record distinguishes them
via `epitope_present_in_construct`:

- the deposited construct **does not contain** the selected epitope — the zero
  is silence about that surface, not evidence against it;
- the epitope **is** present and the binder engages elsewhere — informative.

Scale matters too. `is_minibinder_scale` flags binders at or below
`MINIBINDER_MAX_RESIDUES` (100). A DARPin or scFv footprint is much larger than
anything a mini-binder campaign will produce, so a large scaffold binding a
different surface is a fact about that scaffold, not about the epitope.
