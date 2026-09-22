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


## Stage 6 - the binders the pipeline actually designed

Stage 6 reads a **BindCraft2 campaign's own output**. It is adaptive, not
hardcoded: point it at one or more `project_folder`s, or let it look under
`results/` and `campaigns/` in the run directory.

```bash
ind2b stage6 --efo-id MONDO_0005233
ind2b stage6 --efo-id MONDO_0005233 --project-folder /path/to/campaign
```

It reads the files BindCraft2 writes - `3_Ranked/!_Ranked.csv` (the record of
what was accepted), `2_Refolded/!_Refolded.csv`, `1_Trajectories/!_Trajectories.csv`,
`summary.csv` and `campaign_metadata.json` - and for each accepted design records
its sequence, every metric column the campaign wrote, its structure file, and
whether it engaged the hotspots the stage 5 specification requested. Both sides of
that comparison are in the written structure's author numbering, which is what the
specification asked for, so they compare directly.

If no campaign output exists, stage 6 says so and the report's **Designed
binders** section states that no binder has been designed. Nothing is
substituted for it.

### What stage 6 refuses to let you misread

These follow the upstream output documentation:

- **Acceptance is not affinity.** An accepted design passed the configured
  computational filters; that establishes neither affinity, specificity nor
  experimental stability.
- **A blank is not a zero.** Missing metrics stay missing; `failed_filters`
  distinguishes "not measured" from a poor score.
- **Rank is not probability.** `rank` is a position in one ranking.
- **Some attempts ran an easier task.** Where `autotuned` names `initial_guess`,
  `target_flexibility`, `validation_model` or a raised `design_recycles`, that
  attempt ran on the desperation ladder - a weaker problem than the campaign
  specified. Those designs are flagged in the record and in the report, because
  their scores are not comparable with the rest.
- **A missing stage folder is informative.** It appears only once its first
  result is written, so its absence means the campaign did not get that far -
  not that the target is undesignable.

## Reference binders - a comparison, not a pipeline stage

`ind2b reference-binders` mines the PDB for binders that **already exist**
against the shortlisted targets (de novo designs, miniproteins, DARPins,
designed scFvs, antibodies) and compares their footprint with the selected
epitope.

It is deliberately **not** part of the pipeline and not wired into `all`:
other people's binders are not this pipeline's output, and a report that
showed them by default would invite exactly that confusion. Run it when you
want the comparison - it is what the committed example report shows, because
no campaign has been run for that example.

Its output lands in `reference_binders.json`, and the report renders it in a
**Reference binders (not designed here)** section, separate from Designed
binders and labelled as reference molecules throughout.

### Reading the overlap correctly

Overlap with an existing binder is evidence that a surface is *bindable*; it
says nothing about a binder not yet made. A zero overlap has two very different
causes, and the record distinguishes them via `epitope_present_in_construct`:
the deposited construct may not contain the selected epitope at all (silence
about that surface), or the epitope is present and the binder engages elsewhere
(informative). `is_minibinder_scale` separates mini-binder analogues from large
designed scaffolds, whose footprints are far larger than a campaign produces.

## 3-D viewers

Structure figures are rotatable, and there are two renderers. Both are
self-contained: nothing is fetched at view time.

| `--viewer` | What you get | Example report size |
| --- | --- | --- |
| `auto` (default) | cartoon rendering via 3Dmol.js where a trimmed structure exists, backbone trace otherwise | ~4.6 MB |
| `rich` | as `auto`, but fails loudly if the vendored library is missing | ~4.6 MB |
| `light` | backbone trace only, no library inlined | ~1.8 MB |

```bash
python scripts/ind2b_interactive.py --run-dir runs/MONDO_0005233 --viewer light
```

**Cartoon rendering** uses 3Dmol.js, vendored at `scripts/vendor/3Dmol-min.js`
(0.54 MB, BSD-3-Clause; `VERSION.json` records the version, source URL and
SHA-256). It is vendored rather than loaded from a CDN so a report renders
reproducibly and views offline. The library is inlined only when a rich viewer
was actually emitted, so a run with nothing to show in 3-D does not carry half a
megabyte for nothing. Target chain is grey, binder purple, epitope blue, hotspots
red with spheres on their Cα.

**The backbone trace** is the fallback and needs no library: a canvas drawing of a
persisted principal-axis frame of the complex's CA coordinates. It is what you see
when the library is absent, when `--viewer light` is used, **or when the browser
provides no WebGL** - in that last case the report says so in place of the cartoon
view rather than showing an empty box. Printing falls back to a flat SVG of the
same geometry.

Neither is a publication figure: a WebGL render is not deterministic across GPUs
and browsers. For figures, open the deposited entry in PyMOL, ChimeraX or Mol*.

### Structures for the viewer

`ind2b reference-binders` trims each projected complex to the target chain plus
the one binder chain whose interface was measured, writing
`viewer_structures/<SYMBOL>_<ENTRY>_<CHAIN>.pdb`. Only that pair: a deposited
entry often holds several copies of the binder, and showing the non-contacting
ones would imply contacts that are not there.

Output is PDB because in-browser viewers parse it most consistently, which is
safe only for small complexes - legacy PDB caps atom serials at 99,999 and chain
ids at one character. `structure_export` refuses to write a file that would
overflow either, rather than emitting one a viewer would misread.
