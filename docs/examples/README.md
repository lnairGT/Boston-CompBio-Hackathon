# Example report

`nsclc_run_report.html` (1.6 MB) is the report this renderer produced for
the non-small cell lung carcinoma run (`MONDO_0005233`), committed so the output
can be reviewed without running the pipeline.

GitHub will not execute the page in the file view — **download it and open it
locally**, or use the "Raw" link. It is self-contained: no network access is
needed to view it.

Reproduce it with:

```bash
ind2b stage0 --indication "non-small cell lung carcinoma"
for s in stage1 stage2 stage3 stage4 stage5; do
  ind2b $s --efo-id MONDO_0005233
done
python scripts/ind2b_interactive.py --run-dir runs/MONDO_0005233 \
    --out docs/examples/nsclc_run_report.html
```

## What this run contains

100 candidate targets from 12,475 disease associations; 33 reachable by an
extracellular binder; 24 with a qualifying experimental complex; 59 accepted
epitopes; 24 BindCraft2 specifications.

**No binder was designed.** The run stops at specifications, and the report
says so in its second section. `design_run` and `candidates` are absent and are
reported absent rather than illustrated.

`run_manifest.json` is the provenance sidecar for that run: every layer is
labelled `real`, having been fetched live from Open Targets, UniProt, RCSB PDB,
ChEMBL and ClinicalTrials.gov. It sits here beside the report as an example of
the label format; the renderer reads it from the top level of a run directory.
