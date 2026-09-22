# BindCraft2 campaign runner

Takes an `ind2b` stage-5 specification to a BindCraft2 campaign on Modal, and leaves the
`project_folder` that stage 6 reads.

```
stage 5  ->  BindCraft2 spec  ->  [this]  ->  campaign output  ->  stage 6
```

```bash
python tools/bc2_runner/run_campaign.py \
    --spec-dir runs/MONDO_0004980/stage5_specs \
    --campaign IL4R_1IAR_B \
    --timeout 1800 [--dry-run]
```

## Why it lives in `tools/` and not as a stage

The pipeline's contract is that it *"never submits a design job — stage 5 writes
specification files and stops"*. That boundary is worth keeping: spec generation is cheap,
offline and reproducible, while design execution is GPU-bound, unbounded in time and bills
someone. This is a consumer of stage 5's output, not a stage of the pipeline, and
`ind2b` does not import it.

## What it fixes, and why each fix exists

Stage 5's output is not directly runnable by BindCraft2. Each of these was a real refusal
on a GPU, not a hypothetical, and each has a test:

**1. A bare `target` name cannot be resolved.** Stage 5 writes
`{"target": "IL4R_1IAR_B"}` in the campaign plus a separate
`settings/target/IL4R_1IAR_B.json`. BC2 resolves a bare name *only* against the presets it
ships and never searches a local target directory:

```
campaign refused:
unknown target 'IL4R_1IAR_B'; this package ships dynorphin_a, hIL2R, hIL7RA, hPD1, hPDL1, mPDL1
```

The runner inlines the `targets` array into the campaign file, which is BC2's documented way
to supply your own protein, and drops the unresolvable name.

**2. `description` is rejected at campaign level.** It is valid in a target file; BC2
validates campaign keys strictly. The text is kept in `PROVENANCE.json` beside the job
instead of in a settings file that refuses it.

**3. `target_path` is ambiguous as written.** Stage 5 records it relative to its own output
root, but BC2 resolves it *from the directory of the JSON file it appears in*. The runner
rewrites it absolute under `--remote-root`, which is unambiguous under every one of BC2's
resolution rules.

All three are filed upstream as an issue against stage 5, since fixing them there would make
the specs runnable without a wrapper. The runner does the conversion so it can consume
today's output unchanged.

## Two guards that are load-bearing

**It asserts `jax.devices()` shows a CUDA device before the campaign starts.** BC2's per-GPU
fan-out log names every card even when jax has silently fallen back to the CPU, because that
line is read off `CUDA_VISIBLE_DEVICES` rather than off the devices jax opened. So the log is
not evidence, and without this check a campaign can burn its whole budget ~100× slower with
nothing in the output saying why.

**It points `BINDCRAFT_AF2_PARAMS` at pre-staged parameters.** Left to itself BC2 lazily
downloads 5.3 GB on first campaign. On an account with no egress at job time that hangs
rather than failing fast.

## Harvesting: the failure that cost two campaigns

Read this before running anything expensive. Two separate campaigns completed, wrote their
output, and had it **silently discarded** — 76 files and 23 files. Nothing in the exit code, the
logs or the campaign's own output indicated a problem.

Two rules, and the second is the one that is easy to get wrong:

**1. `project_folder` must sit under `./out/`.** That is the only directory collected from a
job. The runner's default already does this; do not point `project_folder` elsewhere.

**2. Do not name an output filter unless you are certain of what it is relative to.** On this
harness the filter is interpreted **relative to the `./out/` root**, so asking for `out/` means
`./out/out/` and matches nothing — every real file is dropped. Omitting the filter returns
everything under `./out/`, which is what you want.

```python
# WRONG — silently keeps nothing, job still exits normally
submit_job(..., outputs=["out/"])

# RIGHT — everything under ./out/ comes back
submit_job(...)
```

The loss surfaces only as a note on the completion record along the lines of
*"N unselected file(s) under ./out/ were not kept"*. If you see that, the files are already gone:
the working directory does not persist between jobs, so a follow-up job cannot retrieve them.

Related: a campaign stopped by its wall-clock bound exits **124** with results already on disk.
A non-zero exit here does not mean "no output".

## What it does not promise

A finished campaign. BindCraft2 runs until it reaches the requested number of accepted
binders, **with no limit on design attempts**, so there is no bounded runtime to promise.
The runner bounds the wall clock, harvests whatever was written, and records
`time_bounded: true` in `run_manifest.json`.

**An absence of accepted designs therefore means the budget ran out — not that the target is
undesignable.** Stage 6 already reads a partial campaign correctly: it discovers campaigns by
their records rather than by folder name, so one that accepted nothing is still found, and a
missing stage folder means the campaign did not get that far.

## Test coverage — read this before trusting it

```bash
python -m pytest tools/bc2_runner -q      # 10 passed, offline
```

**Covered:** the staging logic — target inlining, key rejection, absolute path rewriting
under a configurable root, structure copying, failing *locally* when a referenced structure
is missing, `resume` defaulting to false, and preservation of stage 5's scientific settings
(binder lengths, design count, trajectory cap, modality). Also validated against two real
stage-5 specs (`IL4R_1IAR_B`, `TNFSF4_2HEV_F`).

**Not covered:** `main()`'s Modal submission and harvesting path has never been executed. The
campaign verified for this PR was submitted through a different mechanism with the same
command string, so the staged spec is confirmed to pass BC2's validation, but these
particular submission lines are unexercised. Treat `--dry-run` output as the trustworthy
part and the submission as needing a first real run.

## Environment

Verified 2026-09-22 against BC2 commit `d5bae16`, `jax 0.11.2`, `[CudaDevice(id=0)]` on an
A100-40GB, with AF2 parameters on a Modal volume. Image and the six documented silent-failure
modes are recorded in the project's Modal provider notes.

**Licensing.** BindCraft2 is source-available and **hosting-restricted** — free for internal
and commercial use on your own infrastructure, but offering it, or a service deriving its core
functionality, to third parties as a hosted or managed offering needs a separate licence from
the Pacesa Lab. Running this yourself is unrestricted; running it on someone else's behalf is
the restricted case.
