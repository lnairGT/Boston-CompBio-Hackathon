#!/usr/bin/env python3
"""Run a BindCraft2 campaign on Modal from an ind2b stage-5 specification.

This sits BETWEEN stage 5 and stage 6 and is deliberately NOT a stage of the pipeline.
ind2b's stated contract is that it "never submits a design job -- stage 5 writes
specification files and stops", and that boundary is worth keeping: spec generation is
cheap, reproducible and offline, while design execution is expensive, GPU-bound and bills
someone.

    stage 5  ->  BindCraft2 spec  ->  [this runner]  ->  campaign output  ->  stage 6

Usage
-----
    python run_campaign.py --spec-dir runs/MONDO_0004980/stage5_specs \\
                           --campaign <NAME> --timeout 3600 [--dry-run]

What it does not promise
------------------------
A finished campaign. BindCraft2 runs until it reaches the requested number of accepted
binders, with **no limit on design attempts**, so there is no bounded runtime to promise.
This runner bounds the wall clock and harvests whatever the campaign wrote, recording that
the run was time-bounded. Stage 6 already reads that correctly -- it discovers campaigns by
their records rather than by folder name, so a campaign that accepted nothing is still
found, and a missing stage folder means the campaign did not get that far rather than that
the target is undesignable.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

# Verified on this account 2026-09-22: BC2 at commit d5bae16, jax 0.11.2 reporting
# [CudaDevice(id=0)] on A100, AF2 params staged on a volume.
IMAGE_ID = "im-uj023rwPkhbFApR2O9hEXa"
AF2_VOLUME = "bindcraft2-af2-params"
AF2_MOUNT = "/weights"
AF2_PARAMS = "/weights/af2/alphafold"


def stage_campaign(spec_dir: Path, campaign: str, workdir: Path,
                   remote_root: str = "/work") -> dict:
    """Copy a stage-5 spec into a layout that satisfies BindCraft2's path resolution.

    BC2's rules (docs/reference.md):

      * ``target_path`` written in a JSON file -> resolved from THAT file's directory
      * ``project_folder``                     -> from the directory the command runs in,
                                                  unless an absolute path is given

    Stage 5 writes ``target_path: "structures/<name>.pdb"`` into
    ``settings/target/<name>.json`` while the structure lands at
    ``<stage5>/structures/<name>.pdb``. Relative to the target JSON's own directory that
    string does not resolve. Rather than depend on which directory BC2 happens to resolve
    from, this rewrites ``target_path`` to an absolute in-job path, which is unambiguous
    under every rule above. Without that, a missing structure would surface well into a
    paid GPU job instead of here.
    """
    spec_dir = spec_dir.resolve()
    camp_src = spec_dir / "campaigns" / f"{campaign}.json"
    targ_src = spec_dir / "settings" / "target" / f"{campaign}.json"
    for p, what in ((camp_src, "campaign"), (targ_src, "target")):
        if not p.exists():
            raise SystemExit(f"no {what} spec at {p}")

    (workdir / "settings" / "target").mkdir(parents=True, exist_ok=True)
    (workdir / "structures").mkdir(parents=True, exist_ok=True)

    targ = json.loads(targ_src.read_text())
    staged = []
    for entry in targ.get("targets") or []:
        rel = entry.get("target_path")
        if not rel:
            continue
        src = (spec_dir / rel).resolve()   # stage 5 records it relative to its own root
        if not src.exists():
            raise SystemExit(
                f"structure referenced by the spec is missing: {rel} (looked in {src}). "
                "Stage 5 records target_path relative to its output root; if that "
                "convention changed this runner needs updating rather than guessing."
            )
        shutil.copy2(src, workdir / "structures" / src.name)
        # Absolute under the root the job will actually unpack into. This was a real
        # failure, not a hypothetical: a hardcoded /work while the job extracted into
        # /work/job produced "target 'X': no file at /work/structures/X.pdb" -- after BC2
        # had already reached the GPU and written a compile cache, so it read as a crashed
        # campaign rather than a bad path.
        entry["target_path"] = f"{remote_root.rstrip('/')}/structures/{src.name}"
        staged.append(src.name)

    (workdir / "settings" / "target" / f"{campaign}.json").write_text(json.dumps(targ, indent=2))

    camp = json.loads(camp_src.read_text())

    # INLINE the target definition into the campaign file.
    #
    # Stage 5 writes two linked files: a campaign with `"target": "<NAME>"` and a separate
    # settings/target/<NAME>.json. BC2 cannot join them -- a bare `target` name is resolved
    # ONLY against the presets the package ships, and a local settings/target directory is
    # not searched. Running stage 5's output unmodified fails with:
    #
    #   campaign refused:
    #   unknown target 'IL4R_1IAR_B'; this package ships dynorphin_a, hIL2R, hIL7RA,
    #   hPD1, hPDL1, mPDL1
    #
    # The documented way to supply your own protein is the inline `targets` array
    # (reference.md: "targets | List of target objects | Supply your own proteins..."), so
    # the target rows are merged in and the unresolvable name dropped.
    if targ.get("targets"):
        camp["targets"] = targ["targets"]
        camp.pop("target", None)
        # Do NOT carry the target file's `description` up to campaign level. BC2 validates
        # campaign keys strictly and refuses it:
        #   unrecognized campaign settings: description (did you mean desperation?)
        # The text is provenance for a human, so it is preserved alongside the job rather
        # than smuggled into a settings file that rejects it.

    camp["project_folder"] = f"{remote_root.rstrip('/')}/out/{camp.get('campaign_name') or campaign}"
    # An existing folder is continued by default; refusing a non-empty one keeps an
    # accidental re-submit from becoming a silent resume.
    camp.setdefault("resume", False)
    (workdir / f"{campaign}.json").write_text(json.dumps(camp, indent=2))

    return {
        "campaign_file": f"{campaign}.json",
        "project_folder": camp["project_folder"],
        "campaign_name": camp.get("campaign_name") or campaign,
        "structures": staged,
        "n_designs": camp.get("number_of_final_designs"),
        "binder_lengths": camp.get("binder_lengths"),
        "hotspots": [e.get("hotspots") for e in (targ.get("targets") or [])],
        "chains": [e.get("chains") for e in (targ.get("targets") or [])],
    }


def build_command(staged: dict) -> str:
    """The job command. Two parts are load-bearing rather than defensive.

    ``BINDCRAFT_AF2_PARAMS`` points at the staged 5.3 GB parameter set; without it BC2
    lazily downloads them on first campaign, which on an account with no job-time egress
    hangs rather than failing fast.

    The ``jax.devices()`` assertion runs BEFORE the campaign and aborts on CPU fallback.
    BC2's own fan-out log names every GPU even when jax has fallen back to the CPU, because
    that line is read off ``CUDA_VISIBLE_DEVICES`` rather than off the devices jax opened.
    So the log is not evidence, and without this check a campaign could burn its entire
    timeout ~100x slower with nothing in the output explaining why.
    """
    return (
        "set -eu; cd /work; "
        "python3 -c \"import jax,sys; d=jax.devices(); print('devices:',d); "
        "sys.exit(0 if any(x.platform=='gpu' for x in d) else 'ABORT: jax fell back to CPU')\"; "
        f"export BINDCRAFT_AF2_PARAMS={AF2_PARAMS}; "
        f"bindcraft design {staged['campaign_file']} 2>&1 | tail -60"
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--spec-dir", required=True, type=Path)
    ap.add_argument("--campaign", required=True)
    ap.add_argument("--timeout", type=int, default=3600,
                    help="wall-clock bound in seconds; a BC2 campaign has no natural end")
    ap.add_argument("--workdir", type=Path, default=Path("bc2_job"))
    ap.add_argument("--remote-root", default="/work",
                    help="absolute directory the job unpacks into. target_path and "
                         "project_folder are written under it, so it must match where the "
                         "files actually land -- a mismatch surfaces as a missing structure "
                         "several seconds into a paid GPU job.")
    ap.add_argument("--dry-run", action="store_true", help="stage and print; submit nothing")
    args = ap.parse_args()

    staged = stage_campaign(args.spec_dir, args.campaign, args.workdir,
                            remote_root=args.remote_root)
    cmd = build_command(staged)
    print(json.dumps({"staged": staged, "command": cmd, "image": IMAGE_ID,
                      "volumes": {AF2_MOUNT: AF2_VOLUME}, "timeout_s": args.timeout}, indent=2))
    if args.dry_run:
        print("\n[dry-run] nothing submitted", file=sys.stderr)
        return 0

    try:
        import modal
    except ImportError:
        raise SystemExit("pip install modal, or submit the printed command yourself")

    app = modal.App.lookup("bc2-campaigns", create_if_missing=True)
    sb = modal.Sandbox.create(
        image=modal.Image.from_id(IMAGE_ID), app=app, gpu="A100",
        volumes={AF2_MOUNT: modal.Volume.from_name(AF2_VOLUME)},
        timeout=args.timeout + 300, workdir="/work",
    )
    try:
        for f in sorted(args.workdir.rglob("*")):
            if f.is_file():
                with sb.open(f"/work/{f.relative_to(args.workdir)}", "wb") as fh:
                    fh.write(f.read_bytes())
        p = sb.exec("bash", "-lc", cmd, timeout=args.timeout)
        print(p.stdout.read())
        err = p.stderr.read()
        if err.strip():
            print("STDERR:", err[-2000:], file=sys.stderr)

        # Harvest whatever exists: a time-bounded campaign is a legitimate partial result.
        listing = sb.exec("bash", "-lc",
                          f"find {staged['project_folder']} -type f 2>/dev/null | head -300")
        files = [l for l in listing.stdout.read().splitlines() if l.strip()]
        out_local = args.workdir / "harvested"
        for rf in files:
            dest = out_local / Path(rf).relative_to(staged["project_folder"])
            dest.parent.mkdir(parents=True, exist_ok=True)
            try:
                with sb.open(rf, "rb") as fh:
                    dest.write_bytes(fh.read())
            except Exception as e:      # a partial write mid-campaign is expected
                print(f"  skip {rf}: {e}")
        (args.workdir / "run_manifest.json").write_text(json.dumps({
            "time_bounded": True, "timeout_s": args.timeout, "staged": staged,
            "harvested_files": len(files),
            "note": ("Wall clock was bounded by this runner. BindCraft2 has no natural "
                     "stopping point, so an absence of accepted designs means the budget "
                     "ran out, NOT that the target is undesignable."),
        }, indent=2))
        print(f"\nharvested {len(files)} file(s) -> {out_local}")
    finally:
        sb.terminate()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
