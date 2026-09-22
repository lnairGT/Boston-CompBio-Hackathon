"""Tests for the stage-5 -> BindCraft2 staging step.

Each test below corresponds to a way a real campaign was refused on a GPU. They exist so
the next person does not rediscover these by burning a job: BC2 validates the spec only
after it has loaded, reached the accelerator and written a compile cache, so a malformed
input reads like a campaign that started and crashed.

Offline: no network, no Modal, no GPU. The submission path in ``main()`` is deliberately
NOT covered here -- see the README.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from run_campaign import build_command, stage_campaign


def make_spec(tmp_path: Path, *, structure: bool = True, name: str = "IL4R_1IAR_B") -> Path:
    """Build a stage-5 spec directory in the exact shape ind2b writes."""
    spec = tmp_path / "stage5_specs"
    (spec / "settings" / "target").mkdir(parents=True)
    (spec / "campaigns").mkdir()
    (spec / "structures").mkdir()

    (spec / "settings" / "target" / f"{name}.json").write_text(json.dumps({
        "description": f"{name} - interface with IL4 from PDB 1IAR chain B, 830 A^2 buried.",
        "targets": [{
            "name": name,
            "target_path": f"structures/{name}.pdb",   # relative to the stage-5 ROOT
            "chains": "B",
            "hotspots": "127,67,69,126,68,13",
        }],
    }))
    (spec / "campaigns" / f"{name}.json").write_text(json.dumps({
        "target": name,                                 # a bare NAME -- BC2 cannot resolve it
        "modality": "binder",
        "campaign_name": name.lower(),
        "binder_lengths": [55, 120],
        "number_of_final_designs": 4,
        "project_folder": f"results/{name.lower()}",
        "max_trajectories": 20,
    }))
    if structure:
        (spec / "structures" / f"{name}.pdb").write_text("ATOM      1  N   MET B   1\nEND\n")
    return spec


def test_target_is_inlined_and_bare_name_dropped(tmp_path):
    """Refusal 1: "unknown target 'X'; this package ships dynorphin_a, hIL2R, ...".

    A bare `target` name resolves only against BC2's shipped presets; a local
    settings/target directory is never searched, so the two stage-5 files never join.
    """
    spec = make_spec(tmp_path)
    stage_campaign(spec, "IL4R_1IAR_B", tmp_path / "job")
    camp = json.loads((tmp_path / "job" / "IL4R_1IAR_B.json").read_text())

    assert "target" not in camp, "the unresolvable bare name must be removed"
    assert len(camp["targets"]) == 1
    assert camp["targets"][0]["name"] == "IL4R_1IAR_B"
    assert camp["targets"][0]["hotspots"] == "127,67,69,126,68,13"
    assert camp["targets"][0]["chains"] == "B"


def test_description_is_not_carried_to_campaign_level(tmp_path):
    """Refusal 2: "unrecognized campaign settings: description (did you mean desperation?)".

    Valid in a target file, rejected at campaign level -- BC2 validates keys strictly.
    """
    spec = make_spec(tmp_path)
    stage_campaign(spec, "IL4R_1IAR_B", tmp_path / "job")
    camp = json.loads((tmp_path / "job" / "IL4R_1IAR_B.json").read_text())
    assert "description" not in camp


def test_target_path_is_absolute_under_the_remote_root(tmp_path):
    """Refusal 3: "target 'X': no file at /work/structures/X.pdb".

    Stage 5 records target_path relative to its own root, but BC2 resolves it from the
    directory of the JSON file it appears in. Absolute removes the ambiguity -- and the
    root must track where the job actually unpacks, which is why it is a parameter.
    """
    spec = make_spec(tmp_path)
    staged = stage_campaign(spec, "IL4R_1IAR_B", tmp_path / "job", remote_root="/scratch/run7")
    camp = json.loads((tmp_path / "job" / "IL4R_1IAR_B.json").read_text())

    tp = camp["targets"][0]["target_path"]
    assert tp == "/scratch/run7/structures/IL4R_1IAR_B.pdb"
    assert staged["project_folder"] == "/scratch/run7/out/il4r_1iar_b"
    # trailing slash on the root must not double up
    staged2 = stage_campaign(spec, "IL4R_1IAR_B", tmp_path / "job2", remote_root="/work/")
    assert "//" not in staged2["project_folder"]


def test_structure_is_copied_next_to_the_campaign(tmp_path):
    spec = make_spec(tmp_path)
    stage_campaign(spec, "IL4R_1IAR_B", tmp_path / "job")
    assert (tmp_path / "job" / "structures" / "IL4R_1IAR_B.pdb").exists()


def test_missing_structure_fails_before_submission(tmp_path):
    """The whole point of staging locally: fail here, not 4 seconds into a GPU job."""
    spec = make_spec(tmp_path, structure=False)
    with pytest.raises(SystemExit) as e:
        stage_campaign(spec, "IL4R_1IAR_B", tmp_path / "job")
    assert "missing" in str(e.value)


def test_missing_spec_files_are_named(tmp_path):
    spec = make_spec(tmp_path)
    with pytest.raises(SystemExit) as e:
        stage_campaign(spec, "NO_SUCH_TARGET", tmp_path / "job")
    assert "campaign" in str(e.value)


def test_resume_defaults_to_false(tmp_path):
    """BC2 continues an existing project_folder by default; an accidental re-submit
    should not become a silent resume."""
    spec = make_spec(tmp_path)
    stage_campaign(spec, "IL4R_1IAR_B", tmp_path / "job")
    camp = json.loads((tmp_path / "job" / "IL4R_1IAR_B.json").read_text())
    assert camp["resume"] is False


def test_campaign_settings_are_preserved(tmp_path):
    """Stage 5's scientific choices must survive staging untouched."""
    spec = make_spec(tmp_path)
    stage_campaign(spec, "IL4R_1IAR_B", tmp_path / "job")
    camp = json.loads((tmp_path / "job" / "IL4R_1IAR_B.json").read_text())
    assert camp["binder_lengths"] == [55, 120]
    assert camp["number_of_final_designs"] == 4
    assert camp["max_trajectories"] == 20
    assert camp["modality"] == "binder"


def test_command_asserts_gpu_before_spending_the_budget(tmp_path):
    """BC2's fan-out log names every GPU even when jax is on the CPU, so the log is not
    evidence. Without this the campaign burns its whole timeout ~100x slower."""
    spec = make_spec(tmp_path)
    staged = stage_campaign(spec, "IL4R_1IAR_B", tmp_path / "job")
    cmd = build_command(staged)
    assert "jax.devices()" in cmd
    assert cmd.index("jax.devices()") < cmd.index("bindcraft design")
    assert "BINDCRAFT_AF2_PARAMS" in cmd


def _write_campaign_output(root: Path, name: str = "il4r_1iar_b") -> Path:
    """A minimal BindCraft2 project_folder, in the layout stage 6 documents."""
    p = root / name
    for sub in ("3_Ranked", "2_Refolded", "1_Trajectories"):
        (p / sub).mkdir(parents=True)
    (p / "campaign_metadata.json").write_text(json.dumps({
        "campaign_name": name, "modality": "binder",
        "targets": [{"name": "IL4R_1IAR_B", "chains": "B",
                     "hotspots": "127,67,69,126,68,13"}],
        "number_of_final_designs": 4,
    }))
    (p / "3_Ranked" / "!_Ranked.csv").write_text(
        "design,rank,binder_sequence,autotuned\nil4r_l55_s1,1,AEELKKLLEEA,\n")
    (p / "2_Refolded" / "!_Refolded.csv").write_text(
        "design,accepted,failed_filters\nil4r_l55_s1,True,\nil4r_l55_s2,False,shape_complementarity\n")
    (p / "1_Trajectories" / "!_Trajectories.csv").write_text(
        "design,terminated\nil4r_l55_s1,\nil4r_l55_s2,clash\n")
    (p / "summary.csv").write_text("metric,value\nn_accepted,1\n")
    return p


def test_project_folder_is_under_out_so_the_harness_keeps_it(tmp_path):
    """A harness constraint that cost a real 28-minute campaign.

    Output is collected from ``./out/`` only. BindCraft2 writes wherever
    ``project_folder`` points, so if that is not under ``out/`` the files are written,
    the job exits, and they are discarded with the per-job working directory -- exit
    code and logs all look normal. Keep the default pointing into ``out/``.
    """
    spec = make_spec(tmp_path)
    staged = stage_campaign(spec, "IL4R_1IAR_B", tmp_path / "job", remote_root="/work")
    assert staged["project_folder"].startswith("/work/out/")


def test_stage6_can_read_a_campaign_folder(tmp_path):
    """The integration claim, tested against ind2b's real reader rather than asserted.

    Skips where stage 6 is not importable -- it arrives in a separate PR.
    """
    stage6 = pytest.importorskip("ind2b.stage6_designs",
                                 reason="stage 6 not present in this checkout")
    folder = _write_campaign_output(tmp_path / "out")

    found = stage6.discover_campaigns([tmp_path / "out"])
    assert [f.name for f in found] == ["il4r_1iar_b"], "campaign discovered by its records"

    rec = stage6.read_campaign(found[0],
                               requested_hotspots=["127", "67", "69", "126", "68", "13"])
    assert rec["campaign_name"] == "il4r_1iar_b"
    assert rec["n_accepted"] == 1
    assert rec["n_attempts"] == 2
    assert rec["n_candidates_scored"] == 2
    assert len(rec["designs"]) == 1


def test_stage6_distinguishes_not_measured_from_measured_zero(tmp_path):
    """A blank is not a zero, and stage 6 enforces that -- this test initially asserted
    the opposite and was wrong.

    A time-bounded campaign is the common case for this runner, so the distinction is
    load-bearing:

      * **no ranked table at all** -> ``n_accepted is None``. The campaign never reached
        the ranking stage. Reporting 0 here would manufacture a measurement, and a reader
        could take it as "nothing passed the filters" when nothing was ever filtered.
      * **an empty ranked table** -> ``n_accepted == 0``. The stage ran and accepted
        nothing, which is a real result.

    Both campaigns are still discovered, so neither disappears from the report.
    """
    stage6 = pytest.importorskip("ind2b.stage6_designs",
                                 reason="stage 6 not present in this checkout")

    never_ranked = tmp_path / "out" / "never_ranked"
    (never_ranked / "1_Trajectories").mkdir(parents=True)
    (never_ranked / "campaign_metadata.json").write_text(
        json.dumps({"campaign_name": "never_ranked"}))
    (never_ranked / "1_Trajectories" / "!_Trajectories.csv").write_text(
        "design,terminated\nd1,clash\n")

    ranked_empty = tmp_path / "out" / "ranked_empty"
    (ranked_empty / "3_Ranked").mkdir(parents=True)
    (ranked_empty / "1_Trajectories").mkdir(parents=True)
    (ranked_empty / "campaign_metadata.json").write_text(
        json.dumps({"campaign_name": "ranked_empty"}))
    (ranked_empty / "3_Ranked" / "!_Ranked.csv").write_text("design,rank,binder_sequence\n")
    (ranked_empty / "1_Trajectories" / "!_Trajectories.csv").write_text(
        "design,terminated\nd1,clash\n")

    found = {f.name for f in stage6.discover_campaigns([tmp_path / "out"])}
    assert found == {"never_ranked", "ranked_empty"}, "both must still be discovered"

    a = stage6.read_campaign(tmp_path / "out" / "never_ranked")
    assert a["n_accepted"] is None, "no ranked table means not measured, not zero"
    assert a["n_attempts"] == 1

    b = stage6.read_campaign(tmp_path / "out" / "ranked_empty")
    assert b["n_accepted"] == 0, "an empty ranked table is a measured zero"


def test_command_points_at_the_staged_parameters(tmp_path):
    """Without the override BC2 lazily downloads 5.3 GB on first campaign, which on an
    account with no job-time egress hangs rather than failing fast."""
    spec = make_spec(tmp_path)
    staged = stage_campaign(spec, "IL4R_1IAR_B", tmp_path / "job")
    assert "/weights/af2/alphafold" in build_command(staged)
