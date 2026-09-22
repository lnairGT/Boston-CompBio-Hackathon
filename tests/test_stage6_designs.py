"""Offline tests for the BindCraft2 design reader.

No campaign can run here, so these build a fixture directory in the layout the
upstream output documentation specifies and check the reader against it. The
fixture exists only in tests: it is never written into a run directory and
never reaches a report.

What these pin is the set of documented ways these tables can be misread:
a blank is not a zero, "not measured" is not a low score, rank is not a
probability, and an attempt that ran on the desperation ladder was given an
easier problem than the campaign asked for.
"""

from __future__ import annotations

import json

import pytest

from ind2b.stage6_designs import (
    ACCEPTANCE_MEANING,
    LADDER_SETTINGS,
    discover_campaigns,
    epitope_engagement,
    parse_interface_residues,
    read_campaign,
)

RANKED_HEADER = (
    "design,hash,rank,trajectory,outcome,failed_filters,autotuned,"
    "Binder_Sequence,Interface_Binder_Residues,Interface_Target_Residues,"
    "i_pDAE,iptm,Binder_RMSD\n"
)


def _campaign(tmp_path, ranked_rows=(), refolded=None, trajectories=None, meta=True):
    folder = tmp_path / "results" / "CD274_campaign"
    (folder / "3_Ranked").mkdir(parents=True)
    if meta:
        (folder / "campaign_metadata.json").write_text(json.dumps({
            "campaign_name": "CD274_4ZQK_A",
            "target": "CD274_4ZQK_A",
            "modality": "binder",
            "binder_lengths": [55, 120],
            "source_revision": "abc1234",
            "resolved": {"validation_model": "model_1"},
        }))
    if ranked_rows:
        (folder / "3_Ranked" / "!_Ranked.csv").write_text(
            RANKED_HEADER + "".join(ranked_rows))
    if refolded is not None:
        (folder / "2_Refolded").mkdir(parents=True, exist_ok=True)
        (folder / "2_Refolded" / "!_Refolded.csv").write_text(refolded)
    if trajectories is not None:
        (folder / "1_Trajectories").mkdir(parents=True, exist_ok=True)
        (folder / "1_Trajectories" / "!_Trajectories.csv").write_text(trajectories)
    return folder


def test_discovers_campaign_by_its_records_not_its_name(tmp_path):
    folder = _campaign(tmp_path)
    found = discover_campaigns([tmp_path / "results"])
    assert folder in found


def test_campaign_with_no_output_is_not_discovered(tmp_path):
    (tmp_path / "results" / "empty").mkdir(parents=True)
    assert discover_campaigns([tmp_path / "results"]) == []


def test_reads_accepted_design_and_its_sequence(tmp_path):
    row = ("d1,h1,1,7,passed,,,QLEEL/GG,3,56,0.21,0.83,1.4\n")
    folder = _campaign(tmp_path, [row])
    rec = read_campaign(folder, requested_hotspots=["56", "58"])
    assert rec["n_accepted"] == 1
    d = rec["designs"][0]
    assert d["binder_length"] == 7      # the chain separator is not a residue
    assert d["binder_chains"] == 2
    assert rec["acceptance_meaning"] == ACCEPTANCE_MEANING


def test_requested_hotspot_engagement_is_measured_not_assumed(tmp_path):
    row = "d1,h1,1,7,passed,,,QLEEL,3,\"56,60,61\",0.21,0.83,1.4\n"
    folder = _campaign(tmp_path, [row])
    rec = read_campaign(folder, requested_hotspots=["56", "58"])
    eng = rec["designs"][0]["vs_requested_epitope"]
    assert eng["comparable"] is True
    assert eng["hotspots_engaged"] == ["56"]
    assert eng["hotspots_missed"] == ["58"]
    assert eng["fraction_hotspots_engaged"] == 0.5


def test_missing_interface_residues_is_not_zero_contacts(tmp_path):
    row = "d1,h1,1,7,passed,,,QLEEL,,,0.21,0.83,1.4\n"
    folder = _campaign(tmp_path, [row])
    rec = read_campaign(folder, requested_hotspots=["56"])
    eng = rec["designs"][0]["vs_requested_epitope"]
    assert eng["comparable"] is False
    assert "not recorded" in eng["reason"]


def test_blank_metric_stays_none_rather_than_zero(tmp_path):
    row = "d1,h1,1,7,passed,,,QLEEL,3,56,,0.83,\n"
    folder = _campaign(tmp_path, [row])
    metrics = read_campaign(folder)["designs"][0]["metrics"]
    assert metrics["i_pDAE"] is None
    assert metrics["Binder_RMSD"] is None
    assert metrics["iptm"] == pytest.approx(0.83)


def test_desperation_ladder_attempt_is_flagged(tmp_path):
    row = ("d1,h1,1,7,passed,,initial_guess=True; design_recycles=6,"
           "QLEEL,3,56,0.21,0.83,1.4\n")
    d = read_campaign(_campaign(tmp_path, [row]))["designs"][0]
    assert d["ran_on_desperation_ladder"] is True
    assert set(d["ladder_settings"]) == {"initial_guess", "design_recycles"}


def test_ordinary_autotuning_is_not_a_ladder_flag(tmp_path):
    row = "d1,h1,1,7,passed,,soft_iterations=60,QLEEL,3,56,0.21,0.83,1.4\n"
    d = read_campaign(_campaign(tmp_path, [row]))["designs"][0]
    assert d["ran_on_desperation_ladder"] is False


def test_ladder_settings_match_documented_names():
    assert set(LADDER_SETTINGS) == {
        "initial_guess", "target_flexibility", "validation_model", "design_recycles"}


def test_campaign_that_accepted_nothing_reports_why(tmp_path):
    folder = _campaign(
        tmp_path,
        refolded="design,outcome,failed_filters\nc1,rejected,\"i_pDAE,Binder_RMSD\"\n"
                 "c2,rejected,i_pDAE\n",
        trajectories="design,hash,terminated\nt1,h1,soft\n")
    rec = read_campaign(folder)
    assert rec["n_accepted"] is None
    assert rec["stages_present"]["ranked"] is False
    assert rec["failed_filter_tally"]["i_pDAE"] == 2
    assert rec["termination_tally"]["soft"] == 1
    assert "not a statement about the target" in rec["designs_note"]


def test_interface_residue_parsing_splits_chains_and_commas():
    assert parse_interface_residues("12,14/20,21") == ["12", "14", "20", "21"]
    assert parse_interface_residues(None) == []
    assert parse_interface_residues("") == []


def test_engagement_without_requested_hotspots_is_not_comparable():
    out = epitope_engagement("56,57", [])
    assert out["comparable"] is False
