"""Offline tests for the structural stages.

These pin the decisions that keep the epitope shortlist honest: which
complexes are admitted, that author numbering is converted rather than
assumed, that hotspots form one compact patch, and that epitope size is
scored as a band.
"""

from __future__ import annotations

import pytest

from ind2b.config import StructureGate
from ind2b.stage3_complexes import classify_entry, gate_entry
from ind2b.stage4_interface import (
    IDEAL_BURIED_HIGH,
    IDEAL_BURIED_LOW,
    auth_to_uniprot,
    designability,
    select_hotspot_patch,
)
from ind2b.stage5_specs import spec_name


def _entity(acc, length, chains, desc):
    return {
        "rcsb_id": "X_1",
        "rcsb_polymer_entity": {"pdbx_description": desc},
        "entity_poly": {"rcsb_sample_sequence_length": length, "type": "polypeptide(L)"},
        "rcsb_polymer_entity_container_identifiers": {
            "auth_asym_ids": chains,
            "asym_ids": chains,
            "reference_sequence_identifiers": (
                [{"database_accession": acc, "database_name": "UniProt"}] if acc else []
            ),
        },
    }


def _entry(entities, method="X-RAY DIFFRACTION", resolution=2.0):
    return {
        "rcsb_id": "TEST",
        "exptl": [{"method": method}],
        "struct": {"title": "test entry"},
        "rcsb_entry_info": {
            "resolution_combined": [resolution] if resolution else None,
            "polymer_entity_count_protein": len(entities),
        },
        "polymer_entities": entities,
    }


TARGET = "P00533"
PARTNERS = {"P01133": {"symbol": "EGF", "accessions": ["P01133"]}}


def test_native_partner_complex_is_admitted():
    entry = _entry([_entity(TARGET, 620, ["A"], "receptor"),
                    _entity("P01133", 53, ["C"], "growth factor")])
    ok, reason = gate_entry(classify_entry(entry, TARGET, PARTNERS), StructureGate())
    assert ok and reason == "admitted"


def test_peptide_fragment_of_target_is_rejected():
    """EGFR appears in the PDB as a 13-residue peptide; that is not an epitope."""
    entry = _entry([_entity(TARGET, 13, ["B"], "EGFR phosphopeptide"),
                    _entity("P01133", 331, ["A"], "growth factor")])
    ok, reason = gate_entry(classify_entry(entry, TARGET, PARTNERS), StructureGate())
    assert not ok
    assert "too_short" in reason


def test_antibody_complex_is_not_a_native_partner():
    entry = _entry([_entity(TARGET, 620, ["A"], "receptor"),
                    _entity(None, 214, ["L"], "Fab DL11 light chain"),
                    _entity(None, 228, ["H"], "Fab DL11 heavy chain")])
    classified = classify_entry(entry, TARGET, PARTNERS)
    roles = {o["role"] for o in classified["other_entities"]}
    assert roles == {"antibody_fragment"}
    ok, reason = gate_entry(classified, StructureGate())
    assert not ok and "no_native_partner" in reason


def test_crystallisation_aid_is_classified_separately():
    entry = _entry([_entity(TARGET, 620, ["A"], "receptor"),
                    _entity("P00698", 129, ["B"], "Lysozyme C")])
    classified = classify_entry(entry, TARGET, PARTNERS)
    assert {o["role"] for o in classified["other_entities"]} == {"crystallisation_aid"}


def test_low_resolution_is_rejected():
    entry = _entry([_entity(TARGET, 620, ["A"], "receptor"),
                    _entity("P01133", 53, ["C"], "growth factor")], resolution=6.0)
    ok, reason = gate_entry(classify_entry(entry, TARGET, PARTNERS), StructureGate())
    assert not ok and "above_max" in reason


def test_missing_resolution_rejected_unless_allowed():
    entry = _entry([_entity(TARGET, 620, ["A"], "receptor"),
                    _entity("P01133", 53, ["C"], "growth factor")],
                   method="SOLUTION NMR", resolution=None)
    classified = classify_entry(entry, TARGET, PARTNERS)
    assert not gate_entry(classified, StructureGate())[0]
    ok, _ = gate_entry(classified, StructureGate(allow_missing_resolution=True))
    assert ok


# --- numbering ------------------------------------------------------------

def test_auth_to_uniprot_offsets_are_applied():
    """4ZQK PD-L1: entity position 1 corresponds to UniProt 18."""
    align = {
        "entities": [
            {
                "entity_id": "4ZQK_1",
                "aligns": [{"accession": "Q9NZQ7",
                             "regions": [{"entity_beg": 1, "ref_beg": 18, "length": 5}]}],
                "chains": {"A": ["18", "19", "20", "21", "22"]},
            }
        ]
    }
    mapping = auth_to_uniprot(align, "Q9NZQ7", "A")
    assert mapping == {"18": 18, "19": 19, "20": 20, "21": 21, "22": 22}


def test_auth_to_uniprot_handles_renumbered_construct():
    """A construct numbered from 1 must still map to real UniProt positions."""
    align = {
        "entities": [
            {
                "entity_id": "X_1",
                "aligns": [{"accession": "P12345",
                             "regions": [{"entity_beg": 1, "ref_beg": 100, "length": 3}]}],
                "chains": {"A": ["1", "2", "3"]},
            }
        ]
    }
    assert auth_to_uniprot(align, "P12345", "A") == {"1": 100, "2": 101, "3": 102}


def test_auth_to_uniprot_returns_empty_for_unknown_chain():
    assert auth_to_uniprot({"entities": []}, "P12345", "A") == {}


# --- hotspots and designability -------------------------------------------

def _res(auth, sasa, coord):
    return {"auth_seq_id": auth, "delta_sasa": sasa, "coord": coord, "is_contact": True}


def test_hotspot_patch_is_compact_not_top_n_by_area():
    """A far-away high-burial residue must not join the patch."""
    candidates = [
        _res("10", 100.0, [0.0, 0.0, 0.0]),
        _res("11", 90.0, [3.0, 0.0, 0.0]),
        _res("12", 80.0, [6.0, 0.0, 0.0]),
        _res("99", 95.0, [60.0, 0.0, 0.0]),   # high burial, opposite face
    ]
    patch = select_hotspot_patch(candidates, 4)
    ids = {c["auth_seq_id"] for c in patch}
    assert "99" not in ids
    assert ids == {"10", "11", "12"}


def test_hotspot_patch_seeds_on_max_burial():
    candidates = [
        _res("1", 20.0, [0.0, 0.0, 0.0]),
        _res("2", 150.0, [30.0, 0.0, 0.0]),
        _res("3", 60.0, [32.0, 0.0, 0.0]),
    ]
    patch = select_hotspot_patch(candidates, 2)
    assert {c["auth_seq_id"] for c in patch} == {"2", "3"}


def test_hotspot_patch_falls_back_without_coordinates():
    candidates = [{"auth_seq_id": "1", "delta_sasa": 50.0, "is_contact": True}]
    assert len(select_hotspot_patch(candidates, 3)) == 1


def test_designability_prefers_the_band_over_larger_interfaces():
    mid = designability((IDEAL_BURIED_LOW + IDEAL_BURIED_HIGH) / 2, 2.0, 6, 6)
    huge = designability(3136.0, 2.0, 6, 6)   # the composite ERBB4 interface
    tiny = designability(120.0, 2.0, 6, 6)
    assert mid["size_score"] == 1.0
    assert huge["designability"] < mid["designability"]
    assert tiny["designability"] < mid["designability"]


def test_designability_rewards_resolution():
    good = designability(800.0, 1.8, 6, 6)
    poor = designability(800.0, 3.9, 6, 6)
    assert good["designability"] > poor["designability"]


def test_spec_name_is_filesystem_safe():
    assert spec_name("CD274", "4ZQK", "A") == "CD274_4ZQK_A"
    assert spec_name("HLA-DRA", "1ICF", "A") == "HLA_DRA_1ICF_A"
