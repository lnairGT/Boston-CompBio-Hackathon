"""Offline tests for stage 2 scoring.

These pin the properties that make the ranking defensible: the accessibility
gate is absolute, the safety penalty does not compound with annotation count,
and a composite is reconstructible from its per-component columns.
"""

from __future__ import annotations

import pytest

from ind2b.config import DEFAULT_PENALTIES, load_weights
from ind2b.stage2_score import components, penalty, score_targets


def _target(**over):
    base = {
        "symbol": "TEST1",
        "ensembl_id": "ENSG00000000001",
        "uniprot": ["P00000"],
        "name": "test target",
        "rank_by_ot_score": 1,
        "ot_overall_score": 0.8,
        "datatype_scores": {
            "genetic_association": 0.6,
            "somatic_mutation": 0.8,
            "literature": 0.9,
            "affected_pathway": 0.5,
            "animal_model": 0.3,
        },
        "surface_accessible": True,
        "accessibility": {
            "mode": "membrane_ectodomain",
            "extracellular_spans": [(20, 240)],
            "reason": None,
        },
        "location_keyword_hint": {"hint": True, "agrees_with_topology": True, "matches": []},
        "tractability": {"AB": ["Approved Drug"]},
        "safety_liabilities": {"count": 0, "events": []},
        "target_class": ["Kinase"],
        "known_drugs": {"count": 10, "has_approved": True, "drug_types": {"Antibody": 2},
                         "has_biologic_precedent": True},
        "chembl": {"n_distinct_mechanisms": 3},
        "trials": {"total_trials": 50, "stopped_trials": [], "example_nct_ids": []},
        "interactions": {"n_partners": 5, "partners": []},
        "expression": {"max_specificity": {"specificity_score": 0.5}},
        "layers_fetched": True,
    }
    base.update(over)
    return base


def test_genetic_component_averages_available_evidence():
    comp = components(_target())
    assert comp["genetic_evidence"] == pytest.approx(0.7)  # mean of 0.6, 0.8


def test_genetic_component_ignores_missing_datatypes():
    comp = components(_target(datatype_scores={"somatic_mutation": 0.9}))
    assert comp["genetic_evidence"] == pytest.approx(0.9)


def test_biologic_precedent_lifts_tractability():
    comp = components(_target(tractability={"AB": ["UniProt loc med conf"]}))
    assert comp["antibody_tractability"] == pytest.approx(0.9)


def test_safety_penalty_is_flat_not_compounding():
    """A target with 21 liabilities must not be punished 21 times over."""
    one, _ = penalty(_target(safety_liabilities={"count": 1}))
    many, _ = penalty(_target(safety_liabilities={"count": 21}))
    assert one == many == pytest.approx(DEFAULT_PENALTIES["safety_liability_any"])


def test_no_penalty_without_liabilities():
    mult, reasons = penalty(_target())
    assert mult == pytest.approx(1.0)
    assert reasons == []


def test_low_specificity_is_penalised():
    mult, reasons = penalty(
        _target(expression={"max_specificity": {"specificity_score": 0.05}})
    )
    assert mult == pytest.approx(DEFAULT_PENALTIES["broad_expression"])
    assert "specificity" in reasons[0]


def test_accessibility_gate_zeroes_the_composite():
    ev = {
        "targets": [
            _target(symbol="REACHABLE"),
            _target(
                symbol="INTRACELLULAR",
                surface_accessible=False,
                accessibility={
                    "mode": None,
                    "extracellular_spans": [],
                    "reason": "membrane_associated_on_cytoplasmic_side",
                },
            ),
        ]
    }
    df, _ = score_targets(ev)
    reach = df[df.symbol == "REACHABLE"].iloc[0]
    intra = df[df.symbol == "INTRACELLULAR"].iloc[0]
    assert reach["composite_score"] > 0
    assert intra["composite_score"] == 0.0
    # the filtered target is still present, with its reason preserved
    assert intra["accessibility_reason"] == "membrane_associated_on_cytoplasmic_side"
    assert reach["rank"] == 1


def test_composite_reconstructs_from_component_columns():
    df, weights = score_targets({"targets": [_target()]})
    row = df.iloc[0]
    rebuilt = sum(row[f"wtd_{k}"] for k in weights) * row["penalty_multiplier"]
    assert rebuilt == pytest.approx(row["composite_score"], abs=1e-4)


def test_weights_file_rejects_unknown_keys(tmp_path):
    bad = tmp_path / "w.json"
    bad.write_text('{"not_a_component": 1.0}')
    with pytest.raises(ValueError):
        load_weights(bad)
