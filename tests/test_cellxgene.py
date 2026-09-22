"""Offline tests for CELLxGENE Discover coverage.

Network is stubbed. What must not drift is the part that decides *what a result
means*: that a term is matched on its ontology id rather than its label, that an
unreachable index is distinguishable from a genuinely absent term, and that a
within-study disease-vs-normal dataset is preferred over a larger disease-only one.
"""

from __future__ import annotations

import pytest

from ind2b.sources import cellxgene as cx


def _ds(did, terms, cells, title="t", tissue="skin of body"):
    return {
        "dataset_id": did,
        "dataset_version_id": did + "-v1",
        "title": title,
        "collection_id": "c1",
        "cell_count": cells,
        "disease": [{"ontology_term_id": t, "label": "irrelevant"} for t in terms],
        "tissue": [{"label": tissue}],
        "assay": [{"label": "10x 3' v3"}],
        "cell_type": [{"label": "keratinocyte"}, {"label": "T cell"}],
    }


INDEX = [
    _ds("ds_paired", ["MONDO:0004980", "PATO:0000461"], 1000),
    _ds("ds_big_disease_only", ["MONDO:0004980"], 9000),
    _ds("ds_other", ["MONDO:0005083"], 500),
]


@pytest.fixture
def stub(monkeypatch):
    monkeypatch.setattr(cx, "_datasets", lambda: INDEX)


def test_matches_on_ontology_id_not_label(stub):
    # The CELLxGENE label for MONDO:0004980 is "atopic eczema", not the input string
    # "atopic dermatitis". Every label in the stub is deliberately "irrelevant".
    assert cx.coverage("MONDO:0004980")["status"] == "covered"


def test_underscore_and_colon_forms_are_the_same_term(stub):
    assert cx.coverage("MONDO_0004980") == cx.coverage("MONDO:0004980")


def test_absent_term_is_a_finding_not_a_failure(stub):
    out = cx.coverage("MONDO:9999999")
    assert out["status"] == "not_in_corpus"
    assert out["n_datasets_scanned"] == len(INDEX)


def test_unreachable_index_is_not_reported_as_absence(monkeypatch):
    monkeypatch.setattr(cx, "_datasets", lambda: None)
    out = cx.coverage("MONDO:0004980")
    assert out["status"] == "source_unavailable"
    assert out["status"] != "not_in_corpus"
    assert "not absent" in out["note"]


def test_every_result_names_its_corpus(stub):
    for term in ("MONDO:0004980", "MONDO:9999999"):
        assert cx.coverage(term)["corpus"] == "discover"


def test_within_study_contrast_beats_a_bigger_disease_only_dataset(stub):
    best = cx.best_contrast_dataset("MONDO:0004980")
    assert best["dataset_id"] == "ds_paired"          # 1,000 cells
    assert best["contrast"] == "within_study_disease_vs_normal"
    # ...even though ds_big_disease_only has nine times as many cells.


def test_disease_only_best_carries_a_between_study_caveat(stub):
    best = cx.best_contrast_dataset("MONDO:0005083")
    assert best["paired_normal_in_same_dataset"] is False
    assert "batch is not separable" in best["caveat"]


def test_no_contrast_dataset_when_not_covered(stub):
    assert cx.best_contrast_dataset("MONDO:9999999") is None