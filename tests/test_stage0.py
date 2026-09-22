"""Offline tests for stage 0 term selection.

These run without network access: the selection rule is pure logic over a
list of ontology hits, which is exactly the part that must not drift.
"""

from __future__ import annotations

import pytest

from ind2b.stage0_resolve import select_disease


HITS = [
    {"id": "EFO_0003060", "name": "non-small cell lung carcinoma"},
    {"id": "MONDO_0005233", "name": "non-small cell lung carcinoma (disease)"},
    {"id": "EFO_0001071", "name": "lung carcinoma"},
]


def test_exact_name_match_wins_over_rank():
    reordered = [HITS[2], HITS[0], HITS[1]]
    hit, rule, ambiguous = select_disease("non-small cell lung carcinoma", reordered)
    assert hit["id"] == "EFO_0003060"
    assert rule == "exact_name_match"
    assert ambiguous is False


def test_case_and_whitespace_insensitive():
    hit, rule, _ = select_disease("  Non-Small Cell Lung Carcinoma ", HITS)
    assert hit["id"] == "EFO_0003060"
    assert rule == "exact_name_match"


def test_token_subset_match_when_no_exact_name():
    hits = [{"id": "EFO_0001071", "name": "lung carcinoma"},
            {"id": "EFO_0003060", "name": "non-small cell lung carcinoma"}]
    hit, rule, ambiguous = select_disease("small cell lung", hits)
    assert hit["id"] == "EFO_0003060"
    assert rule == "all_query_tokens_in_name"
    assert ambiguous is False


def test_falls_back_to_rank_and_flags_ambiguity():
    hits = [{"id": "EFO_0000001", "name": "alpha"}, {"id": "EFO_0000002", "name": "beta"}]
    hit, rule, ambiguous = select_disease("gamma", hits)
    assert hit["id"] == "EFO_0000001"
    assert rule == "top_search_rank"
    assert ambiguous is True


def test_single_rank_hit_is_not_ambiguous():
    hits = [{"id": "EFO_0000001", "name": "alpha"}]
    _, rule, ambiguous = select_disease("gamma", hits)
    assert rule == "top_search_rank"
    assert ambiguous is False


def test_empty_hits_raises():
    with pytest.raises(ValueError):
        select_disease("anything", [])
