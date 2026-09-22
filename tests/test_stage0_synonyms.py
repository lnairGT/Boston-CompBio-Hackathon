"""Offline tests for exact-synonym term selection.

The bug these pin down: a disease's everyday name is often a MONDO *synonym*
rather than its preferred label, and matching query tokens against labels alone
can then select a different disease with high confidence. The real case is
"atopic dermatitis", whose intended term is labelled "atopic eczema" - sharing
no token with the query - while a rare Mendelian term contains the query
verbatim.

These run without network access: the synonym lookup is injected.
"""

from __future__ import annotations

from ind2b.stage0_resolve import select_disease


# Ranked exactly as Open Targets returns them for "atopic dermatitis": the
# intended disease is the top hit, and the rare Mendelian term is third.
AD_HITS = [
    {"id": "MONDO_0004980", "name": "atopic eczema"},
    {"id": "MONDO_0002406", "name": "dermatitis"},
    {"id": "MONDO_0054697", "name": "immunodeficiency 11b with atopic dermatitis"},
    {"id": "MONDO_0005202", "name": "atopic IgE-mediated allergic disorder"},
]

# Abridged from the live API. "atopic dermatitis" is an exact synonym of the
# intended term and of nothing else here.
AD_SYNONYMS = {
    "MONDO_0004980": {
        "hasExactSynonym": [
            "ATOD", "allergic dermatitis", "atopic dermatitis", "atopic eczema",
            "dermatitis, atopic", "eczema", "eczematous dermatitis",
        ],
        "hasRelatedSynonym": [],
    },
    "MONDO_0002406": {"hasExactSynonym": ["dermatitis"]},
    "MONDO_0054697": {"hasExactSynonym": ["immunodeficiency 11b with atopic dermatitis"]},
    "MONDO_0005202": {"hasExactSynonym": []},
}


def _syn(ids):
    return {i: AD_SYNONYMS.get(i, {}) for i in ids}


def test_synonym_match_beats_token_match_on_another_term():
    """The regression: without synonyms this returns the rare immunodeficiency."""
    hit, rule, ambiguous = select_disease("atopic dermatitis", AD_HITS, synonyms_of=_syn)
    assert hit["id"] == "MONDO_0004980"
    assert rule == "exact_synonym_match"
    assert ambiguous is False


def test_without_synonyms_the_old_precedence_is_unchanged():
    """Omitting the lookup must not change existing behaviour."""
    hit, rule, _ = select_disease("atopic dermatitis", AD_HITS)
    assert hit["id"] == "MONDO_0054697"
    assert rule == "all_query_tokens_in_name"


def test_exact_name_still_wins_over_synonym():
    """A label match is more specific evidence than a synonym match."""
    hit, rule, _ = select_disease("dermatitis", AD_HITS, synonyms_of=_syn)
    assert hit["id"] == "MONDO_0002406"
    assert rule == "exact_name_match"


def test_broad_and_narrow_synonyms_are_not_identity():
    """Broad/narrow synonyms name a different disease and must not match."""
    hits = [{"id": "X_1", "name": "some disease"}]
    syn = {"X_1": {"hasBroadSynonym": ["skin disease"],
                   "hasNarrowSynonym": ["severe skin disease"],
                   "hasExactSynonym": []}}
    hit, rule, _ = select_disease("skin disease", hits, synonyms_of=lambda ids: syn)
    assert rule == "top_search_rank"  # fell through; the broad synonym was ignored


def test_synonym_lookup_failure_falls_back_rather_than_raising():
    """A synonym-service outage must not take term resolution down with it."""
    def boom(ids):
        raise RuntimeError("OLS is down")

    hit, rule, _ = select_disease("atopic dermatitis", AD_HITS, synonyms_of=boom)
    assert rule == "all_query_tokens_in_name"  # degraded, not crashed


def test_no_hits_still_raises():
    import pytest

    with pytest.raises(ValueError):
        select_disease("atopic dermatitis", [], synonyms_of=_syn)
