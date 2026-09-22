"""Offline tests for stage 6 binder collection.

These pin the distinctions that keep the binder section honest: what counts as
a made binder, what counts as mini-binder scale, and that an overlap of zero is
reported differently when the deposited construct never contained the epitope.
"""

from __future__ import annotations

import pytest

from ind2b.stage6_binders import (
    MINIBINDER_MAX_RESIDUES,
    _binder_kind,
    compare_to_epitope,
)


def _entity(desc, accessions=(), role="other", length=60):
    return {"description": desc, "accessions": list(accessions),
            "role": role, "length": length, "chains": ["B"]}


def test_designed_binder_recognised_by_description():
    assert _binder_kind(_entity("PD-L1 de novo binder")) == "designed"
    assert _binder_kind(_entity("Designed TrkA-binding miniprotein")) == "designed"
    assert _binder_kind(_entity("Designed Ankyrin Repeat Protein 9_29")) == "designed"


def test_antibody_recognised_by_role():
    assert _binder_kind(_entity("Fab heavy chain", role="antibody_fragment")) == "antibody"


def test_natural_protein_is_not_a_made_binder():
    """A chain with a UniProt accession exists in nature: partner, not binder."""
    assert _binder_kind(_entity("Programmed cell death 1", accessions=["Q15116"])) is None
    # even if its name happens to contain a trigger word
    assert _binder_kind(
        _entity("designed variant of natural ligand", accessions=["P01234"])) is None


def test_unremarkable_chain_is_not_a_binder():
    assert _binder_kind(_entity("Lysozyme")) is None


def test_minibinder_scale_threshold_is_explicit():
    assert MINIBINDER_MAX_RESIDUES == 100


def test_overlap_counts_shared_positions_and_hotspots():
    cmp_ = compare_to_epitope([1, 2, 3, 4], [3, 4, 5], [3, 9])
    assert cmp_["comparable"] is True
    assert cmp_["n_shared"] == 2
    assert cmp_["shared_positions"] == [3, 4]
    assert cmp_["fraction_of_hotspots_covered"] == 0.5
    assert cmp_["hotspots_covered"] == [3]
    assert cmp_["hotspots_missed"] == [9]


def test_jaccard_is_symmetric_intersection_over_union():
    cmp_ = compare_to_epitope([1, 2], [2, 3])
    assert cmp_["jaccard"] == pytest.approx(1 / 3, abs=1e-3)


def test_unmapped_footprint_is_not_comparable():
    """An empty side means no comparison, not a zero overlap."""
    assert compare_to_epitope([], [1, 2], [1])["comparable"] is False
    assert compare_to_epitope([1, 2], [], [])["comparable"] is False


def test_no_hotspots_gives_none_not_zero():
    cmp_ = compare_to_epitope([1, 2], [1, 2], [])
    assert cmp_["fraction_of_hotspots_covered"] is None
