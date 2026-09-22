"""Offline tests for trimming complexes down to what a viewer should show.

The guards matter more than the happy path: legacy PDB silently cannot
represent large atom counts or multi-character chain ids, and a viewer reading
a truncated file would show the wrong molecule without complaining.
"""

from __future__ import annotations

import pytest

gemmi = pytest.importorskip("gemmi")

from ind2b.structure_export import (
    PDB_MAX_ATOMS,
    StructureTooLarge,
    trim_complex,
)


def _structure(chain_specs, with_water=False):
    """Build a tiny structure: {chain_name: n_residues}."""
    st = gemmi.Structure()
    st.spacegroup_hm = "P 1"
    st.cell = gemmi.UnitCell(50, 50, 50, 90, 90, 90)
    model = gemmi.Model("1")
    for name, n_res in chain_specs.items():
        chain = gemmi.Chain(name)
        for i in range(1, n_res + 1):
            res = gemmi.Residue()
            res.name = "ALA"
            res.seqid = gemmi.SeqId(i, " ")
            for atom_name, dx in (("N", 0.0), ("CA", 1.5), ("C", 3.0)):
                atom = gemmi.Atom()
                atom.name = atom_name
                atom.element = gemmi.Element("C")
                atom.pos = gemmi.Position(i * 4.0 + dx, 0, 0)
                res.add_atom(atom)
            chain.add_residue(res)
        model.add_chain(chain)
    if with_water:
        w = gemmi.Chain("W")
        res = gemmi.Residue()
        res.name = "HOH"
        res.seqid = gemmi.SeqId(1, " ")
        atom = gemmi.Atom()
        atom.name = "O"
        atom.element = gemmi.Element("O")
        atom.pos = gemmi.Position(9, 9, 9)
        res.add_atom(atom)
        w.add_residue(res)
        model.add_chain(w)
    st.add_model(model)
    st.setup_entities()
    return st


def _write_cif(st, path):
    st.make_mmcif_document().write_file(str(path))
    return path


def test_keeps_only_requested_chains(tmp_path):
    src = _write_cif(_structure({"A": 4, "B": 4, "C": 4}), tmp_path / "in.cif")
    info = trim_complex(src, ["A", "C"], tmp_path / "out.pdb")
    assert sorted(info["chains_written"]) == ["A", "C"]
    assert info["n_atoms"] == 24          # 8 residues x 3 atoms
    text = (tmp_path / "out.pdb").read_text()
    assert " B " not in text


def test_requested_chain_that_is_absent_is_reported_not_silently_dropped(tmp_path):
    src = _write_cif(_structure({"A": 3}), tmp_path / "in.cif")
    info = trim_complex(src, ["A", "Z"], tmp_path / "out.pdb")
    assert info["chains_missing"] == ["Z"]


def test_solvent_is_dropped_by_default(tmp_path):
    src = _write_cif(_structure({"A": 3}, with_water=True), tmp_path / "in.cif")
    info = trim_complex(src, ["A", "W"], tmp_path / "out.pdb")
    assert "HOH" not in (tmp_path / "out.pdb").read_text()
    assert info["n_atoms"] == 9


def test_multi_character_chain_id_is_refused(tmp_path):
    """Legacy PDB has one column for the chain id; two characters cannot fit."""
    src = _write_cif(_structure({"AA": 3}), tmp_path / "in.cif")
    with pytest.raises(StructureTooLarge, match="multi-character chain id"):
        trim_complex(src, ["AA"], tmp_path / "out.pdb")


def test_atom_limit_is_enforced_before_writing(tmp_path, monkeypatch):
    """The guard fires rather than emitting a file with overflowed serials."""
    monkeypatch.setattr("ind2b.structure_export.PDB_MAX_ATOMS", 10)
    src = _write_cif(_structure({"A": 8}), tmp_path / "in.cif")
    dest = tmp_path / "out.pdb"
    with pytest.raises(StructureTooLarge, match="above the legacy PDB limit"):
        trim_complex(src, ["A"], dest)
    assert not dest.exists()


def test_documented_pdb_atom_limit():
    assert PDB_MAX_ATOMS == 99_999
