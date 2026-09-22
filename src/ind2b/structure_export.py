"""Trim complexes down to the chains a report should display.

A deposited entry usually carries more than the pair being shown: several
copies of the binder in the asymmetric unit, crystallisation additives,
waters. Embedding all of it in a report is misleading as well as heavy - the
non-contacting copies imply contacts that are not there - so the viewer gets
exactly the target chain and the one binder chain whose interface was
measured.

The output is PDB rather than mmCIF because it is what in-browser viewers
parse most consistently. That is safe only for small complexes: legacy PDB
caps atom serials at 99,999 and chain ids at one character, so this module
refuses to write a file that would silently overflow those fields.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Sequence

log = logging.getLogger("ind2b.structure_export")

# Legacy PDB field limits. Beyond these the format cannot represent the
# structure without renumbering, and a viewer would read the wrong thing.
PDB_MAX_ATOMS = 99_999
PDB_MAX_CHAIN_ID_LEN = 1


class StructureTooLarge(ValueError):
    """The trimmed complex cannot be written as legacy PDB without corruption."""


def trim_complex(
    source: Path,
    keep_chains: Sequence[str],
    dest: Path,
    *,
    drop_solvent: bool = True,
) -> dict:
    """Write ``source`` restricted to ``keep_chains`` as a PDB file.

    Returns a record of what was written, including the atom count, so a
    caller can report the real size rather than an assumed one.
    """
    import gemmi

    source, dest = Path(source), Path(dest)
    st = gemmi.read_structure(str(source))
    st.setup_entities()
    keep = {c for c in keep_chains if c}
    model = st[0]
    for name in [ch.name for ch in model if ch.name not in keep]:
        model.remove_chain(name)
    if drop_solvent:
        st.remove_ligands_and_waters()
    st.remove_alternative_conformations()
    st.remove_empty_chains()
    st.setup_entities()

    model = st[0]
    chains = [ch.name for ch in model]
    missing = sorted(keep - set(chains))
    n_atoms = sum(len(res) for ch in model for res in ch)
    if n_atoms > PDB_MAX_ATOMS:
        raise StructureTooLarge(
            f"{source.name} trimmed to {n_atoms} atoms, above the legacy PDB "
            f"limit of {PDB_MAX_ATOMS}; write mmCIF for this complex instead"
        )
    too_long = [c for c in chains if len(c) > PDB_MAX_CHAIN_ID_LEN]
    if too_long:
        raise StructureTooLarge(
            f"{source.name} has multi-character chain id(s) {too_long}, which "
            f"legacy PDB cannot represent; write mmCIF instead"
        )

    dest.parent.mkdir(parents=True, exist_ok=True)
    st.write_pdb(str(dest))
    return {
        "path": str(dest),
        "chains_written": chains,
        "chains_missing": missing,
        "n_atoms": n_atoms,
        "bytes": dest.stat().st_size,
        "source": str(source),
    }


def export_for_viewer(
    records: Iterable[dict],
    run_dir: Path,
    *,
    subdir: str = "viewer_structures",
) -> int:
    """Trim each binder record's complex and note the result on the record.

    Each record needs ``entry_id``, ``target_chain``, ``footprint`` and
    ``projection``; records without a resolvable structure are left alone
    rather than given a path that does not exist.
    """
    run_dir = Path(run_dir)
    out_dir = run_dir / subdir
    written = 0
    for symbol, rec in records:
        fp = rec.get("footprint") or {}
        src = run_dir / str(fp.get("structure_file") or "")
        binder_chain = fp.get("binder_chain_used")
        if not fp.get("structure_file") or not src.exists() or not binder_chain:
            continue
        dest = out_dir / f'{symbol}_{rec["entry_id"]}_{binder_chain}.pdb'
        try:
            info = trim_complex(src, [rec.get("target_chain"), binder_chain], dest)
        except StructureTooLarge as exc:
            rec["viewer_structure_error"] = str(exc)
            log.warning("%s %s: %s", symbol, rec["entry_id"], exc)
            continue
        except Exception as exc:  # noqa: BLE001
            rec["viewer_structure_error"] = f"{type(exc).__name__}: {exc}"
            log.warning("%s %s: trim failed: %s", symbol, rec["entry_id"], exc)
            continue
        rec["viewer_structure"] = f'{subdir}/{dest.name}'
        rec["viewer_structure_info"] = {
            k: info[k] for k in ("chains_written", "n_atoms", "bytes")
        }
        written += 1
    return written
