"""Deterministic global pairwise alignment, used to project an epitope region from
the target onto a paralogue and read identity over *that region only*.

Needleman-Wunsch with affine gaps (Gotoh), BLOSUM62, BLAST protein defaults
(gap open -11, gap extend -1). Pure Python so the lane carries no alignment
dependency; O(n*m) time and memory, which is acceptable for single-protein
comparisons (a 900x900 pair is ~0.8 M cells).

The alignment is *global*, which is the right choice here: we need a
position-to-position map across the whole chain in order to say which paralogue
residue corresponds to a given target residue. Identity read off a global
alignment of two distantly related proteins is a conservative (low) estimate.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .blosum62 import score as blosum_score

GAP_OPEN = -11
GAP_EXTEND = -1

_M, _IX, _IY = 0, 1, 2  # match / gap-in-b / gap-in-a


@dataclass
class Alignment:
    """A global alignment of `a` (query/target) against `b` (subject/paralogue)."""

    a_aligned: str
    b_aligned: str
    score: float
    # 1-based position of `a` -> 1-based position of `b`; None where `a` aligns to a gap.
    a_to_b: dict[int, int | None] = field(default_factory=dict)
    aligned_columns: int = 0
    identical_columns: int = 0

    @property
    def percent_identity(self) -> float | None:
        """Identity over aligned (non-gap/non-gap) columns, in percent."""
        if self.aligned_columns == 0:
            return None
        return 100.0 * self.identical_columns / self.aligned_columns


def align_global(a: str, b: str, *, gap_open: int = GAP_OPEN, gap_extend: int = GAP_EXTEND) -> Alignment:
    """Gotoh affine-gap global alignment. Deterministic: ties resolved M > Ix > Iy."""
    a = (a or "").upper()
    b = (b or "").upper()
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return Alignment(a_aligned=a or "-" * m, b_aligned=b or "-" * n, score=0.0, a_to_b={}, aligned_columns=0)

    neg = float("-inf")
    # Score rows
    prev_M = [neg] * (m + 1)
    prev_Ix = [neg] * (m + 1)
    prev_Iy = [neg] * (m + 1)
    prev_M[0] = 0.0
    for j in range(1, m + 1):
        prev_Iy[j] = gap_open + (j - 1) * gap_extend
    # Traceback: one byte per cell per state
    tb_M = [bytearray(m + 1) for _ in range(n + 1)]
    tb_Ix = [bytearray(m + 1) for _ in range(n + 1)]
    tb_Iy = [bytearray(m + 1) for _ in range(n + 1)]
    for j in range(1, m + 1):
        tb_Iy[0][j] = _IY if j > 1 else _M

    for i in range(1, n + 1):
        ai = a[i - 1]
        cur_M = [neg] * (m + 1)
        cur_Ix = [neg] * (m + 1)
        cur_Iy = [neg] * (m + 1)
        cur_Ix[0] = gap_open + (i - 1) * gap_extend
        tb_Ix[i][0] = _IX if i > 1 else _M
        rM, rIx, rIy = tb_M[i], tb_Ix[i], tb_Iy[i]
        for j in range(1, m + 1):
            s = blosum_score(ai, b[j - 1])
            # M: diagonal from best of the three
            dM, dIx, dIy = prev_M[j - 1], prev_Ix[j - 1], prev_Iy[j - 1]
            best, st = dM, _M
            if dIx > best:
                best, st = dIx, _IX
            if dIy > best:
                best, st = dIy, _IY
            cur_M[j] = best + s
            rM[j] = st
            # Ix: gap in b (consume a) -- come from row above
            openv, extv = prev_M[j] + gap_open, prev_Ix[j] + gap_extend
            if openv >= extv:
                cur_Ix[j], rIx[j] = openv, _M
            else:
                cur_Ix[j], rIx[j] = extv, _IX
            # Iy: gap in a (consume b) -- come from the left in this row
            openv, extv = cur_M[j - 1] + gap_open, cur_Iy[j - 1] + gap_extend
            if openv >= extv:
                cur_Iy[j], rIy[j] = openv, _M
            else:
                cur_Iy[j], rIy[j] = extv, _IY
        prev_M, prev_Ix, prev_Iy = cur_M, cur_Ix, cur_Iy

    # Choose terminal state
    end_scores = [prev_M[m], prev_Ix[m], prev_Iy[m]]
    state = max(range(3), key=lambda k: (end_scores[k], -k))
    best_score = end_scores[state]

    out_a: list[str] = []
    out_b: list[str] = []
    a_to_b: dict[int, int | None] = {}
    i, j = n, m
    while i > 0 or j > 0:
        if i == 0:
            out_a.append("-"); out_b.append(b[j - 1]); j -= 1; continue
        if j == 0:
            out_a.append(a[i - 1]); out_b.append("-"); a_to_b[i] = None; i -= 1; continue
        if state == _M:
            nxt = tb_M[i][j]
            out_a.append(a[i - 1]); out_b.append(b[j - 1]); a_to_b[i] = j
            i -= 1; j -= 1
        elif state == _IX:
            nxt = tb_Ix[i][j]
            out_a.append(a[i - 1]); out_b.append("-"); a_to_b[i] = None
            i -= 1
        else:
            nxt = tb_Iy[i][j]
            out_a.append("-"); out_b.append(b[j - 1])
            j -= 1
        state = nxt

    aa = "".join(reversed(out_a))
    bb = "".join(reversed(out_b))
    cols = sum(1 for x, y in zip(aa, bb) if x != "-" and y != "-")
    ident = sum(1 for x, y in zip(aa, bb) if x != "-" and y != "-" and x == y)
    return Alignment(
        a_aligned=aa,
        b_aligned=bb,
        score=float(best_score),
        a_to_b=a_to_b,
        aligned_columns=cols,
        identical_columns=ident,
    )


def region_identity(aln: Alignment, a_positions: list[int]) -> dict:
    """Identity of `b` against `a` restricted to the given 1-based positions of `a`.

    Returns counts as well as the percentage, because a percentage over 4 aligned
    positions means something very different from one over 40. Positions of `a`
    that align to a gap in `b` are reported separately and are NOT counted as
    matches or mismatches -- they are missing information.
    """
    per_position = []
    matched = aligned = gapped = 0
    for p in sorted(set(int(x) for x in a_positions)):
        if p < 1 or p > len(aln.a_aligned.replace("-", "")):
            per_position.append({"target_position": p, "subject_position": None, "status": "out_of_range"})
            continue
        q = aln.a_to_b.get(p)
        if q is None:
            gapped += 1
            per_position.append({"target_position": p, "subject_position": None, "status": "gap_in_subject"})
            continue
        # residue letters
        a_res = _residue_at(aln.a_aligned, p)
        b_res = _residue_at(aln.b_aligned, q)
        aligned += 1
        same = a_res == b_res
        matched += 1 if same else 0
        sim = blosum_score(a_res, b_res) > 0
        per_position.append(
            {
                "target_position": p,
                "target_residue": a_res,
                "subject_position": q,
                "subject_residue": b_res,
                "identical": same,
                "similar_blosum62": bool(sim),
                "status": "aligned",
            }
        )
    pct = (100.0 * matched / aligned) if aligned else None
    return {
        "n_requested": len(set(int(x) for x in a_positions)),
        "n_aligned": aligned,
        "n_gapped_in_subject": gapped,
        "n_identical": matched,
        "percent_identity_region": pct,
        "per_position": per_position,
    }


def _residue_at(aligned_seq: str, ungapped_index: int) -> str:
    """Residue at a 1-based ungapped index of an aligned string."""
    k = 0
    for ch in aligned_seq:
        if ch == "-":
            continue
        k += 1
        if k == ungapped_index:
            return ch
    return "X"
