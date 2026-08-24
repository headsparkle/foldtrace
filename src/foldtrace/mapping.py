"""Stage 2 of foldtrace: superpose a candidate on the reference and read active-site state.

The reading is *order-aware*. Rather than taking the residue whose CA happens to sit
nearest each reference catalytic atom (which can grab a spatially-close but
sequence-unrelated residue), foldtrace uses the TM-align residue correspondence, which
preserves sequence order, and only then checks that the corresponding residue's CA lies
within a distance cutoff of the reference position. A call is emitted only when fold
correspondence and spatial placement agree; otherwise the site is reported as
``unresolved`` with an explicit reason rather than dropped.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from tmtools import tm_align

from .io import Structure, Site

# active-site states
RETAINED = "retained"
ALTERED = "altered"
LOST = "lost"
UNRESOLVED = "unresolved"


@dataclass
class SiteCall:
    label: str
    ref_resnum: int
    obs_resnum: int | None
    obs_res: str            # one-letter, or "-" for a gap
    ca_offset: float | None  # CA-CA distance (Angstrom) after superposition
    state: str
    reason: str = ""         # populated for unresolved/altered/lost calls


@dataclass
class CandidateResult:
    name: str
    tm_norm_ref: float
    rmsd: float
    fold_ok: bool
    sites: list[SiteCall] = field(default_factory=list)
    verdict: str = UNRESOLVED
    state_reason: str = ""   # explanation for the verdict (any state), "" when self-evident


def _correspondence(seq_x: str, seq_y: str, rn1: list[int], rn2: list[int]) -> dict[int, tuple[int, str, int]]:
    """From aligned reference/candidate sequences, map ref_resnum -> (cand_resnum, cand_1letter, cand_index)."""
    mapping: dict[int, tuple[int, str, int]] = {}
    i = j = 0
    for a, b in zip(seq_x, seq_y):
        if a != "-" and b != "-":
            mapping[rn1[i]] = (rn2[j], b, j)
        if a != "-":
            i += 1
        if b != "-":
            j += 1
    return mapping


def _apply_transform(coords: np.ndarray, u: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Apply the tmtools rigid-body transform (rotates/translates chain 1 into chain 2's frame)."""
    return coords @ np.asarray(u).T + np.asarray(t)


def map_candidate(ref: Structure, cand: Structure, sites: list[Site],
                  candidate_name: str = "candidate",
                  offset_threshold: float = 4.0,
                  tm_gate: float = 0.5) -> CandidateResult:
    """Superpose ``cand`` on ``ref`` with TM-align and read each site's active-site state.

    Parameters
    ----------
    offset_threshold : max CA-CA distance (Angstrom) between a reference catalytic
        position and its order-aware corresponding candidate residue for a call to be
        trusted; beyond it the site is ``unresolved`` (default 4.0).
    tm_gate : minimum TM-score (normalised by the reference length) for the fold to be
        treated as the same; below it ``fold_ok`` is False and every site is
        ``unresolved`` (default 0.5).

    Each site is classified retained / altered / lost / unresolved: ``expected`` residues
    are retained, ``altered`` residues are altered (changed but related chemistry), a gap
    or an out-of-tolerance offset is unresolved (with a reason), and anything else is lost.
    """
    res = tm_align(ref.ca, cand.ca, ref.seq, cand.seq)
    tm_norm = float(res.tm_norm_chain1)
    fold_ok = tm_norm >= tm_gate

    corr = _correspondence(res.seqxA, res.seqyA, ref.resnums, cand.resnums)
    ref_on_cand = _apply_transform(ref.ca, res.u, res.t)
    ref_on_cand_by_num = {rn: ref_on_cand[i] for i, rn in enumerate(ref.resnums)}

    calls: list[SiteCall] = []
    for site in sites:
        if not fold_ok:
            calls.append(SiteCall(site.label, site.resnum, None, "-", None, UNRESOLVED,
                                  f"fold below TM gate (TM {tm_norm:.2f} < {tm_gate:.2f})"))
            continue
        if site.resnum not in corr:
            calls.append(SiteCall(site.label, site.resnum, None, "-", None, UNRESOLVED,
                                  "reference position aligned to a gap"))
            continue
        cand_resnum, cand_res, cand_idx = corr[site.resnum]
        offset = float(np.linalg.norm(ref_on_cand_by_num[site.resnum] - cand.ca[cand_idx]))
        if offset > offset_threshold:
            state, reason = UNRESOLVED, f"CA offset {offset:.2f} A > {offset_threshold:.2f} A"
        elif cand_res in site.expected:
            state, reason = RETAINED, ""
        elif cand_res in site.altered:
            state, reason = ALTERED, f"{cand_res} in altered set"
        else:
            state, reason = LOST, f"{cand_res} not in expected {{{','.join(sorted(site.expected))}}}"
        calls.append(SiteCall(site.label, site.resnum, cand_resnum, cand_res, offset, state, reason))

    # verdict from the key sites (fall back to all sites if none flagged key)
    key_labels = {s.label for s in sites if s.key}
    driving = [c for c in calls if c.label in key_labels] or calls

    if not fold_ok:
        verdict = UNRESOLVED
        reason = f"fold below TM gate (TM {tm_norm:.2f} < {tm_gate:.2f})"
    elif any(c.state == UNRESOLVED for c in driving):
        verdict = UNRESOLVED
        bad = next(c for c in driving if c.state == UNRESOLVED)
        reason = f"{bad.label}: {bad.reason}"
    elif any(c.state == LOST for c in driving):
        verdict = LOST
        reason = "; ".join(f"{c.label} {c.reason}" for c in driving if c.state == LOST)
    elif any(c.state == ALTERED for c in driving):
        verdict = ALTERED
        reason = "; ".join(f"{c.label} {c.reason}" for c in driving if c.state == ALTERED)
    else:
        verdict = RETAINED
        # a positive, useful reason: the key catalytic residues confirmed in place
        reason = "; ".join(f"{c.label} {c.obs_res}{c.obs_resnum} retained" for c in driving)

    return CandidateResult(candidate_name, tm_norm, float(res.rmsd), fold_ok, calls, verdict, reason)
