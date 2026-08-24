"""Core of foldtrace: superpose a candidate on the reference and read the active-site state.

The reading is *order-aware*. Rather than taking the residue whose CA happens to sit
nearest each reference catalytic atom (which can grab a spatially-close but
sequence-unrelated residue), foldtrace uses the TM-align residue correspondence, which
preserves sequence order, and only then checks that the corresponding residue's CA lies
within a distance cutoff of the reference position. A call is emitted only when fold
correspondence and spatial placement agree.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from tmtools import tm_align

from .io import Structure, Site

# active-site states
RETAINED = "retained"
LOST = "lost"
UNRESOLVED = "unresolved"


@dataclass
class SiteCall:
    label: str
    ref_resnum: int
    obs_resnum: int | None
    obs_res: str          # one-letter, or "-" for a gap
    ca_offset: float | None
    state: str


@dataclass
class CandidateResult:
    name: str
    tm_norm_ref: float
    rmsd: float
    fold_ok: bool
    sites: list[SiteCall] = field(default_factory=list)
    verdict: str = UNRESOLVED


def _correspondence(seq_x: str, seq_y: str, rn1: list[int], rn2: list[int],
                    names2: list[str]) -> dict[int, tuple[int, str, int]]:
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
        trusted; beyond it the site is ``unresolved``.
    tm_gate : minimum TM-score (normalised by the reference length) for the fold to be
        considered the same; below it ``fold_ok`` is False and the verdict is ``unresolved``.
    """
    res = tm_align(ref.ca, cand.ca, ref.seq, cand.seq)
    tm_norm = float(res.tm_norm_chain1)
    fold_ok = tm_norm >= tm_gate

    corr = _correspondence(res.seqxA, res.seqyA, ref.resnums, cand.resnums, cand.resnames)
    ref_ca_by_num = {rn: ref.ca[i] for i, rn in enumerate(ref.resnums)}
    ref_on_cand = _apply_transform(ref.ca, res.u, res.t)
    ref_on_cand_by_num = {rn: ref_on_cand[i] for i, rn in enumerate(ref.resnums)}

    calls: list[SiteCall] = []
    for site in sites:
        if site.resnum not in corr:
            calls.append(SiteCall(site.label, site.resnum, None, "-", None, UNRESOLVED))
            continue
        cand_resnum, cand_res, cand_idx = corr[site.resnum]
        offset = float(np.linalg.norm(ref_on_cand_by_num[site.resnum] - cand.ca[cand_idx]))
        if not fold_ok or offset > offset_threshold:
            state = UNRESOLVED
        elif cand_res in site.expected:
            state = RETAINED
        else:
            state = LOST
        calls.append(SiteCall(site.label, site.resnum, cand_resnum, cand_res, offset, state))

    # verdict from the key sites (fall back to all sites if none flagged key)
    key_calls = [c for c in calls if any(s.key and s.label == c.label for s in sites)]
    driving = key_calls or calls
    if not fold_ok:
        verdict = UNRESOLVED
    elif any(c.state == UNRESOLVED for c in driving):
        verdict = UNRESOLVED
    elif all(c.state == RETAINED for c in driving):
        verdict = RETAINED
    else:
        verdict = LOST

    return CandidateResult(candidate_name, tm_norm, float(res.rmsd), fold_ok, calls, verdict)
