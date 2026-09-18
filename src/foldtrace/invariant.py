"""Unexpected invariants: positions nobody declared that refuse to vary.

The declared-site map answers a question you already asked. This module asks the
complementary one -- across everything that shares the fold, what *else* is held
in place, especially in the members where the declared chemistry is gone. The
GFP-fold case is the template: the search was seeded on a chromophore tyrosine
that turned out to be absent, and the family-defining signal was five buried
residues that had never been declared and were invariant in every member.

Two stages, split because the first is expensive and the second is not:

``profile``
    Superpose every candidate on the reference and record the observed residue
    and CA offset at *every* reference position, not only the declared ones.
    One pass, reusable.

``analyse``
    Per position, per state group, compute how conserved it is, and contrast
    that with a background of positions far from the site in the same proteins.
    Rank undeclared near-site positions that hold while the declared chemistry
    does not.

Three guards against the obvious ways this goes wrong:

* A position is only counted where the order-aware correspondence resolves it,
  so a poorly aligned region cannot masquerade as variable.
* Conservation is computed over sequence clusters when a clustering is supplied,
  because a hit set is not a sample of independent proteins and an expanded
  clade will otherwise read as invariance.
* Every finding carries its percentile against the within-set background. If the
  background is also highly conserved -- a small or closely related hit set --
  the percentile says so immediately, and nothing here is significant on its own.

Requires numpy and tmtools, both already runtime dependencies.
"""
from __future__ import annotations

import csv
import math
from collections import Counter
from dataclasses import dataclass, field, asdict

import numpy as np

from .io import Structure, Site
from .mapping import _correspondence, _apply_transform, RETAINED, ALTERED, LOST, UNRESOLVED

PROFILE_COLUMNS = ["candidate", "ref_resnum", "obs_res", "obs_resnum", "offsetA", "resolved"]
FINDING_COLUMNS = ["ref_resnum", "ref_res", "score", "zone", "distA", "modal", "conservation",
                   "coverage", "n", "bg_percentile", "contrast", "entropy_bits", "reason"]


class InvariantError(Exception):
    """Raised for malformed or insufficient input."""


# ------------------------------------------------------------------ profiling
@dataclass
class PositionObs:
    candidate: str
    ref_resnum: int
    obs_res: str
    obs_resnum: int | None
    offset: float | None
    resolved: bool


def profile_candidate(ref: Structure, cand: Structure, candidate_name: str,
                      offset_threshold: float = 4.0, tm_gate: float = 0.5,
                      ) -> tuple[list[PositionObs], float, bool]:
    """Record the observed residue at every reference position, order-aware.

    Returns the per-position observations, the reference-normalised TM-score and
    the fold gate result. A candidate below the gate returns unresolved rows
    rather than nothing, so the number of rows out is independent of the answer.
    """
    from tmtools import tm_align

    res = tm_align(ref.ca, cand.ca, ref.seq, cand.seq)
    tm_norm = float(res.tm_norm_chain1)
    fold_ok = tm_norm >= tm_gate

    if not fold_ok:
        return ([PositionObs(candidate_name, rn, "-", None, None, False)
                 for rn in ref.resnums], tm_norm, False)

    corr = _correspondence(res.seqxA, res.seqyA, ref.resnums, cand.resnums)
    ref_on_cand = _apply_transform(ref.ca, res.u, res.t)
    ref_on_cand_by_num = {rn: ref_on_cand[i] for i, rn in enumerate(ref.resnums)}

    out: list[PositionObs] = []
    for rn in ref.resnums:
        if rn not in corr:
            out.append(PositionObs(candidate_name, rn, "-", None, None, False))
            continue
        cand_resnum, cand_res, cand_idx = corr[rn]
        offset = float(np.linalg.norm(ref_on_cand_by_num[rn] - cand.ca[cand_idx]))
        out.append(PositionObs(candidate_name, rn, cand_res, cand_resnum, offset,
                               offset <= offset_threshold))
    return out, tm_norm, True


def write_profile(obs: list[PositionObs], out, delimiter: str = "\t") -> None:
    close = False
    if hasattr(out, "write"):
        fh = out
    else:
        fh, close = open(out, "w"), True
    try:
        fh.write(delimiter.join(PROFILE_COLUMNS) + "\n")
        for o in obs:
            fh.write(delimiter.join([
                o.candidate, str(o.ref_resnum), o.obs_res,
                str(o.obs_resnum) if o.obs_resnum is not None else "NA",
                f"{o.offset:.2f}" if o.offset is not None else "NA",
                str(o.resolved).lower()]) + "\n")
    finally:
        if close:
            fh.close()


def read_profile(path: str) -> list[PositionObs]:
    obs: list[PositionObs] = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r or (r.get("candidate") or "").startswith("#"):
                continue
            off = (r.get("offsetA") or "NA").strip()
            num = (r.get("obs_resnum") or "NA").strip()
            obs.append(PositionObs(
                r["candidate"].strip(), int(r["ref_resnum"]), (r.get("obs_res") or "-").strip(),
                None if num in ("NA", "") else int(num),
                None if off in ("NA", "") else float(off),
                (r.get("resolved") or "").strip().lower() == "true"))
    if not obs:
        raise InvariantError(f"no profile rows in {path}")
    return obs


# ------------------------------------------------------------------- geometry
def site_centroid(ref: Structure, sites: list[Site]) -> np.ndarray:
    idx = {rn: i for i, rn in enumerate(ref.resnums)}
    pts = [ref.ca[idx[s.resnum]] for s in sites if s.resnum in idx]
    if not pts:
        raise InvariantError("no declared site position is present in the reference structure")
    return np.mean(np.asarray(pts), axis=0)


def distances_to_site(ref: Structure, centroid: np.ndarray) -> dict[int, float]:
    return {rn: float(np.linalg.norm(ref.ca[i] - centroid)) for i, rn in enumerate(ref.resnums)}


def zone_of(dist: float, near: float, background: float) -> str:
    if dist <= near:
        return "site"
    if dist >= background:
        return "background"
    return "mid"


# ----------------------------------------------------------------- statistics
@dataclass
class PositionStats:
    ref_resnum: int
    ref_res: str
    dist: float
    zone: str
    declared: bool
    modal: str = "-"
    conservation: float = 0.0
    coverage: float = 0.0
    n: int = 0
    entropy: float = 0.0
    contrast: float | None = None
    bg_percentile: float | None = None
    score: float = 0.0
    reason: str = ""
    counts: dict = field(default_factory=dict)


def _conservation(residues: list[str]) -> tuple[str, float, float, dict]:
    """Modal residue, its frequency, Shannon entropy in bits, and the counts."""
    counts = Counter(r for r in residues if r and r not in ("-", "X"))
    if not counts:
        return "-", 0.0, 0.0, {}
    total = sum(counts.values())
    modal, top = counts.most_common(1)[0]
    ent = -sum((c / total) * math.log2(c / total) for c in counts.values())
    return modal, top / total, ent, dict(counts)


def _collapse_to_clusters(per_candidate: dict[str, str],
                          clusters: dict[str, str]) -> list[str]:
    """One residue per sequence cluster: the cluster's own modal residue.

    A hit set is not a sample of independent proteins. Without this an expanded
    clade reads as invariance, which is the same non-independence that forced
    cluster reduction and patristic matching in the companion analyses.
    """
    by_cluster: dict[str, list[str]] = {}
    for cand, res in per_candidate.items():
        by_cluster.setdefault(clusters.get(cand, cand), []).append(res)
    out = []
    for members in by_cluster.values():
        counts = Counter(m for m in members if m not in ("-", "X"))
        if counts:
            out.append(counts.most_common(1)[0][0])
    return out


def analyse(profile: list[PositionObs], ref: Structure, sites: list[Site],
            groups: dict[str, str], focus: str = LOST, compare: str | None = RETAINED,
            clusters: dict[str, str] | None = None,
            near: float = 12.0, background: float = 20.0,
            min_coverage: float = 0.60, min_conservation: float = 0.90,
            ) -> tuple[list[PositionStats], dict]:
    """Per-position conservation in the focus group, against a within-set background.

    ``groups`` maps candidate -> state (typically the verdict column of a calls
    table). ``focus`` is the group in which an invariant is unexpected -- by
    default the members that lost the declared chemistry.
    """
    declared = {s.resnum for s in sites}
    ref_res = {rn: ref.seq[i] for i, rn in enumerate(ref.resnums)}
    dists = distances_to_site(ref, site_centroid(ref, sites))

    focus_members = {c for c, g in groups.items() if g == focus}
    if not focus_members:
        raise InvariantError(f"no candidates in focus group {focus!r}")
    compare_members = {c for c, g in groups.items() if g == compare} if compare else set()

    # gather observations per position
    seen_focus: dict[int, dict[str, str]] = {}
    seen_compare: dict[int, dict[str, str]] = {}
    attempted_focus: dict[int, set[str]] = {}
    for o in profile:
        if o.candidate in focus_members:
            attempted_focus.setdefault(o.ref_resnum, set()).add(o.candidate)
            if o.resolved:
                seen_focus.setdefault(o.ref_resnum, {})[o.candidate] = o.obs_res
        elif o.candidate in compare_members and o.resolved:
            seen_compare.setdefault(o.ref_resnum, {})[o.candidate] = o.obs_res

    stats: list[PositionStats] = []
    for rn in ref.resnums:
        obs = seen_focus.get(rn, {})
        attempted = len(attempted_focus.get(rn, set())) or len(focus_members)
        residues = (_collapse_to_clusters(obs, clusters) if clusters
                    else [r for r in obs.values() if r not in ("-", "X")])
        modal, cons, ent, counts = _conservation(residues)
        st = PositionStats(
            ref_resnum=rn, ref_res=ref_res.get(rn, "X"), dist=dists[rn],
            zone=zone_of(dists[rn], near, background), declared=rn in declared,
            modal=modal, conservation=cons, coverage=len(obs) / attempted if attempted else 0.0,
            n=len(residues), entropy=ent, counts=counts)
        if seen_compare.get(rn):
            comp = (_collapse_to_clusters(seen_compare[rn], clusters) if clusters
                    else [r for r in seen_compare[rn].values() if r not in ("-", "X")])
            _, ccons, _, _ = _conservation(comp)
            st.contrast = cons - ccons
        stats.append(st)

    # within-set background: positions far from the site with adequate coverage
    bg = sorted(s.conservation for s in stats
                if s.zone == "background" and s.coverage >= min_coverage)
    summary = {
        "focus_group": focus, "focus_n": len(focus_members),
        "compare_group": compare or "-", "compare_n": len(compare_members),
        "clusters": len({clusters[c] for c in clusters if c in focus_members}) if clusters else 0,
        "background_positions": len(bg),
        "background_median_conservation": f"{bg[len(bg) // 2]:.3f}" if bg else "NA",
        "background_p90_conservation": f"{bg[int(0.9 * (len(bg) - 1))]:.3f}" if bg else "NA",
    }

    for s in stats:
        s.bg_percentile = (sum(1 for b in bg if b < s.conservation) / len(bg)) if bg else None
        s.score, s.reason = _score(s, min_coverage, min_conservation, near)

    # how much of the site zone came back as invariant. A high fraction usually
    # means the hit set is too close or too small to discriminate, not that the
    # whole site is a conserved device.
    site_tested = [s for s in stats if s.zone == "site" and not s.declared
                   and s.coverage >= min_coverage]
    reported = [s for s in site_tested if s.score > 0]
    summary["site_positions_tested"] = len(site_tested)
    summary["site_positions_reported"] = len(reported)
    summary["site_fraction_reported"] = (
        f"{len(reported) / len(site_tested):.2f}" if site_tested else "NA")
    return stats, summary


def _score(s: PositionStats, min_coverage: float, min_conservation: float,
           near: float) -> tuple[float, str]:
    if s.declared:
        return 0.0, "declared position; this is the chemistry you asked about"
    if s.zone != "site":
        return 0.0, f"{s.zone} zone ({s.dist:.1f} A from the site centroid)"
    if s.coverage < min_coverage:
        return 0.0, (f"coverage {s.coverage:.2f} < {min_coverage:.2f}; too few focus-group "
                     "members resolve this position to judge it")
    if s.n < 3:
        return 0.0, f"only {s.n} independent observations"
    if s.conservation < min_conservation:
        return 0.0, f"conservation {s.conservation:.2f} < {min_conservation:.2f}"

    prox = max(0.0, min(1.0, 1.0 - s.dist / near))
    pct = s.bg_percentile if s.bg_percentile is not None else 0.5
    score = 0.45 * s.conservation + 0.35 * pct + 0.20 * prox
    reason = (f"{s.modal} in {s.conservation:.0%} of {s.n} independent observations, "
              f"{s.dist:.1f} A from the site centroid, above {pct:.0%} of background positions")
    if s.contrast is not None and s.contrast > 0.05:
        reason += f"; more conserved here than in the comparison group (+{s.contrast:.2f})"
    return score, reason


def rank(stats: list[PositionStats], min_score: float = 0.0) -> list[PositionStats]:
    keep = [s for s in stats if s.score > min_score]
    return sorted(keep, key=lambda s: (-s.score, s.ref_resnum))


def write_findings(stats: list[PositionStats], out, delimiter: str = "\t") -> None:
    close = False
    if hasattr(out, "write"):
        fh = out
    else:
        fh, close = open(out, "w"), True
    try:
        fh.write(delimiter.join(FINDING_COLUMNS) + "\n")
        for s in stats:
            fh.write(delimiter.join([
                str(s.ref_resnum), s.ref_res, f"{s.score:.3f}", s.zone, f"{s.dist:.1f}",
                s.modal, f"{s.conservation:.3f}", f"{s.coverage:.3f}", str(s.n),
                f"{s.bg_percentile:.3f}" if s.bg_percentile is not None else "NA",
                f"{s.contrast:+.3f}" if s.contrast is not None else "NA",
                f"{s.entropy:.2f}", s.reason.replace(delimiter, " ")]) + "\n")
    finally:
        if close:
            fh.close()


def emit_features(stats: list[PositionStats], top: int = 1) -> list[str]:
    """Suggest a `prime adjudicate` call for the recovered position set."""
    hits = rank(stats)[:12]
    if not hits:
        return []
    listing = ", ".join(f"{s.modal}{s.ref_resnum}" for s in hits[:6])
    noun = "position" if len(hits) == 1 else "positions"
    return [
        f"foldtrace prime adjudicate --register register.json \\\n"
        f"  --feature {f'{len(hits)} undeclared near-site {noun} invariant in the focus group ({listing})'!r} \\\n"
        f"  --surfaced-by unexpected_invariant --status TODO --permitted TODO"]


# ------------------------------------------------------------------- loading
def load_groups(path: str, key: str = "candidate_id", value: str = "verdict") -> dict[str, str]:
    """Group labels from a calls table (candidate -> verdict) or any two-column TSV."""
    out: dict[str, str] = {}
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise InvariantError(f"no rows in {path}")
    kcol = key if key in rows[0] else next(
        (c for c in ("candidate_id", "candidate", "target") if c in rows[0]), None)
    vcol = value if value in rows[0] else None
    if kcol is None or vcol is None:
        raise InvariantError(f"{path}: need a candidate column and a {value!r} column")
    for r in rows:
        k = (r.get(kcol) or "").strip()
        if k:
            out[k] = (r.get(vcol) or UNRESOLVED).strip().lower()
    return out


def load_clusters(path: str) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = (r.get("candidate") or r.get("candidate_id") or r.get("target") or "").strip()
            c = (r.get("cluster") or "").strip()
            if k and c:
                out[k] = c
    if not out:
        raise InvariantError(f"no candidate/cluster pairs in {path}")
    return out
