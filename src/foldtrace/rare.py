"""Rare residues: chemistry carried by a small minority on a common scaffold.

The mirror image of the invariant detector. There, a position refuses to vary
and the question is what is being held. Here a position varies freely and the
question is whether a handful of members carry something at it that the rest do
not -- the His287 case, present in 15 of 1000 beta-propellers and the difference
between a lactonase and an enzyme that breaks a P-F bond.

This is the distribution annotation handles worst. A rare activity on a common
scaffold is almost never in the protein's name, so neither a text search nor a
profile search will separate the carriers from the family around them, while
reading the position directly does it in one pass.

Two outputs. Per position, the rare residue and who carries it. Then, across
positions, the members carrying several rare residues at once -- a complete site
rather than a single lucky substitution, which is the difference between 15
carriers and the 9 with the full device.

Consumes the table from ``foldtrace profile``. Stdlib plus numpy.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field

from .invariant import (PositionObs, InvariantError, site_centroid, distances_to_site,
                        zone_of, _collapse_to_clusters)
from .io import Structure, Site

# residues able to carry acid/base, nucleophile, metal-ligand or cation chemistry
CHEMISTRY_RESIDUES = frozenset("HCDEKRSTY W".replace(" ", ""))

RARE_COLUMNS = ["ref_resnum", "residue", "n_carriers", "fraction", "distA", "zone",
                "score", "modal", "modal_frac", "carriers", "reason"]
COMBO_COLUMNS = ["n_positions", "positions", "n_members", "members", "score", "reason"]


@dataclass
class RareFinding:
    ref_resnum: int
    residue: str
    n_carriers: int
    fraction: float
    dist: float
    zone: str
    modal: str
    modal_frac: float = 0.0
    carriers: list[str] = field(default_factory=list)
    score: float = 0.0
    reason: str = ""


@dataclass
class ComboFinding:
    positions: list[int]
    members: list[str]
    score: float = 0.0
    reason: str = ""


def analyse_rare(profile: list[PositionObs], ref: Structure, sites: list[Site],
                 near: float = 12.0, background: float = 20.0,
                 residues: frozenset[str] = CHEMISTRY_RESIDUES,
                 min_count: int = 3, max_fraction: float = 0.10,
                 min_modal: float = 0.50,
                 clusters: dict[str, str] | None = None,
                 include_declared: bool = False) -> tuple[list[RareFinding], dict]:
    """Rare chemistry-capable residues at positions near the declared site.

    ``max_fraction`` is what makes a residue rare: present in at most this share
    of the members that resolve the position. ``min_count`` stops a single model
    from becoming a finding. When ``clusters`` is supplied the carrier count is
    the number of distinct sequence clusters carrying the residue, so an
    expanded clade does not manufacture a recurrent gain.

    ``min_modal`` is the one that matters most in practice. A rare residue is
    only informative where the position is otherwise consistent -- a catalytic
    histidine standing out against an invariant glycine. At a position whose
    composition is diffuse, every residue is rare and none of them means
    anything, so positions below this modal frequency are not tested at all.
    """
    declared = {s.resnum for s in sites}
    dists = distances_to_site(ref, site_centroid(ref, sites))
    ref_res = {rn: ref.seq[i] for i, rn in enumerate(ref.resnums)}

    seen: dict[int, dict[str, str]] = defaultdict(dict)
    for o in profile:
        if o.resolved and o.obs_res not in ("-", "X"):
            seen[o.ref_resnum][o.candidate] = o.obs_res

    findings: list[RareFinding] = []
    for rn, per_cand in seen.items():
        if rn in declared and not include_declared:
            continue
        zone = zone_of(dists[rn], near, background)
        if zone != "site":
            continue
        # Cluster handling differs from the invariant detector on purpose. There,
        # a cluster votes with its modal residue. Here that would erase the
        # finding: a rare residue is a minority by definition and would be
        # outvoted inside its own cluster. So a cluster *carries* a residue if
        # any of its members does, while the denominator stays the cluster
        # count -- one clade still gets one vote, but a minority is not silenced.
        pool = (_collapse_to_clusters(per_cand, clusters) if clusters
                else list(per_cand.values()))
        if not pool:
            continue
        units = ({clusters.get(c, c) for c in per_cand} if clusters else set(per_cand))
        total = len(units)
        if not total:
            continue
        modal_counts = Counter(pool)
        modal, modal_n = modal_counts.most_common(1)[0]
        modal_frac = modal_n / len(pool)
        if modal_frac < min_modal:
            continue  # diffuse position: rarity here carries no information
        by_residue: dict[str, set[str]] = {}
        for cand, r in per_cand.items():
            by_residue.setdefault(r, set()).add(clusters.get(cand, cand) if clusters else cand)
        for aa, unit_set in by_residue.items():
            n = len(unit_set)
            if aa not in residues or aa == modal:
                continue
            frac = n / total
            if n < min_count or frac > max_fraction:
                continue
            carriers = sorted(c for c, r in per_cand.items() if r == aa)
            prox = max(0.0, min(1.0, 1.0 - dists[rn] / near))
            rarity = 1.0 - frac / max_fraction
            score = (0.35 * rarity + 0.30 * modal_frac + 0.20 * prox
                     + 0.15 * min(1.0, n / 10))
            findings.append(RareFinding(
                rn, aa, n, frac, dists[rn], zone, modal, modal_frac, carriers, score,
                f"{aa} at {n} of {total} ({frac:.1%}) where the position is otherwise "
                f"{modal} in {modal_frac:.0%} of members; {dists[rn]:.1f} A from the "
                "site centroid"))

    summary = {
        "positions_tested": sum(1 for rn in seen if zone_of(dists[rn], near, background) == "site"
                                and (include_declared or rn not in declared)),
        "rare_calls": len(findings),
        "max_fraction": f"{max_fraction:.2f}",
        "min_count": min_count,
        "min_modal": f"{min_modal:.2f}",
        "clustered": bool(clusters),
    }
    return sorted(findings, key=lambda f: (-f.score, f.ref_resnum)), summary


def combinations(findings: list[RareFinding], min_positions: int = 2) -> list[ComboFinding]:
    """Members carrying a rare residue at several positions at once.

    One rare substitution is a substitution. Several in the same protein, at
    positions that share a site, is a candidate complete device -- and a much
    shorter list to test.
    """
    by_member: dict[str, list[int]] = defaultdict(list)
    for f in findings:
        for c in f.carriers:
            by_member[c].append(f.ref_resnum)

    groups: dict[tuple[int, ...], list[str]] = defaultdict(list)
    for member, positions in by_member.items():
        if len(positions) >= min_positions:
            groups[tuple(sorted(positions))].append(member)

    out: list[ComboFinding] = []
    for positions, members in groups.items():
        score = min(1.0, 0.40 + 0.15 * len(positions) + 0.05 * len(members))
        out.append(ComboFinding(
            list(positions), sorted(members), score,
            f"{len(members)} member(s) carry rare residues at {len(positions)} positions "
            f"({', '.join(str(p) for p in positions)}); a candidate complete site rather "
            "than an isolated substitution"))
    return sorted(out, key=lambda c: (-len(c.positions), -len(c.members)))


# ------------------------------------------------------------------- outputs
def _open(out):
    if hasattr(out, "write"):
        return out, False
    return open(out, "w"), True


def write_rare(findings: list[RareFinding], out, delimiter: str = "\t") -> None:
    fh, close = _open(out)
    try:
        fh.write(delimiter.join(RARE_COLUMNS) + "\n")
        for f in findings:
            fh.write(delimiter.join([
                str(f.ref_resnum), f.residue, str(f.n_carriers), f"{f.fraction:.4f}",
                f"{f.dist:.1f}", f.zone, f"{f.score:.3f}", f.modal, f"{f.modal_frac:.3f}",
                ",".join(f.carriers[:20]), f.reason.replace(delimiter, " ")]) + "\n")
    finally:
        if close:
            fh.close()


def write_combos(combos: list[ComboFinding], out, delimiter: str = "\t") -> None:
    fh, close = _open(out)
    try:
        fh.write(delimiter.join(COMBO_COLUMNS) + "\n")
        for c in combos:
            fh.write(delimiter.join([
                str(len(c.positions)), ",".join(str(p) for p in c.positions),
                str(len(c.members)), ",".join(c.members[:20]), f"{c.score:.3f}",
                c.reason.replace(delimiter, " ")]) + "\n")
    finally:
        if close:
            fh.close()


def emit_features(findings: list[RareFinding], combos: list[ComboFinding]) -> list[str]:
    lines: list[str] = []
    if findings:
        top = findings[0]
        lines.append(
            f"foldtrace prime adjudicate --register register.json \\\n"
            f"  --feature {f'{top.residue} present at reference position {top.ref_resnum} in only {top.n_carriers} members ({top.fraction:.1%}), usually {top.modal}'!r} \\\n"
            f"  --surfaced-by rare_residue --status TODO --permitted TODO")
    if combos:
        c = combos[0]
        lines.append(
            f"foldtrace prime adjudicate --register register.json \\\n"
            f"  --feature {f'{len(c.members)} members carry rare residues at {len(c.positions)} site positions simultaneously'!r} \\\n"
            f"  --surfaced-by rare_residue --status TODO --permitted TODO")
    return lines
