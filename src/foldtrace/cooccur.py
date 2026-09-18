"""Co-occurrence: features that travel with a fold more often than they should.

Nothing in the declared-site map would have found TIR-FGE. The signal was that
an unrelated domain kept turning up attached to the hits, and somebody noticed.
This module does the noticing: given any per-target feature table -- domain
architecture, genomic-neighbourhood family, taxon, host, whatever was
annotated -- it counts how often each feature appears in a focus group against a
background and ranks the excess.

Two warnings are built into the output rather than left to the writer.

Enumerating features from the same dataset you then test them in is selection,
and the resulting p-value is exploratory however small it is. Every finding says
so unless a held-out background is supplied with ``--background``, and the
summary records which mode was used. This is the caveat the HEPN neighbourhood
analysis had to carry in prose three times.

Taxonomic clustering imitates enrichment. If the carriers are one genus, the
feature may be tracking ancestry rather than the fold, so each finding reports
how many distinct values of a grouping column its carriers span when one is
given.

Exact Fisher test implemented on ``math.comb``; no scipy.
"""
from __future__ import annotations

import csv
import math
from collections import defaultdict
from dataclasses import dataclass, field

FEATURE_COLUMNS = ["feature_type", "feature", "focus_n", "focus_total", "focus_frac",
                   "bg_n", "bg_total", "bg_frac", "odds_ratio", "p_value", "q_value",
                   "spread", "score", "evidence_mode", "carriers", "reason"]


class CooccurError(Exception):
    """Raised for malformed or insufficient input."""


@dataclass
class FeatureFinding:
    feature_type: str
    feature: str
    focus_n: int
    focus_total: int
    bg_n: int
    bg_total: int
    odds_ratio: float
    p_value: float
    q_value: float
    evidence_mode: str
    spread: int | None = None
    carriers: list[str] = field(default_factory=list)
    score: float = 0.0
    reason: str = ""

    @property
    def focus_frac(self) -> float:
        return self.focus_n / self.focus_total if self.focus_total else 0.0

    @property
    def bg_frac(self) -> float:
        return self.bg_n / self.bg_total if self.bg_total else 0.0


# --------------------------------------------------------------- statistics
def _hypergeom_pmf(k: int, K: int, n: int, N: int) -> float:
    """P(X = k) drawing n from N with K successes."""
    if k < 0 or k > K or n - k < 0 or n - k > N - K:
        return 0.0
    return math.comb(K, k) * math.comb(N - K, n - k) / math.comb(N, n)


def fisher_exact(a: int, b: int, c: int, d: int,
                 alternative: str = "greater") -> tuple[float, float]:
    """Fisher exact test on [[a, b], [c, d]]; returns (odds ratio, p).

    Table: a = focus with feature, b = focus without, c = background with,
    d = background without. The odds ratio uses a Haldane-Anscombe 0.5
    correction so a zero cell does not produce an infinity.

    ``alternative`` defaults to ``greater``. This module only ever reports
    features *enriched* in the focus group -- ``min_odds`` discards the other
    direction before the test is read -- so charging a two-sided penalty for a
    one-sided question is inconsistent, and on the sets this is used for it
    costs real signal: two independently built controls put a planted
    co-occurrence at p = 0.056 and p = 0.081 two-sided, both just outside a
    calibrated cutoff, and at 0.028 and 0.045 one-sided. Pass ``two-sided``
    when depletion matters as much as enrichment.
    """
    N, K, n = a + b + c + d, a + c, a + b
    if N == 0 or K == 0 or n == 0 or K == N or n == N:
        return 1.0, 1.0
    orat = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    lo, hi = max(0, n - (N - K)), min(K, n)
    if alternative == "greater":
        p = sum(_hypergeom_pmf(k, K, n, N) for k in range(a, hi + 1))
    elif alternative == "less":
        p = sum(_hypergeom_pmf(k, K, n, N) for k in range(lo, a + 1))
    elif alternative == "two-sided":
        observed = _hypergeom_pmf(a, K, n, N)
        p = sum(_hypergeom_pmf(k, K, n, N) for k in range(lo, hi + 1)
                if _hypergeom_pmf(k, K, n, N) <= observed * (1 + 1e-9))
    else:
        raise CooccurError(f"unknown alternative {alternative!r}")
    return orat, min(1.0, p)


def benjamini_hochberg(pvals: list[float]) -> list[float]:
    """BH-adjusted q values, order preserved."""
    m = len(pvals)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        idx = m - rank + 1
        prev = min(prev, pvals[i] * m / idx)
        q[i] = prev
    return q


# ------------------------------------------------------------------ loading
def load_features(path: str) -> dict[str, list[tuple[str, str]]]:
    """Long-format feature table: ``target``, ``feature_type``, ``feature``."""
    out: dict[str, list[tuple[str, str]]] = defaultdict(list)
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise CooccurError(f"no rows in {path}")
    for r in rows:
        t = (r.get("target") or r.get("candidate") or r.get("candidate_id") or "").strip()
        ft = (r.get("feature_type") or "feature").strip()
        f = (r.get("feature") or "").strip()
        if t and f:
            out[t].append((ft, f))
    if not out:
        raise CooccurError(f"{path}: need target and feature columns")
    return dict(out)


def load_column(path: str, column: str) -> dict[str, str]:
    """A single per-target column, e.g. genus, used for the spread check."""
    out: dict[str, str] = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            t = (r.get("target") or r.get("candidate") or r.get("candidate_id") or "").strip()
            v = (r.get(column) or "").strip()
            if t and v:
                out[t] = v
    return out


# ------------------------------------------------------------------ analysis
def analyse(features: dict[str, list[tuple[str, str]]], focus: set[str],
            background: set[str] | None = None,
            spread_by: dict[str, str] | None = None,
            min_count: int = 3, min_odds: float = 2.0, max_q: float = 0.02,
            alternative: str = "greater",
            held_out: bool = False) -> tuple[list[FeatureFinding], dict]:
    """Rank features enriched in ``focus`` against ``background``.

    With no explicit background the remaining annotated targets are used, which
    is a within-set comparison and is labelled as such.

    ``max_q`` is the multiple-testing gate. Every feature in the table is tested,
    so on a set of any size a handful will clear an odds-ratio threshold by
    chance; a random split of this control produced four such features before
    this gate existed. Benjamini-Hochberg is applied across all tested features
    and anything above ``max_q`` is dropped rather than merely down-ranked.
    """
    if not focus:
        raise CooccurError("focus set is empty")
    annotated = set(features)
    focus = focus & annotated
    bg = (background & annotated) if background is not None else (annotated - focus)
    if not focus:
        raise CooccurError("no focus target carries any feature")
    mode = "held_out_background" if held_out and background is not None else "within_dataset"

    carriers: dict[tuple[str, str], set[str]] = defaultdict(set)
    for target, feats in features.items():
        for ft, f in feats:
            carriers[(ft, f)].add(target)

    raw: list[FeatureFinding] = []
    for (ft, f), members in carriers.items():
        a = len(members & focus)
        c = len(members & bg)
        if a < min_count:
            continue
        b, d = len(focus) - a, len(bg) - c
        orat, p = fisher_exact(a, b, c, d, alternative=alternative)
        if orat < min_odds:
            continue
        found = sorted(members & focus)
        spread = len({spread_by[m] for m in found if m in spread_by}) if spread_by else None
        raw.append(FeatureFinding(ft, f, a, len(focus), c, len(bg), orat, p, 1.0, mode,
                                  spread, found))

    qs = benjamini_hochberg([f.p_value for f in raw])
    tested = len(raw)
    for f, q in zip(raw, qs):
        f.q_value = q
        prevalence = f.focus_frac
        sig = 1.0 - min(1.0, q / max(max_q, 1e-9))
        f.score = 0.40 * prevalence + 0.35 * sig + 0.25 * min(1.0, math.log10(f.odds_ratio) / 1.5)
        bits = [f"{f.feature} in {f.focus_n}/{f.focus_total} of the focus group "
                f"against {f.bg_n}/{f.bg_total} background "
                f"(odds ratio {f.odds_ratio:.1f}, p = {f.p_value:.2g}, q = {q:.2g})"]
        if f.spread is not None:
            bits.append(f"carriers span {f.spread} distinct grouping values"
                        + ("; a single group means this may track ancestry rather than the fold"
                           if f.spread <= 1 else ""))
        if f.evidence_mode == "within_dataset":
            bits.append("features were enumerated from this dataset, so this is exploratory "
                        "and the p value is not a test of a prespecified hypothesis")
        f.reason = "; ".join(bits)

    kept = [f for f in raw if f.q_value <= max_q]
    suppressed = [f for f in raw if f.q_value > max_q and f.q_value <= min(1.0, max_q * 5)]

    summary = {
        "focus_n": len(focus), "background_n": len(bg),
        "features_considered": len(carriers), "features_tested": tested,
        "features_reported": len(kept), "max_q": f"{max_q:.3f}",
        "alternative": alternative,
        "features_suppressed_near_cutoff": len(suppressed),
        "evidence_mode": mode,
        "note": ("held-out background supplied" if mode == "held_out_background"
                 else "WITHIN-DATASET comparison: every result here is exploratory"),
    }
    if suppressed:
        # Silently dropping a near-miss is the same error as reporting noise.
        # On a small set the q gate costs real signal: a feature at 8/24 against
        # 1/17 sits just the wrong side of it.
        summary["suppressed_detail"] = "; ".join(
            f"{f.feature} (q = {f.q_value:.3g}, {f.focus_n}/{f.focus_total} vs "
            f"{f.bg_n}/{f.bg_total})" for f in sorted(suppressed, key=lambda x: x.q_value)[:5])
        summary["suppressed_note"] = (
            "these fell just outside --max-q; on a small feature table that is a power "
            "limit, not evidence of absence. Re-run with a looser --max-q to see them.")
    return sorted(kept, key=lambda f: (-f.score, f.feature)), summary


# ------------------------------------------------------------------- outputs
def write(findings: list[FeatureFinding], out, delimiter: str = "\t") -> None:
    close = False
    if hasattr(out, "write"):
        fh = out
    else:
        fh, close = open(out, "w"), True
    try:
        fh.write(delimiter.join(FEATURE_COLUMNS) + "\n")
        for f in findings:
            fh.write(delimiter.join([
                f.feature_type, f.feature, str(f.focus_n), str(f.focus_total),
                f"{f.focus_frac:.3f}", str(f.bg_n), str(f.bg_total), f"{f.bg_frac:.3f}",
                f"{f.odds_ratio:.2f}", f"{f.p_value:.3g}", f"{f.q_value:.3g}",
                str(f.spread) if f.spread is not None else "NA",
                f"{f.score:.3f}", f.evidence_mode, ",".join(f.carriers[:20]),
                f.reason.replace(delimiter, " ")]) + "\n")
    finally:
        if close:
            fh.close()


def emit_features(findings: list[FeatureFinding], top: int = 2) -> list[str]:
    lines = []
    for f in findings[:top]:
        lines.append(
            f"foldtrace prime adjudicate --register register.json \\\n"
            f"  --feature {f'{f.feature} co-occurs with the fold in {f.focus_n} of {f.focus_total} members (odds ratio {f.odds_ratio:.1f})'!r} \\\n"
            f"  --surfaced-by cooccurrence --status TODO --permitted TODO")
    return lines
