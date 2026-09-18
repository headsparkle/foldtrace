"""Discordance detection: rank disagreements between independent assertions.

Every annotation-error and hidden-biology result in the companion studies came
from two records that should have agreed and did not -- a product string against
a mapped fold, a sequence-derived state against a three-dimensional one, a
monomer against an oligomer, a predicted model against deposited coordinates.
This module makes that search systematic instead of incidental.

Three detectors, all offline:

``annotation``
    A mapped fold and site state against what the record's name and
    cross-references claim. Finds names that assert a different family, intact
    sites carrying no family annotation at all, and -- demoted on purpose --
    records whose name is a known alias of the query family with the
    cross-reference simply missing.

``paired``
    The same targets in two calls tables. Any pair of channels works: two seeds,
    sequence-mapped against structure-mapped, monomer against oligomer,
    predicted against experimental.

``knockout``
    Within one table: the signature of a point knockout on an otherwise intact
    fold, which is how engineered constructs resident in a natural-sequence
    resource announce themselves.

A finding is a ranked row with an explicit reason, in the same style as
``state_reason``. Nothing here decides what a discordance means; it decides what
is worth looking at first.

Stdlib only.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field, asdict

RETAINED, ALTERED, LOST, UNRESOLVED = "retained", "altered", "lost", "unresolved"
_PRESENT = (RETAINED, ALTERED)

# vocabulary classes recognised in a --vocab file
VOCAB_CLASSES = ("family_name", "family_accession", "alias", "foreign")

# residues a point knockout usually installs
_KNOCKOUT_RESIDUES = frozenset("ASG")


class DiscordError(Exception):
    """Raised for malformed inputs."""


@dataclass
class Finding:
    kind: str
    score: float
    target: str
    observation: str
    reason: str
    evidence: dict = field(default_factory=dict)

    def row(self) -> dict:
        d = asdict(self)
        d["score"] = f"{self.score:.3f}"
        d["evidence"] = "; ".join(f"{k}={v}" for k, v in sorted(self.evidence.items()))
        return d


FINDING_COLUMNS = ["kind", "score", "target", "observation", "reason", "evidence"]


# ------------------------------------------------------------------- loading
def _read_tsv(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r and not (list(r.values())[0] or "").startswith("#")]
    if not rows:
        raise DiscordError(f"no rows in {path}")
    return rows


def _target_of(row: dict) -> str:
    for key in ("candidate_id", "candidate", "target", "accession"):
        v = (row.get(key) or "").strip()
        if v:
            return v
    raise DiscordError("no candidate_id/candidate/target column found")


def site_labels(rows: list[dict]) -> list[str]:
    """Site labels present in a calls table, from the ``<label>_state`` columns."""
    return [k[: -len("_state")] for k in rows[0] if k.endswith("_state")]


@dataclass
class Calls:
    """A calls table keyed by target."""
    by_target: dict[str, dict]
    labels: list[str]

    @property
    def targets(self) -> set[str]:
        return set(self.by_target)


def load_calls(path: str) -> Calls:
    rows = _read_tsv(path)
    return Calls({_target_of(r): r for r in rows}, site_labels(rows))


def load_annotations(path: str) -> dict[str, dict]:
    """Annotation metadata per target: name, pfam, clan, taxon, lineage (all optional)."""
    return {_target_of(r): r for r in _read_tsv(path)}


def load_vocab(path: str) -> dict[str, list[tuple[str, str]]]:
    """Read a vocabulary file: ``class``, ``term``, optional ``note``.

    ``family_name``       a name substring meaning the record already says it is in the query family
    ``family_accession``  a Pfam/clan accession meaning the same
    ``alias``             a name that is a *known* alias of the query family; a missing
                          cross-reference here is an annotation gap, not a discovery, and is
                          scored down accordingly
    ``foreign``           a name asserting a different family; note may be its label
    """
    vocab: dict[str, list[tuple[str, str]]] = {c: [] for c in VOCAB_CLASSES}
    for r in _read_tsv(path):
        cls = (r.get("class") or "").strip()
        term = (r.get("term") or "").strip()
        if cls not in VOCAB_CLASSES:
            raise DiscordError(f"unknown vocab class {cls!r}; expected one of {VOCAB_CLASSES}")
        if term:
            vocab[cls].append((term.lower(), (r.get("note") or "").strip()))
    return vocab


# -------------------------------------------------------------------- helpers
def _f(row: dict, key: str) -> float | None:
    v = (row.get(key) or "").strip()
    if v in ("", "NA", "None"):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def _tm(row: dict) -> float | None:
    return _f(row, "tmalign_tm_norm_ref") or _f(row, "tm_norm_ref")


def _fold_ok(row: dict) -> bool:
    return (row.get("fold_ok") or "").strip().lower() == "true"


def _verdict(row: dict) -> str:
    return (row.get("verdict") or UNRESOLVED).strip().lower()


def _hits(text: str, terms: list[tuple[str, str]]) -> list[tuple[str, str]]:
    low = (text or "").lower()
    return [(t, note) for t, note in terms if t and t in low]


def _accession_hits(row: dict, terms: list[tuple[str, str]]) -> list[tuple[str, str]]:
    blob = " ".join((row.get(k) or "") for k in ("pfam", "clan", "interpro", "xref")).lower()
    return [(t, note) for t, note in terms if t and t in blob]


def _confidence(tm: float | None) -> float:
    """Map fold similarity onto 0-1. A discordance on a marginal fold is worth less."""
    if tm is None:
        return 0.5
    return max(0.0, min(1.0, (tm - 0.5) / 0.4))


# ------------------------------------------------------------ 1. annotation
def scan_annotation(calls: Calls, annos: dict[str, dict],
                    vocab: dict[str, list[tuple[str, str]]]) -> list[Finding]:
    """Mapped fold and site state against what the record's own metadata claims."""
    out: list[Finding] = []
    for target, row in calls.by_target.items():
        if not _fold_ok(row):
            continue
        verdict = _verdict(row)
        if verdict not in _PRESENT:
            continue
        anno = annos.get(target, {})
        name = (anno.get("name") or anno.get("protein_name") or "").strip()
        tm = _tm(row)
        conf = _confidence(tm)
        ev = {"tm": f"{tm:.3f}" if tm is not None else "NA", "verdict": verdict,
              "name": name or "(none)"}

        foreign = _hits(name, vocab.get("foreign", []))
        fam_name = _hits(name, vocab.get("family_name", []))
        fam_acc = _accession_hits(anno, vocab.get("family_accession", []))
        alias = _hits(name, vocab.get("alias", []))

        if foreign and not fam_acc:
            term, label = foreign[0]
            out.append(Finding(
                "name_contradicts_fold", 0.60 + 0.35 * conf, target,
                f"name asserts {label or term!r} but the {verdict} query site is mapped on this fold",
                f"name matched foreign term {term!r}; no query-family accession on the record; "
                f"fold_ok with TM {ev['tm']}", ev))
            continue

        if fam_acc or fam_name:
            continue  # the record already says what the structure says

        if alias:
            term, note = alias[0]
            out.append(Finding(
                "annotation_gap", 0.10 + 0.10 * conf, target,
                f"name {term!r} is a known alias of the query family; cross-reference absent",
                f"alias {term!r}{' -> ' + note if note else ''}; recovering these is incomplete "
                "per-record cross-referencing, not family discovery", ev))
            continue

        out.append(Finding(
            "unannotated_intact_site", 0.45 + 0.35 * conf, target,
            f"{verdict} query site on the mapped fold, with no family name or accession",
            f"no family term, accession or known alias matched; fold_ok with TM {ev['tm']}", ev))
    return out


# ---------------------------------------------------------------- 2. paired
def scan_paired(a: Calls, b: Calls, label_a: str = "A", label_b: str = "B") -> list[Finding]:
    """Same target, two channels, different answer."""
    out: list[Finding] = []
    shared = sorted(a.targets & b.targets)
    for target in shared:
        ra, rb = a.by_target[target], b.by_target[target]
        va, vb = _verdict(ra), _verdict(rb)
        fa, fb = _fold_ok(ra), _fold_ok(rb)
        conf = max(_confidence(_tm(ra)), _confidence(_tm(rb)))
        ev = {label_a: va, label_b: vb, "tm_" + label_a: (f"{_tm(ra):.3f}" if _tm(ra) else "NA"),
              "tm_" + label_b: (f"{_tm(rb):.3f}" if _tm(rb) else "NA")}

        if fa != fb:
            out.append(Finding(
                "fold_disagreement", 0.50 + 0.25 * conf, target,
                f"fold_ok {label_a}={str(fa).lower()} {label_b}={str(fb).lower()}",
                "one channel places this target in the fold and the other does not", ev))

        if va != vb:
            hard = {va, vb} <= {RETAINED, ALTERED, LOST}
            kind = "state_flip" if hard else "resolution_disagreement"
            score = (0.70 + 0.30 * conf) if hard else (0.35 + 0.25 * conf)
            out.append(Finding(
                kind, score, target, f"{label_a} says {va}, {label_b} says {vb}",
                ("a hard state flip between two independent channels; at most one can be right"
                 if hard else
                 "one channel resolved this target and the other abstained; check the "
                 "abstention reason before treating either as the answer"),
                dict(ev, reason_a=(ra.get("state_reason") or "")[:120],
                     reason_b=(rb.get("state_reason") or "")[:120])))

        for label in sorted(set(a.labels) & set(b.labels)):
            sa = (ra.get(f"{label}_state") or "").strip().lower()
            sb = (rb.get(f"{label}_state") or "").strip().lower()
            if not (sa and sb) or sa == sb or not {sa, sb} <= {RETAINED, ALTERED, LOST}:
                continue
            if (sa, sb) == (va, vb):
                # this site flip is what the verdict flip already reports; emitting both
                # would double-count one disagreement. A site flip the verdict masks --
                # two sites moving in opposite directions -- still gets its own row.
                continue
            if True:
                out.append(Finding(
                    "site_flip", 0.55 + 0.30 * conf, target,
                    f"site {label}: {label_a}={sa}, {label_b}={sb}",
                    f"per-site disagreement at {label} between two independent channels",
                    dict(ev, site=label,
                         obs_a=ra.get(f"{label}_obs", ""), obs_b=rb.get(f"{label}_obs", ""))))
    return out


# -------------------------------------------------------------- 3. knockout
_OBS_RE = re.compile(r"^([A-Za-z])(\d+)$")


def scan_knockout(calls: Calls, min_tm: float = 0.85,
                  max_prevalence: float = 0.05) -> list[Finding]:
    """The signature of a point knockout on an otherwise intact fold.

    One site lost to Ala, Ser or Gly, every other site retained, and a close
    superposition. In a natural-sequence resource this is usually a deposited
    catalysis-deficient construct whose sequence inherited the engineered
    residue -- simultaneously an inadvertent positive control and an annotation
    error. It is also, occasionally, a real single-residue loss, which is why
    this is a flag and not a verdict.

    Prevalence is what separates the two, and it is checked here rather than
    left to the reader. Engineered constructs are one-offs: two in three hundred
    alpha/beta-hydrolase hits. A fold in which a fifth of all members show the
    same single-ligand loss is describing its own biology, not a lab artefact.
    Above ``max_prevalence`` the individual calls are replaced by one
    ``recurrent_single_site_loss`` row, which is the more useful finding: it
    names a position this fold routinely gives up.
    """
    out: list[Finding] = []
    if len(calls.labels) < 2:
        return out
    eligible = sum(1 for r in calls.by_target.values()
                   if _fold_ok(r) and (_tm(r) or 0.0) >= min_tm)
    for target, row in calls.by_target.items():
        tm = _tm(row)
        if not _fold_ok(row) or tm is None or tm < min_tm:
            continue
        states = {lab: (row.get(f"{lab}_state") or "").strip().lower() for lab in calls.labels}
        lost = [lab for lab, s in states.items() if s == LOST]
        if len(lost) != 1:
            continue
        if any(s != RETAINED for lab, s in states.items() if lab != lost[0]):
            continue
        lab = lost[0]
        m = _OBS_RE.match((row.get(f"{lab}_obs") or "").strip())
        if not m or m.group(1).upper() not in _KNOCKOUT_RESIDUES:
            continue
        obs = m.group(1).upper()
        out.append(Finding(
            "point_knockout_signature", 0.55 + 0.40 * _confidence(tm), target,
            f"{lab} lost to {obs}{m.group(2)}; all other sites retained at TM {tm:.3f}",
            f"single-site loss to {obs} on an otherwise intact site is the signature of an "
            "engineered construct resident in a natural-sequence record; check the source entry "
            "before treating it as a natural loss",
            {"tm": f"{tm:.3f}", "site": lab, "observed": f"{obs}{m.group(2)}"}))

    # prevalence check: a signature this common is the fold's biology, not a construct
    by_site: dict[str, list[Finding]] = {}
    for f in out:
        by_site.setdefault(f.evidence["site"], []).append(f)
    kept: list[Finding] = []
    for lab, group in by_site.items():
        frac = len(group) / eligible if eligible else 0.0
        # a prevalence needs a denominator; on a handful of targets it means nothing
        if eligible < 20 or frac <= max_prevalence:
            kept.extend(group)
            continue
        examples = ", ".join(f.target for f in sorted(group, key=lambda x: x.target)[:3])
        kept.append(Finding(
            "recurrent_single_site_loss", 0.55 + 0.25 * min(1.0, frac / 0.5),
            f"{len(group)} targets",
            f"site {lab} is singly lost to Ala/Ser/Gly in {len(group)} of {eligible} "
            f"({frac:.0%}) while the other sites are retained",
            "too common to be engineered constructs; this is a position the fold routinely "
            "gives up while keeping the rest of the site, which is a statement about the "
            f"family rather than about any one record (e.g. {examples})",
            {"site": lab, "n": len(group), "eligible": eligible, "fraction": f"{frac:.3f}"}))
    return kept


# ------------------------------------------------------- 4. sequence identity
def scan_identity(profile_rows, annotations: dict[str, dict],
                  min_identity: float = 0.98, organism_key: str = "organism",
                  sample: int = 24) -> list[Finding]:
    """Near-identical sequences attributed to different organisms.

    A natural-sequence resource should not contain the same protein under two
    unrelated taxa. When it does, the usual cause is a reporter construct
    carried into an assembly: an EGFP or mKate2 sequence filed under whatever
    organism was being engineered. Finding these by eye is how the problem is
    normally caught, and it does not scale.

    Candidates are bucketed on a sample of positions before any pairwise
    comparison, so this stays usable on a full hit set while still catching
    matches short of exact identity -- the published cases sit at 99.6%, not
    100%.

    Takes the row objects from ``foldtrace profile``.
    """
    seqs: dict[str, dict[int, str]] = {}
    for o in profile_rows:
        if getattr(o, "resolved", False) and o.obs_res not in ("-", "X"):
            seqs.setdefault(o.candidate, {})[o.ref_resnum] = o.obs_res
    if len(seqs) < 2:
        return []

    positions = sorted({rn for d in seqs.values() for rn in d})
    if not positions:
        return []
    # Several independent bands rather than one probe set: a single substitution
    # inside a sample would otherwise put two 99%-identical records in different
    # buckets and hide exactly the case this exists for.
    bands, per_band = 4, max(3, sample // 4)
    probes = [positions[b::bands][:per_band] for b in range(bands)]

    buckets: dict[tuple, list[str]] = {}
    for cand, d in seqs.items():
        for bi, probe in enumerate(probes):
            if probe:
                buckets.setdefault((bi,) + tuple(d.get(p, "-") for p in probe), []).append(cand)

    def identity(a: str, b: str) -> float:
        da, db = seqs[a], seqs[b]
        shared = da.keys() & db.keys()
        if not shared:
            return 0.0
        return sum(da[p] == db[p] for p in shared) / len(shared)

    out: list[Finding] = []
    seen: set[frozenset] = set()
    for members in buckets.values():
        if len(members) < 2:
            continue
        members = sorted(members)
        group = [members[0]]
        for cand in members[1:]:
            if identity(members[0], cand) >= min_identity:
                group.append(cand)
        if len(group) < 2 or frozenset(group) in seen:
            continue
        seen.add(frozenset(group))
        orgs = {(annotations.get(c, {}).get(organism_key) or "").strip()
                for c in group}
        orgs.discard("")
        if len(orgs) < 2:
            continue
        out.append(Finding(
            "identical_sequence_across_organisms",
            0.70 + 0.25 * min(1.0, len(orgs) / 4), ", ".join(group[:6]),
            f"{len(group)} records are at least {min_identity:.0%} identical across mapped "
            f"positions but are attributed to {len(orgs)} organisms: {', '.join(sorted(orgs)[:5])}",
            "the same sequence under unrelated taxa is the signature of a reporter construct "
            "carried into an assembly; check these against a reference set of engineered "
            "proteins before treating any of them as a natural record",
            {"n_records": len(group), "n_organisms": len(orgs),
             "organisms": "|".join(sorted(orgs)[:5])}))
    return out


# ------------------------------------------------------------------ assembly
def rank(findings: list[Finding]) -> list[Finding]:
    return sorted(findings, key=lambda f: (-f.score, f.kind, f.target))


def summarise(findings: list[Finding]) -> list[str]:
    counts: dict[str, int] = {}
    for f in findings:
        counts[f.kind] = counts.get(f.kind, 0) + 1
    return [f"{n:>6}  {k}" for k, n in sorted(counts.items(), key=lambda kv: -kv[1])]


def write(findings: list[Finding], out, delimiter: str = "\t") -> None:
    close = False
    if hasattr(out, "write"):
        fh = out
    else:
        fh, close = open(out, "w"), True
    try:
        fh.write(delimiter.join(FINDING_COLUMNS) + "\n")
        for f in findings:
            r = f.row()
            fh.write(delimiter.join(str(r[c]).replace(delimiter, " ") for c in FINDING_COLUMNS) + "\n")
    finally:
        if close:
            fh.close()


# Kinds that describe a known annotation gap rather than a candidate finding.
# Emitting these as features would put a non-novelty result in front of the gate.
_NOT_A_FEATURE = {"annotation_gap"}


def emit_features(findings: list[Finding], top: int = 5) -> list[str]:
    """Suggest `foldtrace prime adjudicate` calls for the leading findings.

    Each kind becomes one feature, because prior art is checked per phenomenon,
    not per protein. Every emitted feature must be adjudicated before the
    novelty gate will pass a draft that claims it.
    """
    groups: dict[str, list[Finding]] = {}
    for f in findings:
        if f.kind in _NOT_A_FEATURE:
            continue
        groups.setdefault(f.kind, []).append(f)
    lines: list[str] = []
    for kind, fs in sorted(groups.items(), key=lambda kv: -max(f.score for f in kv[1]))[:top]:
        ex = ", ".join(f.target for f in rank(fs)[:3])
        noun = "target" if len(fs) == 1 else "targets"
        lines.append(
            f"foldtrace prime adjudicate --register register.json \\\n"
            f"  --feature {f'{len(fs)} {noun} showing {kind} (e.g. {ex})'!r} \\\n"
            f"  --surfaced-by discordance --status TODO --permitted TODO")
    return lines
