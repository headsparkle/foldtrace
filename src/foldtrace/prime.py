"""Prior-art register: what was known about a fold before the search was run.

The register has two jobs. Site definitions are built from cited experimental
evidence rather than recall, and a novelty claim cannot be drafted until the
feature it concerns has been adjudicated against a record that predates the
results.

Two invariants are enforced here rather than left to discipline:

1. Nothing in this module can mark a claim ``verified``. Promotion happens only
   in :func:`apply_verification`, from a verification file produced by an
   external retrieval tool. A register whose rows were written but never
   retrieved is worse than no register, because it licenses exactly the novelty
   claim it exists to block.

2. The register records what is *known*, never what is expected. There is no
   field for a hypothesis, a target chemistry, or a candidate ranking, and the
   register is advisory to :mod:`foldtrace.pipeline` only -- it must not change
   retrieval, thresholds, scoring or state calls.

Stdlib only; no new runtime dependencies.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone

REGISTER_VERSION = "0.1"

CLAIM_TYPES = (
    "mechanism",
    "residue_evidence",
    "described_family",
    "described_architecture",
    "characterized_pseudoenzyme",
    "existing_structure",
    "annotation_caveat",
)

EVIDENCE_STRENGTHS = ("experimental", "structural", "computational", "annotation")

VERIFICATION_STATES = ("unverified", "verified", "failed")

SURFACED_BY = (
    "cooccurrence",
    "discordance",
    "rare_residue",
    "unexpected_invariant",
    "stratification",
    "manual",
)

PRIOR_ART_STATUS = (
    "already_described",
    "partially_described",
    "no_prior_description_found",
)

_RESIDUE_RE = re.compile(r"^([A-Z])(\d+)$")


class PrimeError(Exception):
    """Raised for malformed registers, rejected freezes and gate violations."""


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ---------------------------------------------------------------------- model
@dataclass
class Claim:
    """One thing the literature already says about this fold."""
    claim_id: str
    claim_type: str
    claim: str
    identifier: str | None = None
    identifier_resolved: bool = False
    source_verified: str = "unverified"
    verification_method: str | None = None
    residues: list[str] = field(default_factory=list)
    evidence_strength: str = "annotation"
    discovered_post_run: bool = False
    trigger: str | None = None
    note: str = ""

    def __post_init__(self) -> None:
        # a JSON null for an optional list is normal in hand-written registers
        if self.residues is None:
            self.residues = []

    def validate(self) -> list[str]:
        p: list[str] = []
        if self.claim_type not in CLAIM_TYPES:
            p.append(f"{self.claim_id}: claim_type {self.claim_type!r} not in {CLAIM_TYPES}")
        if self.evidence_strength not in EVIDENCE_STRENGTHS:
            p.append(f"{self.claim_id}: evidence_strength {self.evidence_strength!r} invalid")
        if self.source_verified not in VERIFICATION_STATES:
            p.append(f"{self.claim_id}: source_verified {self.source_verified!r} invalid")
        if not self.claim.strip():
            p.append(f"{self.claim_id}: empty claim text")
        if self.discovered_post_run and not (self.trigger or "").strip():
            p.append(f"{self.claim_id}: post-run claim has no trigger")
        for r in self.residues:
            if not _RESIDUE_RE.match(r):
                p.append(f"{self.claim_id}: residue {r!r} not in <one-letter><number> form, e.g. H102")
        return p


@dataclass
class Adjudication:
    """One emergent feature, checked against the frozen record."""
    feature_id: str
    feature: str
    surfaced_by: str
    prior_art_status: str
    permitted_claim: str
    retrieval_date: str | None = None
    described_component: str | None = None
    undescribed_component: str | None = None
    supporting_claims: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.supporting_claims is None:
            self.supporting_claims = []

    def validate(self) -> list[str]:
        p: list[str] = []
        if self.surfaced_by not in SURFACED_BY:
            p.append(f"{self.feature_id}: surfaced_by {self.surfaced_by!r} invalid")
        if self.prior_art_status not in PRIOR_ART_STATUS:
            p.append(f"{self.feature_id}: prior_art_status {self.prior_art_status!r} invalid")
        if self.prior_art_status == "partially_described":
            # The Kibby case: forcing both halves at adjudication time is what
            # produces "121 compact fusions undescribed, 44 already published"
            # in the first draft rather than in revision.
            if not (self.described_component or "").strip():
                p.append(f"{self.feature_id}: partially_described requires described_component")
            if not (self.undescribed_component or "").strip():
                p.append(f"{self.feature_id}: partially_described requires undescribed_component")
        if not self.permitted_claim.strip():
            p.append(f"{self.feature_id}: empty permitted_claim")
        return p


@dataclass
class Register:
    fold_id: str
    reference_structure: str
    reference_source: str | None = None
    scope_statement: str = ""
    novelty_criteria: list[str] = field(default_factory=list)
    non_novelty_criteria: list[str] = field(default_factory=list)
    retrieval_date: str | None = None
    retrieval_tools: list[str] = field(default_factory=list)
    register_version: str = REGISTER_VERSION
    frozen_at: str | None = None
    sha256: str | None = None
    claims: list[Claim] = field(default_factory=list)
    adjudications: list[Adjudication] = field(default_factory=list)

    def __post_init__(self) -> None:
        for name in ("novelty_criteria", "non_novelty_criteria", "retrieval_tools",
                     "claims", "adjudications"):
            if getattr(self, name) is None:
                setattr(self, name, [])

    # ------------------------------------------------------------- accessors
    def claim(self, claim_id: str) -> Claim:
        for c in self.claims:
            if c.claim_id == claim_id:
                return c
        raise PrimeError(f"no claim {claim_id!r} in register")

    def next_claim_id(self) -> str:
        n = 0
        for c in self.claims:
            m = re.match(r"^PA-(\d+)$", c.claim_id)
            if m:
                n = max(n, int(m.group(1)))
        return f"PA-{n + 1:03d}"

    def next_feature_id(self) -> str:
        n = 0
        for a in self.adjudications:
            m = re.match(r"^EF-(\d+)$", a.feature_id)
            if m:
                n = max(n, int(m.group(1)))
        return f"EF-{n + 1:03d}"

    def pre_run_claims(self) -> list[Claim]:
        return [c for c in self.claims if not c.discovered_post_run]


# ------------------------------------------------------------------ load/save
def new_register(fold_id: str, reference_structure: str, scope: str = "",
                 tools: list[str] | None = None) -> Register:
    return Register(fold_id=fold_id, reference_structure=reference_structure,
                    scope_statement=scope, retrieval_tools=list(tools or []),
                    retrieval_date=_now())


def load(path: str) -> Register:
    with open(path) as fh:
        raw = json.load(fh)
    raw.pop("_provenance_note", None)
    claims = [Claim(**c) for c in raw.pop("claims", [])]
    adjs = [Adjudication(**a) for a in raw.pop("adjudications", [])]
    known = {f for f in Register.__dataclass_fields__ if f not in ("claims", "adjudications")}
    unknown = set(raw) - known
    if unknown:
        raise PrimeError(f"unknown register fields: {sorted(unknown)}")
    return Register(claims=claims, adjudications=adjs, **raw)


def save(reg: Register, path: str) -> None:
    d = asdict(reg)
    with open(path, "w") as fh:
        json.dump(d, fh, indent=2, sort_keys=False)
        fh.write("\n")


# --------------------------------------------------------------------- hashing
def canonical_body(reg: Register) -> str:
    """Canonical JSON of the pre-run body: the part the freeze hash covers.

    Appended adjudications and post-run claims sit outside the hash and carry
    their own timestamps, so the record of what was known *before* the run stays
    verifiable after the run has added to it.
    """
    body = {
        "fold_id": reg.fold_id,
        "reference_structure": reg.reference_structure,
        "reference_source": reg.reference_source,
        "scope_statement": reg.scope_statement,
        "novelty_criteria": reg.novelty_criteria,
        "non_novelty_criteria": reg.non_novelty_criteria,
        "retrieval_date": reg.retrieval_date,
        "retrieval_tools": reg.retrieval_tools,
        "register_version": reg.register_version,
        "claims": [asdict(c) for c in reg.pre_run_claims()],
    }
    return json.dumps(body, sort_keys=True, separators=(",", ":"))


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------- verification
def apply_verification(reg: Register, records: list[dict]) -> tuple[int, int]:
    """Promote claims from a verification file produced by an external tool.

    Each record needs ``claim_id`` and ``status`` (``verified`` or ``failed``),
    plus ``method`` describing how the identifier was resolved. This module
    never fetches anything itself, which is what keeps the package dependency-
    free and makes the retrieval tool pluggable -- and what makes it impossible
    for the register writer to also be the register verifier.
    """
    ok = failed = 0
    for rec in records:
        cid = rec.get("claim_id")
        status = rec.get("status")
        if status not in ("verified", "failed"):
            raise PrimeError(f"{cid}: verification status must be verified|failed, got {status!r}")
        c = reg.claim(cid)
        method = (rec.get("method") or "").strip()
        if status == "verified" and not method:
            raise PrimeError(f"{cid}: cannot verify without a method string")
        c.source_verified = status
        c.verification_method = method or None
        c.identifier_resolved = status == "verified"
        if rec.get("identifier"):
            c.identifier = rec["identifier"]
        ok += status == "verified"
        failed += status == "failed"
    return ok, failed


# ----------------------------------------------------------------- validation
def validate(reg: Register, sites_path: str | None = None) -> list[str]:
    """Return every reason this register is not fit to freeze. Empty means fit."""
    problems: list[str] = []
    seen: set[str] = set()
    for c in reg.claims:
        if c.claim_id in seen:
            problems.append(f"duplicate claim_id {c.claim_id}")
        seen.add(c.claim_id)
        problems.extend(c.validate())
    fseen: set[str] = set()
    for a in reg.adjudications:
        if a.feature_id in fseen:
            problems.append(f"duplicate feature_id {a.feature_id}")
        fseen.add(a.feature_id)
        problems.extend(a.validate())
        for cid in a.supporting_claims:
            if cid not in seen:
                problems.append(f"{a.feature_id}: supporting claim {cid} not in register")

    if not reg.novelty_criteria:
        problems.append("novelty_criteria is empty; it must be written before the run")
    if not reg.non_novelty_criteria:
        problems.append(
            "non_novelty_criteria is empty; state what would NOT count as new "
            "however it is recovered (this is the field that failed for DUF86/PF01934)")
    if not reg.scope_statement.strip():
        problems.append("scope_statement is empty")

    pre = reg.pre_run_claims()
    if not pre:
        problems.append("register has no pre-run claims")
    for c in pre:
        if c.source_verified == "unverified":
            problems.append(f"{c.claim_id}: unverified; every pre-run claim must be "
                            "verified or explicitly marked failed")
        elif c.source_verified == "failed":
            # A documented retrieval failure is provenance worth keeping: it records
            # that the claim was checked and could not be confirmed. It only blocks
            # the freeze if a declared site rests on it.
            if sites_path and c.claim_type == "residue_evidence" and c.residues:
                try:
                    declared = {rn for _, rn in _declared_positions(sites_path)}
                except (ValueError, FileNotFoundError):
                    declared = set()
                for r in c.residues:
                    m = _RESIDUE_RE.match(r)
                    if m and int(m.group(2)) in declared:
                        problems.append(
                            f"{c.claim_id}: verification failed but site {r} is declared; "
                            "either verify the source or remove the position")
    if not any(c.source_verified == "verified" for c in pre):
        problems.append("no pre-run claim was verified; the register records nothing confirmed")

    if sites_path:
        problems.extend(_unbacked_sites(reg, sites_path))
    return problems


def _declared_positions(sites_path: str) -> list[tuple[str, int]]:
    from .io import parse_sites
    return [(s.label, s.resnum) for s in parse_sites(sites_path)]


def _unbacked_sites(reg: Register, sites_path: str) -> list[str]:
    try:
        declared = _declared_positions(sites_path)
    except (ValueError, FileNotFoundError) as exc:
        return [f"sites file {sites_path}: {exc} "
                "(an empty skeleton means no verified residue_evidence claims yet)"]
    backed: set[int] = set()
    for c in reg.claims:
        if c.claim_type != "residue_evidence":
            continue
        for r in c.residues:
            m = _RESIDUE_RE.match(r)
            if m:
                backed.add(int(m.group(2)))
    out = []
    for label, resnum in declared:
        if resnum not in backed:
            out.append(f"site {label} (ref {resnum}) has no backing residue_evidence claim")
    return out


# --------------------------------------------------------------------- freeze
def freeze(reg: Register, sites_path: str, extra_paths: list[str] | None = None,
           out_path: str = "FREEZE.md") -> str:
    """Validate, stamp and write FREEZE.md. Raises PrimeError on any problem."""
    problems = validate(reg, sites_path=sites_path)
    if problems:
        raise PrimeError("register not fit to freeze:\n  - " + "\n  - ".join(problems))

    reg.frozen_at = _now()
    reg.sha256 = sha256_text(canonical_body(reg))

    lines = [
        "# FREEZE",
        "",
        f"Locked: {reg.frozen_at}",
        f"Fold: {reg.fold_id}   Reference: {reg.reference_structure}",
        "",
        "The prior-art register hash covers the pre-run body only: header, novelty and",
        "non-novelty criteria, and all claims with discovered_post_run = false.",
        "Adjudications and post-run claims are appended after this lock and carry their",
        "own timestamps.",
        "",
        "## Hashes",
        "",
        "| File | SHA256 |",
        "| --- | --- |",
        f"| register (canonical pre-run body) | `{reg.sha256}` |",
        f"| {os.path.basename(sites_path)} | `{sha256_file(sites_path)}` |",
    ]
    for p in extra_paths or []:
        lines.append(f"| {os.path.basename(p)} | `{sha256_file(p)}` |")

    lines += ["", "## Criteria locked before the run", "", "**Would count as new**", ""]
    lines += [f"- {c}" for c in reg.novelty_criteria]
    lines += ["", "**Would not count as new, however recovered**", ""]
    lines += [f"- {c}" for c in reg.non_novelty_criteria]
    lines += ["", f"## Pre-run claims ({len(reg.pre_run_claims())}, all verified)", ""]
    lines += ["| id | type | identifier | claim |", "| --- | --- | --- | --- |"]
    for c in reg.pre_run_claims():
        lines.append(f"| {c.claim_id} | {c.claim_type} | {c.identifier or '-'} | {c.claim} |")
    text = "\n".join(lines) + "\n"
    with open(out_path, "w") as fh:
        fh.write(text)
    return text


# ------------------------------------------------------------ sites skeleton
SITES_HEADER = "label\tref_resnum\texpected\trole\tkey\taltered"


def sites_skeleton(reg: Register) -> str:
    """Emit a sites.tsv skeleton from verified residue_evidence claims.

    Proposed, not authoritative: the human accepts or edits. The altered column
    is left empty on purpose and flagged, because declaring one residue where a
    family tolerates two is the failure mode whose signature is a spuriously
    clean loss.
    """
    rows = [
        "# generated by `foldtrace prime sites` -- review before use",
        f"# fold: {reg.fold_id}   reference: {reg.reference_structure}",
        "#",
        "# CHECK BEFORE FREEZING: does any position tolerate a second residue?",
        "# An acid that accepts Asp or Glu, an aromatic that accepts Tyr or Phe,",
        "# must list both in `expected` or the map will report a clean, wrong loss.",
        "#",
        SITES_HEADER,
    ]
    seen: set[int] = set()
    for c in reg.claims:
        if c.claim_type != "residue_evidence" or c.source_verified != "verified":
            continue
        for r in c.residues:
            m = _RESIDUE_RE.match(r)
            if not m or int(m.group(2)) in seen:
                continue
            letter, num = m.group(1), int(m.group(2))
            seen.add(num)
            rows.append(f"{letter}{num}\t{num}\t{letter}\tTODO({c.claim_id})\t1\t")
    if not seen:
        rows.append("# no verified residue_evidence claims yet; run `prime verify` first")
    return "\n".join(rows) + "\n"


# ---------------------------------------------------------------- adjudication
def add_adjudication(reg: Register, feature: str, surfaced_by: str,
                     prior_art_status: str, permitted_claim: str,
                     described: str | None = None, undescribed: str | None = None,
                     supporting: list[str] | None = None) -> Adjudication:
    if reg.frozen_at is None:
        raise PrimeError("cannot adjudicate against an unfrozen register; run `prime freeze` first")
    a = Adjudication(
        feature_id=reg.next_feature_id(), feature=feature, surfaced_by=surfaced_by,
        prior_art_status=prior_art_status, permitted_claim=permitted_claim,
        retrieval_date=_now(), described_component=described,
        undescribed_component=undescribed, supporting_claims=list(supporting or []))
    problems = a.validate()
    for cid in a.supporting_claims:
        reg.claim(cid)  # raises if absent
    if problems:
        raise PrimeError("adjudication rejected:\n  - " + "\n  - ".join(problems))
    reg.adjudications.append(a)
    return a


# ----------------------------------------------------------------- novelty gate
NOVELTY_PATTERNS = [
    re.compile(p, re.I) for p in (
        r"\bnovel\b", r"\bnew (family|architecture|class|fold|system|lineage|enzyme)",
        r"\bpreviously (unknown|unrecognis|unrecogniz|undescribed)",
        r"\bfirst (report|description|example)", r"\bhitherto\b", r"\bunprecedented\b",
        r"\bwe discover", r"\bdiscovery of a\b",
    )
]


_STOP = {"which", "their", "there", "these", "those", "other", "where", "while",
         "among", "records", "members", "within", "already", "search", "level",
         "family", "families", "protein", "proteins", "component", "recovery"}


def _keywords(text: str, limit: int = 6) -> set[str]:
    # allow digits after the first letter so accessions (duf86, upf0331, pf01934)
    # count as content words -- they are usually the most specific term present
    words = [w for w in re.findall(r"[a-z][a-z0-9]{4,}", text.lower()) if w not in _STOP]
    seen: list[str] = []
    for w in words:
        if w not in seen:
            seen.append(w)
    return set(seen[:limit])


def _sentences(text: str) -> list[str]:
    parts = re.split(r"(?<=[.!?])\s+|\n{2,}", text)
    return [p.strip() for p in parts if p.strip()]


def novelty_gate(reg: Register, text: str) -> list[str]:
    """Flag novelty language in draft text that no adjudication supports.

    Scoped to the sentence. A novelty claim about one feature must not be
    flagged against a different, already-described feature elsewhere in the
    draft -- a gate that cries wolf is a gate the author learns to skip.

    Three checks: novelty language with no adjudications at all; novelty
    language in a sentence that names an ``already_described`` feature; and a
    ``partially_described`` feature whose undescribed component the draft never
    names. Wording is left to the author.
    """
    violations: list[str] = []
    if reg.frozen_at is None:
        violations.append("register is not frozen; no novelty claim can be adjudicated")

    sentences = _sentences(text)
    flagged = [s for s in sentences if any(p.search(s) for p in NOVELTY_PATTERNS)]
    if not flagged:
        return violations

    if not reg.adjudications:
        hits = sorted({m.group(0) for p in NOVELTY_PATTERNS for m in p.finditer(text)})
        violations.append(
            f"novelty language present ({', '.join(hits)}) but the register "
            "contains no adjudications")
        return violations

    def _anchor(a: Adjudication) -> set[str]:
        """Terms that count as naming this feature in the draft.

        Built from the feature and its two components, not from
        permitted_claim: that field is prose and its generic words
        ("search", "recovered", "known") match almost any sentence.
        """
        blob = " ".join(filter(None, [a.feature, a.described_component,
                                      a.undescribed_component]))
        return _keywords(blob, limit=25)

    def _names(a: Adjudication, sent: str) -> bool:
        # two shared content words, not one: a single generic overlap is noise,
        # and a gate that cries wolf is a gate the author learns to skip
        return len(_anchor(a) & _keywords(sent, limit=40)) >= 2

    described = [a for a in reg.adjudications if a.prior_art_status == "already_described"]
    for sent in flagged:
        for a in described:
            if _names(a, sent):
                violations.append(
                    f"{a.feature_id} is already_described but novelty language attaches to it "
                    f"here: \"{sent.strip()}\". Permitted: {a.permitted_claim}")
        if not any(_names(a, sent) for a in reg.adjudications):
            violations.append(
                f"novelty language in \"{sent.strip()}\" names no adjudicated feature")

    for a in reg.adjudications:
        if a.prior_art_status == "partially_described" and a.undescribed_component:
            key = _keywords(a.undescribed_component, limit=4)
            if key and not (key & _keywords(text, limit=400)):
                violations.append(
                    f"{a.feature_id} is partially_described but the draft does not name its "
                    f"undescribed component. Permitted: {a.permitted_claim}")
    return violations
