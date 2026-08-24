"""Stage 1 of foldtrace: fold search with Foldseek.

``foldtrace search`` wraps Foldseek's ``easy-search``: it runs a structural search of a
query structure against a Foldseek database (e.g. the AlphaFold/UniProt50 set, ``afdb50``),
parses the tabular output, applies optional thresholds, and writes a ranked machine-readable
table that preserves every Foldseek target identifier the mapping stage needs.

Foldseek itself is an external dependency (https://github.com/steineggerlab/foldseek); it is
invoked as a subprocess. When Foldseek is not installed, an already-produced Foldseek table
can be parsed instead with ``--from-foldseek`` / :func:`parse_foldseek`, which is also how the
offline tests and the bundled end-to-end example run.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass, asdict

# Foldseek easy-search output columns foldtrace requests, in order.
FOLDSEEK_FORMAT_FIELDS = [
    "query", "target", "fident", "alnlen", "qstart", "qend", "tstart", "tend",
    "evalue", "bits", "prob", "qcov", "tcov", "alntmscore",
]
FOLDSEEK_FORMAT = ",".join(FOLDSEEK_FORMAT_FIELDS)

# Output columns of a foldtrace search table.
HIT_COLUMNS = [
    "rank", "query", "target", "prob", "evalue", "bits", "fident",
    "alnlen", "qcov", "tcov", "alntmscore",
]


@dataclass
class Hit:
    query: str
    target: str
    prob: float | None
    evalue: float | None
    bits: float | None
    fident: float | None       # Foldseek sequence identity (fraction 0-1, or percent if >1 in source)
    alnlen: int | None
    qcov: float | None
    tcov: float | None
    alntmscore: float | None
    rank: int = 0


def _num(x):
    x = (x or "").strip()
    if x in ("", "NA", "nan", "None", "*"):
        return None
    try:
        return float(x)
    except ValueError:
        return None


def parse_foldseek(path: str, fields: list[str] | None = None) -> list[Hit]:
    """Parse a Foldseek easy-search table (TSV) into Hit objects.

    ``fields`` is the column order of the file; it defaults to :data:`FOLDSEEK_FORMAT_FIELDS`.
    A one-line ``#``-prefixed header, if present, is ignored. ``fident`` reported as a
    percentage (>1) is normalised to a 0-1 fraction so thresholds are unambiguous.
    """
    fields = fields or FOLDSEEK_FORMAT_FIELDS
    hits: list[Hit] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            row = dict(zip(fields, parts))
            fident = _num(row.get("fident"))
            if fident is not None and fident > 1.0:
                fident = fident / 100.0
            hits.append(Hit(
                query=row.get("query", "").strip(),
                target=row.get("target", "").strip(),
                prob=_num(row.get("prob")),
                evalue=_num(row.get("evalue")),
                bits=_num(row.get("bits")),
                fident=fident,
                alnlen=int(_num(row.get("alnlen"))) if _num(row.get("alnlen")) is not None else None,
                qcov=_num(row.get("qcov")),
                tcov=_num(row.get("tcov")),
                alntmscore=_num(row.get("alntmscore")),
            ))
    return hits


def filter_and_rank(hits: list[Hit], min_prob: float = 0.0, min_coverage: float = 0.0,
                    min_identity: float = 0.0, max_hits: int | None = None) -> list[Hit]:
    """Filter hits by thresholds and rank them.

    Thresholds are permissive by default (0.0) so nothing is silently dropped unless the
    user asks. ``min_coverage`` is applied to query coverage (``qcov``); ``min_identity``
    to Foldseek ``fident`` (0-1 fraction). Hits missing a value for an *active* threshold
    (threshold > 0) are excluded, because they cannot be shown to pass it. Ranking is by
    ``prob`` then ``alntmscore`` then ``bits`` (descending), with Foldseek's own order as
    the final tiebreak.
    """
    kept = []
    for h in hits:
        if min_prob > 0.0 and (h.prob is None or h.prob < min_prob):
            continue
        if min_coverage > 0.0 and (h.qcov is None or h.qcov < min_coverage):
            continue
        if min_identity > 0.0 and (h.fident is None or h.fident < min_identity):
            continue
        kept.append(h)
    order = list(range(len(kept)))
    kept_sorted = sorted(
        zip(order, kept),
        key=lambda ok: (
            -(ok[1].prob if ok[1].prob is not None else -1.0),
            -(ok[1].alntmscore if ok[1].alntmscore is not None else -1.0),
            -(ok[1].bits if ok[1].bits is not None else -1.0),
            ok[0],
        ),
    )
    ranked = [h for _, h in kept_sorted]
    if max_hits is not None:
        ranked = ranked[:max_hits]
    for i, h in enumerate(ranked, start=1):
        h.rank = i
    return ranked


def run_foldseek(query: str, database: str, tmp_dir: str, max_hits: int = 1000,
                 foldseek_bin: str = "foldseek", extra_args: list[str] | None = None) -> str:
    """Run ``foldseek easy-search`` and return the path to the raw result table.

    Requires the Foldseek binary on PATH (or an explicit ``foldseek_bin``). ``--format-output``
    is fixed to :data:`FOLDSEEK_FORMAT` so parsing is deterministic; ``--max-seqs`` is set to
    ``max_hits``. Raises ``FileNotFoundError`` if Foldseek is not available and
    ``subprocess.CalledProcessError`` if it fails.
    """
    if shutil.which(foldseek_bin) is None and not os.path.exists(foldseek_bin):
        raise FileNotFoundError(
            f"Foldseek binary '{foldseek_bin}' not found. Install Foldseek "
            "(https://github.com/steineggerlab/foldseek), or supply an existing Foldseek "
            "table with --from-foldseek.")
    os.makedirs(tmp_dir, exist_ok=True)
    raw = os.path.join(tmp_dir, "foldseek_raw.tsv")
    cmd = [
        foldseek_bin, "easy-search", query, database, raw, os.path.join(tmp_dir, "fs_tmp"),
        "--format-output", FOLDSEEK_FORMAT, "--max-seqs", str(max_hits),
    ]
    if extra_args:
        cmd += list(extra_args)
    subprocess.run(cmd, check=True)
    return raw


def write_hits_tsv(hits: list[Hit], path, delimiter: str = "\t") -> None:
    """Write ranked hits as a foldtrace search table (columns :data:`HIT_COLUMNS`)."""
    def fmt(v):
        if v is None:
            return "NA"
        if isinstance(v, float):
            return f"{v:.4g}"
        return str(v)
    close = False
    if hasattr(path, "write"):
        fh = path
    else:
        fh = open(path, "w")
        close = True
    try:
        fh.write(delimiter.join(HIT_COLUMNS) + "\n")
        for h in hits:
            d = asdict(h)
            fh.write(delimiter.join(fmt(d[c]) for c in HIT_COLUMNS) + "\n")
    finally:
        if close:
            fh.close()


def search(query: str, database: str | None, out, tmp_dir: str = "foldtrace_tmp",
           max_hits: int = 1000, min_prob: float = 0.0, min_coverage: float = 0.0,
           min_identity: float = 0.0, foldseek_bin: str = "foldseek",
           from_foldseek: str | None = None, delimiter: str = "\t") -> list[Hit]:
    """End-to-end stage 1: run (or parse) Foldseek, filter, rank, and write the hit table."""
    if from_foldseek is not None:
        raw = from_foldseek
    else:
        if not database:
            raise ValueError("a Foldseek --database is required unless --from-foldseek is given")
        raw = run_foldseek(query, database, tmp_dir, max_hits=max_hits, foldseek_bin=foldseek_bin)
    hits = parse_foldseek(raw)
    ranked = filter_and_rank(hits, min_prob=min_prob, min_coverage=min_coverage,
                             min_identity=min_identity, max_hits=max_hits)
    if out is not None:
        write_hits_tsv(ranked, out, delimiter=delimiter)
    return ranked
