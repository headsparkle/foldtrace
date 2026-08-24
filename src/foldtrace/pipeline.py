"""Top-level ``foldtrace run``: stage 1 (search) followed by stage 2 (active-site mapping).

``run`` obtains a ranked set of candidate hits (by running Foldseek, by reading a Foldseek
table, by reading a foldtrace ``search`` table, or from explicit candidate files), resolves a
structure for each hit, maps it with :func:`foldtrace.mapping.map_candidate`, and writes one
merged record per candidate that carries the Foldseek search metrics, the TM-align score, the
per-site residue/offset/state, the overall verdict, and an explicit reason for any unresolved
call. Candidates whose structure cannot be found are reported as unresolved, never dropped.
"""
from __future__ import annotations

import glob
import os
import urllib.request
from dataclasses import dataclass

from .io import load_structure, parse_sites, Site
from .mapping import map_candidate, CandidateResult, UNRESOLVED
from .search import Hit, HIT_COLUMNS, parse_foldseek, filter_and_rank, search as run_search

_STRUCTURE_EXTS = (".pdb", ".cif", ".mmcif")


def read_search_table(path: str) -> list[Hit]:
    """Read a foldtrace ``search`` table (header = :data:`foldtrace.search.HIT_COLUMNS`)."""
    hits: list[Hit] = []
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = [p.strip() for p in parts]
                continue
            row = dict(zip(header, parts))

            def num(k):
                v = (row.get(k) or "").strip()
                if v in ("", "NA", "None"):
                    return None
                try:
                    return float(v)
                except ValueError:
                    return None
            hits.append(Hit(
                query=row.get("query", "").strip(), target=row.get("target", "").strip(),
                prob=num("prob"), evalue=num("evalue"), bits=num("bits"), fident=num("fident"),
                alnlen=int(num("alnlen")) if num("alnlen") is not None else None,
                qcov=num("qcov"), tcov=num("tcov"), alntmscore=num("alntmscore"),
                rank=int(num("rank")) if num("rank") is not None else 0,
            ))
    return hits


def resolve_structure_path(target: str, candidates_dir: str | None,
                           fetch_dir: str | None = None, fetch: bool = False) -> str | None:
    """Find a structure file for a hit target: look in ``candidates_dir`` (``{target}.pdb`` etc.),
    then optionally download the AlphaFold model into ``fetch_dir``. Returns a path or None."""
    if candidates_dir:
        for ext in _STRUCTURE_EXTS:
            p = os.path.join(candidates_dir, target + ext)
            if os.path.exists(p):
                return p
    if fetch and fetch_dir:
        os.makedirs(fetch_dir, exist_ok=True)
        dest = os.path.join(fetch_dir, target + ".pdb")
        if os.path.exists(dest):
            return dest
        try:
            api = f"https://alphafold.ebi.ac.uk/api/prediction/{target}"
            import json
            meta = json.loads(urllib.request.urlopen(api, timeout=30).read().decode())
            url = meta[0]["pdbUrl"]
            urllib.request.urlretrieve(url, dest)
            return dest
        except Exception:
            return None
    return None


@dataclass
class RunRecord:
    hit: Hit | None
    result: CandidateResult | None
    reason: str = ""   # populated when there is no mapping result (e.g. structure not found)


def run(reference: str, sites_path: str, *,
        database: str | None = None, query: str | None = None,
        from_foldseek: str | None = None, hits_table: str | None = None,
        candidates: list[str] | None = None, candidates_dir: str | None = None,
        fetch: bool = False, fetch_dir: str | None = None,
        max_hits: int = 1000, min_prob: float = 0.0, min_coverage: float = 0.0,
        min_identity: float = 0.0, foldseek_bin: str = "foldseek",
        offset_threshold: float = 4.0, tm_gate: float = 0.5,
        ref_chain: str | None = None, tmp_dir: str = "foldtrace_tmp") -> tuple[list[Site], list[RunRecord]]:
    """Execute the two-stage workflow and return (sites, per-candidate RunRecords)."""
    ref = load_structure(reference, chain=ref_chain)
    sites = parse_sites(sites_path)

    # ---- Stage 1: obtain ranked hits
    if candidates:
        expanded: list[str] = []
        for c in candidates:
            g = sorted(glob.glob(c))
            expanded.extend(g if g else [c])
        hits = [Hit(query=query or os.path.basename(reference),
                    target=os.path.splitext(os.path.basename(p))[0],
                    prob=None, evalue=None, bits=None, fident=None, alnlen=None,
                    qcov=None, tcov=None, alntmscore=None, rank=i)
                for i, p in enumerate(expanded, start=1)]
        direct_paths = {os.path.splitext(os.path.basename(p))[0]: p for p in expanded}
    elif hits_table:
        hits = filter_and_rank(read_search_table(hits_table), min_prob, min_coverage, min_identity, max_hits)
        direct_paths = {}
    elif from_foldseek:
        hits = filter_and_rank(parse_foldseek(from_foldseek), min_prob, min_coverage, min_identity, max_hits)
        direct_paths = {}
    elif database:
        hits = run_search(query or reference, database, out=None, tmp_dir=tmp_dir, max_hits=max_hits,
                          min_prob=min_prob, min_coverage=min_coverage, min_identity=min_identity,
                          foldseek_bin=foldseek_bin)
        direct_paths = {}
    else:
        raise ValueError("run needs one of: --database, --from-foldseek, --hits, or --candidates")

    # ---- Stage 2: resolve structure + map each hit
    records: list[RunRecord] = []
    for h in hits:
        path = direct_paths.get(h.target) or resolve_structure_path(
            h.target, candidates_dir, fetch_dir=fetch_dir, fetch=fetch)
        if path is None:
            records.append(RunRecord(hit=h, result=None, reason="structure_not_found"))
            continue
        try:
            cand = load_structure(path)
        except Exception as exc:
            records.append(RunRecord(hit=h, result=None, reason=f"load_error: {exc}"))
            continue
        result = map_candidate(ref, cand, sites, candidate_name=h.target,
                               offset_threshold=offset_threshold, tm_gate=tm_gate)
        records.append(RunRecord(hit=h, result=result))
    return sites, records


def result_columns(sites: list[Site]) -> list[str]:
    cols = ["candidate_id", "rank", "foldseek_prob", "foldseek_evalue", "foldseek_bits",
            "foldseek_fident", "foldseek_qcov", "foldseek_alntmscore",
            "tmalign_tm_norm_ref", "tmalign_rmsd", "fold_ok"]
    for s in sites:
        cols += [f"{s.label}_obs", f"{s.label}_offsetA", f"{s.label}_state"]
    cols += ["verdict", "state_reason"]
    return cols


def write_results(sites: list[Site], records: list[RunRecord], out, delimiter: str = "\t") -> None:
    cols = result_columns(sites)

    def fmt(v):
        if v is None:
            return "NA"
        if isinstance(v, float):
            return f"{v:.4g}"
        return str(v)

    close = False
    if hasattr(out, "write"):
        fh = out
    else:
        fh = open(out, "w")
        close = True
    try:
        fh.write(delimiter.join(cols) + "\n")
        for rec in records:
            h, r = rec.hit, rec.result
            row = {c: "NA" for c in cols}
            row["candidate_id"] = h.target if h else (r.name if r else "?")
            if h:
                row.update({"rank": h.rank, "foldseek_prob": fmt(h.prob),
                            "foldseek_evalue": fmt(h.evalue), "foldseek_bits": fmt(h.bits),
                            "foldseek_fident": fmt(h.fident), "foldseek_qcov": fmt(h.qcov),
                            "foldseek_alntmscore": fmt(h.alntmscore)})
            if r is None:
                row["verdict"] = UNRESOLVED
                row["state_reason"] = rec.reason
            else:
                row.update({"tmalign_tm_norm_ref": fmt(r.tm_norm_ref), "tmalign_rmsd": fmt(r.rmsd),
                            "fold_ok": str(r.fold_ok).lower(), "verdict": r.verdict,
                            "state_reason": r.state_reason})
                for c in r.sites:
                    obs = f"{c.obs_res}{c.obs_resnum}" if c.obs_resnum is not None else "-"
                    row[f"{c.label}_obs"] = obs
                    row[f"{c.label}_offsetA"] = fmt(c.ca_offset)
                    row[f"{c.label}_state"] = c.state
            fh.write(delimiter.join(fmt(row[c]) for c in cols) + "\n")
    finally:
        if close:
            fh.close()
