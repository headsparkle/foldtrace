"""Command-line interface for foldtrace (search, map, run)."""
from __future__ import annotations

import argparse
import glob
import os
import sys

from . import __version__
from .io import load_structure, parse_sites
from .mapping import map_candidate
from . import search as search_mod
from . import pipeline


def _expand(paths: list[str]) -> list[str]:
    out: list[str] = []
    for p in paths:
        hits = sorted(glob.glob(p))
        out.extend(hits if hits else [p])
    return out


def _delim(fmt: str) -> str:
    return "," if fmt == "csv" else "\t"


# --------------------------------------------------------------------- search
def cmd_search(args: argparse.Namespace) -> int:
    out = open(args.out, "w") if args.out else sys.stdout
    try:
        hits = search_mod.search(
            query=args.query, database=args.database, out=out, tmp_dir=args.tmp_dir,
            max_hits=args.max_hits, min_prob=args.min_prob, min_coverage=args.min_coverage,
            min_identity=args.min_identity, foldseek_bin=args.foldseek,
            from_foldseek=args.from_foldseek, delimiter=_delim(args.format))
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    finally:
        if args.out:
            out.close()
    if args.out:
        print(f"wrote {args.out}: {len(hits)} hits")
    return 0


# ------------------------------------------------------------------------ map
def cmd_map(args: argparse.Namespace) -> int:
    ref = load_structure(args.reference, chain=args.ref_chain)
    sites = parse_sites(args.sites)
    candidates = _expand(args.candidates)
    if not candidates:
        print("error: no candidate structures given", file=sys.stderr)
        return 2

    delim = _delim(args.format)
    header = ["candidate", "tm_norm_ref", "rmsd", "fold_ok"]
    for s in sites:
        header += [f"{s.label}_obs", f"{s.label}_offsetA", f"{s.label}_state"]
    header += ["verdict", "state_reason"]

    rows = []
    for path in candidates:
        name = os.path.splitext(os.path.basename(path))[0]
        try:
            cand = load_structure(path, chain=args.cand_chain)
        except Exception as exc:
            print(f"warning: skipping {name}: {exc}", file=sys.stderr)
            continue
        r = map_candidate(ref, cand, sites, candidate_name=name,
                          offset_threshold=args.offset_threshold, tm_gate=args.tm_gate)
        row = [r.name, f"{r.tm_norm_ref:.3f}", f"{r.rmsd:.2f}", str(r.fold_ok).lower()]
        for c in r.sites:
            obs = f"{c.obs_res}{c.obs_resnum}" if c.obs_resnum is not None else "-"
            off = f"{c.ca_offset:.2f}" if c.ca_offset is not None else "NA"
            row += [obs, off, c.state]
        row += [r.verdict, r.state_reason]
        rows.append(row)

    out = open(args.out, "w") if args.out else sys.stdout
    try:
        out.write(delim.join(header) + "\n")
        for row in rows:
            out.write(delim.join(row) + "\n")
    finally:
        if args.out:
            out.close()
    if args.out:
        n = {}
        for row in rows:
            n[row[-2]] = n.get(row[-2], 0) + 1
        summary = ", ".join(f"{v} {k}" for k, v in sorted(n.items()))
        print(f"wrote {args.out}: {len(rows)} candidates ({summary})")
    return 0


# ------------------------------------------------------------------------ run
def cmd_run(args: argparse.Namespace) -> int:
    try:
        sites, records = pipeline.run(
            reference=args.reference, sites_path=args.sites,
            database=args.database, query=args.query, from_foldseek=args.from_foldseek,
            hits_table=args.hits, candidates=_expand(args.candidates) if args.candidates else None,
            candidates_dir=args.candidates_dir, fetch=args.fetch, fetch_dir=args.fetch_dir,
            max_hits=args.max_hits, min_prob=args.min_prob, min_coverage=args.min_coverage,
            min_identity=args.min_identity, foldseek_bin=args.foldseek,
            offset_threshold=args.offset_threshold, tm_gate=args.tm_gate,
            ref_chain=args.ref_chain, tmp_dir=args.tmp_dir)
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    out = open(args.out, "w") if args.out else sys.stdout
    try:
        pipeline.write_results(sites, records, out, delimiter=_delim(args.format))
    finally:
        if args.out:
            out.close()
    if args.out:
        n = {}
        for rec in records:
            v = rec.result.verdict if rec.result else "unresolved"
            n[v] = n.get(v, 0) + 1
        summary = ", ".join(f"{v} {k}" for k, v in sorted(n.items()))
        print(f"wrote {args.out}: {len(records)} candidates ({summary})")
    return 0


# --------------------------------------------------------------------- parser
def _add_threshold_opts(p):
    p.add_argument("--max-hits", type=int, default=1000, help="max hits to keep (default: 1000)")
    p.add_argument("--min-prob", type=float, default=0.0, help="min Foldseek match probability (default: 0.0 = off)")
    p.add_argument("--min-coverage", type=float, default=0.0, help="min query coverage 0-1 (default: 0.0 = off)")
    p.add_argument("--min-identity", type=float, default=0.0, help="min Foldseek identity 0-1 (default: 0.0 = off)")


def _add_map_opts(p):
    p.add_argument("--offset-threshold", type=float, default=4.0,
                   help="max CA-CA offset (A) for a trusted call; beyond it a site is unresolved (default: 4.0)")
    p.add_argument("--tm-gate", type=float, default=0.5,
                   help="min reference-normalised TM-score for same-fold (default: 0.5)")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="foldtrace",
        description="Structure-first search (stage 1) and active-site mapping (stage 2): "
                    "find remote homologs of a fold and read their catalytic state from predicted structures.")
    p.add_argument("--version", action="version", version=f"foldtrace {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    # search
    s = sub.add_parser("search", help="stage 1: fold search with Foldseek -> ranked hit table")
    s.add_argument("--query", help="query structure (PDB/mmCIF) for Foldseek")
    s.add_argument("--database", help="Foldseek database (e.g. afdb50 for AlphaFold/UniProt50)")
    s.add_argument("--from-foldseek", help="parse an existing Foldseek table instead of running Foldseek")
    s.add_argument("--out", help="output hit table (default: stdout)")
    s.add_argument("--format", choices=["tsv", "csv"], default="tsv", help="output delimiter (default: tsv)")
    s.add_argument("--foldseek", default="foldseek", help="Foldseek binary (default: foldseek on PATH)")
    s.add_argument("--tmp-dir", default="foldtrace_tmp", help="scratch dir for Foldseek (default: foldtrace_tmp)")
    _add_threshold_opts(s)
    s.set_defaults(func=cmd_search)

    # map
    m = sub.add_parser("map", help="stage 2: read active-site state of candidate structures")
    m.add_argument("--reference", required=True, help="reference structure (PDB/mmCIF)")
    m.add_argument("--sites", required=True, help="catalytic-sites TSV")
    m.add_argument("--candidates", required=True, nargs="+", help="candidate structures (files or globs)")
    m.add_argument("--out", help="output TSV/CSV (default: stdout)")
    m.add_argument("--format", choices=["tsv", "csv"], default="tsv", help="output delimiter (default: tsv)")
    m.add_argument("--ref-chain", default=None, help="reference chain id (default: first)")
    m.add_argument("--cand-chain", default=None, help="candidate chain id (default: first)")
    _add_map_opts(m)
    m.set_defaults(func=cmd_map)

    # run
    r = sub.add_parser("run", help="stages 1+2: search then map (end-to-end)")
    r.add_argument("--reference", required=True, help="reference structure with known catalytic residues")
    r.add_argument("--sites", required=True, help="catalytic-sites TSV")
    r.add_argument("--query", help="query structure for the search (defaults to --reference)")
    src = r.add_argument_group("stage-1 source (choose one)")
    src.add_argument("--database", help="Foldseek database; runs the search")
    src.add_argument("--from-foldseek", help="existing Foldseek table")
    src.add_argument("--hits", help="a foldtrace search table")
    src.add_argument("--candidates", nargs="+", help="explicit candidate structures (skip search)")
    r.add_argument("--candidates-dir", help="directory of hit structures named {target}.pdb")
    r.add_argument("--fetch", action="store_true", help="download missing AlphaFold models for hit targets")
    r.add_argument("--fetch-dir", default="foldtrace_models", help="where --fetch saves models")
    r.add_argument("--out", help="output TSV/CSV (default: stdout)")
    r.add_argument("--format", choices=["tsv", "csv"], default="tsv", help="output delimiter (default: tsv)")
    r.add_argument("--ref-chain", default=None, help="reference chain id (default: first)")
    r.add_argument("--foldseek", default="foldseek", help="Foldseek binary (default: foldseek on PATH)")
    r.add_argument("--tmp-dir", default="foldtrace_tmp", help="scratch dir for Foldseek")
    _add_threshold_opts(r)
    _add_map_opts(r)
    r.set_defaults(func=cmd_run)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
