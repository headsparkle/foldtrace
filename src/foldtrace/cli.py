"""Command-line interface for foldtrace."""
from __future__ import annotations

import argparse
import glob
import os
import sys

from . import __version__
from .io import load_structure, parse_sites
from .mapping import map_candidate


def _expand(paths: list[str]) -> list[str]:
    out: list[str] = []
    for p in paths:
        hits = sorted(glob.glob(p))
        out.extend(hits if hits else [p])
    return out


def cmd_map(args: argparse.Namespace) -> int:
    ref = load_structure(args.reference, chain=args.ref_chain)
    sites = parse_sites(args.sites)
    candidates = _expand(args.candidates)
    if not candidates:
        print("error: no candidate structures given", file=sys.stderr)
        return 2

    header = ["candidate", "tm_norm_ref", "rmsd", "fold_ok"]
    for s in sites:
        header += [f"{s.label}_obs", f"{s.label}_offsetA", f"{s.label}_state"]
    header += ["verdict"]

    rows = []
    for path in candidates:
        name = os.path.splitext(os.path.basename(path))[0]
        try:
            cand = load_structure(path, chain=args.cand_chain)
        except Exception as exc:  # keep going through a batch
            print(f"warning: skipping {name}: {exc}", file=sys.stderr)
            continue
        r = map_candidate(ref, cand, sites, candidate_name=name,
                          offset_threshold=args.offset_threshold, tm_gate=args.tm_gate)
        row = [r.name, f"{r.tm_norm_ref:.3f}", f"{r.rmsd:.2f}", str(r.fold_ok).lower()]
        for c in r.sites:
            obs = f"{c.obs_res}{c.obs_resnum}" if c.obs_resnum is not None else "-"
            off = f"{c.ca_offset:.2f}" if c.ca_offset is not None else "-"
            row += [obs, off, c.state]
        row += [r.verdict]
        rows.append(row)

    out = open(args.out, "w") if args.out else sys.stdout
    try:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")
    finally:
        if args.out:
            out.close()
    if args.out:
        n = {"retained": 0, "lost": 0, "unresolved": 0}
        for row in rows:
            n[row[-1]] = n.get(row[-1], 0) + 1
        print(f"wrote {args.out}: {len(rows)} candidates "
              f"({n.get('retained',0)} retained, {n.get('lost',0)} lost, {n.get('unresolved',0)} unresolved)")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="foldtrace",
        description="Structure-first active-site mapping: read the catalytic state of a "
                    "fold from predicted structures.")
    p.add_argument("--version", action="version", version=f"foldtrace {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    m = sub.add_parser("map", help="read active-site state of candidate structures against a reference")
    m.add_argument("--reference", required=True, help="reference structure (PDB/mmCIF)")
    m.add_argument("--sites", required=True, help="catalytic-sites TSV (see docs)")
    m.add_argument("--candidates", required=True, nargs="+",
                   help="candidate structures (files or globs)")
    m.add_argument("--out", help="output TSV (default: stdout)")
    m.add_argument("--ref-chain", default=None, help="reference chain id (default: first)")
    m.add_argument("--cand-chain", default=None, help="candidate chain id (default: first)")
    m.add_argument("--offset-threshold", type=float, default=4.0,
                   help="max CA-CA offset in Angstrom for a trusted call (default: 4.0)")
    m.add_argument("--tm-gate", type=float, default=0.5,
                   help="minimum reference-normalised TM-score for same-fold (default: 0.5)")
    m.set_defaults(func=cmd_map)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
