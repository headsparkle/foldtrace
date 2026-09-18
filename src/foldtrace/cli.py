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
        if not args.no_advice:
            from . import advise as advise_mod
            directory = os.path.dirname(os.path.abspath(args.out)) or "."
            text = advise_mod.next_steps(directory, brief=True)
            if text.strip():
                print("\n--- not looked at yet " + "-" * 44)
                print(text, end="")
    return 0


def cmd_next(args: argparse.Namespace) -> int:
    from . import advise as advise_mod
    sys.stdout.write(advise_mod.next_steps(args.dir, brief=args.brief))
    return 0


# --------------------------------------------------------------------- guided
def _guided_load(path: str):
    from .guided import GuidedProject
    if not os.path.exists(path):
        raise FileNotFoundError(f"no guided project at {path} (run 'foldtrace guided init' first)")
    return GuidedProject.load(path)


def cmd_guided(args: argparse.Namespace) -> int:
    from .guided import GuidedProject, GuidedError

    def _dump(obj):
        import json
        print(json.dumps(obj, indent=2))

    try:
        if args.guided_cmd == "init":
            proj = GuidedProject(args.name, args.reference, args.sites, args.candidate,
                                 tm_gate=args.tm_gate, offset_threshold=args.offset_threshold)
            proj.save(args.project)
            print(f"created guided project '{args.name}' at {args.project}; next: observe")
            return 0

        proj = _guided_load(args.project)

        if args.guided_cmd == "observe":
            _dump(proj.observe(args.notes or ""))
        elif args.guided_cmd == "predict":
            preds = {}
            for item in args.set or []:
                if "=" not in item:
                    print(f"error: --set expects label=state, got {item!r}", file=sys.stderr)
                    return 2
                lab, state = item.split("=", 1)
                preds[lab.strip()] = state.strip()
            _dump(proj.predict(preds, rationale=args.rationale or "", verdict=args.verdict))
        elif args.guided_cmd == "compute":
            r = proj.compute()
            print(f"verdict: {r.verdict} (TM {r.tm_norm_ref:.2f}, RMSD {r.rmsd:.2f} A) -- prediction locked")
        elif args.guided_cmd == "compare":
            _dump(proj.compare())
        elif args.guided_cmd == "decide":
            _dump(proj.decide(args.conclusion, next_action=args.next or ""))
        elif args.guided_cmd == "status":
            _dump(proj.status())
        elif args.guided_cmd == "report":
            out = open(args.out, "w") if args.out else sys.stdout
            try:
                out.write(proj.report())
            finally:
                if args.out:
                    out.close()
            if args.out:
                print(f"wrote {args.out}")
            return 0
        else:  # pragma: no cover - argparse guards this
            print("error: unknown guided subcommand", file=sys.stderr)
            return 2
    except (GuidedError, FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.guided_cmd not in ("report",):
        proj.save(args.project)
    return 0



# ---------------------------------------------------------------------- prime
def cmd_prime(args: argparse.Namespace) -> int:
    import json
    from . import prime as prime_mod
    from .prime import PrimeError

    try:
        if args.prime_cmd == "init":
            reg = prime_mod.new_register(args.fold, args.reference,
                                         scope=args.scope or "",
                                         tools=args.tool or [])
            if os.path.exists(args.register) and not args.force:
                print(f"error: {args.register} exists (use --force)", file=sys.stderr)
                return 2
            prime_mod.save(reg, args.register)
            print(f"created register {args.register} for fold '{args.fold}'; "
                  "next: add claims, then 'prime verify'")
            return 0

        reg = prime_mod.load(args.register)

        if args.prime_cmd == "add-claim":
            c = prime_mod.Claim(
                claim_id=args.id or reg.next_claim_id(), claim_type=args.type,
                claim=args.claim, identifier=args.identifier,
                residues=args.residue or [], evidence_strength=args.evidence,
                discovered_post_run=args.post_run, trigger=args.trigger,
                note=args.note or "")
            problems = c.validate()
            if problems:
                print("error: " + "; ".join(problems), file=sys.stderr)
                return 2
            reg.claims.append(c)
            print(f"added {c.claim_id} (unverified)")

        elif args.prime_cmd == "criteria":
            reg.novelty_criteria += args.novel or []
            reg.non_novelty_criteria += args.not_novel or []
            print(f"novelty: {len(reg.novelty_criteria)}, "
                  f"non-novelty: {len(reg.non_novelty_criteria)}")

        elif args.prime_cmd == "verify":
            with open(args.evidence) as fh:
                records = [json.loads(l) for l in fh if l.strip()]
            ok, failed = prime_mod.apply_verification(reg, records)
            print(f"verified {ok}, failed {failed}")
            if args.fail_on_unresolved:
                left = [c.claim_id for c in reg.pre_run_claims()
                        if c.source_verified != "verified"]
                if left:
                    prime_mod.save(reg, args.register)
                    print("error: unverified pre-run claims: " + ", ".join(left),
                          file=sys.stderr)
                    return 2

        elif args.prime_cmd == "validate":
            problems = prime_mod.validate(reg, sites_path=args.sites)
            if problems:
                for p_ in problems:
                    print(f"- {p_}")
                return 1
            print("register is fit to freeze")
            return 0

        elif args.prime_cmd == "sites":
            text = prime_mod.sites_skeleton(reg)
            if args.out:
                with open(args.out, "w") as fh:
                    fh.write(text)
                print(f"wrote {args.out} -- review the altered column before freezing")
            else:
                sys.stdout.write(text)
            return 0

        elif args.prime_cmd == "freeze":
            prime_mod.freeze(reg, args.sites, extra_paths=args.also or [],
                             out_path=args.out)
            prime_mod.save(reg, args.register)
            print(f"locked {reg.frozen_at}; wrote {args.out}")
            return 0

        elif args.prime_cmd == "adjudicate":
            a = prime_mod.add_adjudication(
                reg, feature=args.feature, surfaced_by=args.surfaced_by,
                prior_art_status=args.status, permitted_claim=args.permitted,
                described=args.described, undescribed=args.undescribed,
                supporting=args.supports or [])
            print(f"added {a.feature_id} ({a.prior_art_status})")

        elif args.prime_cmd == "gate":
            with open(args.draft) as fh:
                text = fh.read()
            violations = prime_mod.novelty_gate(reg, text)
            if violations:
                for v in violations:
                    print(f"- {v}")
                return 1
            print("no unadjudicated novelty language found")
            return 0

        elif args.prime_cmd == "export":
            cols = ["claim_id", "claim_type", "identifier", "source_verified",
                    "evidence_strength", "discovered_post_run", "residues", "claim"]
            lines = ["\t".join(cols)]
            for c in reg.claims:
                lines.append("\t".join([
                    c.claim_id, c.claim_type, c.identifier or "-", c.source_verified,
                    c.evidence_strength, str(c.discovered_post_run).lower(),
                    ",".join(c.residues), c.claim]))
            text = "\n".join(lines) + "\n"
            if args.out:
                with open(args.out, "w") as fh:
                    fh.write(text)
                print(f"wrote {args.out}: {len(reg.claims)} claims")
            else:
                sys.stdout.write(text)
            return 0

        else:  # pragma: no cover - argparse guards this
            print("error: unknown prime subcommand", file=sys.stderr)
            return 2

    except (PrimeError, FileNotFoundError, ValueError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    prime_mod.save(reg, args.register)
    return 0



# ------------------------------------------------------------------- discord
def cmd_discord(args: argparse.Namespace) -> int:
    from . import discord as dis
    from .discord import DiscordError

    try:
        calls = dis.load_calls(args.calls)
        findings = []

        if args.annotations:
            if not args.vocab:
                print("error: --annotations needs --vocab", file=sys.stderr)
                return 2
            findings += dis.scan_annotation(calls, dis.load_annotations(args.annotations),
                                            dis.load_vocab(args.vocab))
        if args.compare:
            other = dis.load_calls(args.compare)
            findings += dis.scan_paired(calls, other, args.label_a, args.label_b)
        if args.profile:
            from . import invariant as inv
            findings += dis.scan_identity(inv.read_profile(args.profile),
                                          dis.load_annotations(args.annotations)
                                          if args.annotations else {},
                                          min_identity=args.min_identity,
                                          organism_key=args.organism_key)
        if not args.no_knockout:
            findings += dis.scan_knockout(calls, min_tm=args.knockout_tm,
                                          max_prevalence=args.knockout_max_prevalence)

        findings = dis.rank(findings)
        if args.min_score:
            findings = [f for f in findings if f.score >= args.min_score]
        if args.kind:
            keep = set(args.kind)
            findings = [f for f in findings if f.kind in keep]
    except (DiscordError, FileNotFoundError, ValueError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    out = open(args.out, "w") if args.out else sys.stdout
    try:
        dis.write(findings, out, delimiter=_delim(args.format))
    finally:
        if args.out:
            out.close()
    if args.out:
        print(f"wrote {args.out}: {len(findings)} findings")
        for line in dis.summarise(findings):
            print(line)
    if args.emit_features:
        print("\n# suggested adjudications (every feature needs one before the gate passes)")
        for line in dis.emit_features(findings):
            print(line)
    return 0



# ----------------------------------------------------------------- invariant
def cmd_profile(args: argparse.Namespace) -> int:
    from . import invariant as inv

    ref = load_structure(args.reference, chain=args.ref_chain)
    candidates = _expand(args.candidates)
    if not candidates:
        print("error: no candidate structures given", file=sys.stderr)
        return 2

    rows, skipped, below = [], 0, 0
    for path in candidates:
        name = os.path.splitext(os.path.basename(path))[0]
        try:
            cand = load_structure(path, chain=args.cand_chain)
        except Exception as exc:
            print(f"warning: skipping {name}: {exc}", file=sys.stderr)
            skipped += 1
            continue
        obs, tm, fold_ok = inv.profile_candidate(
            ref, cand, name, offset_threshold=args.offset_threshold, tm_gate=args.tm_gate)
        below += not fold_ok
        rows.extend(obs)

    out = open(args.out, "w") if args.out else sys.stdout
    try:
        inv.write_profile(rows, out, delimiter=_delim(args.format))
    finally:
        if args.out:
            out.close()
    if args.out:
        n = len(candidates) - skipped
        print(f"wrote {args.out}: {len(rows)} rows, {n} candidates "
              f"({below} below the TM gate, {skipped} unreadable)")
    return 0


def cmd_invariant(args: argparse.Namespace) -> int:
    from . import invariant as inv
    from .invariant import InvariantError

    try:
        ref = load_structure(args.reference, chain=args.ref_chain)
        sites = parse_sites(args.sites)
        profile = inv.read_profile(args.profile)
        groups = inv.load_groups(args.groups)
        clusters = inv.load_clusters(args.clusters) if args.clusters else None
        stats, summary = inv.analyse(
            profile, ref, sites, groups, focus=args.focus, compare=args.compare,
            clusters=clusters, near=args.near, background=args.background,
            min_coverage=args.min_coverage, min_conservation=args.min_conservation)
    except (InvariantError, FileNotFoundError, ValueError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    findings = inv.rank(stats) if not args.all else sorted(
        stats, key=lambda s: (-s.score, s.ref_resnum))

    out = open(args.out, "w") if args.out else sys.stdout
    try:
        inv.write_findings(findings, out, delimiter=_delim(args.format))
    finally:
        if args.out:
            out.close()
    if args.out:
        print(f"wrote {args.out}: {len(findings)} positions")
        for k, v in summary.items():
            print(f"  {k}: {v}")
        if summary["background_positions"] < 10:
            print("  warning: too few background positions to calibrate; "
                  "percentiles are not meaningful", file=sys.stderr)
        frac = summary.get("site_fraction_reported", "NA")
        if frac != "NA" and float(frac) > 0.5:
            print(f"  warning: {frac} of the tested site zone came back invariant; "
                  "the hit set is probably too close or too small to discriminate",
                  file=sys.stderr)
    if args.emit_features:
        print("\n# suggested adjudication (needed before the gate passes)")
        for line in inv.emit_features(findings):
            print(line)
    return 0



# ---------------------------------------------------------------------- rare
def cmd_rare(args: argparse.Namespace) -> int:
    from . import invariant as inv
    from . import rare as rare_mod
    from .invariant import InvariantError

    try:
        ref = load_structure(args.reference, chain=args.ref_chain)
        sites = parse_sites(args.sites)
        profile = inv.read_profile(args.profile)
        clusters = inv.load_clusters(args.clusters) if args.clusters else None
        findings, summary = rare_mod.analyse_rare(
            profile, ref, sites, near=args.near, background=args.background,
            residues=frozenset(args.residues.upper()), min_count=args.min_count,
            max_fraction=args.max_fraction, min_modal=args.min_modal, clusters=clusters,
            include_declared=args.include_declared)
        combos = rare_mod.combinations(findings, min_positions=args.min_positions)
    except (InvariantError, FileNotFoundError, ValueError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    out = open(args.out, "w") if args.out else sys.stdout
    try:
        rare_mod.write_rare(findings, out, delimiter=_delim(args.format))
    finally:
        if args.out:
            out.close()
    if args.out:
        print(f"wrote {args.out}: {len(findings)} rare calls")
        for k, v in summary.items():
            print(f"  {k}: {v}")
    if args.combos:
        rare_mod.write_combos(combos, args.combos, delimiter=_delim(args.format))
        print(f"wrote {args.combos}: {len(combos)} combinations")
    if args.emit_features:
        print("\n# suggested adjudications (needed before the gate passes)")
        for line in rare_mod.emit_features(findings, combos):
            print(line)
    return 0


# ------------------------------------------------------------------ cooccur
def cmd_cooccur(args: argparse.Namespace) -> int:
    from . import cooccur as co
    from .cooccur import CooccurError

    try:
        features = co.load_features(args.features)
        focus = _read_target_set(args.focus_targets)
        bg = _read_target_set(args.background_targets) if args.background_targets else None
        spread = co.load_column(args.features, args.spread_by) if args.spread_by else None
        findings, summary = co.analyse(
            features, focus, background=bg, spread_by=spread,
            min_count=args.min_count, min_odds=args.min_odds, max_q=args.max_q,
            alternative=args.alternative, held_out=args.held_out)
    except (CooccurError, FileNotFoundError, ValueError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    out = open(args.out, "w") if args.out else sys.stdout
    try:
        co.write(findings, out, delimiter=_delim(args.format))
    finally:
        if args.out:
            out.close()
    if args.out:
        print(f"wrote {args.out}: {len(findings)} features")
        for k, v in summary.items():
            print(f"  {k}: {v}")
    if args.emit_features:
        print("\n# suggested adjudications (needed before the gate passes)")
        for line in co.emit_features(findings):
            print(line)
    return 0


def _read_target_set(path: str) -> set[str]:
    """Targets from a one-per-line list or the first column of a TSV."""
    out = set()
    with open(path) as fh:
        for i, line in enumerate(fh):
            tok = line.strip().split("\t")[0].strip()
            if not tok or tok.startswith("#"):
                continue
            if i == 0 and tok in ("target", "candidate", "candidate_id"):
                continue
            out.add(tok)
    if not out:
        raise ValueError(f"no targets in {path}")
    return out



# ------------------------------------------------------------------- annotate
def cmd_annotate(args: argparse.Namespace) -> int:
    import json as _json
    from . import annotate as ann
    from .annotate import AnnotateError

    try:
        if args.probe:
            # Observe one real response before making a thousand requests on the
            # strength of an assumed shape.
            for label, url in (("UniProtKB", ann.UNIPROT.format(acc=args.probe)),
                               ("UniParc", ann.UNIPARC_BY_ACC.format(acc=args.probe))):
                payload, status = ann.fetch_json(url)
                print(f"\n=== {label}: {url}\n--- status: {status}")
                if payload is not None:
                    text = _json.dumps(payload, indent=2)
                    print(text[:args.probe_chars])
                    if len(text) > args.probe_chars:
                        print(f"... [{len(text) - args.probe_chars} more characters]")
                    if label == "UniProtKB" and status == "ok":
                        parsed = ann.parse_uniprot(args.probe, payload)
                        print("--- parsed by foldtrace:")
                        for k, v in zip(ann.ANNOTATION_COLUMNS, parsed.row()):
                            print(f"      {k:10} {v!r}")
                        print("    If any field above is empty but present in the JSON, the "
                              "parser needs fixing before a bulk run.")
            return 0

        accessions = ann.read_accessions(args.accessions)
        rows, man = ann.annotate(accessions, cache_dir=args.cache, clans=not args.no_clans,
                                 pause=args.pause)
        man.output = args.out
        man.output_sha256 = ann.write_annotations(rows, args.out)
        ann.write_manifest(man, args.manifest)
    except (AnnotateError, FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print(f"wrote {args.out}: {len(rows)} rows for {man.requested} accessions")
    for status, n in sorted(man.counts.items(), key=lambda kv: -kv[1]):
        print(f"  {status}: {n}")
    print(f"wrote {args.manifest} (sha256 {man.output_sha256[:16]}...) -- "
          f"pass it to `prime freeze --also {args.manifest}`")
    bad = sum(n for s, n in man.counts.items() if s not in ("ok", "withdrawn_recovered"))
    if bad:
        print(f"  note: {bad} accessions have a non-ok status; these are retrieval "
              "outcomes, not biological absences", file=sys.stderr)
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
    r.add_argument("--no-advice", action="store_true",
                   help="suppress the 'not looked at yet' note after the run")
    r.set_defaults(func=cmd_run)

    # guided: a paced, checkpointed investigation of one candidate (a JSON journal)
    g = sub.add_parser("guided", help="paced observe->predict->compute->compare->decide workflow")
    g.set_defaults(func=cmd_guided)
    gsub = g.add_subparsers(dest="guided_cmd", required=True)

    gi = gsub.add_parser("init", help="create a project journal")
    gi.add_argument("--project", required=True, help="journal file to create (JSON)")
    gi.add_argument("--name", required=True, help="short label for the investigation")
    gi.add_argument("--reference", required=True, help="reference structure (PDB/mmCIF)")
    gi.add_argument("--sites", required=True, help="prespecified catalytic-sites TSV")
    gi.add_argument("--candidate", required=True, help="the one candidate structure to study")
    _add_map_opts(gi)

    for name, help_ in (("observe", "record observations; print the briefing"),
                        ("predict", "commit a state per site (label=state)"),
                        ("compute", "run FOLDTRACE mapping; locks the prediction"),
                        ("compare", "score the prediction against the result"),
                        ("decide", "record a conclusion and next step"),
                        ("status", "show which stages are done"),
                        ("report", "print a Markdown write-up")):
        gp = gsub.add_parser(name, help=help_)
        gp.add_argument("--project", required=True, help="journal file (JSON)")
        if name == "observe":
            gp.add_argument("--notes", help="free-text observations")
        if name == "predict":
            gp.add_argument("--set", action="append", metavar="LABEL=STATE",
                            help="a site prediction, e.g. Glu642=lost (repeatable)")
            gp.add_argument("--rationale", help="reasoning for the prediction")
            gp.add_argument("--verdict", choices=["retained", "altered", "lost"],
                            help="optional overall verdict prediction")
        if name == "decide":
            gp.add_argument("--conclusion", required=True, help="what you conclude")
            gp.add_argument("--next", help="the next experiment/step")
        if name == "report":
            gp.add_argument("--out", help="write to a file instead of stdout")

    # prime: the prior-art register -- what was known before the search was run
    pr = sub.add_parser("prime", help="prior-art register: build, verify, freeze, adjudicate")
    pr.set_defaults(func=cmd_prime)
    psub = pr.add_subparsers(dest="prime_cmd", required=True)

    pi = psub.add_parser("init", help="create an empty register for a fold")
    pi.add_argument("--register", required=True, help="register file to create (JSON)")
    pi.add_argument("--fold", required=True, help="short fold id, e.g. hepn")
    pi.add_argument("--reference", required=True, help="reference structure, e.g. 5YEP:B")
    pi.add_argument("--scope", help="what the prior-art retrieval was asked to cover")
    pi.add_argument("--tool", action="append", help="retrieval tool used (repeatable)")
    pi.add_argument("--force", action="store_true", help="overwrite an existing register")

    pa = psub.add_parser("add-claim", help="add one prior-art claim (always unverified)")
    pa.add_argument("--register", required=True)
    pa.add_argument("--type", required=True, help="claim type (see prime.CLAIM_TYPES)")
    pa.add_argument("--claim", required=True, help="one sentence: what is known")
    pa.add_argument("--identifier", help="DOI, PMID, PDB id or Pfam accession")
    pa.add_argument("--residue", action="append", metavar="H102",
                    help="reference position this claim supports (repeatable)")
    pa.add_argument("--evidence", default="annotation",
                    help="experimental|structural|computational|annotation")
    pa.add_argument("--id", help="claim id (default: next PA-nnn)")
    pa.add_argument("--post-run", action="store_true", help="retrieved after the run")
    pa.add_argument("--trigger", help="for post-run claims: the feature that prompted retrieval")
    pa.add_argument("--note", help="free text")

    pc = psub.add_parser("criteria", help="append novelty / non-novelty criteria")
    pc.add_argument("--register", required=True)
    pc.add_argument("--novel", action="append", help="what would count as new (repeatable)")
    pc.add_argument("--not-novel", action="append",
                    help="what would NOT count as new however recovered (repeatable)")

    pv = psub.add_parser("verify", help="promote claims from an external verification file")
    pv.add_argument("--register", required=True)
    pv.add_argument("--evidence", required=True,
                    help="JSONL from the retrieval tool: claim_id, status, method")
    pv.add_argument("--fail-on-unresolved", action="store_true",
                    help="exit non-zero if any pre-run claim is still unverified")

    pd = psub.add_parser("validate", help="list every reason the register cannot be frozen")
    pd.add_argument("--register", required=True)
    pd.add_argument("--sites", help="also check that each declared site has backing evidence")

    ps = psub.add_parser("sites", help="emit a sites.tsv skeleton from verified evidence")
    ps.add_argument("--register", required=True)
    ps.add_argument("--out", help="write to a file instead of stdout")

    pf = psub.add_parser("freeze", help="validate, hash and lock; writes FREEZE.md")
    pf.add_argument("--register", required=True)
    pf.add_argument("--sites", required=True)
    pf.add_argument("--also", action="append", help="extra file to hash (repeatable)")
    pf.add_argument("--out", default="FREEZE.md")

    pj = psub.add_parser("adjudicate", help="record prior art for one emergent feature")
    pj.add_argument("--register", required=True)
    pj.add_argument("--feature", required=True, help="what the run surfaced")
    pj.add_argument("--surfaced-by", required=True, help="see prime.SURFACED_BY")
    pj.add_argument("--status", required=True, help="see prime.PRIOR_ART_STATUS")
    pj.add_argument("--permitted", required=True,
                    help="the strongest sentence the evidence supports")
    pj.add_argument("--described", help="required when partially_described")
    pj.add_argument("--undescribed", help="required when partially_described")
    pj.add_argument("--supports", action="append", metavar="PA-001",
                    help="backing claim id (repeatable)")

    pg = psub.add_parser("gate", help="refuse novelty language no adjudication supports")
    pg.add_argument("--register", required=True)
    pg.add_argument("--draft", required=True, help="draft text file to check")

    pe = psub.add_parser("export", help="flat TSV of the claim table")
    pe.add_argument("--register", required=True)
    pe.add_argument("--out", help="write to a file instead of stdout")

    # discord: rank disagreements between independent assertions about a target
    d = sub.add_parser("discord", help="rank discordances: annotation vs fold, channel vs channel")
    d.add_argument("--calls", required=True, help="a foldtrace map/run calls table")
    d.add_argument("--annotations", help="per-target metadata TSV (name, pfam, clan, taxon)")
    d.add_argument("--vocab", help="vocabulary TSV: class, term, note (see discord.load_vocab)")
    d.add_argument("--compare", help="second calls table for paired comparison")
    d.add_argument("--label-a", default="A", help="label for --calls (default: A)")
    d.add_argument("--label-b", default="B", help="label for --compare (default: B)")
    d.add_argument("--no-knockout", action="store_true", help="skip the point-knockout scan")
    d.add_argument("--knockout-tm", type=float, default=0.85,
                   help="min TM for the knockout signature (default: 0.85)")
    d.add_argument("--profile", help="a `foldtrace profile` table; enables the "
                                     "identical-sequence-across-organisms scan")
    d.add_argument("--min-identity", type=float, default=0.98,
                   help="identity across mapped positions for the contamination scan "
                        "(default: 0.98)")
    d.add_argument("--organism-key", default="organism",
                   help="column in --annotations holding the organism (default: organism)")
    d.add_argument("--knockout-max-prevalence", type=float, default=0.05,
                   help="above this share of eligible targets the knockout signature is "
                        "reported once as a recurrent single-site loss (default: 0.05)")
    d.add_argument("--min-score", type=float, default=0.0, help="drop findings below this score")
    d.add_argument("--kind", action="append", help="keep only this finding kind (repeatable)")
    d.add_argument("--emit-features", action="store_true",
                   help="print suggested 'prime adjudicate' calls for the leading findings")
    d.add_argument("--out", help="output TSV/CSV (default: stdout)")
    d.add_argument("--format", choices=["tsv", "csv"], default="tsv", help="output delimiter")
    d.set_defaults(func=cmd_discord)

    # profile: observed residue at EVERY reference position, for every candidate
    pp = sub.add_parser("profile", help="record the residue at every reference position (stage A)")
    pp.add_argument("--reference", required=True, help="reference structure (PDB/mmCIF)")
    pp.add_argument("--candidates", required=True, nargs="+", help="candidate structures or globs")
    pp.add_argument("--out", help="output profile TSV (default: stdout)")
    pp.add_argument("--format", choices=["tsv", "csv"], default="tsv")
    pp.add_argument("--ref-chain", default=None)
    pp.add_argument("--cand-chain", default=None)
    _add_map_opts(pp)
    pp.set_defaults(func=cmd_profile)

    # invariant: undeclared positions that hold while the declared chemistry does not
    iv = sub.add_parser("invariant", help="find undeclared near-site positions that do not vary")
    iv.add_argument("--profile", required=True, help="table from `foldtrace profile`")
    iv.add_argument("--reference", required=True, help="reference structure (PDB/mmCIF)")
    iv.add_argument("--sites", required=True, help="declared sites TSV (defines the site centre)")
    iv.add_argument("--groups", required=True,
                    help="calls table giving each candidate a verdict, or a candidate/verdict TSV")
    iv.add_argument("--clusters", help="candidate/cluster TSV; conservation is then computed "
                                       "per sequence cluster rather than per protein")
    iv.add_argument("--focus", default="lost",
                    help="group in which an invariant is unexpected (default: lost)")
    iv.add_argument("--compare", default="retained",
                    help="group to contrast against, or '' for none (default: retained)")
    iv.add_argument("--near", type=float, default=12.0,
                    help="site-zone radius from the declared-site centroid, A (default: 12)")
    iv.add_argument("--background", type=float, default=20.0,
                    help="background zone starts at this distance, A (default: 20)")
    iv.add_argument("--min-coverage", type=float, default=0.60,
                    help="min fraction of focus members resolving a position (default: 0.60)")
    iv.add_argument("--min-conservation", type=float, default=0.90,
                    help="min modal-residue frequency to report (default: 0.90)")
    iv.add_argument("--all", action="store_true", help="report every position, not only findings")
    iv.add_argument("--emit-features", action="store_true",
                    help="print a suggested 'prime adjudicate' call")
    iv.add_argument("--ref-chain", default=None)
    iv.add_argument("--out", help="output TSV/CSV (default: stdout)")
    iv.add_argument("--format", choices=["tsv", "csv"], default="tsv")
    iv.set_defaults(func=cmd_invariant)

    # rare: chemistry carried by a small minority at a near-site position
    rr = sub.add_parser("rare", help="find rare chemistry-capable residues near the site")
    rr.add_argument("--profile", required=True, help="table from `foldtrace profile`")
    rr.add_argument("--reference", required=True)
    rr.add_argument("--sites", required=True, help="declared sites TSV (defines the site centre)")
    rr.add_argument("--clusters", help="candidate/cluster TSV; carriers counted per cluster")
    rr.add_argument("--residues", default="HCDEKRSTYW",
                    help="residues treated as chemistry-capable (default: HCDEKRSTYW)")
    rr.add_argument("--min-count", type=int, default=3, help="min carriers (default: 3)")
    rr.add_argument("--max-fraction", type=float, default=0.10,
                    help="max share of members for a residue to count as rare (default: 0.10)")
    rr.add_argument("--min-modal", type=float, default=0.50,
                    help="min modal-residue frequency for a position to be testable; "
                         "below it the position is diffuse and rarity is uninformative "
                         "(default: 0.50)")
    rr.add_argument("--min-positions", type=int, default=2,
                    help="min rare positions per member for a combination (default: 2)")
    rr.add_argument("--near", type=float, default=12.0)
    rr.add_argument("--background", type=float, default=20.0)
    rr.add_argument("--include-declared", action="store_true",
                    help="also test the declared positions")
    rr.add_argument("--combos", help="write the multi-position combinations here")
    rr.add_argument("--emit-features", action="store_true")
    rr.add_argument("--ref-chain", default=None)
    rr.add_argument("--out", help="output TSV/CSV (default: stdout)")
    rr.add_argument("--format", choices=["tsv", "csv"], default="tsv")
    rr.set_defaults(func=cmd_rare)

    # cooccur: features that travel with the fold more often than background
    cc = sub.add_parser("cooccur", help="rank features enriched among the hits")
    cc.add_argument("--features", required=True,
                    help="long TSV: target, feature_type, feature")
    cc.add_argument("--focus-targets", required=True, help="target list or TSV")
    cc.add_argument("--background-targets", help="held-out background target list")
    cc.add_argument("--held-out", action="store_true",
                    help="assert the background is independent of feature selection")
    cc.add_argument("--spread-by", help="a column in --features to count carrier spread over "
                                        "(e.g. genus), to expose ancestry-tracking features")
    cc.add_argument("--min-count", type=int, default=3, help="min carriers (default: 3)")
    cc.add_argument("--min-odds", type=float, default=2.0, help="min odds ratio (default: 2.0)")
    cc.add_argument("--max-q", type=float, default=0.02,
                    help="Benjamini-Hochberg q cutoff across all tested features; calibrated "
                         "against random splits of a signal-free set to keep at least 90%% of "
                         "null splits empty (default: 0.02)")
    cc.add_argument("--alternative", choices=["greater", "two-sided", "less"], default="greater",
                    help="Fisher alternative; 'greater' matches the directional --min-odds "
                         "filter and is the default")
    cc.add_argument("--emit-features", action="store_true")
    cc.add_argument("--out", help="output TSV/CSV (default: stdout)")
    cc.add_argument("--format", choices=["tsv", "csv"], default="tsv")
    cc.set_defaults(func=cmd_cooccur)

    # annotate: the fetch layer -- the only stage that touches the network
    an = sub.add_parser("annotate", help="build annotations.tsv from UniProt/InterPro")
    an.add_argument("--accessions", help="one accession per line, or a TSV first column")
    an.add_argument("--out", default="annotations.tsv", help="output TSV")
    an.add_argument("--manifest", default="annotations_manifest.json",
                    help="provenance manifest for `prime freeze --also`")
    an.add_argument("--cache", help="directory for cached responses; makes a re-run free")
    an.add_argument("--no-clans", action="store_true", help="skip the InterPro clan lookup")
    an.add_argument("--pause", type=float, default=0.34,
                    help="seconds between requests (default: 0.34)")
    an.add_argument("--probe", metavar="ACCESSION",
                    help="fetch one record, print the raw JSON and how it parses, and stop")
    an.add_argument("--probe-chars", type=int, default=4000,
                    help="how much raw JSON to print with --probe (default: 4000)")
    an.set_defaults(func=cmd_annotate)

    # next: what this working directory has not been looked at yet
    nx = sub.add_parser("next", help="what has not been looked at yet in this directory")
    nx.add_argument("--dir", default=".", help="working directory (default: .)")
    nx.add_argument("--brief", action="store_true", help="show only the top two")
    nx.set_defaults(func=cmd_next)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except BrokenPipeError:
        # a downstream `head` closed the pipe; that is not an error
        try:
            sys.stdout.close()
        except Exception:
            pass
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
