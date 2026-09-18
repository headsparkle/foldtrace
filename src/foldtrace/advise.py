"""What has not been looked at yet.

A Kaede run once put four members of a previously described family at ranks 37
to 100 and said nothing about them. Everything needed to surface them was
present -- the models were on disk, the members were in the calls table, and the
detector that finds shared undeclared positions was installed -- but nothing ran,
because nothing told the user it could.

That is the failure this module exists to stop. The detectors are most useful
exactly when you do not know what you are looking for, which is precisely when
you will not think to invoke them.

It deliberately does not orchestrate. Each stage still writes a file you inspect
before the next, because the judgement in between is where the errors get
caught. What it removes is the need to already know what the next stage is:
:func:`advise` reads a working directory, works out which stages have run, and
prints the next commands with the reason each one matters for *this* set of
results.

Stdlib only; reads files, runs nothing.
"""
from __future__ import annotations

import csv
import glob
import json
import os
from dataclasses import dataclass, field


@dataclass
class Advice:
    """One thing worth doing next, and why it matters here."""
    title: str
    why: str
    commands: list[str] = field(default_factory=list)
    urgency: int = 2          # 1 highest


@dataclass
class State:
    """What exists in a working directory."""
    calls: str | None = None
    profile: str | None = None
    annotations: str | None = None
    clusters: str | None = None
    register: str | None = None
    freeze: str | None = None
    models: int = 0
    reference: str | None = None
    sites: str | None = None
    detector_outputs: list[str] = field(default_factory=list)
    verdicts: dict = field(default_factory=dict)
    n_candidates: int = 0


def _first(directory: str, *names: str) -> str | None:
    for n in names:
        p = os.path.join(directory, n)
        if os.path.exists(p):
            return p
    return None


def _read_verdicts(path: str) -> tuple[dict, int]:
    counts: dict[str, int] = {}
    total = 0
    try:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                v = (row.get("verdict") or "unresolved").strip().lower()
                counts[v] = counts.get(v, 0) + 1
                total += 1
    except (OSError, csv.Error):
        pass
    return counts, total


def inspect(directory: str = ".") -> State:
    s = State()
    s.calls = _first(directory, "calls.tsv", "calls.csv")
    s.profile = _first(directory, "profile.tsv")
    s.annotations = _first(directory, "annotations.tsv")
    s.clusters = _first(directory, "clusters.tsv")
    s.register = _first(directory, "register.json")
    s.freeze = _first(directory, "FREEZE.md")
    s.sites = _first(directory, "sites.tsv")
    s.reference = next(iter(sorted(glob.glob(os.path.join(directory, "*.pdb")))), None)
    for pattern in ("models/*.pdb", "models/*.cif"):
        s.models += len(glob.glob(os.path.join(directory, pattern)))
    for name in ("discord.tsv", "inv.tsv", "invariant.tsv", "rare.tsv", "cooccur.tsv"):
        hit = _first(directory, name)
        if hit:
            s.detector_outputs.append(os.path.basename(hit))
    if s.calls:
        s.verdicts, s.n_candidates = _read_verdicts(s.calls)
    return s


def _lost(state: State) -> int:
    return sum(n for v, n in state.verdicts.items() if v in ("lost", "altered"))


def advise(state: State, directory: str = ".") -> list[Advice]:
    """What to do next, most consequential first."""
    out: list[Advice] = []
    ref = os.path.basename(state.reference) if state.reference else "reference.pdb"
    sites = os.path.basename(state.sites) if state.sites else "sites.tsv"

    if not state.calls:
        out.append(Advice(
            "Map the declared chemistry",
            "Nothing has been mapped yet. `run` searches, fetches models and reads the "
            "declared positions in one pass.",
            [f"foldtrace run --reference {ref} --sites {sites} "
             "--from-foldseek hits.tsv --fetch --fetch-dir models --out calls.tsv"], 1))
        return out

    lost = _lost(state)

    # 1. the omission that motivated this module
    if not state.profile:
        why = ("`run` reads only the positions you declared. It cannot tell you what the "
               "candidates share at positions you did not declare -- which is how a buried "
               "network, or any undeclared shared chemistry, stays invisible.")
        if lost:
            why += (f" {lost} of {state.n_candidates} candidates here lost or altered the "
                    "declared site; whatever those members hold in common is computable "
                    "and has not been computed.")
        out.append(Advice(
            "Profile every position, not only the declared ones", why,
            [f"foldtrace profile --reference {ref} --candidates 'models/*.pdb' "
             "--out profile.tsv"],
            1 if lost else 2))
    else:
        ran = set(state.detector_outputs)
        if not any(n.startswith("inv") for n in ran):
            out.append(Advice(
                "Ask what the losing members hold in common",
                f"{lost} candidates lost or altered the declared site. `invariant` reports "
                "positions that do not vary among them but do vary elsewhere -- the step "
                "that finds a device nobody declared. Sweep the shell: a site-spanning "
                "network will not fit inside an arbitrary 12 A ball.",
                [f"foldtrace invariant --profile profile.tsv --reference {ref} "
                 f"--sites {sites} --groups calls.tsv "
                 + ("--clusters clusters.tsv " if state.clusters else "")
                 + "--focus lost --compare retained --near 12 --out inv_12.tsv",
                 f"foldtrace invariant --profile profile.tsv --reference {ref} "
                 f"--sites {sites} --groups calls.tsv "
                 + ("--clusters clusters.tsv " if state.clusters else "")
                 + "--focus lost --compare retained --near 18 --background 26 "
                   "--out inv_18.tsv"], 1))
        if not any(n.startswith("rare") for n in ran):
            out.append(Advice(
                "Look for chemistry a minority carries",
                "The inverse question: a residue present in a few members at a position "
                "that is otherwise consistent. This is the distribution annotation handles "
                "worst, and the threshold is worth sweeping -- a real case has sat 0.008 "
                "the wrong side of the default.",
                [f"foldtrace rare --profile profile.tsv --reference {ref} --sites {sites} "
                 + ("--clusters clusters.tsv " if state.clusters else "")
                 + "--max-fraction 0.10 --out rare_10.tsv --combos combos.tsv",
                 f"foldtrace rare --profile profile.tsv --reference {ref} --sites {sites} "
                 + ("--clusters clusters.tsv " if state.clusters else "")
                 + "--max-fraction 0.15 --out rare_15.tsv"], 2))

    # 2. contamination: cheap, and it invalidates findings rather than adding them
    if state.profile and "discord.tsv" not in state.detector_outputs:
        out.append(Advice(
            "Check for records that contradict each other",
            "Includes the contamination scan: sequences that are near-identical but "
            "attributed to different organisms. An engineered construct deposited under a "
            "natural accession passes every site check with a flawless intact site, so "
            "this does not add findings -- it removes false ones.",
            ["foldtrace discord --calls calls.tsv --profile profile.tsv "
             + (f"--annotations {os.path.basename(state.annotations)} --vocab vocab.tsv "
                if state.annotations else "")
             + "--min-identity 0.95 --out discord.tsv"], 1))

    # 3. a second reference answers a different question
    if state.calls and lost:
        out.append(Advice(
            "Consider a second reference",
            "A verdict is relative to the chemistry you declared. A candidate can read "
            "`lost` against one reference and `retained` against another because the two "
            "ask different questions. Mapping the same candidates against a second "
            "reference and comparing the two tables turns that into a ranked list of "
            "candidates the references disagree about.",
            ["foldtrace run --reference other.pdb --sites other_sites.tsv "
             "--candidates 'models/*.pdb' --out calls_other.tsv",
             "foldtrace discord --calls calls.tsv --compare calls_other.tsv "
             "--label-a first --label-b second --out ref_discord.tsv"], 3))

    # 4. supporting inputs
    if not state.clusters:
        out.append(Advice(
            "Build a sequence clustering",
            "Without `--clusters` every detector treats near-identical proteins as "
            "independent observations, and an expanded clade reads as conservation. "
            "Two columns: candidate, cluster.",
            ["# e.g. MMseqs2 at 50% identity, 80% coverage, then write candidate/cluster"], 2))
    if not state.annotations:
        out.append(Advice(
            "Fetch the annotations",
            "Names, Pfam accessions, clan membership and organism. Needed by `discord` and "
            "by anything asking what travels with the fold.",
            ["cut -f1 calls.tsv | tail -n +2 > accessions.txt",
             "foldtrace annotate --accessions accessions.txt --cache .cache "
             "--out annotations.tsv --manifest annotations_manifest.json"], 2))

    # 5. the record
    if not state.register:
        out.append(Advice(
            "Record what was already known",
            "Nothing yet states what would count as new for this fold. Written afterwards, "
            "that judgement is made knowing the answer.",
            ["foldtrace prime init --register register.json --fold <name> "
             f"--reference {ref}"], 3))
    elif not state.freeze:
        out.append(Advice(
            "Freeze the register",
            "The register exists but is not locked, so there is no timestamped record of "
            "what was known before the results.",
            [f"foldtrace prime freeze --register register.json --sites {sites}"], 2))
    elif state.detector_outputs:
        out.append(Advice(
            "Adjudicate before writing",
            "Each detector finding needs its prior art recorded before it appears in a "
            "draft. `prime gate` will refuse novelty language that no adjudication "
            "supports.",
            ["foldtrace discord --calls calls.tsv --profile profile.tsv --emit-features",
             "foldtrace prime gate --register register.json --draft draft.md"], 3))

    return sorted(out, key=lambda a: a.urgency)


def format_advice(items: list[Advice], state: State, brief: bool = False) -> str:
    if not items:
        return "Nothing obviously left undone in this directory."
    lines: list[str] = []
    if state.calls and state.verdicts:
        counts = ", ".join(f"{n} {v}" for v, n in sorted(state.verdicts.items()))
        lines.append(f"{state.n_candidates} candidates mapped ({counts})"
                     + (f", {state.models} models on disk" if state.models else ""))
        lines.append("")
    shown = items[:2] if brief else items
    for i, a in enumerate(shown, start=1):
        lines.append(f"{i}. {a.title}")
        lines.append(f"   {a.why}")
        for c in a.commands:
            lines.append(f"     $ {c}")
        lines.append("")
    if brief and len(items) > len(shown):
        lines.append(f"   ({len(items) - len(shown)} more: foldtrace next)")
    return "\n".join(lines).rstrip() + "\n"


def next_steps(directory: str = ".", brief: bool = False) -> str:
    state = inspect(directory)
    return format_advice(advise(state, directory), state, brief=brief)
