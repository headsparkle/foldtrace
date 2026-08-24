# foldtrace

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22081699.svg)](https://doi.org/10.5281/zenodo.22081699)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Structure-first search and active-site mapping: find remote homologs of a fold and read their catalytic state from predicted structures.**

Author: **Marc Zimmer**, Department of Chemistry, **Connecticut College** — ORCID [0000-0001-8460-9064](https://orcid.org/0000-0001-8460-9064).
Concept DOI: [10.5281/zenodo.22081699](https://doi.org/10.5281/zenodo.22081699). Licence: MIT.

Sequence homology assigns a protein the function of its closest annotated relative. That fails
when a familiar fold is retained after the chemistry inside it has been altered, repurposed, or
lost. foldtrace implements the two-stage workflow of Zimmer, *Structure-first search and
active-site mapping across protein folds*:

1. **`foldtrace search`** — searches a Foldseek database (e.g. the AlphaFold/UniProt50 set) by
   fold and returns a ranked hit table.
2. **`foldtrace map`** — superposes each candidate on a reference with TM-align and reads every
   chemistry-defining position as **retained**, **altered**, **lost**, or **unresolved**.
3. **`foldtrace run`** — does search then mapping end-to-end.

The mapping is *order-aware*: rather than taking the residue whose Cα sits nearest a reference
catalytic atom (which can grab a spatially close but sequence-unrelated residue), foldtrace uses
the TM-align residue correspondence, which preserves sequence order, then checks that the
corresponding residue lies within a distance cutoff of the reference position. A call is emitted
only when fold correspondence and spatial placement agree; otherwise the site is reported as
`unresolved` **with an explicit reason, never dropped**.

## Install

```bash
git clone https://github.com/headsparkle/foldtrace
cd foldtrace
pip install -e .
```

Python dependencies are pinned in `requirements.txt` (NumPy, Biopython 1.82, tmtools 0.3.0);
TM-align runs in-process via `tmtools`, no external binary needed. **Foldseek** is required only
for the live `search` stage (https://github.com/steineggerlab/foldseek); everything else, and the
worked example, runs without it. Python >= 3.9.

## Quick start (end-to-end, ~1 s)

```bash
bash examples/tir/run_end_to_end.sh
```

This runs `foldtrace run` on the SARM1 TIR NADase site across four bundled AlphaFold candidates.
Stage 1 is supplied by `examples/tir/example_search_table.tsv`, a **demonstration fixture** — four
hand-picked candidates chosen to span active-site states, *not* an unbiased live AlphaFold/Foldseek
search (run `foldtrace search --database afdb50` for that). The example reproduces the calls in the
manuscript: A0A933QTC5 = **lost** (catalytic Glu → Thr; pocket Trp → Leu), the other three =
**retained** (including Q9FHM1, a TIR NADase site currently annotated as a CD38-like enzyme).

## Commands

### `foldtrace search` (stage 1)

```bash
foldtrace search --query reference/SARM1_TIR_6O0R_A.pdb --database afdb50 --out hits.tsv
# offline: parse an existing Foldseek table instead of running Foldseek
foldtrace search --from-foldseek foldseek_raw.tsv --out hits.tsv
```

**Inputs:** `--query` structure (PDB/mmCIF) and a Foldseek `--database` (e.g. `afdb50`), or
`--from-foldseek <table>`. **Output:** a ranked TSV/CSV with columns `rank, query, target, prob,
evalue, bits, fident, alnlen, qcov, tcov, alntmscore` — every Foldseek target id preserved for
the mapping stage. Foldseek is invoked with a fixed `--format-output` so parsing is deterministic.
**Options:** `--max-hits` (default 1000), `--min-prob`, `--min-coverage`, `--min-identity`
(all 0 = off by default), `--foldseek <binary>`, `--format tsv|csv`.

### `foldtrace map` (stage 2)

```bash
foldtrace map --reference ref.pdb --sites sites.tsv --candidates models/*.pdb --out calls.tsv
```

**Inputs:** `--reference` structure, a `--sites` TSV, and `--candidates` structures (files/globs).
**Output:** per candidate — `tm_norm_ref, rmsd, fold_ok`, then `{site}_obs, {site}_offsetA,
{site}_state` for each site, then `verdict` and `unresolved_reason`.
**Sites TSV** (tab-separated, header required):

```
label    ref_resnum    expected    role               key    altered
Tyr568   568           Y,F         pocket_aromatic    0
Trp638   638           W           nicotinamide_stack 0
Glu642   642           E           catalytic          1
```

`expected` = residues counting as *retained* (comma-separated); optional `altered` = residues
counting as *altered* (changed but related chemistry); `key=1` marks the residue that drives the
verdict. **Options:** `--offset-threshold` (default **4.0 Å**), `--tm-gate` (default **0.5**),
`--ref-chain`, `--cand-chain`, `--format tsv|csv`.

### `foldtrace run` (stages 1 + 2)

```bash
# from a live Foldseek search, fetching hit models from AFDB
foldtrace run --reference ref.pdb --sites sites.tsv --database afdb50 --fetch --out calls.tsv
# from a precomputed search table + local structures (offline)
foldtrace run --reference ref.pdb --sites sites.tsv --hits hits.tsv --candidates-dir models/ --out calls.tsv
```

Stage-1 source is one of `--database` (runs Foldseek), `--from-foldseek`, `--hits` (a foldtrace
search table), or `--candidates` (skip search). Hit structures are resolved from
`--candidates-dir` (files named `{target}.pdb`) or downloaded with `--fetch`. **Output** merges
the search metrics and the mapping result per candidate: `candidate_id, rank, foldseek_prob,
foldseek_evalue, foldseek_bits, foldseek_fident, foldseek_qcov, foldseek_alntmscore,
tmalign_tm_norm_ref, tmalign_rmsd, fold_ok, {site}_obs/_offsetA/_state ..., verdict,
unresolved_reason`. A hit whose structure cannot be found is reported with verdict `unresolved`
and reason `structure_not_found`.

## Default thresholds

| stage | parameter | default | meaning |
|---|---|---|---|
| search | `--max-hits` | 1000 | hits kept (matches the manuscript afdb50 runs) |
| search | `--min-prob` / `--min-coverage` / `--min-identity` | 0.0 (off) | optional hit filters; nothing dropped unless set |
| map | `--offset-threshold` | 4.0 Å | max Cα-Cα offset for a trusted call; beyond it → unresolved |
| map | `--tm-gate` | 0.5 | min reference-normalised TM-score for same-fold; below → unresolved |

## Tests

```bash
pip install -e ".[test]"
pytest
```

## Citing

Cite the concept DOI [10.5281/zenodo.22081699](https://doi.org/10.5281/zenodo.22081699) and see
`CITATION.cff`. TM-align: Zhang & Skolnick, *Nucleic Acids Res.* 33, 2302 (2005); `tmtools`:
https://github.com/jvkersch/tmtools; Foldseek: van Kempen et al., *Nat. Biotechnol.* 42, 243 (2024).

## License

MIT (see `LICENSE`).
