# Changelog

## 0.2.0 (2026-08-24)

Completes the two-stage workflow described in the Nature Methods manuscript. Archived at
Zenodo, concept DOI [10.5281/zenodo.22081699](https://doi.org/10.5281/zenodo.22081699).

- **`foldtrace search`** (new, stage 1): wraps Foldseek `easy-search` against a Foldseek
  database (e.g. `afdb50` for AlphaFold/UniProt50), parses the output, applies optional
  hit-count / probability / coverage / identity thresholds, and writes a ranked, machine-
  readable hit table that preserves every Foldseek target identifier. Foldseek runs in a
  subprocess; an existing Foldseek table can be parsed with `--from-foldseek`.
- **`foldtrace run`** (new): end-to-end stage 1 + stage 2. Obtains hits (by running Foldseek,
  or from `--from-foldseek` / `--hits` / explicit `--candidates`), resolves a structure per
  hit, maps it, and writes one merged record carrying the Foldseek metrics, the TM-align
  score, per-site residue/offset/state, the verdict, and an explicit `unresolved_reason`.
  Candidates whose structure is missing are reported as unresolved, never dropped.
- **`foldtrace map`** (unchanged interface, backward compatible): now also emits an
  `unresolved_reason` column and supports an optional `altered` state.
- New **`altered`** active-site state (four states total: retained / altered / lost /
  unresolved) via an optional `altered` column in the sites TSV.
- `--format tsv|csv` on `map` and `run`; documented default thresholds
  (offset 4.0 A, TM gate 0.5; search thresholds off by default).
- Worked end-to-end example (`examples/tir/run_end_to_end.sh`) and expanded test suite
  (search parsing/ranking, end-to-end run, unresolved reporting, all four states).

## 0.1.0 (2026-08-24)

First release. Implements the active-site-mapping stage of structure-first search and
active-site mapping.

- `foldtrace map`: superpose candidate structures on a reference with TM-align (order-aware
  residue correspondence via tmtools) and read each chemistry-defining position as
  retained / lost / unresolved, with a per-candidate verdict.
- Worked TIR demo (SARM1 NADase site, four AlphaFold candidates, ~1 s) with checked-in
  expected output and an end-to-end test suite.
- Pinned dependencies (NumPy, Biopython 1.82, tmtools 0.3.0); MIT licence; CI on
  Python 3.9/3.11/3.12.

Planned for 0.2: a native `foldtrace search` wrapper over Foldseek, and an `altered` state
for conservative-but-changed chemistry.
