# Changelog

## 0.1.0 (2026-08-24)

First release. Implements the active-site-mapping stage of structure-first search and
active-site mapping.

- `foldsite map`: superpose candidate structures on a reference with TM-align (order-aware
  residue correspondence via tmtools) and read each chemistry-defining position as
  retained / lost / unresolved, with a per-candidate verdict.
- Worked TIR demo (SARM1 NADase site, four AlphaFold candidates, ~1 s) with checked-in
  expected output and an end-to-end test suite.
- Pinned dependencies (NumPy, Biopython 1.82, tmtools 0.3.0); MIT licence; CI on
  Python 3.9/3.11/3.12.

Planned for 0.2: a native `foldsite search` wrapper over Foldseek, and an `altered` state
for conservative-but-changed chemistry.
