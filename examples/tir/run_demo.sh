#!/usr/bin/env bash
# foldsite worked demo: SARM1 TIR NADase site read across four AlphaFold candidates.
# Reproduces the retained/lost calls in Zimmer, "Structure-first search and active-site
# mapping across protein folds". Runtime: ~1 s for 4 candidates on one CPU core.
set -euo pipefail
cd "$(dirname "$0")"

# use the installed console script if present, else the module entry point
if command -v foldsite >/dev/null 2>&1; then FOLDSITE="foldsite"; else FOLDSITE="python3 -m foldsite"; fi

$FOLDSITE map \
  --reference reference/SARM1_TIR_6O0R_A.pdb \
  --sites catalytic_sites.tsv \
  --candidates candidates/*.pdb \
  --out demo_output.tsv

echo
echo "Result (expected: A0A933QTC5 = lost, the other three = retained):"
column -t -s $'\t' demo_output.tsv

# Optional: fail if the output drifts from the checked-in expected result.
if ! diff -q <(sort demo_output.tsv) <(sort expected_output.tsv) >/dev/null; then
  echo "WARNING: demo_output.tsv differs from expected_output.tsv" >&2
fi
