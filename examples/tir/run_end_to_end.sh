#!/usr/bin/env bash
# foldtrace end-to-end demo (stages 1 + 2) on the SARM1/TIR example.
#
# Stage 1 (search) here is supplied as a precomputed foldtrace search table
# (example_search_table.tsv) so the demo runs offline; with Foldseek installed the same
# table is produced by:
#   foldtrace search --query reference/SARM1_TIR_6O0R_A.pdb --database afdb50 --out example_search_table.tsv
# Stage 2 (map) superposes each hit on the SARM1 reference and reads the NADase site.
# Runtime: ~1 s for the four candidates on one CPU core.
set -euo pipefail
cd "$(dirname "$0")"
if command -v foldtrace >/dev/null 2>&1; then FOLDTRACE="foldtrace"; else FOLDTRACE="python3 -m foldtrace"; fi

$FOLDTRACE run \
  --reference reference/SARM1_TIR_6O0R_A.pdb \
  --sites catalytic_sites.tsv \
  --hits example_search_table.tsv \
  --candidates-dir candidates \
  --out end_to_end_output.tsv

echo
echo "End-to-end result (search metrics + active-site state; A0A933QTC5 = lost):"
column -t -s $'\t' end_to_end_output.tsv
