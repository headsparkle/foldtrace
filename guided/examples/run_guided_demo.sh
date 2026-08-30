#!/usr/bin/env bash
# FOLDTRACE Guided demo: one full observe->predict->compute->compare->decide pass
# on the bundled TIR example (offline; no Foldseek, no network).
#
# Run from the repo root:  bash guided/examples/run_guided_demo.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

FT() { python3 -m foldtrace "$@"; }   # falls back to the module if the console script isn't on PATH
command -v foldtrace >/dev/null 2>&1 && FT() { foldtrace "$@"; }

J="$(mktemp -d)/tir_demo.json"
REF="examples/tir/reference/SARM1_TIR_6O0R_A.pdb"
SITES="examples/tir/catalytic_sites.tsv"
CAND="examples/tir/candidates/A0A933QTC5.pdb"   # this candidate has LOST the catalytic Glu642

echo "== init =="
FT guided init --project "$J" --name tir-demo --reference "$REF" --sites "$SITES" --candidate "$CAND"

echo "== observe =="
FT guided observe --project "$J" --notes "SARM1 TIR NADase; Glu642 is the catalytic key."

echo "== predict (deliberately wrong on Trp638 to show a 'surprise') =="
FT guided predict --project "$J" \
  --set Tyr568=retained --set Trp638=retained --set Glu642=lost \
  --rationale "Good fold match, but the catalytic Glu looks substituted." --verdict lost

echo "== compute (locks the prediction) =="
FT guided compute --project "$J"

echo "== compare =="
FT guided compare --project "$J"

echo "== decide =="
FT guided decide --project "$J" \
  --conclusion "Candidate has lost the NADase Glu642; not a functional homolog." \
  --next "Screen the next hit A0A953THP3 (predicted retained)."

echo "== report =="
FT guided report --project "$J"
