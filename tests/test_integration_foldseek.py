"""Conditional integration test that invokes the real Foldseek executable.

Skips gracefully when Foldseek is not installed, so the suite passes everywhere; when Foldseek
is on PATH it runs a genuine ``foldtrace search`` of the SARM1 TIR reference against the bundled
candidate structures (used directly as the target set) and checks the known targets come back.
"""
import os
import shutil

import pytest

from foldtrace.search import search

HERE = os.path.dirname(__file__)
TIR = os.path.join(HERE, "..", "examples", "tir")
REF = os.path.join(TIR, "reference", "SARM1_TIR_6O0R_A.pdb")
CANDS = os.path.join(TIR, "candidates")

foldseek = shutil.which("foldseek")


@pytest.mark.skipif(foldseek is None, reason="foldseek not installed")
def test_real_foldseek_search_recovers_known_targets(tmp_path):
    # foldseek easy-search accepts a directory of structures as the target set
    hits = search(query=REF, database=CANDS, out=str(tmp_path / "hits.tsv"),
                  tmp_dir=str(tmp_path / "fs"), max_hits=50, foldseek_bin=foldseek)
    targets = {h.target for h in hits}
    # at least one of the bundled TIR candidates must be recovered structurally
    known = {"A0A953THP3", "A0A933QTC5", "A0A7Y7TLD0", "Q9FHM1"}
    assert targets & known, f"expected some of {known}, got {targets}"
    # every returned hit must carry the Foldseek metrics the mapping stage needs
    for h in hits:
        assert h.alntmscore is not None and h.target
