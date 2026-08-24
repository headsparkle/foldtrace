"""End-to-end test on the bundled TIR demo: the known states must be reproduced."""
import os

import pytest

from foldsite.io import load_structure, parse_sites
from foldsite.mapping import map_candidate, RETAINED, LOST

HERE = os.path.dirname(__file__)
TIR = os.path.join(HERE, "..", "examples", "tir")
REF = os.path.join(TIR, "reference", "SARM1_TIR_6O0R_A.pdb")
SITES = os.path.join(TIR, "catalytic_sites.tsv")


def _run(acc):
    ref = load_structure(REF)
    sites = parse_sites(SITES)
    cand = load_structure(os.path.join(TIR, "candidates", f"{acc}.pdb"))
    return map_candidate(ref, cand, sites, candidate_name=acc)


@pytest.mark.parametrize("acc,expected", [
    ("A0A953THP3", RETAINED),   # full SARM1-like triad
    ("A0A7Y7TLD0", RETAINED),   # full triad
    ("Q9FHM1", RETAINED),       # TIR NADase site annotated as CD38-like
    ("A0A933QTC5", LOST),       # catalytic Glu naturally lost (Thr)
])
def test_verdict(acc, expected):
    r = _run(acc)
    assert r.fold_ok, f"{acc} did not clear the TM-score fold gate"
    assert r.verdict == expected, f"{acc}: got {r.verdict}, expected {expected}"


def test_glu_lost_reads_thr_and_leu():
    """The catalytic-lost member should read a non-Glu at 642 and a non-Trp at 638."""
    r = _run("A0A933QTC5")
    by = {c.label: c for c in r.sites}
    assert by["Glu642"].state == LOST and by["Glu642"].obs_res != "E"
    assert by["Trp638"].state == LOST and by["Trp638"].obs_res != "W"


def test_offsets_are_small_for_retained():
    """Order-aware correspondence should place retained catalytic residues within ~1.5 A."""
    r = _run("A0A953THP3")
    for c in r.sites:
        assert c.ca_offset is not None and c.ca_offset < 2.0
