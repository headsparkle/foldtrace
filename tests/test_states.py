"""Explicit coverage of all four active-site states: retained / altered / lost / unresolved."""
import os

from foldtrace.io import load_structure, Site
from foldtrace.mapping import map_candidate, RETAINED, ALTERED, LOST, UNRESOLVED

HERE = os.path.dirname(__file__)
TIR = os.path.join(HERE, "..", "examples", "tir")
REF = os.path.join(TIR, "reference", "SARM1_TIR_6O0R_A.pdb")


def _ref():
    return load_structure(REF)


def _cand(acc):
    return load_structure(os.path.join(TIR, "candidates", f"{acc}.pdb"))


def test_retained_and_lost_on_glu642():
    ref = _ref()
    site = [Site("Glu642", 642, frozenset({"E"}), "catalytic", True)]
    assert map_candidate(ref, _cand("A0A953THP3"), site).verdict == RETAINED
    assert map_candidate(ref, _cand("A0A933QTC5"), site).verdict == LOST


def test_altered_state_when_residue_in_altered_set():
    # A0A933QTC5 has Thr at the Glu642 position; declaring T as 'altered' yields an altered call
    ref = _ref()
    site = [Site("Glu642", 642, frozenset({"E"}), "catalytic", True, altered=frozenset({"T"}))]
    r = map_candidate(ref, _cand("A0A933QTC5"), site)
    assert r.verdict == ALTERED
    assert r.sites[0].state == ALTERED


def test_unresolved_by_offset_threshold():
    # an impossibly tight offset threshold forces an unresolved call with a reason
    ref = _ref()
    site = [Site("Glu642", 642, frozenset({"E"}), "catalytic", True)]
    r = map_candidate(ref, _cand("A0A953THP3"), site, offset_threshold=0.0)
    assert r.verdict == UNRESOLVED
    assert "offset" in r.sites[0].reason


def test_unresolved_by_fold_gate():
    ref = _ref()
    site = [Site("Glu642", 642, frozenset({"E"}), "catalytic", True)]
    r = map_candidate(ref, _cand("A0A953THP3"), site, tm_gate=0.99)
    assert r.verdict == UNRESOLVED
    assert not r.fold_ok
    assert "TM gate" in r.state_reason
