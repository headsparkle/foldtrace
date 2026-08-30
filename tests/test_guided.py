"""FOLDTRACE Guided: stage ordering, the prediction lock, scoring, and persistence.

Uses the bundled TIR example so the tests run offline (no Foldseek, no network).
"""
import os

import pytest

from foldtrace.guided import (
    GuidedProject, CheckpointError, PredictionLockedError,
    STAGE_IDS, STAGES, TM_GATE, OFFSET_A,
)
from foldtrace.mapping import RETAINED, LOST

HERE = os.path.dirname(__file__)
TIR = os.path.join(HERE, "..", "examples", "tir")
REF = os.path.join(TIR, "reference", "SARM1_TIR_6O0R_A.pdb")
SITES = os.path.join(TIR, "catalytic_sites.tsv")
LOST_CAND = os.path.join(TIR, "candidates", "A0A933QTC5.pdb")   # Glu642 -> lost
KEPT_CAND = os.path.join(TIR, "candidates", "A0A953THP3.pdb")   # Glu642 -> retained


def _project(cand=LOST_CAND, name="t"):
    return GuidedProject(name, REF, SITES, cand)


def test_constants_and_stage_order():
    assert STAGE_IDS == ("observe", "predict", "compute", "compare", "decide")
    assert tuple(s.id for s in STAGES) == STAGE_IDS
    assert (TM_GATE, OFFSET_A) == (0.5, 4.0)


def test_checkpoints_enforce_order():
    p = _project()
    with pytest.raises(CheckpointError):
        p.predict({l: RETAINED for l in p.site_labels})   # before observe
    p.observe()
    with pytest.raises(CheckpointError):
        p.compute()                                        # before predict
    p.predict({l: RETAINED for l in p.site_labels})
    with pytest.raises(CheckpointError):
        p.compare()                                        # before compute


def test_prediction_can_be_revised_until_compute_then_locks():
    p = _project()
    p.observe()
    p.predict({"Tyr568": RETAINED, "Trp638": RETAINED, "Glu642": RETAINED})
    # revision allowed before compute
    rec = p.predict({"Tyr568": RETAINED, "Trp638": RETAINED, "Glu642": LOST})
    assert rec["site_predictions"]["Glu642"] == LOST
    p.compute()
    with pytest.raises(PredictionLockedError):
        p.predict({l: RETAINED for l in p.site_labels})


def test_predict_validates_labels_and_states():
    p = _project()
    p.observe()
    with pytest.raises(ValueError):
        p.predict({"Glu642": RETAINED})                    # missing sites
    with pytest.raises(ValueError):
        p.predict({"Tyr568": RETAINED, "Trp638": RETAINED,
                   "Glu642": RETAINED, "Bogus999": LOST})  # unknown label
    with pytest.raises(ValueError):
        p.predict({"Tyr568": RETAINED, "Trp638": RETAINED, "Glu642": "unresolved"})  # not predictable


def test_compare_scores_hits_and_surprises():
    p = _project(cand=LOST_CAND)
    p.observe()
    # correct on Glu642 (lost), wrong on Trp638 (predict retained, truly lost)
    p.predict({"Tyr568": RETAINED, "Trp638": RETAINED, "Glu642": LOST})
    result = p.compute()
    assert result.verdict == LOST
    sc = p.compare()
    assert sc["computed_verdict"] == LOST
    assert "Trp638" in sc["surprises"]
    assert sc["n_correct"] == sc["n_scored"] - len(sc["surprises"])


def test_full_workflow_completes_and_reports():
    p = _project(cand=KEPT_CAND)
    p.observe("Glu642 is the catalytic key.")
    p.predict({"Tyr568": RETAINED, "Trp638": RETAINED, "Glu642": RETAINED},
              rationale="Strong fold match; expect the site intact.", verdict=RETAINED)
    p.compute()
    p.compare()
    p.decide("Functional NADase candidate.", next_action="Order the construct.")
    st = p.status()
    assert st["complete"] and st["next"] is None
    assert "## Decision" in p.report()


def test_unresolved_sites_are_not_scored_against_the_student(tmp_path):
    # an impossibly tight offset makes every site unresolved; the student is not penalised
    p = GuidedProject("u", REF, SITES, LOST_CAND, offset_threshold=0.0)
    p.observe()
    p.predict({"Tyr568": RETAINED, "Trp638": RETAINED, "Glu642": LOST})
    p.compute()
    sc = p.compare()
    assert sc["n_scored"] == 0
    assert sc["n_unresolved"] == len(p.site_labels)


def test_cli_end_to_end(tmp_path, capsys):
    from foldtrace.cli import main
    j = str(tmp_path / "j.json")
    assert main(["guided", "init", "--project", j, "--name", "c",
                 "--reference", REF, "--sites", SITES, "--candidate", LOST_CAND]) == 0
    assert main(["guided", "observe", "--project", j]) == 0
    assert main(["guided", "predict", "--project", j,
                 "--set", "Tyr568=retained", "--set", "Trp638=retained",
                 "--set", "Glu642=lost"]) == 0
    assert main(["guided", "compute", "--project", j]) == 0
    # the lock is enforced through the CLI too
    assert main(["guided", "predict", "--project", j,
                 "--set", "Tyr568=retained", "--set", "Trp638=retained",
                 "--set", "Glu642=retained"]) == 2
    assert main(["guided", "compare", "--project", j]) == 0
    assert main(["guided", "decide", "--project", j, "--conclusion", "lost"]) == 0
    capsys.readouterr()
    assert main(["guided", "report", "--project", j]) == 0
    assert "FOLDTRACE Guided" in capsys.readouterr().out


def test_save_and_load_roundtrip(tmp_path):
    p = _project(cand=LOST_CAND, name="roundtrip")
    p.observe()
    p.predict({"Tyr568": RETAINED, "Trp638": RETAINED, "Glu642": LOST})
    p.compute()
    path = tmp_path / "journal.json"
    p.save(str(path))
    q = GuidedProject.load(str(path))
    assert q.name == "roundtrip"
    assert q.status()["prediction_locked"] is True
    # continues from where it left off
    q.compare()
    assert q.status()["done"][:4] == ["observe", "predict", "compute", "compare"]
    # locked prediction survives reload
    with pytest.raises(PredictionLockedError):
        q.predict({l: RETAINED for l in q.site_labels})
