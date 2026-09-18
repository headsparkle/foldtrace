"""Tests for the 'what has not been looked at yet' advisory.

The case that motivated it: a run where members of a known family sit in the
calls table with a lost verdict, the models are on disk, and nothing has told
the user that what those members share is computable.
"""
from __future__ import annotations

import os

from foldtrace import advise


def _calls(tmp_path, rows, name="calls.tsv"):
    p = tmp_path / name
    lines = ["candidate\ttmalign_tm_norm_ref\tverdict\tstate_reason"]
    for cand, verdict in rows:
        lines.append(f"{cand}\t0.90\t{verdict}\t")
    p.write_text("\n".join(lines) + "\n")
    return p


def _kaede_like(tmp_path):
    """A run that has mapped, with a substantial lost group, and nothing else."""
    (tmp_path / "kaede.pdb").write_text("ATOM\n")
    (tmp_path / "sites.tsv").write_text(
        "label\tref_resnum\texpected\trole\tkey\taltered\nH65\t65\tH\tx\t1\t\n")
    models = tmp_path / "models"
    models.mkdir()
    for i in range(5):
        (models / f"m{i}.pdb").write_text("ATOM\n")
    rows = [(f"c{i:03d}", "retained" if i % 3 else "lost") for i in range(60)]
    _calls(tmp_path, rows)
    return tmp_path


# ------------------------------------------------------------------- inspect
def test_inspect_finds_what_exists(tmp_path):
    d = _kaede_like(tmp_path)
    s = advise.inspect(str(d))
    assert s.calls and s.sites and s.reference
    assert s.models == 5
    assert s.n_candidates == 60
    assert s.verdicts["lost"] == 20
    assert s.profile is None and s.clusters is None


def test_inspect_of_an_empty_directory_is_harmless(tmp_path):
    s = advise.inspect(str(tmp_path))
    assert s.calls is None and s.n_candidates == 0
    assert advise.advise(s, str(tmp_path))[0].title.startswith("Map")


# ------------------------------------------------- the motivating omission
def test_profiling_is_the_first_advice_after_a_run_with_losses(tmp_path):
    """The Kaede case: four family members sat at ranks 37-100 and nothing said so."""
    d = _kaede_like(tmp_path)
    items = advise.advise(advise.inspect(str(d)), str(d))
    top = items[0]
    assert "Profile every position" in top.title
    assert "did not declare" in top.why
    assert "20 of 60" in top.why, "the advice names the size of the unexamined group"
    assert any("foldtrace profile" in c for c in top.commands)


def test_advice_moves_on_once_profiling_is_done(tmp_path):
    d = _kaede_like(tmp_path)
    (d / "profile.tsv").write_text("candidate\tref_resnum\tobs_res\tresolved\n")
    titles = [a.title for a in advise.advise(advise.inspect(str(d)), str(d))]
    assert not any("Profile every position" in t for t in titles)
    assert any("hold in common" in t for t in titles)


def test_the_invariant_advice_sweeps_the_shell(tmp_path):
    d = _kaede_like(tmp_path)
    (d / "profile.tsv").write_text("candidate\tref_resnum\tobs_res\tresolved\n")
    inv = next(a for a in advise.advise(advise.inspect(str(d)), str(d))
               if "hold in common" in a.title)
    assert len(inv.commands) == 2
    assert any("--near 12" in c for c in inv.commands)
    assert any("--near 18" in c for c in inv.commands)


def test_clusters_are_used_in_the_suggested_command_when_present(tmp_path):
    d = _kaede_like(tmp_path)
    (d / "profile.tsv").write_text("candidate\tref_resnum\tobs_res\tresolved\n")
    (d / "clusters.tsv").write_text("candidate\tcluster\nc000\tcl0\n")
    inv = next(a for a in advise.advise(advise.inspect(str(d)), str(d))
               if "hold in common" in a.title)
    assert all("--clusters clusters.tsv" in c for c in inv.commands)


def test_missing_clusters_are_flagged_with_the_reason(tmp_path):
    d = _kaede_like(tmp_path)
    item = next(a for a in advise.advise(advise.inspect(str(d)), str(d))
                if "clustering" in a.title)
    assert "expanded clade" in item.why


# -------------------------------------------------------------- other advice
def test_contamination_check_is_advised_and_framed_as_removing_findings(tmp_path):
    d = _kaede_like(tmp_path)
    (d / "profile.tsv").write_text("candidate\tref_resnum\tobs_res\tresolved\n")
    item = next(a for a in advise.advise(advise.inspect(str(d)), str(d))
                if "contradict" in a.title)
    assert "removes false ones" in item.why
    assert item.urgency == 1


def test_a_second_reference_is_suggested_when_there_are_losses(tmp_path):
    d = _kaede_like(tmp_path)
    item = next(a for a in advise.advise(advise.inspect(str(d)), str(d))
                if "second reference" in a.title)
    assert "different questions" in item.why


def test_register_advice_tracks_its_state(tmp_path):
    d = _kaede_like(tmp_path)
    titles = lambda: [a.title for a in advise.advise(advise.inspect(str(d)), str(d))]
    assert any("already known" in t for t in titles())
    (d / "register.json").write_text("{}")
    assert any("Freeze the register" in t for t in titles())
    (d / "FREEZE.md").write_text("# FREEZE")
    assert not any("Freeze the register" in t for t in titles())


def test_detector_outputs_trigger_the_adjudication_reminder(tmp_path):
    d = _kaede_like(tmp_path)
    (d / "profile.tsv").write_text("candidate\tref_resnum\tobs_res\tresolved\n")
    (d / "register.json").write_text("{}")
    (d / "FREEZE.md").write_text("# FREEZE")
    (d / "inv.tsv").write_text("ref_resnum\n")
    assert any("Adjudicate" in a.title
               for a in advise.advise(advise.inspect(str(d)), str(d)))


# ----------------------------------------------------------------- formatting
def test_formatting_leads_with_the_verdict_counts(tmp_path):
    d = _kaede_like(tmp_path)
    text = advise.next_steps(str(d))
    assert "60 candidates" in text and "20 lost" in text
    assert "$ foldtrace profile" in text


def test_brief_shows_two_and_says_how_many_remain(tmp_path):
    d = _kaede_like(tmp_path)
    text = advise.next_steps(str(d), brief=True)
    assert text.count("\n1. ") == 1 and text.count("\n2. ") == 1
    assert "3." not in text
    assert "more: foldtrace next" in text


def test_advice_never_raises_on_a_malformed_calls_file(tmp_path):
    (tmp_path / "calls.tsv").write_text("not\ta\tvalid\theader\n\x00\n")
    text = advise.next_steps(str(tmp_path))
    assert isinstance(text, str) and text
