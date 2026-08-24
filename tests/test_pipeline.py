"""End-to-end tests for ``foldtrace run`` (stage 1 table -> stage 2 mapping), offline."""
import os

from foldtrace import cli
from foldtrace.io import parse_sites
from foldtrace.pipeline import run, write_results, result_columns

HERE = os.path.dirname(__file__)
TIR = os.path.join(HERE, "..", "examples", "tir")
REF = os.path.join(TIR, "reference", "SARM1_TIR_6O0R_A.pdb")
SITES = os.path.join(TIR, "catalytic_sites.tsv")
HITS = os.path.join(TIR, "example_search_table.tsv")
CANDS = os.path.join(TIR, "candidates")


def _run():
    return run(REF, SITES, hits_table=HITS, candidates_dir=CANDS)


def test_end_to_end_verdicts_match_known_states():
    sites, records = _run()
    verdict = {r.hit.target: (r.result.verdict if r.result else "unresolved") for r in records}
    assert verdict["A0A953THP3"] == "retained"
    assert verdict["A0A7Y7TLD0"] == "retained"
    assert verdict["Q9FHM1"] == "retained"
    assert verdict["A0A933QTC5"] == "lost"


def test_run_merges_foldseek_metrics_and_mapping():
    sites, records = _run()
    rec = next(r for r in records if r.hit.target == "A0A933QTC5")
    # foldseek search metric carried through
    assert rec.hit.alntmscore is not None
    # mapping produced a TM-align score and a reason for the lost verdict
    assert rec.result.tm_norm_ref > 0.5
    assert "not in expected" in rec.result.state_reason


def test_unresolved_structure_not_found_is_reported_not_dropped(tmp_path):
    # a hits table naming a target with no structure file must yield an unresolved row, not a drop
    hits = tmp_path / "hits.tsv"
    hits.write_text("rank\tquery\ttarget\tprob\tevalue\tbits\tfident\talnlen\tqcov\ttcov\talntmscore\n"
                    "1\tQ\tGHOST999\tNA\tNA\tNA\t0.2\t100\t0.8\t0.8\t0.7\n")
    sites, records = run(REF, SITES, hits_table=str(hits), candidates_dir=CANDS)
    assert len(records) == 1
    assert records[0].result is None
    assert records[0].reason == "structure_not_found"


def test_run_output_columns_and_cli(tmp_path):
    out = tmp_path / "run.tsv"
    rc = cli.main(["run", "--reference", REF, "--sites", SITES, "--hits", HITS,
                   "--candidates-dir", CANDS, "--out", str(out)])
    assert rc == 0
    lines = out.read_text().strip().splitlines()
    sites = parse_sites(SITES)
    assert lines[0].split("\t") == result_columns(sites)
    assert len(lines) == 5  # header + 4 candidates
    # state_reason column present and populated for the lost candidate
    assert any("not in expected" in ln for ln in lines)


def test_run_candidates_mode_without_search(tmp_path):
    # stage-1 skipped: map explicit candidate files directly
    out = tmp_path / "run2.tsv"
    rc = cli.main(["run", "--reference", REF, "--sites", SITES,
                   "--candidates", os.path.join(CANDS, "A0A933QTC5.pdb"), "--out", str(out)])
    assert rc == 0
    assert "lost" in out.read_text()


def test_run_requires_a_stage1_source():
    import pytest
    with pytest.raises(ValueError):
        run(REF, SITES)
