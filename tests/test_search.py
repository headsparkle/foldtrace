"""Tests for stage 1: Foldseek parsing, filtering/ranking, and the search CLI (offline)."""
import os

from foldtrace.search import (
    parse_foldseek, filter_and_rank, write_hits_tsv, HIT_COLUMNS,
    FOLDSEEK_FORMAT_FIELDS,
)
from foldtrace import cli

HERE = os.path.dirname(__file__)
SAMPLE = os.path.join(HERE, "data", "sample_foldseek.m8")


def test_parse_foldseek_reads_columns_and_normalizes_percent_identity():
    hits = parse_foldseek(SAMPLE)
    assert len(hits) == 3
    a = {h.target: h for h in hits}
    assert a["targetA"].prob == 0.99 and a["targetA"].alntmscore == 0.88
    # fident given as 25.0 (percent) must be normalised to 0.25
    assert abs(a["targetB"].fident - 0.25) < 1e-9
    assert a["targetA"].fident == 0.42


def test_parse_foldseek_strips_webserver_description(tmp_path):
    # The Foldseek webserver returns afdb50 ids with a trailing free-text
    # description; parse_foldseek must keep only the id token (reported bug:
    # otherwise the description pollutes the id and everything zeroes out).
    fields = FOLDSEEK_FORMAT_FIELDS
    vals = {
        "query": "queryX chain description",
        "target": "AF-A0A8C0MM32-F1-model_v6 Carbonic anhydrase",
        "prob": "0.9", "evalue": "1e-9", "bits": "500", "fident": "0.98",
        "qcov": "0.9", "tcov": "0.9", "alntmscore": "0.85",
    }
    line = "\t".join(vals.get(f, "0") for f in fields)
    p = tmp_path / "webserver.m8"
    p.write_text(line + "\n")
    h = parse_foldseek(str(p))[0]
    assert h.target == "AF-A0A8C0MM32-F1-model_v6"
    assert h.query == "queryX"
    assert " " not in h.target and " " not in h.query


def test_filter_and_rank_thresholds_and_order():
    hits = parse_foldseek(SAMPLE)
    # permissive: keep all, ranked by prob desc
    ranked = filter_and_rank(hits)
    assert [h.target for h in ranked] == ["targetA", "targetB", "targetC"]
    assert [h.rank for h in ranked] == [1, 2, 3]
    # probability threshold drops the low-prob hit
    kept = filter_and_rank(hits, min_prob=0.8)
    assert [h.target for h in kept] == ["targetA"]
    # coverage + identity thresholds
    assert {h.target for h in filter_and_rank(hits, min_coverage=0.5)} == {"targetA", "targetB"}
    assert {h.target for h in filter_and_rank(hits, min_identity=0.3)} == {"targetA"}


def test_max_hits_caps_results():
    hits = parse_foldseek(SAMPLE)
    assert len(filter_and_rank(hits, max_hits=2)) == 2


def test_write_hits_roundtrip(tmp_path):
    hits = filter_and_rank(parse_foldseek(SAMPLE))
    out = tmp_path / "hits.tsv"
    write_hits_tsv(hits, str(out))
    lines = out.read_text().strip().splitlines()
    assert lines[0].split("\t") == HIT_COLUMNS
    assert len(lines) == 4  # header + 3


def test_search_cli_from_foldseek(tmp_path):
    out = tmp_path / "search.tsv"
    rc = cli.main(["search", "--from-foldseek", SAMPLE, "--out", str(out)])
    assert rc == 0
    assert out.exists()
    assert "targetA" in out.read_text()


def test_search_cli_requires_database_or_foldseek(capsys):
    # no database and no --from-foldseek -> graceful error, exit 2
    rc = cli.main(["search", "--query", "x.pdb"])
    assert rc == 2
