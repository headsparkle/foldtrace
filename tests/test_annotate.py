"""Tests for the fetch layer.

The network half cannot be tested here and is confined to one function. The
parsing and provenance half is pure, and is tested against payloads shaped like
the documented responses -- with the shapes deliberately varied, because the
point of these tests is that an unexpected shape degrades a row rather than
ending a run.
"""
from __future__ import annotations

import json

import pytest

from foldtrace import annotate as ann
from foldtrace.annotate import AnnotateError


UNIPROT_OK = {
    "uniProtkbId": "HEPN_SHEON",
    "entryType": "UniProtKB reviewed (Swiss-Prot)",
    "proteinDescription": {"recommendedName": {"fullName": {"value": "HEPN domain protein"}}},
    "organism": {"scientificName": "Shewanella oneidensis", "taxonId": 70863},
    "sequence": {"length": 133},
    "uniProtKBCrossReferences": [
        {"database": "Pfam", "id": "PF01934"},
        {"database": "InterPro", "id": "IPR002711"},
        {"database": "Pfam", "id": "PF13157"},
    ],
}

UNIPARC_OK = {
    "results": [{
        "sequence": {"length": 241},
        "uniParcCrossReferences": [
            {"database": "EMBL", "proteinName": "phosphotriesterase-related protein",
             "organism": {"scientificName": "Sulfolobus tengchongensis", "taxonId": 207809}},
        ],
    }],
}


# ------------------------------------------------------------- uniprot parser
def test_parses_a_complete_record():
    a = ann.parse_uniprot("Q8EE45", UNIPROT_OK)
    assert a.name == "HEPN domain protein"
    assert a.organism == "Shewanella oneidensis" and a.taxon_id == "70863"
    assert a.length == "133" and a.status == "ok"
    assert a.pfam == "PF01934;PF13157", "multiple Pfam accessions are kept, in order"


def test_submission_name_is_used_when_there_is_no_recommended_name():
    payload = dict(UNIPROT_OK)
    payload["proteinDescription"] = {
        "submissionNames": [{"fullName": {"value": "Uncharacterized protein"}}]}
    assert ann.parse_uniprot("X", payload).name == "Uncharacterized protein"


def test_falls_back_to_the_entry_id_when_no_name_is_present():
    payload = {k: v for k, v in UNIPROT_OK.items() if k != "proteinDescription"}
    assert ann.parse_uniprot("X", payload).name == "HEPN_SHEON"


def test_a_record_with_no_pfam_is_not_an_error():
    payload = {k: v for k, v in UNIPROT_OK.items() if k != "uniProtKBCrossReferences"}
    a = ann.parse_uniprot("X", payload)
    assert a.pfam == "" and a.status == "ok"


def test_missing_fields_degrade_the_row_rather_than_raising():
    """A shape change must cost one field, not the run."""
    a = ann.parse_uniprot("X", {"organism": {}})
    assert a.status == "ok"
    assert a.name == "" and a.organism == "" and a.length == "" and a.pfam == ""


def test_a_non_dict_payload_is_a_parse_error_not_a_crash():
    assert ann.parse_uniprot("X", ["unexpected"]).status == "parse_error"


def test_an_inactive_entry_is_marked_withdrawn():
    payload = dict(UNIPROT_OK)
    payload["entryType"] = "Inactive"
    assert ann.parse_uniprot("X", payload).status == "withdrawn"


# ------------------------------------------------------------- uniparc parser
def test_uniparc_recovers_organism_and_name():
    a = ann.parse_uniparc("UPI0001", UNIPARC_OK)
    assert a.status == "withdrawn_recovered"
    assert a.organism == "Sulfolobus tengchongensis"
    assert a.name == "phosphotriesterase-related protein"
    assert a.length == "241" and a.source == "uniparc"


def test_uniparc_with_no_results_is_not_found():
    assert ann.parse_uniparc("X", {"results": []}).status == "not_found"


# --------------------------------------------------------------- clan parser
@pytest.mark.parametrize("payload,expected", [
    ({"metadata": {"set_info": {"accession": "CL0291"}}}, "CL0291"),
    ({"metadata": {"sets": [{"accession": "CL0291"}]}}, "CL0291"),
    ({"metadata": {"sets": ["CL0291"]}}, "CL0291"),
    ({"metadata": {}}, ""),
    ({}, ""),
    ("not a dict", ""),
])
def test_clan_parser_tolerates_shape_variation(payload, expected):
    assert ann.parse_interpro_clan(payload) == expected


# ------------------------------------------------------------------ pipeline
def _fetcher(table):
    """An injected fetcher: {url_substring: (payload, status)}."""
    def fetch(url, **kw):
        for key, value in table.items():
            if key in url:
                return value
        return None, "not_found"
    return fetch


def test_every_accession_gets_a_row_whatever_happens():
    """Rows out must equal accessions in; a failed lookup is not an absence."""
    fetch = _fetcher({"Q8EE45": (UNIPROT_OK, "ok"),
                      "PF01934": ({"metadata": {"set_info": {"accession": "CL0291"}}}, "ok")})
    rows, man = ann.annotate(["Q8EE45", "MISSING1", "MISSING2"], fetcher=fetch)
    assert [r.target for r in rows] == ["Q8EE45", "MISSING1", "MISSING2"]
    assert man.requested == 3
    assert man.counts["ok"] == 1


def test_clan_is_attached_from_interpro():
    fetch = _fetcher({"Q8EE45": (UNIPROT_OK, "ok"),
                      "PF01934": ({"metadata": {"set_info": {"accession": "CL0291"}}}, "ok"),
                      "PF13157": ({"metadata": {"set_info": {"accession": "CL0291"}}}, "ok")})
    rows, _ = ann.annotate(["Q8EE45"], fetcher=fetch)
    assert rows[0].clan == "CL0291", "duplicate clans collapse to one"


def test_a_withdrawn_record_is_recovered_from_uniparc():
    withdrawn = dict(UNIPROT_OK); withdrawn["entryType"] = "Inactive"
    fetch = _fetcher({"uniprotkb": (withdrawn, "ok"), "uniparc": (UNIPARC_OK, "ok"),
                      "PF01934": ({"metadata": {}}, "ok"), "PF13157": ({"metadata": {}}, "ok")})
    rows, man = ann.annotate(["A0AAX4L3M9"], fetcher=fetch)
    assert rows[0].status == "withdrawn_recovered"
    assert rows[0].organism == "Sulfolobus tengchongensis"
    assert rows[0].pfam == "PF01934;PF13157", "Pfam from the inactive record is carried over"
    assert man.counts["withdrawn_recovered"] == 1


def test_a_record_missing_everywhere_stays_not_found():
    rows, man = ann.annotate(["NOPE"], fetcher=_fetcher({}))
    assert rows[0].status == "not_found" and man.counts["not_found"] == 1


def test_http_errors_are_recorded_not_raised():
    rows, man = ann.annotate(["X"], fetcher=_fetcher({"X": (None, "http_error:503")}))
    assert rows[0].status == "http_error:503"
    assert "http_error:503" in man.counts


def test_no_accessions_is_an_error():
    with pytest.raises(AnnotateError):
        ann.annotate([], fetcher=_fetcher({}))


# ------------------------------------------------------------------- outputs
def test_manifest_records_endpoints_counts_and_hash(tmp_path):
    fetch = _fetcher({"Q8EE45": (UNIPROT_OK, "ok"), "PF": ({"metadata": {}}, "ok")})
    rows, man = ann.annotate(["Q8EE45", "MISSING"], fetcher=fetch)
    out = str(tmp_path / "annotations.tsv")
    man.output = out
    man.output_sha256 = ann.write_annotations(rows, out)
    mpath = str(tmp_path / "manifest.json")
    ann.write_manifest(man, mpath)

    text = open(out).read().strip().split("\n")
    assert text[0].split("\t") == ann.ANNOTATION_COLUMNS
    assert len(text) == 3

    m = json.load(open(mpath))
    assert m["requested"] == 2
    assert any("rest.uniprot.org" in e for e in m["endpoints"])
    assert len(m["output_sha256"]) == 64
    assert "never a biological absence" in m["note"]


def test_the_output_hash_changes_with_the_content(tmp_path):
    a = [ann.Annotation(target="X", name="one", status="ok")]
    b = [ann.Annotation(target="X", name="two", status="ok")]
    h1 = ann.write_annotations(a, str(tmp_path / "a.tsv"))
    h2 = ann.write_annotations(b, str(tmp_path / "b.tsv"))
    assert h1 != h2


def test_cache_removes_the_second_request(tmp_path):
    calls = []

    def counting(url, **kw):
        calls.append(url)
        return (UNIPROT_OK, "ok") if "uniprotkb" in url else ({"metadata": {}}, "ok")

    cache = str(tmp_path / "cache")
    ann.annotate(["Q8EE45"], cache_dir=cache, fetcher=counting)
    first = len(calls)
    ann.annotate(["Q8EE45"], cache_dir=cache, fetcher=counting)
    assert len(calls) == first, "a cached re-run makes no requests"


def test_accession_list_reading(tmp_path):
    p = tmp_path / "acc.txt"
    p.write_text("target\nQ8EE45\n# a comment\n\nP12345\nQ8EE45\n")
    assert ann.read_accessions(str(p)) == ["Q8EE45", "P12345"]


def test_accession_list_accepts_a_tsv_first_column(tmp_path):
    p = tmp_path / "calls.tsv"
    p.write_text("candidate\tverdict\nA0A852UQU6\tlost\nA0A7Y5GXD8\tretained\n")
    assert ann.read_accessions(str(p)) == ["A0A852UQU6", "A0A7Y5GXD8"]
