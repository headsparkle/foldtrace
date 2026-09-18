"""Tests for the unexpected-invariant detector.

Fixtures are CA-only helices so the geometry is controlled: positions 10-20 sit
inside the site zone around the declared centroid, positions 1-5 and 36-40 sit
in the background zone. Candidates share the reference coordinates and differ
only in sequence, so correspondence is exact and the statistics are what is
under test.
"""
from __future__ import annotations

import math

import pytest

from foldtrace import invariant as inv
from foldtrace.invariant import InvariantError
from foldtrace.io import load_structure, parse_sites

N_RES = 40
AA3 = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
       "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
       "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
       "W": "TRP", "Y": "TYR"}


def _helix_coords(n=N_RES, radius=8.0, rise=1.5, per_turn=3.6):
    out = []
    for i in range(n):
        a = 2 * math.pi * i / per_turn
        out.append((radius * math.cos(a), radius * math.sin(a), rise * i))
    return out


def _line_coords(n=N_RES, step=3.8):
    return [(step * i, 0.0, 0.0) for i in range(n)]


def _write_pdb(path, seq, coords):
    lines = []
    for i, (aa, (x, y, z)) in enumerate(zip(seq, coords), start=1):
        lines.append(
            f"ATOM  {i:>5}  CA  {AA3[aa]} A{i:>4}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00 90.00           C")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")
    return str(path)


POOL = "ADEFGIKLMNQSTVY"


def _seq(base="A", **at):
    """A sequence of `base`, with {position(1-based): residue} overrides."""
    s = [base] * N_RES
    for pos, aa in at.items():
        s[int(pos) - 1] = aa
    return "".join(s)


def _varied(seed, **at):
    """A sequence that differs at every position, then the given overrides.

    Real scaffolds vary. An all-alanine background would make every position an
    "invariant" and tell us nothing about the detector.
    """
    s = [POOL[(seed * 7 + 3 * i) % len(POOL)] for i in range(N_RES)]
    for pos, aa in at.items():
        s[int(pos) - 1] = aa
    return "".join(s)


@pytest.fixture
def ref_and_sites(tmp_path):
    # declared chemistry at 12 and 18; centroid sits between them
    seq = _seq("A", **{"12": "H", "18": "R"})
    ref = _write_pdb(tmp_path / "ref.pdb", seq, _helix_coords())
    sites = tmp_path / "sites.tsv"
    sites.write_text("label\tref_resnum\texpected\trole\tkey\taltered\n"
                     "H12\t12\tH\tcatalytic\t1\t\n"
                     "R18\t18\tR\tpositioning\t1\t\n")
    return load_structure(ref), parse_sites(str(sites))


def _candidate(tmp_path, name, seq, coords=None):
    return _write_pdb(tmp_path / f"{name}.pdb", seq, coords or _helix_coords())


# ------------------------------------------------------------------ profiling
def test_profile_covers_every_reference_position(tmp_path, ref_and_sites):
    ref, _ = ref_and_sites
    cand = load_structure(_candidate(tmp_path, "c1", _seq("A", **{"12": "H", "18": "R"})))
    obs, tm, fold_ok = inv.profile_candidate(ref, cand, "c1")
    assert fold_ok and tm > 0.99
    assert len(obs) == N_RES
    assert {o.ref_resnum for o in obs} == set(range(1, N_RES + 1))
    assert all(o.resolved for o in obs)


def test_below_gate_returns_rows_not_silence(tmp_path, ref_and_sites):
    """Rows out must not depend on the answer."""
    ref, _ = ref_and_sites
    cand = load_structure(_candidate(tmp_path, "far", _seq("A"), _line_coords()))
    obs, tm, fold_ok = inv.profile_candidate(ref, cand, "far")
    assert not fold_ok and tm < 0.5
    assert len(obs) == N_RES and not any(o.resolved for o in obs)


def test_profile_round_trip(tmp_path, ref_and_sites):
    ref, _ = ref_and_sites
    cand = load_structure(_candidate(tmp_path, "c1", _seq("A", **{"12": "H"})))
    obs, _, _ = inv.profile_candidate(ref, cand, "c1")
    p = tmp_path / "profile.tsv"
    inv.write_profile(obs, str(p))
    back = inv.read_profile(str(p))
    assert len(back) == len(obs)
    assert back[11].obs_res == "H" and back[11].resolved


# ------------------------------------------------------------------- geometry
def test_zones_follow_distance_from_the_declared_centroid(tmp_path, ref_and_sites):
    ref, sites = ref_and_sites
    d = inv.distances_to_site(ref, inv.site_centroid(ref, sites))
    assert inv.zone_of(d[15], 12.0, 20.0) == "site"
    assert inv.zone_of(d[1], 12.0, 20.0) == "background"
    assert inv.zone_of(d[40], 12.0, 20.0) == "background"


def test_centroid_requires_a_declared_position_present(tmp_path, ref_and_sites):
    ref, _ = ref_and_sites
    from foldtrace.io import Site
    bogus = [Site("X999", 999, frozenset("H"), "catalytic", True)]
    with pytest.raises(InvariantError, match="no declared site position"):
        inv.site_centroid(ref, bogus)


# ------------------------------------------------------------------- analysis
def _profile_set(tmp_path, ref, specs):
    """specs: {name: (sequence, group)}. Returns (profile rows, groups)."""
    rows, groups = [], {}
    for name, (seq, group) in specs.items():
        cand = load_structure(_candidate(tmp_path, name, seq))
        obs, _, _ = inv.profile_candidate(ref, cand, name)
        rows.extend(obs)
        groups[name] = group
    return rows, groups


def _lost_set(n=10, network=True, conserved_background=False):
    """Members that lost the declared His but keep an undeclared Trp at 15."""
    specs = {}
    for i in range(n):
        over = {"12": "K", "18": "R"}                      # declared chemistry gone
        over["15"] = "W" if network else POOL[i % len(POOL)]
        if conserved_background:
            over["2"] = over["3"] = "C"                    # invariant, but far from the site
        specs[f"lost{i}"] = (_varied(i, **over), "lost")
    return specs


def test_finds_the_undeclared_invariant(tmp_path, ref_and_sites):
    ref, sites = ref_and_sites
    specs = _lost_set()
    profile, groups = _profile_set(tmp_path, ref, specs)
    stats, summary = inv.analyse(profile, ref, sites, groups, compare=None)
    found = inv.rank(stats)
    assert found, "expected at least one finding"
    top = found[0]
    assert top.ref_resnum == 15 and top.modal == "W"
    assert top.conservation == 1.0
    assert "independent observations" in top.reason
    assert summary["focus_n"] == 10


def test_variable_near_site_position_is_not_reported(tmp_path, ref_and_sites):
    ref, sites = ref_and_sites
    profile, groups = _profile_set(tmp_path, ref, _lost_set())
    stats, _ = inv.analyse(profile, ref, sites, groups, compare=None)
    by_num = {s.ref_resnum: s for s in stats}
    variable = [n for n, st in by_num.items()
                if st.zone == "site" and not st.declared and st.conservation < 0.9]
    assert variable, "the fixture must contain varied site-zone positions"
    for n in variable:
        assert by_num[n].score == 0.0
        assert "conservation" in by_num[n].reason


def test_declared_positions_are_excluded_with_a_reason(tmp_path, ref_and_sites):
    ref, sites = ref_and_sites
    profile, groups = _profile_set(tmp_path, ref, _lost_set())
    stats, _ = inv.analyse(profile, ref, sites, groups, compare=None)
    by_num = {s.ref_resnum: s for s in stats}
    for n in (12, 18):
        assert by_num[n].score == 0.0
        assert "declared position" in by_num[n].reason
    assert 12 not in {s.ref_resnum for s in inv.rank(stats)}


def test_background_positions_are_excluded(tmp_path, ref_and_sites):
    """An invariant far from the site is scaffold, not a candidate device."""
    ref, sites = ref_and_sites
    specs = _lost_set(conserved_background=True)
    profile, groups = _profile_set(tmp_path, ref, specs)
    stats, _ = inv.analyse(profile, ref, sites, groups, compare=None)
    by_num = {s.ref_resnum: s for s in stats}
    assert by_num[2].conservation == 1.0 and by_num[2].modal == "C"
    assert by_num[2].score == 0.0 and "background zone" in by_num[2].reason


def test_contrast_against_the_comparison_group(tmp_path, ref_and_sites):
    ref, sites = ref_and_sites
    specs = _lost_set()
    pool = "ADEFGIKLMNQ"
    for i in range(6):  # retained members vary at 15
        specs[f"kept{i}"] = (_varied(100 + i, **{"12": "H", "18": "R", "15": pool[i]}),
                             "retained")
    profile, groups = _profile_set(tmp_path, ref, specs)
    stats, summary = inv.analyse(profile, ref, sites, groups, focus="lost", compare="retained")
    top = inv.rank(stats)[0]
    assert top.ref_resnum == 15
    assert top.contrast is not None and top.contrast > 0.5
    assert "more conserved here" in top.reason
    assert summary["compare_n"] == 6


def test_percentile_reports_an_uninformative_background(tmp_path, ref_and_sites):
    """If the background is also invariant, the percentile says so."""
    ref, sites = ref_and_sites
    profile, groups = _profile_set(tmp_path, ref, _lost_set(conserved_background=True))
    stats, _ = inv.analyse(profile, ref, sites, groups, compare=None)
    by_num = {s.ref_resnum: s for s in stats}
    # two background positions are also invariant, so the percentile is not 1.0
    assert by_num[15].bg_percentile is not None and by_num[15].bg_percentile < 1.0


def test_empty_focus_group_is_an_error(tmp_path, ref_and_sites):
    ref, sites = ref_and_sites
    profile, groups = _profile_set(tmp_path, ref, _lost_set())
    with pytest.raises(InvariantError, match="focus group"):
        inv.analyse(profile, ref, sites, groups, focus="altered")


# ---------------------------------------------------------- non-independence
def test_clustering_collapses_an_expanded_clade(tmp_path, ref_and_sites):
    """Eight near-identical members must not vote eight times.

    Position 15 is Trp in eight members of one cluster and something else in two
    others. Per protein that reads as 80% conserved; per sequence cluster it is
    one cluster of three, which is the honest number.
    """
    ref, sites = ref_and_sites
    specs = {}
    for i in range(8):
        specs[f"clade{i}"] = (_varied(i, **{"12": "K", "18": "R", "15": "W"}), "lost")
    specs["out1"] = (_varied(50, **{"12": "K", "18": "R", "15": "D"}), "lost")
    specs["out2"] = (_varied(60, **{"12": "K", "18": "R", "15": "E"}), "lost")
    profile, groups = _profile_set(tmp_path, ref, specs)

    raw, _ = inv.analyse(profile, ref, sites, groups, compare=None, min_conservation=0.75)
    assert {s.ref_resnum for s in inv.rank(raw)} == {15}

    clusters = {f"clade{i}": "c1" for i in range(8)}
    clusters.update({"out1": "c2", "out2": "c3"})
    clustered, summary = inv.analyse(profile, ref, sites, groups, compare=None,
                                     clusters=clusters, min_conservation=0.75)
    by_num = {s.ref_resnum: s for s in clustered}
    assert by_num[15].n == 3
    assert by_num[15].conservation == pytest.approx(1 / 3)
    assert by_num[15].score == 0.0
    assert summary["clusters"] == 3


def test_low_coverage_position_is_not_judged(tmp_path, ref_and_sites):
    ref, sites = ref_and_sites
    specs = _lost_set(n=8)
    profile, groups = _profile_set(tmp_path, ref, specs)
    # simulate an unresolvable position in most members
    for o in profile:
        if o.ref_resnum == 15 and o.candidate not in ("lost0", "lost1"):
            o.resolved = False
    stats, _ = inv.analyse(profile, ref, sites, groups, compare=None)
    by_num = {s.ref_resnum: s for s in stats}
    assert by_num[15].score == 0.0 and "coverage" in by_num[15].reason


# ------------------------------------------------------------------- outputs
def test_groups_load_from_a_calls_table(tmp_path):
    p = tmp_path / "calls.tsv"
    p.write_text("candidate_id\ttmalign_tm_norm_ref\tverdict\n"
                 "P1\t0.9\tlost\nP2\t0.9\tRETAINED\n")
    g = inv.load_groups(str(p))
    assert g == {"P1": "lost", "P2": "retained"}


def test_clusters_require_pairs(tmp_path):
    p = tmp_path / "cl.tsv"
    p.write_text("candidate\tnot_cluster\nP1\tx\n")
    with pytest.raises(InvariantError):
        inv.load_clusters(str(p))


def test_findings_write_and_emit(tmp_path, ref_and_sites):
    ref, sites = ref_and_sites
    profile, groups = _profile_set(tmp_path, ref, _lost_set())
    stats, _ = inv.analyse(profile, ref, sites, groups, compare=None)
    found = inv.rank(stats)
    out = tmp_path / "found.tsv"
    inv.write_findings(found, str(out))
    lines = out.read_text().strip().split("\n")
    assert lines[0].split("\t") == inv.FINDING_COLUMNS
    feat = inv.emit_features(found)
    assert feat and "--surfaced-by unexpected_invariant" in feat[0]
    assert "W15" in feat[0]


def test_summary_reports_site_zone_saturation(tmp_path, ref_and_sites):
    """A flooded result must say so rather than look like a discovery."""
    ref, sites = ref_and_sites
    varied, groups = _profile_set(tmp_path, ref, _lost_set())
    _, lean = inv.analyse(varied, ref, sites, groups, compare=None)
    assert float(lean["site_fraction_reported"]) < 0.3

    flat = {f"flat{i}": (_seq("A", **{"12": "K", "18": "R"}), "lost") for i in range(10)}
    profile, groups = _profile_set(tmp_path, ref, flat)
    _, sat = inv.analyse(profile, ref, sites, groups, compare=None)
    assert float(sat["site_fraction_reported"]) > 0.9
