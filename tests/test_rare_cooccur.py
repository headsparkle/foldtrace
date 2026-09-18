"""Tests for the rare-residue and co-occurrence detectors.

Fixtures mirror the corpus cases: His287 present in a handful of beta-propeller
members on a scaffold that is otherwise Gly, and an FGE domain that keeps
turning up attached to TIR hits.
"""
from __future__ import annotations

import math

import pytest

from foldtrace import rare, cooccur as co
from foldtrace.cooccur import CooccurError
from foldtrace.invariant import PositionObs, InvariantError
from foldtrace.io import load_structure, parse_sites

from test_invariant import _write_pdb, _helix_coords, _varied, _seq, N_RES  # noqa: F401


@pytest.fixture
def ref_and_sites(tmp_path):
    seq = _seq("A", **{"12": "H", "18": "R"})
    ref = _write_pdb(tmp_path / "ref.pdb", seq, _helix_coords())
    sites = tmp_path / "sites.tsv"
    sites.write_text("label\tref_resnum\texpected\trole\tkey\taltered\n"
                     "H12\t12\tH\tcatalytic\t1\t\n"
                     "R18\t18\tR\tpositioning\t1\t\n")
    return load_structure(ref), parse_sites(str(sites))


def _obs(candidates_at_position):
    """Build a profile: {ref_resnum: {candidate: residue}}."""
    rows = []
    for rn, per_cand in candidates_at_position.items():
        for cand, res in per_cand.items():
            rows.append(PositionObs(cand, rn, res, rn, 0.5, True))
    return rows


# ------------------------------------------------------------ rare residues
def test_rare_residue_on_a_common_scaffold(ref_and_sites):
    """His287: the marking residue is the one almost nobody has."""
    ref, sites = ref_and_sites
    members = [f"m{i:02d}" for i in range(40)]
    at15 = {m: ("H" if i < 4 else "G") for i, m in enumerate(members)}
    f, summary = rare.analyse_rare(_obs({15: at15}), ref, sites)
    assert len(f) == 1
    top = f[0]
    assert top.ref_resnum == 15 and top.residue == "H"
    assert top.n_carriers == 4 and top.modal == "G"
    assert top.carriers == ["m00", "m01", "m02", "m03"]
    assert "otherwise G in 90% of members" in top.reason
    assert summary["rare_calls"] == 1


def test_a_common_residue_is_not_rare(ref_and_sites):
    ref, sites = ref_and_sites
    members = [f"m{i:02d}" for i in range(20)]
    at15 = {m: ("H" if i < 9 else "G") for i, m in enumerate(members)}
    f, _ = rare.analyse_rare(_obs({15: at15}), ref, sites)
    assert f == []


def test_a_single_carrier_is_not_a_finding(ref_and_sites):
    ref, sites = ref_and_sites
    at15 = {f"m{i:02d}": ("H" if i == 0 else "G") for i in range(40)}
    assert rare.analyse_rare(_obs({15: at15}), ref, sites)[0] == []


def test_non_chemistry_residues_are_ignored(ref_and_sites):
    """A rare leucine is not a candidate active site."""
    ref, sites = ref_and_sites
    at15 = {f"m{i:02d}": ("L" if i < 4 else "G") for i in range(40)}
    assert rare.analyse_rare(_obs({15: at15}), ref, sites)[0] == []


def test_declared_positions_are_skipped_unless_asked(ref_and_sites):
    ref, sites = ref_and_sites
    at12 = {f"m{i:02d}": ("C" if i < 3 else "G") for i in range(40)}
    assert rare.analyse_rare(_obs({12: at12}), ref, sites)[0] == []
    inc, _ = rare.analyse_rare(_obs({12: at12}), ref, sites, include_declared=True)
    assert len(inc) == 1 and inc[0].ref_resnum == 12


def test_background_zone_positions_are_skipped(ref_and_sites):
    ref, sites = ref_and_sites
    at2 = {f"m{i:02d}": ("H" if i < 4 else "G") for i in range(40)}
    assert rare.analyse_rare(_obs({2: at2}), ref, sites)[0] == []


def test_clustering_stops_one_clade_faking_a_recurrent_gain(ref_and_sites):
    ref, sites = ref_and_sites
    members = [f"m{i:02d}" for i in range(40)]
    at15 = {m: ("H" if i < 4 else "G") for i, m in enumerate(members)}
    clusters = {m: ("one" if i < 4 else f"c{i}") for i, m in enumerate(members)}
    f, _ = rare.analyse_rare(_obs({15: at15}), ref, sites, clusters=clusters, min_count=2)
    assert f == [], "four members of a single cluster are one observation, not four"


def test_a_diffuse_position_is_not_tested(ref_and_sites):
    """Where every residue is rare, rarity means nothing."""
    ref, sites = ref_and_sites
    pool = "ADEFGIKLMNQSTVWY"
    members = [f"m{i:02d}" for i in range(40)]
    at15 = {m: ("H" if i < 4 else pool[i % len(pool)]) for i, m in enumerate(members)}
    f, _ = rare.analyse_rare(_obs({15: at15}), ref, sites)
    assert f == []


def test_modal_dominance_raises_the_score(ref_and_sites):
    ref, sites = ref_and_sites
    members = [f"m{i:02d}" for i in range(40)]
    clean = {m: ("H" if i < 4 else "G") for i, m in enumerate(members)}
    mixed = {m: ("H" if i < 4 else "S" if i < 12 else "G") for i, m in enumerate(members)}
    a, _ = rare.analyse_rare(_obs({15: clean}), ref, sites)
    b, _ = rare.analyse_rare(_obs({16: mixed}), ref, sites)
    assert a[0].score > b[0].score
    assert a[0].modal_frac > b[0].modal_frac


def test_a_minority_inside_a_cluster_is_not_erased(ref_and_sites):
    """Clustering must stop a clade voting twice without silencing rare residues.

    Each carrier is the only member of its own three-protein cluster to carry
    the residue. Collapsing clusters to a modal residue would delete all four.
    """
    ref, sites = ref_and_sites
    members = [f"m{i:02d}" for i in range(36)]
    at15 = {m: ("C" if i % 9 == 0 else "G") for i, m in enumerate(members)}
    clusters = {m: f"cl{i // 3}" for i, m in enumerate(members)}
    f, _ = rare.analyse_rare(_obs({15: at15}), ref, sites, clusters=clusters,
                             max_fraction=0.40)
    assert len(f) == 1
    assert f[0].residue == "C"
    assert f[0].n_carriers == 4, "four distinct clusters carry it"
    assert len(f[0].carriers) == 4


def test_one_clade_still_counts_once(ref_and_sites):
    ref, sites = ref_and_sites
    members = [f"m{i:02d}" for i in range(36)]
    at15 = {m: ("C" if i < 4 else "G") for i, m in enumerate(members)}
    clusters = {m: ("one" if i < 4 else f"cl{i}") for i, m in enumerate(members)}
    f, _ = rare.analyse_rare(_obs({15: at15}), ref, sites, clusters=clusters, min_count=2)
    assert f == [], "four members of one cluster are one observation"


def test_rarer_residues_rank_higher(ref_and_sites):
    ref, sites = ref_and_sites
    members = [f"m{i:02d}" for i in range(100)]
    prof = _obs({
        15: {m: ("H" if i < 3 else "G") for i, m in enumerate(members)},
        16: {m: ("C" if i < 9 else "G") for i, m in enumerate(members)},
    })
    f, _ = rare.analyse_rare(prof, ref, sites)
    assert [x.ref_resnum for x in f] == [15, 16]


# ------------------------------------------------------------- combinations
def test_members_carrying_several_rare_residues(ref_and_sites):
    """15 carriers, 9 with the complete site: the combination is the short list."""
    ref, sites = ref_and_sites
    members = [f"m{i:02d}" for i in range(60)]
    prof = _obs({
        15: {m: ("H" if i < 9 else "G") for i, m in enumerate(members)},
        16: {m: ("D" if i < 5 else "G") for i, m in enumerate(members)},
        14: {m: ("C" if i < 5 else "G") for i, m in enumerate(members)},
    })
    f, _ = rare.analyse_rare(prof, ref, sites, max_fraction=0.20)
    combos = rare.combinations(f)
    assert combos
    best = combos[0]
    assert best.positions == [14, 15, 16]
    assert len(best.members) == 5
    assert "complete site" in best.reason


def test_no_combination_when_carriers_do_not_overlap(ref_and_sites):
    ref, sites = ref_and_sites
    members = [f"m{i:02d}" for i in range(60)]
    prof = _obs({
        15: {m: ("H" if i < 4 else "G") for i, m in enumerate(members)},
        16: {m: ("D" if 30 <= i < 34 else "G") for i, m in enumerate(members)},
    })
    f, _ = rare.analyse_rare(prof, ref, sites)
    assert len(f) == 2 and rare.combinations(f) == []


def test_rare_emits_features(ref_and_sites):
    ref, sites = ref_and_sites
    at15 = {f"m{i:02d}": ("H" if i < 4 else "G") for i in range(40)}
    f, _ = rare.analyse_rare(_obs({15: at15}), ref, sites)
    lines = rare.emit_features(f, [])
    assert lines and "--surfaced-by rare_residue" in lines[0]


# ------------------------------------------------------------------ fisher
def test_fisher_matches_a_hand_computed_table():
    orat, p = co.fisher_exact(10, 2, 2, 10)
    assert p < 0.01 and orat > 10


def test_one_sided_is_the_default_and_is_more_powerful():
    """The directional --min-odds filter makes a two-sided penalty inconsistent."""
    _, one = co.fisher_exact(8, 16, 1, 16)
    _, two = co.fisher_exact(8, 16, 1, 16, alternative="two-sided")
    assert one < two
    assert co.fisher_exact(8, 16, 1, 16)[1] == one, "greater is the default"


def test_depletion_is_not_reported_as_enrichment():
    _, p = co.fisher_exact(1, 23, 10, 7)
    assert p > 0.9, "a depleted feature must not get a small one-sided p"


def test_fisher_on_an_uninformative_table():
    _, p = co.fisher_exact(5, 5, 5, 5, alternative="two-sided")
    assert p == pytest.approx(1.0, abs=1e-9)


def test_fisher_handles_a_zero_cell():
    orat, p = co.fisher_exact(8, 0, 0, 8)
    assert math.isfinite(orat) and p < 0.01


def test_benjamini_hochberg_is_monotone():
    q = co.benjamini_hochberg([0.001, 0.01, 0.2, 0.5])
    assert q == sorted(q) and all(a <= b for a, b in zip(q, [0.004, 0.02, 0.267, 0.5]))


# ------------------------------------------------------------------ cooccur
def _features(hits_with_fge=9, hits=100, bg_with_fge=1, bg=400, genus_spread=True):
    feats: dict[str, list[tuple[str, str]]] = {}
    focus = set()
    for i in range(hits):
        t = f"hit{i:03d}"
        focus.add(t)
        feats[t] = [("domain", "TIR")]
        if i < hits_with_fge:
            feats[t].append(("domain", "FGE/SUMF1"))
    for i in range(bg):
        t = f"bg{i:03d}"
        feats[t] = [("domain", "other")]
        if i < bg_with_fge:
            feats[t].append(("domain", "FGE/SUMF1"))
    spread = {}
    for i in range(hits):
        spread[f"hit{i:03d}"] = f"genus{i % 11}" if genus_spread else "genusA"
    return feats, focus, spread


def test_cooccurring_domain_is_found():
    feats, focus, spread = _features()
    f, summary = co.analyse(feats, focus, spread_by=spread)
    names = [x.feature for x in f]
    assert "FGE/SUMF1" in names
    hit = next(x for x in f if x.feature == "FGE/SUMF1")
    assert hit.focus_n == 9 and hit.bg_n == 1
    assert hit.odds_ratio > 5 and hit.p_value < 0.01
    assert hit.spread == 9


def test_within_dataset_mode_is_labelled_exploratory():
    feats, focus, _ = _features()
    f, summary = co.analyse(feats, focus)
    assert summary["evidence_mode"] == "within_dataset"
    assert "exploratory" in summary["note"]
    assert "exploratory" in f[0].reason


def test_held_out_background_changes_the_label():
    feats, focus, _ = _features()
    bg = {t for t in feats if t.startswith("bg")}
    f, summary = co.analyse(feats, focus, background=bg, held_out=True)
    assert summary["evidence_mode"] == "held_out_background"
    assert "exploratory" not in f[0].reason


def test_single_genus_carriers_are_flagged_as_possible_ancestry():
    feats, focus, spread = _features(genus_spread=False)
    f, _ = co.analyse(feats, focus, spread_by=spread)
    hit = next(x for x in f if x.feature == "FGE/SUMF1")
    assert hit.spread == 1
    assert "track ancestry" in hit.reason


def test_rare_features_below_min_count_are_dropped():
    feats, focus, _ = _features(hits_with_fge=2)
    f, _ = co.analyse(feats, focus)
    assert "FGE/SUMF1" not in [x.feature for x in f]


def test_unenriched_feature_is_dropped():
    feats, focus, _ = _features(hits_with_fge=10, bg_with_fge=40)
    f, _ = co.analyse(feats, focus)
    assert "FGE/SUMF1" not in [x.feature for x in f]


def _null_split_counts(trials=25, max_q=0.02, seed=5):
    """How many features survive on random splits of a set with no real signal."""
    import random
    rng = random.Random(seed)
    feats = {}
    genera = [f"genus{i}" for i in range(12)]
    for i in range(256):
        t = f"n{i:03d}"
        feats[t] = [("domain", "NTN"), ("genus", rng.choice(genera))]
        if rng.random() < 0.06:
            feats[t].append(("domain", rng.choice(["TPR", "HEAT", "WD40"])))
    names = sorted(feats)
    counts = []
    for _ in range(trials):
        rng.shuffle(names)
        f, _ = co.analyse(feats, set(names[:128]), background=set(names[128:]),
                          held_out=True, max_q=max_q)
        counts.append(len(f))
    return counts


def test_null_split_false_positive_rate_is_low():
    """A random split should usually return nothing.

    It will not always. With only a dozen or so features in the table, BH has
    little to correct over, and a borderline q near the cutoff is noise rather
    than a finding. The gate is calibration, not a guarantee -- so this asserts
    a rate, which is the honest claim.
    """
    counts = _null_split_counts()
    assert sum(1 for c in counts if c == 0) >= 0.8 * len(counts)
    assert sum(counts) / len(counts) < 0.5
    # the looser cutoff is measurably worse, which is why it is not the default
    loose = _null_split_counts(max_q=0.10)
    assert sum(loose) > sum(counts)


def test_the_q_gate_is_what_removes_them():
    gated = sum(_null_split_counts(trials=10))
    ungated = sum(_null_split_counts(trials=10, max_q=1.0))
    assert ungated > gated


def test_near_miss_features_are_reported_as_suppressed_not_dropped():
    """A small set costs real signal at the gate; that must be visible."""
    feats, focus, _ = _features(hits_with_fge=8, hits=24, bg_with_fge=1, bg=17)
    f, summary = co.analyse(feats, focus)
    assert "FGE/SUMF1" not in [x.feature for x in f]
    assert summary["features_suppressed_near_cutoff"] == 1
    assert "FGE/SUMF1" in summary["suppressed_detail"]
    assert "power limit" in summary["suppressed_note"]
    loose, _ = co.analyse(feats, focus, max_q=0.20)
    assert "FGE/SUMF1" in [x.feature for x in loose]


def test_q_value_is_reported():
    feats, focus, _ = _features()
    f, _ = co.analyse(feats, focus)
    assert f[0].q_value <= 0.02


def test_empty_focus_is_an_error():
    feats, _, _ = _features()
    with pytest.raises(CooccurError):
        co.analyse(feats, set())


def test_cooccur_emits_features():
    feats, focus, _ = _features()
    f, _ = co.analyse(feats, focus)
    lines = co.emit_features(f)
    assert lines and "--surfaced-by cooccurrence" in lines[0]


def test_feature_table_round_trip(tmp_path):
    p = tmp_path / "feat.tsv"
    p.write_text("target\tfeature_type\tfeature\nT1\tdomain\tFGE\nT1\tdomain\tTIR\nT2\tdomain\tTIR\n")
    feats = co.load_features(str(p))
    assert feats["T1"] == [("domain", "FGE"), ("domain", "TIR")]


def test_feature_table_needs_the_right_columns(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("a\tb\n1\t2\n")
    with pytest.raises(CooccurError):
        co.load_features(str(p))
