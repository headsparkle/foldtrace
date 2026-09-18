"""Tests for the discordance detector.

The fixtures are the real cases: a TIR fold whose record says CD38, a HEPN fold
whose record says DUF86, sequence-mapped against structure-mapped glutamate
calls, and a PETase-fold entry carrying an engineered Ser-to-Ala.
"""
from __future__ import annotations

import pytest

from foldtrace import discord as dis
from foldtrace.discord import DiscordError


def _calls(tmp_path, rows, name="calls.tsv", labels=("Glu642",)):
    cols = ["candidate_id", "tmalign_tm_norm_ref", "fold_ok"]
    for lab in labels:
        cols += [f"{lab}_obs", f"{lab}_offsetA", f"{lab}_state"]
    cols += ["verdict", "state_reason"]
    p = tmp_path / name
    lines = ["\t".join(cols)]
    for r in rows:
        lines.append("\t".join(str(r.get(c, "NA")) for c in cols))
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def _annos(tmp_path, rows, name="annos.tsv"):
    cols = ["target", "name", "pfam", "clan"]
    p = tmp_path / name
    lines = ["\t".join(cols)]
    for r in rows:
        lines.append("\t".join(str(r.get(c, "")) for c in cols))
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def _vocab(tmp_path, rows, name="vocab.tsv"):
    p = tmp_path / name
    lines = ["class\tterm\tnote"]
    for cls, term, note in rows:
        lines.append(f"{cls}\t{term}\t{note}")
    p.write_text("\n".join(lines) + "\n")
    return str(p)


TIR_VOCAB = [
    ("family_name", "tir", ""),
    ("family_accession", "CL0173", ""),
    ("alias", "duf4062", "PF13271, a TIR-clan family"),
    ("alias", "sefir", "PF08357, a TIR-clan family"),
    ("foreign", "adp-ribosyl cyclase", "CD38/CD157-like ectoenzyme"),
]


# ------------------------------------------------------------------ loading
def test_site_labels_and_target_column(tmp_path):
    p = _calls(tmp_path, [{"candidate_id": "P1", "verdict": "retained"}])
    c = dis.load_calls(p)
    assert c.labels == ["Glu642"] and c.targets == {"P1"}


def test_missing_target_column_is_an_error(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("foo\tbar\n1\t2\n")
    with pytest.raises(DiscordError):
        dis.load_calls(str(p))


def test_unknown_vocab_class_is_an_error(tmp_path):
    with pytest.raises(DiscordError, match="unknown vocab class"):
        dis.load_vocab(_vocab(tmp_path, [("nonsense", "x", "")]))


# --------------------------------------------------------------- annotation
def _scan(tmp_path, calls_rows, anno_rows, vocab_rows=TIR_VOCAB):
    return dis.scan_annotation(
        dis.load_calls(_calls(tmp_path, calls_rows)),
        dis.load_annotations(_annos(tmp_path, anno_rows)),
        dis.load_vocab(_vocab(tmp_path, vocab_rows)))


def test_name_contradicting_the_fold_outranks_everything(tmp_path):
    """The CD38 case: intact SARM1-like site, record says ADP-ribosyl cyclase."""
    f = _scan(tmp_path,
              [{"candidate_id": "Q9FHM1", "tmalign_tm_norm_ref": 0.82, "fold_ok": "true",
                "Glu642_state": "retained", "verdict": "retained"}],
              [{"target": "Q9FHM1",
                "name": "ADP-ribosyl cyclase / cyclic ADP-ribose hydrolase"}])
    assert len(f) == 1
    assert f[0].kind == "name_contradicts_fold"
    assert "CD38" in f[0].observation
    assert f[0].score > 0.7


def test_known_alias_is_demoted_not_promoted(tmp_path):
    """DUF4062/SEFIR, and by the same logic DUF86: an annotation gap, not a find."""
    f = _scan(tmp_path,
              [{"candidate_id": "A1", "tmalign_tm_norm_ref": 0.90, "fold_ok": "true",
                "Glu642_state": "retained", "verdict": "retained"}],
              [{"target": "A1", "name": "DUF4062 domain-containing protein"}])
    assert f[0].kind == "annotation_gap"
    assert f[0].score < 0.3
    assert "not family discovery" in f[0].reason


def test_recognised_family_records_produce_nothing(tmp_path):
    f = _scan(tmp_path,
              [{"candidate_id": "A1", "tmalign_tm_norm_ref": 0.9, "fold_ok": "true",
                "Glu642_state": "retained", "verdict": "retained"}],
              [{"target": "A1", "name": "TIR domain protein", "clan": "CL0173"}])
    assert f == []


def test_accession_beats_a_foreign_name(tmp_path):
    """A record carrying the family accession is not contradicted by its name."""
    f = _scan(tmp_path,
              [{"candidate_id": "A1", "tmalign_tm_norm_ref": 0.9, "fold_ok": "true",
                "Glu642_state": "retained", "verdict": "retained"}],
              [{"target": "A1", "name": "ADP-ribosyl cyclase", "clan": "CL0173"}])
    assert f == []


def test_unannotated_intact_site_is_reported(tmp_path):
    f = _scan(tmp_path,
              [{"candidate_id": "A1", "tmalign_tm_norm_ref": 0.88, "fold_ok": "true",
                "Glu642_state": "retained", "verdict": "retained"}],
              [{"target": "A1", "name": "Uncharacterized protein"}])
    assert f[0].kind == "unannotated_intact_site"


def test_lost_and_failed_folds_are_not_annotation_findings(tmp_path):
    f = _scan(tmp_path,
              [{"candidate_id": "A1", "tmalign_tm_norm_ref": 0.9, "fold_ok": "true",
                "Glu642_state": "lost", "verdict": "lost"},
               {"candidate_id": "A2", "tmalign_tm_norm_ref": 0.3, "fold_ok": "false",
                "Glu642_state": "unresolved", "verdict": "unresolved"}],
              [{"target": "A1", "name": "Uncharacterized protein"},
               {"target": "A2", "name": "Uncharacterized protein"}])
    assert f == []


def test_marginal_folds_score_below_confident_ones(tmp_path):
    rows = [{"candidate_id": "HI", "tmalign_tm_norm_ref": 0.95, "fold_ok": "true",
             "Glu642_state": "retained", "verdict": "retained"},
            {"candidate_id": "LO", "tmalign_tm_norm_ref": 0.55, "fold_ok": "true",
             "Glu642_state": "retained", "verdict": "retained"}]
    annos = [{"target": t, "name": "ADP-ribosyl cyclase"} for t in ("HI", "LO")]
    f = dis.rank(_scan(tmp_path, rows, annos))
    assert f[0].target == "HI" and f[1].target == "LO"


# -------------------------------------------------------------------- paired
def test_state_flip_between_channels(tmp_path):
    """21 of 23 sequence-derived glutamate losses were overturned by superposition."""
    a = dis.load_calls(_calls(tmp_path, [
        {"candidate_id": "P1", "tmalign_tm_norm_ref": 0.80, "fold_ok": "true",
         "Glu642_obs": "K640", "Glu642_state": "lost", "verdict": "lost"}], name="seq.tsv"))
    b = dis.load_calls(_calls(tmp_path, [
        {"candidate_id": "P1", "tmalign_tm_norm_ref": 0.80, "fold_ok": "true",
         "Glu642_obs": "E644", "Glu642_state": "retained", "verdict": "retained"}], name="str.tsv"))
    f = dis.rank(dis.scan_paired(a, b, "sequence", "structure"))
    kinds = [x.kind for x in f]
    # the verdict flip already reports this site; it must not be counted twice
    assert kinds == ["state_flip"]
    top = f[0]
    assert top.score > 0.7
    assert "at most one can be right" in top.reason


def test_resolution_disagreement_scores_below_a_hard_flip(tmp_path):
    a = dis.load_calls(_calls(tmp_path, [
        {"candidate_id": "P1", "tmalign_tm_norm_ref": 0.9, "fold_ok": "true",
         "Glu642_state": "retained", "verdict": "retained"}], name="a.tsv"))
    b = dis.load_calls(_calls(tmp_path, [
        {"candidate_id": "P1", "tmalign_tm_norm_ref": 0.9, "fold_ok": "true",
         "Glu642_state": "unresolved", "verdict": "unresolved",
         "state_reason": "CA offset 5.10 A > 4.00 A"}], name="b.tsv"))
    f = [x for x in dis.scan_paired(a, b) if x.kind == "resolution_disagreement"]
    assert f and f[0].score < 0.7
    assert "abstention reason" in f[0].reason


def test_monomer_to_dimer_flip_is_caught(tmp_path):
    """A0A852UQU6 holds its shift in the dimer; A0A7Y5GXD8 does not."""
    mono = dis.load_calls(_calls(tmp_path, [
        {"candidate_id": "A0A852UQU6", "tmalign_tm_norm_ref": 0.9, "fold_ok": "true",
         "H102_state": "lost", "verdict": "lost"},
        {"candidate_id": "A0A7Y5GXD8", "tmalign_tm_norm_ref": 0.9, "fold_ok": "true",
         "H102_state": "lost", "verdict": "lost"}], name="mono.tsv", labels=("H102",)))
    dimer = dis.load_calls(_calls(tmp_path, [
        {"candidate_id": "A0A852UQU6", "tmalign_tm_norm_ref": 0.9, "fold_ok": "true",
         "H102_state": "lost", "verdict": "lost"},
        {"candidate_id": "A0A7Y5GXD8", "tmalign_tm_norm_ref": 0.9, "fold_ok": "true",
         "H102_state": "retained", "verdict": "retained"}], name="dimer.tsv", labels=("H102",)))
    f = dis.scan_paired(mono, dimer, "monomer", "dimer")
    assert {x.target for x in f} == {"A0A7Y5GXD8"}


def test_targets_in_only_one_channel_are_ignored(tmp_path):
    a = dis.load_calls(_calls(tmp_path, [
        {"candidate_id": "P1", "verdict": "retained", "fold_ok": "true"}], name="a.tsv"))
    b = dis.load_calls(_calls(tmp_path, [
        {"candidate_id": "P2", "verdict": "lost", "fold_ok": "true"}], name="b.tsv"))
    assert dis.scan_paired(a, b) == []


# ------------------------------------------------------------------ knockout
KO_LABELS = ("Ser160", "Asp206", "His237")


def test_point_knockout_signature(tmp_path):
    """A0AA82WPD4: a deposited S-to-A construct resident as a natural entry."""
    c = dis.load_calls(_calls(tmp_path, [
        {"candidate_id": "A0AA82WPD4", "tmalign_tm_norm_ref": 0.97, "fold_ok": "true",
         "Ser160_obs": "A160", "Ser160_state": "lost",
         "Asp206_obs": "D206", "Asp206_state": "retained",
         "His237_obs": "H237", "His237_state": "retained",
         "verdict": "lost"}], labels=KO_LABELS))
    f = dis.scan_knockout(c)
    assert len(f) == 1 and f[0].kind == "point_knockout_signature"
    assert "engineered construct" in f[0].reason


def test_a_common_single_site_loss_is_the_folds_biology_not_a_construct(tmp_path):
    """Cupins lose the acidic ligand constantly; that is not 19 lab constructs."""
    rows = []
    for i in range(40):
        engineered = i < 19
        rows.append({
            "candidate_id": f"cup{i:03d}", "tmalign_tm_norm_ref": 0.95, "fold_ok": "true",
            "Ser160_obs": ("A160" if engineered else "S160"),
            "Ser160_state": ("lost" if engineered else "retained"),
            "Asp206_obs": "D206", "Asp206_state": "retained",
            "His237_obs": "H237", "His237_state": "retained",
            "verdict": "lost" if engineered else "retained"})
    c = dis.load_calls(_calls(tmp_path, rows, labels=KO_LABELS))
    f = dis.scan_knockout(c)
    assert [x.kind for x in f] == ["recurrent_single_site_loss"]
    assert "19 of 40" in f[0].observation
    assert "too common to be engineered" in f[0].reason


def test_prevalence_needs_a_denominator(tmp_path):
    """On a handful of targets a fraction means nothing, so the rule stands down."""
    rows = [{"candidate_id": f"p{i}", "tmalign_tm_norm_ref": 0.95, "fold_ok": "true",
             "Ser160_obs": "A160", "Ser160_state": "lost",
             "Asp206_obs": "D206", "Asp206_state": "retained",
             "His237_obs": "H237", "His237_state": "retained", "verdict": "lost"}
            for i in range(3)]
    c = dis.load_calls(_calls(tmp_path, rows, labels=KO_LABELS))
    assert {x.kind for x in dis.scan_knockout(c)} == {"point_knockout_signature"}


def test_a_rare_single_site_loss_still_reads_as_a_construct(tmp_path):
    rows = []
    for i in range(40):
        engineered = i < 2
        rows.append({
            "candidate_id": f"pet{i:03d}", "tmalign_tm_norm_ref": 0.95, "fold_ok": "true",
            "Ser160_obs": ("A160" if engineered else "S160"),
            "Ser160_state": ("lost" if engineered else "retained"),
            "Asp206_obs": "D206", "Asp206_state": "retained",
            "His237_obs": "H237", "His237_state": "retained",
            "verdict": "lost" if engineered else "retained"})
    c = dis.load_calls(_calls(tmp_path, rows, labels=KO_LABELS))
    f = dis.scan_knockout(c)
    assert {x.kind for x in f} == {"point_knockout_signature"}
    assert len(f) == 2


def test_natural_double_loss_is_not_a_knockout(tmp_path):
    """A0A2U8FN02 lost nucleophile and base: a real non-catalytic module."""
    c = dis.load_calls(_calls(tmp_path, [
        {"candidate_id": "A0A2U8FN02", "tmalign_tm_norm_ref": 0.92, "fold_ok": "true",
         "Ser160_obs": "A160", "Ser160_state": "lost",
         "Asp206_obs": "N206", "Asp206_state": "lost",
         "His237_obs": "H237", "His237_state": "retained",
         "verdict": "lost"}], labels=KO_LABELS))
    assert dis.scan_knockout(c) == []


def test_loss_to_an_unusual_residue_is_not_a_knockout(tmp_path):
    c = dis.load_calls(_calls(tmp_path, [
        {"candidate_id": "P1", "tmalign_tm_norm_ref": 0.97, "fold_ok": "true",
         "Ser160_obs": "W160", "Ser160_state": "lost",
         "Asp206_obs": "D206", "Asp206_state": "retained",
         "His237_obs": "H237", "His237_state": "retained",
         "verdict": "lost"}], labels=KO_LABELS))
    assert dis.scan_knockout(c) == []


def test_knockout_requires_a_close_superposition(tmp_path):
    c = dis.load_calls(_calls(tmp_path, [
        {"candidate_id": "P1", "tmalign_tm_norm_ref": 0.62, "fold_ok": "true",
         "Ser160_obs": "A160", "Ser160_state": "lost",
         "Asp206_obs": "D206", "Asp206_state": "retained",
         "His237_obs": "H237", "His237_state": "retained",
         "verdict": "lost"}], labels=KO_LABELS))
    assert dis.scan_knockout(c) == []


# -------------------------------------------------------- sequence identity
def _prof(rows):
    """rows: {candidate: sequence string}, positions numbered from 1."""
    from foldtrace.invariant import PositionObs
    out = []
    for cand, seq in rows.items():
        for i, aa in enumerate(seq, start=1):
            out.append(PositionObs(cand, i, aa, i, 0.5, True))
    return out


def test_reporter_construct_across_organisms_is_flagged():
    """The Paper 1 case: EGFP carried into three unrelated assemblies."""
    egfp = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTG"
    prof = _prof({"a1": egfp, "a2": egfp, "a3": egfp,
                  "nat": "MQRTYKLLNDVWAGCHPEFSMIRQTAVNGKLDEHWPYSCTMFAGQRLDVKN"})
    annos = {"a1": {"organism": "Escherichia"}, "a2": {"organism": "Vibrio"},
             "a3": {"organism": "Pseudomonas"}, "nat": {"organism": "Acropora"}}
    f = dis.scan_identity(prof, annos)
    assert len(f) == 1
    assert f[0].kind == "identical_sequence_across_organisms"
    assert f[0].evidence["n_organisms"] == 3
    assert "reporter construct" in f[0].reason


def test_near_identity_short_of_exact_is_still_caught():
    """The published contamination sat at 99.6% identity, not 100%."""
    base = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTG"
    variant = base[:20] + "A" + base[21:]
    prof = _prof({"a1": base, "a2": variant})
    annos = {"a1": {"organism": "Escherichia"}, "a2": {"organism": "Vibrio"}}
    assert dis.scan_identity(prof, annos)


def test_identical_records_from_one_organism_are_not_flagged():
    """Paralogues or duplicate records in one genome are not contamination."""
    seq = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTG"
    prof = _prof({"a1": seq, "a2": seq})
    annos = {"a1": {"organism": "Escherichia"}, "a2": {"organism": "Escherichia"}}
    assert dis.scan_identity(prof, annos) == []


def test_unrelated_sequences_are_not_flagged():
    prof = _prof({"a1": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTG",
                  "a2": "WQRTYKLLNDVWAGCHPEFSMIRQTAVNGKLDEHWPYSCTMFAGQRLDVKN"})
    annos = {"a1": {"organism": "Escherichia"}, "a2": {"organism": "Acropora"}}
    assert dis.scan_identity(prof, annos) == []


# ------------------------------------------------------------------ assembly
def test_findings_are_ranked_and_writable(tmp_path):
    f = [dis.Finding("a", 0.2, "T1", "o", "r"), dis.Finding("b", 0.9, "T2", "o", "r")]
    assert [x.target for x in dis.rank(f)] == ["T2", "T1"]
    out = tmp_path / "f.tsv"
    dis.write(dis.rank(f), str(out))
    lines = out.read_text().strip().split("\n")
    assert lines[0].split("\t") == dis.FINDING_COLUMNS and len(lines) == 3


def test_annotation_gaps_are_never_emitted_as_features():
    f = [dis.Finding("annotation_gap", 0.2, "T1", "o", "r"),
         dis.Finding("name_contradicts_fold", 0.9, "T2", "o", "r")]
    lines = "\n".join(dis.emit_features(f))
    assert "name_contradicts_fold" in lines
    assert "annotation_gap" not in lines
    assert "--status TODO" in lines


def test_a_site_flip_the_verdict_masks_is_still_reported(tmp_path):
    """Two sites moving in opposite directions leave the verdict unchanged."""
    labs = ("S1", "S2")
    a = dis.load_calls(_calls(tmp_path, [
        {"candidate_id": "P1", "tmalign_tm_norm_ref": 0.9, "fold_ok": "true",
         "S1_state": "retained", "S2_state": "lost", "verdict": "lost"}],
        name="a.tsv", labels=labs))
    b = dis.load_calls(_calls(tmp_path, [
        {"candidate_id": "P1", "tmalign_tm_norm_ref": 0.9, "fold_ok": "true",
         "S1_state": "lost", "S2_state": "retained", "verdict": "lost"}],
        name="b.tsv", labels=labs))
    f = dis.scan_paired(a, b)
    assert {x.kind for x in f} == {"site_flip"}
    assert {x.evidence["site"] for x in f} == {"S1", "S2"}


def test_feature_suggestions_are_grammatical():
    one = dis.emit_features([dis.Finding("state_flip", 0.9, "T1", "o", "r")])
    many = dis.emit_features([dis.Finding("state_flip", 0.9, f"T{i}", "o", "r")
                              for i in range(3)])
    assert "1 target showing" in one[0]
    assert "3 targets showing" in many[0]
