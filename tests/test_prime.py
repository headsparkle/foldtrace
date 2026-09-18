"""Tests for the prior-art register.

The two invariants worth testing are negative ones: nothing can write a
verified claim except the verification step, and the register cannot change a
single state call.
"""
from __future__ import annotations

import json
import os

import pytest

from foldtrace import prime
from foldtrace.prime import PrimeError


def _reg(tmp_path):
    reg = prime.new_register("hepn", "5YEP:B", scope="HEPN ribonuclease fold",
                             tools=["pubmed"])
    reg.novelty_criteria = ["A HEPN-fold family not in clan CL0291 at family level."]
    reg.non_novelty_criteria = [
        "Recovery of DUF86/UPF0331: PF01934 is a HepT-like HEPN family inside CL0291."]
    reg.claims.append(prime.Claim(
        claim_id="PA-001", claim_type="residue_evidence",
        claim="Arg97 positions the scissile phosphate; His102 is the general base.",
        identifier="10.1074/jbc.RA118.002421", residues=["R97", "H102"],
        evidence_strength="experimental"))
    reg.claims.append(prime.Claim(
        claim_id="PA-002", claim_type="annotation_caveat",
        claim="PF01934 (DUF86) is a HEPN-clan family whose cross-reference is often absent.",
        identifier="PF01934", evidence_strength="annotation"))
    return reg


def _verify_all(reg):
    prime.apply_verification(reg, [
        {"claim_id": c.claim_id, "status": "verified", "method": "doi resolved"}
        for c in reg.pre_run_claims()])


def _sites(tmp_path, nums=(97, 102)):
    p = tmp_path / "sites.tsv"
    rows = ["label\tref_resnum\texpected\trole\tkey\taltered"]
    for n in nums:
        rows.append(f"X{n}\t{n}\tH\tcatalytic\t1\t")
    p.write_text("\n".join(rows) + "\n")
    return str(p)


# ------------------------------------------------------------------ invariants
def test_new_claims_are_always_unverified():
    c = prime.Claim(claim_id="PA-001", claim_type="mechanism", claim="x")
    assert c.source_verified == "unverified"
    assert c.identifier_resolved is False


def test_verification_requires_a_method(tmp_path):
    reg = _reg(tmp_path)
    with pytest.raises(PrimeError, match="without a method"):
        prime.apply_verification(reg, [{"claim_id": "PA-001", "status": "verified"}])


def test_freeze_rejects_unverified_claims(tmp_path):
    reg = _reg(tmp_path)
    with pytest.raises(PrimeError, match="must be verified"):
        prime.freeze(reg, _sites(tmp_path), out_path=str(tmp_path / "FREEZE.md"))


def test_a_documented_verification_failure_does_not_block_freeze(tmp_path):
    """Trying to verify and failing is provenance, not an error -- unless a site rests on it."""
    reg = _reg(tmp_path)
    prime.apply_verification(reg, [
        {"claim_id": "PA-001", "status": "verified", "method": "doi resolved"},
        {"claim_id": "PA-002", "status": "failed", "method": "no resolvable identifier found"}])
    assert prime.validate(reg, sites_path=_sites(tmp_path)) == []


def test_a_failed_claim_backing_a_declared_site_blocks_freeze(tmp_path):
    reg = _reg(tmp_path)
    prime.apply_verification(reg, [
        {"claim_id": "PA-001", "status": "failed", "method": "numbering could not be confirmed"},
        {"claim_id": "PA-002", "status": "verified", "method": "InterPro retrieved"}])
    problems = prime.validate(reg, sites_path=_sites(tmp_path))
    assert any("site R97 is declared" in p or "site H102 is declared" in p for p in problems)


def test_freeze_rejects_a_register_with_nothing_verified(tmp_path):
    reg = _reg(tmp_path)
    prime.apply_verification(reg, [
        {"claim_id": c.claim_id, "status": "failed", "method": "not found"}
        for c in reg.pre_run_claims()])
    problems = prime.validate(reg, sites_path=_sites(tmp_path))
    assert any("nothing confirmed" in p for p in problems)


def test_freeze_rejects_empty_non_novelty_criteria(tmp_path):
    reg = _reg(tmp_path)
    _verify_all(reg)
    reg.non_novelty_criteria = []
    problems = prime.validate(reg, sites_path=_sites(tmp_path))
    assert any("non_novelty_criteria" in p for p in problems)


def test_freeze_rejects_site_without_backing_evidence(tmp_path):
    reg = _reg(tmp_path)
    _verify_all(reg)
    sites = _sites(tmp_path, nums=(97, 102, 555))
    problems = prime.validate(reg, sites_path=sites)
    assert any("555" in p for p in problems)


def test_freeze_writes_hashes_and_stamps(tmp_path):
    reg = _reg(tmp_path)
    _verify_all(reg)
    out = str(tmp_path / "FREEZE.md")
    text = prime.freeze(reg, _sites(tmp_path), out_path=out)
    assert reg.frozen_at and reg.sha256
    assert reg.sha256 in text
    assert "Would not count as new" in text
    assert os.path.exists(out)


def test_post_run_claims_do_not_change_the_pre_run_hash(tmp_path):
    reg = _reg(tmp_path)
    _verify_all(reg)
    before = prime.sha256_text(prime.canonical_body(reg))
    reg.claims.append(prime.Claim(
        claim_id="PA-003", claim_type="described_architecture",
        claim="Fused MNT-HEPN systems are regulated by oligoAMPylation.",
        discovered_post_run=True, trigger="EF-002 fused architecture"))
    assert prime.sha256_text(prime.canonical_body(reg)) == before


# --------------------------------------------------------------- sites skeleton
def test_sites_skeleton_only_uses_verified_evidence(tmp_path):
    reg = _reg(tmp_path)
    assert "no verified residue_evidence" in prime.sites_skeleton(reg)
    _verify_all(reg)
    text = prime.sites_skeleton(reg)
    assert "R97\t97\tR" in text and "H102\t102\tH" in text
    # the Asp/Glu warning must survive, it is the point of the comment block
    assert "Asp or Glu" in text


# ---------------------------------------------------------------- adjudication
def test_cannot_adjudicate_before_freeze(tmp_path):
    reg = _reg(tmp_path)
    _verify_all(reg)
    with pytest.raises(PrimeError, match="unfrozen"):
        prime.add_adjudication(reg, "feature", "cooccurrence",
                               "no_prior_description_found", "claim")


def test_partially_described_requires_both_components(tmp_path):
    reg = _reg(tmp_path)
    _verify_all(reg)
    prime.freeze(reg, _sites(tmp_path), out_path=str(tmp_path / "FREEZE.md"))
    with pytest.raises(PrimeError, match="undescribed_component"):
        prime.add_adjudication(
            reg, feature="TIR-FGE architecture", surfaced_by="cooccurrence",
            prior_art_status="partially_described",
            permitted_claim="x", described="44 members published as bNLRs")


def test_supporting_claim_must_exist(tmp_path):
    reg = _reg(tmp_path)
    _verify_all(reg)
    prime.freeze(reg, _sites(tmp_path), out_path=str(tmp_path / "FREEZE.md"))
    with pytest.raises(PrimeError, match="PA-999"):
        prime.add_adjudication(reg, "f", "manual", "already_described", "c",
                               supporting=["PA-999"])


# ---------------------------------------------------------------- novelty gate
def test_gate_blocks_novelty_language_with_no_adjudication(tmp_path):
    reg = _reg(tmp_path)
    _verify_all(reg)
    prime.freeze(reg, _sites(tmp_path), out_path=str(tmp_path / "FREEZE.md"))
    v = prime.novelty_gate(reg, "We report a new family of HEPN proteins.")
    assert v and "no adjudications" in v[0]


def test_gate_blocks_novelty_attached_to_an_already_described_feature(tmp_path):
    reg = _reg(tmp_path)
    _verify_all(reg)
    prime.freeze(reg, _sites(tmp_path), out_path=str(tmp_path / "FREEZE.md"))
    prime.add_adjudication(
        reg, feature="750 of 836 hits are DUF86/UPF0331", surfaced_by="discordance",
        prior_art_status="already_described",
        permitted_claim="Recovery of a known family with missing per-record cross-reference.",
        supporting=["PA-002"])
    v = prime.novelty_gate(
        reg, "The 750 DUF86/UPF0331 records are a previously unrecognised HEPN family.")
    assert v and "already_described" in v[0]


def test_gate_requires_the_undescribed_component_to_be_named(tmp_path):
    reg = _reg(tmp_path)
    _verify_all(reg)
    prime.freeze(reg, _sites(tmp_path), out_path=str(tmp_path / "FREEZE.md"))
    prime.add_adjudication(
        reg, feature="TIR-FGE", surfaced_by="cooccurrence",
        prior_art_status="partially_described",
        described="Receptor-embedded members are published as bacterial NACHT proteins.",
        undescribed="The compact two-domain fusions lacking a nucleotide-binding sensor.",
        permitted_claim="The compact fusions are the undescribed component.")
    bad = prime.novelty_gate(reg, "We describe a novel TIR-FGE defence architecture.")
    assert bad
    good = prime.novelty_gate(
        reg, "The compact fusions, lacking a nucleotide-binding sensor, are novel; "
             "the receptor-embedded members are not.")
    assert not good


# ---------------------------------------------------------------- round-trip
def test_round_trip_preserves_everything(tmp_path):
    reg = _reg(tmp_path)
    _verify_all(reg)
    prime.freeze(reg, _sites(tmp_path), out_path=str(tmp_path / "FREEZE.md"))
    prime.add_adjudication(reg, "f", "rare_residue", "no_prior_description_found", "c")
    p = str(tmp_path / "register.json")
    prime.save(reg, p)
    back = prime.load(p)
    assert back.sha256 == reg.sha256
    assert prime.canonical_body(back) == prime.canonical_body(reg)
    assert len(back.adjudications) == 1


def test_json_nulls_for_optional_lists_are_tolerated(tmp_path):
    """Hand-written registers carry `null` where a list is absent."""
    p = tmp_path / "reg.json"
    p.write_text(json.dumps({
        "fold_id": "hepn", "reference_structure": "5YEP:B",
        "scope_statement": "s", "novelty_criteria": ["a"], "non_novelty_criteria": ["b"],
        "retrieval_tools": None,
        "claims": [{"claim_id": "PA-001", "claim_type": "mechanism", "claim": "x",
                    "residues": None, "source_verified": "verified",
                    "verification_method": "m", "identifier_resolved": True}],
        "adjudications": []}))
    reg = prime.load(str(p))
    assert reg.claims[0].residues == []
    assert reg.retrieval_tools == []
    assert prime.validate(reg) == []  # must not raise


def test_load_rejects_unknown_fields(tmp_path):
    p = tmp_path / "bad.json"
    p.write_text(json.dumps({"fold_id": "x", "reference_structure": "y",
                             "expected_outcome": "a new family"}))
    with pytest.raises(PrimeError, match="unknown register fields"):
        prime.load(str(p))


def test_gate_flags_novelty_that_names_no_adjudicated_feature(tmp_path):
    reg = _reg(tmp_path)
    _verify_all(reg)
    prime.freeze(reg, _sites(tmp_path), out_path=str(tmp_path / "FREEZE.md"))
    prime.add_adjudication(
        reg, feature="750 of 836 hits are DUF86/UPF0331", surfaced_by="discordance",
        prior_art_status="already_described",
        permitted_claim="Recovery of a known family with a missing cross-reference.")
    v = prime.novelty_gate(reg, "We report an unprecedented phosphatase cassette.")
    assert v and "names no adjudicated feature" in v[0]


def test_gate_is_scoped_to_the_sentence(tmp_path):
    """A novel claim about one feature must not be flagged against another."""
    reg = _reg(tmp_path)
    _verify_all(reg)
    prime.freeze(reg, _sites(tmp_path), out_path=str(tmp_path / "FREEZE.md"))
    prime.add_adjudication(
        reg, feature="750 of 836 hits are DUF86/UPF0331", surfaced_by="discordance",
        prior_art_status="already_described",
        permitted_claim="Recovery of a known family with a missing cross-reference.")
    prime.add_adjudication(
        reg, feature="fused HEPN-MNT records annotated as nucleotidyltransferase",
        surfaced_by="cooccurrence", prior_art_status="partially_described",
        described="Fused MNT-HEPN systems and their regulation are established.",
        undescribed="Four recovered fusions have lost the catalytic histidine.",
        permitted_claim="The histidine-lost fusions are the undescribed part.")
    text = ("The 750 DUF86/UPF0331 records belong to PF01934 at family level. "
            "Among the fused nucleotidyltransferase-annotated records, four have "
            "lost the catalytic histidine, which is the novel component.")
    assert not prime.novelty_gate(reg, text)
