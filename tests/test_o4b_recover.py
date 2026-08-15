"""O4b recovery contract tests.

The recovery republishes a verdict the campaign computed but never wrote.
The danger is that it becomes a back door for editing a result: so these
tests pin that it recomputes the gates from the preserved statistics with
the frozen functions, that it refuses to publish anything that does not
reproduce the preserved verdict bit-for-bit, that it verifies the frozen
decision surface, and that it never draws, resamples, or changes a gate.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import o4b_recover as rec  # noqa: E402

_PREREG = _REPO / "docs" / "prereg"
_FREEZE_SHA = "715865abc684224785e71a3130e17c50db35f947"
_RESERVATION_OBJ = "46acee340bc247511546964b2925953721d5bb59"


def _incident() -> dict:
    return json.loads((_PREREG / "p14_o4b_incident.json")
                      .read_text(encoding="utf-8"))


# ------------------------------------------- the recomputation is faithful

def test_recover_reproduces_the_preserved_verdicts_bit_for_bit():
    result = rec.recover()
    pres = _incident()["preserved"]["completed_gates"]

    for k in ("status", "n", "v_s1_lo", "v_s1_hi",
              "identified_discrepancy", "band_abs"):
        assert result["g1"][k] == pres["g1"][k], f"g1.{k}"
    for k in ("status", "n", "leaking_points", "leak_upper_abs",
              "budget_abs"):
        assert result["g2"][k] == pres["g2"][k], f"g2.{k}"
    assert result["g1"]["status"] == "concordant"
    assert result["g2"]["status"] == "concordant"


def test_g1_interval_is_recomputed_from_only_n_mean_var():
    """The independence claim: G1's number comes from the preserved
    sufficient statistics through the frozen empirical-Bernstein formula,
    nothing else."""

    ck = json.loads((_PREREG / "p14_o4b_checkpoint.json")
                    .read_text(encoding="utf-8"))
    stats = {"n": ck["statistics"]["n"],
             "mean_z": ck["statistics"]["mean_z"],
             "var_z": ck["statistics"]["var_z"]}
    g1 = rec.recompute_g1(stats)
    assert g1 == rec.recover()["g1"]


def test_g2_leak_lower_is_carried_not_recomputed_and_labelled():
    result = rec.recover()
    pres = _incident()["preserved"]["completed_gates"]["g2"]
    assert result["g2"]["leak_lower_abs"] == pres["leak_lower_abs"]
    assert "NOT independently" in result["g2"]["leak_lower_abs_provenance"]
    assert "5c" in result["g2"]["leak_lower_abs_provenance"]


# ------------------------------------------------------- the integrity gate

def test_the_integrity_gate_refuses_a_tampered_statistic(monkeypatch):
    """If the recomputation and the preserved verdict disagree, nothing
    is published -- the recovery is not a channel for editing a result."""

    real = rec.recompute_g1

    def tampered(stats):
        out = real(stats)
        out["v_s1_hi"] = out["v_s1_hi"] + 1.0     # a value that cannot be
        return out                                # reproduced from record

    monkeypatch.setattr(rec, "recompute_g1", tampered)
    with pytest.raises(rec.RecoveryError, match="disagree"):
        rec.recover()


def test_a_mismatched_status_is_caught():
    with pytest.raises(rec.RecoveryError, match="g2.status"):
        rec._assert_reproduces(
            "g2", {"status": "discordant"}, {"status": "concordant"},
            ("status",))


# --------------------------------------------------- the frozen surface

def test_verify_frozen_surface_passes_on_the_frozen_tree():
    manifest = rec.verify_frozen_surface()
    assert "files" in manifest


def test_a_drifted_decision_file_refuses_recovery(monkeypatch):
    """Every executed-manifest file but the ledger must match byte for
    byte; a drifted empirical_bernstein would silently change the
    verdict, so it aborts instead."""

    real = rec.run._sha256
    target = _REPO / "experiments" / "oracle" / "empirical_bernstein.py"

    def drifted(path):
        if Path(path) == target:
            return "0" * 64
        return real(path)

    monkeypatch.setattr(rec.run, "_sha256", drifted)
    with pytest.raises(rec.RecoveryError, match="drifted"):
        rec.verify_frozen_surface()


def test_recovery_refuses_if_an_o4b_seed_is_still_fresh(monkeypatch):
    monkeypatch.setattr(rec.ledger, "FRESH_PROBE_SCALARS",
                        {"o4b_g1_audit": 40_000_401})
    with pytest.raises(rec.RecoveryError, match="FRESH"):
        rec.verify_frozen_surface()


def test_a_later_campaigns_fresh_seed_does_not_block_recovery(monkeypatch):
    """FRESH_PROBE_SCALARS is the program-wide active ledger, so a future
    campaign freezing its own scalar is a normal state. Only the O4b
    seeds' own retirement is the recovery's business (review PR #77 R2)."""

    monkeypatch.setattr(rec.ledger, "FRESH_PROBE_SCALARS",
                        {"o9_future_audit": 40_000_999})
    # the O4b seeds are still OBSERVED, so the surface check passes
    assert rec.verify_frozen_surface()["files"]


# ------------------------------------------------ provenance of the result

def test_it_is_marked_a_recovery_and_draws_nothing():
    r = rec.recover()
    assert r["kind"] == "results"
    assert r["run_kind"] == "recovered_completed_campaign"
    rc = r["recovery"]
    assert rc["reproduced_preserved_verdicts_bit_exact"] is True
    assert rc["no_new_seed"] is True
    assert rc["no_solver_call"] is True
    assert rc["no_resampling"] is True
    assert rc["no_gate_change"] is True
    assert "reproduction" in rc["seeds_status"]


def test_the_reservation_and_freeze_provenance():
    r = rec.recover()
    assert r["freeze_sha"] == _FREEZE_SHA
    assert r["reservation"]["object"] == _RESERVATION_OBJ
    assert r["reservation"]["ref"] == "refs/o4b/reservation"
    assert r["reservation"]["seeds_spent"] is True
    assert r["reservation"]["verified_at_recovery"] is True
    assert r["seeds"] == {"o4b_g1_audit": 40_000_401,
                          "o4b_g2_leakage": 40_000_411}


def test_recover_refuses_a_record_from_a_different_freeze(monkeypatch):
    real = rec._load

    def wrong(path):
        d = real(path)
        if path == rec._INCIDENT:
            d = dict(d, freeze_sha="0" * 40)
        return d

    monkeypatch.setattr(rec, "_load", wrong)
    with pytest.raises(rec.RecoveryError, match="freeze_sha"):
        rec.recover()


# ----------------------------------- authenticating the preserved inputs

def test_a_tampered_preserved_blob_is_caught(monkeypatch):
    """The blob pin: any edit to the incident or checkpoint changes a
    sha256 checked here, before a single value is read from it."""

    real = rec.run._sha256

    def altered(path):
        if Path(path) == rec._INCIDENT:
            return "f" * 64
        return real(path)

    monkeypatch.setattr(rec.run, "_sha256", altered)
    with pytest.raises(rec.RecoveryError, match="sha256"):
        rec._authenticate()


def test_a_single_file_g2_edit_diverges_from_the_second_copy(monkeypatch):
    """The cross-check: the incident and the checkpoint each carry the G2
    block, so editing one alone -- even self-consistently -- disagrees
    with the other and is caught before recomputation (review PR #77
    R1)."""

    real = rec._load

    def tamper_incident_g2(path):
        d = real(path)
        if path == rec._INCIDENT:
            d = json.loads(json.dumps(d))         # deep copy
            g2 = d["preserved"]["completed_gates"]["g2"]
            g2["n"] = g2["n"] + 1                 # self-consistent-looking
        return d

    monkeypatch.setattr(rec, "_load", tamper_incident_g2)
    with pytest.raises(rec.RecoveryError, match="G2 block"):
        rec._authenticate()


def test_the_pinned_blob_hashes_match_the_committed_records():
    """The pins are the committed record, so recovery reads exactly what
    PR #76 preserved."""

    assert rec.run._sha256(rec._INCIDENT) == rec.INCIDENT_SHA256
    assert rec.run._sha256(rec._CHECKPOINT) == rec.CHECKPOINT_SHA256


# ------------------------------------- the published artifact, once it exists

def test_the_published_result_matches_the_recomputation():
    """If the result has been published, its scientific content is
    exactly what recover() produces (environment is call-time, so it is
    not compared)."""

    path = _PREREG / "p14_o4b_results.json"
    if not path.exists():
        pytest.skip("result not yet published")
    published = json.loads(path.read_text(encoding="utf-8"))
    r = rec.recover()
    for key in ("kind", "run_kind", "freeze_sha", "manifest_digest",
                "seeds", "g1", "g2", "reservation", "recovery", "order"):
        assert published[key] == r[key], key
