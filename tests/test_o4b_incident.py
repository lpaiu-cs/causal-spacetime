"""O4b abort-record contract tests.

Unlike the O4 abort, this one is a publication-wiring defect, not a
scientific fail-closed stop: both gates completed WITH a status and the
run stopped only at the exit provenance check. So the record here
deliberately PRESERVES the gate statistics -- the recovery audit reads
them -- and the tests pin the opposite of O4's negative contract: the
estimates are present, they are flagged non-verdict, and the top-level
verdict stays null, so nothing downstream can promote a preserved
partial into a published result. The lineage must stay verifiable
against the freeze commit and the reservation the run actually claimed.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

_PREREG = _REPO / "docs" / "prereg"
_INCIDENT = _PREREG / "p14_o4b_incident.json"
_CHECKPOINT = _PREREG / "p14_o4b_checkpoint.json"
_INCIDENT_MD = _PREREG / "p14_o4b_incident.md"
_EXECUTED = _PREREG / "p14_o4b_executed_freeze_manifest.json"

_FREEZE_REV = "715865abc684224785e71a3130e17c50db35f947"
_RESERVATION_OBJ = "46acee340bc247511546964b2925953721d5bb59"
_MANIFEST_DIGEST = (
    "cec650b9391af0fc11e4b6bb94455cdbc1a037c18ecec83149d0d4693e7d7be2")


def _incident() -> dict:
    return json.loads(_INCIDENT.read_text(encoding="utf-8"))


def _checkpoint() -> dict:
    return json.loads(_CHECKPOINT.read_text(encoding="utf-8"))


# --------------------------------------------- the negative contract

def test_the_verdict_is_null_and_the_outcome_is_abort():
    """The whole point: no scientific verdict is published, even though
    the gates ran to completion."""

    rec = _incident()
    assert rec["kind"] == "incident"
    assert rec["outcome"] == "ABORT"
    assert rec["verdict"] is None
    assert rec["stage"] == "g2"
    # The CAMPAIGN run published no verdict -- that is what the incident
    # records. A result file may nonetheless exist, published by the
    # separate RECOVERY (a reproduction of the completed run's preserved
    # statistics), and the two coexist by run_kind: the campaign aborted,
    # the recovery republished. What must never appear is a campaign-kind
    # verdict beside this incident.
    result = _PREREG / "p14_o4b_results.json"
    if result.exists():
        published = json.loads(result.read_text(encoding="utf-8"))
        assert published["run_kind"] == "recovered_completed_campaign"
        assert published["kind"] == "results"


def test_stage_g2_is_the_pre_publish_label_not_an_incomplete_gate():
    """Runner-record defect 5b, pinned. The stop happened AFTER run_g2
    wrote g2_complete, at the pre-publish re-read that the runner wraps
    in _staged("g2"). So `stage == "g2"` marks where the run stopped,
    not what completed -- and the same record proves G2 completed. A
    recovery audit keying on `stage` alone would wrongly discard the
    G2 sample; this test and the incident doc forbid that reading."""

    pres = _incident()["preserved"]
    assert _incident()["stage"] == "g2"                 # where it stopped
    assert pres["completed_gates"]["g2"]["status"] == "concordant"
    assert pres["completed_gates"]["g2"]["n"] == 1_072_696   # G2 is done
    assert pres["completed_gates"]["g1"]["status"] == "concordant"


def test_the_preserved_gate_statistics_are_flagged_non_verdict():
    """The distinction the recovery audit rests on: the numbers are
    here, but they are not a status. If a later edit drops the
    non-verdict flags while keeping the numbers, this fails."""

    pres = _incident()["preserved"]
    # the estimates ARE preserved -- this is the recovery material
    assert pres["completed_gates"]["g1"]["status"] == "concordant"
    assert pres["completed_gates"]["g2"]["status"] == "concordant"
    assert pres["g1_partial"]["n"] == 26_200_000
    # ...and every one of them is explicitly disclaimed
    assert "not_a_verdict" in pres["completed_gates_are_not_a_verdict"] \
        or "did not publish" in pres["completed_gates_are_not_a_verdict"]
    assert "preserved partial" in pres["is_not_a_verdict"]
    assert _incident()["verdict"] is None


def test_the_stop_was_wiring_not_science_and_not_a_cap():
    """Neither a fail-closed scientific stop nor a cap: the budget shows
    both caps unbound, and the reason names the internal defect."""

    rec = _incident()
    assert "no claim object to verify" in rec["termination_reason"]
    b = rec["preserved"]["budget"]
    assert b["completed"] < b["max_calls"]           # call cap unbound
    assert b["wall_s"] < b["max_wall_s"]             # wall cap unbound


# ------------------------------------------------------- the lineage

def test_the_reservation_was_claimed_and_the_seeds_are_spent():
    rec = _incident()
    assert rec["freeze_sha"] == _FREEZE_REV
    assert rec["reservation_claimed"] is True
    assert rec["reservation_object"] == _RESERVATION_OBJ
    assert rec["reservation_uncertainty"] is None
    assert rec["seeds_spent"] is True
    assert rec["seeds"] == {"o4b_g1_audit": 40_000_401,
                            "o4b_g2_leakage": 40_000_411}
    probe = subprocess.run(["git", "cat-file", "-e", _FREEZE_REV],
                           cwd=_REPO, capture_output=True)
    assert probe.returncode == 0, "the freeze commit must be reachable"


def test_the_executed_freeze_manifest_is_the_frozen_blob_verbatim():
    """What the campaign actually verified against, preserved even as
    the live manifest moves on with the ledger (the S4/S5/O4 split)."""

    blob = subprocess.run(
        ["git", "show", f"{_FREEZE_REV}:docs/prereg/"
         f"p14_o4b_freeze_manifest.json"],
        cwd=_REPO, capture_output=True, check=True).stdout
    assert _EXECUTED.read_bytes() == blob
    # and its digest is exactly the one the run recorded it ran against
    assert hashlib.sha256(_EXECUTED.read_bytes()).hexdigest() \
        == _MANIFEST_DIGEST
    assert _incident()["manifest_digest"] == _MANIFEST_DIGEST


def test_the_checkpoint_is_a_partial_non_verdict_at_g2_complete():
    cp = _checkpoint()
    assert cp["kind"] == "checkpoint"
    assert cp["stage"] == "g2_complete"
    assert cp["partial"] is True
    assert cp["non_verdict"] is True
    assert cp["freeze_sha"] == _FREEZE_REV


def test_the_g2_complete_checkpoint_rng_is_the_g1_stream():
    """Runner-record defect 5c, pinned. `Campaign.checkpoint` always
    serialises the G1 audit seed and generator; run_g2 uses a separate,
    unserialised generator. So this `g2_complete` record's seed/rng is
    the G1 stream (40,000,401), NOT G2's (40,000,411). A recovery that
    replays G2 from this seed would draw the wrong sample -- the
    incident doc names this and the seeds block keeps both recoverable."""

    cp = _checkpoint()
    assert cp["seed"] == 40_000_401                      # the G1 seed
    assert cp["rng_position"]["bit_generator"] == "PCG64"
    # both seeds are recoverable from the incident, so which is which
    # is never lost even though the checkpoint carries only one
    assert _incident()["seeds"]["o4b_g2_leakage"] == 40_000_411


# --------------------------------------------------------- the ledger

def test_the_drawn_scalars_are_retired_under_their_functional_names():
    """Spent, in OBSERVED, and refused as a fresh allocation -- but NOT
    renamed `aborted`, because a recovery may yet make them a verdict's
    provenance."""

    import probe_seed_ledger as ledger

    assert ledger.OBSERVED_PROBE_SCALARS["o4b_g1_audit"] == 40_000_401
    assert ledger.OBSERVED_PROBE_SCALARS["o4b_g2_leakage"] == 40_000_411
    assert 40_000_401 in ledger.spent_scalars()
    assert 40_000_411 in ledger.spent_scalars()
    # per-seed, not global-empty: FRESH is the program-wide active
    # ledger, and a later campaign's own allocation (e.g. the O5
    # ambiguity pilot) is a normal state (the PR #77 R2 lesson)
    assert "o4b_g1_audit" not in ledger.FRESH_PROBE_SCALARS
    assert "o4b_g2_leakage" not in ledger.FRESH_PROBE_SCALARS
    for name in ("o4b_g1_audit", "o4b_g2_leakage"):
        with pytest.raises(KeyError):
            ledger.assert_fresh_scalar(name)


def test_the_never_drawn_scalar_stays_unspent():
    import probe_seed_ledger as ledger

    assert 40_000_301 not in ledger.spent_scalars()
    assert 40_000_301 not in ledger.FRESH_PROBE_SCALARS.values()


# ------------------------------------------------ the human narrative

def test_the_document_names_all_defects_and_refuses_to_grade():
    doc = _INCIDENT_MD.read_text(encoding="utf-8")
    # the return-value defect (§2)
    assert "returns `None`" in doc or "returned `None`" in doc
    assert "no claim object to verify" in doc
    # the integration-test gap (§3)
    assert "returning" in doc and "stub" in doc
    # the runner-record provenance defects (§5a-5c)
    assert 'stage: "g3a"' in doc or "reserved_by_stage" in doc     # 5a
    assert "read `completed_gates`, not `stage`" in doc            # 5b
    assert "must not read\n`40,000,401` as the G2 seed" in doc \
        or "must not read `40,000,401` as the G2 seed" in doc      # 5c
    # the JSON is preserved verbatim; corrections live in the doc
    assert "preserved verbatim" in doc
    # it grades nothing
    assert "grades nothing" in doc.lower() or "no scientific" in doc
    assert "reproduction" in doc                      # recovery framing
