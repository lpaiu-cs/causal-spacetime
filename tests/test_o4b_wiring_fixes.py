"""The O4b abort's wiring defects, fixed and pinned.

The campaign completed both gates and then failed on a one-line exit
provenance bug: main's claim callback stashed the reservation object and
returned None, so Campaign.run carried None to the pre-publish re-verify.
Every prior success-path test drove Campaign.run with a RETURNING stub,
so the real main -> claim -> Campaign.run wiring was never run end to end.
These tests close that gap and pin the related record-provenance fixes:
the budget stage label, the incident failure_point, and the checkpoint's
named RNG stream.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import o4_sizing as sz  # noqa: E402
import o4b_checkpoint as ck  # noqa: E402
import o4b_stages as stages  # noqa: E402
import o4b_volume_audit as run  # noqa: E402

_SEEDS = {"o4b_g1_audit": 40_000_401, "o4b_g2_leakage": 40_000_411}


def _cfg(**over) -> dict:
    small = dict(run.FROZEN)
    small.update(n_g1=64, n_g2=16, n_avail=4, k_g3b=2, chunk=8)
    small.update(over)
    return small


def _draw(rng, k, r_lo, r_hi, psi):
    u = rng.random(k)
    return r_lo + u * (r_hi - r_lo), np.full(k, 0.5 * psi)


def _state(r, theta, tol, gap=6.0):
    half = 0.5 * (sz.DT - gap)
    return {"t1": half, "err1": 1e-9, "t2": half, "err2": 1e-9,
            "theta": theta, "dpsi": theta, "recovery_shift": 0.0,
            "family1": "one-turn", "family2": "one-turn"}


def _ell(r, theta, tol, gap=6.0):
    half = 0.5 * (sz.DT - gap)
    ell = sz.DT - 2 * half
    return (ell if ell > 0.0 else 0.0), half, half


def _campaign(tmp_path, **kw):
    return run.Campaign(_cfg(), _SEEDS, "a" * 40, "b" * 64,
                        checkpoint_path=tmp_path / "ck.json",
                        draw=_draw, state_of=_state, ell=_ell, **kw)


# ============================ item 1 + 2: the claim object wiring, via main

def test_main_wires_the_claim_object_through_to_the_publish_check(
        tmp_path, monkeypatch):
    """THE REGRESSION. main()'s real claim closure must return the
    reservation object so Campaign.run carries it to the pre-publish
    re-verify. This drives the actual main path -- not a returning stub
    -- so a claim that stashed and returned None (the abort) is caught."""

    claimed_obj = "d" * 40
    seen = {}

    monkeypatch.setattr(run, "preflight", lambda rev: {
        "seeds": _SEEDS, "git": {"rev": "a" * 40}, "manifest": {}})
    monkeypatch.setattr(run, "_sha256", lambda p: "b" * 64)
    monkeypatch.setattr(run.reservation, "claim",
                        lambda payload: claimed_obj)

    def verify_still_held(obj):
        seen["verified"] = obj                 # what the check received
        return obj

    monkeypatch.setattr(run.reservation, "verify_still_held",
                        verify_still_held)
    monkeypatch.setattr(run.Campaign, "run_g3a",
                        lambda self: {"passed": True, "conditions": {}})
    monkeypatch.setattr(run.Campaign, "run_g3b",
                        lambda self: {"availability": "complete"})
    monkeypatch.setattr(run.Campaign, "run_g1",
                        lambda self: {"status": "concordant", "n": 1})
    monkeypatch.setattr(run.Campaign, "run_g2",
                        lambda self: {"status": "concordant", "n": 1})
    monkeypatch.setattr(run, "_ARTIFACT", tmp_path / "results.json")
    monkeypatch.setattr(run, "_INCIDENT", tmp_path / "incident.json")
    monkeypatch.setattr(run, "_WRITE_ONCE",
                        (tmp_path / "results.json",
                         tmp_path / "incident.json"))
    monkeypatch.setattr(sys, "argv",
                        ["prog", "--freeze-rev", "a" * 40])

    run.main()

    # the object claim() returned actually reached the pre-publish check
    assert seen["verified"] == claimed_obj
    # a result was published carrying it, and NO incident was filed
    result = json.loads((tmp_path / "results.json").read_text("utf-8"))
    assert result["kind"] == "results"
    assert result["reservation"]["object"] == claimed_obj
    assert not (tmp_path / "incident.json").exists()


def test_the_real_claim_closure_returns_the_object_it_stashes(
        tmp_path, monkeypatch):
    """Directly: the closure main builds records the object AND returns
    it. The abort was a closure that did only the first."""

    monkeypatch.setattr(run, "preflight", lambda rev: {
        "seeds": _SEEDS, "git": {"rev": "a" * 40}, "manifest": {}})
    monkeypatch.setattr(run, "_sha256", lambda p: "b" * 64)
    monkeypatch.setattr(run.reservation, "claim", lambda payload: "e" * 40)

    captured = {}

    def fake_run(self, on_g3a_passed=None, on_before_publish=None):
        captured["returned"] = on_g3a_passed()      # the real claim()
        return {"kind": "results", "reservation": {"object": "e" * 40},
                "g1": {"status": "concordant"},
                "g2": {"status": "concordant"}}

    monkeypatch.setattr(run.Campaign, "run", fake_run)
    monkeypatch.setattr(run, "publish_write_once",
                        lambda *a, **k: True)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", "a" * 40])

    run.main()
    assert captured["returned"] == "e" * 40         # not None


# ==================================== item 3: the budget stage label moves

def test_the_budget_stage_label_follows_the_run(tmp_path, monkeypatch):
    """The budget's stage label is advanced at each stage, so calls are
    attributed to the stage that spent them rather than all to g3a
    (incident defect 5a). Recorded directly, since the stub gates here
    do not call the solver and so reserve nothing to compare by count."""

    campaign = _campaign(tmp_path)
    monkeypatch.setattr(campaign, "_check_contract", lambda *a, **k: None)

    entered = []
    real_enter = campaign.budget.enter

    def record(stage):
        entered.append(stage)
        real_enter(stage)

    monkeypatch.setattr(campaign.budget, "enter", record)
    campaign.run()

    assert entered[0] == "g3a"                 # starts where it started
    assert {"g3b", "g1", "g2"} <= set(entered)  # and moves on, per stage
    assert campaign.budget.stage == "g2"       # ended in the last stage


def test_the_budget_label_charges_reserved_calls_to_the_current_stage():
    """The mechanism under the fix: enter() then reserve() attributes to
    the entered stage, not to g3a."""

    budget = run.o4b_budget.Budget(max_calls=10_000, max_wall_s=1_000.0)
    budget.reserve(5)                          # starts at g3a
    budget.enter("g1")
    budget.reserve(7)
    by_stage = budget.state()["reserved_by_stage"]
    assert by_stage["g3a"] == 5
    assert by_stage["g1"] == 7


# ============================== item 4: the incident names the failure_point

def test_a_pre_publish_failure_is_labelled_pre_publish_not_g2(
        tmp_path, monkeypatch):
    """The stop is under stage g2 but is NOT a claim that G2 is
    unfinished; failure_point records the actual step (incident
    defect 5b)."""

    campaign = _campaign(tmp_path)
    monkeypatch.setattr(campaign, "_check_contract", lambda *a, **k: None)

    def boom(obj):
        raise RuntimeError("exit check blew up")

    with pytest.raises(stages.StageFailure) as caught:
        campaign.run(on_g3a_passed=lambda: "c" * 40,
                     on_before_publish=boom)
    assert caught.value.stage == "g2"
    assert caught.value.failure_point == "pre_publish_verify"
    record = stages.incident(caught.value, {})
    assert record["failure_point"] == "pre_publish_verify"
    assert record["stage"] == "g2"


def test_a_publication_failure_is_labelled_artifact_publication(
        tmp_path, monkeypatch):
    """A post-gate failure -- both gates finished and only the publish
    step failed -- is recorded as `artifact_publication`, not defaulted
    to a `g2` gate failure that would contradict `detail.at` (review
    PR #78 R1)."""

    monkeypatch.setattr(run, "preflight", lambda rev: {
        "seeds": _SEEDS, "git": {"rev": "a" * 40}, "manifest": {}})
    monkeypatch.setattr(run, "_sha256", lambda p: "b" * 64)
    monkeypatch.setattr(run.reservation, "claim", lambda payload: "d" * 40)
    monkeypatch.setattr(run.reservation, "verify_still_held",
                        lambda obj: obj)
    monkeypatch.setattr(run.Campaign, "run_g3a",
                        lambda self: {"passed": True, "conditions": {}})
    monkeypatch.setattr(run.Campaign, "run_g3b",
                        lambda self: {"availability": "complete"})
    monkeypatch.setattr(run.Campaign, "run_g1",
                        lambda self: {"status": "concordant", "n": 1})
    monkeypatch.setattr(run.Campaign, "run_g2",
                        lambda self: {"status": "concordant", "n": 1})
    incident_path = tmp_path / "incident.json"
    monkeypatch.setattr(run, "_ARTIFACT", tmp_path / "results.json")
    monkeypatch.setattr(run, "_INCIDENT", incident_path)
    monkeypatch.setattr(run, "_WRITE_ONCE",
                        (tmp_path / "results.json", incident_path))

    real_publish = run.publish_write_once

    def boom_on_artifact(path, payload, receipt=None):
        if path == run._ARTIFACT:               # publish fails, uncommitted
            raise OSError("disk full")
        return real_publish(path, payload, receipt)   # the incident write

    monkeypatch.setattr(run, "publish_write_once", boom_on_artifact)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", "a" * 40])

    with pytest.raises(SystemExit):
        run.main()

    record = json.loads(incident_path.read_text("utf-8"))
    assert record["failure_point"] == "artifact_publication"
    assert record["detail"]["at"] == "artifact publication"
    assert record["stage"] == "g2"              # the run-stage is still g2
    assert record["verdict"] is None


def test_a_reservation_claim_failure_is_labelled_as_such(
        tmp_path, monkeypatch):
    campaign = _campaign(tmp_path)

    def boom():
        raise RuntimeError("push failed")

    with pytest.raises(stages.StageFailure) as caught:
        campaign.run(on_g3a_passed=boom)
    assert caught.value.stage == "g3b"
    assert caught.value.failure_point == "reservation_claim"


def test_an_ordinary_stage_failure_points_at_its_stage(tmp_path,
                                                       monkeypatch):
    """failure_point defaults to the stage, so a real G1 failure is not
    mislabelled."""

    campaign = _campaign(tmp_path)
    monkeypatch.setattr(campaign, "_check_contract", lambda *a, **k: None)

    def blow_up_g1():
        raise RuntimeError("g1 numerical failure")

    monkeypatch.setattr(campaign, "run_g1", blow_up_g1)
    with pytest.raises(stages.StageFailure) as caught:
        campaign.run(on_g3a_passed=lambda: "c" * 40)
    assert caught.value.stage == "g1"
    assert caught.value.failure_point == "g1"


# ============================ item 5: the checkpoint names its RNG stream

def test_the_runner_checkpoint_names_the_g1_stream(tmp_path, monkeypatch):
    """Every checkpoint the runner writes serialises the G1 stream and
    says so, so a g2_complete record is not read as carrying G2's
    position (incident defect 5c)."""

    campaign = _campaign(tmp_path)
    monkeypatch.setattr(campaign, "_check_contract", lambda *a, **k: None)
    campaign.run()
    record = ck.read(tmp_path / "ck.json")
    assert record["rng_stream"] == "o4b_g1_audit"
    assert record["seed"] == 40_000_401
