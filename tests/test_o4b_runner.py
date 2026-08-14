"""The O4b runner: the frozen order, and what survives a failure.

The campaign itself is 26.2M points and is not approved to run. These
tests drive the same `Campaign` object with a tiny frozen config and a
stub solver, so the ORDER, the accumulation rule, the checkpoints and
the preservation rules are exercised end to end without touching a
campaign seed or the real solver.
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
    """Deterministic and cheap; the real sampler is exercised by the
    O4 suite and is imported, not reimplemented, by the runner."""

    u = rng.random(k)
    return r_lo + u * (r_hi - r_lo), np.full(k, 0.5 * psi)


def _state(r, theta, tol, gap=6.0):
    """A cluster with a comfortable window, so eligibility holds."""

    half = 0.5 * (sz.DT - gap)
    return {"t1": half, "err1": 1e-9, "t2": half, "err2": 1e-9,
            "family1": "one-turn", "family2": "one-turn"}


def _campaign(tmp_path, cfg=None, state_of=_state, **kw):
    return run.Campaign(cfg or _cfg(), _SEEDS, "a" * 40, "b" * 64,
                        checkpoint_path=tmp_path / "ck.json",
                        draw=_draw, state_of=state_of, **kw)


# ------------------------------------------------------ the ordering

def test_the_generator_does_not_exist_until_g3a_has_passed(tmp_path,
                                                           monkeypatch):
    """A G3a failure has spent nothing, and it is the ordering that
    guarantees it rather than a promise in a docstring."""

    campaign = _campaign(tmp_path)
    assert campaign.rng is None

    monkeypatch.setattr(
        run.g3a, "run_preflight",
        lambda tol, eta: {"passed": False,
                          "failed_conditions": ["tri_state_rows"],
                          "conditions": {}})
    with pytest.raises(stages.StageFailure) as caught:
        campaign.run()
    assert caught.value.stage == "g3a"
    assert caught.value.outcome == "INVALID"
    assert caught.value.preserved["fresh_seed_touched"] is False
    assert campaign.rng is None


def test_the_stages_run_in_the_frozen_order(tmp_path, monkeypatch):
    seen = []
    campaign = _campaign(tmp_path)
    for stage in ("run_g3a", "run_g3b", "run_g1", "run_g2"):
        real = getattr(campaign, stage)

        def wrapped(*a, _s=stage, _r=real, **k):
            seen.append(_s)
            return _r(*a, **k)

        monkeypatch.setattr(campaign, stage, wrapped)
    monkeypatch.setattr(campaign, "_check_contract",
                        lambda *a, **k: None)
    campaign.run()
    assert seen == ["run_g3a", "run_g3b", "run_g1", "run_g2"]


# ------------------------------------------- the accumulation rule

def test_every_prefix_point_reaches_the_g1_accumulator(tmp_path,
                                                       monkeypatch):
    """G3b consumes the front of the stream. Points it rejected must
    still be in `V_S1`'s sample."""

    campaign = _campaign(tmp_path)
    monkeypatch.setattr(campaign, "_check_contract",
                        lambda *a, **k: None)
    campaign.run_g3a()
    campaign.run_g3b()
    assert campaign.g1.n == campaign.scanned > 0
    assert campaign.prefix.report()["accumulated_points"] == (
        campaign.g1.n)


def test_ineligible_points_are_accumulated_too(tmp_path, monkeypatch):
    """A window too narrow to probe is still a sample point. Keeping
    only the eligible ones estimates a conditional."""

    def narrow(r, theta, tol):
        return _state(r, theta, tol, gap=-1.0)   # L < 0 -> Z = 0

    campaign = _campaign(tmp_path, state_of=narrow)
    campaign.run_g3a()
    with pytest.raises(stages.StageFailure) as caught:
        campaign.run_g3b()
    assert caught.value.outcome == "INCONCLUSIVE"
    # nothing was eligible, and every point is still in the estimator
    assert campaign.prefix.eligible == 0
    assert campaign.g1.n == campaign.scanned > 0
    assert caught.value.preserved["g1_partial"]["n"] == campaign.g1.n


# ------------------------------------------------------ checkpoints

def test_each_stage_boundary_and_each_g1_chunk_checkpoints(
        tmp_path, monkeypatch):
    written = []
    real = ck.write

    def watching(path, stage, payload):
        written.append(stage)
        return real(path, stage, payload)

    monkeypatch.setattr(ck, "write", watching)
    campaign = _campaign(tmp_path)
    monkeypatch.setattr(campaign, "_check_contract",
                        lambda *a, **k: None)
    campaign.run()
    assert written[0] == "g3b"
    assert written[-1] == "g2_complete"
    assert "g1_complete" in written
    assert written.count("g1_chunk") >= 1

    record = json.loads((tmp_path / "ck.json").read_text(
        encoding="utf-8"))
    assert record["partial"] is True and record["non_verdict"] is True
    assert record["freeze_sha"] == "a" * 40
    assert record["seed"] == 40_000_401
    assert record["rng_position"]["bit_generator"] == "PCG64"


def test_the_checkpoint_is_not_one_of_the_write_once_outputs():
    """It is overwritten by design, and it is a different file from
    both the results and the incident, so a partial cannot be read as
    either."""

    assert run._CHECKPOINT not in run._WRITE_ONCE
    assert set(run._WRITE_ONCE) == {run._ARTIFACT, run._INCIDENT}
    assert run._ARTIFACT != run._INCIDENT


# -------------------------------------------- failure preservation

def test_a_g3b_contract_failure_keeps_the_g1_prefix(tmp_path,
                                                    monkeypatch):
    campaign = _campaign(tmp_path)

    def mismatch(r, theta, state, verdict):
        campaign.mismatches.append({"r": r, "probe": "inside"})

    monkeypatch.setattr(campaign, "_check_contract", mismatch)
    campaign.run_g3a()
    with pytest.raises(stages.StageFailure) as caught:
        campaign.run_g3b()
    failure = caught.value
    assert failure.stage == "g3b"
    assert failure.outcome == "INVALID"
    assert failure.preserved["g1_partial"]["n"] > 0
    assert "no gate has a status" in failure.preserved["is_not_a_verdict"]

    record = stages.incident(failure, {"freeze_sha": "a" * 40})
    assert record["verdict"] is None
    assert record["preserved"]["g1_partial"]["n"] > 0


def test_the_stream_is_continued_not_restarted_after_g3b(tmp_path,
                                                         monkeypatch):
    """The prefix already consumed part of the stream; a fresh
    generator would re-draw points the estimator has taken."""

    campaign = _campaign(tmp_path)
    monkeypatch.setattr(campaign, "_check_contract",
                        lambda *a, **k: None)
    campaign.run_g3a()
    campaign.run_g3b()
    before = campaign.rng.bit_generator.state["state"]["state"]
    campaign.run_g1()
    after = campaign.rng.bit_generator.state["state"]["state"]
    assert after != before
    assert campaign.g1.n == campaign.cfg["n_g1"]


# ---------------------------------------------- the freeze envelope

def test_the_run_demands_the_exact_approved_sha():
    state = {"rev": "a" * 40, "dirty": False}
    run.verify_rev("preflight", "A" * 40, state)          # case-fold
    with pytest.raises(SystemExit, match="40-hex commit"):
        run.verify_rev("preflight", "abc", state)
    with pytest.raises(SystemExit, match="the approval does not name"):
        run.verify_rev("preflight", "b" * 40, state)


def test_the_o4_reservation_is_retained_and_is_not_this_ref():
    import o4b_reservation as res

    assert res.REF == "refs/o4b/reservation"
    assert res.RETAINED_REF == "refs/o4/reservation"
    assert res.RETAINED_OBJECT == (
        "c4da1626463e6a6505374813cf3f56d6b429c209")
    source = (_REPO / "experiments" / "oracle"
              / "o4b_reservation.py").read_text(encoding="utf-8")
    # the retained ref is only ever READ
    assert 'REMOTE, f"{obj}:{RETAINED_REF}"' not in source
    assert "RETAINED_REF" in source


def test_a_missing_o4_reservation_stops_the_run(monkeypatch):
    import o4b_reservation as res

    monkeypatch.setattr(res, "_ls_remote", lambda ref: None)
    with pytest.raises(SystemExit, match="is GONE"):
        res.verify_o4_ref_retained()
    monkeypatch.setattr(res, "_ls_remote", lambda ref: "d" * 40)
    with pytest.raises(SystemExit, match="has been rewritten"):
        res.verify_o4_ref_retained()


def test_the_frozen_config_is_inherited_not_re_frozen():
    assert run.FROZEN["tau"] == sz.TAU
    assert run.FROZEN["n_g1"] == sz.N_G1
    import o4_g3_redesign as g3
    assert run.FROZEN["eta"] == g3.ETA
    assert run.FROZEN["n_avail"] == g3.N_AVAIL
    assert run.FROZEN["k_g3b"] == g3.K_G3B


def test_the_committed_manifest_is_the_one_this_tree_produces():
    """Regenerable, not a one-off script's output. A surface file that
    quietly stopped being pinned would otherwise be invisible."""

    committed = json.loads(run._MANIFEST.read_text(encoding="utf-8"))
    assert committed == run.build_manifest()
    assert set(committed["files"]) == set(run.PROTOCOL_SURFACE)


def test_the_manifest_pins_the_inherited_o4_modules_and_the_instrument():
    """O4b reads the estimand, oracle, tau and sampler from the O4
    freeze, so an edit to either is a different campaign."""

    pinned = set(run.PROTOCOL_SURFACE)
    for rel in ("experiments/oracle/o4_sizing.py",
                "experiments/oracle/o4_volume_audit.py",
                "experiments/positive_control/s1_schwarzschild_cost.py",
                "experiments/positive_control/probe_seed_ledger.py",
                "docs/prereg/p14_o4b_sizing.json",
                "docs/prereg/p14_o3_volume.json"):
        assert rel in pinned, rel
    # every O4b module the run actually executes
    for rel in run.PROTOCOL_SURFACE:
        assert (_REPO / rel).exists(), rel
    o4b_modules = {p.name for p in (_REPO / "experiments"
                                    / "oracle").glob("o4b_*.py")}
    assert o4b_modules == {Path(r).name for r in pinned
                           if Path(r).name.startswith("o4b_")}


def test_the_manifest_admits_it_cannot_certify_itself():
    committed = json.loads(run._MANIFEST.read_text(encoding="utf-8"))
    assert "cannot certify itself" in committed["note"]
    assert "exact approved freeze SHA" in committed["note"]
    assert committed["environment"] == run.environment()


def test_the_sampler_is_imported_rather_than_copied():
    """A second copy could drift into being a different instrument
    while still passing every digest check."""

    source = (_REPO / "experiments" / "oracle"
              / "o4b_volume_audit.py").read_text(encoding="utf-8")
    assert "o4._draw" in source
    assert "np.cbrt" not in source          # not reimplemented here
