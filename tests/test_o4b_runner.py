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
    """The PREDICATE's view: a cluster with a comfortable window, so
    eligibility holds."""

    half = 0.5 * (sz.DT - gap)
    return {"t1": half, "err1": 1e-9, "t2": half, "err2": 1e-9,
            "theta": theta, "dpsi": theta, "recovery_shift": 0.0,
            "family1": "one-turn", "family2": "one-turn"}


def _ell(r, theta, tol, gap=6.0):
    """The ESTIMAND's view, at the drawn theta. Deliberately a
    different function from `_state`: the runner must not substitute
    one coordinate for the other."""

    half = 0.5 * (sz.DT - gap)
    ell = sz.DT - 2 * half
    return (ell if ell > 0.0 else 0.0), half, half


def _campaign(tmp_path, cfg=None, state_of=_state, ell=_ell, **kw):
    return run.Campaign(cfg or _cfg(), _SEEDS, "a" * 40, "b" * 64,
                        checkpoint_path=tmp_path / "ck.json",
                        draw=_draw, state_of=state_of, ell=ell, **kw)


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

    def narrow_ell(r, theta, tol):
        return _ell(r, theta, tol, gap=-1.0)

    campaign = _campaign(tmp_path, state_of=narrow, ell=narrow_ell)
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


# ------------------------------------------------- the two coordinates

def test_the_predicate_state_is_built_at_the_recovered_dpsi():
    """The abort's third non-adaptive margin: a probe placed from
    `T(theta)` and judged at `t_min(acos(cos theta))` is placed in one
    coordinate and read in another. The recovery error is not bounded
    in ulps near zero -- it reaches 4.14e-13 at theta = 1e-5, the same
    order as eta."""

    import o4_g3_redesign as g3
    import s1_schwarzschild_cost as s1

    theta = 1e-5
    state = run.solver_state(14.0, theta, s1.DEFAULT_TOL)
    assert state["dpsi"] == g3.wrapper_dpsi(theta) != theta
    assert abs(state["recovery_shift"]) > 0.0

    # the shift is the same order as eta, which is why placing in one
    # coordinate and reading in the other is not a rounding detail
    assert abs(state["recovery_shift"]) == pytest.approx(4.137e-13,
                                                         rel=1e-3)
    assert abs(state["recovery_shift"]) > 0.4 * g3.ETA

    # and the state IS the solver at the recovered angle
    assert state["t1"] == s1.flight_time(
        sz.R_IN, 14.0, state["dpsi"], s1.M, s1.DEFAULT_TOL)[0]
    assert state["t2"] == s1.flight_time(
        14.0, sz.R_OUT, state["dpsi"], s1.M, s1.DEFAULT_TOL)[0]


def test_g1_uses_the_estimand_and_g3b_uses_the_predicate_state(
        tmp_path, monkeypatch):
    """Two coordinates, two functions. Fixing `solver_state` alone
    would have moved the G1 sample off the frozen estimand instead."""

    called = {"ell": 0, "state": 0}

    def counting_ell(r, theta, tol):
        called["ell"] += 1
        return _ell(r, theta, tol)

    def counting_state(r, theta, tol):
        called["state"] += 1
        return _state(r, theta, tol)

    campaign = _campaign(tmp_path, ell=counting_ell,
                         state_of=counting_state)
    monkeypatch.setattr(campaign, "_check_contract",
                        lambda *a, **k: None)
    campaign.run_g3a()
    campaign.run_g3b()
    # every drawn point costs an estimand evaluation; only candidates
    # cost a predicate-coordinate one
    assert called["ell"] == campaign.scanned
    assert called["state"] == campaign.prefix.candidates + (
        campaign.prefix.scan_candidates)
    assert called["state"] <= called["ell"]


# ------------------------------------------------ the probe contract

def test_every_probe_checks_both_legs():
    """The design fixes both legs of all three probes, and the cost
    model froze six predicate calls per fully-testable cluster on that
    basis. Checking only the starved leg spends four and lets a defect
    on the other one through."""

    assert set(run.Campaign.CONTRACT) == {"inside", "outside_above",
                                          "outside_below"}
    for wanted in run.Campaign.CONTRACT.values():
        assert set(wanted) == {"p_to_x", "x_to_q"}
    assert run.Campaign.CONTRACT["outside_above"] == {"p_to_x": True,
                                                      "x_to_q": False}
    assert run.Campaign.CONTRACT["outside_below"] == {"p_to_x": False,
                                                      "x_to_q": True}
    calls = 2 * len(run.Campaign.CONTRACT)
    import o4b_sizing as sizing
    assert calls == sizing.PROBE_CALLS_PER_CLUSTER == 6


def test_a_defect_on_the_unstarved_leg_is_a_mismatch(tmp_path,
                                                     monkeypatch):
    """outside-above starves `x -> q`; a wrapper returning `None` on
    the other leg used never to be looked at."""

    campaign = _campaign(tmp_path)
    monkeypatch.setattr(
        campaign, "_legs",
        lambda r, th, t_x: {"p_to_x": None, "x_to_q": False})
    campaign._check_contract(
        13.0, 0.1, _state(13.0, 0.1, 1e-8),
        {"probes": {"inside": {"reached": True, "t_x": 1.0},
                    "outside_above": {"reached": True, "dt": 1.0},
                    "outside_below": {"reached": True, "dt": 1.0}}})
    legs = {(m["probe"], m["leg"]) for m in campaign.mismatches}
    assert ("outside_above", "p_to_x") in legs
    assert all(m["got"] == "None" or m["want"] is not False
               for m in campaign.mismatches
               if m["leg"] == "p_to_x")


# ------------------------------------------------------- the cap path

def test_a_fired_cap_becomes_an_inconclusive_incident(tmp_path,
                                                      monkeypatch):
    """Not a `CapReached` escaping the entry point. If it fires during
    G3b -- before the first checkpoint exists -- the seeds are already
    spent and the reservation forbids re-running them, so everything
    accumulated would vanish with the process."""

    import o4b_budget as bud

    campaign = _campaign(tmp_path)
    campaign.run_g3a()

    def fires():
        raise bud.CapReached("max-wall", 8_000_000, 7_999_999,
                             86_401.0, 1)

    with pytest.raises(stages.StageFailure) as caught:
        campaign._staged("g3b", fires)
    failure = caught.value
    assert failure.stage == "g3b"
    assert failure.outcome == "INCONCLUSIVE"
    assert failure.detail["cap"]["reason"] == "max-wall"
    assert "g1_partial" in failure.preserved
    assert failure.preserved["budget"]["reserved"] >= 0
    # the incident is the record here; see
    # `test_a_g3b_cap_writes_no_checkpoint_stage_it_did_not_reach`
    # for why no checkpoint is written at this point


def test_a_cap_is_inconclusive_and_never_invalid():
    """The instrument said nothing wrong; the run ran out. `INVALID`
    would indict the wrapper for a budget decision."""

    source = (_REPO / "experiments" / "oracle"
              / "o4b_volume_audit.py").read_text(encoding="utf-8")
    assert "a cap is not a contract failure" in source.lower()
    assert "CapReached" in source


# ------------------------------------------------ the reservation nonce

def test_two_attempts_in_one_second_claim_different_objects(
        monkeypatch):
    """Git's commit timestamp is one-second resolution and the payload
    and identity are fixed, so without a nonce two runs on one host
    could build the SAME commit -- and then the second push is
    `everything up-to-date` and both pass the `held() == obj` check."""

    import o4b_reservation as res

    messages = []
    monkeypatch.setattr(res, "verify_o4_ref_retained", lambda: "c" * 40)
    monkeypatch.setattr(res, "held", lambda: None)
    monkeypatch.setattr(res, "_make_commit",
                        lambda m: messages.append(m) or "obj")
    monkeypatch.setattr(res.subprocess, "run",
                        lambda *a, **k: type("R", (), {
                            "returncode": 0, "stdout": "",
                            "stderr": ""})())
    # free before the push, taken after it
    pushed = {"n": 0}

    def held():
        pushed["n"] += 1
        return None if pushed["n"] % 2 == 1 else "obj"

    monkeypatch.setattr(res, "held", held)

    payload = {"campaign": "o4b", "freeze_rev": "a" * 40}
    res.claim(dict(payload))
    res.claim(dict(payload))
    assert len(messages) == 2
    assert messages[0] != messages[1]
    assert all("attempt_nonce" in m for m in messages)


# --------------------------------------- the batched draw's boundary

def test_the_chunk_position_is_recorded_at_both_ends(tmp_path,
                                                     monkeypatch):
    """`rng_position` is the END of the last chunk drawn; the
    accumulator holds only the consumed prefix of it. Resuming from
    the position alone would skip the unconsumed tail and change the
    frozen stream."""

    campaign = _campaign(tmp_path)
    monkeypatch.setattr(campaign, "_check_contract",
                        lambda *a, **k: None)
    campaign.run_g3a()
    campaign.run_g3b()

    chunk = ck.read(tmp_path / "ck.json")["statistics"]["chunk"]
    assert chunk["state_before_draw"]["bit_generator"] == "PCG64"
    assert chunk["state_before_draw"] != (
        ck.read(tmp_path / "ck.json")["rng_position"])
    assert 0 < chunk["consumed"] <= chunk["size"]
    assert "skip the unconsumed" in chunk["why"]


def test_the_scan_stops_at_the_stopping_time_not_the_chunk_end(
        tmp_path, monkeypatch):
    """Otherwise a mismatch found among points AFTER the sample the
    rule defined could turn an already complete G3b INVALID."""

    campaign = _campaign(tmp_path, cfg=_cfg(n_avail=2, k_g3b=2,
                                            chunk=64))
    checked = []
    monkeypatch.setattr(
        campaign, "_check_contract",
        lambda r, th, st, v: checked.append(r))
    campaign.run_g3a()
    campaign.run_g3b()
    # exactly the stopping sample, not the whole 64-point buffer
    assert campaign.scanned == 2
    assert len(checked) == 2
    assert campaign.chunk_consumed == 2 < campaign.chunk_size == 64


# ------------------------------------------ every exit leaves a record

def test_a_solver_systemexit_becomes_an_abort_incident(tmp_path):
    """The fail-closed solver paths raise `SystemExit`, which is a
    BaseException and slips past `except Exception` entirely. Without
    conversion the process dies with the seeds spent and no incident
    -- the observability defect this freeze exists to close."""

    campaign = _campaign(tmp_path)
    campaign.run_g3a()

    def fails():
        raise SystemExit("fail-closed: flight_time returned nan")

    with pytest.raises(stages.StageFailure) as caught:
        campaign._staged("g3b", fails)
    failure = caught.value
    assert failure.outcome == "ABORT"
    assert failure.detail["exception"]["type"] == "SystemExit"
    assert "fail-closed" in failure.detail["exception"]["message"]
    assert "g1_partial" in failure.preserved


def test_an_ordinary_exception_becomes_an_abort_incident(tmp_path):
    campaign = _campaign(tmp_path)
    campaign.run_g3a()
    with pytest.raises(stages.StageFailure) as caught:
        campaign._staged("g1", lambda: 1 / 0)
    assert caught.value.outcome == "ABORT"
    assert caught.value.detail["exception"]["type"] == (
        "ZeroDivisionError")


def test_a_stage_failure_passes_through_unchanged(tmp_path):
    """It is already the right shape; re-wrapping it as ABORT would
    turn an INVALID contract failure into an unexplained crash."""

    campaign = _campaign(tmp_path)
    original = stages.StageFailure("g3b", "INVALID", "mismatch")

    def raises():
        raise original

    with pytest.raises(stages.StageFailure) as caught:
        campaign._staged("g3b", raises)
    assert caught.value is original


def test_a_keyboard_interrupt_is_not_swallowed(tmp_path):
    campaign = _campaign(tmp_path)
    with pytest.raises(KeyboardInterrupt):
        campaign._staged("g1", lambda: (_ for _ in ()).throw(
            KeyboardInterrupt()))


# ------------------------------- the cap must not delete a result

def test_a_g2_cap_does_not_overwrite_the_completed_g1_checkpoint(
        tmp_path, monkeypatch):
    """Entering G2 means `g1_complete` is on disk with G1's status and
    interval, and this exit publishes no results artifact. Overwriting
    it would delete a G1 result the run had already established."""

    import o4b_budget as bud

    campaign = _campaign(tmp_path)
    monkeypatch.setattr(campaign, "_check_contract",
                        lambda *a, **k: None)
    campaign.run_g3a()
    campaign.run_g3b()
    campaign.run_g1()
    assert ck.read(tmp_path / "ck.json")["stage"] == "g1_complete"
    before = ck.read(tmp_path / "ck.json")

    def fires():
        raise bud.CapReached("max-calls", 8 * 10 ** 7, 8 * 10 ** 7,
                             100.0, 1)

    with pytest.raises(stages.StageFailure) as caught:
        campaign._staged("g2", fires)
    assert ck.read(tmp_path / "ck.json") == before
    assert caught.value.preserved["g1_partial"]["n"] > 0


def test_a_g3b_cap_writes_no_checkpoint_stage_it_did_not_reach(
        tmp_path):
    """A stop inside G3b has reached no checkpoint stage; writing
    `g3b` would claim G3b finished. The incident is the record."""

    import o4b_budget as bud

    campaign = _campaign(tmp_path)
    campaign.run_g3a()
    with pytest.raises(stages.StageFailure) as caught:
        campaign._staged("g3b", lambda: (_ for _ in ()).throw(
            bud.CapReached("max-wall", 1, 1, 9.0, 1)))
    assert not (tmp_path / "ck.json").exists()
    assert caught.value.preserved["availability"]["candidates"] == 0


# ------------------------------- the tail G3b drew but did not use

def test_g1_continues_from_the_tail_g3b_left(tmp_path, monkeypatch):
    """The RNG has already advanced past that tail. Discarding it
    would make G1 skip a stretch of the frozen stream and integrate a
    different sample -- with `n_avail=2, k_g3b=2, chunk=64`, points
    3..64 would simply vanish."""

    campaign = _campaign(tmp_path, cfg=_cfg(n_avail=2, k_g3b=2,
                                            chunk=64, n_g1=64))
    monkeypatch.setattr(campaign, "_check_contract",
                        lambda *a, **k: None)
    campaign.run_g3a()
    campaign.run_g3b()
    assert campaign.scanned == 2
    assert len(campaign.pending) == 62      # drawn, not yet used

    tail = list(campaign.pending)
    campaign.run_g1()
    assert campaign.pending == []
    assert campaign.g1.n == 64
    # and the same stream: G3b's 2 + the 62 it handed over + 0 new
    assert campaign.scanned == 64
    assert tail[0] is not None


def test_the_whole_frozen_prefix_is_the_same_points_either_way(
        tmp_path, monkeypatch):
    """The estimator's sample must not depend on where G3b happened to
    stop. Run the same seed with two different stopping points and the
    G1 sample has to be identical."""

    seen = {}

    def record(key):
        def ell(r, theta, tol):
            seen.setdefault(key, []).append((r, theta))
            return _ell(r, theta, tol)
        return ell

    for key, k_g3b in (("early", 2), ("late", 6)):
        campaign = _campaign(
            tmp_path, cfg=_cfg(n_avail=k_g3b, k_g3b=k_g3b, chunk=64,
                               n_g1=64), ell=record(key))
        monkeypatch.setattr(campaign, "_check_contract",
                            lambda *a, **k: None)
        campaign.run_g3a()
        campaign.run_g3b()
        campaign.run_g1()
        assert campaign.g1.n == 64
    assert seen["early"][:64] == seen["late"][:64]


# --------------------------- one point, one commit boundary

def test_a_cap_inside_judge_leaves_the_cursor_and_g1_agreeing(
        tmp_path):
    """`observe` adds `z` and advances the cursor together, so an
    incident can never report a G1 sample containing a point the
    cursor says was never reached."""

    import o4b_budget as bud

    campaign = _campaign(tmp_path, cfg=_cfg(n_avail=8, k_g3b=8))
    campaign.run_g3a()

    def exploding(r, theta, tol):
        if campaign.g1.n >= 3:
            raise bud.CapReached("max-calls", 10, 10, 1.0, 1)
        return _state(r, theta, tol)

    campaign.state_of = exploding
    with pytest.raises(stages.StageFailure) as caught:
        campaign._staged("g3b", campaign.run_g3b)
    preserved = caught.value.preserved
    assert preserved["g1_partial"]["n"] == preserved["scanned_points"]
    assert campaign.g1.n == campaign.scanned


def test_an_unjudged_candidate_is_not_counted(tmp_path):
    """Otherwise a cap at the `N_avail`-th candidate could report
    `complete = True` and publish availability rates over a candidate
    that was never judged."""

    import o4b_budget as bud

    campaign = _campaign(tmp_path, cfg=_cfg(n_avail=3, k_g3b=3))
    campaign.run_g3a()

    def exploding(r, theta, tol):
        if campaign.prefix.candidates == 2:
            raise bud.CapReached("max-calls", 10, 10, 1.0, 1)
        return _state(r, theta, tol)

    campaign.state_of = exploding
    with pytest.raises(stages.StageFailure) as caught:
        campaign._staged("g3b", campaign.run_g3b)
    report = caught.value.preserved["availability"]
    assert report["candidates"] == 2          # not 3
    assert report["complete"] is False
    assert "rates_withheld" in report


# ------------------ a contract failure outranks a later interruption

def test_a_mismatch_seen_before_a_cap_is_INVALID_not_INCONCLUSIVE(
        tmp_path):
    """`INCONCLUSIVE` indicts nothing; `INVALID` indicts the
    instrument. Having seen it disagree, the second is the true
    sentence and the cap only decided when the run stopped."""

    import o4b_budget as bud

    campaign = _campaign(tmp_path)
    campaign.run_g3a()
    campaign.mismatches.append({"probe": "inside", "leg": "p_to_x"})

    def fires():
        raise bud.CapReached("max-calls", 8 * 10 ** 7, 8 * 10 ** 7,
                             100.0, 1)

    with pytest.raises(stages.StageFailure) as caught:
        campaign._staged("g3b", fires)
    failure = caught.value
    assert failure.outcome == "INVALID"
    assert failure.detail["mismatches"]
    assert failure.detail["stopped_by"]["type"] == "CapReached"
    assert "observed first" in failure.detail["why_not_inconclusive"]


def test_a_mismatch_seen_before_a_crash_is_also_INVALID(tmp_path):
    campaign = _campaign(tmp_path)
    campaign.run_g3a()
    campaign.mismatches.append({"probe": "outside_above"})
    with pytest.raises(stages.StageFailure) as caught:
        campaign._staged("g3b", lambda: 1 / 0)
    assert caught.value.outcome == "INVALID"
    assert caught.value.detail["stopped_by"]["type"] == (
        "ZeroDivisionError")
