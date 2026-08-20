"""S6 M=3.0 Poisson-count campaign freeze contract tests.

The frozen ruled values re-derived from the runner's OWN 96-bit rule,
cross-checked against the frozen P14 P2 engine; the acceptance window
proven equivalent to the verdict rule at U = 0; the power sentence's
numbers pinned; the import gate, the SHA gates, the preflight
refusals, the reservation ordering and retained set, the chunked
atomic checkpoints, and the REAL `main()` wiring end to end -- the
O4b lesson: no returning stubs standing in for the actual wiring.
NO RESULTS: execution waits on the integrated approval."""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import p14_probe_p2 as p2  # noqa: E402
import s6_m30_count as camp  # noqa: E402

_SHA = "a" * 40


# ------------------------------------------------ the frozen values

def test_frozen_configuration_is_the_ruled_one():
    assert camp.FROZEN == {
        "a_intensity": 440.0,
        "tau": 0.025,
        "alpha": 0.05,
        "v_lo": 121.22502058308072,
        "v_hi": 122.44328048247328,
        "v_ref": 121.834150532777,
        "u_power_ceiling": 30,
        "seed_name": "s6_m30_count",
        "tol": 1e-08,
        "chunk": 65_536,
        "scan_call_cap": 11_000_000,
        "max_calls": 11_006_466,
        "max_wall_s": 86_400.0,
    }
    assert camp.G3A_PREFLIGHT_CALLS == 6_466
    assert camp.FROZEN["max_calls"] == (
        camp.FROZEN["scan_call_cap"] + camp.G3A_PREFLIGHT_CALLS)
    assert camp.ACCEPTANCE == (52_986, 54_220)


def test_the_frozen_endpoints_are_the_o3p_artifact_verbatim():
    """The import gate's claim, asserted independently: the frozen
    literals equal the serialized artifact values, and the pilot
    artifact certifies the frozen U ceiling."""

    o3p = json.loads((_REPO / "docs" / "prereg"
                      / "p14_s6_m30_volume.json").read_text(
                          encoding="utf-8"))
    assert o3p["result"]["v_lo"] == camp.FROZEN["v_lo"]
    assert o3p["result"]["v_hi"] == camp.FROZEN["v_hi"]
    assert o3p["result"]["status"] == "target-met"
    assert (o3p["v_ref_rung_recommendation"]["value"]
            == camp.FROZEN["v_ref"])
    pilot = json.loads((_REPO / "docs" / "prereg"
                        / "p14_s6_m30_pilot.json").read_text(
                            encoding="utf-8"))
    assert pilot["decision"]["verdict"] == "FEASIBLE"
    assert pilot["decision"]["u_max"] == camp.FROZEN["u_power_ceiling"]


def test_the_seed_is_fresh_and_the_pilot_stream_stays_retired():
    import probe_seed_ledger as ledger

    assert "s6_m30_count" not in ledger.FRESH_PROBE_SCALARS
    assert ledger.OBSERVED_PROBE_SCALARS["s6_m30_count"] == 40_000_511
    assert ledger.replay_scalar("s6_m30_count") == 40_000_511
    assert ledger.OBSERVED_PROBE_SCALARS["s6_m30_pilot"] == 40_000_501
    assert 40_000_301 not in ledger.spent_scalars()
    assert 40_000_301 not in ledger.FRESH_PROBE_SCALARS.values()


def test_the_acceptance_window_is_the_verdict_rule_at_zero_u():
    """decide(k, 0) is CONCORDANT exactly on the frozen window --
    the sizing's integers and the runner's 96-bit rule are one."""

    k_lo, k_hi = camp.ACCEPTANCE
    assert camp.decide(k_lo - 1, 0)["verdict"] == "INCONCLUSIVE"
    assert camp.decide(k_lo, 0)["verdict"] == "CONCORDANT"
    assert camp.decide(k_hi, 0)["verdict"] == "CONCORDANT"
    assert camp.decide(k_hi + 1, 0)["verdict"] == "INCONCLUSIVE"


def test_ambiguity_shrinks_the_window_from_the_top_exactly():
    """With u ambiguous points the upper bound reads K_certain + u,
    so concordance needs K_certain + u <= k_hi -- checked at the
    exact integer boundary."""

    k_lo, k_hi = camp.ACCEPTANCE
    assert camp.decide(k_hi - 30, 30)["verdict"] == "CONCORDANT"
    assert camp.decide(k_hi - 29, 30)["verdict"] == "INCONCLUSIVE"
    assert camp.decide(k_lo, 30)["verdict"] == "CONCORDANT"


def test_the_three_way_verdict_is_exhaustive_and_disjoint():
    cases = {
        (53_500, 5): "CONCORDANT",
        (40_000, 0): "DISCORDANT",     # far below, disjoint
        (70_000, 0): "DISCORDANT",     # far above, disjoint
        (0, 0): "DISCORDANT",
        (52_300, 0): "INCONCLUSIVE",   # overlaps without containment
        (52_904, 0): "INCONCLUSIVE",
    }
    for (k, u), want in cases.items():
        assert camp.decide(k, u)["verdict"] == want, (k, u)


def test_decide_publishes_its_whole_arithmetic():
    d = camp.decide(53_500, 5)
    assert d["k_certain"] == 53_500 and d["u_ambiguous"] == 5
    assert 0 < d["garwood_low_counts"] < d["garwood_up_counts"]
    # the division happens at 96 bits BEFORE serialization, so the
    # published volumes agree with a double division only to 1 ulp
    assert math.isclose(d["c_lo_volume"],
                        d["garwood_low_counts"] / 440.0,
                        rel_tol=1e-15)
    assert math.isclose(d["c_hi_volume"],
                        d["garwood_up_counts"] / 440.0,
                        rel_tol=1e-15)
    assert math.isclose(d["band"], 0.025 * 121.834150532777,
                        rel_tol=1e-15)
    assert d["d_lo"] < d["d_hi"]
    assert d["oracle"]["v_lo"] == camp.FROZEN["v_lo"]
    assert "general in (K_certain, U_amb)" in d["rule"]
    with pytest.raises(ValueError):
        camp.decide(-1, 0)
    with pytest.raises(ValueError):
        camp.decide(0, -1)


def test_the_garwood_bounds_have_their_defining_coverage():
    """The frozen rule's own bounds satisfy their exact definitions,
    evaluated with the frozen 96-bit CDF at the bound."""

    for k in (1, 5, 52_986, 54_220):
        lam = camp.garwood_low(k, 0.05)
        # P(X >= k; lam) == alpha/2 at the lower bound
        assert abs(float(1 - camp.pois_cdf(k - 1, lam)) - 0.025) < 1e-9
        lam = camp.garwood_up(k, 0.05)
        assert abs(float(camp.pois_cdf(k, lam)) - 0.025) < 1e-9
    assert float(camp.garwood_low(0, 0.05)) == 0.0


def test_the_p2_engine_re_derives_the_acceptance_integers():
    """The frozen repo P2 engine (independent double-precision
    implementation) re-derives the SAME acceptance integers from the
    same definition -- the design audit's cross-engine agreement,
    kept alive as a contract."""

    a = camp.FROZEN["a_intensity"]
    band = camp.FROZEN["tau"] * camp.FROZEN["v_ref"]
    need_l = a * (camp.FROZEN["v_hi"] - band)
    cap_u = a * (camp.FROZEN["v_lo"] + band)
    k_lo, k_hi = camp.ACCEPTANCE
    lo_of = lambda k: p2.poisson_mean_ci(k, level=0.95)[0]  # noqa: E731
    up_of = lambda k: p2.poisson_mean_ci(k, level=0.95)[1]  # noqa: E731
    assert lo_of(k_lo) >= need_l > lo_of(k_lo - 1)
    assert up_of(k_hi) <= cap_u < up_of(k_hi + 1)


def test_the_power_sentence_numbers_reproduce():
    """Endpoint powers, the U = 30 window power, the U floor
    crossing, the Chernoff exponent and the joint bound -- from the
    frozen P2 CDF, pinned to the frozen 6-decimal quotes."""

    a = camp.FROZEN["a_intensity"]
    k_lo, k_hi = camp.ACCEPTANCE

    def window_power(lo, hi, v):
        lam = a * v
        return p2._poisson_cdf(hi, lam) - p2._poisson_cdf(lo - 1, lam)

    p_lo = window_power(k_lo, k_hi, camp.FROZEN["v_lo"])
    p_hi = window_power(k_lo, k_hi, camp.FROZEN["v_hi"])
    assert abs(p_lo - 0.937123) < 5e-7
    assert abs(p_hi - 0.931493) < 5e-7
    assert p_hi < p_lo                      # the worst endpoint
    p30 = window_power(k_lo + 30, k_hi - 30, camp.FROZEN["v_hi"])
    assert abs(p30 - 0.912741) < 5e-7
    p47 = window_power(k_lo + 47, k_hi - 47, camp.FROZEN["v_hi"])
    p48 = window_power(k_lo + 48, k_hi - 48, camp.FROZEN["v_hi"])
    assert p47 >= 0.90 > p48
    # exact Chernoff exponent on P(2N > 11M), N ~ Poisson(A*SCALE)
    lam_n = a * camp._RUNG["scale"]
    m = camp.FROZEN["scan_call_cap"] / 2.0
    u = (m - lam_n) / lam_n
    expo = lam_n * ((1 + u) * math.log1p(u) - u)
    assert 81_307 < expo < 81_308
    joint = p30 - 0.001 - 0.010 - math.exp(-700)   # exp(-expo) < 1e-300
    assert joint >= 0.901741 - 1e-6


# ------------------------------------------------ measure identities

def test_the_box_measure_identity_is_exact():
    g = camp._RUNG
    b = (2.0 * math.pi / 3.0 * (g["r_hi"] ** 3 - g["r_lo"] ** 3)
         * (1.0 - math.cos(g["psi_max"])))
    assert abs(b - g["scale"] / g["dt"]) <= 2e-13 * b


class _FixedU:
    def __init__(self, u, col):
        self._u, self._col = u, col

    def random(self, shape):
        out = np.zeros(shape)
        out[:, self._col] = self._u
        return out


def test_the_sampler_is_the_pilot_sampler_bit_for_bit():
    """The campaign sprinkles the SAME box measure the pilot scanned
    -- same inverse-CDF identities, and identical rows for identical
    uniforms."""

    import s6_m30_amb_pilot as pilot

    us = np.linspace(0.0, 1.0, 101)
    for col in range(4):
        ours = camp.draw_points(_FixedU(us, col), 101)
        theirs = pilot.draw_events(_FixedU(us, col), 101)
        assert (ours == theirs).all()
    ev = camp.draw_points(_FixedU(us, 0), 101)
    r3 = ev[:, 1] ** 3
    lo3 = camp._RUNG["r_lo"] ** 3
    hi3 = camp._RUNG["r_hi"] ** 3
    want = lo3 + us * (hi3 - lo3)
    assert np.allclose(r3, want, rtol=1e-12, atol=0.0)


def test_points_land_inside_the_box():
    rng = np.random.default_rng(40_000_311)   # the retired smoke stream
    ev = camp.draw_points(rng, 4096)
    t, r, psi, phi = ev[:, 0], ev[:, 1], ev[:, 2], ev[:, 3]
    assert (0 <= t).all() and (t <= camp._RUNG["dt"]).all()
    assert (camp._RUNG["r_lo"] <= r).all() and (r <= camp._RUNG["r_hi"]).all()
    assert (0 <= psi).all() and (psi <= camp._RUNG["psi_max"]).all()
    assert (0 <= phi).all() and (phi < 2 * math.pi).all()


def test_point_membership_is_tristate_with_early_exit(monkeypatch):
    """('in'|'out'|'ambiguous', calls) for every leg combination --
    and leg1 False never asks leg2, which is exactly why the
    campaign's ambiguity is A FORTIORI inside the pilot's."""

    def with_legs(*legs):
        seq = iter(legs)
        monkeypatch.setattr(camp.s1, "causal_relation",
                            lambda a, b, m, tol: next(seq))
        return camp.point_membership(1.0, 15.0, 0.1, 1e-8)

    assert with_legs(None) == ("ambiguous", 1)
    assert with_legs(False) == ("out", 1)
    assert with_legs(True, None) == ("ambiguous", 2)
    assert with_legs(True, False) == ("out", 2)
    assert with_legs(True, True) == ("in", 2)


# ------------------------------------------------ gates and refusals

@pytest.mark.parametrize("bad", ["a1f2063", "b" * 39, "g" * 40, ""])
def test_short_or_malformed_freeze_rev_is_refused(bad):
    with pytest.raises(SystemExit, match="full 40-hex"):
        camp.require_full_sha(bad)


def test_missing_freeze_rev_is_refused_by_the_cli():
    script = str(_REPO / "experiments" / "oracle"
                 / "s6_m30_count.py")
    out = subprocess.run([sys.executable, script, "--preflight"],
                         capture_output=True, text=True)
    assert out.returncode != 0
    assert "full 40-hex" in (out.stderr + out.stdout)


def _pass_static(monkeypatch, tmp_path, rev=_SHA):
    # the real stream is retired (results commit); the gate under
    # test is never the ledger here, so stub the allocation -- the
    # ledger's own refusal has its own test below
    monkeypatch.setattr(camp, "assert_fresh_scalar",
                        lambda name: 40_000_511)
    monkeypatch.setattr(camp, "verify_freeze", lambda stage: {})
    monkeypatch.setattr(camp, "_git_state",
                        lambda: {"rev": rev, "dirty": False,
                                 "dirt": []})
    monkeypatch.setattr(camp, "_ARTIFACT", tmp_path / "a.json")
    monkeypatch.setattr(camp, "_INCIDENT", tmp_path / "i.json")
    monkeypatch.setattr(camp, "_CHECKPOINT", tmp_path / "c.json")
    monkeypatch.setattr(camp.reservation, "verify_retained",
                        lambda: None)
    monkeypatch.setattr(camp.reservation, "held", lambda: None)
    monkeypatch.setattr(camp.reservation, "probe_namespace",
                        lambda: None)


def test_preflight_passes_only_on_the_exact_sha(monkeypatch, tmp_path):
    _pass_static(monkeypatch, tmp_path)
    assert camp.preflight(_SHA)["seed"] == 40_000_511
    with pytest.raises(SystemExit, match="exact\\s+approved commit"):
        camp.preflight("b" * 40)


def test_preflight_refuses_leftover_artifacts(monkeypatch, tmp_path):
    _pass_static(monkeypatch, tmp_path)
    (tmp_path / "i.json").write_text("{}", encoding="utf-8")
    with pytest.raises(SystemExit, match="write-once"):
        camp.preflight(_SHA)
    (tmp_path / "i.json").unlink()
    (tmp_path / "c.json").write_text("{}", encoding="utf-8")
    with pytest.raises(SystemExit, match="earlier attempt"):
        camp.preflight(_SHA)


def test_preflight_refuses_a_held_reservation(monkeypatch, tmp_path):
    _pass_static(monkeypatch, tmp_path)
    monkeypatch.setattr(camp.reservation, "held", lambda: "x" * 40)
    with pytest.raises(SystemExit, match="already held"):
        camp.preflight(_SHA)


# ------------------------------------------------ the reservation

def test_the_retained_set_names_all_seven_prior_claims():
    import s6_m30_count_reservation as res

    assert res.REF == "refs/s6m30/reservation"
    assert dict(res.RETAINED) == {
        "refs/o4/reservation":
            "c4da1626463e6a6505374813cf3f56d6b429c209",
        "refs/o4b/reservation":
            "46acee340bc247511546964b2925953721d5bb59",
        "refs/o5pilot/reservation":
            "f5223636bc7c60603b18ecae3f474dd3dc37146f",
        "refs/o5/reservation":
            "b5d9d6131644e7c8c07477c0b676696aaf2d1ce1",
        "refs/s6m14pilot/reservation":
            "3fc1b28bdf8014245078f2ac3fd08faa17ea1f5b",
        "refs/s6m14/reservation":
            "a6c8fb74cbe7910d05f228368bef20843f3d8ff9",
        "refs/s6m18pilot/reservation":
            "f695f96c03e5a0781428f4a2ddb77097ff16e94d",
        "refs/s6m18/reservation":
            "ba9c0e7fc21f9382d530cd429e32118d01c62233",
        "refs/s6m30pilot/reservation":
            "b1b10cd4a0a9fe064f8b2561d8f5dc5777a3242a",
    }


def test_a_cleared_or_rewritten_retained_ref_refuses(monkeypatch):
    import s6_m30_count_reservation as res

    monkeypatch.setattr(res.base, "_ls_remote", lambda ref: None)
    with pytest.raises(SystemExit, match="GONE"):
        res.verify_retained()
    monkeypatch.setattr(res.base, "_ls_remote", lambda ref: "e" * 40)
    with pytest.raises(SystemExit, match="rewritten"):
        res.verify_retained()


def test_make_commit_refuses_cleanly_when_git_fails(monkeypatch):
    """The PR #85 pre-push contract, built in from the start here."""

    import s6_m30_count_reservation as res

    class _Fail:
        returncode = 1
        stdout = ""
        stderr = "boom"

    monkeypatch.setattr(res.subprocess, "run", lambda *a, **k: _Fail())
    with pytest.raises(SystemExit, match="hash-object"):
        res._make_commit("x")

    def no_git(*a, **k):
        raise FileNotFoundError("git is not available")

    monkeypatch.setattr(res.subprocess, "run", no_git)
    with pytest.raises(SystemExit, match="could not run"):
        res._make_commit("x")


# ------------------------------------------------ checkpoints

def test_checkpoints_enforce_their_required_keys(monkeypatch,
                                                 tmp_path):
    monkeypatch.setattr(camp, "_CHECKPOINT", tmp_path / "ck.json")
    with pytest.raises(ValueError, match="missing"):
        camp.write_checkpoint({"freeze_sha": "x"})


def test_checkpoints_are_stamped_non_verdict_by_the_writer(
        monkeypatch, tmp_path):
    ck = tmp_path / "ck.json"
    monkeypatch.setattr(camp, "_CHECKPOINT", ck)
    payload = {kk: 0 for kk in camp._CK_REQUIRED}
    camp.write_checkpoint(payload)
    rec = json.loads(ck.read_text(encoding="utf-8"))
    assert rec["partial"] is True
    assert rec["non_verdict"] is True
    assert rec["stage"] == "chunk"


# ------------------------------------------------ the real main()

def _tiny(monkeypatch, tmp_path, ambiguous_every=0, out_every=0):
    """Drive the ACTUAL main() with a tiny frozen intensity and a stub
    predicate; nothing else is stubbed away from the real wiring.
    (The stub never calls flight_time, so the meter charges nothing
    -- the caps stay unbound, which these tests do not exercise.)"""

    small = dict(camp.FROZEN, a_intensity=64.0 / camp._RUNG["scale"],
                 chunk=16,
                 scan_call_cap=1_000, max_calls=1_000,
                 max_wall_s=600.0)
    monkeypatch.setattr(camp, "FROZEN", small)
    monkeypatch.setattr(camp, "G3A_PREFLIGHT_CALLS", 0)
    monkeypatch.setattr(camp.g3a, "run_preflight",
                        lambda tol, m=None, geometry=None:
                        {"passed": True})
    _pass_static(monkeypatch, tmp_path)
    monkeypatch.setattr(camp.reservation, "claim",
                        lambda payload: "d" * 40)
    monkeypatch.setattr(camp.reservation, "verify_still_held",
                        lambda obj: obj)
    calls = {"n": 0}

    def stub_relation(a, b, m, tol):
        calls["n"] += 1
        if ambiguous_every and calls["n"] % ambiguous_every == 0:
            return None
        if out_every and calls["n"] % out_every == 0:
            return False
        return True

    monkeypatch.setattr(camp.s1, "causal_relation", stub_relation)
    return small, calls


def test_main_publishes_with_the_general_rule(monkeypatch, tmp_path):
    small, _ = _tiny(monkeypatch, tmp_path, ambiguous_every=13)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    camp.main()
    art = json.loads((tmp_path / "a.json").read_text(encoding="utf-8"))
    scan = art["scan"]
    assert scan["points_done"] == scan["n_total"] > 0
    assert scan["u_ambiguous"] > 0          # the stub produced some
    # the published decision is EXACTLY the general rule applied to
    # the observed tallies under the run's own frozen config
    assert art["decision"] == camp.decide(
        scan["k_certain"], scan["u_ambiguous"], small)
    assert art["reservation"]["object"] == "d" * 40
    assert art["reservation"]["verified_at_exit"] is True
    ck = json.loads((tmp_path / "c.json").read_text(encoding="utf-8"))
    assert ck["non_verdict"] is True
    assert ck["n_total"] == scan["n_total"]
    assert ck["points_done"] == scan["points_done"]


def test_the_sprinkled_n_is_poisson_from_the_stream(monkeypatch,
                                                    tmp_path):
    """N is drawn ONCE from the same seeded stream, before any point
    -- the artifact's n_total must equal an independent replay of the
    frozen draw."""

    small, _ = _tiny(monkeypatch, tmp_path)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    camp.main()
    art = json.loads((tmp_path / "a.json").read_text(encoding="utf-8"))
    rng = np.random.default_rng(40_000_511)
    want = int(rng.poisson(small["a_intensity"] * camp._RUNG["scale"]))
    assert art["scan"]["n_total"] == want


def test_a_mid_chunk_failure_preserves_the_last_observation_point(
        monkeypatch, tmp_path):
    """The incident carries the ACTUAL last completed observation,
    not the last successful chunk (review PR #87 R1): the stub dies
    on point 21's first leg, so with chunk = 16 the incident must say
    20 points -- four PAST the chunk boundary -- and n_total must be
    the drawn N, never null."""

    small, calls = _tiny(monkeypatch, tmp_path)

    def dying_relation(a, b, m, tol):
        calls["n"] += 1
        if calls["n"] > 40:                   # 2 calls per point
            raise RuntimeError("solver died mid-scan")
        return True

    monkeypatch.setattr(camp.s1, "causal_relation", dying_relation)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(SystemExit, match="incident written"):
        camp.main()
    inc = json.loads((tmp_path / "i.json").read_text(encoding="utf-8"))
    assert inc["failure_point"] == "scan"
    assert inc["verdict"] is None
    assert inc["seed_spent"] is True
    pre = inc["preserved"]
    assert pre["points_done"] == 20 > small["chunk"]
    assert pre["k_certain"] == 20 and pre["calls"] == 40
    rng = np.random.default_rng(40_000_511)
    assert pre["n_total"] == int(
        rng.poisson(small["a_intensity"] * camp._RUNG["scale"]))
    # the chunk-granular checkpoint holds the last COMPLETE chunk
    ck = json.loads((tmp_path / "c.json").read_text(encoding="utf-8"))
    assert ck["points_done"] == small["chunk"]
    assert not (tmp_path / "a.json").exists()


def test_a_first_chunk_failure_still_preserves_the_drawn_n(
        monkeypatch, tmp_path):
    """Even before ANY chunk completes, the incident must carry the
    drawn N and the in-flight tallies -- the first-chunk corner the
    review named."""

    small, calls = _tiny(monkeypatch, tmp_path)

    def dying_relation(a, b, m, tol):
        calls["n"] += 1
        if calls["n"] > 6:                    # dies inside chunk 1
            raise RuntimeError("solver died early")
        return True

    monkeypatch.setattr(camp.s1, "causal_relation", dying_relation)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(SystemExit, match="incident written"):
        camp.main()
    inc = json.loads((tmp_path / "i.json").read_text(encoding="utf-8"))
    pre = inc["preserved"]
    assert pre["n_total"] is not None and pre["n_total"] > 0
    assert pre["points_done"] == 3
    assert not (tmp_path / "c.json").exists()   # no chunk completed


def test_a_race_losing_clean_refusal_files_no_incident(monkeypatch,
                                                       tmp_path):
    _tiny(monkeypatch, tmp_path)

    def losing_claim(payload):
        raise SystemExit("reservation: refs/s6m30/reservation is "
                         "already held by someone")

    monkeypatch.setattr(camp.reservation, "claim", losing_claim)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(SystemExit, match="already held"):
        camp.main()
    assert not (tmp_path / "i.json").exists()
    assert not (tmp_path / "a.json").exists()


def test_an_uncertain_claim_files_a_fail_closed_incident(monkeypatch,
                                                         tmp_path):
    _tiny(monkeypatch, tmp_path)

    def uncertain_claim(payload):
        raise camp.reservation.ClaimUncertain(
            "the push did not return", pushed=None, obj="e" * 40,
            detail="socket")

    monkeypatch.setattr(camp.reservation, "claim", uncertain_claim)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(SystemExit, match="incident written"):
        camp.main()
    inc = json.loads((tmp_path / "i.json").read_text(encoding="utf-8"))
    assert inc["failure_point"] == "reservation_claim"
    assert inc["reservation_claimed"] == "uncertain"
    assert inc["seed_spent"] is True
    rec = inc["reservation_uncertainty"]
    assert rec["seeds_must_be_treated_as"] == "spent"


def test_publish_time_reverify_uncertainty_names_publish(
        monkeypatch, tmp_path):
    """The pilot's PR #85 contract, built in from the start: the
    stage names the failure point, the record lands under its own
    key, and the confirmed claim stays confirmed."""

    small, _ = _tiny(monkeypatch, tmp_path)

    def uncertain(obj):
        raise camp.reservation.ClaimUncertain(
            "the ref could not be re-read before publication",
            pushed=None, obj=obj, detail="network failed")

    monkeypatch.setattr(camp.reservation, "verify_still_held",
                        uncertain)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(SystemExit, match="incident written"):
        camp.main()
    assert not (tmp_path / "a.json").exists()
    inc = json.loads((tmp_path / "i.json").read_text(encoding="utf-8"))
    assert inc["failure_point"] == "publish"
    assert inc["reservation_claimed"] is True
    assert inc["reservation_object"] == "d" * 40
    assert inc["reservation_uncertainty"] is None
    rec = inc["publish_reverify_uncertainty"]
    assert rec["seeds_must_be_treated_as"] == "spent"
    assert inc["preserved"]["points_done"] == inc["preserved"]["n_total"]


def test_an_exit_rebind_refuses_publication(monkeypatch, tmp_path):
    """HEAD moved mid-run: the exit-lineage gate refuses and the
    incident names the publish stage."""

    _tiny(monkeypatch, tmp_path)
    states = iter([{"rev": _SHA, "dirty": False, "dirt": []},
                   {"rev": "f" * 40, "dirty": False, "dirt": []}])
    monkeypatch.setattr(camp, "_git_state", lambda: next(states))
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(SystemExit, match="incident written"):
        camp.main()
    inc = json.loads((tmp_path / "i.json").read_text(encoding="utf-8"))
    assert inc["failure_point"] == "publish"
    assert "exit lineage" in inc["termination_reason"]
    assert not (tmp_path / "a.json").exists()


def test_a_broken_pipe_after_publication_files_no_incident(
        monkeypatch, tmp_path):
    _tiny(monkeypatch, tmp_path)
    real_print = print

    def broken_print(*args, **kwargs):
        if args and str(args[0]).startswith("result:"):
            raise BrokenPipeError("stdout gone")
        real_print(*args, **kwargs)

    monkeypatch.setattr("builtins.print", broken_print)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(BrokenPipeError):
        camp.main()
    assert (tmp_path / "a.json").exists()       # published
    assert not (tmp_path / "i.json").exists()   # and no incident


def test_a_write_once_refusal_on_our_own_artifact_files_no_incident(
        monkeypatch, tmp_path):
    _tiny(monkeypatch, tmp_path)
    (tmp_path / "a.json").write_text(json.dumps(
        {"reservation": {"object": "d" * 40}}), encoding="utf-8")
    monkeypatch.setattr(
        camp, "preflight",
        lambda rev: {"manifest": {}, "seed": 40_000_511,
                     "git": {"rev": _SHA, "dirty": False,
                             "dirt": []}})
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(SystemExit, match="was published"):
        camp.main()
    assert not (tmp_path / "i.json").exists()


def test_a_write_once_refusal_on_a_foreign_artifact_files_an_incident(
        monkeypatch, tmp_path):
    _tiny(monkeypatch, tmp_path)
    (tmp_path / "a.json").write_text(json.dumps(
        {"reservation": {"object": "f" * 40}}), encoding="utf-8")
    monkeypatch.setattr(
        camp, "preflight",
        lambda rev: {"manifest": {}, "seed": 40_000_511,
                     "git": {"rev": _SHA, "dirty": False,
                             "dirt": []}})
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(SystemExit, match="incident written"):
        camp.main()
    inc = json.loads((tmp_path / "i.json").read_text(encoding="utf-8"))
    assert inc["failure_point"] == "publish"
    assert inc["seed_spent"] is True


# ------------------------------------------------ manifest and doc

def test_the_committed_manifest_is_the_tree_it_describes():
    committed = json.loads(camp._MANIFEST.read_text(encoding="utf-8"))
    assert committed == camp.build_manifest()
    assert set(committed["files"]) == set(camp.PROTOCOL_SURFACE)
    # the artifacts the frozen literals quote are pinned
    for rel in ("docs/prereg/p14_s6_m30_volume.json",
                "docs/prereg/p14_o3_volume.json",
                "docs/prereg/p14_s6_m30_pilot.json"):
        assert rel in committed["files"]


def test_the_doc_freezes_the_discipline():
    doc = (_REPO / "docs" / "prereg"
           / "p14_s6_m30_count.md").read_text(encoding="utf-8")
    flat = " ".join(doc.split())
    assert "no auto-raise, ever" in flat
    assert "stops this rung for a PI ruling" in flat
    assert "general in (K_certain, U_amb)" in flat.lower() or \
        "general in (k_certain, u_amb)" in flat.lower()
    assert "Wall-cap risk is NOT included" in flat
    assert "never twice" in flat
    assert "a fortiori narrower" in flat
    assert "V_op − V_true" in flat or "V_op - V_true" in flat


def test_the_rung_literals_match_the_one_exact_path():
    import s6_rungs as s6

    assert camp._RUNG == {k: s6.rung_geometry(3.0)[k]
                          for k in camp._RUNG}
    assert camp._P[1] == s6.R_IN == 12.0
    assert camp._Q[1] == s6.R_OUT == 18.0
    assert camp._Q[0] == camp._RUNG["dt"]


def test_the_published_rule_string_is_unprimed():
    d = camp.decide(53_500, 5)
    assert "V'" not in d["rule"] and "V_ref'" not in d["rule"]
    assert "O5" not in d["rule"]
