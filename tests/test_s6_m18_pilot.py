"""S6 M=1.8 ambiguity-pilot freeze contract tests.

The frozen ruled values, the general-k decision rule cross-checked by
two independent exact engines plus the frozen P2 engine, the sampler's
measure identities, the approved-SHA gates, the preflight refusals,
the reservation ordering, the chunked atomic checkpoints, and the
REAL `main()` wiring end to end -- success, mid-run rebind refusal,
and the incident path (the O4b lesson: no returning stubs standing in
for the actual wiring)."""

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

import s6_m18_amb_pilot as pilot  # noqa: E402

_SHA = "a" * 40


# ------------------------------------------------ the frozen values

def test_frozen_configuration_is_the_ruled_one():
    assert pilot.FROZEN == {
        "n_events": 3_940_846,
        "u_max": 30,
        "alpha_pilot": 0.01,
        "tail_budget": 0.001,
        "a_provisional": 720.0,
        "seed_name": "s6_m18_pilot",
        "tol": 1e-08,
        "chunk": 65_536,
        "max_calls": 7_888_153,
        "max_wall_s": 86_400.0,
    }
    assert pilot.G3A_PREFLIGHT_CALLS == 6_461
    assert pilot.FROZEN["max_calls"] == (
        2 * pilot.FROZEN["n_events"] + pilot.G3A_PREFLIGHT_CALLS)


def test_the_seed_is_fresh_and_301_untouched():
    import probe_seed_ledger as ledger

    assert "s6_m18_pilot" not in ledger.FRESH_PROBE_SCALARS
    assert ledger.OBSERVED_PROBE_SCALARS["s6_m18_pilot"] == 40_000_481
    assert ledger.replay_scalar("s6_m18_pilot") == 40_000_481
    assert 40_000_301 not in ledger.spent_scalars()
    assert 40_000_301 not in ledger.FRESH_PROBE_SCALARS.values()


# ------------------------------------------------ measure identities

def test_the_box_measure_identity_is_exact():
    """B_spatial = 2pi/3 (R_HI^3 - R_LO^3)(1 - cos PSI_MAX) must equal
    SCALE/DT to the last bit of the frozen constants -- the sampler
    and the campaign estimand share one measure."""

    g = pilot._RUNG
    b = (2.0 * math.pi / 3.0 * (g["r_hi"] ** 3 - g["r_lo"] ** 3)
         * (1.0 - math.cos(g["psi_max"])))
    assert abs(b - g["scale"] / g["dt"]) <= 2e-13 * b
    import s6_rungs as s6
    assert g == {k: s6.rung_geometry(1.8)[k] for k in g}


class _FixedU:
    """A 'generator' returning fixed uniforms in one chosen column."""

    def __init__(self, u, col):
        self._u, self._col = u, col

    def random(self, shape):
        out = np.zeros(shape)
        out[:, self._col] = self._u
        return out


def test_the_r_sampler_inverts_the_r_squared_cdf():
    """r(u)^3 is EXACTLY linear in u between R_LO^3 and R_HI^3 -- the
    inverse-CDF identity of the r^2 dr density, checked on a grid."""

    us = np.linspace(0.0, 1.0, 101)
    ev = pilot.draw_events(_FixedU(us, 0), 101)
    r3 = ev[:, 1] ** 3
    lo3 = pilot._RUNG["r_lo"] ** 3
    hi3 = pilot._RUNG["r_hi"] ** 3
    want = lo3 + us * (hi3 - lo3)
    assert np.allclose(r3, want, rtol=1e-12, atol=0.0)
    # and cos(psi) is exactly linear between cos(PSI_MAX) and 1
    ev2 = pilot.draw_events(_FixedU(us, 1), 101)
    cpsi = np.cos(ev2[:, 2])
    cmin = math.cos(pilot._RUNG["psi_max"])
    want2 = cmin + us * (1.0 - cmin)
    assert np.allclose(cpsi, want2, rtol=1e-12, atol=1e-15)


def test_events_land_inside_the_box():
    rng = np.random.default_rng(40_000_311)   # the retired smoke stream
    ev = pilot.draw_events(rng, 4096)
    t, r, psi, phi = ev[:, 0], ev[:, 1], ev[:, 2], ev[:, 3]
    assert (0 <= t).all() and (t <= pilot._RUNG["dt"]).all()
    assert (pilot._RUNG["r_lo"] <= r).all() and (r <= pilot._RUNG["r_hi"]).all()
    assert (0 <= psi).all() and (psi <= pilot._RUNG["psi_max"]).all()
    assert (0 <= phi).all() and (phi < 2 * math.pi).all()


# ------------------------------------------------ the frozen rule

def test_the_rule_is_general_in_k_and_matches_the_reference_table():
    """The frozen reference evaluations, from the doc -- but the RULE
    computes them; nothing is hardcoded to k <= 1."""

    expect = {
        0: (11.5092, 1.474929e-06, "FEASIBLE"),
        1: (16.5905, 9.999959e-04, "FEASIBLE"),
        2: (21.0081, 2.426052e-02, "INCONCLUSIVE"),
        3: (25.1047, 1.414954e-01, "INCONCLUSIVE"),
    }
    for k, (lam, tail, verdict) in expect.items():
        d = pilot.decide(k)
        assert d["verdict"] == verdict
        assert abs(d["lambda_u"] - lam) < 5e-4
        # 6-significant-digit reference values -> relative tolerance
        assert abs(d["tail_p_u_gt_u_max"] - tail) < 5e-7 * tail + 1e-12
    # far k keeps working and keeps failing feasibility
    far = pilot.decide(1_000)
    assert far["verdict"] == "INCONCLUSIVE"
    assert far["lambda_u"] > 1_000


def test_the_verdict_path_is_cross_checked_by_two_other_engines():
    """The RUNTIME verdict path is 96-bit gmpy2 (review PR #83 R1:
    the k=1 boundary sits 4.11e-9 of tail below the budget, and a
    platform-libm double engine could flip it). The cross-checks are
    an independent double-precision log-space engine -- agreement
    above its ~1e-8 lgamma-cancellation floor -- and the frozen P14
    P2 gamma engine for the Poisson tail."""

    import math as m

    import p14_probe_p2 as p2

    n, alpha = pilot.FROZEN["n_events"], pilot.FROZEN["alpha_pilot"]

    def cp_double(k):
        if k == 0:
            return 1.0 - alpha ** (1.0 / n)
        lo, hi = k / n, 1.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            s, lp, lq = 0.0, m.log(mid), m.log1p(-mid)
            for i in range(k + 1):
                s += m.exp(m.lgamma(n + 1) - m.lgamma(i + 1)
                           - m.lgamma(n - i + 1) + i * lp
                           + (n - i) * lq)
            if s > alpha:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    for k in (0, 1, 2, 5):
        pu = pilot.cp_upper(k, n, alpha)          # the 96-bit path
        pu_d = cp_double(k)
        assert abs(pu - pu_d) < 1e-7 * pu + 1e-300
        lam = pilot.FROZEN["a_provisional"] * pilot._RUNG["scale"] * pu
        t1 = pilot.pois_tail_gt(pilot.FROZEN["u_max"], lam)
        t2 = 1.0 - p2._poisson_cdf(pilot.FROZEN["u_max"], lam)
        assert abs(t1 - t2) < 1e-11
        # verdict-level agreement between the engines is exact
        lam_d = pilot.FROZEN["a_provisional"] * pilot._RUNG["scale"] * pu_d
        t_d = pilot.pois_tail_gt(pilot.FROZEN["u_max"], lam_d)
        assert ((t1 <= pilot.FROZEN["tail_budget"])
                == (t_d <= pilot.FROZEN["tail_budget"]))


def test_the_k1_boundary_margin_is_the_mpfr_determined_value():
    """The margin the review computed, reproduced by the runtime path:
    k=1 sits 4.1062e-09 of tail BELOW the budget. Pinned tightly --
    the 96-bit path is host-independent (gmpy2/MPFR are in the
    environment lock; platform libm is not on this path)."""

    margin = (pilot.FROZEN["tail_budget"]
              - pilot.decide(1)["tail_p_u_gt_u_max"])
    assert abs(margin - 4.1061750e-09) < 1e-15
    assert pilot.decide(1)["verdict"] == "FEASIBLE"
    assert pilot.decide(2)["verdict"] == "INCONCLUSIVE"


# ------------------------------------------------ gates and refusals

@pytest.mark.parametrize("bad", ["30d286a", "a" * 39, "g" * 40, ""])
def test_short_or_malformed_freeze_rev_is_refused(bad):
    with pytest.raises(SystemExit, match="full 40-hex"):
        pilot.require_full_sha(bad)


def test_missing_freeze_rev_is_refused_by_the_cli():
    script = str(_REPO / "experiments" / "oracle" / "s6_m18_amb_pilot.py")
    out = subprocess.run([sys.executable, script, "--preflight"],
                         capture_output=True, text=True)
    assert out.returncode != 0
    assert "full 40-hex" in (out.stderr + out.stdout)


def _pass_static(monkeypatch, tmp_path, rev=_SHA):
    # the real stream is retired (results commit); the gate under
    # test is never the ledger here
    monkeypatch.setattr(pilot, "assert_fresh_scalar",
                        lambda name: 40_000_481)
    monkeypatch.setattr(pilot, "verify_freeze", lambda stage: {})
    monkeypatch.setattr(pilot, "_git_state",
                        lambda: {"rev": rev, "dirty": False,
                                 "dirt": []})
    monkeypatch.setattr(pilot, "_ARTIFACT", tmp_path / "a.json")
    monkeypatch.setattr(pilot, "_INCIDENT", tmp_path / "i.json")
    monkeypatch.setattr(pilot, "_CHECKPOINT", tmp_path / "c.json")
    monkeypatch.setattr(pilot.reservation, "verify_retained",
                        lambda: None)
    monkeypatch.setattr(pilot.reservation, "held", lambda: None)
    monkeypatch.setattr(pilot.reservation, "probe_namespace",
                        lambda: None)


def test_preflight_passes_only_on_the_exact_sha(monkeypatch, tmp_path):
    _pass_static(monkeypatch, tmp_path)
    assert pilot.preflight(_SHA)["seed"] == 40_000_481
    with pytest.raises(SystemExit, match="exact\\s+approved commit"):
        pilot.preflight("b" * 40)


def test_preflight_refuses_leftover_artifacts(monkeypatch, tmp_path):
    _pass_static(monkeypatch, tmp_path)
    (tmp_path / "i.json").write_text("{}", encoding="utf-8")
    with pytest.raises(SystemExit, match="write-once"):
        pilot.preflight(_SHA)
    (tmp_path / "i.json").unlink()
    (tmp_path / "c.json").write_text("{}", encoding="utf-8")
    with pytest.raises(SystemExit, match="earlier attempt"):
        pilot.preflight(_SHA)


def test_preflight_refuses_a_held_reservation(monkeypatch, tmp_path):
    _pass_static(monkeypatch, tmp_path)
    monkeypatch.setattr(pilot.reservation, "held", lambda: "x" * 40)
    with pytest.raises(SystemExit, match="already held"):
        pilot.preflight(_SHA)


# ------------------------------------------------ checkpoints

def test_checkpoints_enforce_their_required_keys(monkeypatch,
                                                 tmp_path):
    monkeypatch.setattr(pilot, "_CHECKPOINT", tmp_path / "ck.json")
    with pytest.raises(ValueError, match="missing"):
        pilot.write_checkpoint({"freeze_sha": "x"})


def test_checkpoints_are_stamped_non_verdict_by_the_writer(
        monkeypatch, tmp_path):
    ck = tmp_path / "ck.json"
    monkeypatch.setattr(pilot, "_CHECKPOINT", ck)
    payload = {kk: 0 for kk in pilot._CK_REQUIRED}
    pilot.write_checkpoint(payload)
    rec = json.loads(ck.read_text(encoding="utf-8"))
    assert rec["partial"] is True
    assert rec["non_verdict"] is True
    assert rec["stage"] == "chunk"


# ------------------------------------------------ the real main()

def _tiny(monkeypatch, tmp_path, ambiguous_every=0):
    """Drive the ACTUAL main() with a tiny frozen config and a stub
    predicate; nothing else is stubbed away from the real wiring.
    (The stub never calls flight_time, so the meter charges nothing
    -- the caps stay unbound, which these tests do not exercise.)"""

    small = dict(pilot.FROZEN, n_events=64, chunk=16,
                 max_calls=2 * 64, max_wall_s=600.0)
    monkeypatch.setattr(pilot, "FROZEN", small)
    monkeypatch.setattr(pilot, "G3A_PREFLIGHT_CALLS", 0)
    monkeypatch.setattr(pilot.g3a, "run_preflight",
                        lambda tol, m=None, geometry=None: {"passed": True})
    _pass_static(monkeypatch, tmp_path)
    monkeypatch.setattr(pilot.reservation, "claim",
                        lambda payload: "d" * 40)
    monkeypatch.setattr(pilot.reservation, "verify_still_held",
                        lambda obj: obj)
    calls = {"n": 0}

    def stub_relation(a, b, m, tol):
        calls["n"] += 1
        if ambiguous_every and calls["n"] % ambiguous_every == 0:
            return None
        return True

    monkeypatch.setattr(pilot.s1, "causal_relation", stub_relation)
    return small, calls


def test_main_publishes_with_the_general_rule(monkeypatch, tmp_path):
    small, _ = _tiny(monkeypatch, tmp_path, ambiguous_every=13)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    pilot.main()
    art = json.loads((tmp_path / "a.json").read_text(encoding="utf-8"))
    k = art["scan"]["k_ambiguous"]
    assert k > 0                     # the stub produced ambiguity
    # the published decision is EXACTLY the general rule applied to
    # the observed k under the run's own frozen config
    want = pilot.decide(k, small)
    assert art["decision"]["p_upper"] == want["p_upper"]
    assert art["decision"]["tail_p_u_gt_u_max"] == \
        want["tail_p_u_gt_u_max"]
    assert art["decision"]["verdict"] == want["verdict"]
    assert art["reservation"]["object"] == "d" * 40
    assert art["code"]["start"] == art["code"]["end"]
    assert art["scan"]["events_done"] == 64
    assert not (tmp_path / "i.json").exists()
    # every chunk left an atomic non-verdict checkpoint behind
    ck = json.loads((tmp_path / "c.json").read_text(encoding="utf-8"))
    assert ck["events_done"] == 64
    assert ck["rng_stream"] == "s6_m18_pilot"
    assert ck["non_verdict"] is True


def test_a_mid_run_rebind_refuses_before_publishing(monkeypatch,
                                                    tmp_path):
    _tiny(monkeypatch, tmp_path)
    # preflight sees the approved SHA; the exit re-read sees another
    seen = {"n": 0}

    def moving_state():
        seen["n"] += 1
        if seen["n"] >= 2:                      # after preflight
            return {"rev": "b" * 40, "dirty": False, "dirt": []}
        return {"rev": _SHA, "dirty": False, "dirt": []}

    monkeypatch.setattr(pilot, "_git_state", moving_state)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(SystemExit, match="incident written"):
        pilot.main()
    inc = json.loads((tmp_path / "i.json").read_text(encoding="utf-8"))
    assert inc["failure_point"] == "publish"
    assert "exit lineage" in inc["termination_reason"]
    assert not (tmp_path / "a.json").exists()


def test_a_scan_failure_files_an_incident_with_preserved_tallies(
        monkeypatch, tmp_path):
    _tiny(monkeypatch, tmp_path)
    boom = {"n": 0}

    def failing(a, b, m, tol):
        boom["n"] += 1
        if boom["n"] > 40:              # chunk 1 (32 calls) completes
            raise RuntimeError("solver blew up")
        return True

    monkeypatch.setattr(pilot.s1, "causal_relation", failing)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(SystemExit, match="incident written"):
        pilot.main()
    inc = json.loads((tmp_path / "i.json").read_text(encoding="utf-8"))
    assert inc["failure_point"] == "scan"
    assert inc["verdict"] is None
    assert inc["seed_spent"] is True            # the claim succeeded
    assert inc["preserved"]["events_done"] == 16  # first chunk landed


def test_a_race_losing_clean_claim_refusal_files_no_incident(
        monkeypatch, tmp_path):
    """Review PR #83 R2: two checkouts pass preflight on an empty ref;
    A claims first; B's claim() sees A's object at held() and refuses
    with a PRE-PUSH SystemExit. B wrote nothing and spent nothing, so
    B may not file the pilot's global incident -- it would sit beside
    A's legitimate result and poison the write-once provenance."""

    _tiny(monkeypatch, tmp_path)

    def lost_race(payload):
        raise SystemExit(
            "reservation: refs/o5pilot/reservation is already held "
            "by " + "e" * 40)

    monkeypatch.setattr(pilot.reservation, "claim", lost_race)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(SystemExit, match="already held"):
        pilot.main()
    assert not (tmp_path / "i.json").exists()
    assert not (tmp_path / "a.json").exists()


def test_an_uncertain_push_still_files_a_fail_closed_incident(
        monkeypatch, tmp_path):
    """The other side of the R2 boundary: once the push was ATTEMPTED
    the outcome is uncertain and the seed possibly spent, so the
    incident must exist and must say so."""

    _tiny(monkeypatch, tmp_path)

    def uncertain(payload):
        raise pilot.reservation.ClaimUncertain(
            "the push did not return", pushed=None, obj="c" * 40,
            detail="interrupted")

    monkeypatch.setattr(pilot.reservation, "claim", uncertain)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(SystemExit, match="incident written"):
        pilot.main()
    inc = json.loads((tmp_path / "i.json").read_text(encoding="utf-8"))
    assert inc["failure_point"] == "reservation_claim"
    assert inc["reservation_claimed"] == "uncertain"
    assert inc["seed_spent"] is True
    assert inc["reservation_uncertainty"][
        "seeds_must_be_treated_as"] == "spent"


# -------------------------------- the publication commit boundary

def test_a_failure_after_the_link_files_no_incident(monkeypatch,
                                                    tmp_path):
    """Review PR #83 R1: an interrupt one line after os.link -- here a
    BrokenPipeError out of the success print -- must NOT file an
    incident beside the published result."""

    _tiny(monkeypatch, tmp_path)

    def broken_print(*args, **kwargs):
        raise BrokenPipeError("stdout gone")

    monkeypatch.setattr("builtins.print", broken_print)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(BrokenPipeError):
        pilot.main()
    assert (tmp_path / "a.json").exists()       # published
    assert not (tmp_path / "i.json").exists()   # and no incident


def test_a_write_once_refusal_on_our_own_artifact_files_no_incident(
        monkeypatch, tmp_path):
    """If the destination already carries THIS run's claim object, the
    refusal is a post-commit fact, not a failure."""

    _tiny(monkeypatch, tmp_path)
    (tmp_path / "a.json").write_text(json.dumps(
        {"reservation": {"object": "d" * 40}}), encoding="utf-8")
    # preflight must not see it, or it refuses before the run
    monkeypatch.setattr(
        pilot, "preflight",
        lambda rev: {"manifest": {}, "seed": 40_000_481,
                     "git": {"rev": _SHA, "dirty": False,
                             "dirt": []}})
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(SystemExit, match="was published"):
        pilot.main()
    assert not (tmp_path / "i.json").exists()


def test_a_write_once_refusal_on_a_foreign_artifact_files_an_incident(
        monkeypatch, tmp_path):
    """A fully scanned run that finds a FOREIGN artifact on its path
    still leaves a record: the seed is spent and neither result nor
    silence may stand in for the incident."""

    _tiny(monkeypatch, tmp_path)
    (tmp_path / "a.json").write_text(json.dumps(
        {"reservation": {"object": "f" * 40}}), encoding="utf-8")
    monkeypatch.setattr(
        pilot, "preflight",
        lambda rev: {"manifest": {}, "seed": 40_000_481,
                     "git": {"rev": _SHA, "dirty": False,
                             "dirt": []}})
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(SystemExit, match="incident written"):
        pilot.main()
    inc = json.loads((tmp_path / "i.json").read_text(encoding="utf-8"))
    assert inc["failure_point"] == "publish"
    assert inc["seed_spent"] is True


# ------------------------------------------------ manifest and doc

def test_the_committed_manifest_is_the_tree_it_describes():
    committed = json.loads(pilot._MANIFEST.read_text(encoding="utf-8"))
    assert committed == pilot.build_manifest()
    assert set(committed["files"]) == set(pilot.PROTOCOL_SURFACE)


def test_the_doc_freezes_the_discipline():
    doc = (_REPO / "docs" / "prereg"
           / "p14_s6_m18_pilot.md").read_text(encoding="utf-8")
    flat = " ".join(doc.split())
    assert "no auto-raise, ever" in flat
    assert "general in k" in flat.lower()
    assert "sized at 96-bit from the ACTUAL rung endpoints" in flat
    assert "Wall-cap risk is NOT included" in flat
    assert "never twice" in flat
    assert "indicts nothing" in flat


def test_publish_time_reverify_uncertainty_names_publish(
        monkeypatch, tmp_path):
    """The O5 PR #85 contract, inherited by the mirror: the stage
    names the failure point, the record lands under its own key, and
    the confirmed claim stays confirmed."""

    small, _ = _tiny(monkeypatch, tmp_path)

    def uncertain(obj):
        raise pilot.reservation.ClaimUncertain(
            "the ref could not be re-read before publication",
            pushed=None, obj=obj, detail="network failed")

    monkeypatch.setattr(pilot.reservation, "verify_still_held",
                        uncertain)
    monkeypatch.setattr(sys, "argv", ["prog", "--freeze-rev", _SHA])
    with pytest.raises(SystemExit, match="incident written"):
        pilot.main()
    assert not (tmp_path / "a.json").exists()
    inc = json.loads((tmp_path / "i.json").read_text(encoding="utf-8"))
    assert inc["failure_point"] == "publish"
    assert inc["reservation_claimed"] is True
    assert inc["reservation_object"] == "d" * 40
    assert inc["reservation_uncertainty"] is None
    assert inc["seed_spent"] is True
    rec = inc["publish_reverify_uncertainty"]
    assert rec["seeds_must_be_treated_as"] == "spent"
    assert rec["push_reported"] == "did not report"
    assert inc["preserved"]["events_done"] == small["n_events"]


def test_the_anchors_come_from_the_ladder_source_of_truth():
    """_P/_Q reference s6_rungs' absolute anchors -- a ladder-anchor
    change cannot silently desync this runner (review PR #95 R1)."""

    import s6_rungs as s6

    assert pilot._P[1] == s6.R_IN == 12.0
    assert pilot._Q[1] == s6.R_OUT == 18.0
    assert pilot._Q[0] == pilot._RUNG["dt"]
    assert pilot._P[0] == 0.0


def test_make_commit_refuses_cleanly_when_git_fails(monkeypatch):
    """The PR #85 pre-push contract holds in THIS reservation copy."""

    import s6_m18_reservation as res

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
