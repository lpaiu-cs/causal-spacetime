"""O5 ambiguity-pilot freeze contract tests.

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

import o4_sizing as sz  # noqa: E402
import o5_amb_pilot as pilot  # noqa: E402

_SHA = "a" * 40


# ------------------------------------------------ the frozen values

def test_frozen_configuration_is_the_ruled_one():
    assert pilot.FROZEN == {
        "n_events": 10_736_965,
        "u_max": 30,
        "alpha_pilot": 0.01,
        "tail_budget": 0.001,
        "a_provisional": 940.0,
        "seed_name": "o5_amb_pilot",
        "tol": 1e-08,
        "chunk": 65_536,
        "max_calls": 21_480_467,
        "max_wall_s": 86_400.0,
    }
    assert pilot.G3A_PREFLIGHT_CALLS == 6_537
    assert pilot.FROZEN["max_calls"] == (
        2 * pilot.FROZEN["n_events"] + pilot.G3A_PREFLIGHT_CALLS)


def test_the_seed_is_fresh_and_301_untouched():
    import probe_seed_ledger as ledger

    assert ledger.assert_fresh_scalar("o5_amb_pilot") == 40_000_441
    assert 40_000_301 not in ledger.spent_scalars()
    assert 40_000_301 not in ledger.FRESH_PROBE_SCALARS.values()


# ------------------------------------------------ measure identities

def test_the_box_measure_identity_is_exact():
    """B_spatial = 2pi/3 (R_HI^3 - R_LO^3)(1 - cos PSI_MAX) must equal
    SCALE/DT to the last bit of the frozen constants -- the sampler
    and the campaign estimand share one measure."""

    b = (2.0 * math.pi / 3.0 * (sz.R_HI ** 3 - sz.R_LO ** 3)
         * (1.0 - math.cos(sz.PSI_MAX)))
    assert b == sz.SCALE / sz.DT


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
    want = sz.R_LO ** 3 + us * (sz.R_HI ** 3 - sz.R_LO ** 3)
    assert np.allclose(r3, want, rtol=1e-12, atol=0.0)
    # and cos(psi) is exactly linear between cos(PSI_MAX) and 1
    ev2 = pilot.draw_events(_FixedU(us, 1), 101)
    cpsi = np.cos(ev2[:, 2])
    want2 = math.cos(sz.PSI_MAX) + us * (1.0 - math.cos(sz.PSI_MAX))
    assert np.allclose(cpsi, want2, rtol=1e-12, atol=1e-15)


def test_events_land_inside_the_box():
    rng = np.random.default_rng(40_000_311)   # the retired smoke stream
    ev = pilot.draw_events(rng, 4096)
    t, r, psi, phi = ev[:, 0], ev[:, 1], ev[:, 2], ev[:, 3]
    assert (0 <= t).all() and (t <= sz.DT).all()
    assert (sz.R_LO <= r).all() and (r <= sz.R_HI).all()
    assert (0 <= psi).all() and (psi <= sz.PSI_MAX).all()
    assert (0 <= phi).all() and (phi < 2 * math.pi).all()


# ------------------------------------------------ the frozen rule

def test_the_rule_is_general_in_k_and_matches_the_reference_table():
    """The frozen reference evaluations, from the doc -- but the RULE
    computes them; nothing is hardcoded to k <= 1."""

    expect = {
        0: (11.5092, 1.474933e-06, "FEASIBLE"),
        1: (16.5905, 9.999992e-04, "FEASIBLE"),
        2: (21.0081, 2.426060e-02, "INCONCLUSIVE"),
        3: (25.1047, 1.414958e-01, "INCONCLUSIVE"),
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


def test_cp_upper_cross_checked_by_the_p2_engine_and_gmpy2():
    """Two independent exact engines reproduce the rule's numbers:
    the frozen P14 P2 gamma engine (Poisson tail via its gamma link)
    and a 96-bit gmpy2 evaluation of both the CP bound and the tail."""

    import gmpy2
    import p14_probe_p2 as p2

    gmpy2.get_context().precision = 96
    n, alpha = pilot.FROZEN["n_events"], pilot.FROZEN["alpha_pilot"]
    for k in (0, 1, 2, 5):
        pu = pilot.cp_upper(k, n, alpha)
        # gmpy2 re-derivation of the CP bound
        if k == 0:
            pu3 = float(1 - gmpy2.exp(
                gmpy2.log(gmpy2.mpfr(alpha)) / n))
        else:
            lo, hi = gmpy2.mpfr(k) / n, gmpy2.mpfr(1)
            for _ in range(300):
                mid = (lo + hi) / 2
                lp, lq = gmpy2.log(mid), gmpy2.log1p(-mid)
                s = gmpy2.mpfr(0)
                for i in range(k + 1):
                    s += gmpy2.exp(
                        gmpy2.lgamma(gmpy2.mpfr(n + 1))[0]
                        - gmpy2.lgamma(gmpy2.mpfr(i + 1))[0]
                        - gmpy2.lgamma(gmpy2.mpfr(n - i + 1))[0]
                        + i * lp + (n - i) * lq)
                if s > alpha:
                    lo = mid
                else:
                    hi = mid
            pu3 = float((lo + hi) / 2)
        # the double engine's floor is the lgamma cancellation at
        # n ~ 1e7 (relative ~1e-8); the agreement check sits above
        # that floor, and the VERDICT-level agreement below is exact
        assert abs(pu - pu3) < 1e-7 * pu + 1e-300
        # Poisson tail via the frozen P2 gamma link
        lam = pilot.FROZEN["a_provisional"] * sz.SCALE * pu
        t1 = pilot.pois_tail_gt(pilot.FROZEN["u_max"], lam)
        t2 = 1.0 - p2._poisson_cdf(pilot.FROZEN["u_max"], lam)
        assert abs(t1 - t2) < 1e-12
        # and the VERDICT the engines imply is identical
        lam3 = pilot.FROZEN["a_provisional"] * sz.SCALE * pu3
        t3 = pilot.pois_tail_gt(pilot.FROZEN["u_max"], lam3)
        assert ((t1 <= pilot.FROZEN["tail_budget"])
                == (t3 <= pilot.FROZEN["tail_budget"]))


def test_the_verdict_boundary_is_where_the_rule_puts_it():
    """k_max at the frozen n is DERIVED (k=1), not assumed: k=1 passes
    by 8e-10 of tail budget and k=2 fails by 23x."""

    assert pilot.decide(1)["tail_p_u_gt_u_max"] <= 0.001
    assert pilot.decide(2)["tail_p_u_gt_u_max"] > 0.001


# ------------------------------------------------ gates and refusals

@pytest.mark.parametrize("bad", ["30d286a", "a" * 39, "g" * 40, ""])
def test_short_or_malformed_freeze_rev_is_refused(bad):
    with pytest.raises(SystemExit, match="full 40-hex"):
        pilot.require_full_sha(bad)


def test_missing_freeze_rev_is_refused_by_the_cli():
    script = str(_REPO / "experiments" / "oracle" / "o5_amb_pilot.py")
    out = subprocess.run([sys.executable, script, "--preflight"],
                         capture_output=True, text=True)
    assert out.returncode != 0
    assert "full 40-hex" in (out.stderr + out.stdout)


def _pass_static(monkeypatch, tmp_path, rev=_SHA):
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
    assert pilot.preflight(_SHA)["seed"] == 40_000_441
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
                        lambda tol: {"passed": True})
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
    assert ck["rng_stream"] == "o5_amb_pilot"
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


# ------------------------------------------------ manifest and doc

def test_the_committed_manifest_is_the_tree_it_describes():
    committed = json.loads(pilot._MANIFEST.read_text(encoding="utf-8"))
    assert committed == pilot.build_manifest()
    assert set(committed["files"]) == set(pilot.PROTOCOL_SURFACE)


def test_the_doc_freezes_the_discipline():
    doc = (_REPO / "docs" / "prereg"
           / "p14_o5_amb_pilot.md").read_text(encoding="utf-8")
    flat = " ".join(doc.split())
    assert "no auto-raise, ever" in flat
    assert "general in k" in flat.lower()
    assert "exact INEQUALITY upper bound" in flat
    assert "Wall-cap risk is NOT included" in flat
    assert "never twice" in flat
    assert "INCONCLUSIVE indicts nothing" in flat
