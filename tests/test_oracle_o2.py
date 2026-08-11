"""PR-O2 contract tests: the certified volume integrator (L6).

The load-bearing check is the flat closed-form control -- the
assembled machinery (pruning, weights with the (rho, psi) Jacobian,
Lipschitz/tangent cells, directed composition) must produce an
interval CONTAINING pi tau^4 / 24 computed as a 200-bit directed
reference. Everything else is fail-closed behavior and soundness of
the closed-form functionals against certified flight times."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import gmpy2
import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import probe_seed_ledger as ledger  # noqa: E402
import volume_oracle as vo  # noqa: E402
from certified_flight_time import (  # noqa: E402
    FROZEN_ANCHORS,
    angle_cost_iv,
    flight_time_certified,
    tortoise_iv,
)
from certified_interval import CertificationError, Iv  # noqa: E402
from volume_oracle import (  # noqa: E402
    OracleConfig,
    assemble,
    mc_diagnostic,
)

_DN200 = gmpy2.context(precision=200, round=gmpy2.RoundDown)
_UP200 = gmpy2.context(precision=200, round=gmpy2.RoundUp)
_PRICE = _REPO / "docs" / "prereg" / "p14_oracle_price.json"
_MODEW = _REPO / "docs" / "prereg" / "p14_oracle_mode_width.json"


def _flat_exact_bounds(dt: float, d: float) -> tuple:
    """200-bit directed enclosure of pi tau^4 / 24 with
    tau^2 = dt^2 - d^2 (the flat Alexandrov-diamond volume)."""

    def val(ctx):
        tau2 = ctx.sub(ctx.mul(gmpy2.mpfr(dt), gmpy2.mpfr(dt)),
                       ctx.mul(gmpy2.mpfr(d), gmpy2.mpfr(d)))
        return ctx.div(ctx.mul(ctx.const_pi(),
                               ctx.mul(tau2, tau2)), gmpy2.mpfr(24))

    return val(_DN200), val(_UP200)


def test_flat_control_contains_the_closed_form():
    """M = 0, anchors (12, 18, 8.5): the certified interval must
    contain pi tau^4 / 24 = pi 36.25^2 / 24, and be nontrivially
    tight at a modest budget."""

    cfg = OracleConfig(12.0, 18.0, 8.5, m=0.0, target_ratio=0.5,
                       n_sub=16, max_calls=500, max_wall_s=240.0,
                       init_rho=6, init_psi=6)
    res = assemble(cfg)
    lo, hi = _flat_exact_bounds(8.5, 6.0)
    v = res["v"]
    assert v.lo <= lo and hi <= v.hi, (float(v.lo), float(v.hi))
    assert res["ratio"] is not None and res["ratio"] < 1.0
    assert v.certainly_gt(Iv(0))


def test_cost_caps_bind_during_the_initial_grid_and_stay_sound():
    """R1: the caps must bind while the INITIAL grid is being
    evaluated, not only inside the refinement loop -- otherwise a
    small cap is blown through before the first target check. The
    zero-call fallback the capped cells fall back to must keep the
    enclosure valid, so the flat closed form must still be inside."""

    cfg = OracleConfig(12.0, 18.0, 8.5, m=0.0, target_ratio=0.001,
                       n_sub=16, max_calls=20, max_wall_s=240.0,
                       init_rho=8, init_psi=8)
    res = assemble(cfg)
    assert res["calls"] <= 22, res["calls"]   # 2 calls per cell
    assert res["modes"]["uncosted"] > 0
    lo, hi = _flat_exact_bounds(8.5, 6.0)
    assert res["v"].lo <= lo and hi <= res["v"].hi


@pytest.mark.parametrize("max_calls,target", [(320, 0.45), (520, 0.38)])
def test_a_crossing_before_the_cap_is_never_misfiled(max_calls,
                                                     target):
    """R1: target accounting happens on a 32-split cadence, so a
    crossing in the final splits before a cap must still be credited
    -- otherwise the ladder records a rung it actually reached as an
    extrapolation. The invariant is exact in both directions."""

    cfg = OracleConfig(12.0, 18.0, 8.0, m=1.0, n_sub=16,
                       max_calls=max_calls, max_wall_s=240.0,
                       init_rho=8, init_psi=8)
    res = assemble(cfg, targets=[target])
    if res["ratio"] is not None and res["ratio"] <= target:
        assert res["status"] == "target-met"
        assert target in res["crossings"]
        assert res["crossings"][target]["ratio"] <= target
    else:
        assert res["status"] == "target-not-met"
        assert target not in res["crossings"]


def test_empty_diamond_returns_exact_zero():
    res = assemble(OracleConfig(12.0, 18.0, 5.0, m=1.0))
    assert res["status"] == "empty-diamond"
    assert float(res["v"].lo) == 0.0 and float(res["v"].hi) == 0.0
    assert res["calls"] == 0
    assert res["termination_reason"] == "empty-diamond"


# ---------------------------------------------------------------------
# termination provenance (O3 freeze-review ruling): every exit path
# must name its reason -- a 24 h frozen run whose caps fire is not
# fail-explicit if all caps collapse into one `target-not-met`
# ---------------------------------------------------------------------


def test_termination_reason_target_met():
    cfg = OracleConfig(12.0, 18.0, 8.0, m=1.0, target_ratio=0.9,
                       n_sub=16, max_calls=4000, max_wall_s=240.0,
                       init_rho=6, init_psi=6)
    res = assemble(cfg)
    assert res["status"] == "target-met"
    assert res["termination_reason"] == "target-met"


def test_termination_reason_max_calls():
    cfg = OracleConfig(12.0, 18.0, 8.0, m=1.0, target_ratio=1e-9,
                       n_sub=16, max_calls=30, max_wall_s=240.0,
                       init_rho=6, init_psi=6)
    res = assemble(cfg)
    assert res["status"] == "target-not-met"
    assert res["termination_reason"] == "max-calls"
    assert res["uncosted_cells"] > 0
    assert res["calls"] >= 30


def test_termination_reason_max_wall():
    cfg = OracleConfig(12.0, 18.0, 8.0, m=1.0, target_ratio=1e-9,
                       n_sub=16, max_calls=10_000, max_wall_s=1e-6,
                       init_rho=6, init_psi=6)
    res = assemble(cfg)
    assert res["status"] == "target-not-met"
    assert res["termination_reason"] == "max-wall"
    assert res["calls"] < 10_000


def test_termination_reason_max_depth_exhausted():
    """max_depth=0 blocks every split: the heap drains with live
    unpruned cells sitting at the depth cap, and the result must say
    the DEPTH was the limit -- the ruling explicitly forbids
    assuming a depth cap is sufficient without the result showing
    when it was not."""

    cfg = OracleConfig(12.0, 18.0, 8.0, m=1.0, target_ratio=1e-9,
                       n_sub=16, max_calls=4000, max_wall_s=240.0,
                       init_rho=6, init_psi=6, max_depth=0)
    res = assemble(cfg)
    assert res["status"] == "target-not-met"
    assert res["termination_reason"] == "max-depth-exhausted"
    assert res["cells_at_max_depth"] > 0
    assert res["max_depth_reached"] == 0


def test_termination_reason_heap_exhausted(monkeypatch):
    """The defensive branch: a drained heap with NO depth-blocked
    live cells. L6c makes this unreachable naturally (a cell meeting
    the diamond never prunes), so the classification logic is pinned
    with a doctored lower bound that prunes everything."""

    monkeypatch.setattr(
        vo, "_closed_form_lower",
        lambda st, c: vo.Iv(1e6))
    cfg = OracleConfig(12.0, 18.0, 8.0, m=1.0, target_ratio=1e-9,
                       n_sub=16, max_calls=4000, max_wall_s=240.0,
                       init_rho=4, init_psi=4)
    res = assemble(cfg)
    assert res["termination_reason"] == "heap-exhausted"
    assert res["cells_at_max_depth"] == 0
    assert res["modes"]["pruned"] == 16


def test_containment_violation_refuses_without_numbers():
    with pytest.raises(CertificationError):
        assemble(OracleConfig(12.0, 18.0, 12.0, m=1.0))


def test_schwarzschild_assembly_is_sound_and_mc_is_not_a_gate():
    """Coarse neighbor-configuration run: certified positive
    interval, tangent cells firing, and the raw/certified width
    decomposition present and internally consistent.

    The MC cross-check is exercised STRUCTURALLY only. Per the
    PR #65 review ruling, overlap must never be a pass/fail
    criterion: a probabilistic comparison can fail against a
    perfectly correct oracle, so gating CI on it contradicts the
    documented \"MC is never a verdict\" contract. Overlap is
    computed and recorded by the runners; nothing asserts it."""

    cfg = OracleConfig(12.0, 18.0, 8.0, m=1.0, target_ratio=0.5,
                       n_sub=16, max_calls=400, max_wall_s=240.0,
                       init_rho=8, init_psi=8)
    res = assemble(cfg)
    v = res["v"]
    assert v.certainly_gt(Iv(0))
    assert res["modes"]["tangent"] > 0
    assert res["modes"]["pruned"] > 0
    # decomposition contract: the certified width never exceeds the
    # raw sum, and the per-mode split decomposes the raw sum
    assert (res["certified_total_after_intersection"]
            <= res["raw_total_before_intersection"] + 1e-12)
    assert isinstance(res["intersection_active"], bool)
    split = sum(res["raw_width_by_mode"].values())
    import math
    assert math.isclose(split * 2.0 * math.pi,
                        res["raw_total_before_intersection"],
                        rel_tol=1e-6)
    mc = mc_diagnostic(cfg, 500, 40_000_271)
    assert mc["se"] > 0 and mc["ci95"][0] < mc["ci95"][1]
    assert mc["n"] == 500 and mc["seed"] == 40_000_271


def test_angle_lower_bound_is_below_certified_flight_times():
    """The L6b' shell-local functional must lower-bound the true
    optical distance: compare against certified T at cell corners."""

    m_i = Iv(1.0)
    st = vo._State(vo.OracleConfig(12.0, 18.0, 8.0), m_i, Iv(8.0),
                   tortoise_iv(Iv(12.0), m_i),
                   tortoise_iv(Iv(18.0), m_i), angle_cost_iv(m_i))
    st.rho_smin = tortoise_iv(Iv(vo.R_MIN), m_i)
    st.w_shell = vo.w_iv(Iv(vo.R_MIN), m_i)
    for psi0, r_corner in ((0.1, 13.0), (0.3, 15.5), (0.6, 12.5)):
        rho_c = float(tortoise_iv(Iv(r_corner), m_i).lo)
        cell = vo._Cell(rho_c, rho_c + 0.4, psi0, psi0 + 0.1, 0)
        d_lo = vo._anchor_distance_lower(st, cell, st.rho_p)
        ft = flight_time_certified(12.0, r_corner, psi0, 1.0,
                                   n_sub=16)
        assert float(d_lo.lo) <= float(ft.t.hi) + 1e-9, (psi0,
                                                         r_corner)


def test_price_artifact_measures_the_neighbor_not_the_frozen():
    """The [TO SIZE] ladder must be measured NEXT DOOR: the frozen
    (12, 18, 8.5) volume stays unobserved until PR-O3."""

    art = json.loads(_PRICE.read_text(encoding="utf-8"))
    na = art["neighbor_anchors"]
    assert (na["r_in"], na["r_out"]) == FROZEN_ANCHORS[:2]
    assert na["dt"] != FROZEN_ANCHORS[2]
    assert "unobserved" in art["scope"]


def test_price_artifact_separates_measured_from_extrapolated():
    """Every reported crossing must be a real crossing (its recorded
    ratio actually at or below the target), and no target may appear
    as both measured and extrapolated -- the artifact must never let
    a fitted number read as a measurement."""

    art = json.loads(_PRICE.read_text(encoding="utf-8"))
    measured = art["measured_crossings"]
    extrap = art["extrapolated_targets"]
    for row in measured:
        assert row["ratio"] <= row["target_ratio"], row
        assert row["v_lo"] <= row["v_hi"]
        assert row["calls"] > 0 and row["cells"] > 0
    m_targets = [r["target_ratio"] for r in measured]
    e_targets = [r["target_ratio"] for r in extrap]
    assert set(m_targets).isdisjoint(e_targets)
    assert m_targets == sorted(m_targets, reverse=True)
    for row in extrap:
        # never a measured-looking key on an extrapolated row, and
        # reachability is explicitly an ESTIMATE (review ruling)
        assert "calls" not in row and "cells" not in row
        if row["estimated_reachable_at_this_n_sub"]:
            assert row["calls_extrapolated"] > 0
            lo, hi = row["calls_extrapolated_range"]
            # the fit-window spread is the model's honest error bar
            # and must bracket the primary number
            assert lo <= row["calls_extrapolated"] <= hi
            assert lo < hi
        else:
            assert "calls_extrapolated" not in row
            assert "floor" in row["reason"]
    if extrap:
        fit = art["convergence_fit"]
        assert fit is not None
        assert fit["log_log_slope"] < 0
        # the extrapolation model must be the FLOOR-AWARE one: a
        # plain power law contradicts the floor diagnostic
        assert "floor" in fit["model"]
        assert (fit["floor_ratio_used"]
                == art["quadrature_floor"]["rows"][0]["floor_ratio"])
        assert len(fit["fit_windows"]) >= 2
        assert "not an execution cap" in fit["note"]


def test_price_artifact_curve_is_monotone_and_consistent():
    """The streamed trace must be non-decreasing in cost and each
    sample's interval must be ordered; the final row must agree with
    the last curve sample."""

    art = json.loads(_PRICE.read_text(encoding="utf-8"))
    curve = art["curve"]
    assert len(curve) >= 3
    for a, b in zip(curve[:-1], curve[1:], strict=True):
        assert b["calls"] >= a["calls"]
        assert b["wall_s"] >= a["wall_s"]
    for s in curve:
        assert s["v_lo"] <= s["v_hi"]
    assert art["final"]["calls"] >= curve[-1]["calls"]
    assert art["final"]["status"] in ("target-met", "target-not-met")


def test_price_artifact_records_the_binding_lever_and_mc_role():
    """The floor table must name which lever binds and stay labelled
    diagnostic; the MC row must stay diagnostic-only with the ledger
    seed that was registered as observed."""

    art = json.loads(_PRICE.read_text(encoding="utf-8"))
    qf = art["quadrature_floor"]
    assert qf["binding_lever_estimate"] in ("quadrature",
                                            "cell-refinement")
    assert "DIAGNOSTIC" in qf["note"]
    # the n_sub=32/64 rows are probe-based models, and the note has
    # to say so rather than let them read as integrator runs
    assert "NOT full-integrator" in qf["note"]
    for row in qf["rows"]:
        assert 0.0 < row["floor_ratio"] < 1.0
        assert row["probes"] > 0
    mc = art["mc_diagnostic"]
    assert mc["seed"] == 40_000_271
    assert mc["seed"] in ledger.spent_scalars()
    assert "diagnostic only" in mc["role"]


def test_price_artifact_floor_falls_quadratically_with_n_sub():
    """The certified arc integrals are a composite MIDPOINT rule, so
    the flight-time width -- and with it the floor -- must fall like
    n_sub^-2. A departure means the quadrature is not behaving as
    the certification claims."""

    art = json.loads(_PRICE.read_text(encoding="utf-8"))
    rows = art["quadrature_floor"]["rows"]
    assert len(rows) >= 2
    for a, b in zip(rows[:-1], rows[1:], strict=True):
        assert b["n_sub"] > a["n_sub"]
        got = (a["mean_flight_time_width"]
               / b["mean_flight_time_width"])
        want = (b["n_sub"] / a["n_sub"]) ** 2
        assert 0.7 * want <= got <= 1.4 * want, (a, b, got, want)


def test_price_artifact_plan_covers_the_frozen_target():
    """Every quadrature setting must carry a verdict on the FROZEN
    0.01 target, with wall clock, so PR-O3 has a configuration
    answer rather than a number to guess."""

    art = json.loads(_PRICE.read_text(encoding="utf-8"))
    plan = art["frozen_target_plan"]
    assert plan["target_ratio"] == 0.01
    floors = {r["n_sub"] for r in art["quadrature_floor"]["rows"]}
    assert {r["n_sub"] for r in plan["rows"]} == floors
    for row in plan["rows"]:
        assert isinstance(row["frozen_target_estimated_reachable"],
                          bool)
        if row["frozen_target_estimated_reachable"]:
            assert row["calls_extrapolated"] > 0
            assert row["wall_s_extrapolated"] > 0
            lo, hi = row["calls_extrapolated_range"]
            assert lo <= row["calls_extrapolated"] <= hi
    # the wording boundary the review froze: recommended plan and
    # model-based expectation, with the mode-width diagnostic as
    # the prerequisite for any O3 budget decision
    assert "RECOMMENDED PLAN" in plan["reading"]
    assert "MODEL-BASED EXPECTATION" in plan["reading"]
    assert "CANDIDATE" in plan["reading"]
    assert "not frozen" in plan["note"]


def test_mode_width_artifact_carries_the_ruled_decomposition():
    """The PR #65 ruling's prerequisite for O3: the bottleneck
    diagnostic must be run on the NEIGHBOR configuration and its
    artifact must carry the exact decomposition fields the ruling
    listed -- raw per-mode split, raw and certified totals,
    intersection flag, per-mode cell counts -- with the raw/
    certified relationship intact."""

    art = json.loads(_MODEW.read_text(encoding="utf-8"))
    na = art["neighbor_anchors"]
    assert (na["r_in"], na["r_out"]) == FROZEN_ANCHORS[:2]
    assert na["dt"] != FROZEN_ANCHORS[2]
    assert "unobserved" in art["scope"]
    fin = art["final"]
    for key in ("raw_width_by_mode", "raw_total_before_intersection",
                "certified_total_after_intersection",
                "intersection_active", "cell_counts_by_mode"):
        assert key in fin, key
    assert (fin["certified_total_after_intersection"]
            <= fin["raw_total_before_intersection"] + 1e-12)
    assert isinstance(fin["intersection_active"], bool)
    shares = fin["raw_width_shares"]
    assert abs(sum(shares.values()) - 1.0) < 1e-9
    assert set(fin["raw_width_by_mode"]) == set(
        fin["cell_counts_by_mode"])
    # the verdict must stay within what one neighbor run supports
    assert art["verdict"].startswith("MEASURED on the neighbor")
    assert "UPPER BOUND" in art["decomposition_note"]
    # the trend is auditable: every curve sample carries the split
    assert len(art["curve"]) >= 3
    for s in art["curve"]:
        assert "raw_width_by_mode" in s and "raw_width" in s
