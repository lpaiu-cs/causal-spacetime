"""Regressions for the P11 metric instrument (prereg v1.6, Section 8).

The load-bearing pins: the chain estimator must equal the shared
definition's (L - 2)/sqrt(N) under the frozen rho = N/2 convention,
the verdict table must be single-valued by precedence, the power
literals must implement the frozen selection rule, and every seed
window must be private and fresh against the full documented list.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p11_metric import (  # noqa: E402
    DELTA_EQ,
    K_PAIRS,
    P11_LADDER,
    PILOT_BLOCKS,
    SCORE_OFFSET,
    SKIP_CAP,
    STAGE_A_BLOCKS,
    STRIDE,
    TAU_BAND,
    VERIFY_BASE,
    VERIFY_COUNT,
    continuum_uv,
    draw_disjoint_pairs,
    eligible_pool,
    power_requirements,
    run_sample,
    score_pair,
    supports_disjoint,
    verdict,
)


def test_continuum_uv_is_the_rank_grid():
    pi = np.array([2, 0, 3, 1])
    u, v = continuum_uv(pi)
    assert np.allclose(u, np.arange(4) / 4.0)
    assert np.allclose(v, pi / 4.0)


def test_tau_chain_equals_shared_definition_with_endpoints_subtracted():
    """Identity permutation: the interval (0, k) contains k - 1
    interior events forming a full chain, so the closed chain length is
    k + 1 and tau_hat must be exactly (k + 1 - 2)/sqrt(N)."""

    n, k = 100, 40
    pi = np.arange(n)
    u, v = continuum_uv(pi)
    scored = score_pair(u, v, 0, k)
    assert scored["chain_length"] == k + 1
    assert scored["tau_chain"] == pytest.approx((k - 1) / np.sqrt(n))
    # and the volume estimator through the shared metrics definition
    assert scored["m_open"] == k - 1
    assert scored["tau_vol"] == pytest.approx(2.0 * np.sqrt((k - 1) / n))
    # truth via the frozen convention
    assert scored["tau_true"] == pytest.approx(2.0 * k / n)


def test_conditioned_interior_mean_uses_rank_gaps():
    n, k = 100, 40
    pi = np.arange(n)
    u, v = continuum_uv(pi)
    scored = score_pair(u, v, 0, k)
    assert scored["m_conditioned"] == pytest.approx(
        (k - 1) * (k - 1) / (n - 2)
    )


def test_eligible_pool_matches_band_by_construction():
    rng = np.random.default_rng(7)
    pi = rng.permutation(300)
    u, v = continuum_uv(pi)
    pool = eligible_pool(u, v)
    assert pool.shape[0] > 0
    taus = 2.0 * np.sqrt(
        (u[pool[:, 1]] - u[pool[:, 0]]) * (v[pool[:, 1]] - v[pool[:, 0]])
    )
    assert np.all(taus >= TAU_BAND[0]) and np.all(taus <= TAU_BAND[1])
    assert np.all(u[pool[:, 1]] > u[pool[:, 0]])
    assert np.all(v[pool[:, 1]] > v[pool[:, 0]])


def test_supports_disjoint_is_box_disjointness():
    a = (0.0, 0.2, 0.0, 0.2)
    assert supports_disjoint(a, (0.3, 0.5, 0.0, 0.2))   # u-disjoint
    assert supports_disjoint(a, (0.0, 0.2, 0.3, 0.5))   # v-disjoint
    assert not supports_disjoint(a, (0.1, 0.3, 0.1, 0.3))  # overlap


def test_draw_rejection_cap_marks_incomplete():
    """A pool of mutually overlapping candidates can never yield six
    disjoint supports; the frozen cap must fire and mark the sample
    incomplete rather than loop."""

    u = np.linspace(0.0, 1.0, 40)
    v = np.linspace(0.0, 1.0, 40)
    pool = np.array([[0, 39 - k] for k in range(10)])  # nested boxes
    pairs, complete, rejections = draw_disjoint_pairs(
        pool, u, v, np.random.default_rng(0)
    )
    assert not complete
    assert len(pairs) == 1
    assert rejections == 9  # pool exhausted before the cap


def test_run_sample_completes_and_scores_at_bottom_rung():
    record, complete = run_sample(600, 190000, single_stream=True)
    assert complete
    assert record["pool_size"] > 0
    assert np.isfinite(record["y"])
    assert record["min_m_open"] > 0


def test_verification_mode_smoke_completeness():
    """20 single-stream samples at the bottom rung: the Section 4 pin
    demands 1998/2000; this smoke demands 20/20 at test speed (a
    failure here means the full pin is hopeless)."""

    complete_count = sum(
        run_sample(600, 190000 + k, single_stream=True)[1]
        for k in range(20)
    )
    assert complete_count == 20


def test_verdict_precedence_is_single_valued():
    assert verdict(-0.05, -0.01, True) == "IMPROVES"      # also in margin
    assert verdict(0.01, 0.05, True) == "DEGRADES"        # also in margin
    assert verdict(-0.05, 0.05, True) == "FLAT-WITHIN-MARGIN"
    assert verdict(-0.10, 0.05, True) == "INCONCLUSIVE"
    assert verdict(-0.05, 0.05, False) == "UNRESOLVED"
    assert verdict(-0.05, -0.01, False) == "IMPROVES"
    assert abs(DELTA_EQ - 0.067) < 1e-12


def _gaussian_pilot(sigma: float, seed: int, n: int = 200) -> np.ndarray:
    y = np.random.default_rng(seed).standard_normal(n)
    y = (y - y.mean()) / y.std(ddof=1)          # exact sample sigma
    return sigma * y


def test_bonett_bound_is_above_the_variance_and_kurtosis_adaptive():
    from p11_metric import bonett_variance_bound

    gaussian = bonett_variance_bound(_gaussian_pilot(0.10, 1))
    assert gaussian["s2"] == pytest.approx(0.01)
    # Gaussian-kurtosis inflation sits near 1.18, never below 1
    assert 1.05 < gaussian["bound"] / gaussian["s2"] < 1.35

    # heavy tails must widen the bound: same variance, Laplace-ish y
    rng = np.random.default_rng(2)
    heavy = rng.laplace(size=200)
    heavy = 0.10 * (heavy - heavy.mean()) / heavy.std(ddof=1)
    heavy_bound = bonett_variance_bound(heavy)
    assert heavy_bound["g4"] > gaussian["g4"]
    assert (heavy_bound["bound"] / heavy_bound["s2"]
            > gaussian["bound"] / gaussian["s2"])


def test_power_selection_rule_branches():
    """v1.8: sizing from summed per-rung Bonett bounds; branch
    structure exercised with exact-sample-sigma Gaussian pilots."""

    # small sigmas -> n_eq drives, flat available
    p = power_requirements(_gaussian_pilot(0.07, 3), _gaussian_pilot(0.07, 4))
    assert p["flat_available"] and p["n_per_rung"] == p["n_eq"]
    assert p["n_eq"] <= 60 and not p["infeasible"]
    # moderate sigmas -> flat unaffordable, superiority sizing
    p = power_requirements(_gaussian_pilot(0.15, 5), _gaussian_pilot(0.15, 6))
    assert not p["flat_available"] and not p["infeasible"]
    assert p["n_per_rung"] == max(p["n_sup"], 12)
    # large sigmas -> superiority breaches the cap -> infeasible
    p = power_requirements(_gaussian_pilot(0.40, 7), _gaussian_pilot(0.40, 8))
    assert p["infeasible"] and p["n_per_rung"] is None
    # variances enter as sums; the bound exceeds the point estimate
    p = power_requirements(_gaussian_pilot(0.05, 9), _gaussian_pilot(0.13, 10))
    assert p["s2"] == pytest.approx(0.05 ** 2 + 0.13 ** 2)
    assert p["s2_90"] > p["s2"]


def test_stage_gate_rejects_stale_prerequisite(tmp_path, monkeypatch):
    """A prerequisite artifact stamped by an older implementation (or
    left behind by a crashed rerun) must be refused, not consumed."""

    import p11_metric as mod

    (tmp_path / "p11_verification_summary.json").write_text(
        '{"code_version": "0ld5tamp", "pin_passed": true}\n',
        encoding="utf-8",
    )
    monkeypatch.setattr(mod, "_diag_code_version", lambda: "abc1234")
    with pytest.raises(SystemExit, match="implementation changed"):
        mod.run_pilot(tmp_path)


def test_preflight_refuses_a_dirty_worktree(monkeypatch, tmp_path):
    import p11_metric as mod

    monkeypatch.setattr(mod, "_diag_code_version",
                        lambda: "abc1234-dirty")
    with pytest.raises(SystemExit, match="dirty"):
        mod.run_verify(tmp_path)


def test_windows_are_private_fresh_and_reserve_sized():
    """Stride-200 privacy, the documented spans, and freshness against
    every seed range the programme has ever used."""

    used = (set(range(0, 10)) | set(range(100, 120)) | set(range(400, 420))
            | set(range(500, 520)) | set(range(820, 960))
            | set(range(1000, 1020)) | set(range(1100, 1160))
            | set(range(9001, 9060)) | set(range(30000, 30380))
            | set(range(40000, 40169)) | set(range(41000, 41060))
            | set(range(43000, 55000)))
    blocks = []
    for n, (base, slots) in {**PILOT_BLOCKS, **STAGE_A_BLOCKS}.items():
        span = set(range(base, base + STRIDE * slots))
        blocks.append((n, base, span))
        assert not (span & used), (n, base)
        assert max(0, 3, 5, 9, 61, 100, SCORE_OFFSET) < STRIDE
    for _, _, a in blocks:
        for _, _, b in blocks:
            assert a is b or not (a & b)
    assert PILOT_BLOCKS[600] == (200000, 220)
    assert PILOT_BLOCKS[2400] == (244000, 220)
    assert STAGE_A_BLOCKS == {600: (288000, 80), 1200: (304000, 80),
                              2400: (320000, 80)}
    experimental_bottom = min(
        base for base, _ in {**PILOT_BLOCKS, **STAGE_A_BLOCKS}.values()
    )
    for n in P11_LADDER:
        verify_span = set(range(VERIFY_BASE[n], VERIFY_BASE[n] + VERIFY_COUNT))
        assert max(verify_span) < experimental_bottom
        assert not (verify_span & used)
    assert SKIP_CAP == 20 and K_PAIRS == 6
