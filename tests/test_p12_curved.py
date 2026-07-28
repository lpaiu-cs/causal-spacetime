"""Regressions for P12 (prereg v1.1a, Section 8).

The load-bearing pins: the truth formula against numerical geodesic
integration, the sprinkler's realized density against Omega^2, the
null-convention identity that v1.0 got wrong, and the frozen
eligibility conditions that v1.1 added after the completeness
collapse.
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

from p12_curved import (  # noqa: E402
    ELL,
    ETA_RANGE,
    P12_LADDER,
    PATCH_VOLUME,
    PILOT_BLOCKS,
    STAGE_A_BLOCKS,
    STRIDE,
    TAU_BAND,
    VERIFY_BASE,
    X_HALF,
    eligible_pairs,
    run_sample,
    sprinkle,
    tau_curved,
)


def _geodesic_proper_time(e1, x1, e2, x2):
    """INDEPENDENT geodesic integration for ds^2 = Omega^2(-deta^2+dx^2)
    with Omega = ell/|eta| -- no de Sitter invariant, no embedding.

    x is cyclic, so k = Omega^2 dx/dtau is conserved; the normalization
    Omega^2(deta^2 - dx^2) = 1 then gives

        dx/deta   = k / sqrt(Omega^2 + k^2),
        dtau/deta = Omega^2 / sqrt(Omega^2 + k^2).

    Shoot on k until the integrated Delta x matches the endpoints, then
    integrate the proper time. This is a solution of the geodesic
    equation, not a straight coordinate path, and the prereg requires
    agreement at 1e-9.
    """

    nodes, weights = np.polynomial.legendre.leggauss(512)
    lo, hi = min(e1, e2), max(e1, e2)
    grid = 0.5 * (hi - lo) * nodes + 0.5 * (hi + lo)
    scale = 0.5 * (hi - lo)
    omega2 = ELL ** 2 / grid ** 2
    target = abs(x2 - x1)

    def dx_of_k(k):
        return float(np.sum(weights * (k / np.sqrt(omega2 + k ** 2))) * scale)

    if target == 0.0:
        k = 0.0
    else:
        k_hi = 1.0
        while dx_of_k(k_hi) < target:
            k_hi *= 2.0
            assert k_hi < 1e12, "no geodesic reaches this separation"
        k_lo = 0.0
        for _ in range(200):
            mid = 0.5 * (k_lo + k_hi)
            if dx_of_k(mid) < target:
                k_lo = mid
            else:
                k_hi = mid
        k = 0.5 * (k_lo + k_hi)
    return float(
        np.sum(weights * (omega2 / np.sqrt(omega2 + k ** 2))) * scale
    )


def test_truth_matches_numerical_geodesic_integration():
    """Prereg Section 3: the closed form must agree with an independent
    integration of the GEODESIC equation to 1e-9 relative, on a frozen
    grid of separations (the first version integrated a straight
    coordinate path at 2e-3, which could have accepted a materially
    wrong truth formula)."""

    cases = [
        (-1.9, 0.0, -1.6, 0.05),
        (-1.8, -0.2, -1.5, -0.1),
        (-1.5, 0.3, -1.2, 0.35),
        (-2.0, 0.0, -1.7, 0.0),
        (-1.95, -0.4, -1.05, 0.4),
    ]
    for e1, x1, e2, x2 in cases:
        closed = float(tau_curved(e1, x1, e2, x2))
        numeric = _geodesic_proper_time(e1, x1, e2, x2)
        assert np.isfinite(closed)
        assert abs(closed - numeric) / closed < 1e-9


def test_truth_reduces_to_the_local_metric_at_short_separation():
    e, de, dx = -1.5, 2e-3, 1e-3
    closed = float(tau_curved(e, 0.0, e + de, dx))
    local = (ELL / abs(e)) * np.sqrt(de ** 2 - dx ** 2)
    assert abs(closed - local) / local < 1e-3


def test_null_convention_identity_the_v1_0_bug():
    """FULL-null here, HALF-null in the P11 routines: the same number.
    v1.0 quoted the half-null formula against full-null coordinates,
    inflating the flat arm by exactly this factor."""

    e1, x1, e2, x2 = -1.8, 0.0, -1.5, 0.1
    dU = (e2 + x2) - (e1 + x1)
    dV = (e2 - x2) - (e1 - x1)
    du, dv = dU / 2.0, dV / 2.0
    assert np.sqrt(dU * dV) == pytest.approx(2.0 * np.sqrt(du * dv))


def test_sprinkler_density_follows_omega_squared():
    """Chi-square on the realized eta profile against ell^2/eta^2."""

    rng = np.random.default_rng(2_000_101)
    eta = np.concatenate([sprinkle(2400, rng)[0] for _ in range(8)])
    edges = np.linspace(ETA_RANGE[0], ETA_RANGE[1], 11)
    observed, _ = np.histogram(eta, bins=edges)
    weights = np.array([
        np.trapezoid(ELL ** 2 / np.linspace(a, b, 2001) ** 2,
                     np.linspace(a, b, 2001))
        for a, b in zip(edges[:-1], edges[1:], strict=True)
    ])
    expected = eta.size * weights / weights.sum()
    chi2 = float(np.sum((observed - expected) ** 2 / expected))
    assert chi2 < 27.9      # 9 dof, p = 0.001


def test_patch_volume_matches_closed_form():
    closed = (ELL ** 2 * (1.0 / abs(ETA_RANGE[1]) - 1.0 / abs(ETA_RANGE[0]))
              * 2.0 * X_HALF)
    assert PATCH_VOLUME == pytest.approx(closed, rel=1e-9)


def test_eligibility_enforces_band_and_box_inside_patch():
    """Both frozen conditions (v1.1): the box-inside rule is what
    stops chains being truncated by the patch edge."""

    rng = np.random.default_rng(2_000_202)
    eta, x = sprinkle(1200, rng)
    pool = eligible_pairs(eta, x)
    assert pool.shape[0] > 0
    a, b = pool[:, 0], pool[:, 1]
    tau = tau_curved(eta[a], x[a], eta[b], x[b])
    assert np.all(tau >= TAU_BAND[0]) and np.all(tau <= TAU_BAND[1])
    u, v = eta + x, eta - x
    assert np.all(u[b] > u[a]) and np.all(v[b] > v[a])
    e1, x1 = (u[b] + v[a]) / 2.0, (u[b] - v[a]) / 2.0
    e2, x2 = (u[a] + v[b]) / 2.0, (u[a] - v[b]) / 2.0
    for e, xx in ((e1, x1), (e2, x2)):
        assert np.all(e >= ETA_RANGE[0] - 1e-12)
        assert np.all(e <= ETA_RANGE[1] + 1e-12)
        assert np.all(np.abs(xx) <= X_HALF + 1e-12)


def test_completeness_smoke_after_the_v1_1_fix():
    """v1.0's constants completed 12-22% of samples; the pin demands
    1998/2000. This smoke demands 20/20 at test speed."""

    done = sum(
        run_sample(600, VERIFY_BASE[600] + k, single_stream=True)[1]
        for k in range(20)
    )
    assert done == 20


def test_sample_scores_curved_truth_not_flat():
    record, complete = run_sample(1200, VERIFY_BASE[1200], single_stream=True)
    assert complete
    assert np.isfinite(record["y"])
    # the estimator must sit closer to the CURVED truth than to the
    # flat template -- the arm that actually involves tau_hat (the
    # first version compared two truths and never touched the
    # estimator, so it could not support this claim at all)
    assert (record["median_relerr_vs_flat_template"]
            > record["median_relerr_curved"])
    # and the band must carry real curvature to compare against
    assert 0.3 < record["median_flat_curved_gap"] < 0.7


def test_windows_are_private_and_fresh_and_clear_of_design_space():
    spans = []
    for base, slots in (list(PILOT_BLOCKS.values())
                        + list(STAGE_A_BLOCKS.values())):
        span = set(range(base, base + STRIDE * slots))
        spans.append(span)
        assert min(span) >= 1_000_000
        assert max(span) < 2_000_000        # design-check space starts here
    for a in spans:
        for b in spans:
            assert a is b or not (a & b)
    for n in P12_LADDER:
        vspan = set(range(VERIFY_BASE[n], VERIFY_BASE[n] + 2000))
        assert max(vspan) < min(min(s) for s in spans)
