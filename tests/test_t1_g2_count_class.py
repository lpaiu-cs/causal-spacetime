"""Regression tests for the G2 count-class close-out.

Sizing note: the two arms need different budgets. The tube harvest is
cheap per realization but its fit scatters, so it gets the full density
grid; the order-only harvest is the expensive one but fits tightly, so
it gets a shorter grid with more trials. Both are chosen to decide the
1/6 separation with several standard errors to spare while keeping this
file to well under a minute. Full-precision numbers live in
docs/theory/t1_g2_count_class_results.json.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = Path(__file__).resolve().parents[1] / "experiments" / "theory"
sys.path.insert(0, str(EXPERIMENT_DIR))

from t1_g2_count_class import (  # noqa: E402
    CANDIDATES,
    RHO_GRID,
    fit_exponent_with_uncertainty,
    measure_arm,
)

SHORT_GRID = (1000, 2000, 4000, 8000, 16000)
KPZ = 1.0 / 3.0
POISSON = 0.5


def test_exponent_fit_matches_the_analytic_ols_answer():
    """The interval machinery against a hand-computable case:
    log-abscissae 0..3 with residuals (+d, -d, -d, +d) are orthogonal to
    both the constant and the regressor, so OLS returns the planted
    slope and exactly those residuals, giving se = d sqrt(2/5)."""

    d, theta = 0.1, KPZ
    log_x = np.array([0.0, 1.0, 2.0, 3.0])
    residual = np.array([d, -d, -d, d])

    fit = fit_exponent_with_uncertainty(
        list(np.exp(log_x)), list(np.exp(theta * log_x + residual))
    )
    assert abs(fit["theta"] - theta) < 1e-12
    assert abs(fit["stderr"] - d * np.sqrt(2.0 / 5.0)) < 1e-12
    assert fit["dof"] == 2
    assert "halves" not in fit  # four points is too few to split


def test_fit_refuses_too_few_densities():
    """Two points determine a line exactly, leaving no residual to build
    an interval from; that must be an explicit error, not a division by
    a zero degree-of-freedom count."""

    with pytest.raises(ValueError, match="at least 3 densities"):
        fit_exponent_with_uncertainty([1.0, 2.0], [1.0, 2.0])


def test_calibration_arm_returns_one_half():
    """The Poisson clock's count exponent is 1/2 by construction. If the
    estimator misses that, every verdict below is meaningless."""

    measured = measure_arm("thinned", rho_grid=RHO_GRID, trials=200)
    fit = measured["fit"]

    assert measured["clock_failures"] == 0
    assert abs(fit["theta"] - POISSON) < 0.04, fit["theta"]
    assert fit["candidates_inside_ci95"]["poisson_rate_1/2"]
    assert not fit["candidates_inside_ci95"]["kpz_like_1/3"]


def test_order_only_chain_counts_are_kpz_not_poisson():
    """A genuine longest chain, free to wander: Tracy-Widom counts."""

    measured = measure_arm("order_only", rho_grid=SHORT_GRID, trials=200)
    fit = measured["fit"]

    assert measured["clock_failures"] == 0
    assert abs(fit["theta"] - KPZ) < abs(fit["theta"] - POISSON)
    assert fit["t_against"]["poisson_rate_1/2"] < -3.0, fit["t_against"]


def test_tube_chain_counts_are_poisson_not_kpz():
    """The shipped tube is ~rho^-1/2 wide against a natural wandering of
    ~rho^-1/6, so it confines the chain far below the excursions KPZ
    scaling needs and the counts revert to Poisson. This is the finding
    that retires reading the tube's -0.317 error exponent as its count
    class."""

    measured = measure_arm("tube", rho_grid=RHO_GRID, trials=150)
    fit = measured["fit"]

    assert abs(fit["theta"] - POISSON) < abs(fit["theta"] - KPZ)
    assert fit["t_against"]["kpz_like_1/3"] > 3.0, fit["t_against"]


def test_the_two_harvest_protocols_land_in_different_classes():
    """The headline: there is no single count class to find. On the same
    densities, the two harvests come out on opposite sides of the gap --
    which is why asking for "the" class never converged."""

    order_only = measure_arm("order_only", rho_grid=SHORT_GRID, trials=120)
    tube = measure_arm("tube", rho_grid=SHORT_GRID, trials=120)

    separation = abs(CANDIDATES["poisson_rate_1/2"] - CANDIDATES["kpz_like_1/3"])
    assert abs(separation - 1.0 / 6.0) < 1e-12

    gap = tube["fit"]["theta"] - order_only["fit"]["theta"]
    assert gap > 0.5 * separation, (order_only["fit"]["theta"],
                                    tube["fit"]["theta"])


def test_widening_the_tube_restores_kpz_counts():
    """The mechanism, not just the endpoints: a tube that keeps up with
    the chain's natural transverse wandering (~rho^-1/6) stops
    suppressing the optimization, and theta returns toward 1/3."""

    grid = RHO_GRID[:6]
    tight = measure_arm(
        "tube", rho_grid=grid, trials=120,
        tube_width_at=lambda rho: 3.0 * rho ** -0.5,
    )
    wide = measure_arm(
        "tube", rho_grid=grid, trials=120,
        tube_width_at=lambda rho: 1.0 * rho ** (-1.0 / 6.0),
    )

    assert wide["fit"]["theta"] < tight["fit"]["theta"]
    assert abs(wide["fit"]["theta"] - KPZ) < abs(tight["fit"]["theta"] - KPZ)


def test_clock_failures_are_counted_rather_than_dropped():
    """A harvest that raises must be recorded, not silently read as a
    short chain -- the same guard the other T1 harnesses apply."""

    measured = measure_arm("order_only", rho_grid=(1000, 2000, 4000), trials=20)

    assert measured["clock_failures"] == 0
    assert all(row["trials"] == 20 for row in measured["rows"])
