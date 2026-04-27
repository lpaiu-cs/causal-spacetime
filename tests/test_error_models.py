from __future__ import annotations

import math

import pytest

from causal_spacetime_lab.error_models import (
    binned_summary,
    binomial_interval_mean_1p1,
    binomial_interval_var_1p1,
    binomial_tau_abs_std_delta_1p1,
    binomial_tau_rel_std_delta_1p1,
    bootstrap_mean_ci,
    interval_volume_1p1,
    poisson_interval_mean_1p1,
    poisson_tau_abs_std_delta_1p1,
    poisson_tau_rel_std_delta_1p1,
)


def test_interval_volume_1p1() -> None:
    assert interval_volume_1p1(2.0) == pytest.approx(2.0)


def test_poisson_interval_mean_1p1() -> None:
    assert poisson_interval_mean_1p1(tau=2.0, rho=10.0) == pytest.approx(20.0)


def test_poisson_tau_delta_standard_deviations() -> None:
    assert poisson_tau_abs_std_delta_1p1(rho=50.0) == pytest.approx(0.1)
    assert poisson_tau_rel_std_delta_1p1(tau=2.0, rho=50.0) == pytest.approx(0.05)


def test_binomial_interval_mean_and_variance_1p1() -> None:
    assert binomial_interval_mean_1p1(
        tau=1.0,
        T=2.0,
        n_support=100,
    ) == pytest.approx(25.0)
    assert binomial_interval_var_1p1(
        tau=1.0,
        T=2.0,
        n_support=100,
    ) == pytest.approx(18.75)


def test_binomial_tau_delta_standard_deviations() -> None:
    assert binomial_tau_abs_std_delta_1p1(
        tau=1.0,
        T=2.0,
        n_support=100,
        rho=50.0,
    ) == pytest.approx(math.sqrt(18.75) / 50.0)
    assert binomial_tau_rel_std_delta_1p1(
        tau=1.0,
        T=2.0,
        n_support=100,
        rho=50.0,
    ) == pytest.approx(math.sqrt(18.75) / 50.0)


def test_binned_summary_on_deterministic_data() -> None:
    rows = binned_summary(
        x=[0.1, 0.2, 0.4],
        y=[1.0, 2.0, 4.0],
        bins=[0.0, 0.3, 0.5],
    )

    assert rows[0]["count"] == 2
    assert rows[0]["mean"] == pytest.approx(1.5)
    assert rows[0]["rmse"] == pytest.approx(math.sqrt(2.5))
    assert rows[1]["count"] == 1
    assert rows[1]["mean"] == pytest.approx(4.0)


def test_bootstrap_mean_ci_reproducible_with_fixed_seed() -> None:
    first = bootstrap_mean_ci([1.0, 2.0, 3.0], seed=123, n_resamples=100)
    second = bootstrap_mean_ci([1.0, 2.0, 3.0], seed=123, n_resamples=100)

    assert first == pytest.approx(second)
    assert first[0] == pytest.approx(2.0)

