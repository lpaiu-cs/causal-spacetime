from __future__ import annotations

import math

import pytest

from causal_spacetime_lab.estimators import (
    bias,
    estimate_tau_from_interval_count_1p1,
    estimate_tau_from_longest_chain_1p1,
    global_density_1p1,
    mae,
    relative_rmse,
    rmse,
)


def test_global_density_formula_1p1() -> None:
    assert global_density_1p1(100, T=2.0) == pytest.approx(50.0)


def test_interval_count_tau_estimator_1p1() -> None:
    assert estimate_tau_from_interval_count_1p1(
        interval_count=50,
        rho=25.0,
    ) == pytest.approx(2.0)


def test_longest_chain_tau_estimator_endpoint_convention() -> None:
    assert estimate_tau_from_longest_chain_1p1(
        chain_length=6,
        rho=8.0,
        chain_counts_endpoints=True,
    ) == pytest.approx(1.0)
    assert estimate_tau_from_longest_chain_1p1(
        chain_length=6,
        rho=8.0,
        chain_counts_endpoints=False,
    ) == pytest.approx(1.5)


def test_error_metric_functions() -> None:
    truth = [1.0, 2.0]
    prediction = [2.0, 2.0]

    assert mae(truth, prediction) == pytest.approx(0.5)
    assert rmse(truth, prediction) == pytest.approx(math.sqrt(0.5))
    assert bias(truth, prediction) == pytest.approx(0.5)
    assert relative_rmse(truth, prediction) == pytest.approx(math.sqrt(0.5))

