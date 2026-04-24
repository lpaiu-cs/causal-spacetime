from __future__ import annotations

import pytest

from causal_spacetime_lab.metrics import (
    estimate_tau_from_interval_count,
    minkowski_tau_1p1,
)


def test_minkowski_tau_returns_correct_timelike_value() -> None:
    assert minkowski_tau_1p1((0.0, 0.0), (5.0, 3.0)) == pytest.approx(4.0)


def test_minkowski_tau_returns_zero_for_null_separation() -> None:
    assert minkowski_tau_1p1((0.0, 0.0), (1.0, 1.0)) == pytest.approx(0.0)


def test_minkowski_tau_raises_for_spacelike_events() -> None:
    with pytest.raises(ValueError, match="spacelike-separated"):
        minkowski_tau_1p1((0.0, 0.0), (1.0, 2.0))


def test_interval_count_estimator_returns_expected_value() -> None:
    assert estimate_tau_from_interval_count(
        count=50,
        rho=25,
        eta_d=0.5,
        d=2,
    ) == pytest.approx(2.0)

