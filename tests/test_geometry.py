from __future__ import annotations

import math

import pytest

from causal_spacetime_lab.estimators import estimate_tau_from_interval_count_1p1
from causal_spacetime_lab.geometry import (
    alexandrov_eta,
    alexandrov_volume,
    estimate_tau_from_interval_count,
)


def test_alexandrov_eta_known_values() -> None:
    assert alexandrov_eta(2) == pytest.approx(0.5)
    assert alexandrov_eta(3) == pytest.approx(math.pi / 12.0)
    assert alexandrov_eta(4) == pytest.approx(math.pi / 24.0)


def test_alexandrov_volume_formula() -> None:
    assert alexandrov_volume(2.0, spacetime_dim=3) == pytest.approx(2 * math.pi / 3)


def test_generalized_tau_estimator_agrees_with_1p1_estimator() -> None:
    generalized = estimate_tau_from_interval_count(
        interval_count=50,
        rho=25.0,
        spacetime_dim=2,
    )
    existing = estimate_tau_from_interval_count_1p1(
        interval_count=50,
        rho=25.0,
    )

    assert generalized == pytest.approx(existing)

