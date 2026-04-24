from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.lorentz import (
    gamma,
    inverse_lorentz_transform_1p1,
    lorentz_transform_1p1,
    measured_rod_length_lab,
    rod_endpoint_events_simultaneous_in_lab,
)


def test_gamma_returns_expected_values() -> None:
    assert gamma(0.0) == pytest.approx(1.0)
    assert gamma(0.6) == pytest.approx(1.25)


def test_gamma_rejects_invalid_beta() -> None:
    with pytest.raises(ValueError, match=r"abs\(beta\) < 1"):
        gamma(1.0)


def test_lorentz_transform_and_inverse_round_trip() -> None:
    event_prime = (2.0, 0.75)
    beta = 0.4

    event_lab = lorentz_transform_1p1(*event_prime, beta=beta)
    round_trip = inverse_lorentz_transform_1p1(*event_lab, beta=beta)

    assert round_trip == pytest.approx(event_prime)


def test_lab_simultaneous_endpoint_events_are_different_rest_frame_times() -> None:
    L0 = 2.0
    beta = 0.5

    left_event, right_event = rod_endpoint_events_simultaneous_in_lab(L0, beta)
    left_prime = inverse_lorentz_transform_1p1(*left_event, beta=beta)
    right_prime = inverse_lorentz_transform_1p1(*right_event, beta=beta)

    assert left_event[0] == pytest.approx(0.0)
    assert right_event[0] == pytest.approx(0.0)
    assert left_prime[1] == pytest.approx(0.0)
    assert right_prime[1] == pytest.approx(L0)
    assert left_prime[0] != pytest.approx(right_prime[0])


@pytest.mark.parametrize("beta", [0.0, 0.3, 0.6, -0.6])
def test_measured_rod_length_lab_matches_l0_over_gamma(beta: float) -> None:
    L0 = 3.0

    measured_length = measured_rod_length_lab(L0, beta)

    assert measured_length == pytest.approx(L0 / gamma(beta))


def test_gamma_accepts_arrays() -> None:
    betas = np.array([0.0, 0.6])

    values = gamma(betas)

    assert values == pytest.approx(np.array([1.0, 1.25]))

