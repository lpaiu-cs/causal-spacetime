from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.conformal import (
    central_observer_proper_time_1p1,
    conformal_volume_density_1p1,
    constant_profile,
    flat_interval_x_width_1p1,
    flat_profile,
    integrate_conformal_interval_volume_1p1,
    linear_time_profile,
    omega_values,
    sinusoidal_time_profile,
    unweighted_interval_coordinate_volume_estimate_1p1,
    weighted_interval_volume_estimate_1p1,
)


def test_flat_profile_omega_is_one() -> None:
    events = np.array([[-1.0, 0.0], [0.2, 0.1]])

    values = omega_values(flat_profile(), events)

    assert values == pytest.approx(np.ones(2))


def test_constant_profile_omega_is_scale() -> None:
    events = np.array([[-1.0, 0.0], [0.2, 0.1]])

    values = omega_values(constant_profile(2.5), events)

    assert values == pytest.approx(np.full(2, 2.5))


def test_invalid_constant_scale_raises() -> None:
    with pytest.raises(ValueError, match="positive"):
        constant_profile(0.0)


def test_linear_time_profile_rejects_nonpositive_risk() -> None:
    with pytest.raises(ValueError, match="abs\\(amplitude\\)"):
        linear_time_profile(1.0, T=2.0)


def test_sinusoidal_time_profile_rejects_abs_amplitude_at_least_one() -> None:
    with pytest.raises(ValueError, match="abs\\(amplitude\\)"):
        sinusoidal_time_profile(-1.0, T=2.0)


def test_conformal_volume_density_is_omega_squared() -> None:
    events = np.array([[0.0, 0.0], [0.5, 0.0]])

    density = conformal_volume_density_1p1(constant_profile(3.0), events)

    assert density == pytest.approx(np.full(2, 9.0))


def test_flat_interval_x_width_centered_interval() -> None:
    t_values = np.array([-1.0, 0.0, 1.0])

    widths = flat_interval_x_width_1p1(
        t_values,
        np.array([-1.0, 0.0]),
        np.array([1.0, 0.0]),
    )

    assert widths == pytest.approx(np.array([0.0, 2.0, 0.0]))


def test_integrate_flat_full_diamond_volume() -> None:
    p = np.array([-1.0, 0.0])
    q = np.array([1.0, 0.0])

    volume = integrate_conformal_interval_volume_1p1(p, q, flat_profile())

    assert volume == pytest.approx(2.0, rel=1e-4)


def test_integrate_constant_scale_volume() -> None:
    p = np.array([-1.0, 0.0])
    q = np.array([1.0, 0.0])

    volume = integrate_conformal_interval_volume_1p1(p, q, constant_profile(1.5))

    assert volume == pytest.approx(1.5**2 * 2.0, rel=1e-4)


def test_central_observer_proper_time_constant_scale() -> None:
    proper_time = central_observer_proper_time_1p1(
        -1.0,
        1.0,
        constant_profile(1.5),
    )

    assert proper_time == pytest.approx(3.0)


def test_weighted_interval_volume_matches_unweighted_for_flat_profile() -> None:
    events = np.array([[-1.0, 0.0], [0.0, 0.0], [1.0, 0.0]])
    causal_matrix = causal_matrix_1p1(events)
    coordinate_density = 10.0

    weighted = weighted_interval_volume_estimate_1p1(
        events,
        causal_matrix,
        0,
        2,
        flat_profile(),
        coordinate_density,
    )
    unweighted = unweighted_interval_coordinate_volume_estimate_1p1(
        causal_matrix,
        0,
        2,
        coordinate_density,
    )

    assert weighted == pytest.approx(unweighted)
