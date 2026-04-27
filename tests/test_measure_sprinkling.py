from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.conformal import (
    conformal_volume_density_1p1,
    constant_profile,
    flat_profile,
    sinusoidal_time_profile,
)
from causal_spacetime_lab.measure_sprinkling import (
    conformal_full_diamond_volume_1p1,
    conformal_time_bin_masses_1p1,
    coordinate_time_bin_volumes_1p1,
    estimate_conformal_weight_max_1p1,
    estimate_time_bin_density_shape_1p1,
    normalize_profile_shape,
    sprinkle_conformal_measure_1p1,
)


def test_conformal_full_diamond_volume_flat() -> None:
    assert conformal_full_diamond_volume_1p1(2.0, flat_profile()) == pytest.approx(
        2.0,
        rel=1e-4,
    )


def test_conformal_full_diamond_volume_constant_scale() -> None:
    alpha = 1.5

    volume = conformal_full_diamond_volume_1p1(2.0, constant_profile(alpha))

    assert volume == pytest.approx(alpha * alpha * 2.0, rel=1e-4)


def test_coordinate_time_bin_volumes_sum_to_global_volume() -> None:
    edges = np.linspace(-1.0, 1.0, 9)

    volumes = coordinate_time_bin_volumes_1p1(2.0, edges)

    assert np.sum(volumes) == pytest.approx(2.0)


def test_conformal_time_bin_masses_sum_to_full_volume() -> None:
    edges = np.linspace(-1.0, 1.0, 9)
    profile = sinusoidal_time_profile(0.3, T=2.0)

    masses = conformal_time_bin_masses_1p1(2.0, profile, edges)
    full_volume = conformal_full_diamond_volume_1p1(2.0, profile)

    assert np.sum(masses) == pytest.approx(full_volume, rel=1e-3)


def test_estimate_conformal_weight_max_bounds_sampled_weights() -> None:
    T = 2.0
    profile = sinusoidal_time_profile(0.3, T=T)
    t_values = np.linspace(-T / 2.0, T / 2.0, 257)
    events = np.column_stack((t_values, np.zeros_like(t_values)))

    estimated_max = estimate_conformal_weight_max_1p1(T, profile)
    sampled_max = float(np.max(conformal_volume_density_1p1(profile, events)))

    assert estimated_max >= sampled_max


def test_sprinkle_conformal_measure_shape_and_diamond_condition() -> None:
    T = 2.0

    events = sprinkle_conformal_measure_1p1(
        50,
        T,
        sinusoidal_time_profile(0.2, T),
        seed=11,
    )

    assert events.shape == (50, 2)
    assert np.all(np.abs(events[:, 1]) <= T / 2.0 - np.abs(events[:, 0]) + 1e-12)


def test_sprinkle_conformal_measure_rejects_invalid_inputs() -> None:
    with pytest.raises(ValueError, match="n must be positive"):
        sprinkle_conformal_measure_1p1(0, 2.0, flat_profile())

    with pytest.raises(ValueError, match="T must be positive"):
        sprinkle_conformal_measure_1p1(10, 0.0, flat_profile())


def test_estimate_time_bin_density_shape_returns_finite_nonempty_bins() -> None:
    events = np.array(
        [
            [-0.75, 0.0],
            [-0.25, 0.0],
            [0.25, 0.0],
            [0.75, 0.0],
        ]
    )
    edges = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])

    densities = estimate_time_bin_density_shape_1p1(events, 2.0, edges)

    assert np.all(np.isfinite(densities))
    assert np.all(densities > 0.0)


def test_normalize_profile_shape_has_mean_one() -> None:
    values = np.array([2.0, 4.0, 6.0])

    normalized = normalize_profile_shape(values)

    assert np.nanmean(normalized) == pytest.approx(1.0)
    assert normalized.tolist() == pytest.approx([0.5, 1.0, 1.5])
