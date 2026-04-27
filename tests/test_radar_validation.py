from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.radar_validation import (
    accessible_fraction,
    radar_distance_error,
    radar_time_error,
    true_stationary_radar_coordinates_1p1,
)


def test_true_stationary_radar_coordinates_returns_t_and_abs_x() -> None:
    events = np.array(
        [
            [0.25, -0.5],
            [-0.75, 0.25],
        ]
    )

    time, distance = true_stationary_radar_coordinates_1p1(events)

    assert time.tolist() == pytest.approx([0.25, -0.75])
    assert distance.tolist() == pytest.approx([0.5, 0.25])


def test_radar_error_helpers() -> None:
    assert radar_time_error([1.0, 2.0], [0.5, 2.5]).tolist() == pytest.approx(
        [0.5, -0.5]
    )
    assert radar_distance_error([1.0, 2.0], [0.25, 3.0]).tolist() == pytest.approx(
        [0.75, -1.0]
    )
    assert accessible_fraction([True, False, True]) == pytest.approx(2 / 3)

