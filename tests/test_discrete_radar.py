from __future__ import annotations

import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.discrete_radar import (
    discrete_radar_coordinates_from_order,
    find_radar_ticks_from_order,
)


def test_find_radar_ticks_from_order_hand_coded_example() -> None:
    events = np.array(
        [
            [-1.0, 0.0],
            [-0.5, 0.0],
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
            [0.0, 0.5],
        ]
    )
    C = causal_matrix_1p1(events)

    ticks = find_radar_ticks_from_order(
        C,
        observer_indices=np.array([0, 1, 2, 3, 4]),
        target_index=5,
        clock_times=np.array([-1.0, -0.5, 0.0, 0.5, 1.0]),
    )

    assert ticks == (-0.5, 0.5)


def test_find_radar_ticks_from_order_returns_none_when_inaccessible() -> None:
    events = np.array(
        [
            [-1.0, 0.0],
            [0.0, 0.0],
            [0.8, 0.0],
        ]
    )
    C = causal_matrix_1p1(events)

    ticks = find_radar_ticks_from_order(
        C,
        observer_indices=np.array([0, 1]),
        target_index=2,
        clock_times=np.array([-1.0, 0.0]),
    )

    assert ticks is None


def test_discrete_radar_coordinates_from_order_deterministic_example() -> None:
    events = np.array(
        [
            [-1.0, 0.0],
            [-0.5, 0.0],
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
            [0.0, 0.5],
        ]
    )
    C = causal_matrix_1p1(events)

    result = discrete_radar_coordinates_from_order(
        C,
        observer_indices=np.array([0, 1, 2, 3, 4]),
        target_indices=np.array([5]),
        clock_times=np.array([-1.0, -0.5, 0.0, 0.5, 1.0]),
    )[0]

    assert result.accessible
    assert result.radar_time == 0.0
    assert result.radar_distance == 0.5

