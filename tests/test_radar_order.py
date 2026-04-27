from __future__ import annotations

import numpy as np

from causal_spacetime_lab.radar_order import (
    earliest_successor_tick_index,
    latest_predecessor_tick_index,
    radar_return_order_from_successor_ticks,
    radar_tick_brackets_from_order,
)


def _causal_matrix() -> np.ndarray:
    matrix = np.zeros((7, 7), dtype=bool)
    matrix[0, 4] = True
    matrix[2, 4] = True
    matrix[4, 3] = True
    matrix[1, 5] = True
    matrix[5, 3] = True
    matrix[6, 3] = True
    return matrix


def test_latest_predecessor_tick_index_deterministic() -> None:
    observers = np.asarray([0, 1, 2, 3], dtype=int)

    assert latest_predecessor_tick_index(_causal_matrix(), observers, 4) == 2


def test_earliest_successor_tick_index_deterministic() -> None:
    observers = np.asarray([0, 1, 2, 3], dtype=int)

    assert earliest_successor_tick_index(_causal_matrix(), observers, 4) == 3


def test_radar_tick_brackets_from_order_deterministic() -> None:
    observers = np.asarray([0, 1, 2, 3], dtype=int)
    targets = np.asarray([4, 5, 6], dtype=int)

    predecessor, successor, accessible = radar_tick_brackets_from_order(
        _causal_matrix(),
        observers,
        targets,
    )

    assert predecessor.tolist() == [2, 1, -1]
    assert successor.tolist() == [3, 3, 3]
    assert accessible.tolist() == [True, True, False]


def test_radar_return_order_from_successor_ticks_deterministic() -> None:
    successors = np.asarray([2, 4, 3], dtype=int)
    accessible = np.asarray([True, True, False])

    order = radar_return_order_from_successor_ticks(successors, accessible)

    assert order.tolist() == [
        [False, True, False],
        [False, False, False],
        [False, False, False],
    ]
