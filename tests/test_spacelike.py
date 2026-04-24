from __future__ import annotations

import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.spacelike import (
    alexandrov_interval_count_matrix,
    common_future_indices,
    common_future_overlap_count,
    common_past_indices,
    common_past_overlap_count,
    is_spacelike_pair,
    minimal_enclosing_alexandrov_interval,
    minimal_enclosing_alexandrov_interval_count,
)


def test_common_overlap_counts_for_hand_coded_spacelike_pair() -> None:
    events = np.array(
        [
            [-2.0, 0.0],
            [0.0, -1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [1.0, 0.0],
            [-1.0, 0.0],
        ]
    )
    C = causal_matrix_1p1(events)

    assert is_spacelike_pair(C, 1, 2)
    assert set(common_future_indices(C, 1, 2).tolist()) == {3, 4}
    assert set(common_past_indices(C, 1, 2).tolist()) == {0, 5}
    assert common_future_overlap_count(C, 1, 2) == 2
    assert common_past_overlap_count(C, 1, 2) == 2


def test_minimal_enclosing_alexandrov_interval_count() -> None:
    events = np.array(
        [
            [-2.0, 0.0],
            [0.0, -1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [1.0, 0.0],
            [-1.0, 0.0],
        ]
    )
    C = causal_matrix_1p1(events)
    interval_counts = alexandrov_interval_count_matrix(C)

    interval = minimal_enclosing_alexandrov_interval(C, 1, 2, interval_counts)

    assert interval is not None
    assert (interval.past_index, interval.future_index, interval.count) == (5, 4, 2)
    assert minimal_enclosing_alexandrov_interval_count(C, 1, 2) == 2


def test_minimal_enclosing_interval_returns_none_without_common_bounds() -> None:
    events = np.array(
        [
            [0.0, -1.0],
            [0.0, 1.0],
        ]
    )
    C = causal_matrix_1p1(events)

    assert minimal_enclosing_alexandrov_interval(C, 0, 1) is None

