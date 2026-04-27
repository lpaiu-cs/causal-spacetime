from __future__ import annotations

import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.sampling import sample_timelike_pairs


def test_sample_timelike_pairs_uses_causal_matrix() -> None:
    events = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [1.0, 2.0],
        ]
    )
    C = causal_matrix_1p1(events)

    pairs = sample_timelike_pairs(events, C, num_pairs=10, seed=3)

    assert pairs
    assert all(C[i, j] for i, j in pairs)
    assert (0, 3) not in pairs


def test_sample_timelike_pairs_filters_by_min_tau_when_requested() -> None:
    events = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
        ]
    )
    C = causal_matrix_1p1(events)

    pairs = sample_timelike_pairs(events, C, num_pairs=10, seed=1, min_tau=1.5)

    assert pairs == [(0, 2)]


def test_sample_timelike_pairs_filters_by_min_interval_count() -> None:
    events = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
        ]
    )
    C = causal_matrix_1p1(events)

    pairs = sample_timelike_pairs(
        events,
        C,
        num_pairs=10,
        seed=1,
        min_interval_count=1,
    )

    assert pairs == [(0, 2)]

