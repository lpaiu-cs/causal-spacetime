from __future__ import annotations

import numpy as np

from causal_spacetime_lab.probes import (
    count_interval_events_for_probe_pair_1p1,
    count_interval_events_for_probe_pairs_1p1,
    sample_probe_timelike_pairs_1p1,
)


def test_probe_pair_sampler_returns_pairs_inside_global_diamond() -> None:
    T = 2.0
    pairs = sample_probe_timelike_pairs_1p1(100, T=T, seed=5, min_tau=0.05)

    assert pairs.shape == (100, 2, 2)
    t = pairs[:, :, 0]
    x = pairs[:, :, 1]
    assert np.all(np.abs(x) <= T / 2.0 - np.abs(t) + 1e-12)


def test_probe_pair_sampler_returns_ordered_timelike_pairs() -> None:
    pairs = sample_probe_timelike_pairs_1p1(100, T=2.0, seed=7, min_tau=0.05)

    dt = pairs[:, 1, 0] - pairs[:, 0, 0]
    dx = pairs[:, 1, 1] - pairs[:, 0, 1]
    assert np.all(dt > 0.0)
    assert np.all(dt * dt >= dx * dx)


def test_count_interval_events_for_probe_pair_hand_coded_example() -> None:
    events = np.array(
        [
            [-0.5, 0.0],
            [0.0, 0.5],
            [0.9, 0.0],
            [0.0, 1.1],
            [-1.0, 0.0],
            [1.0, 0.0],
        ]
    )
    p = np.array([-1.0, 0.0])
    q = np.array([1.0, 0.0])

    count = count_interval_events_for_probe_pair_1p1(events, p, q)
    vectorized_counts = count_interval_events_for_probe_pairs_1p1(
        events,
        np.array([[p, q]]),
    )

    assert count == 3
    assert vectorized_counts.tolist() == [3]

