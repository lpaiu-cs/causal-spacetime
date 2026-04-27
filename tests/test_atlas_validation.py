from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.atlas_validation import (
    chart_interval_disagreement,
    minkowski_interval_squared_from_coords,
    sample_event_pairs,
)


def test_sample_event_pairs_returns_no_self_pairs() -> None:
    pairs = sample_event_pairs(np.array([0, 1, 2]), num_pairs=20, seed=3)

    assert pairs.shape == (20, 2)
    assert np.all(pairs[:, 0] != pairs[:, 1])


def test_minkowski_interval_squared_from_coords_hand_coded_pairs() -> None:
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 0.5]])
    pairs = np.array([[0, 1], [0, 2], [1, 2]])

    intervals = minkowski_interval_squared_from_coords(coords, pairs)

    assert intervals == pytest.approx(np.array([1.0, 0.75, -0.25]))


def test_chart_interval_disagreement_identical_charts_is_zero() -> None:
    coords = np.array(
        [
            [0.0, 0.0],
            [0.2, 0.1],
            [0.4, -0.1],
            [-0.3, 0.05],
        ]
    )
    accessible = np.array([True, True, True, True])

    summary = chart_interval_disagreement(
        coords,
        coords,
        accessible,
        num_pairs=30,
        seed=5,
    )

    assert summary["pair_count"] == 30.0
    assert summary["interval_rmse"] == pytest.approx(0.0)
    assert summary["interval_mae"] == pytest.approx(0.0)
    assert summary["interval_bias"] == pytest.approx(0.0)
