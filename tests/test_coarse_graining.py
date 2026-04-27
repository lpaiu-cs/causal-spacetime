from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.coarse_graining import (
    rescaled_density_after_thinning,
    thin_events,
    thinning_summary_stats,
)


def test_thin_events_keep_probability_one_returns_all_events() -> None:
    events = np.array([[0.0, 0.0], [0.2, 0.1]])

    thinned = thin_events(events, keep_probability=1.0, seed=3)

    assert thinned == pytest.approx(events)


def test_thin_events_keep_probability_zero_preserves_column_count() -> None:
    events = np.array([[0.0, 0.0], [0.2, 0.1]])

    thinned = thin_events(events, keep_probability=0.0, seed=3)

    assert thinned.shape == (0, 2)


def test_thin_events_rejects_invalid_probabilities() -> None:
    events = np.array([[0.0, 0.0]])

    with pytest.raises(ValueError, match="keep_probability"):
        thin_events(events, keep_probability=-0.1)

    with pytest.raises(ValueError, match="keep_probability"):
        thin_events(events, keep_probability=1.1)


def test_rescaled_density_after_thinning() -> None:
    assert rescaled_density_after_thinning(50.0, 0.25) == pytest.approx(12.5)


def test_thinning_summary_stats_deterministic_counts() -> None:
    stats = thinning_summary_stats(
        original_count=100,
        thinned_count=37,
        keep_probability=0.4,
    )

    assert stats["original_count"] == 100.0
    assert stats["thinned_count"] == 37.0
    assert stats["expected_count"] == pytest.approx(40.0)
    assert stats["count_ratio"] == pytest.approx(0.37)
    assert stats["expected_ratio"] == pytest.approx(0.4)
    assert stats["count_ratio_error"] == pytest.approx(-0.03)
