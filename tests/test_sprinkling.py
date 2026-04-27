from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.sprinkling import (
    sprinkle_1p1_causal_diamond,
    sprinkle_1p1_forward_cone,
    sprinkle_minkowski_causal_diamond,
)


def test_sprinkled_events_satisfy_diamond_condition() -> None:
    T = 2.5
    events = sprinkle_1p1_causal_diamond(1000, T=T, seed=7)

    assert events.shape == (1000, 2)
    t = events[:, 0]
    x = events[:, 1]
    assert np.all(np.abs(x) <= T / 2.0 - np.abs(t) + 1e-12)


def test_sprinkle_rejects_invalid_inputs() -> None:
    with pytest.raises(ValueError, match="n must be non-negative"):
        sprinkle_1p1_causal_diamond(-1)

    with pytest.raises(ValueError, match="T must be positive"):
        sprinkle_1p1_causal_diamond(1, T=0.0)


def test_sprinkled_events_satisfy_forward_cone_condition() -> None:
    T = 5.0
    events = sprinkle_1p1_forward_cone(1000, T=T, seed=19)

    assert events.shape == (1000, 2)
    t = events[:, 0]
    x = events[:, 1]
    assert np.all(t >= 0.0)
    assert np.all(t <= T)
    assert np.all(np.abs(x) <= t + 1e-12)


def test_generalized_sprinkling_shape_and_diamond_condition() -> None:
    T = 2.0
    events = sprinkle_minkowski_causal_diamond(
        500,
        spacetime_dim=4,
        T=T,
        seed=23,
    )

    assert events.shape == (500, 4)
    t = events[:, 0]
    spatial_norm = np.linalg.norm(events[:, 1:], axis=1)
    assert np.all(spatial_norm <= T / 2.0 - np.abs(t) + 1e-12)
