from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.observer import (
    make_stationary_observer_chain_1p1,
    observer_chain_indices,
    validate_observer_chain_clock_order,
)


def test_make_stationary_observer_chain_shape_and_default_x() -> None:
    observer_events, clock_times = make_stationary_observer_chain_1p1(
        T=2.0,
        num_ticks=5,
    )

    assert observer_events.shape == (5, 2)
    assert clock_times.shape == (5,)
    assert np.all(observer_events[:, 1] == 0.0)
    assert observer_events[:, 0].tolist() == pytest.approx(clock_times)


def test_observer_clock_labels_are_strictly_increasing() -> None:
    _, clock_times = make_stationary_observer_chain_1p1(T=2.0, num_ticks=5)

    validated = validate_observer_chain_clock_order(clock_times)

    assert np.all(np.diff(validated) > 0.0)


def test_validate_observer_clock_order_rejects_repeated_label() -> None:
    with pytest.raises(ValueError, match="strictly increasing"):
        validate_observer_chain_clock_order([0.0, 0.0, 1.0])


def test_observer_chain_indices() -> None:
    assert observer_chain_indices(3, 4).tolist() == [3, 4, 5, 6]

