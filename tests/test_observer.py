from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.observer import (
    inertial_chain_inside_diamond_mask,
    make_inertial_observer_chain_1p1,
    make_stationary_observer_chain_1p1,
    observer_chain_indices,
    safe_tau_range_for_inertial_chain_1p1,
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


def test_make_inertial_observer_chain_shape_and_clock_labels() -> None:
    observer_events, clock_times = make_inertial_observer_chain_1p1(
        beta=0.4,
        num_ticks=6,
        tau_min=-0.5,
        tau_max=0.5,
        x_prime=0.1,
    )

    assert observer_events.shape == (6, 2)
    assert clock_times.shape == (6,)
    assert clock_times[0] == pytest.approx(-0.5)
    assert clock_times[-1] == pytest.approx(0.5)
    assert np.all(np.diff(clock_times) > 0.0)


def test_make_inertial_observer_chain_rejects_invalid_beta() -> None:
    with pytest.raises(ValueError, match=r"abs\(beta\) < 1"):
        make_inertial_observer_chain_1p1(
            beta=1.0,
            num_ticks=4,
            tau_min=-0.5,
            tau_max=0.5,
        )


def test_safe_tau_range_for_inertial_chain_lies_inside_diamond() -> None:
    tau_min, tau_max = safe_tau_range_for_inertial_chain_1p1(
        T=2.0,
        beta=0.5,
        x_prime=0.12,
    )

    observer_events, _ = make_inertial_observer_chain_1p1(
        beta=0.5,
        num_ticks=32,
        tau_min=tau_min,
        tau_max=tau_max,
        x_prime=0.12,
    )

    assert np.all(inertial_chain_inside_diamond_mask(observer_events, T=2.0))
