from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.rindler import (
    analytic_rindler_radar_coordinates_1p1,
    analytic_rindler_radar_ticks_1p1,
    left_rindler_wedge_mask,
    make_rindler_observer_chain_1p1,
    right_rindler_wedge_mask,
    rindler_observer_event_1p1,
    safe_tau_range_for_rindler_chain_1p1,
)


def test_rindler_observer_event_lies_on_expected_hyperbola() -> None:
    acceleration = 2.0
    event = rindler_observer_event_1p1(0.3, acceleration)

    interval = event[1] ** 2 - event[0] ** 2

    assert interval == pytest.approx(1.0 / acceleration**2)


def test_make_rindler_observer_chain_shape_and_clocks() -> None:
    events, clocks = make_rindler_observer_chain_1p1(
        acceleration=1.5,
        num_ticks=5,
        tau_min=-0.4,
        tau_max=0.4,
    )

    assert events.shape == (5, 2)
    assert clocks.shape == (5,)
    assert np.all(np.diff(clocks) > 0.0)


def test_right_rindler_wedge_mask_hand_coded_events() -> None:
    events = np.array([[0.0, 1.0], [0.5, 0.6], [0.5, 0.4], [0.0, -1.0]])

    mask = right_rindler_wedge_mask(events)

    assert mask.tolist() == [True, True, False, False]


def test_left_rindler_wedge_mask_hand_coded_events() -> None:
    events = np.array([[0.0, -1.0], [0.5, -0.6], [0.5, -0.4], [0.0, 1.0]])

    mask = left_rindler_wedge_mask(events)

    assert mask.tolist() == [True, True, False, False]


def test_analytic_rindler_radar_ticks_hand_coded_right_wedge() -> None:
    events = np.array([[0.0, 1.0], [0.0, 2.0]])

    tau_minus, tau_plus, accessible = analytic_rindler_radar_ticks_1p1(
        events,
        acceleration=1.0,
    )

    assert accessible.tolist() == [True, True]
    assert tau_minus[0] == pytest.approx(0.0)
    assert tau_plus[0] == pytest.approx(0.0)
    assert tau_minus[1] == pytest.approx(-np.log(2.0))
    assert tau_plus[1] == pytest.approx(np.log(2.0))


def test_analytic_rindler_radar_coordinates_nan_outside_wedge() -> None:
    events = np.array([[0.0, 1.0], [0.0, -1.0]])

    radar_time, radar_distance, accessible = analytic_rindler_radar_coordinates_1p1(
        events,
        acceleration=1.0,
    )

    assert accessible.tolist() == [True, False]
    assert np.isfinite(radar_time[0])
    assert np.isfinite(radar_distance[0])
    assert np.isnan(radar_time[1])
    assert np.isnan(radar_distance[1])


def test_safe_tau_range_for_rindler_chain_lies_inside_global_diamond() -> None:
    tau_min, tau_max = safe_tau_range_for_rindler_chain_1p1(
        T_global=4.0,
        acceleration=2.0,
    )
    events, _ = make_rindler_observer_chain_1p1(
        acceleration=2.0,
        num_ticks=32,
        tau_min=tau_min,
        tau_max=tau_max,
    )

    assert np.all(np.abs(events[:, 1]) <= 2.0 - np.abs(events[:, 0]) + 1e-12)
