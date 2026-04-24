from __future__ import annotations

from math import sqrt

import numpy as np
import pytest

from causal_spacetime_lab.radar import (
    inertial_radar_coordinates_1p1,
    lorentz_boost_to_observer_rest_frame_1p1,
    radar_coordinates_1p1,
    stationary_radar_coordinates_1p1,
)


@pytest.mark.parametrize(
    "event",
    [
        (0.0, 0.0),
        (3.0, 1.25),
        (-2.0, -0.5),
    ],
)
def test_stationary_radar_coordinates_recover_t_and_abs_x(
    event: tuple[float, float],
) -> None:
    coordinates = stationary_radar_coordinates_1p1(event)

    assert coordinates.radar_time == pytest.approx(event[0])
    assert coordinates.radar_distance == pytest.approx(abs(event[1]))
    assert coordinates.tau_minus == pytest.approx(event[0] - abs(event[1]))
    assert coordinates.tau_plus == pytest.approx(event[0] + abs(event[1]))


def test_default_radar_coordinates_are_stationary_coordinates() -> None:
    event = np.array([4.0, -1.5])

    coordinates = radar_coordinates_1p1(event)

    assert coordinates.radar_time == pytest.approx(4.0)
    assert coordinates.radar_distance == pytest.approx(1.5)


def test_moving_inertial_observer_uses_lorentz_rest_frame_coordinates() -> None:
    observer_velocity = 0.6
    gamma = 1.0 / sqrt(1.0 - observer_velocity**2)
    rest_frame_event = np.array([2.0, 0.5])
    lab_frame_event = np.array(
        [
            gamma
            * (rest_frame_event[0] + observer_velocity * rest_frame_event[1]),
            gamma
            * (rest_frame_event[1] + observer_velocity * rest_frame_event[0]),
        ]
    )

    boosted_event = lorentz_boost_to_observer_rest_frame_1p1(
        lab_frame_event,
        observer_velocity,
    )
    coordinates = inertial_radar_coordinates_1p1(
        lab_frame_event,
        observer_velocity=observer_velocity,
    )

    assert boosted_event == pytest.approx(rest_frame_event)
    assert coordinates.radar_time == pytest.approx(rest_frame_event[0])
    assert coordinates.radar_distance == pytest.approx(abs(rest_frame_event[1]))
    assert coordinates.tau_minus == pytest.approx(1.5)
    assert coordinates.tau_plus == pytest.approx(2.5)


def test_radar_coordinates_reject_superluminal_observer_velocity() -> None:
    with pytest.raises(ValueError, match=r"abs\(v\) < 1"):
        inertial_radar_coordinates_1p1((1.0, 0.0), observer_velocity=1.0)

