"""Radar coordinate decomposition in 1+1D special relativity."""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class RadarCoordinates:
    """Emission/reception times and the resulting radar decomposition."""

    tau_minus: float
    tau_plus: float
    radar_time: float
    radar_distance: float


def _as_event(event: ArrayLike) -> NDArray[np.float64]:
    array = np.asarray(event, dtype=float)
    if array.shape != (2,):
        raise ValueError("event must have shape (2,), representing (t, x)")
    return array


def _validate_velocity(velocity: float) -> float:
    v = float(velocity)
    if abs(v) >= 1.0:
        raise ValueError("observer_velocity must satisfy abs(v) < 1 in units c = 1")
    return v


def lorentz_boost_to_observer_rest_frame_1p1(
    event: ArrayLike,
    observer_velocity: float,
) -> NDArray[np.float64]:
    """Transform a lab-frame event into an inertial observer's rest frame.

    The observer moves with constant velocity ``observer_velocity`` along the
    lab-frame x-axis. With natural units ``c = 1``, the passive boost is:

    ``t' = gamma * (t - v x)``, ``x' = gamma * (x - v t)``.
    """

    t, x = _as_event(event)
    v = _validate_velocity(observer_velocity)
    gamma = 1.0 / sqrt(1.0 - v * v)
    return np.array([gamma * (t - v * x), gamma * (x - v * t)])


def stationary_radar_coordinates_1p1(event: ArrayLike) -> RadarCoordinates:
    """Return radar coordinates for the stationary worldline ``O(tau)=(tau, 0)``."""

    t, x = _as_event(event)
    spatial_separation = abs(float(x))
    tau_minus = float(t) - spatial_separation
    tau_plus = float(t) + spatial_separation
    radar_time = 0.5 * (tau_plus + tau_minus)
    radar_distance = 0.5 * (tau_plus - tau_minus)
    return RadarCoordinates(
        tau_minus=tau_minus,
        tau_plus=tau_plus,
        radar_time=radar_time,
        radar_distance=radar_distance,
    )


def inertial_radar_coordinates_1p1(
    event: ArrayLike,
    observer_velocity: float = 0.0,
) -> RadarCoordinates:
    """Return radar coordinates for a moving inertial observer.

    This is an operational coordinate protocol: first transform the event into
    the inertial observer's rest frame, then apply the stationary radar
    decomposition there.
    """

    rest_frame_event = lorentz_boost_to_observer_rest_frame_1p1(
        event,
        observer_velocity,
    )
    return stationary_radar_coordinates_1p1(rest_frame_event)


def radar_coordinates_1p1(
    event: ArrayLike,
    observer_velocity: float = 0.0,
) -> RadarCoordinates:
    """Alias for inertial radar coordinates in 1+1D Minkowski spacetime."""

    return inertial_radar_coordinates_1p1(event, observer_velocity)

