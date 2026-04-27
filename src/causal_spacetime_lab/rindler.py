"""Rindler observer utilities for flat 1+1D validation experiments."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.observer import inertial_chain_inside_diamond_mask


def _validate_acceleration(acceleration: float) -> float:
    value = float(acceleration)
    if value <= 0.0:
        raise ValueError("acceleration must be positive")
    return value


def _validate_direction(direction: int) -> int:
    value = int(direction)
    if value not in {-1, 1}:
        raise ValueError("direction must be -1 or 1")
    return value


def _validate_origin(origin: tuple[float, float]) -> tuple[float, float]:
    if len(origin) != 2:
        raise ValueError("horizon_origin must contain two values")
    return float(origin[0]), float(origin[1])


def _as_event_array(events: ArrayLike) -> NDArray[np.float64]:
    array = np.asarray(events, dtype=float)
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError("events must have shape (n, 2)")
    return array


def rindler_observer_event_1p1(
    proper_time: float,
    acceleration: float,
    horizon_origin: tuple[float, float] = (0.0, 0.0),
    direction: int = 1,
) -> NDArray[np.float64]:
    """Return one uniformly accelerated observer event in 1+1D Minkowski."""

    accel = _validate_acceleration(acceleration)
    wedge_direction = _validate_direction(direction)
    origin_t, origin_x = _validate_origin(horizon_origin)
    tau = float(proper_time)
    xi = 1.0 / accel
    return np.asarray(
        [
            origin_t + xi * np.sinh(accel * tau),
            origin_x + wedge_direction * xi * np.cosh(accel * tau),
        ],
        dtype=np.float64,
    )


def make_rindler_observer_chain_1p1(
    acceleration: float,
    num_ticks: int,
    tau_min: float,
    tau_max: float,
    horizon_origin: tuple[float, float] = (0.0, 0.0),
    direction: int = 1,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return a supplied Rindler observer chain and proper-time labels."""

    _validate_acceleration(acceleration)
    _validate_direction(direction)
    ticks = int(num_ticks)
    if ticks < 2:
        raise ValueError("num_ticks must be at least 2")
    tau_min_value = float(tau_min)
    tau_max_value = float(tau_max)
    if tau_max_value <= tau_min_value:
        raise ValueError("tau_max must be greater than tau_min")

    clock_times = np.linspace(tau_min_value, tau_max_value, ticks, dtype=np.float64)
    events = np.vstack(
        [
            rindler_observer_event_1p1(
                tau,
                acceleration,
                horizon_origin=horizon_origin,
                direction=direction,
            )
            for tau in clock_times
        ]
    )
    return events.astype(np.float64, copy=False), clock_times


def right_rindler_wedge_mask(
    events: ArrayLike,
    horizon_origin: tuple[float, float] = (0.0, 0.0),
) -> NDArray[np.bool_]:
    """Return mask for the right Rindler wedge ``x - x0 > |t - t0|``."""

    event_array = _as_event_array(events)
    origin_t, origin_x = _validate_origin(horizon_origin)
    t_rel = event_array[:, 0] - origin_t
    x_rel = event_array[:, 1] - origin_x
    return x_rel > np.abs(t_rel)


def left_rindler_wedge_mask(
    events: ArrayLike,
    horizon_origin: tuple[float, float] = (0.0, 0.0),
) -> NDArray[np.bool_]:
    """Return mask for the left Rindler wedge ``x - x0 < -|t - t0|``."""

    event_array = _as_event_array(events)
    origin_t, origin_x = _validate_origin(horizon_origin)
    t_rel = event_array[:, 0] - origin_t
    x_rel = event_array[:, 1] - origin_x
    return x_rel < -np.abs(t_rel)


def rindler_wedge_mask(
    events: ArrayLike,
    horizon_origin: tuple[float, float] = (0.0, 0.0),
    direction: int = 1,
) -> NDArray[np.bool_]:
    """Return the Rindler wedge mask for ``direction`` +1 or -1."""

    wedge_direction = _validate_direction(direction)
    if wedge_direction == 1:
        return right_rindler_wedge_mask(events, horizon_origin)
    return left_rindler_wedge_mask(events, horizon_origin)


def analytic_rindler_radar_ticks_1p1(
    events: ArrayLike,
    acceleration: float,
    horizon_origin: tuple[float, float] = (0.0, 0.0),
    direction: int = 1,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.bool_]]:
    """Return analytic Rindler radar ticks for validation only.

    For ``direction=-1`` the spatial coordinate is reflected around the horizon
    origin and the right-wedge formula is applied in that reflected coordinate.
    """

    event_array = _as_event_array(events)
    accel = _validate_acceleration(acceleration)
    wedge_direction = _validate_direction(direction)
    origin_t, origin_x = _validate_origin(horizon_origin)

    t_rel = event_array[:, 0] - origin_t
    directed_x = wedge_direction * (event_array[:, 1] - origin_x)
    u = t_rel - directed_x
    v = t_rel + directed_x
    accessible = (u < 0.0) & (v > 0.0)

    tau_minus = np.full(event_array.shape[0], np.nan, dtype=np.float64)
    tau_plus = np.full(event_array.shape[0], np.nan, dtype=np.float64)
    if np.any(accessible):
        tau_u = -(1.0 / accel) * np.log(-accel * u[accessible])
        tau_v = (1.0 / accel) * np.log(accel * v[accessible])
        tau_minus[accessible] = np.minimum(tau_u, tau_v)
        tau_plus[accessible] = np.maximum(tau_u, tau_v)

    return tau_minus, tau_plus, accessible


def analytic_rindler_radar_coordinates_1p1(
    events: ArrayLike,
    acceleration: float,
    horizon_origin: tuple[float, float] = (0.0, 0.0),
    direction: int = 1,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.bool_]]:
    """Return analytic Rindler radar coordinates for validation only."""

    tau_minus, tau_plus, accessible = analytic_rindler_radar_ticks_1p1(
        events,
        acceleration,
        horizon_origin=horizon_origin,
        direction=direction,
    )
    radar_time = np.full_like(tau_minus, np.nan)
    radar_distance = np.full_like(tau_plus, np.nan)
    radar_time[accessible] = 0.5 * (tau_plus[accessible] + tau_minus[accessible])
    radar_distance[accessible] = 0.5 * (
        tau_plus[accessible] - tau_minus[accessible]
    )
    return radar_time, radar_distance, accessible


def safe_tau_range_for_rindler_chain_1p1(
    T_global: float,
    acceleration: float,
    horizon_origin: tuple[float, float] = (0.0, 0.0),
    direction: int = 1,
    margin: float = 0.98,
    samples: int = 20_001,
) -> tuple[float, float]:
    """Return a conservative Rindler proper-time range inside the diamond."""

    if T_global <= 0.0:
        raise ValueError("T_global must be positive")
    accel = _validate_acceleration(acceleration)
    _validate_direction(direction)
    if not 0.0 < margin <= 1.0:
        raise ValueError("margin must satisfy 0 < margin <= 1")
    sample_count = int(samples)
    if sample_count < 3:
        raise ValueError("samples must be at least 3")

    origin_t, origin_x = _validate_origin(horizon_origin)
    broad_radius = (
        np.log1p(accel * (T_global + abs(origin_t) + abs(origin_x) + 1.0)) / accel
    )
    broad_radius = max(broad_radius + 1.0 / accel, T_global / 2.0)
    candidates = np.linspace(-broad_radius, broad_radius, sample_count)
    events = np.vstack(
        [
            rindler_observer_event_1p1(
                tau,
                accel,
                horizon_origin=(origin_t, origin_x),
                direction=direction,
            )
            for tau in candidates
        ]
    )
    allowed = inertial_chain_inside_diamond_mask(events, T_global)
    if not np.any(allowed):
        raise ValueError("no safe tau range exists for this chain in the diamond")

    allowed_tau = candidates[allowed]
    tau_min = float(np.min(allowed_tau))
    tau_max = float(np.max(allowed_tau))
    midpoint = 0.5 * (tau_min + tau_max)
    half_width = 0.5 * (tau_max - tau_min) * float(margin)
    if half_width <= 0.0:
        raise ValueError("safe tau range has zero width")
    return midpoint - half_width, midpoint + half_width
