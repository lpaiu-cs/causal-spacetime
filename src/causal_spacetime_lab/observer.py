"""Observer-chain utilities for operational reconstruction protocols."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.lorentz import gamma, lorentz_transform_1p1


def make_stationary_observer_chain_1p1(
    T: float,
    num_ticks: int,
    x: float = 0.0,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return a supplied stationary observer chain ``O(tau) = (tau, x)``."""

    if T <= 0:
        raise ValueError("T must be positive")
    ticks = int(num_ticks)
    if ticks < 2:
        raise ValueError("num_ticks must be at least 2")

    clock_times = np.linspace(-T / 2.0, T / 2.0, ticks, dtype=np.float64)
    observer_events = np.column_stack(
        (clock_times, np.full(ticks, float(x), dtype=np.float64))
    )
    return observer_events, clock_times


def observer_chain_indices(start_index: int, num_ticks: int) -> NDArray[np.int_]:
    """Return contiguous indices for observer-chain events in a combined array."""

    start = int(start_index)
    ticks = int(num_ticks)
    if start < 0:
        raise ValueError("start_index must be non-negative")
    if ticks < 0:
        raise ValueError("num_ticks must be non-negative")
    return np.arange(start, start + ticks, dtype=int)


def validate_observer_chain_clock_order(clock_times: ArrayLike) -> NDArray[np.float64]:
    """Validate and return strictly increasing observer clock labels."""

    clocks = np.asarray(clock_times, dtype=float)
    if clocks.ndim != 1:
        raise ValueError("clock_times must be one-dimensional")
    if clocks.size < 2:
        raise ValueError("clock_times must contain at least two ticks")
    if not np.all(np.diff(clocks) > 0.0):
        raise ValueError("clock_times must be strictly increasing")
    return clocks.astype(np.float64, copy=False)


def make_inertial_observer_chain_1p1(
    beta: float,
    num_ticks: int,
    tau_min: float,
    tau_max: float,
    x_prime: float = 0.0,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return an inertial observer or beacon chain transformed to lab frame."""

    beta_value = float(beta)
    if abs(beta_value) >= 1.0:
        raise ValueError("beta must satisfy abs(beta) < 1 in units c = 1")
    ticks = int(num_ticks)
    if ticks < 2:
        raise ValueError("num_ticks must be at least 2")
    tau_min_value = float(tau_min)
    tau_max_value = float(tau_max)
    if tau_max_value <= tau_min_value:
        raise ValueError("tau_max must be greater than tau_min")

    clock_times = np.linspace(tau_min_value, tau_max_value, ticks, dtype=np.float64)
    x_prime_values = np.full(ticks, float(x_prime), dtype=np.float64)
    t_lab, x_lab = lorentz_transform_1p1(clock_times, x_prime_values, beta_value)
    observer_events = np.column_stack((t_lab, x_lab)).astype(np.float64, copy=False)
    return observer_events, clock_times


def inertial_chain_inside_diamond_mask(
    events: ArrayLike,
    T: float,
) -> NDArray[np.bool_]:
    """Return mask for 1+1D events inside the global causal diamond."""

    if T <= 0:
        raise ValueError("T must be positive")
    event_array = np.asarray(events, dtype=float)
    if event_array.ndim != 2 or event_array.shape[1] != 2:
        raise ValueError("events must have shape (n, 2), with columns (t, x)")
    return np.abs(event_array[:, 1]) <= T / 2.0 - np.abs(event_array[:, 0]) + 1e-12


def safe_tau_range_for_inertial_chain_1p1(
    T: float,
    beta: float,
    x_prime: float = 0.0,
    margin: float = 0.98,
) -> tuple[float, float]:
    """Return a conservative tau range whose chain lies inside the diamond."""

    if T <= 0:
        raise ValueError("T must be positive")
    beta_value = float(beta)
    if abs(beta_value) >= 1.0:
        raise ValueError("beta must satisfy abs(beta) < 1 in units c = 1")
    if not 0.0 < margin <= 1.0:
        raise ValueError("margin must satisfy 0 < margin <= 1")

    gamma_value = float(gamma(beta_value))
    center = -beta_value * float(x_prime)
    broad_radius = T / max(2.0 * gamma_value * (1.0 + abs(beta_value)), 1e-12)
    broad_radius += abs(float(x_prime)) + T
    candidates = np.linspace(center - broad_radius, center + broad_radius, 20_001)
    x_prime_values = np.full(candidates.size, float(x_prime), dtype=np.float64)
    t_lab, x_lab = lorentz_transform_1p1(candidates, x_prime_values, beta_value)
    events = np.column_stack((t_lab, x_lab))
    allowed = inertial_chain_inside_diamond_mask(events, T)
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

