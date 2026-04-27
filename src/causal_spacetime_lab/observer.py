"""Observer-chain utilities for operational reconstruction protocols."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


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

