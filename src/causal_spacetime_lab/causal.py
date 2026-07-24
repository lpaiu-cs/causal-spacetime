"""Causal-order utilities for 1+1D Minkowski events."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

#: Interval tolerance of the causal predicate: an event pair counts as
#: causal when ``dt > 0`` and ``dt^2 - ||dx||^2 >= -DEFAULT_CAUSAL_ATOL``.
#: Named because the null-inclusive convention is load-bearing (T1
#: Section 2) and code that evaluates the same relation without going
#: through this module -- vectorized fast paths, interval-membership
#: tests -- must use the SAME tolerance or it silently disagrees on
#: near-null pairs.
DEFAULT_CAUSAL_ATOL = 1e-12


def _as_event_array(events: ArrayLike) -> NDArray[np.float64]:
    array = np.asarray(events, dtype=float)
    if array.ndim != 2 or array.shape[1] < 2:
        raise ValueError("events must have shape (n, D), D >= 2")
    return array


def causal_matrix_minkowski(
    events: ArrayLike,
    *,
    atol: float = DEFAULT_CAUSAL_ATOL,
) -> NDArray[np.bool_]:
    """Return the Minkowski causal precedence matrix for events ``(t, x...)``.

    ``C[i, j]`` is true iff event ``i`` causally precedes event ``j``:

    ``dt > 0`` and ``dt^2 >= ||dx_vec||^2``.
    """

    event_array = _as_event_array(events)
    t = event_array[:, 0]

    dt = t[None, :] - t[:, None]
    spatial_distance_squared = np.zeros_like(dt)
    for axis in range(1, event_array.shape[1]):
        dx = event_array[None, :, axis] - event_array[:, None, axis]
        spatial_distance_squared += dx * dx
    interval_squared = dt * dt - spatial_distance_squared
    return (dt > 0.0) & (interval_squared >= -atol)


def causal_matrix_1p1(
    events: ArrayLike,
    *,
    atol: float = DEFAULT_CAUSAL_ATOL,
) -> NDArray[np.bool_]:
    """Return the 1+1D causal precedence matrix for events ``(t, x)``."""

    event_array = np.asarray(events, dtype=float)
    if event_array.ndim != 2 or event_array.shape[1] != 2:
        raise ValueError("events must have shape (n, 2), with columns (t, x)")
    return causal_matrix_minkowski(event_array, atol=atol)
