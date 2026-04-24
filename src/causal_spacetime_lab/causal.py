"""Causal-order utilities for 1+1D Minkowski events."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _as_event_array(events: ArrayLike) -> NDArray[np.float64]:
    array = np.asarray(events, dtype=float)
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError("events must have shape (n, 2), with columns (t, x)")
    return array


def causal_matrix_1p1(
    events: ArrayLike,
    *,
    atol: float = 1e-12,
) -> NDArray[np.bool_]:
    """Return the 1+1D causal precedence matrix for events ``(t, x)``.

    ``C[i, j]`` is true iff event ``i`` causally precedes event ``j``:

    ``dt > 0`` and ``dt^2 >= dx^2``.
    """

    event_array = _as_event_array(events)
    t = event_array[:, 0]
    x = event_array[:, 1]

    dt = t[None, :] - t[:, None]
    dx = x[None, :] - x[:, None]
    interval_squared = dt * dt - dx * dx
    return (dt > 0.0) & (interval_squared >= -atol)

