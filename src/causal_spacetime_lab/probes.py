"""Independent probe-pair utilities for 1+1D causal-set validation."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond


def _as_event_array(events: ArrayLike) -> NDArray[np.float64]:
    array = np.asarray(events, dtype=float)
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError("events must have shape (n, 2), with columns (t, x)")
    return array


def _as_event(event: ArrayLike, name: str) -> NDArray[np.float64]:
    array = np.asarray(event, dtype=float)
    if array.shape != (2,):
        raise ValueError(f"{name} must have shape (2,), representing (t, x)")
    return array


def _validate_tau_bounds(
    min_tau: float | None,
    max_tau: float | None,
) -> tuple[float | None, float | None]:
    lower = None if min_tau is None else float(min_tau)
    upper = None if max_tau is None else float(max_tau)
    if lower is not None and lower < 0.0:
        raise ValueError("min_tau must be non-negative")
    if upper is not None and upper <= 0.0:
        raise ValueError("max_tau must be positive")
    if lower is not None and upper is not None and lower > upper:
        raise ValueError("min_tau must be less than or equal to max_tau")
    return lower, upper


def sample_probe_timelike_pairs_1p1(
    num_pairs: int,
    T: float = 2.0,
    seed: int | None = None,
    min_tau: float | None = None,
    max_tau: float | None = None,
    max_attempts: int = 100_000,
) -> NDArray[np.float64]:
    """Sample independent ordered timelike probe endpoint pairs in a diamond."""

    requested = int(num_pairs)
    if requested < 0:
        raise ValueError("num_pairs must be non-negative")
    if T <= 0:
        raise ValueError("T must be positive")
    if max_attempts <= 0:
        raise ValueError("max_attempts must be positive")
    lower_tau, upper_tau = _validate_tau_bounds(min_tau, max_tau)
    if requested == 0:
        return np.empty((0, 2, 2), dtype=np.float64)

    rng = np.random.default_rng(seed)
    accepted: list[NDArray[np.float64]] = []
    attempts = 0
    while len(accepted) < requested and attempts < max_attempts:
        remaining = requested - len(accepted)
        batch_size = min(max(4 * remaining, 128), max_attempts - attempts)
        candidates = sprinkle_1p1_causal_diamond(2 * batch_size, T=T, seed=rng)
        first = candidates[:batch_size]
        second = candidates[batch_size:]
        attempts += batch_size

        swap = first[:, 0] > second[:, 0]
        p = first.copy()
        q = second.copy()
        p[swap] = second[swap]
        q[swap] = first[swap]

        dt = q[:, 0] - p[:, 0]
        dx = q[:, 1] - p[:, 1]
        tau_squared = dt * dt - dx * dx
        keep = (dt > 0.0) & (tau_squared >= 0.0)
        if lower_tau is not None:
            keep &= tau_squared >= lower_tau * lower_tau
        if upper_tau is not None:
            keep &= tau_squared <= upper_tau * upper_tau

        for pair in np.stack((p[keep], q[keep]), axis=1):
            accepted.append(pair.astype(np.float64, copy=False))
            if len(accepted) == requested:
                break

    if len(accepted) < requested:
        raise RuntimeError(
            "could not sample the requested number of probe timelike pairs "
            f"within {max_attempts} attempts"
        )

    return np.stack(accepted, axis=0)


def _precedes_many(
    source: NDArray[np.float64],
    targets: NDArray[np.float64],
) -> NDArray[np.bool_]:
    dt = targets[:, 0] - source[0]
    dx = targets[:, 1] - source[1]
    return (dt > 0.0) & (dt * dt >= dx * dx)


def count_interval_events_for_probe_pair_1p1(
    events: ArrayLike,
    p: ArrayLike,
    q: ArrayLike,
) -> int:
    """Count support events in the Alexandrov interval between probe endpoints."""

    event_array = _as_event_array(events)
    p_event = _as_event(p, "p")
    q_event = _as_event(q, "q")
    dt_from_p = event_array[:, 0] - p_event[0]
    dx_from_p = event_array[:, 1] - p_event[1]
    dt_to_q = q_event[0] - event_array[:, 0]
    dx_to_q = q_event[1] - event_array[:, 1]
    inside = (
        (dt_from_p > 0.0)
        & (dt_from_p * dt_from_p >= dx_from_p * dx_from_p)
        & (dt_to_q > 0.0)
        & (dt_to_q * dt_to_q >= dx_to_q * dx_to_q)
    )
    return int(np.count_nonzero(inside))


def count_interval_events_for_probe_pairs_1p1(
    events: ArrayLike,
    probe_pairs: ArrayLike,
) -> NDArray[np.int64]:
    """Count support events inside each independent probe Alexandrov interval."""

    event_array = _as_event_array(events)
    pairs = np.asarray(probe_pairs, dtype=float)
    if pairs.ndim != 3 or pairs.shape[1:] != (2, 2):
        raise ValueError("probe_pairs must have shape (num_pairs, 2, 2)")

    counts = np.empty(pairs.shape[0], dtype=np.int64)
    for index, (p, q) in enumerate(pairs):
        dt_from_p = event_array[:, 0] - p[0]
        dx_from_p = event_array[:, 1] - p[1]
        dt_to_q = q[0] - event_array[:, 0]
        dx_to_q = q[1] - event_array[:, 1]
        inside = (
            (dt_from_p > 0.0)
            & (dt_from_p * dt_from_p >= dx_from_p * dx_from_p)
            & (dt_to_q > 0.0)
            & (dt_to_q * dt_to_q >= dx_to_q * dx_to_q)
        )
        counts[index] = np.count_nonzero(inside)
    return counts
