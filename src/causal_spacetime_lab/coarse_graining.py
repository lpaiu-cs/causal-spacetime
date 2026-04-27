"""Coarse-graining and thinning utilities for finite causal-set experiments."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _as_event_array(events: ArrayLike) -> NDArray[np.float64]:
    array = np.asarray(events, dtype=float)
    if array.ndim != 2:
        raise ValueError("events must have shape (n, D)")
    return array


def thin_events(
    events: ArrayLike,
    keep_probability: float,
    seed: int | None = None,
) -> NDArray[np.float64]:
    """Randomly thin an event set with independent keep probability."""

    event_array = _as_event_array(events)
    probability = float(keep_probability)
    if not 0.0 <= probability <= 1.0:
        raise ValueError("keep_probability must satisfy 0 <= p <= 1")
    if probability == 1.0:
        return event_array.astype(np.float64, copy=True)
    if probability == 0.0:
        return np.empty((0, event_array.shape[1]), dtype=np.float64)
    rng = np.random.default_rng(seed)
    keep = rng.uniform(0.0, 1.0, event_array.shape[0]) < probability
    return event_array[keep].astype(np.float64, copy=True)


def rescaled_density_after_thinning(
    original_density: float,
    keep_probability: float,
) -> float:
    """Return expected density after independent random thinning."""

    density = float(original_density)
    probability = float(keep_probability)
    if density <= 0.0:
        raise ValueError("original_density must be positive")
    if not 0.0 <= probability <= 1.0:
        raise ValueError("keep_probability must satisfy 0 <= p <= 1")
    return probability * density


def thinning_summary_stats(
    original_count: int,
    thinned_count: int,
    keep_probability: float,
) -> dict[str, float]:
    """Return simple realized-versus-expected thinning statistics."""

    original = int(original_count)
    thinned = int(thinned_count)
    probability = float(keep_probability)
    if original < 0:
        raise ValueError("original_count must be non-negative")
    if thinned < 0:
        raise ValueError("thinned_count must be non-negative")
    if thinned > original:
        raise ValueError("thinned_count cannot exceed original_count")
    if not 0.0 <= probability <= 1.0:
        raise ValueError("keep_probability must satisfy 0 <= p <= 1")
    expected_count = probability * original
    count_ratio = thinned / original if original > 0 else float("nan")
    return {
        "original_count": float(original),
        "thinned_count": float(thinned),
        "expected_count": float(expected_count),
        "count_ratio": float(count_ratio),
        "expected_ratio": probability,
        "count_ratio_error": float(count_ratio - probability)
        if np.isfinite(count_ratio)
        else float("nan"),
    }
