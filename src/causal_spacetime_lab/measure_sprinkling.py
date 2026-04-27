"""Measure-encoded sprinkling utilities for 1+1D conformal toy models."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.conformal import (
    ConformalProfile,
    conformal_volume_density_1p1,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond


def _validate_t(T: float) -> float:
    value = float(T)
    if value <= 0.0:
        raise ValueError("T must be positive")
    return value


def _diamond_width(t_values: NDArray[np.float64], T: float) -> NDArray[np.float64]:
    return np.maximum(0.0, 2.0 * (T / 2.0 - np.abs(t_values)))


def conformal_full_diamond_volume_1p1(
    T: float,
    profile: ConformalProfile,
    num_t: int = 4096,
) -> float:
    """Integrate full-diamond conformal volume for time-only profiles."""

    duration = _validate_t(T)
    sample_count = int(num_t)
    if sample_count < 2:
        raise ValueError("num_t must be at least 2")
    t_values = np.linspace(-duration / 2.0, duration / 2.0, sample_count)
    events = np.column_stack((t_values, np.zeros(sample_count, dtype=float)))
    weights = conformal_volume_density_1p1(profile, events)
    return float(np.trapezoid(_diamond_width(t_values, duration) * weights, t_values))


def _coordinate_antiderivative(t: NDArray[np.float64], T: float) -> NDArray[np.float64]:
    return np.where(t < 0.0, T * t + t * t, T * t - t * t)


def coordinate_time_bin_volumes_1p1(
    T: float,
    bin_edges: ArrayLike,
) -> NDArray[np.float64]:
    """Return coordinate volume of each time bin in the 1+1D diamond."""

    duration = _validate_t(T)
    edges = np.asarray(bin_edges, dtype=float)
    if edges.ndim != 1 or edges.size < 2:
        raise ValueError("bin_edges must be one-dimensional with at least two values")
    if not np.all(np.diff(edges) > 0.0):
        raise ValueError("bin_edges must be strictly increasing")

    lower = np.clip(edges[:-1], -duration / 2.0, duration / 2.0)
    upper = np.clip(edges[1:], -duration / 2.0, duration / 2.0)
    volumes = _coordinate_antiderivative(upper, duration) - _coordinate_antiderivative(
        lower,
        duration,
    )
    volumes[upper <= lower] = 0.0
    return volumes.astype(np.float64, copy=False)


def conformal_time_bin_masses_1p1(
    T: float,
    profile: ConformalProfile,
    bin_edges: ArrayLike,
    num_t_per_bin: int = 512,
) -> NDArray[np.float64]:
    """Return conformal physical mass of each time bin for time-only profiles."""

    duration = _validate_t(T)
    edges = np.asarray(bin_edges, dtype=float)
    if edges.ndim != 1 or edges.size < 2:
        raise ValueError("bin_edges must be one-dimensional with at least two values")
    if not np.all(np.diff(edges) > 0.0):
        raise ValueError("bin_edges must be strictly increasing")
    sample_count = int(num_t_per_bin)
    if sample_count < 2:
        raise ValueError("num_t_per_bin must be at least 2")

    masses = np.zeros(edges.size - 1, dtype=np.float64)
    for index, (raw_lower, raw_upper) in enumerate(
        zip(edges[:-1], edges[1:], strict=True)
    ):
        lower = max(float(raw_lower), -duration / 2.0)
        upper = min(float(raw_upper), duration / 2.0)
        if upper <= lower:
            continue
        t_values = np.linspace(lower, upper, sample_count)
        events = np.column_stack((t_values, np.zeros(sample_count, dtype=float)))
        weights = conformal_volume_density_1p1(profile, events)
        masses[index] = np.trapezoid(
            _diamond_width(t_values, duration) * weights,
            t_values,
        )
    return masses


def estimate_conformal_weight_max_1p1(
    T: float,
    profile: ConformalProfile,
    samples: int = 4096,
    safety_factor: float = 1.01,
) -> float:
    """Return a conservative sampled maximum of ``Omega^2`` inside the diamond."""

    duration = _validate_t(T)
    sample_count = int(samples)
    if sample_count < 2:
        raise ValueError("samples must be at least 2")
    factor = float(safety_factor)
    if factor < 1.0:
        raise ValueError("safety_factor must be at least 1")
    t_values = np.linspace(-duration / 2.0, duration / 2.0, sample_count)
    events = np.column_stack((t_values, np.zeros(sample_count, dtype=float)))
    return float(np.max(conformal_volume_density_1p1(profile, events)) * factor)


def sprinkle_conformal_measure_1p1(
    n: int,
    T: float,
    profile: ConformalProfile,
    seed: int | np.random.Generator | None = None,
    batch_size: int = 4096,
    max_attempts: int = 10_000_000,
) -> NDArray[np.float64]:
    """Sprinkle events uniformly with respect to conformal physical volume."""

    count = int(n)
    if count <= 0:
        raise ValueError("n must be positive")
    duration = _validate_t(T)
    batch = int(batch_size)
    if batch <= 0:
        raise ValueError("batch_size must be positive")
    attempts_limit = int(max_attempts)
    if attempts_limit < count:
        raise ValueError("max_attempts must be at least n")

    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    weight_max = estimate_conformal_weight_max_1p1(duration, profile)
    accepted: list[NDArray[np.float64]] = []
    accepted_count = 0
    attempts = 0
    while accepted_count < count and attempts < attempts_limit:
        proposal_count = min(batch, attempts_limit - attempts)
        candidates = sprinkle_1p1_causal_diamond(proposal_count, T=duration, seed=rng)
        weights = conformal_volume_density_1p1(profile, candidates)
        keep = rng.uniform(0.0, 1.0, proposal_count) <= weights / weight_max
        kept = candidates[keep]
        if kept.size:
            accepted.append(kept)
            accepted_count += kept.shape[0]
        attempts += proposal_count

    if accepted_count < count:
        raise RuntimeError("could not sample requested events within max_attempts")
    events = np.vstack(accepted)[:count].astype(np.float64, copy=False)
    inside = np.abs(events[:, 1]) <= duration / 2.0 - np.abs(events[:, 0]) + 1e-12
    if not np.all(inside):
        raise RuntimeError("generated event outside causal diamond")
    return events


def estimate_time_bin_density_shape_1p1(
    events: ArrayLike,
    T: float,
    bin_edges: ArrayLike,
) -> NDArray[np.float64]:
    """Estimate coordinate density in each time bin as count/bin volume."""

    event_array = np.asarray(events, dtype=float)
    if event_array.ndim != 2 or event_array.shape[1] != 2:
        raise ValueError("events must have shape (n, 2)")
    volumes = coordinate_time_bin_volumes_1p1(T, bin_edges)
    counts, _ = np.histogram(event_array[:, 0], bins=np.asarray(bin_edges, dtype=float))
    densities = np.full(volumes.shape, np.nan, dtype=np.float64)
    positive = volumes > 0.0
    densities[positive] = counts[positive] / volumes[positive]
    return densities


def normalize_profile_shape(values: ArrayLike) -> NDArray[np.float64]:
    """Normalize finite positive values by their mean."""

    array = np.asarray(values, dtype=float)
    finite_positive = np.isfinite(array) & (array > 0.0)
    if not np.any(finite_positive):
        raise ValueError("values must contain at least one finite positive entry")
    mean_value = float(np.mean(array[finite_positive]))
    normalized = np.full(array.shape, np.nan, dtype=np.float64)
    normalized[finite_positive] = array[finite_positive] / mean_value
    return normalized
