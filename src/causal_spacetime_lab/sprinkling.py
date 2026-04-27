"""Uniform sprinkling routines for simple causal-set experiments."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def _sample_points_in_ball(
    radii: NDArray[np.float64],
    spatial_dim: int,
    rng: np.random.Generator,
) -> NDArray[np.float64]:
    if spatial_dim <= 0:
        return np.empty((radii.size, 0), dtype=np.float64)

    directions = rng.normal(size=(radii.size, spatial_dim))
    norms = np.linalg.norm(directions, axis=1)
    zero_norms = norms == 0.0
    while np.any(zero_norms):
        directions[zero_norms] = rng.normal(
            size=(np.count_nonzero(zero_norms), spatial_dim)
        )
        norms = np.linalg.norm(directions, axis=1)
        zero_norms = norms == 0.0

    unit_directions = directions / norms[:, None]
    radial_fractions = rng.uniform(0.0, 1.0, size=radii.size) ** (1.0 / spatial_dim)
    return unit_directions * (radii * radial_fractions)[:, None]


def sprinkle_minkowski_causal_diamond(
    n: int,
    spacetime_dim: int,
    T: float = 2.0,
    seed: int | np.random.Generator | None = None,
) -> NDArray[np.float64]:
    """Sample events uniformly inside a D-dimensional Minkowski causal diamond.

    The interval endpoints are ``(-T/2, 0, ..., 0)`` and ``(T/2, 0, ..., 0)``.
    Events satisfy ``||x_vec|| <= T/2 - abs(t)``.
    """

    if n < 0:
        raise ValueError("n must be non-negative")
    if spacetime_dim < 2:
        raise ValueError("spacetime_dim must be at least 2")
    if T <= 0:
        raise ValueError("T must be positive")

    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    if n == 0:
        return np.empty((0, spacetime_dim), dtype=np.float64)

    spatial_dim = spacetime_dim - 1
    half_T = T / 2.0
    accepted_t: list[NDArray[np.float64]] = []
    accepted_count = 0

    while accepted_count < n:
        remaining = n - accepted_count
        batch_size = max(4 * remaining, 256)
        t_candidates = rng.uniform(-half_T, half_T, size=batch_size)
        radii = half_T - np.abs(t_candidates)
        accept_probability = (radii / half_T) ** spatial_dim
        keep = rng.uniform(0.0, 1.0, size=batch_size) <= accept_probability
        kept_t = t_candidates[keep]
        if kept_t.size == 0:
            continue
        accepted_t.append(kept_t)
        accepted_count += kept_t.size

    t = np.concatenate(accepted_t)[:n]
    radii = half_T - np.abs(t)
    spatial = _sample_points_in_ball(radii, spatial_dim, rng)
    return np.column_stack((t, spatial)).astype(np.float64, copy=False)


def sprinkle_1p1_causal_diamond(
    n: int,
    T: float = 2.0,
    seed: int | np.random.Generator | None = None,
) -> NDArray[np.float64]:
    """Sample events uniformly inside a 1+1D Minkowski causal diamond."""

    return sprinkle_minkowski_causal_diamond(n, spacetime_dim=2, T=T, seed=seed)


def sprinkle_1p1_forward_cone(
    n: int,
    T: float,
    seed: int | np.random.Generator | None = None,
) -> NDArray[np.float64]:
    """Sample events uniformly in ``0 <= t <= T`` and ``abs(x) <= t``."""

    if n < 0:
        raise ValueError("n must be non-negative")
    if T <= 0:
        raise ValueError("T must be positive")

    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    t = T * np.sqrt(rng.uniform(0.0, 1.0, size=n))
    x = rng.uniform(-t, t)
    return np.column_stack((t, x)).astype(np.float64, copy=False)
