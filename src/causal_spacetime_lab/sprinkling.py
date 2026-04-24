"""Uniform sprinkling routines for simple causal-set experiments."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def sprinkle_1p1_causal_diamond(
    n: int,
    T: float = 2.0,
    seed: int | np.random.Generator | None = None,
) -> NDArray[np.float64]:
    """Sample events uniformly inside a 1+1D Minkowski causal diamond.

    The diamond has boundary events ``(-T/2, 0)`` and ``(T/2, 0)`` and is
    equivalently described by ``abs(x) <= T/2 - abs(t)``.

    Sampling is performed in null coordinates ``u = t + x`` and ``v = t - x``.
    The causal diamond maps to the square ``[-T/2, T/2]^2`` in ``(u, v)``,
    making uniform sampling straightforward.
    """

    if n < 0:
        raise ValueError("n must be non-negative")
    if T <= 0:
        raise ValueError("T must be positive")

    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    lower = -T / 2.0
    upper = T / 2.0
    u = rng.uniform(lower, upper, size=n)
    v = rng.uniform(lower, upper, size=n)

    t = 0.5 * (u + v)
    x = 0.5 * (u - v)
    return np.column_stack((t, x)).astype(np.float64, copy=False)


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

