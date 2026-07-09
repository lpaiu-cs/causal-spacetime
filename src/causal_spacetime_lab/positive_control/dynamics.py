"""Geometry-free causal-set dynamics for the P3 emergence experiment.

Transitive percolation (the simplest classical-sequential-growth dynamics,
Rideout-Sorkin 2000): elements are grown in sequence and each new element is
linked to each earlier element independently with probability p; the relation
is then transitively closed. The dynamics never refers to any embedding
geometry. See docs/prereg/p3_dynamics_emergence.md.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.positive_control.rewire import transitive_closure


def transitive_percolation(
    n: int, p: float, seed: int
) -> tuple[NDArray[np.bool_], NDArray[np.int_]]:
    """Return a transitive-percolation causal order and its growth index.

    ``p`` is the pre-closure link probability between earlier and later growth
    elements. Growth order (index order) is a valid linear extension of the
    resulting order and is returned as the intrinsic time label.
    """

    if n <= 0:
        raise ValueError("n must be positive")
    if not 0.0 <= p <= 1.0:
        raise ValueError("p must be in [0, 1]")
    rng = np.random.default_rng(seed)
    links = np.triu(rng.uniform(size=(n, n)) < p, k=1)
    closed = transitive_closure(links)
    np.fill_diagonal(closed, False)
    return closed, np.arange(n, dtype=int)


def relation_density(causal: NDArray[np.bool_], times: NDArray[np.float64]) -> float:
    """Relation density among time-ordered pairs (matches the P1/P2 convention)."""

    ordered = times[:, None] < times[None, :]
    ordered_count = int(np.sum(ordered))
    if ordered_count == 0:
        raise ValueError("no time-ordered pairs")
    return float(np.sum(causal & ordered) / ordered_count)
