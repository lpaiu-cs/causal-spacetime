from __future__ import annotations

import numpy as np

from causal_spacetime_lab.order_constraints import (
    quadruplet_constraints_from_order_values,
    quadruplet_constraints_from_pair_distance_order,
    quadruplet_constraints_from_successor_ticks,
    scalar_order_constraints_from_values,
)
from causal_spacetime_lab.ordinal_embedding import quadruplet_violation_rate


def test_quadruplet_constraints_from_order_values_are_valid() -> None:
    values = np.asarray([0.0, 1.0, 3.0, 7.0])

    constraints = quadruplet_constraints_from_order_values(values, 20, seed=1)

    assert constraints.shape == (20, 4)
    assert quadruplet_violation_rate(values[:, None], constraints) == 0.0


def test_quadruplet_constraints_from_pair_distance_order_are_valid() -> None:
    pair_indices = np.asarray([[0, 1], [0, 2], [0, 3]], dtype=int)
    pair_values = np.asarray([1.0, 2.0, 4.0])

    constraints = quadruplet_constraints_from_pair_distance_order(
        pair_indices,
        pair_values,
        10,
        seed=2,
    )

    assert constraints.shape == (10, 4)
    assert np.all(constraints >= 0)


def test_scalar_order_constraints_from_values_are_valid() -> None:
    values = np.asarray([0.0, 2.0, 5.0])

    constraints = scalar_order_constraints_from_values(values, 8, seed=3)

    assert constraints.shape == (8, 2)
    assert np.all(values[constraints[:, 0]] < values[constraints[:, 1]])


def test_quadruplet_constraints_from_successor_ticks_returns_scalar_order() -> None:
    ticks = np.asarray([4, 2, 7, -1])
    accessible = np.asarray([True, True, True, False])

    constraints = quadruplet_constraints_from_successor_ticks(
        ticks,
        accessible,
        8,
        seed=4,
    )

    assert constraints.shape == (8, 2)
    assert np.all(ticks[constraints[:, 0]] < ticks[constraints[:, 1]])
