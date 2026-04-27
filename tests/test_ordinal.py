from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.ordinal import (
    order_inversion_rate,
    order_matrix_from_values,
    pair_distance_order_inversion_rate,
    pair_distance_values_1d,
)


def test_order_matrix_from_values_strict_and_nonstrict() -> None:
    values = np.asarray([2.0, 1.0, 2.0])

    strict = order_matrix_from_values(values, strict=True)
    nonstrict = order_matrix_from_values(values, strict=False)

    assert strict.tolist() == [
        [False, False, False],
        [True, False, True],
        [False, False, False],
    ]
    assert np.all(np.diag(nonstrict))


def test_order_inversion_rate_identical_values_is_zero() -> None:
    values = np.asarray([1.0, 2.0, 4.0])

    assert order_inversion_rate(values, values) == pytest.approx(0.0)


def test_order_inversion_rate_positive_scaling_is_zero() -> None:
    values = np.asarray([1.0, 2.0, 4.0])

    assert order_inversion_rate(values, 10.0 * values) == pytest.approx(0.0)


def test_order_inversion_rate_reversed_values_is_one() -> None:
    values = np.asarray([1.0, 2.0, 4.0])

    assert order_inversion_rate(values, values[::-1]) == pytest.approx(1.0)


def test_pair_distance_values_1d_hand_coded() -> None:
    coords = np.asarray([0.0, 2.0, 5.0])
    pairs = np.asarray([[0, 1], [1, 2], [0, 2]])

    distances = pair_distance_values_1d(coords, pairs)

    assert distances.tolist() == pytest.approx([2.0, 3.0, 5.0])


def test_pair_distance_order_inversion_rate_hand_coded() -> None:
    true_coords = np.asarray([0.0, 1.0, 3.0])
    estimated_coords = np.asarray([0.0, 3.0, 4.0])
    pairs = np.asarray([[0, 1], [1, 2], [0, 2]])

    inversion = pair_distance_order_inversion_rate(
        true_coords,
        estimated_coords,
        pairs,
    )

    assert inversion > 0.0
