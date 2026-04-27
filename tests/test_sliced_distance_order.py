from __future__ import annotations

import numpy as np

from causal_spacetime_lab.sliced_distance_order import (
    pair_slice_labels,
    quadruplet_constraints_from_sliced_pair_distances,
    sliced_pair_distance_order_inversion_rate,
)


def test_pair_slice_labels() -> None:
    pairs = np.asarray([[0, 1], [0, 2], [2, 3]])
    labels = np.asarray([0, 0, 1, 1])

    result = pair_slice_labels(pairs, labels)

    assert np.array_equal(result, [0, -1, 1])


def test_sliced_pair_distance_order_inversion_identical_is_zero() -> None:
    coords = np.asarray([0.0, 1.0, 3.0, 6.0])
    pairs = np.asarray([[0, 1], [0, 2], [0, 3]])

    error = sliced_pair_distance_order_inversion_rate(coords, coords, pairs)

    assert error == 0.0


def test_sliced_pair_distance_order_inversion_translated_is_zero() -> None:
    coords = np.asarray([0.0, 1.0, 3.0, 6.0])
    translated = coords + 10.0
    pairs = np.asarray([[0, 1], [0, 2], [0, 3]])

    error = sliced_pair_distance_order_inversion_rate(coords, translated, pairs)

    assert error == 0.0


def test_quadruplet_constraints_from_sliced_pair_distances_respects_slices() -> None:
    pairs = np.asarray([[0, 1], [0, 2], [2, 3], [2, 4]])
    values = np.asarray([1.0, 3.0, 2.0, 5.0])
    pair_labels = np.asarray([0, 0, 1, 1])

    constraints = quadruplet_constraints_from_sliced_pair_distances(
        pairs,
        values,
        pair_labels,
        20,
        seed=3,
    )

    assert constraints.shape == (20, 4)
    original_pairs = {
        tuple(pair): label for pair, label in zip(pairs, pair_labels, strict=True)
    }
    for row in constraints:
        assert original_pairs[tuple(row[:2])] == original_pairs[tuple(row[2:])]
