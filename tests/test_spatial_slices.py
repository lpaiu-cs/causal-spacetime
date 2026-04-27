from __future__ import annotations

import numpy as np

from causal_spacetime_lab.spatial_slices import (
    assign_slices_from_radar_time_rank,
    filter_pairs_by_same_slice,
    radar_distance_rank_from_tick_brackets,
    radar_time_rank_from_tick_brackets,
    same_slice_unordered_pairs,
    slice_sizes,
)


def test_radar_time_rank_from_tick_brackets() -> None:
    predecessor = np.asarray([1, 2, -1])
    successor = np.asarray([3, 5, -1])
    accessible = np.asarray([True, True, False])

    ranks = radar_time_rank_from_tick_brackets(predecessor, successor, accessible)

    assert np.array_equal(ranks, [4, 7, -1])


def test_radar_distance_rank_from_tick_brackets() -> None:
    predecessor = np.asarray([1, 2, -1])
    successor = np.asarray([3, 5, -1])
    accessible = np.asarray([True, True, False])

    ranks = radar_distance_rank_from_tick_brackets(predecessor, successor, accessible)

    assert np.array_equal(ranks, [2, 3, -1])


def test_assign_slices_from_radar_time_rank() -> None:
    ranks = np.asarray([0, 1, 2, 5, -1])
    accessible = np.asarray([True, True, True, True, False])

    labels = assign_slices_from_radar_time_rank(ranks, accessible, bin_width=2)

    assert np.array_equal(labels, [0, 0, 1, 2, -1])


def test_slice_sizes_ignores_inaccessible_label() -> None:
    labels = np.asarray([0, 0, 1, -1, 1, 1])

    assert slice_sizes(labels) == {0: 2, 1: 3}


def test_same_slice_unordered_pairs() -> None:
    labels = np.asarray([0, 0, 1, 1, 1, -1])

    pairs = same_slice_unordered_pairs(labels, max_pairs_per_slice=2, seed=1)

    assert pairs.shape == (3, 2)
    assert np.all(labels[pairs[:, 0]] == labels[pairs[:, 1]])
    assert np.all(pairs[:, 0] < pairs[:, 1])


def test_filter_pairs_by_same_slice() -> None:
    labels = np.asarray([0, 0, 1, -1])
    pairs = np.asarray([[0, 1], [0, 2], [2, 3]])

    filtered = filter_pairs_by_same_slice(pairs, labels)

    assert np.array_equal(filtered, [[0, 1]])
