from __future__ import annotations

from math import factorial

import numpy as np
import pytest

from causal_spacetime_lab.identity_constraints import (
    best_persistence_matchings_with_fixed_points,
    filter_permutations_by_fixed_points,
    generate_partial_identity_constraints,
)
from causal_spacetime_lab.persistence_definability import (
    object_identity_without_persistence,
    pair_distance_history_without_persistence,
)
from causal_spacetime_lab.persistence_history import (
    build_identity_tracks_from_adjacent_matchings,
    infer_adjacent_matchings,
    relational_history_from_inferred_tracks,
    track_consistency_error,
)
from causal_spacetime_lab.persistence_matching import (
    all_permutations,
    apply_permutation_to_positions,
    best_persistence_matchings_by_relational_order,
    matching_accuracy,
    matching_ambiguity_gap,
    signature_disagreement_for_permutation,
)
from causal_spacetime_lab.persistence_scenarios import (
    generate_unlabeled_crossing_history_1d,
    generate_unlabeled_small_motion_history,
    generate_unlabeled_static_history,
)


def test_object_identity_without_persistence_returns_undefined() -> None:
    report = object_identity_without_persistence()

    assert report.defined is False
    assert "undefined" in report.reason


def test_pair_distance_history_without_persistence_returns_undefined() -> None:
    report = pair_distance_history_without_persistence()

    assert report.defined is False
    assert "requires" in report.reason


def test_all_permutations_returns_correct_number_for_n3() -> None:
    permutations = all_permutations(3)

    assert permutations.shape == (factorial(3), 3)


def test_all_permutations_rejects_n_larger_than_max() -> None:
    with pytest.raises(ValueError):
        all_permutations(4, max_n=3)


def test_apply_permutation_to_positions_deterministic_example() -> None:
    values = apply_permutation_to_positions(
        np.asarray([10.0, 20.0, 30.0]),
        np.asarray([2, 0, 1]),
    )

    assert np.array_equal(values, [30.0, 10.0, 20.0])


def test_signature_disagreement_for_permutation_identity_example() -> None:
    positions = np.asarray([0.0, 1.0, 3.0])

    value = signature_disagreement_for_permutation(
        positions,
        positions,
        np.asarray([0, 1, 2]),
    )

    assert value == pytest.approx(0.0)


def test_best_persistence_matchings_returns_sorted_costs() -> None:
    matches = best_persistence_matchings_by_relational_order(
        np.asarray([0.0, 1.0, 3.0]),
        np.asarray([1.0, 3.0, 0.0]),
        top_k=4,
    )

    costs = [float(row["cost"]) for row in matches]
    assert costs == sorted(costs)


def test_matching_ambiguity_gap() -> None:
    gap = matching_ambiguity_gap([{"cost": 0.1}, {"cost": 0.4}])

    assert gap == pytest.approx(0.3)


def test_matching_accuracy() -> None:
    accuracy = matching_accuracy(np.asarray([0, 2, 1]), np.asarray([0, 1, 1]))

    assert accuracy == pytest.approx(2 / 3)


def test_filter_permutations_by_fixed_points() -> None:
    filtered = filter_permutations_by_fixed_points(all_permutations(3), {0: 2})

    assert filtered.shape[0] == 2
    assert np.all(filtered[:, 0] == 2)


def test_best_persistence_matchings_with_fixed_points() -> None:
    matches = best_persistence_matchings_with_fixed_points(
        np.asarray([0.0, 1.0, 3.0]),
        np.asarray([3.0, 0.0, 1.0]),
        {0: 1},
        top_k=3,
    )

    assert matches
    assert all(int(np.asarray(row["permutation"])[0]) == 1 for row in matches)


def test_generate_partial_identity_constraints() -> None:
    fixed = generate_partial_identity_constraints(
        np.asarray([2, 0, 1, 3]),
        0.5,
        seed=2,
    )

    assert len(fixed) == 2
    assert all(value in {0, 1, 2, 3} for value in fixed.values())


def test_infer_adjacent_matchings() -> None:
    positions = {0: np.asarray([0.0, 1.0]), 1: np.asarray([1.0, 0.0])}

    inferred = infer_adjacent_matchings(positions, top_k=2)

    assert set(inferred) == {(0, 1)}
    assert inferred[(0, 1)]


def test_build_identity_tracks_from_adjacent_matchings() -> None:
    tracks = build_identity_tracks_from_adjacent_matchings(
        {(0, 1): np.asarray([2, 0, 1])},
        np.asarray([10, 20, 30]),
    )

    assert tracks[0] == {10: 0, 20: 1, 30: 2}
    assert tracks[1] == {10: 2, 20: 0, 30: 1}


def test_relational_history_from_inferred_tracks() -> None:
    positions = {0: np.asarray([0.0, 1.0]), 1: np.asarray([2.0, 3.0])}
    tracks = {0: {7: 0, 8: 1}, 1: {7: 1, 8: 0}}

    history = relational_history_from_inferred_tracks(positions, tracks)

    assert history[1][7] == pytest.approx(3.0)
    assert history[1][8] == pytest.approx(2.0)


def test_track_consistency_error() -> None:
    inferred = {0: {0: 0, 1: 1}, 1: {0: 1, 1: 0}}
    truth = {0: {0: 0, 1: 1}, 1: {0: 0, 1: 1}}

    assert track_consistency_error(inferred, truth) == pytest.approx(0.5)


def test_generate_unlabeled_static_history() -> None:
    positions, tracks = generate_unlabeled_static_history(
        3,
        np.asarray([0.0, 1.0, 3.0]),
        seed=1,
    )

    assert set(positions) == {0, 1, 2}
    assert all(values.shape == (3,) for values in positions.values())
    assert set(tracks[0]) == {0, 1, 2}


def test_generate_unlabeled_small_motion_history() -> None:
    positions, tracks = generate_unlabeled_small_motion_history(
        2,
        np.asarray([0.0, 1.0]),
        np.asarray([0.1, -0.1]),
        random_per_slice_permutation=False,
    )

    assert positions[1][0] == pytest.approx(0.1)
    assert tracks[1] == {0: 0, 1: 1}


def test_generate_unlabeled_crossing_history_1d() -> None:
    positions, tracks = generate_unlabeled_crossing_history_1d(
        3,
        np.asarray([-1.0, 0.0, 1.0]),
        (0, 2),
        random_per_slice_permutation=False,
    )

    assert positions[0][0] == pytest.approx(-1.0)
    assert positions[2][0] == pytest.approx(1.0)
    assert tracks[0] == {0: 0, 1: 1, 2: 2}
