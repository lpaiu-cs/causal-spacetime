from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.predicate_definability import predicate_requirement_table
from causal_spacetime_lab.relational_evolution import (
    apply_per_slice_affine_position_gauge,
    compare_histories_order_disagreement,
    pair_distance_order_signature,
    pair_order_comparison_matrix,
    relational_change_rate_between_slices,
    relational_shape_history,
    signature_order_disagreement,
    unordered_object_pairs,
)
from causal_spacetime_lab.relational_scenarios import (
    add_position_noise_to_history,
    generate_expanding_configuration_history,
    generate_shear_or_reordering_history_1d,
    generate_static_configuration_history,
)


def test_predicate_requirement_table_contains_expected_predicates() -> None:
    table = predicate_requirement_table()
    by_name = {str(row["name"]): row for row in table}

    assert "pair_distance_order_history" in by_name
    assert "coordinate_velocity" in by_name
    assert by_name["pair_distance_order_history"]["requires_persistence"] is True
    assert by_name["pair_distance_order_history"]["requires_transport"] is False
    assert by_name["coordinate_velocity"]["requires_transport"] is True


def test_unordered_object_pairs_deterministic_example() -> None:
    pairs = unordered_object_pairs(np.asarray([2, 0, 1]))

    assert np.array_equal(pairs, [[0, 1], [0, 2], [1, 2]])


def test_pair_distance_order_signature_deterministic_example() -> None:
    signature = pair_distance_order_signature({0: 0.0, 1: 1.0, 2: 3.0})

    assert np.array_equal(signature.object_pairs, [[0, 1], [0, 2], [1, 2]])
    assert np.array_equal(signature.pair_distance_values, [1.0, 3.0, 2.0])
    assert np.array_equal(signature.pair_order_ranks, [0.0, 2.0, 1.0])


def test_pair_order_comparison_matrix_deterministic_example() -> None:
    signature = pair_distance_order_signature({0: 0.0, 1: 1.0, 2: 3.0})

    matrix = pair_order_comparison_matrix(signature)

    assert matrix[0, 1] == -1
    assert matrix[1, 0] == 1
    assert matrix[0, 0] == 0


def test_signature_order_disagreement_identical_is_zero() -> None:
    signature = pair_distance_order_signature({0: 0.0, 1: 1.0, 2: 3.0})

    assert signature_order_disagreement(signature, signature) == pytest.approx(0.0)


def test_relational_shape_history_returns_one_signature_per_slice() -> None:
    history = generate_static_configuration_history(
        np.asarray([0, 1]),
        np.asarray([0.0, 1.0, 3.0]),
    )

    signatures = relational_shape_history(history)

    assert set(signatures) == {0, 1}
    assert all(
        signature.object_pairs.shape == (3, 2)
        for signature in signatures.values()
    )


def test_relational_change_rate_static_history_is_zero() -> None:
    history = relational_shape_history(
        generate_static_configuration_history(
            np.asarray([0, 1, 2]),
            np.asarray([0.0, 1.0, 3.0]),
        )
    )

    rows = relational_change_rate_between_slices(history)

    assert rows
    assert all(row["order_disagreement"] == pytest.approx(0.0) for row in rows)


def test_generate_expanding_configuration_preserves_pair_distance_order() -> None:
    static = relational_shape_history(
        generate_static_configuration_history(
            np.asarray([0, 1, 2]),
            np.asarray([0.0, 1.0, 3.0]),
        )
    )
    expanding = relational_shape_history(
        generate_expanding_configuration_history(
            np.asarray([0, 1, 2]),
            np.asarray([0.0, 1.0, 3.0]),
            np.asarray([1.0, 2.0, 4.0]),
        )
    )

    assert compare_histories_order_disagreement(static, expanding) == pytest.approx(0.0)


def test_generate_shear_or_reordering_can_produce_nonzero_change() -> None:
    history = relational_shape_history(
        generate_shear_or_reordering_history_1d(
            np.asarray([0, 1, 2]),
            np.asarray([0.0, 1.0, 3.0]),
            moving_object_id=1,
            displacement_by_slice=np.asarray([0.0, 1.5, 3.0]),
        )
    )

    rows = relational_change_rate_between_slices(history)

    assert any(row["order_disagreement"] > 0.0 for row in rows)


def test_apply_per_slice_affine_position_gauge_preserves_relational_history() -> None:
    positions = generate_shear_or_reordering_history_1d(
        np.asarray([0, 1, 2]),
        np.asarray([0.0, 1.0, 3.0]),
        moving_object_id=1,
        displacement_by_slice=np.asarray([0.0, 1.5, 3.0]),
    )
    base = relational_shape_history(positions)
    gauged = relational_shape_history(
        apply_per_slice_affine_position_gauge(positions, seed=11)
    )

    assert compare_histories_order_disagreement(base, gauged) == pytest.approx(0.0)


def test_add_position_noise_to_history_preserves_shape_and_keys() -> None:
    history = generate_static_configuration_history(
        np.asarray([0, 1]),
        np.asarray([0.0, 1.0, 3.0]),
    )

    noisy = add_position_noise_to_history(history, noise_scale=0.01, seed=7)

    assert set(noisy) == set(history)
    assert all(set(noisy[label]) == set(history[label]) for label in history)
