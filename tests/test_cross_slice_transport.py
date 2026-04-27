from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.cross_slice import (
    same_direction_without_transport,
    same_position_without_transport,
    spatial_evolution_without_transport,
    velocity_without_transport,
)
from causal_spacetime_lab.persistence import (
    finite_difference_velocity_by_object,
    generate_persistent_object_events_1p1,
)
from causal_spacetime_lab.slice_anchors import (
    anchor_indices_by_slice,
    anchor_reference_positions_by_slice,
    generate_stationary_anchor_events_1p1,
)
from causal_spacetime_lab.slice_constraint_graph import (
    component_labels_from_constraints,
    constraint_cross_slice_fraction,
    constraint_point_graph_components,
)
from causal_spacetime_lab.slice_local_embedding import (
    assemble_slice_embeddings_with_nan,
    constraints_for_slice,
    fit_slice_local_ordinal_embeddings,
    remap_constraints_to_local_indices,
)
from causal_spacetime_lab.slice_local_validation import (
    summarize_slice_local_validation,
    validate_slice_local_embeddings_against_true_positions,
)
from causal_spacetime_lab.slice_transport import (
    SliceTransportRule,
    align_slices_by_anchor_points,
    apply_random_slice_gauge_transforms,
    apply_slice_transport_rule,
    fit_1d_affine_from_anchors,
    same_position_under_transport,
)


def test_without_transport_predicates_are_undefined() -> None:
    reports = [
        same_position_without_transport(),
        same_direction_without_transport(),
        velocity_without_transport(),
        spatial_evolution_without_transport(),
    ]

    assert all(not report.defined for report in reports)
    assert all(report.value is None for report in reports)
    assert all("undefined" in report.reason.lower() for report in reports)


def test_constraint_point_graph_components_deterministic() -> None:
    constraints = np.asarray([[0, 1, 1, 2], [3, 4, 4, 5]])

    components = constraint_point_graph_components(7, constraints)

    assert [component.tolist() for component in components] == [[0, 1, 2], [3, 4, 5]]


def test_component_labels_from_constraints() -> None:
    constraints = np.asarray([[0, 1, 1, 2], [3, 4, 4, 5]])

    labels = component_labels_from_constraints(7, constraints)

    assert np.array_equal(labels, [0, 0, 0, 1, 1, 1, -1])


def test_constraint_cross_slice_fraction() -> None:
    labels = np.asarray([0, 0, 0, 1, 1, 1])
    same_slice = np.asarray([[0, 1, 1, 2], [3, 4, 4, 5]])
    cross_slice = np.asarray([[0, 1, 3, 4], [3, 4, 4, 5]])

    assert constraint_cross_slice_fraction(same_slice, labels) == pytest.approx(0.0)
    assert constraint_cross_slice_fraction(cross_slice, labels) == pytest.approx(0.5)


def test_constraints_for_slice_and_remap() -> None:
    labels = np.asarray([0, 0, 0, 1])
    constraints = np.asarray([[0, 1, 1, 2], [0, 1, 2, 3]])

    local = constraints_for_slice(constraints, labels, 0)
    remapped, mapping = remap_constraints_to_local_indices(local, np.asarray([0, 1, 2]))

    assert np.array_equal(local, [[0, 1, 1, 2]])
    assert np.array_equal(remapped, [[0, 1, 1, 2]])
    assert mapping == {0: 0, 1: 1, 2: 2}


def test_fit_slice_local_ordinal_embeddings_tiny_example() -> None:
    labels = np.asarray([0, 0, 0, 0, 1])
    constraints = np.asarray(
        [
            [0, 1, 0, 2],
            [0, 1, 0, 3],
            [1, 2, 0, 3],
            [2, 3, 0, 3],
        ]
    )

    embeddings = fit_slice_local_ordinal_embeddings(
        5,
        labels,
        constraints,
        min_points_per_slice=3,
        min_constraints_per_slice=2,
        steps=10,
        restarts=1,
        seed=2,
    )

    assert set(embeddings) == {0}
    assert np.asarray(embeddings[0]["embedding"]).shape == (4, 1)


def test_assemble_slice_embeddings_with_nan() -> None:
    slice_embeddings = {
        2: {
            "embedding": np.asarray([[0.0], [1.0]]),
            "global_indices": np.asarray([1, 3]),
        }
    }

    assembled = assemble_slice_embeddings_with_nan(5, slice_embeddings)

    assert assembled.shape == (5, 1)
    assert np.isnan(assembled[0, 0])
    assert assembled[1, 0] == pytest.approx(0.0)
    assert assembled[3, 0] == pytest.approx(1.0)


def test_validate_slice_local_embeddings_against_true_positions() -> None:
    slice_embeddings = {
        0: {
            "embedding": np.asarray([[0.0], [1.0], [2.0]]),
            "global_indices": np.asarray([0, 1, 2]),
            "constraints": np.empty((0, 4), dtype=int),
        }
    }

    rows = validate_slice_local_embeddings_against_true_positions(
        slice_embeddings,
        np.asarray([10.0, 12.0, 14.0]),
    )
    summary = summarize_slice_local_validation(rows)

    assert rows[0]["local_rmse"] < 1e-12
    assert summary["slice_count"] == 1.0


def test_apply_slice_transport_rule_deterministic() -> None:
    rule = SliceTransportRule(
        name="test",
        description="test rule",
        slice_labels=np.asarray([0, 0, 1]),
        scale_by_slice={0: 2.0, 1: 3.0},
        shift_by_slice={0: 1.0, 1: -1.0},
        reflection_by_slice={0: 1.0, 1: -1.0},
    )

    transported = apply_slice_transport_rule(
        np.asarray([0.0, 1.0, 2.0]),
        np.asarray([0, 0, 1]),
        rule,
    )

    assert np.array_equal(transported, [1.0, 3.0, -7.0])


def test_apply_random_slice_gauge_transforms_preserves_shape() -> None:
    coords = np.asarray([0.0, 1.0, 2.0, 3.0])
    labels = np.asarray([0, 0, 1, 1])

    transformed = apply_random_slice_gauge_transforms(coords, labels, seed=5)

    assert transformed.shape == coords.shape
    assert np.all(np.isfinite(transformed))


def test_same_position_under_transport_requires_covered_slices() -> None:
    rule = SliceTransportRule(
        name="partial",
        description="covers one slice",
        slice_labels=np.asarray([0, 1]),
        scale_by_slice={0: 1.0},
        shift_by_slice={0: 0.0},
        reflection_by_slice={0: 1.0},
    )

    covered = same_position_under_transport(1.0, 0, 1.0, 0, rule)
    missing = same_position_under_transport(1.0, 0, 1.0, 1, rule)

    assert covered.defined
    assert covered.value is True
    assert not missing.defined


def test_fit_1d_affine_from_anchors_recovers_map() -> None:
    scale, offset, rmse = fit_1d_affine_from_anchors(
        np.asarray([0.0, 1.0, 2.0]),
        np.asarray([1.0, 3.0, 5.0]),
    )

    assert scale == pytest.approx(2.0)
    assert offset == pytest.approx(1.0)
    assert rmse < 1e-12


def test_align_slices_by_anchor_points_deterministic() -> None:
    coords = np.asarray([0.0, 1.0, 5.0, 6.0])
    labels = np.asarray([0, 0, 1, 1])

    aligned = align_slices_by_anchor_points(
        coords,
        labels,
        {0: np.asarray([0, 1]), 1: np.asarray([2, 3])},
        {0: np.asarray([10.0, 12.0]), 1: np.asarray([20.0, 22.0])},
    )

    assert np.allclose(aligned, [10.0, 12.0, 20.0, 22.0])


def test_generate_stationary_anchor_events_1p1_shape_and_labels() -> None:
    events, slice_ids, anchor_ids = generate_stationary_anchor_events_1p1(
        np.asarray([0.0, 1.0]),
        np.asarray([-1.0, 1.0]),
    )

    assert events.shape == (4, 2)
    assert np.array_equal(slice_ids, [0, 0, 1, 1])
    assert np.array_equal(anchor_ids, [0, 1, 0, 1])
    assert anchor_indices_by_slice(np.arange(4), slice_ids)[1].tolist() == [2, 3]
    refs = anchor_reference_positions_by_slice(
        anchor_ids,
        np.asarray([-1.0, 1.0]),
        slice_ids,
    )
    assert np.array_equal(refs[0], [-1.0, 1.0])


def test_generate_persistent_object_events_and_velocity() -> None:
    events, slice_ids, object_ids = generate_persistent_object_events_1p1(
        np.asarray([0.0, 1.0, 2.0]),
        np.asarray([1.0, -1.0]),
        np.asarray([0.5, -0.25]),
    )

    velocities = finite_difference_velocity_by_object(
        events[:, 1],
        slice_ids,
        object_ids,
        {0: 0.0, 1: 1.0, 2: 2.0},
    )

    assert events.shape == (6, 2)
    assert set(object_ids) == {0, 1}
    assert velocities[0] == pytest.approx(0.5)
    assert velocities[1] == pytest.approx(-0.25)
