"""Exact sanity checks for cross-slice transport utilities."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.cross_slice import (
    same_position_without_transport,
    velocity_without_transport,
)
from causal_spacetime_lab.persistence import (
    finite_difference_velocity_by_object,
    generate_persistent_object_events_1p1,
)
from causal_spacetime_lab.slice_anchors import generate_stationary_anchor_events_1p1
from causal_spacetime_lab.slice_constraint_graph import (
    component_labels_from_constraints,
    constraint_cross_slice_fraction,
    constraint_point_graph_components,
)
from causal_spacetime_lab.slice_local_embedding import (
    constraints_for_slice,
    remap_constraints_to_local_indices,
)
from causal_spacetime_lab.slice_transport import (
    SliceTransportRule,
    align_slices_by_anchor_points,
    apply_random_slice_gauge_transforms,
    fit_1d_affine_from_anchors,
)

DEFAULT_OUTPUT = Path("outputs/data/cross_slice_transport_exact_sanity.csv")


def _row(check: str, passed: bool, value: float) -> dict[str, float | str]:
    return {"check": check, "passed": float(passed), "value": float(value)}


def run_experiment() -> list[dict[str, float | str]]:
    """Run deterministic exact checks."""

    constraints = np.asarray([[0, 1, 1, 2], [3, 4, 4, 5]], dtype=int)
    labels = np.asarray([0, 0, 0, 1, 1, 1])
    components = constraint_point_graph_components(6, constraints)
    component_labels = component_labels_from_constraints(6, constraints)
    local_constraints = constraints_for_slice(constraints, labels, 0)
    remapped, mapping = remap_constraints_to_local_indices(local_constraints, [0, 1, 2])
    coords = np.asarray([0.0, 1.0, 2.0, -1.0, 0.0, 1.0])
    random_coords = apply_random_slice_gauge_transforms(coords, labels, seed=1)
    scale, offset, rmse = fit_1d_affine_from_anchors(
        np.asarray([0.0, 2.0]),
        np.asarray([1.0, 5.0]),
    )
    anchor_indices = {0: np.asarray([0, 1]), 1: np.asarray([3, 4])}
    anchor_refs = {0: np.asarray([10.0, 12.0]), 1: np.asarray([20.0, 22.0])}
    aligned = align_slices_by_anchor_points(coords, labels, anchor_indices, anchor_refs)
    anchor_events, anchor_slice_ids, anchor_ids = generate_stationary_anchor_events_1p1(
        np.asarray([0.0, 1.0]),
        np.asarray([-0.5, 0.5]),
    )
    object_events, object_slice_ids, object_ids = generate_persistent_object_events_1p1(
        np.asarray([0.0, 1.0, 2.0]),
        np.asarray([0.0, 1.0]),
        np.asarray([0.5, -0.25]),
    )
    velocities = finite_difference_velocity_by_object(
        object_events[:, 1],
        object_slice_ids,
        object_ids,
        {0: 0.0, 1: 1.0, 2: 2.0},
    )
    rule = SliceTransportRule(
        name="identity",
        description="identity",
        slice_labels=labels,
        scale_by_slice={0: 1.0, 1: 1.0},
        shift_by_slice={0: 0.0, 1: 0.0},
        reflection_by_slice={0: 1.0, 1: 1.0},
    )
    _ = rule
    return [
        _row(
            "same_position_without_transport",
            not same_position_without_transport().defined,
            0.0,
        ),
        _row(
            "velocity_without_transport",
            not velocity_without_transport().defined,
            0.0,
        ),
        _row("constraint_components", len(components) == 2, len(components)),
        _row(
            "component_labels",
            np.array_equal(component_labels, [0, 0, 0, 1, 1, 1]),
            component_labels[0],
        ),
        _row(
            "cross_slice_fraction_same",
            constraint_cross_slice_fraction(constraints, labels) == 0.0,
            0.0,
        ),
        _row(
            "constraints_for_slice",
            local_constraints.shape[0] == 1,
            local_constraints.shape[0],
        ),
        _row(
            "remap_constraints",
            np.array_equal(remapped, [[0, 1, 1, 2]]) and mapping[2] == 2,
            remapped.shape[0],
        ),
        _row(
            "random_gauge_shape",
            random_coords.shape == coords.shape,
            random_coords.size,
        ),
        _row(
            "affine_fit",
            abs(scale - 2.0) < 1e-12
            and abs(offset - 1.0) < 1e-12
            and rmse < 1e-12,
            rmse,
        ),
        _row(
            "anchor_alignment",
            np.allclose(aligned[[0, 1, 3, 4]], [10.0, 12.0, 20.0, 22.0]),
            aligned[0],
        ),
        _row(
            "anchor_generation",
            anchor_events.shape == (4, 2)
            and anchor_slice_ids.size == 4
            and anchor_ids.size == 4,
            anchor_events.shape[0],
        ),
        _row(
            "persistent_generation",
            object_events.shape == (6, 2) and set(object_ids) == {0, 1},
            object_events.shape[0],
        ),
        _row(
            "finite_velocity",
            abs(velocities[0] - 0.5) < 1e-12
            and abs(velocities[1] + 0.25) < 1e-12,
            velocities[0],
        ),
    ]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_path: Path = DEFAULT_OUTPUT,
) -> Path:
    """Write exact sanity rows."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote cross-slice transport exact sanity results: {output_path}")


if __name__ == "__main__":
    main()
