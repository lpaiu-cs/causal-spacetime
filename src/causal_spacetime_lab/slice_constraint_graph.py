"""Constraint graph diagnostics for slice-local order data."""

from __future__ import annotations

from collections import defaultdict

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _as_constraints(constraints: ArrayLike) -> NDArray[np.int_]:
    array = np.asarray(constraints, dtype=int)
    if array.ndim != 2 or array.shape[1] != 4:
        raise ValueError("constraints must have shape (m, 4)")
    if np.any(array < 0):
        raise IndexError("constraints contain negative indices")
    return array


class _UnionFind:
    def __init__(self, n: int) -> None:
        self.parent = np.arange(n, dtype=int)
        self.touched = np.zeros(n, dtype=bool)

    def find(self, item: int) -> int:
        root = int(item)
        while self.parent[root] != root:
            root = int(self.parent[root])
        while self.parent[item] != item:
            parent = int(self.parent[item])
            self.parent[item] = root
            item = parent
        return root

    def union(self, a: int, b: int) -> None:
        root_a = self.find(a)
        root_b = self.find(b)
        if root_a != root_b:
            self.parent[root_b] = root_a


def constraint_point_graph_components(
    n_points: int,
    constraints: ArrayLike,
) -> list[NDArray[np.int_]]:
    """Return connected point components induced by quadruplet constraints."""

    n = int(n_points)
    if n < 0:
        raise ValueError("n_points must be nonnegative")
    constraint_array = _as_constraints(constraints)
    if constraint_array.size and np.any(constraint_array >= n):
        raise IndexError("constraints contain an index out of range")
    graph = _UnionFind(n)
    for row in constraint_array:
        unique = np.unique(row)
        graph.touched[unique] = True
        for value in unique[1:]:
            graph.union(int(unique[0]), int(value))

    groups: dict[int, list[int]] = defaultdict(list)
    for index in np.flatnonzero(graph.touched):
        groups[graph.find(int(index))].append(int(index))
    return [np.asarray(values, dtype=int) for values in groups.values()]


def component_labels_from_constraints(
    n_points: int,
    constraints: ArrayLike,
) -> NDArray[np.int_]:
    """Return component labels for constrained points and ``-1`` for isolated."""

    labels = np.full(int(n_points), -1, dtype=int)
    for label, component in enumerate(
        constraint_point_graph_components(n_points, constraints)
    ):
        labels[component] = label
    return labels


def slice_component_summary(
    slice_labels: ArrayLike,
    constraints: ArrayLike,
    n_points: int,
) -> dict[str, float]:
    """Summarize constraint components and slice coverage."""

    labels = np.asarray(slice_labels, dtype=int)
    if labels.ndim != 1 or labels.shape[0] != int(n_points):
        raise ValueError("slice_labels must have shape (n_points,)")
    components = constraint_point_graph_components(n_points, constraints)
    component_labels = component_labels_from_constraints(n_points, constraints)
    nonnegative_slices = sorted(set(int(value) for value in labels if value >= 0))
    components_per_slice: list[int] = []
    points_per_slice: list[int] = []
    for slice_label in nonnegative_slices:
        in_slice = labels == slice_label
        points_per_slice.append(int(np.count_nonzero(in_slice)))
        components_here = {
            int(value) for value in component_labels[in_slice] if value >= 0
        }
        components_per_slice.append(len(components_here))
    largest = max((component.size for component in components), default=0)
    return {
        "point_count": float(n_points),
        "component_count": float(len(components)),
        "largest_component_size": float(largest),
        "isolated_count": float(np.count_nonzero(component_labels < 0)),
        "slice_count": float(len(nonnegative_slices)),
        "mean_points_per_slice": float(np.mean(points_per_slice))
        if points_per_slice
        else float("nan"),
        "mean_components_per_slice": float(np.mean(components_per_slice))
        if components_per_slice
        else float("nan"),
    }


def constraint_cross_slice_fraction(
    constraints: ArrayLike,
    slice_labels: ArrayLike,
) -> float:
    """Return fraction of constraints whose point IDs span multiple slices."""

    constraint_array = _as_constraints(constraints)
    labels = np.asarray(slice_labels, dtype=int)
    if labels.ndim != 1:
        raise ValueError("slice_labels must be one-dimensional")
    if constraint_array.size == 0:
        return float("nan")
    if np.any(constraint_array >= labels.size):
        raise IndexError("constraints contain an index out of range")
    cross_slice = 0
    for row in constraint_array:
        row_labels = labels[row]
        nonnegative = row_labels[row_labels >= 0]
        if np.unique(nonnegative).size > 1:
            cross_slice += 1
    return float(cross_slice / constraint_array.shape[0])
