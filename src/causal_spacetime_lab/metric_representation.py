"""Finite diagnostics for metric representations of distance-order data."""

from __future__ import annotations

from collections import deque

import numpy as np
from numpy.typing import ArrayLike, NDArray


def unordered_pair_indices(n: int) -> NDArray[np.int_]:
    """Return all unordered index pairs ``(i, j)`` with ``i < j``."""

    count = int(n)
    if count < 0:
        raise ValueError("n must be nonnegative")
    pairs = [(i, j) for i in range(count) for j in range(i + 1, count)]
    return np.asarray(pairs, dtype=int).reshape((-1, 2))


def distance_order_constraints_from_values(
    values: ArrayLike,
    tolerance: float = 0.0,
) -> NDArray[np.int_]:
    """Return strict ordinal constraints between pair-distance value indices."""

    distances = np.asarray(values, dtype=float)
    if distances.ndim != 1:
        raise ValueError("values must be one-dimensional")
    tol = float(tolerance)
    if tol < 0.0:
        raise ValueError("tolerance must be nonnegative")
    constraints: list[tuple[int, int]] = []
    for left in range(distances.size):
        for right in range(left + 1, distances.size):
            if distances[left] + tol < distances[right]:
                constraints.append((left, right))
            elif distances[right] + tol < distances[left]:
                constraints.append((right, left))
    return np.asarray(constraints, dtype=int).reshape((-1, 2))


def _validate_ordered_pairs(
    num_items: int,
    ordered_pairs: ArrayLike,
) -> NDArray[np.int_]:
    count = int(num_items)
    if count < 0:
        raise ValueError("num_items must be nonnegative")
    pairs = np.asarray(ordered_pairs, dtype=int)
    if pairs.size == 0:
        return np.empty((0, 2), dtype=int)
    if pairs.ndim != 2 or pairs.shape[1] != 2:
        raise ValueError("ordered_pairs must have shape (m, 2)")
    if np.any((pairs < 0) | (pairs >= count)):
        raise IndexError("ordered_pairs contain an index out of range")
    return pairs


def _topological_order(
    num_items: int,
    ordered_pairs: NDArray[np.int_],
) -> tuple[list[int], list[list[int]]]:
    outgoing: list[list[int]] = [[] for _ in range(num_items)]
    indegree = np.zeros(num_items, dtype=int)
    for source, target in ordered_pairs:
        outgoing[int(source)].append(int(target))
        indegree[int(target)] += 1

    queue: deque[int] = deque(int(i) for i in np.flatnonzero(indegree == 0))
    order: list[int] = []
    while queue:
        node = queue.popleft()
        order.append(node)
        for target in outgoing[node]:
            indegree[target] -= 1
            if indegree[target] == 0:
                queue.append(target)
    return order, outgoing


def has_order_cycle(num_items: int, ordered_pairs: ArrayLike) -> bool:
    """Return whether strict order constraints contain a directed cycle."""

    count = int(num_items)
    pairs = _validate_ordered_pairs(count, ordered_pairs)
    order, _ = _topological_order(count, pairs)
    return len(order) != count


def topological_rank_representation(
    num_items: int,
    ordered_pairs: ArrayLike,
) -> NDArray[np.float64]:
    """Assign ordinal ranks satisfying acyclic strict order constraints."""

    count = int(num_items)
    pairs = _validate_ordered_pairs(count, ordered_pairs)
    order, outgoing = _topological_order(count, pairs)
    if len(order) != count:
        raise ValueError("ordered_pairs contain a cycle")
    ranks = np.zeros(count, dtype=np.float64)
    for node in order:
        for target in outgoing[node]:
            ranks[target] = max(ranks[target], ranks[node] + 1.0)
    return ranks


def _distance_matrix_from_unordered(
    distances: NDArray[np.float64],
    n_points: int,
) -> NDArray[np.float64]:
    expected = n_points * (n_points - 1) // 2
    if distances.ndim != 1 or distances.size != expected:
        raise ValueError("distances length must be n_points * (n_points - 1) / 2")
    matrix = np.zeros((n_points, n_points), dtype=np.float64)
    pairs = unordered_pair_indices(n_points)
    matrix[pairs[:, 0], pairs[:, 1]] = distances
    matrix[pairs[:, 1], pairs[:, 0]] = distances
    return matrix


def triangle_inequality_violations(
    distances: ArrayLike,
    n_points: int,
) -> int:
    """Count triangle inequality violations in unordered-pair distances."""

    count = int(n_points)
    if count < 0:
        raise ValueError("n_points must be nonnegative")
    values = np.asarray(distances, dtype=float)
    matrix = _distance_matrix_from_unordered(values, count)
    violations = 0
    tol = 1e-12
    for i in range(count):
        for j in range(i + 1, count):
            for k in range(j + 1, count):
                dij = matrix[i, j]
                dik = matrix[i, k]
                djk = matrix[j, k]
                violations += int(dij > dik + djk + tol)
                violations += int(dik > dij + djk + tol)
                violations += int(djk > dij + dik + tol)
    return violations


def euclidean_gram_from_distance_matrix(
    distance_matrix: ArrayLike,
) -> NDArray[np.float64]:
    """Return the classical double-centered Gram matrix for distances."""

    distances = np.asarray(distance_matrix, dtype=float)
    if distances.ndim != 2 or distances.shape[0] != distances.shape[1]:
        raise ValueError("distance_matrix must be square")
    n = distances.shape[0]
    centering = np.eye(n) - np.ones((n, n), dtype=np.float64) / n
    return -0.5 * centering @ (distances * distances) @ centering


def euclidean_embedding_diagnostics(
    distance_matrix: ArrayLike,
    tol: float = 1e-8,
) -> dict[str, float | bool]:
    """Return lightweight diagnostics for a Euclidean distance candidate."""

    distances = np.asarray(distance_matrix, dtype=float)
    if distances.ndim != 2 or distances.shape[0] != distances.shape[1]:
        raise ValueError("distance_matrix must be square")
    tolerance = float(tol)
    if tolerance < 0.0:
        raise ValueError("tol must be nonnegative")
    symmetric = bool(np.allclose(distances, distances.T, atol=tolerance, rtol=0.0))
    zero_diagonal = bool(np.allclose(np.diag(distances), 0.0, atol=tolerance))
    gram = euclidean_gram_from_distance_matrix(distances)
    eigenvalues = np.linalg.eigvalsh(gram) if symmetric else np.asarray([np.nan])
    min_eigenvalue = float(np.min(eigenvalues))
    positive_count = int(np.count_nonzero(eigenvalues > tolerance))
    candidate = symmetric and zero_diagonal and min_eigenvalue >= -tolerance
    return {
        "is_symmetric": symmetric,
        "has_zero_diagonal": zero_diagonal,
        "min_gram_eigenvalue": min_eigenvalue,
        "positive_eigenvalue_count": float(positive_count),
        "is_euclidean_candidate": bool(candidate),
    }
