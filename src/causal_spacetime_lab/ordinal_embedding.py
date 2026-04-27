"""Ordinal embedding diagnostics for effective metric representations."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _as_coord_array(coords: ArrayLike, name: str = "coords") -> NDArray[np.float64]:
    array = np.asarray(coords, dtype=float)
    if array.ndim != 2:
        raise ValueError(f"{name} must have shape (n, d)")
    return array


def _as_constraint_array(constraints: ArrayLike) -> NDArray[np.int_]:
    array = np.asarray(constraints, dtype=int)
    if array.ndim != 2 or array.shape[1] != 4:
        raise ValueError("constraints must have shape (m, 4)")
    if np.any(array < 0):
        raise IndexError("constraints contain negative indices")
    return array


def squared_distance_matrix(coords: ArrayLike) -> NDArray[np.float64]:
    """Return squared Euclidean distance matrix for coordinates."""

    points = _as_coord_array(coords)
    diffs = points[:, None, :] - points[None, :, :]
    return np.sum(diffs * diffs, axis=2).astype(np.float64, copy=False)


def sample_quadruplet_constraints_from_distance_values(
    distance_matrix: ArrayLike,
    num_constraints: int,
    seed: int | None = None,
    tolerance: float = 0.0,
) -> NDArray[np.int_]:
    """Sample constraints ``(i, j, k, l)`` meaning ``d(i,j) < d(k,l)``."""

    distances = np.asarray(distance_matrix, dtype=float)
    if distances.ndim != 2 or distances.shape[0] != distances.shape[1]:
        raise ValueError("distance_matrix must be square")
    count = int(num_constraints)
    if count < 0:
        raise ValueError("num_constraints must be nonnegative")
    tol = float(tolerance)
    if tol < 0.0:
        raise ValueError("tolerance must be nonnegative")
    if count == 0:
        return np.empty((0, 4), dtype=int)

    n = distances.shape[0]
    rng = np.random.default_rng(seed)
    accepted: list[NDArray[np.int_]] = []
    accepted_count = 0
    attempts = 0
    max_attempts = max(100_000, 100 * count)
    while accepted_count < count and attempts < max_attempts:
        batch = min(max(4 * (count - accepted_count), 256), max_attempts - attempts)
        first = rng.integers(0, n, size=(batch, 2))
        second = rng.integers(0, n, size=(batch, 2))
        valid_pairs = (first[:, 0] != first[:, 1]) & (second[:, 0] != second[:, 1])
        dij = distances[first[:, 0], first[:, 1]]
        dkl = distances[second[:, 0], second[:, 1]]
        left_smaller = dij + tol < dkl
        right_smaller = dkl + tol < dij

        rows = np.empty((batch, 4), dtype=int)
        rows[:, 0:2] = first
        rows[:, 2:4] = second
        swapped = rows.copy()
        swapped[:, 0:2] = second
        swapped[:, 2:4] = first

        selected = np.vstack(
            (
                rows[valid_pairs & left_smaller],
                swapped[valid_pairs & right_smaller],
            )
        )
        if selected.size:
            accepted.append(selected)
            accepted_count += selected.shape[0]
        attempts += batch

    if accepted_count == 0:
        return np.empty((0, 4), dtype=int)
    return np.vstack(accepted)[:count].astype(int, copy=False)


def sample_quadruplet_constraints_from_coords(
    coords: ArrayLike,
    num_constraints: int,
    seed: int | None = None,
    tolerance: float = 0.0,
) -> NDArray[np.int_]:
    """Sample distance-order constraints from hidden validation coordinates."""

    return sample_quadruplet_constraints_from_distance_values(
        squared_distance_matrix(coords),
        num_constraints,
        seed=seed,
        tolerance=tolerance,
    )


def quadruplet_violation_rate(
    coords: ArrayLike,
    constraints: ArrayLike,
    tolerance: float = 0.0,
) -> float:
    """Return fraction of quadruplet constraints violated by coordinates."""

    points = _as_coord_array(coords)
    constraint_array = _as_constraint_array(constraints)
    if constraint_array.size == 0:
        return float("nan")
    if np.any(constraint_array >= points.shape[0]):
        raise IndexError("constraints contain an index out of range")
    tol = float(tolerance)
    if tol < 0.0:
        raise ValueError("tolerance must be nonnegative")
    distances = squared_distance_matrix(points)
    lhs = distances[constraint_array[:, 0], constraint_array[:, 1]]
    rhs = distances[constraint_array[:, 2], constraint_array[:, 3]]
    return float(np.mean(lhs + tol >= rhs))


def quadruplet_hinge_loss(
    coords: ArrayLike,
    constraints: ArrayLike,
    margin: float = 1e-3,
) -> float:
    """Return mean squared quadruplet hinge loss."""

    points = _as_coord_array(coords)
    constraint_array = _as_constraint_array(constraints)
    if constraint_array.size == 0:
        return 0.0
    if np.any(constraint_array >= points.shape[0]):
        raise IndexError("constraints contain an index out of range")
    distances = squared_distance_matrix(points)
    lhs = distances[constraint_array[:, 0], constraint_array[:, 1]]
    rhs = distances[constraint_array[:, 2], constraint_array[:, 3]]
    active = np.maximum(0.0, float(margin) + lhs - rhs)
    return float(np.mean(active * active))


def _normalize_embedding(coords: NDArray[np.float64]) -> NDArray[np.float64]:
    normalized = coords - np.mean(coords, axis=0, keepdims=True)
    rms = float(np.sqrt(np.mean(normalized * normalized)))
    if rms > 0.0:
        normalized = normalized / rms
    return normalized


def _batch_gradient(
    coords: NDArray[np.float64],
    constraints: NDArray[np.int_],
    margin: float,
) -> tuple[NDArray[np.float64], float]:
    if constraints.size == 0:
        return np.zeros_like(coords), 0.0
    i = constraints[:, 0]
    j = constraints[:, 1]
    k = constraints[:, 2]
    ell = constraints[:, 3]
    vij = coords[i] - coords[j]
    vkl = coords[k] - coords[ell]
    dij = np.sum(vij * vij, axis=1)
    dkl = np.sum(vkl * vkl, axis=1)
    hinge = float(margin) + dij - dkl
    active = hinge > 0.0
    if not np.any(active):
        return np.zeros_like(coords), 0.0

    active_hinge = hinge[active]
    active_vij = vij[active]
    active_vkl = vkl[active]
    grad = np.zeros_like(coords)
    scale = 4.0 * active_hinge[:, None] / constraints.shape[0]
    np.add.at(grad, i[active], scale * active_vij)
    np.add.at(grad, j[active], -scale * active_vij)
    np.add.at(grad, k[active], -scale * active_vkl)
    np.add.at(grad, ell[active], scale * active_vkl)
    loss = float(np.mean(active_hinge * active_hinge))
    return grad, loss


def fit_ordinal_embedding_gradient_descent(
    n_points: int,
    embedding_dim: int,
    constraints: ArrayLike,
    *,
    steps: int = 2000,
    learning_rate: float = 0.01,
    margin: float = 1e-3,
    seed: int | None = None,
    restarts: int = 5,
    batch_size: int | None = None,
) -> tuple[NDArray[np.float64], dict[str, float]]:
    """Fit a finite ordinal embedding with simple hinge-loss gradient descent."""

    n = int(n_points)
    dim = int(embedding_dim)
    if n <= 0:
        raise ValueError("n_points must be positive")
    if dim <= 0:
        raise ValueError("embedding_dim must be positive")
    step_count = int(steps)
    if step_count <= 0:
        raise ValueError("steps must be positive")
    lr = float(learning_rate)
    if lr <= 0.0:
        raise ValueError("learning_rate must be positive")
    restart_count = int(restarts)
    if restart_count <= 0:
        raise ValueError("restarts must be positive")

    constraint_array = _as_constraint_array(constraints)
    if constraint_array.size and np.any(constraint_array >= n):
        raise IndexError("constraints contain an index out of range")
    rng = np.random.default_rng(seed)
    batch = (
        max(1, constraint_array.shape[0]) if batch_size is None else int(batch_size)
    )
    if batch <= 0:
        raise ValueError("batch_size must be positive")

    best_coords: NDArray[np.float64] | None = None
    best_loss = float("inf")
    best_violation = float("nan")
    for _ in range(restart_count):
        coords = _normalize_embedding(rng.normal(size=(n, dim)))
        for _step in range(step_count):
            if constraint_array.shape[0] > batch:
                selected = rng.choice(
                    constraint_array.shape[0],
                    size=batch,
                    replace=False,
                )
                batch_constraints = constraint_array[selected]
            else:
                batch_constraints = constraint_array
            grad, _ = _batch_gradient(coords, batch_constraints, margin)
            coords = _normalize_embedding(coords - lr * grad)
        final_loss = quadruplet_hinge_loss(coords, constraint_array, margin=margin)
        final_violation = quadruplet_violation_rate(coords, constraint_array)
        if final_loss < best_loss:
            best_loss = final_loss
            best_violation = final_violation
            best_coords = coords.copy()

    if best_coords is None:
        raise RuntimeError("ordinal embedding failed to produce coordinates")
    diagnostics = {
        "final_loss": float(best_loss),
        "violation_rate": float(best_violation),
        "steps": float(step_count),
        "restarts": float(restart_count),
        "constraint_count": float(constraint_array.shape[0]),
    }
    return best_coords, diagnostics


def procrustes_align(
    estimated: ArrayLike,
    reference: ArrayLike,
    allow_scaling: bool = True,
) -> tuple[NDArray[np.float64], dict[str, float]]:
    """Align estimated coordinates to a reference by orthogonal Procrustes."""

    estimate = _as_coord_array(estimated, "estimated")
    target = _as_coord_array(reference, "reference")
    if estimate.shape != target.shape:
        raise ValueError("estimated and reference must have the same shape")
    estimate_mean = np.mean(estimate, axis=0, keepdims=True)
    target_mean = np.mean(target, axis=0, keepdims=True)
    x = estimate - estimate_mean
    y = target - target_mean
    u, singular_values, vt = np.linalg.svd(x.T @ y, full_matrices=False)
    rotation = u @ vt
    scale = 1.0
    if allow_scaling:
        denom = float(np.sum(x * x))
        if denom > 0.0:
            scale = float(np.sum(singular_values) / denom)
    aligned = scale * x @ rotation + target_mean
    rmse = float(np.sqrt(np.mean((aligned - target) ** 2)))
    translation = target_mean - scale * estimate_mean @ rotation
    diagnostics = {
        "rmse": rmse,
        "scale": float(scale),
        "translation_norm": float(np.linalg.norm(translation)),
    }
    return aligned.astype(np.float64, copy=False), diagnostics


def embedding_distance_order_error(
    estimated: ArrayLike,
    reference: ArrayLike,
    num_pair_comparisons: int = 10_000,
    seed: int | None = None,
) -> float:
    """Sample pair-pair comparisons and return distance-order error."""

    estimate = _as_coord_array(estimated, "estimated")
    target = _as_coord_array(reference, "reference")
    if estimate.shape[0] != target.shape[0]:
        raise ValueError("estimated and reference must have the same point count")
    count = int(num_pair_comparisons)
    if count <= 0:
        raise ValueError("num_pair_comparisons must be positive")
    rng = np.random.default_rng(seed)
    n = target.shape[0]
    first = rng.integers(0, n, size=(count, 2))
    second = rng.integers(0, n, size=(count, 2))
    valid = (first[:, 0] != first[:, 1]) & (second[:, 0] != second[:, 1])
    first = first[valid]
    second = second[valid]
    ref_distances = squared_distance_matrix(target)
    est_distances = squared_distance_matrix(estimate)
    true_delta = ref_distances[first[:, 0], first[:, 1]] - ref_distances[
        second[:, 0],
        second[:, 1],
    ]
    estimated_delta = est_distances[first[:, 0], first[:, 1]] - est_distances[
        second[:, 0],
        second[:, 1],
    ]
    true_sign = np.sign(true_delta)
    estimated_sign = np.sign(estimated_delta)
    comparable = true_sign != 0.0
    if not np.any(comparable):
        return float("nan")
    return float(np.mean(true_sign[comparable] != estimated_sign[comparable]))


def ordinal_embedding_dimension_curve(
    n_points: int,
    constraints: ArrayLike,
    candidate_dims: list[int],
    *,
    steps: int = 1500,
    seed: int | None = None,
    restarts: int = 3,
    learning_rate: float = 0.01,
    batch_size: int | None = None,
) -> list[dict[str, float]]:
    """Fit candidate embedding dimensions and report finite ordinal stress."""

    rows: list[dict[str, float]] = []
    for index, dim in enumerate(candidate_dims):
        embedding, diagnostics = fit_ordinal_embedding_gradient_descent(
            n_points,
            int(dim),
            constraints,
            steps=steps,
            learning_rate=learning_rate,
            seed=None if seed is None else seed + 10_000 * index,
            restarts=restarts,
            batch_size=batch_size,
        )
        rows.append(
            {
                "embedding_dim": float(dim),
                "final_loss": diagnostics["final_loss"],
                "violation_rate": quadruplet_violation_rate(embedding, constraints),
                "constraint_count": diagnostics["constraint_count"],
            }
        )
    return rows
