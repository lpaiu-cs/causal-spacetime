"""Slice-local ordinal embedding utilities."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.ordinal_embedding import (
    fit_ordinal_embedding_gradient_descent,
)


def _as_constraints(constraints: ArrayLike) -> NDArray[np.int_]:
    array = np.asarray(constraints, dtype=int)
    if array.ndim != 2 or array.shape[1] != 4:
        raise ValueError("constraints must have shape (m, 4)")
    return array


def constraints_for_slice(
    constraints: ArrayLike,
    slice_labels: ArrayLike,
    slice_label: int,
) -> NDArray[np.int_]:
    """Return constraints whose four point IDs all have the requested slice."""

    constraint_array = _as_constraints(constraints)
    labels = np.asarray(slice_labels, dtype=int)
    if labels.ndim != 1:
        raise ValueError("slice_labels must be one-dimensional")
    if constraint_array.size == 0:
        return np.empty((0, 4), dtype=int)
    if np.any((constraint_array < 0) | (constraint_array >= labels.size)):
        raise IndexError("constraints contain an index out of range")
    keep = np.all(labels[constraint_array] == int(slice_label), axis=1)
    return constraint_array[keep].astype(int, copy=False)


def remap_constraints_to_local_indices(
    constraints: ArrayLike,
    global_indices: ArrayLike,
) -> tuple[NDArray[np.int_], dict[int, int]]:
    """Map global point IDs in constraints to local indices."""

    constraint_array = _as_constraints(constraints)
    indices = np.asarray(global_indices, dtype=int)
    if indices.ndim != 1:
        raise ValueError("global_indices must be one-dimensional")
    mapping = {int(global_id): local for local, global_id in enumerate(indices)}
    remapped = np.empty_like(constraint_array)
    for row_index, row in enumerate(constraint_array):
        for col_index, value in enumerate(row):
            try:
                remapped[row_index, col_index] = mapping[int(value)]
            except KeyError as exc:
                msg = "constraint contains point outside global_indices"
                raise KeyError(msg) from exc
    return remapped, mapping


def fit_slice_local_ordinal_embeddings(
    n_points: int,
    slice_labels: ArrayLike,
    constraints: ArrayLike,
    embedding_dim: int = 1,
    min_points_per_slice: int = 4,
    min_constraints_per_slice: int = 20,
    steps: int = 600,
    restarts: int = 2,
    learning_rate: float = 0.1,
    seed: int | None = None,
) -> dict[int, dict[str, object]]:
    """Fit independent ordinal embeddings within each slice."""

    labels = np.asarray(slice_labels, dtype=int)
    if labels.ndim != 1 or labels.shape[0] != int(n_points):
        raise ValueError("slice_labels must have shape (n_points,)")
    constraint_array = _as_constraints(constraints)
    result: dict[int, dict[str, object]] = {}
    for slice_label in sorted(set(int(value) for value in labels if value >= 0)):
        global_indices = np.flatnonzero(labels == slice_label)
        if global_indices.size < int(min_points_per_slice):
            continue
        slice_constraints = constraints_for_slice(
            constraint_array,
            labels,
            slice_label,
        )
        if slice_constraints.shape[0] < int(min_constraints_per_slice):
            continue
        local_constraints, mapping = remap_constraints_to_local_indices(
            slice_constraints,
            global_indices,
        )
        embedding, diagnostics = fit_ordinal_embedding_gradient_descent(
            global_indices.size,
            int(embedding_dim),
            local_constraints,
            steps=steps,
            restarts=restarts,
            learning_rate=learning_rate,
            seed=None if seed is None else seed + 10_000 * slice_label,
            batch_size=min(2048, max(1, local_constraints.shape[0])),
        )
        result[slice_label] = {
            "embedding": embedding,
            "global_indices": global_indices,
            "constraints": local_constraints,
            "global_constraints": slice_constraints,
            "mapping": mapping,
            "diagnostics": diagnostics,
        }
    return result


def assemble_slice_embeddings_with_nan(
    n_points: int,
    slice_embeddings: dict[int, dict[str, object]],
) -> NDArray[np.float64]:
    """Assemble slice-local 1D embeddings into one vector with NaNs for gaps."""

    assembled = np.full((int(n_points), 1), np.nan, dtype=float)
    for payload in slice_embeddings.values():
        embedding = np.asarray(payload["embedding"], dtype=float)
        indices = np.asarray(payload["global_indices"], dtype=int)
        if embedding.ndim == 1:
            embedding = embedding[:, None]
        if embedding.shape[0] != indices.size:
            raise ValueError("embedding and global_indices size mismatch")
        assembled[indices, 0] = embedding[:, 0]
    return assembled
