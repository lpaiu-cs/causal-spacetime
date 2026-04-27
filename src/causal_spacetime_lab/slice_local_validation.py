"""Validation helpers for slice-local embeddings."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from causal_spacetime_lab.ordinal_embedding import (
    embedding_distance_order_error,
    procrustes_align,
)


def validate_slice_local_embeddings_against_true_positions(
    slice_embeddings: dict[int, dict[str, object]],
    true_positions: ArrayLike,
) -> list[dict[str, float]]:
    """Validate slice-local embeddings against hidden positions per slice."""

    truth = np.asarray(true_positions, dtype=float)
    if truth.ndim == 2 and truth.shape[1] == 1:
        truth = truth[:, 0]
    if truth.ndim != 1:
        raise ValueError("true_positions must be one-dimensional or shape (n, 1)")
    rows: list[dict[str, float]] = []
    for slice_label, payload in sorted(slice_embeddings.items()):
        embedding = np.asarray(payload["embedding"], dtype=float)
        if embedding.ndim == 1:
            embedding = embedding[:, None]
        indices = np.asarray(payload["global_indices"], dtype=int)
        local_truth = truth[indices, None]
        aligned, diagnostics = procrustes_align(embedding, local_truth)
        order_error = embedding_distance_order_error(
            aligned,
            local_truth,
            seed=17 + int(slice_label),
        )
        constraints = np.asarray(payload["constraints"], dtype=int)
        rows.append(
            {
                "slice_label": float(slice_label),
                "point_count": float(indices.size),
                "constraint_count": float(constraints.shape[0]),
                "local_rmse": diagnostics["rmse"],
                "local_order_error": order_error,
            }
        )
    return rows


def summarize_slice_local_validation(
    rows: list[dict[str, float]],
) -> dict[str, float]:
    """Summarize slice-local validation rows."""

    if not rows:
        return {
            "slice_count": 0.0,
            "mean_local_rmse": float("nan"),
            "median_local_rmse": float("nan"),
            "mean_local_order_error": float("nan"),
            "median_local_order_error": float("nan"),
            "total_points_validated": 0.0,
        }
    rmse = np.asarray([row["local_rmse"] for row in rows], dtype=float)
    order = np.asarray([row["local_order_error"] for row in rows], dtype=float)
    points = np.asarray([row["point_count"] for row in rows], dtype=float)
    return {
        "slice_count": float(len(rows)),
        "mean_local_rmse": float(np.nanmean(rmse)),
        "median_local_rmse": float(np.nanmedian(rmse)),
        "mean_local_order_error": float(np.nanmean(order)),
        "median_local_order_error": float(np.nanmedian(order)),
        "total_points_validated": float(np.nansum(points)),
    }
