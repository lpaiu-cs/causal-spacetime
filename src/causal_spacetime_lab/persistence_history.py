"""Utilities for histories under inferred persistence hypotheses."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.persistence_matching import (
    best_persistence_matchings_by_relational_order,
)


def infer_adjacent_matchings(
    positions_by_slice_unlabeled: dict[int, ArrayLike],
    *,
    max_n: int = 8,
    tolerance: float = 0.0,
    top_k: int = 3,
) -> dict[tuple[int, int], list[dict[str, object]]]:
    """Infer top persistence hypotheses between consecutive slice labels."""

    labels = sorted(positions_by_slice_unlabeled)
    result: dict[tuple[int, int], list[dict[str, object]]] = {}
    for left, right in zip(labels[:-1], labels[1:], strict=False):
        result[(int(left), int(right))] = (
            best_persistence_matchings_by_relational_order(
                positions_by_slice_unlabeled[left],
                positions_by_slice_unlabeled[right],
                max_n=max_n,
                tolerance=tolerance,
                top_k=top_k,
            )
        )
    return result


def build_identity_tracks_from_adjacent_matchings(
    matchings: dict[tuple[int, int], ArrayLike],
    initial_object_ids: ArrayLike,
) -> dict[int, dict[int, int]]:
    """Build object-id -> local-index tracks from adjacent permutations."""

    object_ids = np.asarray(initial_object_ids, dtype=int)
    if object_ids.ndim != 1:
        raise ValueError("initial_object_ids must be one-dimensional")
    if not matchings:
        return {0: {int(object_id): i for i, object_id in enumerate(object_ids)}}
    first_slice = min(left for left, _ in matchings)
    tracks: dict[int, dict[int, int]] = {
        int(first_slice): {
            int(object_id): int(local_index)
            for local_index, object_id in enumerate(object_ids)
        }
    }
    for left, right in sorted(matchings):
        if left not in tracks:
            raise ValueError("adjacent matchings must form a connected chain")
        permutation = np.asarray(matchings[(left, right)], dtype=int)
        tracks[int(right)] = {
            object_id: int(permutation[local_index])
            for object_id, local_index in tracks[int(left)].items()
        }
    return tracks


def relational_history_from_inferred_tracks(
    positions_by_slice_unlabeled: dict[int, ArrayLike],
    tracks: dict[int, dict[int, int]],
) -> dict[int, dict[int, float]]:
    """Convert unlabeled positions and tracks to persistent-object histories."""

    output: dict[int, dict[int, float]] = {}
    for slice_label, object_map in sorted(tracks.items()):
        positions = np.asarray(
            positions_by_slice_unlabeled[int(slice_label)],
            dtype=float,
        )
        output[int(slice_label)] = {}
        for object_id, local_index in object_map.items():
            if local_index < 0 or local_index >= positions.size:
                raise IndexError("track local index is out of range")
            output[int(slice_label)][int(object_id)] = float(
                positions[int(local_index)]
            )
    return output


def track_consistency_error(
    inferred_tracks: dict[int, dict[int, int]],
    true_tracks: dict[int, dict[int, int]],
) -> float:
    """Return the fraction of mismatched object-to-local-index assignments."""

    total = 0
    mismatches = 0
    for slice_label in sorted(set(inferred_tracks) & set(true_tracks)):
        inferred = inferred_tracks[slice_label]
        truth = true_tracks[slice_label]
        for object_id in sorted(set(inferred) & set(truth)):
            total += 1
            mismatches += int(inferred[object_id] != truth[object_id])
    return float(mismatches / total) if total else 0.0


def adjacent_true_permutation_from_tracks(
    true_tracks: dict[int, dict[int, int]],
    left: int,
    right: int,
) -> NDArray[np.int_]:
    """Return validation permutation from true tracks for adjacent slices."""

    left_track = true_tracks[int(left)]
    right_track = true_tracks[int(right)]
    n = len(left_track)
    permutation = np.empty(n, dtype=int)
    for object_id, left_index in left_track.items():
        permutation[int(left_index)] = int(right_track[int(object_id)])
    return permutation
