"""Exploratory spacelike-distance proxies for finite causal sets."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class EnclosingInterval:
    """A candidate Alexandrov interval enclosing a spacelike pair."""

    past_index: int
    future_index: int
    count: int


def _as_square_bool_matrix(matrix: ArrayLike) -> NDArray[np.bool_]:
    array = np.asarray(matrix, dtype=bool)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError("C must be a square boolean matrix")
    return array


def _validate_index(index: int, n: int, name: str) -> int:
    integer_index = int(index)
    if integer_index < 0 or integer_index >= n:
        raise IndexError(f"{name} index out of range")
    return integer_index


def is_spacelike_pair(C: ArrayLike, i: int, j: int) -> bool:
    """Return true when neither event causally precedes the other."""

    causal_matrix = _as_square_bool_matrix(C)
    n = causal_matrix.shape[0]
    i_index = _validate_index(i, n, "i")
    j_index = _validate_index(j, n, "j")
    if i_index == j_index:
        return False
    return not bool(causal_matrix[i_index, j_index] or causal_matrix[j_index, i_index])


def common_future_indices(C: ArrayLike, i: int, j: int) -> NDArray[np.int_]:
    """Return indices of events in the common causal future of ``i`` and ``j``."""

    causal_matrix = _as_square_bool_matrix(C)
    n = causal_matrix.shape[0]
    i_index = _validate_index(i, n, "i")
    j_index = _validate_index(j, n, "j")
    return np.flatnonzero(causal_matrix[i_index] & causal_matrix[j_index])


def common_past_indices(C: ArrayLike, i: int, j: int) -> NDArray[np.int_]:
    """Return indices of events in the common causal past of ``i`` and ``j``."""

    causal_matrix = _as_square_bool_matrix(C)
    n = causal_matrix.shape[0]
    i_index = _validate_index(i, n, "i")
    j_index = _validate_index(j, n, "j")
    return np.flatnonzero(causal_matrix[:, i_index] & causal_matrix[:, j_index])


def common_future_overlap_count(C: ArrayLike, i: int, j: int) -> int:
    """Return the common-future overlap count for a pair of events."""

    return int(common_future_indices(C, i, j).size)


def common_past_overlap_count(C: ArrayLike, i: int, j: int) -> int:
    """Return the common-past overlap count for a pair of events."""

    return int(common_past_indices(C, i, j).size)


def alexandrov_interval_count_matrix(C: ArrayLike) -> NDArray[np.int64]:
    """Return all pairwise Alexandrov interval cardinalities.

    Entry ``[i, j]`` counts events ``k`` such that ``i`` precedes ``k`` and
    ``k`` precedes ``j``.
    """

    causal_matrix = _as_square_bool_matrix(C).astype(np.int64)
    return causal_matrix @ causal_matrix


def minimal_enclosing_alexandrov_interval(
    C: ArrayLike,
    i: int,
    j: int,
    interval_counts: ArrayLike | None = None,
) -> EnclosingInterval | None:
    """Find a minimal sampled Alexandrov interval enclosing a spacelike pair.

    The returned count is the number of sampled events inside the enclosing
    interval, excluding the enclosing endpoints. It includes the queried pair
    when both events lie inside the interval.
    """

    causal_matrix = _as_square_bool_matrix(C)
    n = causal_matrix.shape[0]
    i_index = _validate_index(i, n, "i")
    j_index = _validate_index(j, n, "j")

    past = common_past_indices(causal_matrix, i_index, j_index)
    future = common_future_indices(causal_matrix, i_index, j_index)
    if past.size == 0 or future.size == 0:
        return None

    if interval_counts is None:
        counts = alexandrov_interval_count_matrix(causal_matrix)
    else:
        counts = np.asarray(interval_counts, dtype=np.int64)
        if counts.shape != causal_matrix.shape:
            raise ValueError("interval_counts must have the same shape as C")

    valid_enclosures = causal_matrix[np.ix_(past, future)]
    if not np.any(valid_enclosures):
        return None

    candidate_counts = np.where(
        valid_enclosures,
        counts[np.ix_(past, future)],
        np.iinfo(np.int64).max,
    )
    flat_index = int(np.argmin(candidate_counts))
    past_position, future_position = np.unravel_index(
        flat_index,
        candidate_counts.shape,
    )

    return EnclosingInterval(
        past_index=int(past[past_position]),
        future_index=int(future[future_position]),
        count=int(candidate_counts[past_position, future_position]),
    )


def minimal_enclosing_alexandrov_interval_count(
    C: ArrayLike,
    i: int,
    j: int,
    interval_counts: ArrayLike | None = None,
) -> int | None:
    """Return only the minimal enclosing interval count for a pair."""

    interval = minimal_enclosing_alexandrov_interval(C, i, j, interval_counts)
    if interval is None:
        return None
    return interval.count
