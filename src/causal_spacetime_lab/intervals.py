"""Alexandrov interval helpers for finite causal orders."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _as_square_bool_matrix(matrix: ArrayLike) -> NDArray[np.bool_]:
    array = np.asarray(matrix, dtype=bool)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError("C must be a square boolean matrix")
    return array


def alexandrov_interval_indices(C: ArrayLike, i: int, j: int) -> NDArray[np.int_]:
    """Return indices ``k`` such that ``i`` precedes ``k`` and ``k`` precedes ``j``."""

    causal_matrix = _as_square_bool_matrix(C)
    n = causal_matrix.shape[0]
    i_index = int(i)
    j_index = int(j)
    if i_index < 0 or i_index >= n:
        raise IndexError("i index out of range")
    if j_index < 0 or j_index >= n:
        raise IndexError("j index out of range")
    return np.flatnonzero(causal_matrix[i_index] & causal_matrix[:, j_index])

