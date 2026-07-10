"""Benincasa-Dowker 2D causal-set action, raw and smeared (nonlocal).

Conventions verified against Dowker & Glaser 2013 (arXiv:1305.2588):

- ``n(x, y)`` is the number of elements strictly between a related pair
  (open-interval cardinality), i.e. the inclusive interval size minus 2.
- Raw 2D action (eq. 1 layer coefficients ``C^(2) = (1, -2, 1)``,
  ``alpha_2 = -2``, ``beta_2 = 4``; action from ``R(x) = -2 l^2 B1(x)``)::

      S / hbar = 2 (N - 2 n_0 + 4 n_1 - 2 n_2)

  where ``n_k`` counts related pairs with exactly ``k`` elements between.
  Self-check: on sprinklings of a flat 2D causal diamond the action
  fluctuates around 0 (the Einstein-Hilbert value for flat space).
- Smeared (nonlocal) family (eqs. 25-26), nonlocality parameter
  ``eps = (l / xi)^2 in (0, 1]``::

      S_eps = 2 eps N - 4 eps^2 sum_{x < y} f2(n(x, y), eps)
      f2(n, eps) = (1-eps)^n [1 - 2 eps n / (1-eps)
                              + eps^2 n (n-1) / (2 (1-eps)^2)]

  ``eps = 1`` reduces exactly to the raw action. Smearing damps the
  per-realization fluctuations on sprinklings (Sorkin damping).

Known structural fact (established during E1/E2 exploration and used by the
P4 prereg): the complete bipartite order on N elements has
``S_eps = eps N (2 - eps N)`` exactly, which lies far below the sprinkled
value whenever ``eps N >> 2`` -- the action alone does not select manifoldlike
orders in the unrestricted poset ensemble.
"""

from __future__ import annotations

import numpy as np


def interval_sizes(causal_matrix: np.ndarray) -> np.ndarray:
    """Open-interval cardinality n(x, y) for every related pair.

    Returns a 1D array with one entry per related pair (order unspecified).
    """
    strict = np.array(causal_matrix, dtype=bool)
    np.fill_diagonal(strict, False)
    between = strict.astype(np.int32) @ strict.astype(np.int32)
    return between[strict]


def abundances(causal_matrix: np.ndarray, kmax: int = 2) -> tuple[int, list[int]]:
    """Return (R, [n_0 .. n_kmax]): relation count and interval abundances."""
    sizes = interval_sizes(causal_matrix)
    return int(sizes.size), [int(np.sum(sizes == k)) for k in range(kmax + 1)]

def bd_action_2d(causal_matrix: np.ndarray) -> float:
    """Raw 2D Benincasa-Dowker action S/hbar = 2(N - 2 n0 + 4 n1 - 2 n2)."""
    n_elements = causal_matrix.shape[0]
    _, (n0, n1, n2) = abundances(causal_matrix, kmax=2)
    return 2.0 * (n_elements - 2 * n0 + 4 * n1 - 2 * n2)


def smearing_f2(sizes: np.ndarray, eps: float) -> np.ndarray:
    """2D smearing function f2(n, eps) of Dowker-Glaser eq. 26.

    At ``eps = 1`` returns the raw layer coefficients (1, -2, 1) on
    n = 0, 1, 2 and 0 beyond, the exact raw-action limit.
    """
    if eps >= 1.0:
        out = np.zeros(sizes.shape, dtype=float)
        out[sizes == 0] = 1.0
        out[sizes == 1] = -2.0
        out[sizes == 2] = 1.0
        return out
    ratio = eps / (1.0 - eps)
    n = sizes.astype(float)
    return (1.0 - eps) ** n * (
        1.0 - 2.0 * ratio * n + 0.5 * ratio * ratio * n * (n - 1.0)
    )


def smeared_action_2d(causal_matrix: np.ndarray, eps: float) -> float:
    """Smeared 2D BD action S_eps = 2 eps N - 4 eps^2 sum f2(n, eps)."""
    sizes = interval_sizes(causal_matrix)
    n_elements = causal_matrix.shape[0]
    return 2.0 * eps * n_elements - 4.0 * eps * eps * float(
        np.sum(smearing_f2(sizes, eps))
    )
