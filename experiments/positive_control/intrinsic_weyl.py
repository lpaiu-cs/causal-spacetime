"""Development estimator: local vacuum Weyl tensor from interval cardinalities.

No spacetime coordinates, profile, or paired flat sample enter this module.
Four-dimensional vacuum and uniform physical-volume sampling are assumptions.
Finite-sample accuracy and automatic landmark selection are not established.
"""

import itertools
from functools import lru_cache

import numpy as np
from scipy.linalg import null_space

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
PAIRS = list(itertools.combinations(range(4), 2))


def weyl_basis():
    """Ten algebraic Weyl tensors: five electric and five magnetic components."""
    symmetric = []
    for diagonal in ((1, -1, 0), (0, 1, -1)):
        symmetric.append(np.diag(diagonal))
    for i, j in itertools.combinations(range(3), 2):
        value = np.zeros((3, 3))
        value[i, j] = value[j, i] = 1
        symmetric.append(value)
    epsilon = np.zeros((3, 3, 3))
    for i, j, k in itertools.permutations(range(3)):
        epsilon[i, j, k] = (i - j) * (j - k) * (k - i) / 2
    result = []
    for channel in range(2):
        for value in symmetric:
            electric = value if channel == 0 else np.zeros((3, 3))
            magnetic = value if channel == 1 else np.zeros((3, 3))
            matrix = np.zeros((6, 6))
            for p, (i, j) in enumerate(PAIRS):
                for q, (k, ell) in enumerate(PAIRS):
                    if i == 0 and k == 0:
                        matrix[p, q] = electric[j - 1, ell - 1]
                    elif i == 0:
                        matrix[p, q] = epsilon[k - 1, ell - 1] @ magnetic[j - 1]
                    elif k == 0:
                        matrix[p, q] = matrix[q, p]
                    else:
                        matrix[p, q] = (
                            (i == k) * electric[j - 1, ell - 1]
                            - (i == ell) * electric[j - 1, k - 1]
                            - (j == k) * electric[i - 1, ell - 1]
                            + (j == ell) * electric[i - 1, k - 1]
                        )
            tensor = np.zeros((4, 4, 4, 4))
            for p, (i, j) in enumerate(PAIRS):
                for q, (k, ell) in enumerate(PAIRS):
                    tensor[i, j, k, ell] = tensor[j, i, ell, k] = matrix[p, q]
                    tensor[j, i, k, ell] = tensor[i, j, ell, k] = -matrix[p, q]
            result.append(tensor)
    return np.array(result)


WEYL = weyl_basis()


@lru_cache(maxsize=1)
def gradient_basis():
    """24 vacuum Weyl gradients, in derivative-major (4,10) component order."""
    rows = []
    for a, b, c in itertools.combinations(range(4), 3):
        for d, e in itertools.combinations(range(4), 2):
            row = np.zeros((4, 10))
            row[a] += WEYL[:, b, c, d, e]
            row[b] += WEYL[:, c, a, d, e]
            row[c] += WEYL[:, a, b, d, e]
            rows.append(row.ravel())
    constraint = np.array(rows)
    basis = null_space(constraint)
    if basis.shape != (40, 24) or np.max(abs(constraint @ basis)) >= 1e-12:
        raise ArithmeticError("Invalid differential Bianchi null space")
    basis.setflags(write=False)
    return basis


def quartic(x, y):
    """Return C(x,y,x,y) for each tensor in the fixed algebraic basis."""
    return np.einsum("a,b,c,d,kabcd->k", x, y, x, y, WEYL)


def residual_design(references, coefficients):
    """Linear curvature response after eliminating the unknown flat embedding."""
    points = coefficients @ references
    rr = np.array([[quartic(x, y) for y in references] for x in references])
    rp = np.array([[quartic(x, y) for x in references] for y in points])
    pairs = list(itertools.combinations_with_replacement(range(len(points)), 2))
    rows = []
    for i, j in pairs:
        rows.append(
            (
                quartic(points[i], points[j])
                - coefficients[j] @ rp[i]
                - coefficients[i] @ rp[j]
                + np.einsum("a,b,abk->k", coefficients[i], coefficients[j], rr)
            )
            / 6
        )
    return np.array(rows), pairs


def reference_frame(gram):
    """Lorentzian embedding, with time axis along the first reference event."""
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    if not (eigenvalues[0] < 0 < eigenvalues[1]):
        raise ValueError("Reference interval Gram matrix lacks Lorentzian signature")
    coordinates = eigenvectors * np.sqrt(np.abs(eigenvalues))
    time = coordinates[0] / np.sqrt(-gram[0, 0])
    basis = [time]
    for value in coordinates[1:]:
        residual = value + (value @ ETA @ time) * time
        for spatial in basis[1:]:
            residual -= (residual @ ETA @ spatial) * spatial
        norm = residual @ ETA @ residual
        if norm <= 1e-14:
            raise ValueError("Reference events do not span three spatial directions")
        basis.append(residual / np.sqrt(norm))
    frame = np.column_stack(basis)
    return coordinates @ np.linalg.inv(frame).T


def fit_from_interval_counts(counts):
    """Fit dimensionless Weyl components in the internally reconstructed frame.

    counts is a symmetric matrix of interval cardinalities (or continuum
    volumes for a separately labelled oracle diagnostic). Zero off-diagonal
    entries are unavailable pair clocks. Event 0 precedes all other events;
    events 1..4 specify a spanning reference frame with all pair clocks known.
    Every probe must have a known pair clock to the origin and four references.
    The returned tensor is scaled by the squared time from event 0 to event 4.
    A low residual is not a finite-sample confidence statement.
    """
    counts = np.asarray(counts, dtype=float)
    if counts.ndim != 2 or counts.shape[0] != counts.shape[1] or len(counts) < 10:
        raise ValueError("At least ten events and a square interval matrix required")
    if not np.all(np.isfinite(counts)) or not np.allclose(counts, counts.T):
        raise ValueError("Counts must be finite and symmetric")
    if np.any(np.diag(counts) != 0) or np.any(counts < 0):
        raise ValueError("Counts must be nonnegative with zero diagonal")
    required = np.ones_like(counts[:5], dtype=bool)
    required[np.arange(5), np.arange(5)] = False
    if np.any(counts[:5][required] <= 0):
        raise ValueError("Positive interval counts required to all five references")
    squared = -np.sqrt(counts / counts[0, 4])
    gram = (squared[0, 1:, None] + squared[0, None, 1:] - squared[1:, 1:]) / 2
    reference_gram = gram[:4, :4]
    references = reference_frame(reference_gram)
    coefficients = np.linalg.solve(reference_gram, gram[:4, 4:]).T
    residual = gram[4:, 4:] - coefficients @ reference_gram @ coefficients.T
    design, pairs = residual_design(references, coefficients)
    available = [i == j or counts[i + 5, j + 5] > 0 for i, j in pairs]
    design = design[available]
    pairs = [pair for pair, keep in zip(pairs, available, strict=True) if keep]
    data = np.array([residual[i, j] for i, j in pairs])
    solution, _, rank, singular = np.linalg.lstsq(design, data, rcond=None)
    if rank != 10:
        raise ValueError(
            f"Weyl tensor unidentifiable for this configuration: rank {rank}"
        )
    tensor = np.einsum("k,kabcd->abcd", solution, WEYL)
    electric = tensor[0, 1:, 0, 1:]
    return {
        "coefficients": solution.tolist(),
        "electric_eigenvalues": np.linalg.eigvalsh(electric).tolist(),
        "electric_norm_squared": float(np.sum(electric**2)),
        "design_rank": int(rank),
        "design_condition": float(singular[0] / singular[-1]),
        "relative_residual": float(
            np.linalg.norm(design @ solution - data) / max(np.linalg.norm(data), 1e-30)
        ),
    }


def fit_from_order(causal, landmarks):
    """Use only the order and selected event IDs, including a reference frame.

    Selection of a statistically adequate, independent landmark set is still
    a research requirement, not supplied by this low-level fitting function.
    """
    causal = np.asarray(causal, dtype=bool)
    ids = np.asarray(landmarks, dtype=int)
    if (
        causal.ndim != 2
        or causal.shape[0] != causal.shape[1]
        or np.any(np.diag(causal))
    ):
        raise ValueError("Expected a strict causal-order matrix")
    if (
        ids.ndim != 1
        or len(set(ids)) != len(ids)
        or np.any(ids < 0)
        or np.any(ids >= len(causal))
    ):
        raise ValueError("Landmark IDs must be distinct and in range")
    if len(ids) < 10 or not np.all(causal[ids[0], ids[1:]]):
        raise ValueError("At least ten landmarks with a common first event required")
    counts = np.zeros((len(ids), len(ids)))
    for i, j in itertools.combinations(range(len(ids)), 2):
        a, b = ids[i], ids[j]
        if causal[b, a]:
            a, b = b, a
        if causal[a, b]:
            counts[i, j] = counts[j, i] = np.count_nonzero(causal[a] & causal[:, b])
    return fit_from_interval_counts(counts)
