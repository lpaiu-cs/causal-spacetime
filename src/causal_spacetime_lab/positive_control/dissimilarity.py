"""Profile dissimilarity, measured margin, and leak-free constraint splits.

PC-V1 Section 6: RMS dissimilarity over common reachable columns, margin
derived from the measured resolution, and a pair-level train/held-out split
so no target pair contributes to both sides.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.positive_control.echo_profiles import EchoProfileMatrix


def parallax_profiles(profiles: EchoProfileMatrix) -> NDArray[np.float64]:
    """Return reference-centered (parallax) bracket-width profiles.

    Each target's bracket widths are centered across the references that
    reach it, removing the shared per-target scalar (a global/temporal
    common mode that any time-respecting order injects) and retaining only
    the cross-observer parallax that carries spatial information. Enforcing
    this operationalizes the response-order underdetermination principle: a
    single shared scalar across observers is not a distance structure.
    Unreachable entries are left as NaN.
    """

    delays = profiles.delay_ranks
    reachable = profiles.reachable
    centered = np.full_like(delays, np.nan, dtype=np.float64)
    for row in range(profiles.target_count):
        columns = reachable[row]
        if not np.any(columns):
            continue
        values = delays[row, columns].astype(np.float64)
        centered[row, columns] = values - values.mean()
    return centered


def profile_dissimilarity_matrix(
    profiles: EchoProfileMatrix,
    min_common_columns: int = 4,
    center_references: bool = True,
) -> NDArray[np.float64]:
    """Return the RMS profile dissimilarity matrix (NaN where undefined).

    With ``center_references`` (the fixed default), the dissimilarity is
    computed over reference-centered parallax profiles so that a shared
    per-target scalar cannot masquerade as spatial structure. Set it False
    only for diagnostics that deliberately inspect the raw common mode.
    """

    if min_common_columns < 1:
        raise ValueError("min_common_columns must be positive")
    reachable = profiles.reachable
    delays = (
        parallax_profiles(profiles)
        if center_references
        else profiles.delay_ranks.astype(np.float64)
    )
    n = profiles.target_count
    dissimilarity = np.full((n, n), np.nan, dtype=np.float64)
    np.fill_diagonal(dissimilarity, 0.0)
    for i in range(n - 1):
        common = reachable[i] & reachable[i + 1 :]
        counts = np.sum(common, axis=1)
        diffs = np.where(common, delays[i] - delays[i + 1 :], 0.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            rms = np.sqrt(np.sum(diffs * diffs, axis=1) / counts)
        rms = np.where(counts >= min_common_columns, rms, np.nan)
        dissimilarity[i, i + 1 :] = rms
        dissimilarity[i + 1 :, i] = rms
    return dissimilarity


def margin_from_probe_quantile(
    dissimilarity: NDArray[np.float64],
    quantile: float = 0.25,
    num_probes: int = 20_000,
    seed: int = 0,
) -> float:
    """Derive the constraint margin from measured |D(i,j) - D(k,l)| probes."""

    if not 0.0 < quantile < 1.0:
        raise ValueError("quantile must be in (0, 1)")
    matrix = np.asarray(dissimilarity, dtype=float)
    n = matrix.shape[0]
    rng = np.random.default_rng(seed)
    first = rng.integers(0, n, size=(num_probes, 2))
    second = rng.integers(0, n, size=(num_probes, 2))
    valid = (first[:, 0] != first[:, 1]) & (second[:, 0] != second[:, 1])
    dij = matrix[first[valid, 0], first[valid, 1]]
    dkl = matrix[second[valid, 0], second[valid, 1]]
    gaps = np.abs(dij - dkl)
    gaps = gaps[np.isfinite(gaps) & (gaps > 0.0)]
    if gaps.size == 0:
        raise ValueError("no positive finite dissimilarity gaps to probe")
    return float(np.quantile(gaps, quantile))


def pair_split_matrix(
    n_targets: int,
    train_fraction: float = 0.8,
    seed: int = 0,
) -> NDArray[np.bool_]:
    """Return a symmetric boolean matrix marking train pairs."""

    if not 0.0 < train_fraction < 1.0:
        raise ValueError("train_fraction must be in (0, 1)")
    rng = np.random.default_rng(seed)
    upper = rng.uniform(size=(n_targets, n_targets)) < train_fraction
    train = np.triu(upper, k=1)
    return train | train.T


@dataclass(frozen=True)
class ConstraintSplit:
    """Train and held-out quadruplet constraints with the margin used."""

    train: NDArray[np.int_]
    heldout: NDArray[np.int_]
    margin: float


def _sample_side_constraints(
    dissimilarity: NDArray[np.float64],
    side_pairs: NDArray[np.bool_],
    count: int,
    margin: float,
    rng: np.random.Generator,
) -> NDArray[np.int_]:
    n = dissimilarity.shape[0]
    accepted: list[NDArray[np.int_]] = []
    accepted_count = 0
    attempts = 0
    max_attempts = max(200_000, 200 * count)
    while accepted_count < count and attempts < max_attempts:
        batch = min(max(4 * (count - accepted_count), 512), max_attempts - attempts)
        first = rng.integers(0, n, size=(batch, 2))
        second = rng.integers(0, n, size=(batch, 2))
        distinct = (first[:, 0] != first[:, 1]) & (second[:, 0] != second[:, 1])
        on_side = (
            side_pairs[first[:, 0], first[:, 1]]
            & side_pairs[second[:, 0], second[:, 1]]
        )
        dij = dissimilarity[first[:, 0], first[:, 1]]
        dkl = dissimilarity[second[:, 0], second[:, 1]]
        finite = np.isfinite(dij) & np.isfinite(dkl)
        usable = distinct & on_side & finite
        left_smaller = usable & (dij + margin < dkl)
        right_smaller = usable & (dkl + margin < dij)

        rows = np.concatenate(
            (
                np.hstack((first[left_smaller], second[left_smaller])),
                np.hstack((second[right_smaller], first[right_smaller])),
            )
        )
        if rows.size:
            accepted.append(rows)
            accepted_count += rows.shape[0]
        attempts += batch

    if accepted_count < count:
        raise ValueError(
            f"could only sample {accepted_count} of {count} constraints; "
            "margin or split leaves too few comparable quadruples"
        )
    return np.vstack(accepted)[:count].astype(int, copy=False)


def build_constraint_split(
    dissimilarity: NDArray[np.float64],
    train_count: int,
    heldout_count: int,
    margin: float,
    train_fraction: float = 0.8,
    seed: int = 0,
) -> ConstraintSplit:
    """Sample leak-free train/held-out constraints per PC-V1 Section 6."""

    n = dissimilarity.shape[0]
    train_pairs = pair_split_matrix(n, train_fraction=train_fraction, seed=seed)
    heldout_pairs = ~train_pairs
    np.fill_diagonal(heldout_pairs, False)
    rng = np.random.default_rng(seed + 17)
    train = _sample_side_constraints(
        dissimilarity, train_pairs, train_count, margin, rng
    )
    heldout = _sample_side_constraints(
        dissimilarity, heldout_pairs, heldout_count, margin, rng
    )
    return ConstraintSplit(train=train, heldout=heldout, margin=float(margin))


def flip_constraint_sides(
    constraints: NDArray[np.int_],
    flip_probability: float = 0.5,
    seed: int = 0,
) -> NDArray[np.int_]:
    """Report-only C-FLIP control: swap constraint sides at random."""

    rng = np.random.default_rng(seed)
    flipped = np.asarray(constraints, dtype=int).copy()
    mask = rng.uniform(size=flipped.shape[0]) < flip_probability
    flipped[mask] = flipped[mask][:, [2, 3, 0, 1]]
    return flipped
