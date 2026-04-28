"""Diagnostics for order-level reference-chain brackets."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from causal_spacetime_lab.ordinal import order_agreement_rate, order_inversion_rate


def bracket_coverage_summary(
    predecessor_positions: ArrayLike,
    successor_positions: ArrayLike,
    accessible_mask: ArrayLike,
    reference_chain_length: int,
) -> dict[str, float]:
    """Summarize two-sided bracket coverage for one reference chain."""

    predecessors = np.asarray(predecessor_positions, dtype=int)
    successors = np.asarray(successor_positions, dtype=int)
    accessible = np.asarray(accessible_mask, dtype=bool)
    if predecessors.shape != successors.shape or predecessors.shape != accessible.shape:
        raise ValueError("predecessor, successor, and mask arrays must match")
    if reference_chain_length < 0:
        raise ValueError("reference_chain_length must be nonnegative")
    has_predecessor = predecessors >= 0
    has_successor = successors >= 0
    widths = successors[accessible] - predecessors[accessible]
    return {
        "target_count": float(predecessors.size),
        "reference_chain_length": float(reference_chain_length),
        "accessible_count": float(np.count_nonzero(accessible)),
        "accessible_fraction": float(np.mean(accessible)) if accessible.size else 0.0,
        "predecessor_only_count": float(
            np.count_nonzero(has_predecessor & ~has_successor)
        ),
        "successor_only_count": float(
            np.count_nonzero(~has_predecessor & has_successor)
        ),
        "neither_count": float(np.count_nonzero(~has_predecessor & ~has_successor)),
        "mean_predecessor_position": float(np.mean(predecessors[has_predecessor]))
        if np.any(has_predecessor)
        else float("nan"),
        "mean_successor_position": float(np.mean(successors[has_successor]))
        if np.any(has_successor)
        else float("nan"),
        "mean_bracket_width_rank": float(np.mean(widths)) if widths.size else 0.0,
        "median_bracket_width_rank": float(np.median(widths)) if widths.size else 0.0,
        "max_bracket_width_rank": float(np.max(widths)) if widths.size else 0.0,
    }


def rank_slice_summary(slice_labels: ArrayLike) -> dict[str, float]:
    """Summarize assigned rank-slice labels."""

    labels = np.asarray(slice_labels, dtype=int)
    assigned = labels >= 0
    assigned_labels = labels[assigned]
    if assigned_labels.size:
        _, counts = np.unique(assigned_labels, return_counts=True)
    else:
        counts = np.empty(0, dtype=int)
    return {
        "slice_count": float(counts.size),
        "assigned_count": float(assigned_labels.size),
        "mean_slice_size": float(np.mean(counts)) if counts.size else 0.0,
        "median_slice_size": float(np.median(counts)) if counts.size else 0.0,
        "max_slice_size": float(np.max(counts)) if counts.size else 0.0,
        "singleton_slice_fraction": float(np.count_nonzero(counts == 1) / counts.size)
        if counts.size
        else 0.0,
    }


def compare_bracket_rank_orders(
    ranks_a: ArrayLike,
    ranks_b: ArrayLike,
    accessible_a: ArrayLike,
    accessible_b: ArrayLike,
) -> dict[str, float]:
    """Compare rank orderings on events accessible to two reference chains."""

    rank_values_a = np.asarray(ranks_a, dtype=float)
    rank_values_b = np.asarray(ranks_b, dtype=float)
    mask_a = np.asarray(accessible_a, dtype=bool)
    mask_b = np.asarray(accessible_b, dtype=bool)
    if (
        rank_values_a.shape != rank_values_b.shape
        or rank_values_a.shape != mask_a.shape
        or rank_values_a.shape != mask_b.shape
    ):
        raise ValueError("rank and accessibility arrays must match")
    common = mask_a & mask_b
    common_count = int(np.count_nonzero(common))
    if common_count < 2:
        inversion = float("nan")
        agreement = float("nan")
    else:
        inversion = order_inversion_rate(
            rank_values_a[common],
            rank_values_b[common],
            ignore_ties=True,
        )
        agreement = order_agreement_rate(
            rank_values_a[common],
            rank_values_b[common],
            ignore_ties=True,
        )
    return {
        "common_accessible_count": float(common_count),
        "common_accessible_fraction": float(common_count / rank_values_a.size)
        if rank_values_a.size
        else 0.0,
        "order_inversion_rate": inversion,
        "order_agreement_rate": agreement,
    }
