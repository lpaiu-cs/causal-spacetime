"""Summary diagnostics for same-emission echo-order ranks."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def _validate_echo_arrays(
    delay_ranks: ArrayLike,
    reachable_mask: ArrayLike,
) -> tuple[np.ndarray, np.ndarray]:
    delays = np.asarray(delay_ranks, dtype=int)
    reachable = np.asarray(reachable_mask, dtype=bool)
    if delays.ndim != 1 or reachable.ndim != 1 or delays.shape != reachable.shape:
        raise ValueError("delay_ranks and reachable_mask must be matching vectors")
    return delays, reachable


def echo_reachability_summary(
    return_positions: ArrayLike,
    delay_ranks: ArrayLike,
    reachable_mask: ArrayLike,
    emission_position: int,
    reference_chain_length: int,
) -> dict[str, float]:
    """Summarize echo reachability and return-rank distribution."""

    returns = np.asarray(return_positions, dtype=int)
    delays, reachable = _validate_echo_arrays(delay_ranks, reachable_mask)
    if returns.ndim != 1 or returns.shape != delays.shape:
        raise ValueError("return_positions must match delay_ranks")
    if reference_chain_length < 0:
        raise ValueError("reference_chain_length must be nonnegative")
    reachable_returns = returns[reachable]
    reachable_delays = delays[reachable]
    if reachable_delays.size:
        _, counts = np.unique(reachable_delays, return_counts=True)
        singleton_delay_fraction = float(np.count_nonzero(counts == 1) / counts.size)
        distinct_delay_rank_count = float(counts.size)
    else:
        singleton_delay_fraction = 0.0
        distinct_delay_rank_count = 0.0
    return {
        "target_count": float(delays.size),
        "reachable_count": float(np.count_nonzero(reachable)),
        "reachable_fraction": float(np.mean(reachable)) if reachable.size else 0.0,
        "emission_position": float(emission_position),
        "reference_chain_length": float(reference_chain_length),
        "mean_return_position": float(np.mean(reachable_returns))
        if reachable_returns.size
        else float("nan"),
        "median_return_position": float(np.median(reachable_returns))
        if reachable_returns.size
        else float("nan"),
        "mean_echo_delay_rank": float(np.mean(reachable_delays))
        if reachable_delays.size
        else float("nan"),
        "median_echo_delay_rank": float(np.median(reachable_delays))
        if reachable_delays.size
        else float("nan"),
        "max_echo_delay_rank": float(np.max(reachable_delays))
        if reachable_delays.size
        else 0.0,
        "singleton_delay_fraction": singleton_delay_fraction,
        "distinct_delay_rank_count": distinct_delay_rank_count,
    }


def echo_delay_histogram(
    delay_ranks: ArrayLike,
    reachable_mask: ArrayLike,
) -> dict[int, int]:
    """Return delay-rank counts for reachable targets."""

    delays, reachable = _validate_echo_arrays(delay_ranks, reachable_mask)
    histogram: dict[int, int] = {}
    if np.any(reachable):
        values, counts = np.unique(delays[reachable], return_counts=True)
        histogram = {
            int(value): int(count) for value, count in zip(values, counts, strict=True)
        }
    return histogram


def echo_order_resolution_summary(
    delay_ranks: ArrayLike,
    reachable_mask: ArrayLike,
) -> dict[str, float]:
    """Summarize how many reachable target pairs are tied or strictly ordered."""

    delays, reachable = _validate_echo_arrays(delay_ranks, reachable_mask)
    reachable_delays = delays[reachable]
    reachable_count = int(reachable_delays.size)
    pair_count = reachable_count * (reachable_count - 1) // 2
    if pair_count:
        _, counts = np.unique(reachable_delays, return_counts=True)
        tied_pairs = int(np.sum(counts * (counts - 1) // 2))
        tied_fraction = float(tied_pairs / pair_count)
    else:
        tied_fraction = float("nan")
    strict_fraction = (
        float("nan") if not np.isfinite(tied_fraction) else 1.0 - tied_fraction
    )
    return {
        "reachable_count": float(reachable_count),
        "comparable_pair_count": float(pair_count),
        "tied_pair_fraction": tied_fraction,
        "strict_order_pair_fraction": strict_fraction,
    }
