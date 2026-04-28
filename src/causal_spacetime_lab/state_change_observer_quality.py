"""Quality diagnostics for observer-like chain candidates."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.state_change import StateChangeNetwork
from causal_spacetime_lab.state_change_observers import (
    ObserverChainCandidate,
    chain_bracketing_mask,
    chain_comparability_mask,
    is_chain,
)
from causal_spacetime_lab.state_change_order import causal_interval_indices


@dataclass(frozen=True)
class ChainQualityReport:
    """Heuristic observer-like chain quality report."""

    name: str
    source: str
    chain_length: int
    event_count: int
    is_valid_chain: bool
    comparable_fraction: float
    bracketed_fraction: float
    mean_chain_gap_interval_size: float
    std_chain_gap_interval_size: float
    chain_gap_interval_cv: float
    max_chain_gap_interval_size: float
    local_system_purity: float


def _interval_sizes_for_adjacent_chain_gaps(
    order_matrix: NDArray[np.bool_],
    chain: NDArray[np.int_],
) -> NDArray[np.float64]:
    sizes: list[float] = []
    for left, right in zip(chain[:-1], chain[1:], strict=False):
        if order_matrix[int(left), int(right)]:
            sizes.append(
                float(causal_interval_indices(order_matrix, int(left), int(right)).size)
            )
    return np.asarray(sizes, dtype=float)


def _local_system_purity(
    network: StateChangeNetwork,
    chain: NDArray[np.int_],
) -> float:
    if chain.size == 0:
        return 0.0
    system_ids = np.asarray(
        [network.events[int(event_id)].system_id for event_id in chain]
    )
    _, counts = np.unique(system_ids, return_counts=True)
    return float(np.max(counts) / chain.size)


def evaluate_chain_candidate(
    network: StateChangeNetwork,
    order_matrix: ArrayLike,
    candidate: ObserverChainCandidate,
) -> ChainQualityReport:
    """Evaluate one candidate observer-like chain.

    Local-system purity uses event metadata, so it is not order-only.
    """

    matrix = np.asarray(order_matrix, dtype=bool)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("order_matrix must be square")
    event_count = matrix.shape[0]
    chain = np.asarray(candidate.chain_event_ids, dtype=int)
    valid = is_chain(matrix, chain)
    comparable = chain_comparability_mask(matrix, chain)
    bracketed = chain_bracketing_mask(matrix, chain)
    chain_set = set(chain.tolist())
    non_chain_count = event_count - len(chain_set)
    interval_sizes = _interval_sizes_for_adjacent_chain_gaps(matrix, chain)
    mean_interval = float(np.mean(interval_sizes)) if interval_sizes.size else 0.0
    std_interval = float(np.std(interval_sizes)) if interval_sizes.size else 0.0
    cv = float(std_interval / mean_interval) if mean_interval > 0.0 else float("nan")
    return ChainQualityReport(
        name=candidate.name,
        source=candidate.source,
        chain_length=int(chain.size),
        event_count=int(event_count),
        is_valid_chain=bool(valid),
        comparable_fraction=float(np.mean(comparable)) if event_count else 0.0,
        bracketed_fraction=float(np.count_nonzero(bracketed) / non_chain_count)
        if non_chain_count
        else 0.0,
        mean_chain_gap_interval_size=mean_interval,
        std_chain_gap_interval_size=std_interval,
        chain_gap_interval_cv=cv,
        max_chain_gap_interval_size=float(np.max(interval_sizes))
        if interval_sizes.size
        else 0.0,
        local_system_purity=_local_system_purity(network, chain),
    )


def rank_chain_candidates(
    reports: list[ChainQualityReport],
    *,
    comparable_weight: float = 0.35,
    bracket_weight: float = 0.45,
    length_weight: float = 0.10,
    regularity_weight: float = 0.10,
) -> list[dict[str, float | str]]:
    """Rank candidates by a heuristic chain quality score.

    The score is not a proof of physical observerhood.
    """

    rows: list[dict[str, float | str]] = []
    for report in reports:
        normalized_length = (
            report.chain_length / report.event_count if report.event_count else 0.0
        )
        cv = report.chain_gap_interval_cv
        regularity = 1.0 / (1.0 + cv) if np.isfinite(cv) else 0.0
        score = (
            comparable_weight * report.comparable_fraction
            + bracket_weight * report.bracketed_fraction
            + length_weight * normalized_length
            + regularity_weight * regularity
        )
        if not report.is_valid_chain:
            score = 0.0
        row = asdict(report)
        row.update(
            {
                "normalized_chain_length": float(normalized_length),
                "regularity_score": float(regularity),
                "score": float(score),
            }
        )
        rows.append(row)
    rows.sort(key=lambda row: float(row["score"]), reverse=True)
    for rank, row in enumerate(rows, start=1):
        row["rank"] = float(rank)
    return rows
