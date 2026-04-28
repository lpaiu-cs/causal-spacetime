"""Ambiguity diagnostics for observer-like chain candidates."""

from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_observers import ObserverChainCandidate


def top_score_gap(ranked_rows: list[dict[str, float | str]]) -> float:
    """Return best score minus second-best score, or NaN if unavailable."""

    if len(ranked_rows) < 2:
        return float("nan")
    return float(ranked_rows[0]["score"]) - float(ranked_rows[1]["score"])


def candidate_overlap_fraction(
    candidate_a: ObserverChainCandidate,
    candidate_b: ObserverChainCandidate,
) -> float:
    """Return Jaccard overlap of two candidate chain event-id sets."""

    set_a = set(np.asarray(candidate_a.chain_event_ids, dtype=int).tolist())
    set_b = set(np.asarray(candidate_b.chain_event_ids, dtype=int).tolist())
    union = set_a | set_b
    if not union:
        return 1.0
    return float(len(set_a & set_b) / len(union))


def chain_candidate_diversity(
    candidates: list[ObserverChainCandidate],
) -> dict[str, float]:
    """Return pairwise overlap diagnostics for candidate chains."""

    overlaps: list[float] = []
    for i, candidate_a in enumerate(candidates):
        for candidate_b in candidates[i + 1 :]:
            overlaps.append(candidate_overlap_fraction(candidate_a, candidate_b))
    return {
        "candidate_count": float(len(candidates)),
        "mean_pairwise_overlap": float(np.mean(overlaps)) if overlaps else float("nan"),
        "min_pairwise_overlap": float(np.min(overlaps)) if overlaps else float("nan"),
        "max_pairwise_overlap": float(np.max(overlaps)) if overlaps else float("nan"),
    }
