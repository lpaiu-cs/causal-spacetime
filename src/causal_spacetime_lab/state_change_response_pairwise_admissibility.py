"""Admissibility diagnostics for pairwise response-profile protocols."""

from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
    pairwise_response_dissimilarity,
    pairwise_response_dissimilarity_matrix,
    pairwise_response_order_inversion_rate,
)
from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile


def _tie_and_strict_fractions(values: np.ndarray) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if finite.size < 2:
        return 0.0, 0.0
    total = 0
    ties = 0
    strict = 0
    for i in range(finite.size - 1):
        for j in range(i + 1, finite.size):
            total += 1
            if finite[i] == finite[j]:
                ties += 1
            else:
                strict += 1
    return float(ties / total), float(strict / total)


def pairwise_protocol_admissibility_report(
    profile: EchoResponseProfile,
    protocol: PairwiseResponseComparisonProtocol,
) -> dict[str, float | str]:
    """Report finite admissibility diagnostics for one protocol."""

    matrix = pairwise_response_dissimilarity_matrix(profile, protocol)
    dissimilarity = pairwise_response_dissimilarity(profile, protocol)
    target_count = profile.target_event_ids.size
    pair_count = dissimilarity.valid_pair_mask.size
    finite_values = dissimilarity.dissimilarity_values[
        dissimilarity.valid_pair_mask
    ]
    tie_fraction, strict_fraction = _tie_and_strict_fractions(finite_values)
    return {
        "protocol_name": protocol.name,
        "method": protocol.method,
        "missing_policy": protocol.missing_policy,
        "target_count": float(target_count),
        "valid_pair_count": float(np.sum(dissimilarity.valid_pair_mask)),
        "valid_pair_fraction": float(np.mean(dissimilarity.valid_pair_mask))
        if pair_count
        else 0.0,
        "finite_value_fraction": float(np.mean(np.isfinite(matrix))),
        "diagonal_zero": float(np.allclose(np.diag(matrix), 0.0, equal_nan=False)),
        "symmetric": float(np.allclose(matrix, matrix.T, equal_nan=True)),
        "mean_dissimilarity": float(np.mean(finite_values))
        if finite_values.size
        else float("nan"),
        "tie_fraction": tie_fraction,
        "strict_pair_order_fraction": strict_fraction,
    }


def _tie_change_fraction(values_a: np.ndarray, values_b: np.ndarray) -> float:
    common = np.isfinite(values_a) & np.isfinite(values_b)
    indices = np.flatnonzero(common)
    if indices.size < 2:
        return float("nan")
    changes: list[bool] = []
    for outer in range(indices.size - 1):
        for inner in range(outer + 1, indices.size):
            a = indices[outer]
            b = indices[inner]
            tie_a = values_a[a] == values_a[b]
            tie_b = values_b[a] == values_b[b]
            changes.append(tie_a != tie_b)
    return float(np.mean(changes)) if changes else float("nan")


def compare_pairwise_protocols(
    profile: EchoResponseProfile,
    protocols: list[PairwiseResponseComparisonProtocol],
) -> list[dict[str, float | str]]:
    """Compare response-comparison order disagreement among protocols."""

    dissimilarities = [
        pairwise_response_dissimilarity(profile, protocol) for protocol in protocols
    ]
    rows: list[dict[str, float | str]] = []
    for i in range(len(dissimilarities) - 1):
        for j in range(i + 1, len(dissimilarities)):
            left = dissimilarities[i]
            right = dissimilarities[j]
            common_valid = left.valid_pair_mask & right.valid_pair_mask
            rows.append(
                {
                    "protocol_a": left.protocol_name,
                    "protocol_b": right.protocol_name,
                    "common_valid_pair_count": float(np.sum(common_valid)),
                    "order_inversion_rate": pairwise_response_order_inversion_rate(
                        left,
                        right,
                    ),
                    "tie_change_fraction": _tie_change_fraction(
                        left.dissimilarity_values,
                        right.dissimilarity_values,
                    ),
                }
            )
    return rows
