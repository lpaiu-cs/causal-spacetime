"""Comparison utilities for echo-response order signatures."""

from __future__ import annotations

from collections import Counter

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.state_change_response_signature import (
    EchoResponseSignature,
)


def _common_target_ids(
    signatures: list[EchoResponseSignature],
) -> NDArray[np.int_]:
    if not signatures:
        return np.empty(0, dtype=int)
    common = set(int(value) for value in signatures[0].target_event_ids)
    for signature in signatures[1:]:
        common &= {int(value) for value in signature.target_event_ids}
    return np.asarray(sorted(common), dtype=int)


def _index_by_target(signature: EchoResponseSignature) -> dict[int, int]:
    return {
        int(target_id): int(index)
        for index, target_id in enumerate(signature.target_event_ids)
    }


def _aligned_sign_matrix(
    signature: EchoResponseSignature,
    target_event_ids: NDArray[np.int_],
) -> tuple[NDArray[np.int_], NDArray[np.bool_]]:
    lookup = _index_by_target(signature)
    indices = np.asarray(
        [lookup[int(target)] for target in target_event_ids],
        dtype=int,
    )
    return (
        signature.order_sign_matrix[np.ix_(indices, indices)],
        signature.reachable_mask[indices],
    )


def _unordered_pair_indices(n_targets: int) -> list[tuple[int, int]]:
    return [
        (outer, inner)
        for outer in range(max(0, n_targets - 1))
        for inner in range(outer + 1, n_targets)
    ]


def _tie_fraction_for_pairs(
    signs: NDArray[np.int_],
    pairs: list[tuple[int, int]],
) -> float:
    if not pairs:
        return 0.0
    return float(np.mean([signs[i, j] == 0 for i, j in pairs]))


def compare_response_signatures(
    baseline: EchoResponseSignature,
    variant: EchoResponseSignature,
) -> dict[str, float]:
    """Compare response-order signatures on common reachable target pairs."""

    common_targets = _common_target_ids([baseline, variant])
    if common_targets.size == 0:
        return {
            "common_target_count": 0.0,
            "common_reachable_count": 0.0,
            "common_reachable_fraction": 0.0,
            "pair_count": 0.0,
            "strict_pair_count_baseline": 0.0,
            "strict_pair_count_variant": 0.0,
            "pair_agreement_fraction": float("nan"),
            "pair_disagreement_fraction": float("nan"),
            "pair_tie_changed_fraction": float("nan"),
            "baseline_tie_fraction": float("nan"),
            "variant_tie_fraction": float("nan"),
        }

    baseline_signs, baseline_reachable = _aligned_sign_matrix(
        baseline,
        common_targets,
    )
    variant_signs, variant_reachable = _aligned_sign_matrix(variant, common_targets)
    common_reachable = baseline_reachable & variant_reachable
    reachable_indices = np.flatnonzero(common_reachable)
    pair_positions = _unordered_pair_indices(reachable_indices.size)
    pairs = [
        (int(reachable_indices[i]), int(reachable_indices[j]))
        for i, j in pair_positions
    ]
    pair_count = len(pairs)
    strict_baseline = [baseline_signs[i, j] != 0 for i, j in pairs]
    strict_variant = [variant_signs[i, j] != 0 for i, j in pairs]
    strict_both = [
        index
        for index, (is_strict_baseline, is_strict_variant) in enumerate(
            zip(strict_baseline, strict_variant, strict=True)
        )
        if is_strict_baseline and is_strict_variant
    ]
    if strict_both:
        agreements = [
            baseline_signs[pairs[index]] == variant_signs[pairs[index]]
            for index in strict_both
        ]
        pair_agreement = float(np.mean(agreements))
        pair_disagreement = float(1.0 - pair_agreement)
    else:
        pair_agreement = float("nan")
        pair_disagreement = float("nan")
    tie_changed = [
        (baseline_signs[i, j] == 0) != (variant_signs[i, j] == 0)
        for i, j in pairs
    ]
    return {
        "common_target_count": float(common_targets.size),
        "common_reachable_count": float(reachable_indices.size),
        "common_reachable_fraction": float(
            reachable_indices.size / common_targets.size
        ),
        "pair_count": float(pair_count),
        "strict_pair_count_baseline": float(np.sum(strict_baseline)),
        "strict_pair_count_variant": float(np.sum(strict_variant)),
        "pair_agreement_fraction": pair_agreement,
        "pair_disagreement_fraction": pair_disagreement,
        "pair_tie_changed_fraction": float(np.mean(tie_changed))
        if tie_changed
        else 0.0,
        "baseline_tie_fraction": _tie_fraction_for_pairs(baseline_signs, pairs),
        "variant_tie_fraction": _tie_fraction_for_pairs(variant_signs, pairs),
    }


def stable_response_order_core(
    signatures: list[EchoResponseSignature],
    *,
    min_agreement_fraction: float = 0.8,
    require_nonzero: bool = True,
) -> dict[str, object]:
    """Return pairwise signs stable across enough protocol variants."""

    if not 0.0 <= min_agreement_fraction <= 1.0:
        raise ValueError("min_agreement_fraction must be in [0, 1]")
    targets = _common_target_ids(signatures)
    n_targets = targets.size
    stable_signs = np.zeros((n_targets, n_targets), dtype=int)
    stable_mask = np.zeros((n_targets, n_targets), dtype=bool)
    if not signatures or n_targets < 2:
        return {
            "target_event_ids": targets,
            "stable_order_sign_matrix": stable_signs,
            "stable_pair_mask": stable_mask,
            "stable_pair_fraction": 0.0,
            "variant_count": float(len(signatures)),
        }

    aligned = [_aligned_sign_matrix(signature, targets)[0] for signature in signatures]
    stable_pair_count = 0
    total_pair_count = 0
    for i, j in _unordered_pair_indices(n_targets):
        total_pair_count += 1
        values = [int(signs[i, j]) for signs in aligned]
        if require_nonzero:
            values = [value for value in values if value != 0]
        if not values:
            continue
        counts = Counter(values)
        sign, count = counts.most_common(1)[0]
        if sign == 0 and require_nonzero:
            continue
        if count / len(values) >= min_agreement_fraction:
            stable_signs[i, j] = int(sign)
            stable_signs[j, i] = -int(sign)
            stable_mask[i, j] = True
            stable_mask[j, i] = True
            stable_pair_count += 1
    return {
        "target_event_ids": targets,
        "stable_order_sign_matrix": stable_signs,
        "stable_pair_mask": stable_mask,
        "stable_pair_fraction": float(stable_pair_count / total_pair_count)
        if total_pair_count
        else 0.0,
        "variant_count": float(len(signatures)),
    }


def consensus_response_order_matrix(
    signatures: list[EchoResponseSignature],
    *,
    min_agreement_fraction: float = 0.5,
) -> NDArray[np.int_]:
    """Return a majority-vote response-order sign matrix over variants."""

    core = stable_response_order_core(
        signatures,
        min_agreement_fraction=min_agreement_fraction,
        require_nonzero=False,
    )
    return np.asarray(core["stable_order_sign_matrix"], dtype=int)


def response_order_cycle_count(order_sign_matrix: NDArray[np.int_]) -> int:
    """Count directed 3-cycles in a response-order sign matrix.

    Sign -1 is interpreted as target i preceding target j in response-rank
    order. This finite diagnostic is not a representability theorem.
    """

    signs = np.asarray(order_sign_matrix, dtype=int)
    if signs.ndim != 2 or signs.shape[0] != signs.shape[1]:
        raise ValueError("order_sign_matrix must be square")
    precedes = signs == -1
    count = 0
    n_targets = signs.shape[0]
    for i in range(n_targets - 2):
        for j in range(i + 1, n_targets - 1):
            for k in range(j + 1, n_targets):
                cycle_forward = precedes[i, j] and precedes[j, k] and precedes[k, i]
                cycle_reverse = precedes[i, k] and precedes[k, j] and precedes[j, i]
                if cycle_forward or cycle_reverse:
                    count += 1
    return int(count)
