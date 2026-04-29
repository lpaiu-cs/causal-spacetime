"""Ordinal echo-response signatures for motif target populations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.state_change_echo_motif_validation import (
    recovered_delay_for_motif,
)
from causal_spacetime_lab.state_change_echo_motifs import EchoMotifRecord


@dataclass(frozen=True)
class EchoResponseSignature:
    """Order signature induced by recovered echo-delay ranks.

    This is an ordinal diagnostic over motif targets, not a metric space.
    Unreachable targets have delay rank -1 and do not contribute strict pair
    order.
    """

    target_event_ids: NDArray[np.int_]
    delay_ranks: NDArray[np.int_]
    reachable_mask: NDArray[np.bool_]
    order_sign_matrix: NDArray[np.int_]
    label: str


def _as_vector(values: ArrayLike, name: str) -> NDArray[np.int_]:
    vector = np.asarray(values, dtype=int)
    if vector.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    return vector


def response_order_sign_matrix(
    delay_ranks: ArrayLike,
    reachable_mask: ArrayLike,
    *,
    tie_value: int = 0,
) -> NDArray[np.int_]:
    """Return pairwise response-order signs from delay ranks.

    Entries are -1 when target i has smaller recovered echo-delay rank than
    target j, 1 when it has larger rank, and tie_value for equal, unresolved,
    or diagonal pairs.
    """

    delays = _as_vector(delay_ranks, "delay_ranks")
    reachable = np.asarray(reachable_mask, dtype=bool)
    if reachable.ndim != 1 or reachable.shape != delays.shape:
        raise ValueError("reachable_mask must match delay_ranks")

    n_targets = delays.size
    signs = np.full((n_targets, n_targets), int(tie_value), dtype=int)
    reachable_pairs = reachable[:, None] & reachable[None, :]
    smaller = delays[:, None] < delays[None, :]
    larger = delays[:, None] > delays[None, :]
    signs[reachable_pairs & smaller] = -1
    signs[reachable_pairs & larger] = 1
    np.fill_diagonal(signs, int(tie_value))
    return signs


def _signature_from_arrays(
    target_event_ids: ArrayLike,
    delay_ranks: ArrayLike,
    reachable_mask: ArrayLike,
    *,
    label: str,
) -> EchoResponseSignature:
    targets = _as_vector(target_event_ids, "target_event_ids")
    delays = _as_vector(delay_ranks, "delay_ranks")
    reachable = np.asarray(reachable_mask, dtype=bool)
    if targets.shape != delays.shape or targets.shape != reachable.shape:
        raise ValueError("targets, delay_ranks, and reachable_mask must match")
    return EchoResponseSignature(
        target_event_ids=targets,
        delay_ranks=delays,
        reachable_mask=reachable,
        order_sign_matrix=response_order_sign_matrix(delays, reachable),
        label=label,
    )


def echo_response_signature_from_motif_rows(
    motif_rows: list[dict[str, float | str]],
    *,
    label: str = "signature",
) -> EchoResponseSignature:
    """Build an echo-response signature from motif diagnostic rows."""

    targets: list[int] = []
    delays: list[int] = []
    reachable: list[bool] = []
    for row in motif_rows:
        targets.append(int(float(row["target_event_id"])))
        recovered = float(row["recovered_delay_rank"])
        is_reachable = bool(np.isfinite(recovered))
        delays.append(int(recovered) if is_reachable else -1)
        reachable.append(is_reachable)
    return _signature_from_arrays(
        targets,
        delays,
        reachable,
        label=label,
    )


def echo_response_signature_from_motifs(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    motifs: list[EchoMotifRecord],
    *,
    label: str = "signature",
) -> EchoResponseSignature:
    """Build a recovered response-order signature for motif targets."""

    targets: list[int] = []
    delays: list[int] = []
    reachable: list[bool] = []
    for motif in motifs:
        recovered = recovered_delay_for_motif(
            order_matrix,
            reference_chain_event_ids,
            motif,
        )
        targets.append(int(motif.target_event_id))
        delays.append(int(recovered) if recovered is not None else -1)
        reachable.append(recovered is not None)
    return _signature_from_arrays(
        targets,
        delays,
        reachable,
        label=label,
    )


def _reachable_pair_values(signature: EchoResponseSignature) -> NDArray[np.int_]:
    reachable_indices = np.flatnonzero(signature.reachable_mask)
    if reachable_indices.size < 2:
        return np.empty(0, dtype=int)
    signs: list[int] = []
    for outer_index, target_i in enumerate(reachable_indices[:-1]):
        for target_j in reachable_indices[outer_index + 1 :]:
            signs.append(int(signature.order_sign_matrix[target_i, target_j]))
    return np.asarray(signs, dtype=int)


def signature_tie_fraction(signature: EchoResponseSignature) -> float:
    """Return the tied-pair fraction among reachable target pairs."""

    signs = _reachable_pair_values(signature)
    if signs.size == 0:
        return 0.0
    return float(np.mean(signs == 0))


def signature_strict_pair_fraction(signature: EchoResponseSignature) -> float:
    """Return the strict nonzero pair fraction among reachable target pairs."""

    signs = _reachable_pair_values(signature)
    if signs.size == 0:
        return 0.0
    return float(np.mean(signs != 0))


def signature_reachable_fraction(signature: EchoResponseSignature) -> float:
    """Return the fraction of targets with recovered echo-delay ranks."""

    if signature.target_event_ids.size == 0:
        return 0.0
    return float(np.mean(signature.reachable_mask))

