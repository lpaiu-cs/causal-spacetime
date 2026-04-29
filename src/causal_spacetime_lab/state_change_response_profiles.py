"""Multi-protocol response-profile diagnostics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.state_change_response_signature import (
    EchoResponseSignature,
)


@dataclass(frozen=True)
class EchoResponseProfile:
    """Delay-rank matrix for shared targets across protocol columns.

    A response profile is richer than one scalar response order, but it is
    still pre-metric. Protocol columns may encode reference choices, emission
    positions, gate values, or variant labels.
    """

    target_event_ids: NDArray[np.int_]
    protocol_labels: list[str]
    delay_rank_matrix: NDArray[np.int_]
    reachable_matrix: NDArray[np.bool_]


def _index_by_target(signature: EchoResponseSignature) -> dict[int, int]:
    return {
        int(target_id): int(index)
        for index, target_id in enumerate(signature.target_event_ids)
    }


def _common_target_ids(signatures: list[EchoResponseSignature]) -> NDArray[np.int_]:
    if not signatures:
        return np.empty(0, dtype=int)
    common = set(int(value) for value in signatures[0].target_event_ids)
    for signature in signatures[1:]:
        common &= {int(value) for value in signature.target_event_ids}
    return np.asarray(sorted(common), dtype=int)


def response_profile_from_signatures(
    signatures: list[EchoResponseSignature],
) -> EchoResponseProfile:
    """Align response signatures into a common-target profile matrix."""

    targets = _common_target_ids(signatures)
    delays = np.full((targets.size, len(signatures)), -1, dtype=int)
    reachable = np.zeros((targets.size, len(signatures)), dtype=bool)
    for column, signature in enumerate(signatures):
        lookup = _index_by_target(signature)
        for row, target in enumerate(targets):
            index = lookup[int(target)]
            delays[row, column] = int(signature.delay_ranks[index])
            reachable[row, column] = bool(signature.reachable_mask[index])
    return EchoResponseProfile(
        target_event_ids=targets,
        protocol_labels=[signature.label for signature in signatures],
        delay_rank_matrix=delays,
        reachable_matrix=reachable,
    )


def response_profile_reachable_fraction(profile: EchoResponseProfile) -> float:
    """Return the fraction of profile entries with reachable delay ranks."""

    if profile.reachable_matrix.size == 0:
        return 0.0
    return float(np.mean(profile.reachable_matrix))


def _row_profile_key(delays: NDArray[np.int_], reachable: NDArray[np.bool_]) -> tuple:
    return tuple(
        int(delay) if bool(is_reachable) else None
        for delay, is_reachable in zip(delays, reachable, strict=True)
    )


def response_profile_separation_fraction(profile: EchoResponseProfile) -> float:
    """Return fraction of target pairs separated by at least one protocol."""

    n_targets = profile.target_event_ids.size
    if n_targets < 2:
        return 0.0
    separated = 0
    total = 0
    for i in range(n_targets - 1):
        for j in range(i + 1, n_targets):
            total += 1
            comparable_columns = (
                profile.reachable_matrix[i] & profile.reachable_matrix[j]
            )
            if np.any(
                comparable_columns
                & (profile.delay_rank_matrix[i] != profile.delay_rank_matrix[j])
            ):
                separated += 1
    return float(separated / total) if total else 0.0


def response_profile_equivalence_classes(
    profile: EchoResponseProfile,
) -> list[NDArray[np.int_]]:
    """Group targets with identical reachable delay profiles."""

    groups: dict[tuple, list[int]] = {}
    for row, target_id in enumerate(profile.target_event_ids):
        key = _row_profile_key(
            profile.delay_rank_matrix[row],
            profile.reachable_matrix[row],
        )
        groups.setdefault(key, []).append(int(target_id))
    return [np.asarray(values, dtype=int) for values in groups.values()]


def compare_response_profiles(
    profile_a: EchoResponseProfile,
    profile_b: EchoResponseProfile,
) -> dict[str, float]:
    """Compare two response profiles on common targets and protocol labels."""

    common_targets = sorted(
        set(int(value) for value in profile_a.target_event_ids)
        & set(int(value) for value in profile_b.target_event_ids)
    )
    common_protocols = [
        label
        for label in profile_a.protocol_labels
        if label in profile_b.protocol_labels
    ]
    if not common_targets or not common_protocols:
        return {
            "common_target_count": float(len(common_targets)),
            "common_protocol_count": float(len(common_protocols)),
            "profile_entry_agreement_fraction": float("nan"),
            "profile_pair_separation_agreement": float("nan"),
            "reachable_agreement_fraction": float("nan"),
        }

    target_lookup_a = {
        int(target): index for index, target in enumerate(profile_a.target_event_ids)
    }
    target_lookup_b = {
        int(target): index for index, target in enumerate(profile_b.target_event_ids)
    }
    protocol_lookup_a = {
        label: index for index, label in enumerate(profile_a.protocol_labels)
    }
    protocol_lookup_b = {
        label: index for index, label in enumerate(profile_b.protocol_labels)
    }
    rows_a = np.asarray([target_lookup_a[target] for target in common_targets])
    rows_b = np.asarray([target_lookup_b[target] for target in common_targets])
    cols_a = np.asarray([protocol_lookup_a[label] for label in common_protocols])
    cols_b = np.asarray([protocol_lookup_b[label] for label in common_protocols])
    delays_a = profile_a.delay_rank_matrix[np.ix_(rows_a, cols_a)]
    delays_b = profile_b.delay_rank_matrix[np.ix_(rows_b, cols_b)]
    reachable_a = profile_a.reachable_matrix[np.ix_(rows_a, cols_a)]
    reachable_b = profile_b.reachable_matrix[np.ix_(rows_b, cols_b)]
    reachable_agreement = reachable_a == reachable_b
    common_reachable = reachable_a & reachable_b
    if np.any(common_reachable):
        entry_agreement = float(
            np.mean(delays_a[common_reachable] == delays_b[common_reachable])
        )
    else:
        entry_agreement = float("nan")

    aligned_a = EchoResponseProfile(
        target_event_ids=np.asarray(common_targets, dtype=int),
        protocol_labels=list(common_protocols),
        delay_rank_matrix=delays_a,
        reachable_matrix=reachable_a,
    )
    aligned_b = EchoResponseProfile(
        target_event_ids=np.asarray(common_targets, dtype=int),
        protocol_labels=list(common_protocols),
        delay_rank_matrix=delays_b,
        reachable_matrix=reachable_b,
    )
    separation_a = _pair_separation_flags(aligned_a)
    separation_b = _pair_separation_flags(aligned_b)
    if separation_a.size:
        separation_agreement = float(np.mean(separation_a == separation_b))
    else:
        separation_agreement = 0.0
    return {
        "common_target_count": float(len(common_targets)),
        "common_protocol_count": float(len(common_protocols)),
        "profile_entry_agreement_fraction": entry_agreement,
        "profile_pair_separation_agreement": separation_agreement,
        "reachable_agreement_fraction": float(np.mean(reachable_agreement)),
    }


def _pair_separation_flags(profile: EchoResponseProfile) -> NDArray[np.bool_]:
    flags: list[bool] = []
    n_targets = profile.target_event_ids.size
    for i in range(n_targets - 1):
        for j in range(i + 1, n_targets):
            comparable_columns = (
                profile.reachable_matrix[i] & profile.reachable_matrix[j]
            )
            flags.append(
                bool(
                    np.any(
                        comparable_columns
                        & (
                            profile.delay_rank_matrix[i]
                            != profile.delay_rank_matrix[j]
                        )
                    )
                )
            )
    return np.asarray(flags, dtype=bool)
