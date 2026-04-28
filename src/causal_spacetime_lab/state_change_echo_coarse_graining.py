"""Coarse-graining utilities for return-spectrum echo diagnostics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.state_change_echo_interference import (
    earliest_echo_delay_from_spectrum,
)
from causal_spacetime_lab.state_change_echo_motifs import EchoMotifRecord


@dataclass(frozen=True)
class EventThinningResult:
    """Result of closure-preserving event thinning."""

    retained_indices: NDArray[np.int_]
    old_to_new: dict[int, int]
    new_to_old: NDArray[np.int_]
    restricted_order_matrix: NDArray[np.bool_]


@dataclass(frozen=True)
class ReferenceSubsamplingResult:
    """Result of reference-chain position subsampling."""

    retained_reference_positions: NDArray[np.int_]
    subsampled_reference_chain: NDArray[np.int_]
    old_position_to_new: dict[int, int]


def _as_index_vector(values: ArrayLike, name: str) -> NDArray[np.int_]:
    indices = np.asarray(values, dtype=int)
    if indices.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    return indices


def protected_indices_for_reference_and_motifs(
    reference_chain_event_ids: ArrayLike,
    motifs: list[EchoMotifRecord],
    *,
    include_motif_paths: bool = False,
) -> NDArray[np.int_]:
    """Return event ids protected from event thinning."""

    protected: list[int] = [
        int(value) for value in _as_index_vector(reference_chain_event_ids, "reference")
    ]
    for motif in motifs:
        protected.append(int(motif.target_event_id))
        if include_motif_paths:
            protected.extend(int(value) for value in motif.outward_event_ids)
            protected.extend(int(value) for value in motif.return_event_ids)
            protected.append(int(motif.reference_emission_event_id))
            protected.append(int(motif.reference_return_event_id))
    return np.asarray(sorted(set(protected)), dtype=int)


def sample_retained_indices(
    event_count: int,
    keep_probability: float,
    protected_indices: ArrayLike,
    seed: int | None = None,
) -> NDArray[np.int_]:
    """Sample retained event ids while always keeping protected ids."""

    if event_count < 0:
        raise ValueError("event_count must be nonnegative")
    if not 0.0 <= keep_probability <= 1.0:
        raise ValueError("keep_probability must be in [0, 1]")
    protected = _as_index_vector(protected_indices, "protected_indices")
    if protected.size and (np.min(protected) < 0 or np.max(protected) >= event_count):
        raise IndexError("protected index outside event range")
    rng = np.random.default_rng(seed)
    keep = rng.random(event_count) < keep_probability
    keep[protected] = True
    return np.flatnonzero(keep).astype(int)


def restrict_transitive_order_to_retained_events(
    order_matrix: ArrayLike,
    retained_indices: ArrayLike,
) -> EventThinningResult:
    """Restrict a transitive closure to retained events.

    This is closure-preserving event thinning: causal reachability among
    retained events is inherited from the original transitive closure. Hidden
    intermediate events are integrated out rather than interpreted as removed
    causal paths.
    """

    matrix = np.asarray(order_matrix, dtype=bool)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("order_matrix must be square")
    retained = np.asarray(sorted(set(_as_index_vector(retained_indices, "retained"))))
    if retained.size and (np.min(retained) < 0 or np.max(retained) >= matrix.shape[0]):
        raise IndexError("retained index outside order matrix")
    old_to_new = {int(old): int(new) for new, old in enumerate(retained)}
    restricted = matrix[np.ix_(retained, retained)]
    return EventThinningResult(
        retained_indices=retained.astype(int),
        old_to_new=old_to_new,
        new_to_old=retained.astype(int),
        restricted_order_matrix=restricted,
    )


def remap_reference_chain(
    reference_chain_event_ids: ArrayLike,
    old_to_new: dict[int, int],
) -> NDArray[np.int_]:
    """Remap a reference chain through an event-thinning index map."""

    chain = _as_index_vector(reference_chain_event_ids, "reference_chain_event_ids")
    missing = [int(value) for value in chain if int(value) not in old_to_new]
    if missing:
        raise KeyError(f"reference event ids missing after thinning: {missing}")
    return np.asarray([old_to_new[int(value)] for value in chain], dtype=int)


def remap_echo_motif_record_for_event_thinning(
    motif: EchoMotifRecord,
    old_to_new: dict[int, int],
) -> EchoMotifRecord | None:
    """Remap a motif after closure-preserving event thinning.

    Returns None if the motif target or required reference emission/return
    event is missing. Intermediate path ids are retained only when present.
    """

    required = (
        motif.target_event_id,
        motif.reference_emission_event_id,
        motif.reference_return_event_id,
    )
    if any(int(value) not in old_to_new for value in required):
        return None

    def remap_optional(values: NDArray[np.int_]) -> NDArray[np.int_]:
        return np.asarray(
            [old_to_new[int(value)] for value in values if int(value) in old_to_new],
            dtype=int,
        )

    return EchoMotifRecord(
        target_event_id=old_to_new[int(motif.target_event_id)],
        emission_position=int(motif.emission_position),
        planted_return_position=int(motif.planted_return_position),
        planted_delay_rank=int(motif.planted_delay_rank),
        outward_event_ids=remap_optional(motif.outward_event_ids),
        return_event_ids=remap_optional(motif.return_event_ids),
        reference_emission_event_id=old_to_new[int(motif.reference_emission_event_id)],
        reference_return_event_id=old_to_new[int(motif.reference_return_event_id)],
        motif_kind=motif.motif_kind,
    )


def spectrum_jaccard(spectrum_a: ArrayLike, spectrum_b: ArrayLike) -> float:
    """Return Jaccard agreement between two finite delay-rank sets."""

    set_a = {int(value) for value in np.asarray(spectrum_a, dtype=int)}
    set_b = {int(value) for value in np.asarray(spectrum_b, dtype=int)}
    if not set_a and not set_b:
        return 1.0
    return float(len(set_a & set_b) / len(set_a | set_b))


def return_spectrum_stability_report(
    baseline_spectrum: ArrayLike,
    coarse_spectrum: ArrayLike,
) -> dict[str, float]:
    """Compare baseline and coarse return spectra."""

    baseline = np.asarray(baseline_spectrum, dtype=int)
    coarse = np.asarray(coarse_spectrum, dtype=int)
    baseline_earliest = earliest_echo_delay_from_spectrum(baseline)
    coarse_earliest = earliest_echo_delay_from_spectrum(coarse)
    shift = (
        float(coarse_earliest - baseline_earliest)
        if baseline_earliest is not None and coarse_earliest is not None
        else float("nan")
    )
    return {
        "baseline_size": float(baseline.size),
        "coarse_size": float(coarse.size),
        "spectrum_jaccard": spectrum_jaccard(baseline, coarse),
        "baseline_earliest": float(baseline_earliest)
        if baseline_earliest is not None
        else float("nan"),
        "coarse_earliest": float(coarse_earliest)
        if coarse_earliest is not None
        else float("nan"),
        "earliest_delay_shift": shift,
        "exact_earliest_preserved": float(baseline_earliest == coarse_earliest),
    }


def subsample_reference_chain_positions(
    reference_chain_event_ids: ArrayLike,
    stride: int,
    *,
    protected_positions: ArrayLike | None = None,
    keep_first_last: bool = True,
) -> ReferenceSubsamplingResult:
    """Subsample reference-chain positions by stride.

    This creates coarse reference ranks. It is not clock calibration.
    """

    if stride < 1:
        raise ValueError("stride must be at least 1")
    chain = _as_index_vector(reference_chain_event_ids, "reference_chain_event_ids")
    positions: set[int] = set(range(0, chain.size, stride))
    if keep_first_last and chain.size:
        positions.add(0)
        positions.add(chain.size - 1)
    if protected_positions is not None:
        protected = _as_index_vector(protected_positions, "protected_positions")
        if protected.size and (
            np.min(protected) < 0 or np.max(protected) >= chain.size
        ):
            raise IndexError("protected reference position outside chain")
        positions.update(int(position) for position in protected)
    retained = np.asarray(sorted(positions), dtype=int)
    old_position_to_new = {
        int(old): int(new) for new, old in enumerate(retained)
    }
    return ReferenceSubsamplingResult(
        retained_reference_positions=retained,
        subsampled_reference_chain=chain[retained],
        old_position_to_new=old_position_to_new,
    )


def coarse_emission_position_for_motif(
    motif: EchoMotifRecord,
    subsampling: ReferenceSubsamplingResult,
) -> int | None:
    """Return the coarse emission position when retained."""

    return subsampling.old_position_to_new.get(int(motif.emission_position))


def expected_coarse_return_position_for_motif(
    motif: EchoMotifRecord,
    subsampling: ReferenceSubsamplingResult,
) -> int | None:
    """Return the first retained coarse return position at or after the plant."""

    retained = subsampling.retained_reference_positions
    mask = retained >= int(motif.planted_return_position)
    if not np.any(mask):
        return None
    original_position = int(retained[np.flatnonzero(mask)[0]])
    return subsampling.old_position_to_new[original_position]


def expected_coarse_delay_rank_for_motif(
    motif: EchoMotifRecord,
    subsampling: ReferenceSubsamplingResult,
) -> int | None:
    """Return the expected coarse delay rank under reference subsampling."""

    emission = coarse_emission_position_for_motif(motif, subsampling)
    returned = expected_coarse_return_position_for_motif(motif, subsampling)
    if emission is None or returned is None:
        return None
    return int(returned - emission)
