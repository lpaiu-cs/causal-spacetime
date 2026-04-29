"""Layered controlled echo-response motif helpers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifRecord,
    EchoMotifSpec,
)
from causal_spacetime_lab.state_change_response_signature import (
    EchoResponseSignature,
    response_order_sign_matrix,
)


@dataclass(frozen=True)
class EchoResponseLayerSpec:
    """Specification for one planted response-rank layer."""

    delay_rank: int
    target_count: int
    outward_steps: int = 1
    return_steps: int = 1
    layer_label: str = ""

    def __post_init__(self) -> None:
        if self.delay_rank < 1:
            raise ValueError("delay_rank must be at least 1")
        if self.target_count < 1:
            raise ValueError("target_count must be at least 1")
        if self.outward_steps < 0:
            raise ValueError("outward_steps must be nonnegative")
        if self.return_steps < 0:
            raise ValueError("return_steps must be nonnegative")


def build_layered_echo_motif_specs(
    reference_chain_event_ids: ArrayLike,
    emission_position: int,
    layer_specs: list[EchoResponseLayerSpec],
    seed: int | None = None,
) -> list[EchoMotifSpec]:
    """Build controlled motif specs grouped by planted delay-rank layers."""

    chain = np.asarray(reference_chain_event_ids, dtype=int)
    if chain.ndim != 1:
        raise ValueError("reference_chain_event_ids must be one-dimensional")
    if emission_position < 0 or emission_position >= chain.size:
        raise IndexError("emission_position is outside the reference chain")
    specs: list[EchoMotifSpec] = []
    for layer in layer_specs:
        if emission_position + layer.delay_rank >= chain.size:
            raise ValueError("layer delay extends beyond the reference chain")
        motif_kind = layer.layer_label or f"layer_delay_{layer.delay_rank}"
        specs.extend(
            EchoMotifSpec(
                emission_position=emission_position,
                planted_delay_rank=layer.delay_rank,
                outward_steps=layer.outward_steps,
                return_steps=layer.return_steps,
                motif_kind=motif_kind,
            )
            for _ in range(layer.target_count)
        )
    rng = np.random.default_rng(seed)
    rng.shuffle(specs)
    return specs


def planted_layer_labels_for_motifs(
    motifs: list[EchoMotifRecord],
) -> dict[int, int]:
    """Return target event id to planted delay-rank labels."""

    return {
        int(motif.target_event_id): int(motif.planted_delay_rank)
        for motif in motifs
    }


def planted_response_signature_from_motifs(
    motifs: list[EchoMotifRecord],
    *,
    label: str = "planted",
) -> EchoResponseSignature:
    """Build the planted ordinal response signature for controlled motifs."""

    target_event_ids = np.asarray(
        [int(motif.target_event_id) for motif in motifs],
        dtype=int,
    )
    delay_ranks = np.asarray(
        [int(motif.planted_delay_rank) for motif in motifs],
        dtype=int,
    )
    reachable_mask = np.ones(target_event_ids.size, dtype=bool)
    return EchoResponseSignature(
        target_event_ids=target_event_ids,
        delay_ranks=delay_ranks,
        reachable_mask=reachable_mask,
        order_sign_matrix=response_order_sign_matrix(delay_ranks, reachable_mask),
        label=label,
    )

