"""Predeclared gated echo-delay protocols."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from causal_spacetime_lab.state_change_echo_motifs import EchoMotifRecord
from causal_spacetime_lab.state_change_echo_spectrum_semantics import (
    full_transitive_return_spectrum,
    gated_echo_delay_from_spectrum,
)
from causal_spacetime_lab.state_change_response_signature import (
    EchoResponseSignature,
    response_order_sign_matrix,
)


def gated_echo_delay_for_motif(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    motif: EchoMotifRecord,
    gate_delay_rank: int,
) -> int | None:
    """Return gated echo delay for one motif.

    Gated echo is a separate predeclared protocol, not a correction to
    D_echo after observing shortcuts.
    """

    spectrum = full_transitive_return_spectrum(
        order_matrix,
        reference_chain_event_ids,
        motif.target_event_id,
        motif.emission_position,
    )
    return gated_echo_delay_from_spectrum(spectrum, gate_delay_rank)


def gated_response_signature_from_motifs(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    motifs: list[EchoMotifRecord],
    gate_delay_rank: int,
    label: str = "gated",
) -> EchoResponseSignature:
    """Build response signature under a predeclared gated-return protocol."""

    target_event_ids: list[int] = []
    delay_ranks: list[int] = []
    reachable: list[bool] = []
    for motif in motifs:
        delay = gated_echo_delay_for_motif(
            order_matrix,
            reference_chain_event_ids,
            motif,
            gate_delay_rank,
        )
        target_event_ids.append(int(motif.target_event_id))
        delay_ranks.append(int(delay) if delay is not None else -1)
        reachable.append(delay is not None)
    delays = np.asarray(delay_ranks, dtype=int)
    mask = np.asarray(reachable, dtype=bool)
    return EchoResponseSignature(
        target_event_ids=np.asarray(target_event_ids, dtype=int),
        delay_ranks=delays,
        reachable_mask=mask,
        order_sign_matrix=response_order_sign_matrix(delays, mask),
        label=label,
    )

