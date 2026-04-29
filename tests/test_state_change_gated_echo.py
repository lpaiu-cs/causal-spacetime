from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    insert_echo_motif,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_gated_echo import (
    gated_echo_delay_for_motif,
    gated_response_signature_from_motifs,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)


def test_gated_echo_delay_for_motif() -> None:
    network, reference = generate_reference_backbone_network(10)
    network, motif = insert_echo_motif(
        network,
        reference,
        EchoMotifSpec(2, 4, return_steps=0),
    )
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))

    assert gated_echo_delay_for_motif(closure, reference, motif, 3) == 4
    assert gated_echo_delay_for_motif(closure, reference, motif, 5) == 5
    assert gated_echo_delay_for_motif(closure, reference, motif, 8) is None


def test_gated_response_signature_from_motifs() -> None:
    network, reference = generate_reference_backbone_network(10)
    network, motif = insert_echo_motif(
        network,
        reference,
        EchoMotifSpec(2, 4, return_steps=0),
    )
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))

    signature = gated_response_signature_from_motifs(
        closure,
        reference,
        [motif],
        gate_delay_rank=3,
    )

    assert np.array_equal(signature.delay_ranks, np.asarray([4]))
    assert np.array_equal(signature.reachable_mask, np.asarray([True]))
