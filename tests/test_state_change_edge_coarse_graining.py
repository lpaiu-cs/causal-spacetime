from __future__ import annotations

from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    insert_echo_motif,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_edge_coarse_graining import (
    compare_recovery_before_after_edge_thinning,
    edge_key,
    protected_edge_keys_for_motifs,
    protected_edge_keys_for_reference_chain,
    protected_source_target_pairs_for_reference_chain,
    thin_immediate_trigger_edges,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    topological_order_from_adjacency,
    transitive_closure_dag,
)


def _motif_case():
    network, reference = generate_reference_backbone_network(8)
    network, motif = insert_echo_motif(
        network,
        reference,
        EchoMotifSpec(2, 3, outward_steps=1, return_steps=1),
    )
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    return network, closure, reference, motif


def test_edge_key() -> None:
    network, _, _, _ = _motif_case()

    key = edge_key(network.trigger_edges[0])

    assert len(key) == 3
    assert key[0] == network.trigger_edges[0].source_event_id


def test_protected_reference_chain_pairs_and_keys() -> None:
    _, _, reference, _ = _motif_case()

    pairs = protected_source_target_pairs_for_reference_chain(reference)
    keys = protected_edge_keys_for_reference_chain(reference)

    assert (int(reference[0]), int(reference[1])) in pairs
    assert (int(reference[0]), int(reference[1]), "local_successor") in keys


def test_protected_edge_keys_for_motifs() -> None:
    _, _, reference, motif = _motif_case()

    pairs = protected_edge_keys_for_motifs([motif])

    assert (int(reference[2]), int(motif.outward_event_ids[0])) in pairs
    assert (int(motif.return_event_ids[0]), int(reference[5])) in pairs


def test_thin_immediate_trigger_edges_preserves_acyclicity() -> None:
    network, _, reference, _ = _motif_case()
    protected = protected_source_target_pairs_for_reference_chain(reference)

    thinned, removed = thin_immediate_trigger_edges(
        network,
        removal_probability=1.0,
        protected_source_target_pairs=protected,
        seed=0,
    )

    order = topological_order_from_adjacency(immediate_trigger_adjacency(thinned))
    assert order.size == len(thinned.events)
    assert removed > 0


def test_compare_recovery_before_after_edge_thinning_detects_missing() -> None:
    network, before, reference, motif = _motif_case()
    thinned, _ = thin_immediate_trigger_edges(network, 1.0, seed=0)
    after = transitive_closure_dag(immediate_trigger_adjacency(thinned))

    rows = compare_recovery_before_after_edge_thinning(
        before,
        after,
        reference,
        [motif],
    )

    assert rows[0]["effect"] == "became_missing"
    assert rows[0]["before_delay_rank"] == 3.0
