from __future__ import annotations

import pytest

from causal_spacetime_lab.state_change_echo_interference import (
    return_spectrum_report_for_motif,
)
from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    insert_echo_motif,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_echo_shortcuts import (
    ShortcutInjectionSpec,
    add_random_acyclic_background_edges,
    classify_added_edge_effects_on_motifs,
    inject_shortcut_returns,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    topological_order_from_adjacency,
    transitive_closure_dag,
)


def _clean_motif_case():
    network, reference = generate_reference_backbone_network(8)
    network, motif = insert_echo_motif(
        network,
        reference,
        EchoMotifSpec(2, 3, outward_steps=1, return_steps=1),
    )
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    return network, closure, reference, motif


def test_shortcut_injection_spec_validation() -> None:
    with pytest.raises(ValueError):
        ShortcutInjectionSpec(probability=-0.1)
    with pytest.raises(ValueError):
        ShortcutInjectionSpec(probability=1.1)
    with pytest.raises(ValueError):
        ShortcutInjectionSpec(probability=0.1, min_depth=0)
    with pytest.raises(ValueError):
        ShortcutInjectionSpec(probability=0.1, min_depth=2, max_depth=1)
    with pytest.raises(ValueError):
        ShortcutInjectionSpec(probability=0.1, mode="not_a_mode")


def test_inject_shortcut_returns_adds_expected_shortcut() -> None:
    network, _, reference, motif = _clean_motif_case()

    updated, records = inject_shortcut_returns(
        network,
        reference,
        [motif],
        ShortcutInjectionSpec(probability=1.0, min_depth=1, max_depth=1),
        seed=0,
    )
    closure = transitive_closure_dag(immediate_trigger_adjacency(updated))
    report = return_spectrum_report_for_motif(closure, reference, motif)

    assert len(records) == 1
    assert records[0].shortcut_delay_rank == 2
    assert report["recovered_delay_rank"] == 2.0
    assert report["early_shortcut"] == 1.0


def test_inject_shortcut_returns_preserves_acyclicity() -> None:
    network, _, reference, motif = _clean_motif_case()

    updated, _ = inject_shortcut_returns(
        network,
        reference,
        [motif],
        ShortcutInjectionSpec(probability=1.0, mode="decoy_path_to_early_reference"),
        seed=0,
    )

    order = topological_order_from_adjacency(immediate_trigger_adjacency(updated))
    assert order.size == len(updated.events)


def test_add_random_acyclic_background_edges_preserves_acyclicity() -> None:
    network, _, _, _ = _clean_motif_case()

    updated, added = add_random_acyclic_background_edges(
        network,
        edge_probability=1.0,
        seed=0,
        max_edges=5,
    )

    order = topological_order_from_adjacency(immediate_trigger_adjacency(updated))
    assert order.size == len(updated.events)
    assert 0 <= added <= 5


def test_classify_added_edge_effects_on_motifs_detects_new_shortcut() -> None:
    network, before, reference, motif = _clean_motif_case()
    updated, _ = inject_shortcut_returns(
        network,
        reference,
        [motif],
        ShortcutInjectionSpec(probability=1.0, min_depth=1, max_depth=1),
        seed=0,
    )
    after = transitive_closure_dag(immediate_trigger_adjacency(updated))

    rows = classify_added_edge_effects_on_motifs(before, after, reference, [motif])

    assert rows[0]["effect"] == "new_shortcut"
    assert rows[0]["after_delay_rank"] == 2.0


def test_classify_added_edge_effects_handles_unchanged_case() -> None:
    _, before, reference, motif = _clean_motif_case()

    rows = classify_added_edge_effects_on_motifs(before, before, reference, [motif])

    assert rows[0]["effect"] == "unchanged"
    assert rows[0]["before_delay_rank"] == rows[0]["after_delay_rank"]
