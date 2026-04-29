"""Protocol-variant helpers for echo-response order signatures."""

from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change import StateChangeNetwork
from causal_spacetime_lab.state_change_echo import echo_delay_rank_for_emission
from causal_spacetime_lab.state_change_echo_coarse_graining import (
    coarse_emission_position_for_motif,
    protected_indices_for_reference_and_motifs,
    remap_echo_motif_record_for_event_thinning,
    remap_reference_chain,
    restrict_transitive_order_to_retained_events,
    sample_retained_indices,
    subsample_reference_chain_positions,
)
from causal_spacetime_lab.state_change_echo_motif_validation import (
    recovered_delay_for_motif,
)
from causal_spacetime_lab.state_change_echo_motifs import EchoMotifRecord
from causal_spacetime_lab.state_change_echo_shortcuts import (
    ShortcutInjectionSpec,
    inject_shortcut_returns,
)
from causal_spacetime_lab.state_change_edge_coarse_graining import (
    protected_edge_keys_for_motifs,
    protected_source_target_pairs_for_reference_chain,
    thin_immediate_trigger_edges,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)
from causal_spacetime_lab.state_change_response_signature import (
    EchoResponseSignature,
    echo_response_signature_from_motifs,
    response_order_sign_matrix,
)


def _closure(network: StateChangeNetwork) -> np.ndarray:
    return transitive_closure_dag(immediate_trigger_adjacency(network))


def _signature_from_targets_delays(
    target_event_ids: list[int],
    delay_values: list[int | None],
    label: str,
) -> EchoResponseSignature:
    targets = np.asarray(target_event_ids, dtype=int)
    reachable = np.asarray([value is not None for value in delay_values], dtype=bool)
    delays = np.asarray(
        [int(value) if value is not None else -1 for value in delay_values],
        dtype=int,
    )
    return EchoResponseSignature(
        target_event_ids=targets,
        delay_ranks=delays,
        reachable_mask=reachable,
        order_sign_matrix=response_order_sign_matrix(delays, reachable),
        label=label,
    )


def signature_after_closure_event_thinning(
    baseline_network: StateChangeNetwork,
    reference_chain_event_ids: np.ndarray,
    motifs: list[EchoMotifRecord],
    *,
    keep_probability: float,
    seed: int | None = None,
    label: str = "closure_event_thinning",
) -> EchoResponseSignature:
    """Build a signature after closure-preserving event thinning."""

    closure = _closure(baseline_network)
    protected = protected_indices_for_reference_and_motifs(
        reference_chain_event_ids,
        motifs,
    )
    retained = sample_retained_indices(
        closure.shape[0],
        keep_probability,
        protected,
        seed=seed,
    )
    thinning = restrict_transitive_order_to_retained_events(closure, retained)
    coarse_reference = remap_reference_chain(
        reference_chain_event_ids,
        thinning.old_to_new,
    )
    targets: list[int] = []
    delays: list[int | None] = []
    for motif in motifs:
        targets.append(int(motif.target_event_id))
        remapped = remap_echo_motif_record_for_event_thinning(
            motif,
            thinning.old_to_new,
        )
        if remapped is None:
            delays.append(None)
            continue
        delays.append(
            recovered_delay_for_motif(
                thinning.restricted_order_matrix,
                coarse_reference,
                remapped,
            )
        )
    return _signature_from_targets_delays(targets, delays, label)


def signature_after_reference_subsampling(
    baseline_network: StateChangeNetwork,
    reference_chain_event_ids: np.ndarray,
    motifs: list[EchoMotifRecord],
    *,
    stride: int,
    label: str = "reference_subsampling",
) -> EchoResponseSignature:
    """Build a signature using a subsampled reference chain."""

    closure = _closure(baseline_network)
    protected_positions = np.asarray(
        [int(motif.emission_position) for motif in motifs],
        dtype=int,
    )
    subsampling = subsample_reference_chain_positions(
        reference_chain_event_ids,
        stride,
        protected_positions=protected_positions,
    )
    targets: list[int] = []
    delays: list[int | None] = []
    for motif in motifs:
        targets.append(int(motif.target_event_id))
        coarse_emission = coarse_emission_position_for_motif(motif, subsampling)
        if coarse_emission is None:
            delays.append(None)
            continue
        delays.append(
            echo_delay_rank_for_emission(
                closure,
                subsampling.subsampled_reference_chain,
                motif.target_event_id,
                coarse_emission,
            )
        )
    return _signature_from_targets_delays(targets, delays, label)


def signature_after_edge_thinning(
    baseline_network: StateChangeNetwork,
    reference_chain_event_ids: np.ndarray,
    motifs: list[EchoMotifRecord],
    *,
    removal_probability: float,
    protect_motif_edges: bool,
    seed: int | None = None,
    label: str = "edge_thinning",
) -> EchoResponseSignature:
    """Build a signature after immediate-edge thinning."""

    protected_pairs = protected_source_target_pairs_for_reference_chain(
        reference_chain_event_ids
    )
    if protect_motif_edges:
        protected_pairs |= protected_edge_keys_for_motifs(motifs)
    thinned_network, _ = thin_immediate_trigger_edges(
        baseline_network,
        removal_probability,
        protected_source_target_pairs=protected_pairs,
        seed=seed,
    )
    return echo_response_signature_from_motifs(
        _closure(thinned_network),
        reference_chain_event_ids,
        motifs,
        label=label,
    )


def signature_after_shortcut_injection(
    baseline_network: StateChangeNetwork,
    reference_chain_event_ids: np.ndarray,
    motifs: list[EchoMotifRecord],
    *,
    probability: float,
    seed: int | None = None,
    label: str = "shortcut_injection",
) -> EchoResponseSignature:
    """Build a signature after controlled shortcut-return injection."""

    network, _ = inject_shortcut_returns(
        baseline_network,
        reference_chain_event_ids,
        motifs,
        ShortcutInjectionSpec(probability=probability),
        seed=seed,
    )
    return echo_response_signature_from_motifs(
        _closure(network),
        reference_chain_event_ids,
        motifs,
        label=label,
    )


def build_protocol_variant_signatures(
    baseline_network: StateChangeNetwork,
    reference_chain_event_ids: np.ndarray,
    motifs: list[EchoMotifRecord],
    variant_specs: list[dict[str, object]],
) -> list[EchoResponseSignature]:
    """Build response signatures for simple shared protocol variants."""

    signatures: list[EchoResponseSignature] = []
    baseline_closure = _closure(baseline_network)
    for index, spec in enumerate(variant_specs):
        kind = str(spec.get("kind", "baseline"))
        label = str(spec.get("label", kind))
        seed_value = spec.get("seed")
        seed = int(seed_value) if seed_value is not None else index
        if kind == "baseline":
            signatures.append(
                echo_response_signature_from_motifs(
                    baseline_closure,
                    reference_chain_event_ids,
                    motifs,
                    label=label,
                )
            )
        elif kind == "closure_event_thinning":
            signatures.append(
                signature_after_closure_event_thinning(
                    baseline_network,
                    reference_chain_event_ids,
                    motifs,
                    keep_probability=float(spec["keep_probability"]),
                    seed=seed,
                    label=label,
                )
            )
        elif kind == "reference_subsampling":
            signatures.append(
                signature_after_reference_subsampling(
                    baseline_network,
                    reference_chain_event_ids,
                    motifs,
                    stride=int(spec["stride"]),
                    label=label,
                )
            )
        elif kind == "edge_thinning":
            signatures.append(
                signature_after_edge_thinning(
                    baseline_network,
                    reference_chain_event_ids,
                    motifs,
                    removal_probability=float(spec["removal_probability"]),
                    protect_motif_edges=bool(spec["protect_motif_edges"]),
                    seed=seed,
                    label=label,
                )
            )
        elif kind == "shortcut_injection":
            signatures.append(
                signature_after_shortcut_injection(
                    baseline_network,
                    reference_chain_event_ids,
                    motifs,
                    probability=float(spec["probability"]),
                    seed=seed,
                    label=label,
                )
            )
        else:
            raise ValueError(f"unknown protocol variant kind: {kind}")
    return signatures

