"""Protocol-invariant v4 response-profile generation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from causal_spacetime_lab.state_change_manifest_v4_execution_spec import (
    V4ProtocolExecutionFamilySpec,
)
from causal_spacetime_lab.state_change_manifest_v4_protocol_mapping import (
    measurement_protocol_for_v4_family,
    v4_protocol_id,
)
from causal_spacetime_lab.state_change_measurement_protocol import (
    MeasurementProtocolSpec,
)
from causal_spacetime_lab.state_change_response_profile_metadata import (
    profile_metadata_from_protocols,
)
from causal_spacetime_lab.state_change_response_profiles import (
    EchoResponseProfile,
    EchoResponseProfileWithMetadata,
    attach_profile_metadata,
)


@dataclass(frozen=True)
class V4ProtocolProfileGenerationConfig:
    """Generation config for one protocol-invariant v4 response profile."""

    family_name: str
    family_kind: str
    manifest_index: int
    measurement_protocol: MeasurementProtocolSpec
    measurement_protocol_id: str
    reference_set_id: str
    reference_chain_ids: list[str]
    reference_length: int
    layer_delay_ranks: list[int]
    targets_per_layer: int
    seed: int
    handoff_provenance_type: str
    profile_resolution_policy: str
    target_inclusion_policy: str
    linked_v3_root_causes: list[str]


def default_v4_protocol_profile_configs(
    specs: list[V4ProtocolExecutionFamilySpec],
    seed: int = 0,
) -> list[V4ProtocolProfileGenerationConfig]:
    """Return one profile config per planned v4 manifest."""

    configs: list[V4ProtocolProfileGenerationConfig] = []
    for spec_index, spec in enumerate(specs):
        protocol = measurement_protocol_for_v4_family(spec)
        for manifest_index in range(spec.planned_manifest_count):
            layer_delays = [3, 5, 8, 12, 17, 23]
            targets_per_layer = 10
            reference_count = 10
            if "retained" in spec.family_name:
                layer_delays = [3, 6, 10, 15, 21, 28]
                reference_count = 12
                targets_per_layer = 11
            if "gated" in spec.family_name:
                layer_delays = [4, 7, 11, 16, 22, 29]
            if "high_resolution" in spec.family_name:
                layer_delays = [2, 4, 7, 11, 16, 22, 29, 37]
                reference_count = 14
                targets_per_layer = 12
            if "combined" in spec.family_name:
                reference_count = 12
                targets_per_layer = 10
            if spec.family_kind in {"failed_control", "report_only"}:
                layer_delays = [4, 7]
                targets_per_layer = 3
                reference_count = 4
            family_seed = int(seed) + 20_000 * spec_index + manifest_index
            reference_set_id = f"{spec.family_name}_refs_{manifest_index}"
            configs.append(
                V4ProtocolProfileGenerationConfig(
                    family_name=spec.family_name,
                    family_kind=spec.family_kind,
                    manifest_index=manifest_index,
                    measurement_protocol=protocol,
                    measurement_protocol_id=v4_protocol_id(protocol),
                    reference_set_id=reference_set_id,
                    reference_chain_ids=[
                        f"{reference_set_id}_r{index}"
                        for index in range(reference_count)
                    ],
                    reference_length=128,
                    layer_delay_ranks=layer_delays,
                    targets_per_layer=targets_per_layer,
                    seed=family_seed,
                    handoff_provenance_type=spec.handoff_provenance_type,
                    profile_resolution_policy=spec.profile_resolution_policy,
                    target_inclusion_policy=spec.target_inclusion_policy,
                    linked_v3_root_causes=list(spec.linked_v3_root_causes),
                )
            )
    return configs


def build_v4_protocol_response_profile(
    config: V4ProtocolProfileGenerationConfig,
) -> EchoResponseProfileWithMetadata:
    """Build a deterministic protocol-invariant v4 response profile."""

    rng = np.random.default_rng(config.seed)
    layer_count = len(config.layer_delay_ranks)
    target_count = layer_count * config.targets_per_layer
    target_ids = (
        np.arange(target_count, dtype=int)
        + 40_000
        + 1000 * int(config.manifest_index)
    )
    column_count = len(config.reference_chain_ids)
    delays = np.empty((target_count, column_count), dtype=int)
    reachable = np.ones((target_count, column_count), dtype=bool)
    high_resolution = "high_resolution" in config.family_name
    for layer_index, base_delay in enumerate(config.layer_delay_ranks):
        for within_layer in range(config.targets_per_layer):
            row = layer_index * config.targets_per_layer + within_layer
            target_offset = within_layer % (5 if high_resolution else 4)
            for column in range(column_count):
                reference_offset = (
                    column + layer_index + config.manifest_index
                ) % (5 if high_resolution else 4)
                delays[row, column] = int(base_delay + target_offset + reference_offset)
                if config.family_kind not in {"failed_control", "report_only"}:
                    if "coverage" not in config.target_inclusion_policy:
                        if (row + column + config.manifest_index) % 31 == 0:
                            reachable[row, column] = False
    if config.family_kind in {"failed_control", "report_only"}:
        reachable[:, :] = False
        reachable[::2, ::2] = True
        delays[~reachable] = -1
    else:
        jitter_upper = 1 if "stability" in config.family_name else 2
        jitter = rng.integers(0, jitter_upper + 1, size=delays.shape)
        delays = delays + jitter
        delays[~reachable] = -1
    profile = EchoResponseProfile(
        target_event_ids=target_ids,
        protocol_labels=list(config.reference_chain_ids),
        delay_rank_matrix=delays.astype(int, copy=False),
        reachable_matrix=reachable.astype(bool, copy=False),
    )
    metadata = profile_metadata_from_protocols(
        profile_family_id=config.family_name,
        protocols=[config.measurement_protocol for _ in config.reference_chain_ids],
        reference_chain_ids=config.reference_chain_ids,
        reference_set_id=config.reference_set_id,
    )
    return attach_profile_metadata(profile, metadata)
