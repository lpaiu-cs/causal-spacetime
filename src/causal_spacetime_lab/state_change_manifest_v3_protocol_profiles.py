"""Protocol-invariant patched v3 response-profile generation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from causal_spacetime_lab.state_change_manifest_v3_protocol_execution_spec import (
    V3ProtocolExecutionFamilySpec,
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
class V3ProtocolProfileGenerationConfig:
    """Generation config for one protocol-invariant v3 response profile."""

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


def default_v3_protocol_profile_configs(
    specs: list[V3ProtocolExecutionFamilySpec],
    seed: int = 0,
) -> list[V3ProtocolProfileGenerationConfig]:
    """Return one profile config per manifest index per patched v3 family."""

    configs: list[V3ProtocolProfileGenerationConfig] = []
    for spec_index, spec in enumerate(specs):
        for manifest_index in range(spec.planned_manifest_count):
            layer_delays = [3, 5, 8, 13, 21]
            targets_per_layer = 8
            reference_count = 8
            if "retained" in spec.family_name:
                reference_count = 9
            if "gated" in spec.family_name:
                layer_delays = [4, 7, 11, 16, 22]
            if "coverage" in spec.family_name:
                targets_per_layer = 9
            if "combined" in spec.family_name:
                reference_count = 10
            if spec.family_kind in {"failed_control", "report_only"}:
                layer_delays = [5, 8]
                targets_per_layer = 3
                reference_count = 4
            family_seed = int(seed) + 10_000 * spec_index + manifest_index
            reference_set_id = f"{spec.family_name}_refs_{manifest_index}"
            configs.append(
                V3ProtocolProfileGenerationConfig(
                    family_name=spec.family_name,
                    family_kind=spec.family_kind,
                    manifest_index=manifest_index,
                    measurement_protocol=spec.measurement_protocol,
                    measurement_protocol_id=spec.measurement_protocol_id,
                    reference_set_id=reference_set_id,
                    reference_chain_ids=[
                        f"{reference_set_id}_r{column}"
                        for column in range(reference_count)
                    ],
                    reference_length=96,
                    layer_delay_ranks=layer_delays,
                    targets_per_layer=targets_per_layer,
                    seed=family_seed,
                    handoff_provenance_type=spec.handoff_provenance_type,
                )
            )
    return configs


def build_v3_protocol_response_profile(
    config: V3ProtocolProfileGenerationConfig,
) -> EchoResponseProfileWithMetadata:
    """Build a deterministic protocol-invariant patched v3 response profile."""

    rng = np.random.default_rng(config.seed)
    layer_count = len(config.layer_delay_ranks)
    target_count = layer_count * config.targets_per_layer
    target_ids = (
        np.arange(target_count, dtype=int) + 30_000 + 1000 * int(config.manifest_index)
    )
    column_count = len(config.reference_chain_ids)
    delays = np.empty((target_count, column_count), dtype=int)
    reachable = np.ones((target_count, column_count), dtype=bool)
    for layer_index, base_delay in enumerate(config.layer_delay_ranks):
        for within_layer in range(config.targets_per_layer):
            row = layer_index * config.targets_per_layer + within_layer
            target_offset = within_layer % 4
            for column in range(column_count):
                reference_offset = (column + layer_index + config.manifest_index) % 4
                delays[row, column] = int(base_delay + target_offset + reference_offset)
                if config.family_kind not in {"failed_control", "report_only"}:
                    if (
                        "coverage" not in config.family_name
                        and (row + column) % 23 == 0
                    ):
                        reachable[row, column] = False
    if config.family_kind in {"failed_control", "report_only"}:
        reachable[:, :] = False
        reachable[::2, ::2] = True
        delays[~reachable] = -1
    else:
        jitter = rng.integers(0, 2, size=delays.shape)
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
