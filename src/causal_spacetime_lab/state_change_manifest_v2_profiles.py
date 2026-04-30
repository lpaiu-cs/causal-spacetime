"""Deterministic v2 response-profile generation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from causal_spacetime_lab.state_change_manifest_v2_spec import V2ManifestFamilySpec
from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile


@dataclass(frozen=True)
class V2ProfileGenerationConfig:
    """Configuration for one v2 response-profile family."""

    family_name: str
    reference_length: int
    emission_positions: list[int]
    layer_delay_ranks: list[int]
    targets_per_layer: int
    protocol_column_count: int
    reachable_enrichment: bool
    rank_resolution_enrichment: bool
    seed: int


def default_v2_profile_configs(
    family_specs: list[V2ManifestFamilySpec],
    seed: int = 0,
) -> list[V2ProfileGenerationConfig]:
    """Return deterministic profile-generation configs for v2 families."""

    configs: list[V2ProfileGenerationConfig] = []
    for index, spec in enumerate(family_specs):
        family = spec.family_name
        protocol_count = 8
        layer_delays = [3, 5, 8, 13, 21]
        targets_per_layer = 8
        reachable = "coverage" in family or "combined" in family
        rank_resolution = "rank_resolution" in family
        if "more_protocol_columns" in family:
            protocol_count = 10
        if rank_resolution:
            layer_delays = [3, 6, 10, 15, 21, 28]
        if "failed_controls" in family:
            protocol_count = 4
            targets_per_layer = 3
            layer_delays = [5, 8]
        configs.append(
            V2ProfileGenerationConfig(
                family_name=family,
                reference_length=96,
                emission_positions=[8, 16, 24, 32],
                layer_delay_ranks=layer_delays,
                targets_per_layer=targets_per_layer,
                protocol_column_count=protocol_count,
                reachable_enrichment=reachable,
                rank_resolution_enrichment=rank_resolution,
                seed=int(seed) + 1000 * index,
            )
        )
    return configs


def build_v2_response_profile(
    config: V2ProfileGenerationConfig,
) -> EchoResponseProfile:
    """Build a deterministic pre-metric v2 response profile."""

    rng = np.random.default_rng(config.seed)
    layer_count = len(config.layer_delay_ranks)
    target_count = layer_count * config.targets_per_layer
    target_ids = np.arange(target_count, dtype=int) + 10_000
    delays = np.empty((target_count, config.protocol_column_count), dtype=int)
    reachable = np.ones_like(delays, dtype=bool)
    for layer_index, base_delay in enumerate(config.layer_delay_ranks):
        for within_layer in range(config.targets_per_layer):
            row = layer_index * config.targets_per_layer + within_layer
            target_offset = within_layer % 3
            for column in range(config.protocol_column_count):
                column_offset = (column + layer_index) % 3
                if config.rank_resolution_enrichment:
                    column_offset += column % 2
                delays[row, column] = int(base_delay + target_offset + column_offset)
                if not config.reachable_enrichment and (row + column) % 17 == 0:
                    reachable[row, column] = False
    if "failed_controls" in config.family_name:
        reachable[:, :] = False
        reachable[::2, ::2] = True
        delays[~reachable] = -1
    else:
        # Small deterministic jitter preserves broad layer order while avoiding
        # excessive ties across protocol columns.
        jitter = rng.integers(0, 2, size=delays.shape)
        delays = delays + jitter
        delays[~reachable] = -1
    return EchoResponseProfile(
        target_event_ids=target_ids,
        protocol_labels=[
            f"{config.family_name}_protocol_{column}"
            for column in range(config.protocol_column_count)
        ],
        delay_rank_matrix=delays.astype(int, copy=False),
        reachable_matrix=reachable.astype(bool, copy=False),
    )
