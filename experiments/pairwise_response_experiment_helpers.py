"""Shared helpers for pairwise response-profile experiments."""

from __future__ import annotations

import numpy as np
from exp115_multi_reference_response_profile_diagnostics import (
    ExperimentConfig as ProfileConfig,
)
from exp115_multi_reference_response_profile_diagnostics import (
    build_synthetic_protocol_signatures,
)

from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
    pairwise_response_dissimilarity,
    response_pair_comparison_constraints,
)
from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile
from causal_spacetime_lab.state_change_response_signature import (
    EchoResponseSignature,
    response_order_sign_matrix,
)


def default_pairwise_protocols() -> list[PairwiseResponseComparisonProtocol]:
    """Return the default admissible pairwise response-comparison protocols."""

    return [
        PairwiseResponseComparisonProtocol("separation", "separation_fraction"),
        PairwiseResponseComparisonProtocol("gap_mean", "rank_gap_mean"),
        PairwiseResponseComparisonProtocol("gap_median", "rank_gap_median"),
        PairwiseResponseComparisonProtocol("mismatch", "reachability_mismatch"),
        PairwiseResponseComparisonProtocol(
            "combined",
            "combined_gap_and_mismatch",
            missing_policy="penalize_mismatch",
        ),
    ]


def profile_from_synthetic_config(
    reference_length: int,
    emission_positions: tuple[int, ...],
    layer_delay_ranks: tuple[int, ...],
    targets_per_layer: int,
    repetition: int,
    seed: int,
) -> EchoResponseProfile:
    """Build a synthetic multi-protocol response profile."""

    from causal_spacetime_lab.state_change_response_profiles import (
        response_profile_from_signatures,
    )

    config = ProfileConfig(
        reference_length=reference_length,
        emission_positions=emission_positions,
        layer_delay_ranks=layer_delay_ranks,
        targets_per_layer=targets_per_layer,
        repetitions=1,
        seed=seed,
    )
    return response_profile_from_signatures(
        build_synthetic_protocol_signatures(config, repetition)
    )


def synthetic_random_profile(
    target_count: int,
    protocol_count: int,
    reachable_probability: float,
    unique_rank_count: int,
    seed: int,
) -> EchoResponseProfile:
    """Generate a synthetic response profile with controlled reachability."""

    if not 0.0 <= reachable_probability <= 1.0:
        raise ValueError("reachable_probability must be in [0, 1]")
    rng = np.random.default_rng(seed)
    reachable = rng.random((target_count, protocol_count)) < reachable_probability
    delays = rng.integers(1, unique_rank_count + 1, size=(target_count, protocol_count))
    delays = delays.astype(int)
    delays[~reachable] = -1
    return EchoResponseProfile(
        target_event_ids=np.arange(target_count, dtype=int),
        protocol_labels=[f"p{index}" for index in range(protocol_count)],
        delay_rank_matrix=delays,
        reachable_matrix=reachable,
    )


def signature_from_profile_column(
    profile: EchoResponseProfile,
    column: int,
    label: str,
) -> EchoResponseSignature:
    """Build a response signature from one profile column."""

    delays = profile.delay_rank_matrix[:, column]
    reachable = profile.reachable_matrix[:, column]
    return EchoResponseSignature(
        target_event_ids=profile.target_event_ids,
        delay_ranks=delays,
        reachable_mask=reachable,
        order_sign_matrix=response_order_sign_matrix(delays, reachable),
        label=label,
    )


def constraint_count(
    profile: EchoResponseProfile,
    protocol: PairwiseResponseComparisonProtocol,
    seed: int,
) -> float:
    """Return the available sampled response-comparison constraint count."""

    dissimilarity = pairwise_response_dissimilarity(profile, protocol)
    constraints = response_pair_comparison_constraints(
        dissimilarity,
        num_constraints=500,
        seed=seed,
    )
    return float(constraints.shape[0])
