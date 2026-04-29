"""Pipeline for building response-comparison handoff manifests."""

from __future__ import annotations

from dataclasses import replace

from causal_spacetime_lab.state_change_response_constraint_coverage import (
    constraint_pool_summary,
)
from causal_spacetime_lab.state_change_response_constraint_null_validation import (
    evaluate_constraint_pool_against_nulls,
)
from causal_spacetime_lab.state_change_response_constraint_pool import (
    build_constraint_pool_from_dissimilarity,
)
from causal_spacetime_lab.state_change_response_constraint_validation import (
    ConstraintValidationGate,
    bootstrap_constraint_stability,
    heldout_protocol_constraint_validation,
)
from causal_spacetime_lab.state_change_response_handoff import (
    ResponseConstraintHandoffManifest,
    build_handoff_validation_summary,
    decide_handoff_eligibility,
    forbidden_interpretations_default,
    manifest_digest,
    manifest_to_jsonable,
    split_constraint_indices,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
    pairwise_response_dissimilarity,
)
from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile


def build_candidate_handoff_manifest(
    profile: EchoResponseProfile,
    comparison_protocol: PairwiseResponseComparisonProtocol,
    gate: ConstraintValidationGate,
    *,
    max_constraints: int = 5000,
    min_margin: float = 0.05,
    train_fraction: float = 0.7,
    constraint_seed: int = 0,
    bootstrap_count: int = 20,
    null_repetitions: int = 5,
    source_label: str = "candidate",
) -> ResponseConstraintHandoffManifest:
    """Build a predeclared handoff manifest without fitting an embedding."""

    dissimilarity = pairwise_response_dissimilarity(profile, comparison_protocol)
    pool = build_constraint_pool_from_dissimilarity(
        dissimilarity,
        max_constraints=max_constraints,
        min_margin=min_margin,
        seed=constraint_seed,
        source_label=source_label,
    )
    coverage = constraint_pool_summary(pool, int(profile.target_event_ids.size))
    heldout = heldout_protocol_constraint_validation(
        profile,
        comparison_protocol,
        train_fraction=train_fraction,
        max_constraints=max_constraints,
        min_margin=min_margin,
        seed=constraint_seed,
    )
    bootstrap = bootstrap_constraint_stability(
        profile,
        comparison_protocol,
        pool,
        bootstrap_count=bootstrap_count,
        seed=constraint_seed,
    )
    nulls = evaluate_constraint_pool_against_nulls(
        profile,
        comparison_protocol,
        pool,
        null_repetitions=null_repetitions,
        seed=constraint_seed,
    )
    summary = build_handoff_validation_summary(
        heldout,
        bootstrap,
        nulls,
        coverage,
        int(pool.constraints.shape[0]),
    )
    decision = decide_handoff_eligibility(summary, gate)
    train_indices, heldout_indices = split_constraint_indices(
        int(pool.constraints.shape[0]),
        train_fraction,
        seed=constraint_seed,
    )
    manifest = ResponseConstraintHandoffManifest(
        manifest_id="",
        created_by_milestone="31",
        profile_label=source_label,
        comparison_protocol_name=comparison_protocol.name,
        comparison_method=comparison_protocol.method,
        missing_policy=comparison_protocol.missing_policy,
        min_common_protocols=comparison_protocol.min_common_protocols,
        min_margin=float(min_margin),
        max_constraints=int(max_constraints),
        constraint_seed=int(constraint_seed),
        train_fraction=float(train_fraction),
        validation_gate_name=gate.name,
        validation_summary=summary,
        handoff_decision=decision,
        target_event_ids=pool.target_event_ids.copy(),
        constraints=pool.constraints.copy(),
        margins=pool.margins.copy(),
        train_constraint_indices=train_indices,
        heldout_constraint_indices=heldout_indices,
        null_baseline_labels=[
            "shuffle_delays",
            "shuffle_reachability",
            "permute_profiles",
            "random_same_marginals",
        ],
        forbidden_interpretations=forbidden_interpretations_default(),
    )
    digest = manifest_digest(manifest_to_jsonable(manifest))
    return replace(manifest, manifest_id=digest)


def select_eligible_manifests(
    manifests: list[ResponseConstraintHandoffManifest],
) -> list[ResponseConstraintHandoffManifest]:
    """Return eligible handoff manifests."""

    return [manifest for manifest in manifests if manifest.handoff_decision.eligible]


def summarize_handoff_manifests(
    manifests: list[ResponseConstraintHandoffManifest],
) -> list[dict[str, float | str]]:
    """Summarize handoff manifests for CSV output."""

    rows: list[dict[str, float | str]] = []
    for manifest in manifests:
        summary = manifest.validation_summary
        rows.append(
            {
                "manifest_id": manifest.manifest_id,
                "protocol_name": manifest.comparison_protocol_name,
                "method": manifest.comparison_method,
                "min_margin": float(manifest.min_margin),
                "eligible": float(manifest.handoff_decision.eligible),
                "failed_reasons": ";".join(manifest.handoff_decision.failed_reasons),
                "constraint_count": float(summary.constraint_count),
                "heldout_agreement": float(summary.heldout_agreement_fraction),
                "bootstrap_agreement": float(
                    summary.bootstrap_mean_agreement_fraction
                ),
                "null_z_score": float(summary.null_z_score),
                "target_coverage": float(summary.target_coverage_fraction),
                "pair_node_coverage": float(summary.pair_node_coverage_fraction),
            }
        )
    return rows
