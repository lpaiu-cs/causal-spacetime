"""Exact sanity checks for response-comparison validation utilities."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_response_constraint_coverage import (
    constraint_pair_node_coverage,
    constraint_target_coverage,
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
    bootstrap_profile_protocol_columns,
    heldout_protocol_constraint_validation,
    split_profile_protocol_columns,
    validation_gate_pass_fail,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
    pairwise_response_dissimilarity,
)
from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile

DEFAULT_OUTPUT_DIR = Path("outputs")


def _profile() -> EchoResponseProfile:
    return EchoResponseProfile(
        target_event_ids=np.asarray([10, 11, 12, 13, 14], dtype=int),
        protocol_labels=["a", "b", "c", "d"],
        delay_rank_matrix=np.asarray(
            [
                [1, 1, 2, 2],
                [2, 2, 3, 3],
                [3, 3, 4, 4],
                [5, 5, 6, 6],
                [8, 8, 9, 9],
            ],
            dtype=int,
        ),
        reachable_matrix=np.ones((5, 4), dtype=bool),
    )


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact validation sanity checks."""

    profile = _profile()
    protocol = PairwiseResponseComparisonProtocol("gap", "rank_gap_mean")
    dissimilarity = pairwise_response_dissimilarity(profile, protocol)
    pool = build_constraint_pool_from_dissimilarity(
        dissimilarity,
        max_constraints=40,
        min_margin=0.0,
        seed=0,
    )
    train, test = split_profile_protocol_columns(profile, 0.5, seed=0)
    heldout = heldout_protocol_constraint_validation(
        profile,
        protocol,
        max_constraints=40,
        seed=0,
    )
    bootstrap_profile = bootstrap_profile_protocol_columns(
        profile,
        sample_fraction=0.5,
        seed=0,
    )
    bootstrap = bootstrap_constraint_stability(
        profile,
        protocol,
        pool,
        bootstrap_count=3,
        seed=0,
    )
    nulls = evaluate_constraint_pool_against_nulls(
        profile,
        protocol,
        pool,
        null_repetitions=2,
        seed=0,
    )
    target_coverage = constraint_target_coverage(pool, target_count=5)
    pair_coverage = constraint_pair_node_coverage(pool, target_count=5)
    passing_gate = validation_gate_pass_fail(
        {
            "constraint_count": 100.0,
            "evaluable_fraction": 1.0,
            "agreement_fraction": 0.9,
            "inversion_fraction": 0.0,
            "tie_or_unresolved_fraction": 0.1,
        },
        ConstraintValidationGate("pass", min_constraint_count=10),
    )
    failing_gate = validation_gate_pass_fail(
        {
            "constraint_count": 1.0,
            "evaluable_fraction": 0.1,
            "agreement_fraction": 0.2,
            "inversion_fraction": 0.8,
            "tie_or_unresolved_fraction": 0.8,
        },
        ConstraintValidationGate("fail", min_constraint_count=10),
    )
    checks = [
        (
            "split_preserves_target_ids",
            np.array_equal(train.target_event_ids, test.target_event_ids),
            train.target_event_ids.tolist(),
        ),
        (
            "heldout_expected_keys",
            "agreement_fraction" in heldout and "train_protocol_count" in heldout,
            sorted(heldout),
        ),
        (
            "bootstrap_shape",
            bootstrap_profile.delay_rank_matrix.shape == (5, 2),
            bootstrap_profile.delay_rank_matrix.shape,
        ),
        (
            "bootstrap_expected_keys",
            "stable_constraint_fraction" in bootstrap,
            sorted(bootstrap),
        ),
        (
            "null_expected_keys",
            "null_z_score" in nulls,
            sorted(nulls),
        ),
        (
            "target_coverage_fields",
            target_coverage["touched_target_count"] > 0,
            target_coverage,
        ),
        (
            "pair_node_coverage_fields",
            pair_coverage["touched_pair_node_count"] > 0,
            pair_coverage,
        ),
        (
            "gate_pass_fail",
            passing_gate["passed"] == 1.0 and failing_gate["passed"] == 0.0,
            f"{passing_gate['passed']}/{failing_gate['passed']}",
        ),
    ]
    return [
        {"check": name, "passed": float(passed), "value": str(value)}
        for name, passed, value in checks
    ]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write exact validation sanity CSV."""

    path = output_dir / "data" / "response_constraint_validation_exact_sanity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote response constraint validation exact sanity: {output_path}")


if __name__ == "__main__":
    main()
