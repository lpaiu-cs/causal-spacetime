"""Exact sanity checks for pairwise response-profile comparison protocols."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
    pairwise_response_dissimilarity,
    pairwise_response_dissimilarity_matrix,
    pairwise_response_order_inversion_rate,
    response_pair_comparison_constraints,
)
from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile

DEFAULT_OUTPUT = Path("outputs/data/pairwise_response_profile_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Pairwise response exact sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / DEFAULT_OUTPUT.name


def _profile() -> EchoResponseProfile:
    delays = np.asarray(
        [
            [1, 2, 3],
            [1, 2, 4],
            [2, 2, 3],
            [1, -1, 3],
        ],
        dtype=int,
    )
    reachable = delays >= 0
    return EchoResponseProfile(
        target_event_ids=np.asarray([10, 11, 12, 13], dtype=int),
        protocol_labels=["p0", "p1", "p2"],
        delay_rank_matrix=delays,
        reachable_matrix=reachable,
    )


def run_experiment() -> list[dict[str, float | str]]:
    """Run deterministic protocol checks."""

    profile = _profile()
    protocols = [
        PairwiseResponseComparisonProtocol("sep", "separation_fraction"),
        PairwiseResponseComparisonProtocol("gap", "rank_gap_mean"),
        PairwiseResponseComparisonProtocol("mismatch", "reachability_mismatch"),
        PairwiseResponseComparisonProtocol(
            "combined",
            "combined_gap_and_mismatch",
            missing_policy="penalize_mismatch",
        ),
    ]
    sep_matrix = pairwise_response_dissimilarity_matrix(profile, protocols[0])
    mismatch_matrix = pairwise_response_dissimilarity_matrix(profile, protocols[2])
    dissimilarity = pairwise_response_dissimilarity(profile, protocols[1])
    constraints = response_pair_comparison_constraints(
        dissimilarity,
        5,
        seed=0,
    )
    return [
        {
            "check": "separation_matrix_symmetric",
            "passed": float(np.allclose(sep_matrix, sep_matrix.T, equal_nan=True)),
            "value": str(sep_matrix),
        },
        {
            "check": "diagonal_zero",
            "passed": float(np.allclose(np.diag(sep_matrix), 0.0)),
            "value": str(np.diag(sep_matrix)),
        },
        {
            "check": "expected_separation_value",
            "passed": float(np.isclose(sep_matrix[0, 1], 1.0 / 3.0)),
            "value": str(sep_matrix[0, 1]),
        },
        {
            "check": "expected_reachability_mismatch",
            "passed": float(np.isclose(mismatch_matrix[0, 3], 1.0 / 3.0)),
            "value": str(mismatch_matrix[0, 3]),
        },
        {
            "check": "constraint_generation",
            "passed": float(constraints.shape[1] == 4 and constraints.shape[0] > 0),
            "value": str(constraints),
        },
        {
            "check": "identical_inversion_zero",
            "passed": float(
                pairwise_response_order_inversion_rate(dissimilarity, dissimilarity)
                == 0.0
            ),
            "value": str(
                pairwise_response_order_inversion_rate(dissimilarity, dissimilarity)
            ),
        },
    ]


def write_outputs(
    rows: list[dict[str, float | str]],
    path: Path = DEFAULT_OUTPUT,
) -> Path:
    """Write exact sanity CSV."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    path = parse_args()
    rows = run_experiment()
    output_path = write_outputs(rows, path)
    print(f"Wrote pairwise response profile exact sanity: {output_path}")


if __name__ == "__main__":
    main()
