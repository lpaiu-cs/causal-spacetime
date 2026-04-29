"""Exact sanity checks for pairwise null baselines and admissibility."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
)
from causal_spacetime_lab.state_change_response_pairwise_admissibility import (
    compare_pairwise_protocols,
    pairwise_protocol_admissibility_report,
)
from causal_spacetime_lab.state_change_response_pairwise_nulls import (
    permute_target_profiles,
    random_profile_with_same_marginals,
    shuffle_profile_delays_within_protocols,
    shuffle_profile_reachability_masks,
)
from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile

DEFAULT_OUTPUT = Path(
    "outputs/data/pairwise_response_null_admissibility_exact_sanity.csv"
)


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Pairwise response null exact sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / DEFAULT_OUTPUT.name


def _profile() -> EchoResponseProfile:
    delays = np.asarray(
        [
            [1, 2, 3],
            [2, -1, 3],
            [3, 2, -1],
            [1, 4, 5],
        ],
        dtype=int,
    )
    return EchoResponseProfile(
        target_event_ids=np.asarray([10, 11, 12, 13], dtype=int),
        protocol_labels=["p0", "p1", "p2"],
        delay_rank_matrix=delays,
        reachable_matrix=delays >= 0,
    )


def _column_delay_multisets(profile: EchoResponseProfile) -> list[list[int]]:
    return [
        sorted(
            int(value)
            for value in profile.delay_rank_matrix[
                profile.reachable_matrix[:, column],
                column,
            ]
        )
        for column in range(profile.delay_rank_matrix.shape[1])
    ]


def _reachable_counts(profile: EchoResponseProfile) -> list[int]:
    return [
        int(np.sum(profile.reachable_matrix[:, column]))
        for column in range(profile.reachable_matrix.shape[1])
    ]


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact null/admissibility checks."""

    profile = _profile()
    nulls = [
        shuffle_profile_delays_within_protocols(profile, seed=0),
        shuffle_profile_reachability_masks(profile, seed=0),
        permute_target_profiles(profile, seed=0),
        random_profile_with_same_marginals(profile, seed=0),
    ]
    protocol = PairwiseResponseComparisonProtocol("gap", "rank_gap_mean")
    report = pairwise_protocol_admissibility_report(profile, protocol)
    comparisons = compare_pairwise_protocols(
        profile,
        [
            protocol,
            PairwiseResponseComparisonProtocol("sep", "separation_fraction"),
        ],
    )
    baseline_multisets = _column_delay_multisets(profile)
    baseline_counts = _reachable_counts(profile)
    return [
        {
            "check": "delay_shuffle_preserves_delay_multisets",
            "passed": float(_column_delay_multisets(nulls[0]) == baseline_multisets),
            "value": str(_column_delay_multisets(nulls[0])),
        },
        {
            "check": "reachability_shuffle_preserves_counts",
            "passed": float(_reachable_counts(nulls[1]) == baseline_counts),
            "value": str(_reachable_counts(nulls[1])),
        },
        {
            "check": "permutation_preserves_column_marginals",
            "passed": float(
                _column_delay_multisets(nulls[2]) == baseline_multisets
                and _reachable_counts(nulls[2]) == baseline_counts
            ),
            "value": str(_column_delay_multisets(nulls[2])),
        },
        {
            "check": "random_same_marginals_preserves_marginals",
            "passed": float(
                _column_delay_multisets(nulls[3]) == baseline_multisets
                and _reachable_counts(nulls[3]) == baseline_counts
            ),
            "value": str(_column_delay_multisets(nulls[3])),
        },
        {
            "check": "admissibility_report_fields",
            "passed": float("valid_pair_fraction" in report and "symmetric" in report),
            "value": str(report),
        },
        {
            "check": "compare_pairwise_protocols_rows",
            "passed": float(len(comparisons) == 1),
            "value": str(comparisons),
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
    print(f"Wrote pairwise response null/admissibility exact sanity: {output_path}")


if __name__ == "__main__":
    main()
