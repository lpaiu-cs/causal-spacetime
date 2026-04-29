"""Exact sanity checks for response profiles and requirement ladder."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_representability_requirements import (
    default_response_representability_ladder,
)
from causal_spacetime_lab.state_change_response_profiles import (
    compare_response_profiles,
    response_profile_equivalence_classes,
    response_profile_from_signatures,
    response_profile_separation_fraction,
)
from causal_spacetime_lab.state_change_response_signature import (
    EchoResponseSignature,
    response_order_sign_matrix,
)

DEFAULT_OUTPUT = Path("outputs/data/response_profile_requirement_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Response profile exact sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / DEFAULT_OUTPUT.name


def _signature(
    targets: list[int],
    delays: list[int],
    label: str,
) -> EchoResponseSignature:
    target_array = np.asarray(targets, dtype=int)
    delay_array = np.asarray(delays, dtype=int)
    reachable = delay_array >= 0
    return EchoResponseSignature(
        target_event_ids=target_array,
        delay_ranks=delay_array,
        reachable_mask=reachable,
        order_sign_matrix=response_order_sign_matrix(delay_array, reachable),
        label=label,
    )


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact sanity checks."""

    signatures = [
        _signature([1, 2, 3, 4], [2, 2, 5, 5], "p0"),
        _signature([2, 3, 4, 5], [1, 3, 3, 9], "p1"),
    ]
    profile = response_profile_from_signatures(signatures)
    classes = response_profile_equivalence_classes(profile)
    comparison = compare_response_profiles(profile, profile)
    ladder = default_response_representability_ladder()
    levels = {entry.level for entry in ladder}
    scalar_level = next(
        entry for entry in ladder if entry.level == "scalar_response_rank"
    )
    return [
        {
            "check": "profile_aligns_common_targets",
            "passed": float(
                np.array_equal(profile.target_event_ids, np.asarray([2, 3, 4]))
            ),
            "value": str(profile.target_event_ids),
        },
        {
            "check": "profile_separation_fraction_positive",
            "passed": float(response_profile_separation_fraction(profile) > 0.0),
            "value": str(response_profile_separation_fraction(profile)),
        },
        {
            "check": "equivalence_classes_exist",
            "passed": float(len(classes) >= 1),
            "value": str([group.tolist() for group in classes]),
        },
        {
            "check": "identical_profile_agreement",
            "passed": float(comparison["profile_entry_agreement_fraction"] == 1.0),
            "value": str(comparison),
        },
        {
            "check": "ladder_has_expected_levels",
            "passed": float(
                {
                    "scalar_response_rank",
                    "multi_reference_response_profile",
                    "calibrated_metric_representation",
                }
                <= levels
            ),
            "value": str(sorted(levels)),
        },
        {
            "check": "scalar_level_denies_target_distance",
            "passed": float(
                "target-target distance order"
                in scalar_level.what_it_does_not_imply
            ),
            "value": scalar_level.what_it_does_not_imply,
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
    print(f"Wrote response profile requirement exact sanity: {output_path}")


if __name__ == "__main__":
    main()
