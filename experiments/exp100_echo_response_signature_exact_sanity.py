"""Exact sanity checks for echo-response order signatures."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_response_signature import (
    echo_response_signature_from_motif_rows,
    response_order_sign_matrix,
    signature_reachable_fraction,
    signature_strict_pair_fraction,
    signature_tie_fraction,
)
from causal_spacetime_lab.state_change_response_signature_comparison import (
    compare_response_signatures,
    response_order_cycle_count,
    stable_response_order_core,
)

DEFAULT_OUTPUT = Path("outputs/data/echo_response_signature_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Echo response signature sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / "echo_response_signature_exact_sanity.csv"


def _row(name: str, passed: bool, value: float | str) -> dict[str, float | str]:
    return {"check": name, "passed": float(passed), "value": value}


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact response-signature checks."""

    rows = [
        {"target_event_id": 10.0, "recovered_delay_rank": 2.0},
        {"target_event_id": 11.0, "recovered_delay_rank": 5.0},
        {"target_event_id": 12.0, "recovered_delay_rank": 5.0},
    ]
    signature = echo_response_signature_from_motif_rows(rows, label="baseline")
    expected_signs = np.asarray(
        [
            [0, -1, -1],
            [1, 0, 0],
            [1, 0, 0],
        ],
        dtype=int,
    )
    variant = echo_response_signature_from_motif_rows(
        [
            {"target_event_id": 10.0, "recovered_delay_rank": 2.0},
            {"target_event_id": 11.0, "recovered_delay_rank": 4.0},
            {"target_event_id": 12.0, "recovered_delay_rank": 6.0},
        ],
        label="variant",
    )
    comparison = compare_response_signatures(signature, variant)
    stable = stable_response_order_core([signature, variant])
    cyclic = np.asarray(
        [
            [0, -1, 1],
            [1, 0, -1],
            [-1, 1, 0],
        ],
        dtype=int,
    )
    return [
        _row("reachable_fraction", signature_reachable_fraction(signature) == 1.0, 1.0),
        _row(
            "tie_fraction",
            signature_tie_fraction(signature) == 1.0 / 3.0,
            signature_tie_fraction(signature),
        ),
        _row(
            "strict_pair_fraction",
            signature_strict_pair_fraction(signature) == 2.0 / 3.0,
            signature_strict_pair_fraction(signature),
        ),
        _row(
            "order_sign_matrix",
            np.array_equal(signature.order_sign_matrix, expected_signs),
            str(signature.order_sign_matrix),
        ),
        _row(
            "direct_sign_matrix",
            np.array_equal(
                response_order_sign_matrix([2, 5, 5], [True, True, True]),
                expected_signs,
            ),
            1.0,
        ),
        _row(
            "comparison_pair_count",
            comparison["pair_count"] == 3.0,
            comparison["pair_count"],
        ),
        _row(
            "stable_core_nonempty",
            float(stable["stable_pair_fraction"]) > 0.0,
            float(stable["stable_pair_fraction"]),
        ),
        _row("cycle_count", response_order_cycle_count(cyclic) == 1, 1.0),
    ]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_path: Path = DEFAULT_OUTPUT,
) -> Path:
    """Write exact sanity rows."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def main() -> None:
    output_path = write_outputs(run_experiment(), parse_args())
    print(f"Wrote echo response signature exact sanity: {output_path}")


if __name__ == "__main__":
    main()

