"""Exact sanity checks for response-order underdetermination."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_response_underdetermination import (
    pairwise_distance_order_disagreement,
    response_order_preserved,
)

DEFAULT_OUTPUT = Path("outputs/data/response_order_underdetermination_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Response-order underdetermination exact sanity."
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / DEFAULT_OUTPUT.name


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact underdetermination checks."""

    delays = np.asarray([1, 2, 3, 4], dtype=int)
    reachable = np.ones(delays.size, dtype=bool)
    layouts = [
        np.asarray([1.0, 2.0, 3.0, 4.0]),
        np.asarray([1.0, 2.0, 100.0, 101.0]),
        np.asarray([1.0, 10.0, 11.0, 100.0]),
    ]
    preserved = [
        response_order_preserved(delays, reachable, layout) for layout in layouts
    ]
    disagreements = [
        pairwise_distance_order_disagreement(layouts[0], layouts[1], reachable),
        pairwise_distance_order_disagreement(layouts[0], layouts[2], reachable),
        pairwise_distance_order_disagreement(layouts[1], layouts[2], reachable),
    ]
    return [
        {
            "check": "all_layouts_preserve_response_order",
            "passed": float(all(preserved)),
            "value": str(preserved),
        },
        {
            "check": "pair_distance_order_disagreement_exists",
            "passed": float(max(disagreements) > 0.0),
            "value": str(disagreements),
        },
        {
            "check": "same_scalar_order_different_pair_orders",
            "passed": float(all(preserved) and max(disagreements) > 0.0),
            "value": str(max(disagreements)),
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
    print(f"Wrote response-order underdetermination exact sanity: {output_path}")


if __name__ == "__main__":
    main()
