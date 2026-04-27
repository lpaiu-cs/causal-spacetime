"""Cross-slice predicates are undefined without transport protocols."""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
from pathlib import Path

from causal_spacetime_lab.cross_slice import (
    same_direction_without_transport,
    same_position_without_transport,
    spatial_evolution_without_transport,
    velocity_without_transport,
)

DEFAULT_OUTPUT = Path("outputs/data/cross_slice_predicate_undefined.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Cross-slice predicates without transport.",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / "cross_slice_predicate_undefined.csv"


def run_experiment() -> list[dict[str, float | str | bool | None]]:
    """Return deterministic undefined predicate reports."""

    reports = [
        same_position_without_transport(),
        same_direction_without_transport(),
        velocity_without_transport(),
        spatial_evolution_without_transport(),
    ]
    return [asdict(report) for report in reports]


def write_outputs(
    rows: list[dict[str, float | str | bool | None]],
    output_path: Path = DEFAULT_OUTPUT,
) -> Path:
    """Write undefined predicate rows."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def main() -> None:
    output_path = parse_args()
    rows = run_experiment()
    output_path = write_outputs(rows, output_path)
    print(f"Wrote cross-slice undefined predicate data: {output_path}")


if __name__ == "__main__":
    main()
