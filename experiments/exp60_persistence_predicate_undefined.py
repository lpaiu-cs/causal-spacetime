"""Persistence-dependent predicates are undefined without identity structure."""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
from pathlib import Path

from causal_spacetime_lab.persistence_definability import (
    object_identity_without_persistence,
    pair_distance_history_without_persistence,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Persistence predicate reports.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return args.output_dir


def run_experiment() -> list[dict[str, float | str | bool | None]]:
    """Return deterministic undefined persistence predicate reports."""

    return [
        asdict(object_identity_without_persistence()),
        asdict(pair_distance_history_without_persistence()),
    ]


def write_outputs(
    rows: list[dict[str, float | str | bool | None]],
    output_dir: Path,
) -> Path:
    """Write persistence predicate rows."""

    output_path = output_dir / "data" / "persistence_predicate_undefined.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def main() -> None:
    output_dir = parse_args()
    rows = run_experiment()
    output_path = write_outputs(rows, output_dir)
    print(f"Wrote persistence predicate data: {output_path}")


if __name__ == "__main__":
    main()
