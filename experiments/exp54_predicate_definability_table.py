"""Predicate-definability table for transport and relational evolution."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from causal_spacetime_lab.predicate_definability import predicate_requirement_table

DEFAULT_OUTPUT_DIR = Path("outputs")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Predicate-definability table.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return args.output_dir


def write_outputs(rows: list[dict[str, str | bool]], output_dir: Path) -> Path:
    """Write predicate-definability table."""

    output_path = output_dir / "data" / "predicate_definability_table.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def run_experiment() -> list[dict[str, str | bool]]:
    """Return predicate-definability rows."""

    return predicate_requirement_table()


def main() -> None:
    output_dir = parse_args()
    rows = run_experiment()
    output_path = write_outputs(rows, output_dir)
    print(f"Wrote predicate definability table: {output_path}")


if __name__ == "__main__":
    main()
