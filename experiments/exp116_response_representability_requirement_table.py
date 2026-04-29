"""Write the response representability requirement ladder."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from causal_spacetime_lab.state_change_representability_requirements import (
    representability_ladder_table,
)

DEFAULT_OUTPUT = Path("outputs/data/response_representability_requirement_table.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Write representability ladder table.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / DEFAULT_OUTPUT.name


def run_experiment() -> list[dict[str, str]]:
    """Return ladder rows."""

    return representability_ladder_table()


def write_outputs(rows: list[dict[str, str]], path: Path = DEFAULT_OUTPUT) -> Path:
    """Write ladder CSV."""

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
    print(f"Wrote response representability requirement table: {output_path}")


if __name__ == "__main__":
    main()
