"""Write fixed cross-family robustness criteria to CSV."""

from __future__ import annotations

import csv
from dataclasses import asdict
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def run_experiment() -> list[dict[str, float | str]]:
    """Return fixed robustness criteria as threshold rows."""

    criteria = default_cross_family_robustness_criteria()
    return [
        {"criterion": key, "value": value}
        for key, value in asdict(criteria).items()
    ]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write criteria table."""

    path = output_dir / "data" / "cross_family_robustness_criteria_table.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=["criterion", "value"])
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    path = write_outputs(run_experiment())
    print(f"Wrote cross-family robustness criteria table: {path}")


if __name__ == "__main__":
    main()
