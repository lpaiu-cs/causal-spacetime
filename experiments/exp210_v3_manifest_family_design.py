"""Export planned-only v3 manifest-family design table."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_v3_design import (
    default_v3_manifest_family_designs,
    v3_manifest_family_design_table,
)


def run_experiment() -> list[dict[str, float | str]]:
    return v3_manifest_family_design_table(default_v3_manifest_family_designs())


def main() -> None:
    rows = run_experiment()
    path = write_csv(
        rows,
        Path("outputs/data/v3_manifest_family_design.csv"),
        [
            "family_name",
            "family_kind",
            "planned_manifest_count",
            "execution_status",
        ],
    )
    print(f"Wrote v3 manifest family design: {path}")


if __name__ == "__main__":
    main()

