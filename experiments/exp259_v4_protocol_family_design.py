"""Export planned-only v4 protocol family design table."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_v4_design import (
    default_v4_protocol_family_designs,
    v4_protocol_family_design_table,
)


def run_experiment() -> list[dict[str, float | str]]:
    """Return planned-only v4 design rows."""

    return v4_protocol_family_design_table(default_v4_protocol_family_designs())


def main() -> None:
    path = write_csv(
        run_experiment(),
        Path("outputs/data/v4_protocol_family_design.csv"),
        ["family_name", "execution_status"],
    )
    print(f"Wrote v4 protocol family design: {path}")


if __name__ == "__main__":
    main()

