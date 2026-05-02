"""Export planned-only v5 family design table."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v4_blocked_v5_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_v5_design import (
    default_v5_protocol_family_designs,
    v5_protocol_family_design_table,
)


def run_experiment() -> list[dict[str, float | str]]:
    """Return planned-only v5 design rows."""

    return v5_protocol_family_design_table(default_v5_protocol_family_designs())


def main() -> None:
    output_dir = Path("outputs")
    path = write_csv(
        run_experiment(),
        data_path(output_dir, "v5_protocol_family_design.csv"),
        ["family_name", "execution_status"],
    )
    print(f"Wrote v5 protocol family design: {path}")


if __name__ == "__main__":
    main()
