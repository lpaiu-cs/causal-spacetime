"""Export planned v2 manifest-family specifications."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_remediation_plan import (
    default_new_manifest_family_specs_v2,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def run_experiment() -> list[dict[str, float | str]]:
    """Return planned v2 manifest-family specification rows."""

    return default_new_manifest_family_specs_v2()


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write planned v2 manifest-family design CSV."""

    return write_csv(
        rows,
        output_dir / "data" / "new_manifest_family_design_v2.csv",
        [
            "family_name",
            "family_kind",
            "profile_column_policy",
            "target_inclusion_policy",
            "pair_node_coverage_policy",
            "null_taxonomy_policy",
            "restart_stability_required",
            "latent_order_stability_required",
            "execution_status",
        ],
    )


def main() -> None:
    path = write_outputs(run_experiment())
    print(f"Wrote new manifest family design v2: {path}")


if __name__ == "__main__":
    main()
