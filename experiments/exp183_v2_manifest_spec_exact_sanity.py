"""Exact sanity checks for v2 manifest family specs."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import (
    DEFAULT_OUTPUT_DIR,
    load_v2_specs_with_fallback,
)

from causal_spacetime_lab.state_change_manifest_v2_spec import (
    mark_v2_family_specs_executed,
)


def run_experiment(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, float | str]]:
    """Run deterministic v2 spec checks."""

    specs = load_v2_specs_with_fallback(output_dir)
    executed = mark_v2_family_specs_executed(specs)
    names = {spec.family_name for spec in executed}
    expected = {
        "rank_gap_more_protocol_columns_v2",
        "rank_gap_coverage_enriched_v2",
        "rank_gap_rank_resolution_enriched_v2",
        "combined_diagnostic_complete_v2",
        "failed_controls_v2",
    }
    checks = [
        ("expected_family_names_present", expected <= names),
        (
            "execution_status_marked_executed",
            all(spec.execution_status == "executed" for spec in executed),
        ),
        (
            "planned_specs_loaded",
            all(
                spec.execution_status in {"planned_only", "executed"}
                for spec in specs
            ),
        ),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write exact sanity rows."""

    return write_csv(
        rows,
        output_dir / "data" / "v2_manifest_spec_exact_sanity.csv",
        ["check", "passed"],
    )


def main() -> None:
    path = write_outputs(run_experiment())
    print(f"Wrote v2 manifest spec sanity: {path}")


if __name__ == "__main__":
    main()
