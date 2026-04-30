"""Aggregate v2 carry-forward evaluation into one report card."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path, read_csv_rows

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 report card."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 carry-forward report card.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Build v2 carry-forward report card rows."""

    decisions = read_csv_rows(
        data_path(config.output_dir, "v2_cross_family_robustness_decisions.csv")
    )
    registry_rows = {
        row["family_name"]: row
        for row in read_csv_rows(
            data_path(config.output_dir, "v2_carry_forward_registry_export.csv")
        )
    }
    plan_rows = read_csv_rows(
        data_path(config.output_dir, "v2_stress_test_handoff_plan.csv")
    )
    stop_rows = [row for row in plan_rows if row.get("row_type") == "stop_condition"]
    stop_status = stop_rows[0].get("decision", "") if stop_rows else ""
    rows: list[dict[str, float | str]] = []
    for decision in decisions:
        registry = registry_rows.get(decision["family_name"], {})
        rows.append(
            {
                "family_name": decision["family_name"],
                "family_kind": decision["family_kind"],
                "decision": decision["decision"],
                "diagnostic_complete": decision.get("diagnostic_complete", ""),
                "mean_heldout_violation": decision.get("mean_heldout_violation", ""),
                "destructive_null_gap": decision.get("destructive_null_gap", ""),
                "restart_std": decision.get("restart_std", ""),
                "latent_order_disagreement": decision.get(
                    "latent_order_disagreement",
                    "",
                ),
                "failed_reasons": decision.get("failed_reasons", ""),
                "allowed_future_use": registry.get("allowed_future_use", "report_only"),
                "stop_condition": stop_status,
                "interpretation_warning": (
                    "carry-forward is stress-test eligibility, not geometry"
                ),
            }
        )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write v2 carry-forward report card."""

    return write_csv(
        rows,
        data_path(output_dir, "v2_carry_forward_report_card.csv"),
        [
            "family_name",
            "family_kind",
            "decision",
            "diagnostic_complete",
            "mean_heldout_violation",
            "destructive_null_gap",
            "restart_std",
            "latent_order_disagreement",
            "failed_reasons",
            "allowed_future_use",
            "stop_condition",
            "interpretation_warning",
        ],
    )


def main() -> None:
    config = parse_args()
    path = write_outputs(run_experiment(config), config.output_dir)
    print(f"Wrote v2 carry-forward report card: {path}")


if __name__ == "__main__":
    main()
