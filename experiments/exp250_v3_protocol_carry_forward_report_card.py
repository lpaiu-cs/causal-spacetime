"""Aggregate v3 protocol carry-forward evaluation into one report card."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v3_protocol_carry_forward_experiment_helpers import data_path, read_csv_rows


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 protocol report card.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    decisions = read_csv_rows(
        data_path(
            config.output_dir, "v3_protocol_cross_family_robustness_decisions.csv"
        )
    )
    registry = {
        row["family_name"]: row
        for row in read_csv_rows(
            data_path(
                config.output_dir, "v3_protocol_carry_forward_registry_export.csv"
            )
        )
    }
    plan_rows = read_csv_rows(
        data_path(config.output_dir, "v3_protocol_stress_test_handoff_plan.csv")
    )
    stop_rows = [row for row in plan_rows if row.get("row_type") == "stop_condition"]
    stop_status = stop_rows[0].get("decision", "") if stop_rows else ""
    rows: list[dict[str, float | str]] = []
    for decision in decisions:
        reg = registry.get(decision["family_name"], {})
        key_metrics = (
            f"heldout={decision.get('mean_heldout_violation', '')};"
            f"null_gap={decision.get('destructive_null_gap', '')};"
            f"restart_std={decision.get('restart_std', '')};"
            f"latent_order={decision.get('latent_order_disagreement', '')}"
        )
        rows.append(
            {
                "family_name": decision["family_name"],
                "family_kind": decision["family_kind"],
                "decision": decision["decision"],
                "diagnostic_complete": decision.get("diagnostic_complete", ""),
                "preconditions_passed": decision.get("preconditions_passed", ""),
                "failed_preconditions": decision.get("failed_preconditions", ""),
                "key_metrics": key_metrics,
                "dominant_handoff_provenance_type": decision.get(
                    "dominant_handoff_provenance_type",
                    "",
                ),
                "allowed_future_use": reg.get("allowed_future_use", "report_only"),
                "stop_condition": stop_status,
                "interpretation_warning": (
                    "carry-forward is stress-test eligibility, not geometry"
                ),
            }
        )
    return rows


def main() -> None:
    config = parse_args()
    path = write_csv(
        run_experiment(config),
        data_path(config.output_dir, "v3_protocol_carry_forward_report_card.csv"),
        ["family_name", "decision", "interpretation_warning"],
    )
    print(f"Wrote v3 protocol carry-forward report card: {path}")


if __name__ == "__main__":
    main()
