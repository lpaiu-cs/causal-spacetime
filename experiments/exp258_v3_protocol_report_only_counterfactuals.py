"""Run report-only v3 protocol counterfactual diagnostics."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v3_protocol_blocked_v4_experiment_helpers import (
    data_path,
    read_csv_rows,
    to_float_rows,
)

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_blocking_analysis import (
    decompose_v3_protocol_blocking_by_family,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_counterfactuals import (
    family_would_remain_blocked_after_single_fix,
    heldout_threshold_counterfactual_report,
    latent_order_threshold_counterfactual_report,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")
    heldout_thresholds: tuple[float, ...] = (0.20, 0.25, 0.30, 0.35)
    latent_order_thresholds: tuple[float, ...] = (0.30, 0.40, 0.50)
    ignored_root_causes: tuple[str, ...] = (
        "heldout_failure",
        "latent_order_instability",
        "null_separation_failure",
    )


def _parse_values(values: list[str]) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 report-only counterfactuals.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument(
        "--heldout-thresholds",
        nargs="+",
        default=["0.20", "0.25", "0.30", "0.35"],
    )
    parser.add_argument(
        "--latent-order-thresholds",
        nargs="+",
        default=["0.30", "0.40", "0.50"],
    )
    parser.add_argument(
        "--ignored-root-causes",
        nargs="+",
        default=[
            "heldout_failure",
            "latent_order_instability",
            "null_separation_failure",
        ],
    )
    args = parser.parse_args()
    return ExperimentConfig(
        output_dir=args.output_dir,
        heldout_thresholds=_parse_values(args.heldout_thresholds),
        latent_order_thresholds=_parse_values(args.latent_order_thresholds),
        ignored_root_causes=tuple(args.ignored_root_causes),
    )


def run_experiment(
    config: ExperimentConfig,
) -> tuple[
    list[dict[str, float | str]],
    list[dict[str, float | str]],
    list[dict[str, float | str]],
]:
    metric_rows = to_float_rows(
        read_csv_rows(
            data_path(
                config.output_dir,
                "v3_protocol_cross_family_robustness_decision_metrics.csv",
            )
        )
    )
    decision_rows = to_float_rows(
        read_csv_rows(
            data_path(
                config.output_dir,
                "v3_protocol_cross_family_robustness_decisions.csv",
            )
        )
    )
    precondition_rows = to_float_rows(
        read_csv_rows(
            data_path(config.output_dir, "v3_protocol_precondition_audit.csv")
        )
    )
    records = decompose_v3_protocol_blocking_by_family(
        metric_rows,
        decision_rows,
        precondition_rows,
        default_cross_family_robustness_criteria(),
    )
    single_fix_rows: list[dict[str, float | str]] = []
    for root_cause in config.ignored_root_causes:
        single_fix_rows.extend(
            family_would_remain_blocked_after_single_fix(
                records,
                ignored_root_cause=root_cause,
            )
        )
    return (
        heldout_threshold_counterfactual_report(
            metric_rows,
            hypothetical_thresholds=list(config.heldout_thresholds),
        ),
        latent_order_threshold_counterfactual_report(
            metric_rows,
            hypothetical_thresholds=list(config.latent_order_thresholds),
        ),
        single_fix_rows,
    )


def main() -> None:
    config = parse_args()
    heldout, latent, single_fix = run_experiment(config)
    paths = [
        write_csv(
            heldout,
            data_path(config.output_dir, "v3_protocol_heldout_counterfactual.csv"),
            ["family_name", "hypothetical_threshold", "output_note"],
        ),
        write_csv(
            latent,
            data_path(config.output_dir, "v3_protocol_latent_order_counterfactual.csv"),
            ["family_name", "hypothetical_threshold", "output_note"],
        ),
        write_csv(
            single_fix,
            data_path(config.output_dir, "v3_protocol_single_fix_counterfactual.csv"),
            ["family_name", "ignored_root_cause", "output_note"],
        ),
    ]
    print(
        "Wrote v3 protocol report-only counterfactuals: "
        + ", ".join(map(str, paths))
    )


if __name__ == "__main__":
    main()
