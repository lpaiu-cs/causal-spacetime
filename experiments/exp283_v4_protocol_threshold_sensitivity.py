"""V4 protocol threshold sensitivity without retuning."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import save_metric_bar
from v4_protocol_carry_forward_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_v4_carry_forward import (
    v4_protocol_metrics_rows_from_bundle,
)
from causal_spacetime_lab.state_change_manifest_v4_outputs import (
    load_v4_protocol_output_bundle,
)
from causal_spacetime_lab.state_change_manifest_v4_preconditions import (
    evaluate_v4_protocol_preconditions,
)
from causal_spacetime_lab.state_change_manifest_v4_threshold_sensitivity import (
    v4_protocol_threshold_sensitivity_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")
    heldout_thresholds: tuple[float, ...] = (0.15, 0.20, 0.25)
    null_gap_thresholds: tuple[float, ...] = (0.05, 0.10, 0.15)
    stricter_pass_thresholds: tuple[float, ...] = (0.25, 0.50, 0.75)


def _parse_float_values(values: list[str]) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V4 protocol threshold sensitivity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument(
        "--heldout-thresholds", nargs="+", default=["0.15", "0.20", "0.25"]
    )
    parser.add_argument(
        "--null-gap-thresholds", nargs="+", default=["0.05", "0.10", "0.15"]
    )
    parser.add_argument(
        "--stricter-pass-thresholds",
        nargs="+",
        default=["0.25", "0.50", "0.75"],
    )
    args = parser.parse_args()
    return ExperimentConfig(
        output_dir=args.output_dir,
        heldout_thresholds=_parse_float_values(args.heldout_thresholds),
        null_gap_thresholds=_parse_float_values(args.null_gap_thresholds),
        stricter_pass_thresholds=_parse_float_values(args.stricter_pass_thresholds),
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    bundle = load_v4_protocol_output_bundle(config.output_dir)
    metrics = v4_protocol_metrics_rows_from_bundle(bundle)
    preconditions = evaluate_v4_protocol_preconditions(
        config.output_dir / "manifests_v4",
        bundle,
    )
    return v4_protocol_threshold_sensitivity_rows(
        metrics,
        preconditions,
        heldout_thresholds=list(config.heldout_thresholds),
        null_gap_thresholds=list(config.null_gap_thresholds),
        stricter_pass_thresholds=list(config.stricter_pass_thresholds),
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_csv(
        rows,
        data_path(config.output_dir, "v4_protocol_threshold_sensitivity.csv"),
        [
            "heldout_threshold",
            "null_gap_threshold",
            "stricter_pass_threshold",
            "carry_forward_count",
            "provisional_count",
            "blocked_count",
            "report_only_count",
            "failed_control_count",
            "note",
        ],
    )
    figure = save_metric_bar(
        rows,
        label_key="heldout_threshold",
        value_key="carry_forward_count",
        path=config.output_dir / "figures" / "v4_protocol_threshold_sensitivity.png",
        ylabel="carry-forward count",
    )
    print(f"Wrote v4 protocol threshold sensitivity: {path}")
    print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
