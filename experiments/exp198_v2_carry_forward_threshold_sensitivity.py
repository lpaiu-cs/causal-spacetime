"""V2 carry-forward threshold sensitivity without retuning."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import save_metric_bar

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v2_carry_forward import (
    decide_v2_family_robustness,
    v2_metrics_rows_from_bundle,
)
from causal_spacetime_lab.state_change_manifest_v2_outputs import (
    load_v2_output_bundle,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for threshold sensitivity."""

    output_dir: Path = DEFAULT_OUTPUT_DIR
    heldout_thresholds: tuple[float, ...] = (0.15, 0.20, 0.25)
    null_gap_thresholds: tuple[float, ...] = (0.05, 0.10, 0.15)
    stricter_pass_thresholds: tuple[float, ...] = (0.25, 0.50, 0.75)


def _parse_float_cli_values(values: list[str]) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 threshold sensitivity.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--heldout-thresholds",
        nargs="+",
        default=["0.15", "0.20", "0.25"],
    )
    parser.add_argument(
        "--null-gap-thresholds",
        nargs="+",
        default=["0.05", "0.10", "0.15"],
    )
    parser.add_argument(
        "--stricter-pass-thresholds",
        nargs="+",
        default=["0.25", "0.50", "0.75"],
    )
    args = parser.parse_args()
    return ExperimentConfig(
        output_dir=args.output_dir,
        heldout_thresholds=_parse_float_cli_values(args.heldout_thresholds),
        null_gap_thresholds=_parse_float_cli_values(args.null_gap_thresholds),
        stricter_pass_thresholds=_parse_float_cli_values(args.stricter_pass_thresholds),
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run sensitivity rows without changing the fixed default criteria."""

    metrics = v2_metrics_rows_from_bundle(load_v2_output_bundle(config.output_dir))
    rows: list[dict[str, float | str]] = []
    default = default_cross_family_robustness_criteria()
    for heldout in config.heldout_thresholds:
        for null_gap in config.null_gap_thresholds:
            for stricter in config.stricter_pass_thresholds:
                criteria = replace(
                    default,
                    max_mean_heldout_violation=float(heldout),
                    min_destructive_null_gap=float(null_gap),
                    min_stricter_threshold_pass_fraction=float(stricter),
                )
                decisions = decide_v2_family_robustness(metrics, criteria)
                counts = {
                    "carry_forward": 0,
                    "provisional": 0,
                    "blocked": 0,
                    "report_only": 0,
                    "failed_control": 0,
                }
                for decision in decisions:
                    counts[decision.decision] += 1
                rows.append(
                    {
                        "heldout_threshold": float(heldout),
                        "null_gap_threshold": float(null_gap),
                        "stricter_pass_threshold": float(stricter),
                        **{
                            key + "_count": float(value)
                            for key, value in counts.items()
                        },
                        "analysis_label": "sensitivity_only_not_threshold_retuning",
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write v2 threshold sensitivity rows."""

    return write_csv(
        rows,
        output_dir / "data" / "v2_carry_forward_threshold_sensitivity.csv",
        [
            "heldout_threshold",
            "null_gap_threshold",
            "stricter_pass_threshold",
            "carry_forward_count",
            "provisional_count",
            "blocked_count",
            "report_only_count",
            "failed_control_count",
            "analysis_label",
        ],
    )


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save v2 threshold sensitivity figure."""

    return [
        save_metric_bar(
            rows,
            label_key="heldout_threshold",
            value_key="carry_forward_count",
            path=output_dir
            / "figures"
            / "v2_carry_forward_threshold_sensitivity.png",
            ylabel="carry-forward count",
        )
    ]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote v2 threshold sensitivity: {path}")
    for figure in figure_paths:
        print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
