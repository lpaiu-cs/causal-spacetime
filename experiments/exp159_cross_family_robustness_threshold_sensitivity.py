"""Sensitivity analysis for fixed cross-family robustness thresholds."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from pathlib import Path

from cross_family_robustness_experiment_helpers import load_metrics_and_decisions
from manifest_family_experiment_helpers import save_bar_figure, write_csv

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for threshold-sensitivity diagnostics."""

    output_dir: Path = DEFAULT_OUTPUT_DIR
    heldout_thresholds: tuple[float, ...] = (0.15, 0.20, 0.25)
    null_gap_thresholds: tuple[float, ...] = (0.05, 0.10, 0.15)
    stricter_pass_thresholds: tuple[float, ...] = (0.25, 0.50, 0.75)


def _float_tuple(values: list[str]) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Cross-family threshold sensitivity.")
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
        heldout_thresholds=_float_tuple(args.heldout_thresholds),
        null_gap_thresholds=_float_tuple(args.null_gap_thresholds),
        stricter_pass_thresholds=_float_tuple(args.stricter_pass_thresholds),
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Sweep predeclared thresholds without selecting a preferred setting."""

    rows: list[dict[str, float | str]] = []
    base = default_cross_family_robustness_criteria()
    for heldout in config.heldout_thresholds:
        for null_gap in config.null_gap_thresholds:
            for stricter in config.stricter_pass_thresholds:
                criteria = replace(
                    base,
                    max_mean_heldout_violation=heldout,
                    min_destructive_null_gap=null_gap,
                    min_stricter_threshold_pass_fraction=stricter,
                )
                _, decisions, missing = load_metrics_and_decisions(
                    config.output_dir,
                    criteria,
                )
                counts: dict[str, float] = {}
                for decision in decisions:
                    counts[decision.decision] = counts.get(decision.decision, 0.0) + 1.0
                rows.append(
                    {
                        "heldout_threshold": heldout,
                        "null_gap_threshold": null_gap,
                        "stricter_pass_threshold": stricter,
                        "carry_forward_count": counts.get("carry_forward", 0.0),
                        "provisional_count": counts.get("provisional", 0.0),
                        "blocked_count": counts.get("blocked", 0.0),
                        "report_only_count": counts.get("report_only", 0.0),
                        "failed_control_count": counts.get("failed_control", 0.0),
                        "missing_inputs": ";".join(missing),
                        "analysis_role": "sensitivity_only_no_threshold_selection",
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write threshold sensitivity CSV."""

    return write_csv(
        rows,
        output_dir / "data" / "cross_family_robustness_threshold_sensitivity.csv",
        [
            "heldout_threshold",
            "null_gap_threshold",
            "stricter_pass_threshold",
            "carry_forward_count",
            "provisional_count",
            "blocked_count",
            "report_only_count",
            "failed_control_count",
            "missing_inputs",
            "analysis_role",
        ],
    )


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save threshold sensitivity figure."""

    labels = [
        f"h={row['heldout_threshold']},n={row['null_gap_threshold']},s={row['stricter_pass_threshold']}"
        for row in rows
    ]
    values = [float(row["carry_forward_count"]) for row in rows]
    return [
        save_bar_figure(
            labels,
            values,
            output_dir / "figures" / "cross_family_threshold_sensitivity.png",
            ylabel="carry-forward count",
        )
    ]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_outputs(rows, config.output_dir)
    figures = save_figures(rows, config.output_dir)
    print(f"Wrote cross-family threshold sensitivity: {path}")
    for figure in figures:
        print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
