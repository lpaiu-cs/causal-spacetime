"""Apply fixed stricter criteria to manifest-family diagnostics."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from exp149_manifest_family_fit_comparison import (
    ExperimentConfig as FitComparisonConfig,
)
from exp149_manifest_family_fit_comparison import run_experiment as run_fit_comparison
from manifest_family_experiment_helpers import save_bar_figure, write_csv

from causal_spacetime_lab.state_change_manifest_family_comparison import (
    compare_manifest_family_against_thresholds,
)
from causal_spacetime_lab.state_change_manifest_family_config import (
    FamilyComparisonConfig,
    default_family_comparison_config,
)
from causal_spacetime_lab.state_change_manifest_family_diagnostics import (
    ManifestFitDiagnosticRow,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for stricter family criteria."""

    manifest_dir: Path = Path("outputs/manifests")
    dims: tuple[int, ...] = (1, 2, 3)
    steps: int = 800
    restarts: int = 3
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Manifest-family stricter criteria.")
    parser.add_argument("--manifest-dir", type=Path, default=Path("outputs/manifests"))
    parser.add_argument("--dims", nargs="+", default=["1", "2", "3"])
    parser.add_argument("--steps", type=int, default=800)
    parser.add_argument("--restarts", type=int, default=3)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        manifest_dir=args.manifest_dir,
        dims=_parse_int_values(args.dims),
        steps=args.steps,
        restarts=args.restarts,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _row_to_dataclass(row: dict[str, float | str]) -> ManifestFitDiagnosticRow:
    return ManifestFitDiagnosticRow(
        manifest_id=str(row["manifest_id"]),
        family_name=str(row["family_name"]),
        family_kind=str(row["family_kind"]),
        eligible=bool(float(row["eligible"])),
        fitted=bool(float(row["fitted"])),
        reason_not_fit=str(row["reason_not_fit"]),
        embedding_dim=int(float(row["embedding_dim"])),
        train_violation_rate=float(row["train_violation_rate"]),
        heldout_violation_rate=float(row["heldout_violation_rate"]),
        generalization_gap=float(row["generalization_gap"]),
        train_hinge_loss=float(row["train_hinge_loss"]),
        heldout_hinge_loss=float(row["heldout_hinge_loss"]),
        target_count=int(float(row["target_count"])),
        train_constraint_count=int(float(row["train_constraint_count"])),
        heldout_constraint_count=int(float(row["heldout_constraint_count"])),
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Apply fixed stricter criteria to family fit diagnostics."""

    fit_rows, _summary_rows, _failed_rows = run_fit_comparison(
        FitComparisonConfig(
            manifest_dir=config.manifest_dir,
            dims=config.dims,
            steps=config.steps,
            restarts=config.restarts,
            seed=config.seed,
            output_dir=config.output_dir,
        )
    )
    row_objects = [_row_to_dataclass(row) for row in fit_rows]
    fixed = default_family_comparison_config()
    comparison_config = FamilyComparisonConfig(
        dims=list(config.dims),
        steps=config.steps,
        restarts=config.restarts,
        learning_rate=fixed.learning_rate,
        null_repetitions=fixed.null_repetitions,
        restart_count=fixed.restart_count,
        heldout_violation_threshold=fixed.heldout_violation_threshold,
        generalization_gap_threshold=fixed.generalization_gap_threshold,
        min_null_separation=fixed.min_null_separation,
        max_restart_std=fixed.max_restart_std,
        max_latent_order_disagreement=fixed.max_latent_order_disagreement,
        seed=config.seed,
    )
    rows = compare_manifest_family_against_thresholds(row_objects, comparison_config)
    for row in rows:
        row["null_separation_pass"] = float("nan")
        row["restart_std_pass"] = float("nan")
        row["latent_order_stability_pass"] = float("nan")
        row["criteria_note"] = "fixed_thresholds_no_retuning"
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write stricter criteria CSV."""

    return write_csv(
        rows,
        output_dir / "data" / "manifest_family_stricter_criteria.csv",
        [
            "family_name",
            "family_kind",
            "embedding_dim",
            "manifest_count",
            "fitted_count",
            "no_fit_count",
            "mean_train_violation",
            "mean_heldout_violation",
            "mean_generalization_gap",
            "median_heldout_violation",
            "best_heldout_violation",
            "worst_heldout_violation",
            "heldout_threshold",
            "generalization_gap_threshold",
            "heldout_pass",
            "generalization_gap_pass",
            "threshold_pass",
            "null_separation_pass",
            "restart_std_pass",
            "latent_order_stability_pass",
            "criteria_note",
        ],
    )


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save stricter criteria pass-rate figure."""

    labels = sorted({str(row["family_name"]) for row in rows})
    values = []
    for label in labels:
        family_values = [
            float(row["threshold_pass"]) for row in rows if row["family_name"] == label
        ]
        values.append(float(np.mean(family_values)) if family_values else 0.0)
    return [
        save_bar_figure(
            labels,
            values,
            output_dir / "figures" / "manifest_family_stricter_pass_rate.png",
            ylabel="stricter criteria pass rate",
        )
    ]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote manifest family stricter criteria: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
