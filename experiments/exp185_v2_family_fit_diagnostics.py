"""Produce v2 family latent ordinal fit diagnostics."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import (
    DEFAULT_OUTPUT_DIR,
    data_path,
    figure_path,
    parse_int_values,
    save_metric_bar,
)

from causal_spacetime_lab.state_change_manifest_v2_diagnostics import (
    compute_v2_manifest_family_fit_rows,
    summarize_v2_fit_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 fit diagnostics."""

    manifest_dir: Path = Path("outputs/manifests_v2")
    dims: tuple[int, ...] = (1, 2, 3)
    steps: int = 600
    restarts: int = 2
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 family fit diagnostics.")
    parser.add_argument(
        "--manifest-dir",
        type=Path,
        default=Path("outputs/manifests_v2"),
    )
    parser.add_argument("--dims", nargs="+", default=["1", "2", "3"])
    parser.add_argument("--steps", type=int, default=600)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        manifest_dir=args.manifest_dir,
        dims=parse_int_values(args.dims),
        steps=args.steps,
        restarts=args.restarts,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    """Run fit diagnostics and family summary."""

    rows = compute_v2_manifest_family_fit_rows(
        config.manifest_dir,
        dims=list(config.dims),
        steps=config.steps,
        restarts=config.restarts,
        seed=config.seed,
    )
    return rows, summarize_v2_fit_rows(rows)


def write_outputs(
    rows: list[dict[str, float | str]],
    summary_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write v2 fit diagnostics."""

    return (
        write_csv(
            rows,
            data_path(output_dir, "v2_manifest_family_fit_comparison.csv"),
            [
                "manifest_id",
                "family_name",
                "family_kind",
                "eligible",
                "fitted",
                "reason_not_fit",
                "embedding_dim",
                "train_violation_rate",
                "heldout_violation_rate",
                "generalization_gap",
                "train_hinge_loss",
                "heldout_hinge_loss",
                "target_count",
                "train_constraint_count",
                "heldout_constraint_count",
            ],
        ),
        write_csv(
            summary_rows,
            data_path(output_dir, "v2_manifest_family_fit_summary.csv"),
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
            ],
        ),
    )


def save_figures(
    summary_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> list[Path]:
    """Save v2 fit summary figures."""

    return [
        save_metric_bar(
            summary_rows,
            label_key="family_name",
            value_key="mean_heldout_violation",
            path=figure_path(output_dir, "v2_manifest_family_heldout_violation.png"),
            ylabel="mean held-out violation",
        ),
        save_metric_bar(
            summary_rows,
            label_key="family_name",
            value_key="mean_generalization_gap",
            path=figure_path(output_dir, "v2_manifest_family_generalization_gap.png"),
            ylabel="mean generalization gap",
        ),
    ]


def main() -> None:
    config = parse_args()
    rows, summary_rows = run_experiment(config)
    paths = write_outputs(rows, summary_rows, config.output_dir)
    figure_paths = save_figures(summary_rows, config.output_dir)
    print("Wrote v2 fit diagnostics: " + ", ".join(str(path) for path in paths))
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
