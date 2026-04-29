"""Compare latent ordinal representation fits across manifest families."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import (
    assignment_lookup,
    assignments_for_datasets,
    load_or_create_family_datasets,
    save_family_metric_figure,
    write_csv,
)

from causal_spacetime_lab.state_change_manifest_family_diagnostics import (
    failed_manifest_accounting_summary,
    fit_diagnostic_row_to_dict,
    manifest_fit_diagnostic_row,
    summarize_family_fit_diagnostics,
)
from causal_spacetime_lab.state_change_manifest_representation import (
    ManifestRepresentationConfig,
    fit_manifest_dimension_curve,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for manifest-family fit comparison."""

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

    parser = argparse.ArgumentParser(description="Manifest-family fit comparison.")
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


def run_experiment(
    config: ExperimentConfig,
) -> tuple[
    list[dict[str, float | str]],
    list[dict[str, float | str]],
    list[dict[str, float | str]],
]:
    """Run family-level latent representation fit comparison."""

    datasets = load_or_create_family_datasets(config.manifest_dir, config.output_dir)
    assignments = assignments_for_datasets(datasets)
    assignment_by_id = assignment_lookup(assignments)
    fit_rows = []
    for dataset_index, dataset in enumerate(datasets):
        base_config = ManifestRepresentationConfig(
            embedding_dim=int(config.dims[0]),
            steps=config.steps,
            restarts=config.restarts,
            learning_rate=0.05,
            seed=config.seed + 100_000 * dataset_index,
        )
        fits = fit_manifest_dimension_curve(dataset, list(config.dims), base_config)
        assignment = assignment_by_id[dataset.manifest_id]
        fit_rows.extend(manifest_fit_diagnostic_row(fit, assignment) for fit in fits)
    row_dicts = [fit_diagnostic_row_to_dict(row) for row in fit_rows]
    summary = summarize_family_fit_diagnostics(fit_rows)
    failed = failed_manifest_accounting_summary(assignments)
    return row_dicts, summary, failed


def write_outputs(
    fit_rows: list[dict[str, float | str]],
    summary_rows: list[dict[str, float | str]],
    failed_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path, Path]:
    """Write family fit comparison CSVs."""

    data_dir = output_dir / "data"
    fit_path = write_csv(
        fit_rows,
        data_dir / "manifest_family_fit_comparison.csv",
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
    )
    summary_path = write_csv(
        summary_rows,
        data_dir / "manifest_family_fit_summary.csv",
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
    )
    failed_path = write_csv(
        failed_rows,
        data_dir / "manifest_family_failed_accounting.csv",
        ["row_type", "family_name", "family_kind", "reason", "count"],
    )
    return fit_path, summary_path, failed_path


def save_figures(
    summary_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> list[Path]:
    """Save family fit comparison figures."""

    return [
        save_family_metric_figure(
            summary_rows,
            metric="mean_heldout_violation",
            path=output_dir / "figures" / "manifest_family_heldout_violation.png",
            ylabel="mean held-out violation",
        ),
        save_family_metric_figure(
            summary_rows,
            metric="mean_generalization_gap",
            path=output_dir / "figures" / "manifest_family_generalization_gap.png",
            ylabel="mean generalization gap",
        ),
    ]


def main() -> None:
    config = parse_args()
    fit_rows, summary_rows, failed_rows = run_experiment(config)
    paths = write_outputs(fit_rows, summary_rows, failed_rows, config.output_dir)
    figure_paths = save_figures(summary_rows, config.output_dir)
    print(
        "Wrote manifest family fit comparison: "
        + ", ".join(str(path) for path in paths)
    )
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
