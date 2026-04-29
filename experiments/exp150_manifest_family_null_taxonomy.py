"""Classify frozen-manifest representation nulls by taxonomy."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path

from manifest_family_experiment_helpers import (
    load_or_create_family_datasets,
    save_bar_figure,
    write_csv,
)

from causal_spacetime_lab.state_change_manifest_family_comparison import (
    summarize_nulls_by_taxonomy,
)
from causal_spacetime_lab.state_change_manifest_null_taxonomy import (
    classify_null_type,
)
from causal_spacetime_lab.state_change_manifest_representation import (
    ManifestRepresentationConfig,
)
from causal_spacetime_lab.state_change_manifest_representation_nulls import (
    evaluate_manifest_representation_nulls,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for null taxonomy comparison."""

    manifest_dir: Path = Path("outputs/manifests")
    embedding_dim: int = 2
    null_repetitions: int = 5
    steps: int = 600
    restarts: int = 2
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Manifest-family null taxonomy.")
    parser.add_argument("--manifest-dir", type=Path, default=Path("outputs/manifests"))
    parser.add_argument("--embedding-dim", type=int, default=2)
    parser.add_argument("--null-repetitions", type=int, default=5)
    parser.add_argument("--steps", type=int, default=600)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        manifest_dir=args.manifest_dir,
        embedding_dim=args.embedding_dim,
        null_repetitions=args.null_repetitions,
        steps=args.steps,
        restarts=args.restarts,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    """Evaluate and classify representation nulls."""

    datasets = [
        dataset
        for dataset in load_or_create_family_datasets(
            config.manifest_dir,
            config.output_dir,
        )
        if dataset.eligible
    ]
    rows: list[dict[str, float | str]] = []
    for dataset_index, dataset in enumerate(datasets):
        fit_config = ManifestRepresentationConfig(
            embedding_dim=config.embedding_dim,
            steps=config.steps,
            restarts=config.restarts,
            seed=config.seed + 100_000 * dataset_index,
        )
        summaries = evaluate_manifest_representation_nulls(
            dataset,
            fit_config,
            null_repetitions=config.null_repetitions,
            seed=config.seed + 10_000 * dataset_index,
        )
        for summary in summaries:
            row = asdict(summary)
            row["taxonomy_class"] = classify_null_type(summary.null_type)
            row["taxonomy_note"] = (
                "target_label_permutation_is_symmetry_control"
                if summary.null_type == "permuted_targets"
                else ""
            )
            rows.append(row)
    return rows, summarize_nulls_by_taxonomy(rows)


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write null taxonomy CSV."""

    return write_csv(
        rows,
        output_dir / "data" / "manifest_family_null_taxonomy.csv",
        [
            "manifest_id",
            "null_type",
            "embedding_dim",
            "repetitions",
            "mean_heldout_violation_rate",
            "std_heldout_violation_rate",
            "best_heldout_violation_rate",
            "structured_heldout_violation_rate",
            "structured_minus_null_mean",
            "taxonomy_class",
            "taxonomy_note",
        ],
    )


def save_figures(
    taxonomy_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> list[Path]:
    """Save null taxonomy figure."""

    labels = [str(row["taxonomy_class"]) for row in taxonomy_rows]
    values = [float(row["mean_null_heldout_violation"]) for row in taxonomy_rows]
    return [
        save_bar_figure(
            labels,
            values,
            output_dir / "figures" / "manifest_null_taxonomy_heldout_violation.png",
            ylabel="mean null held-out violation",
        )
    ]


def main() -> None:
    config = parse_args()
    rows, taxonomy_rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(taxonomy_rows, config.output_dir)
    print(f"Wrote manifest family null taxonomy: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
