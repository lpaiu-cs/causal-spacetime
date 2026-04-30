"""Produce v2 null taxonomy diagnostics."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import (
    DEFAULT_OUTPUT_DIR,
    data_path,
    figure_path,
    save_metric_bar,
)

from causal_spacetime_lab.state_change_manifest_family_comparison import (
    summarize_nulls_by_taxonomy,
)
from causal_spacetime_lab.state_change_manifest_v2_diagnostics import (
    compute_v2_null_taxonomy_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 null taxonomy diagnostics."""

    manifest_dir: Path = Path("outputs/manifests_v2")
    embedding_dim: int = 2
    null_repetitions: int = 5
    steps: int = 600
    restarts: int = 2
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 null taxonomy diagnostics.")
    parser.add_argument(
        "--manifest-dir",
        type=Path,
        default=Path("outputs/manifests_v2"),
    )
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


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run v2 representation null taxonomy diagnostics."""

    return compute_v2_null_taxonomy_rows(
        config.manifest_dir,
        embedding_dim=config.embedding_dim,
        null_repetitions=config.null_repetitions,
        steps=config.steps,
        restarts=config.restarts,
        seed=config.seed,
    )


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write v2 null taxonomy diagnostics."""

    return write_csv(
        rows,
        data_path(output_dir, "v2_manifest_family_null_taxonomy.csv"),
        [
            "manifest_id",
            "family_name",
            "family_kind",
            "null_type",
            "taxonomy_class",
            "embedding_dim",
            "repetitions",
            "mean_heldout_violation_rate",
            "std_heldout_violation_rate",
            "best_heldout_violation_rate",
            "structured_heldout_violation_rate",
            "structured_minus_null_mean",
            "taxonomy_note",
        ],
    )


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> list[Path]:
    """Save v2 null taxonomy figures."""

    taxonomy_rows = summarize_nulls_by_taxonomy(rows)
    return [
        save_metric_bar(
            taxonomy_rows,
            label_key="taxonomy_class",
            value_key="mean_null_heldout_violation",
            path=figure_path(
                output_dir,
                "v2_manifest_null_taxonomy_heldout_violation.png",
            ),
            ylabel="mean null held-out violation",
        )
    ]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    output_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote v2 null taxonomy diagnostics: {output_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
