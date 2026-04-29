"""Restart stability diagnostics for frozen-manifest representation fits."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_representation_experiment_helpers import (
    ensure_manifests,
    save_simple_bar,
    write_csv,
)

from causal_spacetime_lab.state_change_manifest_representation import (
    ManifestRepresentationConfig,
)
from causal_spacetime_lab.state_change_manifest_representation_stability import (
    fit_manifest_restarts,
    heldout_violation_stability_summary,
    pairwise_latent_order_stability,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for restart stability diagnostics."""

    manifest_dir: Path = Path("outputs/manifests")
    embedding_dim: int = 2
    restart_count: int = 10
    steps: int = 600
    learning_rate: float = 0.05
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Frozen-manifest fit stability.")
    parser.add_argument("--manifest-dir", type=Path, default=Path("outputs/manifests"))
    parser.add_argument("--embedding-dim", type=int, default=2)
    parser.add_argument("--restart-count", type=int, default=10)
    parser.add_argument("--steps", type=int, default=600)
    parser.add_argument("--learning-rate", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        manifest_dir=args.manifest_dir,
        embedding_dim=args.embedding_dim,
        restart_count=args.restart_count,
        steps=args.steps,
        learning_rate=args.learning_rate,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run restart stability diagnostics for eligible manifests."""

    rows: list[dict[str, float | str]] = []
    loaded = ensure_manifests(config.manifest_dir, output_dir=config.output_dir)
    datasets = [dataset for dataset in loaded if dataset.eligible]
    for dataset_index, dataset in enumerate(datasets):
        fit_config = ManifestRepresentationConfig(
            embedding_dim=config.embedding_dim,
            steps=config.steps,
            restarts=1,
            learning_rate=config.learning_rate,
            seed=config.seed + 100_000 * dataset_index,
        )
        fits = fit_manifest_restarts(
            dataset,
            fit_config,
            restart_count=config.restart_count,
        )
        heldout = heldout_violation_stability_summary(fits)
        latent = pairwise_latent_order_stability(
            fits,
            seed=config.seed + 10_000 * dataset_index,
        )
        rows.append(
            {
                "manifest_id": dataset.manifest_id,
                "embedding_dim": float(config.embedding_dim),
                "restart_count": float(config.restart_count),
                **heldout,
                **latent,
            }
        )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write stability CSV."""

    return write_csv(
        rows,
        output_dir / "data" / "frozen_manifest_fit_stability.csv",
        [
            "manifest_id",
            "embedding_dim",
            "restart_count",
            "fit_count",
            "mean_heldout_violation_rate",
            "std_heldout_violation_rate",
            "min_heldout_violation_rate",
            "max_heldout_violation_rate",
            "fit_pair_count",
            "mean_pair_order_disagreement",
            "max_pair_order_disagreement",
        ],
    )


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save restart stability figures."""

    labels = [str(row["manifest_id"])[:12] for row in rows]
    paths = []
    paths.append(
        save_simple_bar(
            labels,
            [float(row["std_heldout_violation_rate"]) for row in rows],
            output_dir / "figures" / "frozen_manifest_restart_stability.png",
            ylabel="held-out violation std",
        )
    )
    paths.append(
        save_simple_bar(
            labels,
            [float(row["mean_pair_order_disagreement"]) for row in rows],
            output_dir / "figures" / "frozen_manifest_latent_order_stability.png",
            ylabel="mean latent pair-order disagreement",
        )
    )
    return paths


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote frozen manifest fit stability: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
