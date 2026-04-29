"""Null baselines for frozen-manifest latent representation fits."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from manifest_representation_experiment_helpers import ensure_manifests, write_csv

from causal_spacetime_lab.state_change_manifest_representation import (
    ManifestRepresentationConfig,
)
from causal_spacetime_lab.state_change_manifest_representation_nulls import (
    evaluate_manifest_representation_nulls,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for frozen-manifest representation nulls."""

    manifest_dir: Path = Path("outputs/manifests")
    embedding_dim: int = 2
    null_repetitions: int = 5
    steps: int = 600
    restarts: int = 2
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Frozen-manifest null baselines.")
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


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Evaluate null representation baselines for eligible manifests."""

    loaded = ensure_manifests(config.manifest_dir, output_dir=config.output_dir)
    datasets = [dataset for dataset in loaded if dataset.eligible]
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
        rows.extend(asdict(summary) for summary in summaries)
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write null baseline CSV."""

    return write_csv(
        rows,
        output_dir / "data" / "frozen_manifest_representation_nulls.csv",
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
        ],
    )


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save null held-out violation figure."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    labels = sorted({str(row["null_type"]) for row in rows})
    values = []
    for label in labels:
        label_values = np.asarray(
            [
                float(row["mean_heldout_violation_rate"])
                for row in rows
                if row["null_type"] == label
            ],
            dtype=float,
        )
        values.append(float(np.nanmean(label_values)) if label_values.size else 0.0)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(labels, values)
    ax.set_ylabel("mean null held-out violation")
    ax.tick_params(axis="x", labelrotation=30)
    ax.grid(True, axis="y", alpha=0.3)
    path = figure_dir / "frozen_manifest_null_heldout_violation.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return [path]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote frozen manifest representation nulls: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
