"""Fit latent ordinal representations to frozen handoff manifests."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from manifest_representation_experiment_helpers import ensure_manifests, write_csv

from causal_spacetime_lab.state_change_manifest_representation import (
    ManifestRepresentationConfig,
    fit_manifest_dimension_curve,
    representation_fit_to_row,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for frozen-manifest ordinal representation diagnostics."""

    manifest_dir: Path = Path("outputs/manifests")
    dims: tuple[int, ...] = (1, 2, 3)
    steps: int = 800
    restarts: int = 3
    learning_rate: float = 0.05
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Frozen-manifest latent ordinal representation diagnostics."
    )
    parser.add_argument("--manifest-dir", type=Path, default=Path("outputs/manifests"))
    parser.add_argument("--dims", nargs="+", default=["1", "2", "3"])
    parser.add_argument("--steps", type=int, default=800)
    parser.add_argument("--restarts", type=int, default=3)
    parser.add_argument("--learning-rate", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        manifest_dir=args.manifest_dir,
        dims=_parse_int_values(args.dims),
        steps=args.steps,
        restarts=args.restarts,
        learning_rate=args.learning_rate,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Fit latent ordinal representation curves for frozen manifests."""

    datasets = ensure_manifests(config.manifest_dir, output_dir=config.output_dir)
    rows: list[dict[str, float | str]] = []
    for dataset_index, dataset in enumerate(datasets):
        base_config = ManifestRepresentationConfig(
            embedding_dim=int(config.dims[0]),
            steps=config.steps,
            restarts=config.restarts,
            learning_rate=config.learning_rate,
            seed=config.seed + 100_000 * dataset_index,
        )
        fits = fit_manifest_dimension_curve(
            dataset,
            list(config.dims),
            base_config,
            fit_ineligible=False,
        )
        for fit in fits:
            row = representation_fit_to_row(fit)
            row["failed_reasons"] = ";".join(dataset.failed_reasons)
            rows.append(row)
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write fit diagnostics CSV."""

    return write_csv(
        rows,
        output_dir / "data" / "frozen_manifest_ordinal_representation.csv",
        [
            "manifest_id",
            "embedding_dim",
            "eligible",
            "fitted",
            "reason_not_fit",
            "train_violation_rate",
            "heldout_violation_rate",
            "train_hinge_loss",
            "heldout_hinge_loss",
            "train_constraint_count",
            "heldout_constraint_count",
            "target_count",
            "representation_kind",
            "failed_reasons",
        ],
    )


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save fit summary figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    fitted = [row for row in rows if float(row["fitted"]) == 1.0]

    dims = sorted({float(row["embedding_dim"]) for row in fitted})
    means = []
    for dim in dims:
        values = np.asarray(
            [
                float(row["heldout_violation_rate"])
                for row in fitted
                if float(row["embedding_dim"]) == dim
            ],
            dtype=float,
        )
        means.append(float(np.nanmean(values)) if values.size else float("nan"))
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(dims, means, marker="o")
    ax.set_xlabel("latent ordinal representation dimension")
    ax.set_ylabel("mean held-out violation")
    ax.grid(True, alpha=0.3)
    path = figure_dir / "frozen_manifest_heldout_violation_by_dim.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(5, 5))
    train = [float(row["train_violation_rate"]) for row in fitted]
    heldout = [float(row["heldout_violation_rate"]) for row in fitted]
    ax.scatter(train, heldout, alpha=0.75)
    ax.set_xlabel("train violation")
    ax.set_ylabel("held-out violation")
    ax.grid(True, alpha=0.3)
    path = figure_dir / "frozen_manifest_train_vs_heldout.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)
    return paths


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote frozen manifest ordinal representation: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
