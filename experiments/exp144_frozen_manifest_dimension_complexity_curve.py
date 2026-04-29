"""Latent-dimension complexity curves for frozen-manifest representation fits."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from manifest_representation_experiment_helpers import (
    ensure_manifests,
    save_grouped_line,
    save_simple_bar,
    write_csv,
)

from causal_spacetime_lab.state_change_manifest_representation import (
    ManifestRepresentationConfig,
    fit_manifest_dimension_curve,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for dimension-complexity diagnostics."""

    manifest_dir: Path = Path("outputs/manifests")
    candidate_dims: tuple[int, ...] = (1, 2, 3, 4, 5)
    steps: int = 600
    restarts: int = 2
    learning_rate: float = 0.05
    complexity_lambda: float = 0.01
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
        description="Frozen-manifest dimension-complexity curves."
    )
    parser.add_argument("--manifest-dir", type=Path, default=Path("outputs/manifests"))
    parser.add_argument(
        "--candidate-dims",
        nargs="+",
        default=["1", "2", "3", "4", "5"],
    )
    parser.add_argument("--steps", type=int, default=600)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--learning-rate", type=float, default=0.05)
    parser.add_argument("--complexity-lambda", type=float, default=0.01)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        manifest_dir=args.manifest_dir,
        candidate_dims=_parse_int_values(args.candidate_dims),
        steps=args.steps,
        restarts=args.restarts,
        learning_rate=args.learning_rate,
        complexity_lambda=args.complexity_lambda,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run representation complexity curves over latent dimensions."""

    loaded = ensure_manifests(config.manifest_dir, output_dir=config.output_dir)
    datasets = [dataset for dataset in loaded if dataset.eligible]
    rows: list[dict[str, float | str]] = []
    for dataset_index, dataset in enumerate(datasets):
        base_config = ManifestRepresentationConfig(
            embedding_dim=int(config.candidate_dims[0]),
            steps=config.steps,
            restarts=config.restarts,
            learning_rate=config.learning_rate,
            seed=config.seed + 100_000 * dataset_index,
        )
        fits = fit_manifest_dimension_curve(
            dataset,
            list(config.candidate_dims),
            base_config,
        )
        scores = [
            float(fit.heldout_violation_rate)
            + float(config.complexity_lambda) * fit.embedding_dim
            for fit in fits
        ]
        finite_scores = np.asarray(scores, dtype=float)
        best_index = int(np.nanargmin(finite_scores)) if finite_scores.size else -1
        best_dim = fits[best_index].embedding_dim if best_index >= 0 else -1
        for fit, score in zip(fits, scores, strict=True):
            rows.append(
                {
                    "manifest_id": dataset.manifest_id,
                    "embedding_dim": float(fit.embedding_dim),
                    "heldout_violation_rate": float(fit.heldout_violation_rate),
                    "train_violation_rate": float(fit.train_violation_rate),
                    "complexity_lambda": float(config.complexity_lambda),
                    "complexity_score": float(score),
                    "selected_dimension": float(best_dim),
                    "is_selected_dimension": float(fit.embedding_dim == best_dim),
                    "interpretation": (
                        "representation_complexity_not_physical_dimension"
                    ),
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write dimension-complexity CSV."""

    return write_csv(
        rows,
        output_dir / "data" / "frozen_manifest_dimension_complexity_curve.csv",
        [
            "manifest_id",
            "embedding_dim",
            "heldout_violation_rate",
            "train_violation_rate",
            "complexity_lambda",
            "complexity_score",
            "selected_dimension",
            "is_selected_dimension",
            "interpretation",
        ],
    )


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save dimension-complexity figures."""

    paths = [
        save_grouped_line(
            rows,
            x_key="embedding_dim",
            y_key="complexity_score",
            group_key="manifest_id",
            path=output_dir / "figures" / "frozen_manifest_complexity_curve.png",
            ylabel="held-out violation + complexity penalty",
        )
    ]
    selected = Counter(
        int(float(row["selected_dimension"]))
        for row in rows
        if float(row["is_selected_dimension"]) == 1.0
    )
    labels = [str(dim) for dim in sorted(selected)]
    values = [float(selected[int(label)]) for label in labels]
    paths.append(
        save_simple_bar(
            labels,
            values,
            output_dir / "figures" / "frozen_manifest_selected_dimension.png",
            ylabel="selected count",
        )
    )
    return paths


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote frozen manifest dimension-complexity curve: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
