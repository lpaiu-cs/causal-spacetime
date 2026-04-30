"""Produce v2 restart and latent-order stability diagnostics."""

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

from causal_spacetime_lab.state_change_manifest_v2_diagnostics import (
    compute_v2_latent_order_stability_rows,
    compute_v2_restart_stability_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 stability diagnostics."""

    manifest_dir: Path = Path("outputs/manifests_v2")
    embedding_dim: int = 2
    restart_count: int = 8
    steps: int = 500
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 restart stability diagnostics.")
    parser.add_argument(
        "--manifest-dir",
        type=Path,
        default=Path("outputs/manifests_v2"),
    )
    parser.add_argument("--embedding-dim", type=int, default=2)
    parser.add_argument("--restart-count", type=int, default=8)
    parser.add_argument("--steps", type=int, default=500)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        manifest_dir=args.manifest_dir,
        embedding_dim=args.embedding_dim,
        restart_count=args.restart_count,
        steps=args.steps,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    """Run restart and latent-order stability diagnostics."""

    restart_rows = compute_v2_restart_stability_rows(
        config.manifest_dir,
        embedding_dim=config.embedding_dim,
        restart_count=config.restart_count,
        steps=config.steps,
        seed=config.seed,
    )
    latent_rows = compute_v2_latent_order_stability_rows(
        config.manifest_dir,
        embedding_dim=config.embedding_dim,
        restart_count=config.restart_count,
        steps=config.steps,
        seed=config.seed,
    )
    return restart_rows, latent_rows


def write_outputs(
    restart_rows: list[dict[str, float | str]],
    latent_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write v2 stability diagnostics."""

    return (
        write_csv(
            restart_rows,
            data_path(output_dir, "v2_manifest_family_restart_stability.csv"),
            [
                "manifest_id",
                "family_name",
                "family_kind",
                "embedding_dim",
                "restart_count",
                "fit_count",
                "mean_heldout_violation_rate",
                "restart_std",
                "min_heldout_violation_rate",
                "max_heldout_violation_rate",
            ],
        ),
        write_csv(
            latent_rows,
            data_path(output_dir, "v2_manifest_family_latent_order_stability.csv"),
            [
                "manifest_id",
                "family_name",
                "family_kind",
                "embedding_dim",
                "restart_count",
                "fit_pair_count",
                "latent_order_disagreement",
                "max_pair_order_disagreement",
            ],
        ),
    )


def save_figures(
    restart_rows: list[dict[str, float | str]],
    latent_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> list[Path]:
    """Save v2 stability figures."""

    return [
        save_metric_bar(
            restart_rows,
            label_key="family_name",
            value_key="restart_std",
            path=figure_path(output_dir, "v2_manifest_restart_std.png"),
            ylabel="restart std",
        ),
        save_metric_bar(
            latent_rows,
            label_key="family_name",
            value_key="latent_order_disagreement",
            path=figure_path(
                output_dir,
                "v2_manifest_latent_order_disagreement.png",
            ),
            ylabel="latent-order disagreement",
        ),
    ]


def main() -> None:
    config = parse_args()
    restart_rows, latent_rows = run_experiment(config)
    paths = write_outputs(restart_rows, latent_rows, config.output_dir)
    figure_paths = save_figures(restart_rows, latent_rows, config.output_dir)
    print("Wrote v2 stability diagnostics: " + ", ".join(str(path) for path in paths))
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
