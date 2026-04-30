"""Restart and latent-order stability for protocol-invariant v3 manifests."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import save_family_metric_figure, write_csv

from causal_spacetime_lab.state_change_manifest_v3_protocol_diagnostics import (
    compute_v3_protocol_latent_order_stability_rows,
    compute_v3_protocol_restart_stability_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    manifest_dir: Path = Path("outputs/manifests_v3")
    embedding_dim: int = 2
    restart_count: int = 3
    steps: int = 300
    seed: int = 0


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 protocol stability diagnostics.")
    parser.add_argument(
        "--manifest-dir", type=Path, default=Path("outputs/manifests_v3")
    )
    parser.add_argument("--embedding-dim", type=int, default=2)
    parser.add_argument("--restart-count", type=int, default=3)
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    return ExperimentConfig(
        args.manifest_dir,
        args.embedding_dim,
        args.restart_count,
        args.steps,
        args.seed,
    )


def main() -> None:
    config = parse_args()
    restart_rows = compute_v3_protocol_restart_stability_rows(
        config.manifest_dir,
        embedding_dim=config.embedding_dim,
        restart_count=config.restart_count,
        steps=config.steps,
        seed=config.seed,
    )
    latent_rows = compute_v3_protocol_latent_order_stability_rows(
        config.manifest_dir,
        embedding_dim=config.embedding_dim,
        restart_count=config.restart_count,
        steps=config.steps,
        seed=config.seed,
    )
    restart_path = write_csv(
        restart_rows,
        Path("outputs/data/v3_protocol_manifest_family_restart_stability.csv"),
        ["manifest_id", "family_name", "restart_std"],
    )
    latent_path = write_csv(
        latent_rows,
        Path("outputs/data/v3_protocol_manifest_family_latent_order_stability.csv"),
        ["manifest_id", "family_name", "latent_order_disagreement"],
    )
    save_family_metric_figure(
        restart_rows,
        metric="restart_std",
        path=Path("outputs/figures/v3_protocol_manifest_restart_std.png"),
        ylabel="Restart std",
    )
    save_family_metric_figure(
        latent_rows,
        metric="latent_order_disagreement",
        path=Path("outputs/figures/v3_protocol_manifest_latent_order_disagreement.png"),
        ylabel="Latent-order disagreement",
    )
    print(f"Wrote v3 protocol stability diagnostics: {restart_path}, {latent_path}")


if __name__ == "__main__":
    main()
