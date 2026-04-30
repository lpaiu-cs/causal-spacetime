"""Null taxonomy diagnostics for preregistered v4 protocol manifests."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import save_family_metric_figure, write_csv

from causal_spacetime_lab.state_change_manifest_v4_diagnostics import (
    compute_v4_protocol_null_taxonomy_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    manifest_dir: Path = Path("outputs/manifests_v4")
    embedding_dim: int = 2
    null_repetitions: int = 2
    steps: int = 300
    restarts: int = 1
    seed: int = 0


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V4 protocol null taxonomy.")
    parser.add_argument(
        "--manifest-dir",
        type=Path,
        default=Path("outputs/manifests_v4"),
    )
    parser.add_argument("--embedding-dim", type=int, default=2)
    parser.add_argument("--null-repetitions", type=int, default=2)
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--restarts", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    return ExperimentConfig(
        args.manifest_dir,
        args.embedding_dim,
        args.null_repetitions,
        args.steps,
        args.restarts,
        args.seed,
    )


def main() -> None:
    config = parse_args()
    rows = compute_v4_protocol_null_taxonomy_rows(
        config.manifest_dir,
        embedding_dim=config.embedding_dim,
        null_repetitions=config.null_repetitions,
        steps=config.steps,
        restarts=config.restarts,
        seed=config.seed,
    )
    path = write_csv(
        rows,
        Path("outputs/data/v4_protocol_manifest_family_null_taxonomy.csv"),
        ["manifest_id", "family_name", "null_type"],
    )
    save_family_metric_figure(
        rows,
        metric="mean_heldout_violation_rate",
        path=Path("outputs/figures/v4_protocol_manifest_null_taxonomy_heldout_violation.png"),
        ylabel="Null heldout violation",
    )
    print(f"Wrote v4 protocol null taxonomy diagnostics: {path}")


if __name__ == "__main__":
    main()
