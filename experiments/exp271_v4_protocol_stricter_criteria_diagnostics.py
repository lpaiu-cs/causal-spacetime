"""Stricter diagnostic criteria rows for v4 protocol manifests."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import save_family_metric_figure, write_csv

from causal_spacetime_lab.state_change_manifest_v4_diagnostics import (
    compute_v4_protocol_stricter_criteria_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    manifest_dir: Path = Path("outputs/manifests_v4")
    dims: list[int] | None = None
    steps: int = 300
    restarts: int = 1
    seed: int = 0


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V4 stricter criteria diagnostics.")
    parser.add_argument(
        "--manifest-dir",
        type=Path,
        default=Path("outputs/manifests_v4"),
    )
    parser.add_argument("--dims", type=int, nargs="+", default=[1, 2])
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--restarts", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    return ExperimentConfig(
        args.manifest_dir,
        args.dims,
        args.steps,
        args.restarts,
        args.seed,
    )


def main() -> None:
    config = parse_args()
    rows = compute_v4_protocol_stricter_criteria_rows(
        config.manifest_dir,
        dims=config.dims or [1, 2],
        steps=config.steps,
        restarts=config.restarts,
        seed=config.seed,
    )
    path = write_csv(
        rows,
        Path("outputs/data/v4_protocol_manifest_family_stricter_criteria.csv"),
        ["family_name", "family_kind", "embedding_dim", "threshold_pass"],
    )
    save_family_metric_figure(
        rows,
        metric="threshold_pass",
        path=Path("outputs/figures/v4_protocol_manifest_family_stricter_pass_rate.png"),
        ylabel="Stricter pass",
    )
    print(f"Wrote v4 protocol stricter criteria diagnostics: {path}")


if __name__ == "__main__":
    main()
