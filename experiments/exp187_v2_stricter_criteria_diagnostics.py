"""Apply fixed stricter diagnostic criteria to v2 fit rows."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import (
    DEFAULT_OUTPUT_DIR,
    data_path,
    figure_path,
    parse_int_values,
    save_metric_bar,
)

from causal_spacetime_lab.state_change_manifest_v2_diagnostics import (
    compute_v2_stricter_criteria_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 stricter criteria diagnostics."""

    manifest_dir: Path = Path("outputs/manifests_v2")
    dims: tuple[int, ...] = (1, 2, 3)
    steps: int = 600
    restarts: int = 2
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 stricter criteria diagnostics.")
    parser.add_argument(
        "--manifest-dir",
        type=Path,
        default=Path("outputs/manifests_v2"),
    )
    parser.add_argument("--dims", nargs="+", default=["1", "2", "3"])
    parser.add_argument("--steps", type=int, default=600)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        manifest_dir=args.manifest_dir,
        dims=parse_int_values(args.dims),
        steps=args.steps,
        restarts=args.restarts,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run fixed stricter criteria diagnostics."""

    return compute_v2_stricter_criteria_rows(
        config.manifest_dir,
        dims=list(config.dims),
        steps=config.steps,
        restarts=config.restarts,
        seed=config.seed,
    )


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write v2 stricter criteria rows."""

    return write_csv(
        rows,
        data_path(output_dir, "v2_manifest_family_stricter_criteria.csv"),
        [
            "family_name",
            "family_kind",
            "embedding_dim",
            "threshold_pass",
            "heldout_violation_threshold",
            "generalization_gap_threshold",
        ],
    )


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> list[Path]:
    """Save v2 stricter criteria figure."""

    return [
        save_metric_bar(
            rows,
            label_key="family_name",
            value_key="threshold_pass",
            path=figure_path(output_dir, "v2_manifest_family_stricter_pass_rate.png"),
            ylabel="fixed stricter pass indicator",
        )
    ]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    output_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote v2 stricter criteria diagnostics: {output_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
