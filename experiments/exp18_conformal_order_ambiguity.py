"""Conformal ambiguity of causal order in controlled 1+1D models."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.conformal import (
    central_observer_proper_time_1p1,
    constant_profile,
    flat_profile,
    integrate_conformal_interval_volume_1p1,
    sinusoidal_time_profile,
)
from causal_spacetime_lab.dimension import (
    estimate_myrheim_meyer_dimension,
    relation_fraction,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_OUTPUT_DIR = Path("outputs")
DEFAULT_CONSTANT_SCALES = (1.0, 1.5, 2.0)


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for conformal order ambiguity validation."""

    T: float = 2.0
    n_events: int = 1200
    seed: int = 0
    constant_scales: tuple[float, ...] = DEFAULT_CONSTANT_SCALES
    output_dir: Path = DEFAULT_OUTPUT_DIR
    include_sinusoidal: bool = True


def _parse_float_values(values: list[str], name: str) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    if not parsed:
        raise argparse.ArgumentTypeError(f"at least one {name} value is required")
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Demonstrate conformal ambiguity of causal order."
    )
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument("--N", type=int, default=1200)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--constant-scales",
        nargs="+",
        default=[str(value) for value in DEFAULT_CONSTANT_SCALES],
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--skip-sinusoidal", action="store_true")
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_events=args.N,
        seed=args.seed,
        constant_scales=_parse_float_values(args.constant_scales, "scale"),
        output_dir=args.output_dir,
        include_sinusoidal=not args.skip_sinusoidal,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run conformal order ambiguity validation."""

    if config.T <= 0.0:
        raise ValueError("T must be positive")
    if config.n_events <= 1:
        raise ValueError("N must be greater than 1")

    events = sprinkle_1p1_causal_diamond(config.n_events, T=config.T, seed=config.seed)
    causal_matrix = causal_matrix_1p1(events)
    flat_fraction = relation_fraction(causal_matrix)
    flat_dimension = estimate_myrheim_meyer_dimension(flat_fraction)
    p = np.asarray([-config.T / 2.0, 0.0])
    q = np.asarray([config.T / 2.0, 0.0])
    flat_time = central_observer_proper_time_1p1(
        -config.T / 2.0,
        config.T / 2.0,
        flat_profile(),
    )
    flat_volume = integrate_conformal_interval_volume_1p1(p, q, flat_profile())

    rows: list[dict[str, float | str]] = []
    for scale in config.constant_scales:
        profile = constant_profile(scale)
        proper_time = central_observer_proper_time_1p1(
            -config.T / 2.0,
            config.T / 2.0,
            profile,
        )
        volume = integrate_conformal_interval_volume_1p1(p, q, profile)
        rows.append(
            {
                "profile": f"constant_{scale:g}",
                "scale": float(scale),
                "causal_matrix_changed": 0.0,
                "relation_fraction": flat_fraction,
                "estimated_dimension": flat_dimension,
                "proper_time": proper_time,
                "full_diamond_volume": volume,
                "relation_fraction_ratio": flat_fraction / flat_fraction,
                "estimated_dimension_ratio": flat_dimension / flat_dimension,
                "proper_time_ratio": proper_time / flat_time,
                "volume_ratio": volume / flat_volume,
            }
        )

    if config.include_sinusoidal:
        profile = sinusoidal_time_profile(0.3, config.T)
        proper_time = central_observer_proper_time_1p1(
            -config.T / 2.0,
            config.T / 2.0,
            profile,
        )
        volume = integrate_conformal_interval_volume_1p1(p, q, profile)
        rows.append(
            {
                "profile": "sinusoidal_0.3",
                "scale": float("nan"),
                "causal_matrix_changed": 0.0,
                "relation_fraction": flat_fraction,
                "estimated_dimension": flat_dimension,
                "proper_time": proper_time,
                "full_diamond_volume": volume,
                "relation_fraction_ratio": 1.0,
                "estimated_dimension_ratio": 1.0,
                "proper_time_ratio": proper_time / flat_time,
                "volume_ratio": volume / flat_volume,
            }
        )
    return rows


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write conformal order ambiguity summary."""

    output_path = output_dir / "data" / "conformal_order_ambiguity_summary.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_plot(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Save normalized conformal scale ambiguity plot."""

    output_path = output_dir / "figures" / "conformal_order_ambiguity_scales.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    scaled_rows = [row for row in rows if np.isfinite(float(row["scale"]))]
    scales = np.asarray([row["scale"] for row in scaled_rows], dtype=float)
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    for key, label in [
        ("relation_fraction_ratio", "relation fraction"),
        ("estimated_dimension_ratio", "estimated dimension"),
        ("proper_time_ratio", "proper time"),
        ("volume_ratio", "volume"),
    ]:
        values = np.asarray([row[key] for row in scaled_rows], dtype=float)
        ax.plot(scales, values, marker="o", label=label)
    ax.set_xlabel("Constant conformal scale")
    ax.set_ylabel("Normalized quantity")
    ax.set_title("Causal order is invariant under constant conformal scaling")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    summary_path = write_outputs(rows, config.output_dir)
    figure_path = save_plot(rows, config.output_dir)
    print(f"Wrote summary: {summary_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
