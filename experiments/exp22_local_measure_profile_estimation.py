"""Local relative measure-shape estimation from physical-volume sprinkling."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.conformal import (
    ConformalProfile,
    constant_profile,
    flat_profile,
    linear_time_profile,
    sinusoidal_time_profile,
)
from causal_spacetime_lab.estimators import bias, mae, rmse
from causal_spacetime_lab.measure_sprinkling import (
    conformal_time_bin_masses_1p1,
    coordinate_time_bin_volumes_1p1,
    estimate_time_bin_density_shape_1p1,
    normalize_profile_shape,
    sprinkle_conformal_measure_1p1,
)

DEFAULT_N_VALUES = (600, 1200, 2400, 4800)
DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for local measure-shape estimation."""

    T: float = 2.0
    n_values: tuple[int, ...] = DEFAULT_N_VALUES
    repetitions: int = 5
    num_bins: int = 16
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str], name: str) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    if not parsed:
        raise argparse.ArgumentTypeError(f"at least one {name} value is required")
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Estimate local relative measure shape from event density."
    )
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument(
        "--n-values",
        nargs="+",
        default=[str(value) for value in DEFAULT_N_VALUES],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--num-bins", type=int, default=16)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_values=_parse_int_values(args.n_values, "N"),
        repetitions=args.repetitions,
        num_bins=args.num_bins,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def default_profiles(T: float) -> tuple[ConformalProfile, ...]:
    """Return default profiles."""

    return (
        flat_profile(),
        constant_profile(1.5),
        sinusoidal_time_profile(0.3, T),
        linear_time_profile(0.3, T),
    )


def _profile_label(profile: ConformalProfile) -> str:
    if profile.name == "constant":
        return f"constant_{profile.parameters['scale']:g}"
    if profile.name == "sinusoidal_time":
        return f"sinusoidal_{profile.parameters['amplitude']:g}"
    if profile.name == "linear_time":
        return f"linear_{profile.parameters['amplitude']:g}"
    return profile.name


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    """Run local measure-shape estimation."""

    if config.T <= 0.0:
        raise ValueError("T must be positive")
    if any(value <= 0 for value in config.n_values):
        raise ValueError("all n_values must be positive")
    if config.repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if config.num_bins < 2:
        raise ValueError("num_bins must be at least 2")

    edges = np.linspace(-config.T / 2.0, config.T / 2.0, config.num_bins + 1)
    coord_volumes = coordinate_time_bin_volumes_1p1(config.T, edges)
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    bin_rows: list[dict[str, float | str]] = []
    summary_rows: list[dict[str, float | str]] = []

    for profile_index, profile in enumerate(default_profiles(config.T)):
        profile_label = _profile_label(profile)
        masses = conformal_time_bin_masses_1p1(config.T, profile, edges)
        true_density_shape = normalize_profile_shape(masses / coord_volumes)
        for n_index, n_events in enumerate(config.n_values):
            estimates: list[np.ndarray] = []
            for repetition in range(config.repetitions):
                events = sprinkle_conformal_measure_1p1(
                    n_events,
                    config.T,
                    profile,
                    seed=(
                        config.seed
                        + 1_000_000 * profile_index
                        + 10_000 * n_index
                        + repetition
                    ),
                )
                estimated_density = estimate_time_bin_density_shape_1p1(
                    events,
                    config.T,
                    edges,
                )
                estimated_shape = normalize_profile_shape(estimated_density)
                estimates.append(estimated_shape)
                for bin_index in range(config.num_bins):
                    bin_rows.append(
                        {
                            "profile": profile_label,
                            "N": float(n_events),
                            "repetition": float(repetition),
                            "bin_index": float(bin_index),
                            "bin_center": float(bin_centers[bin_index]),
                            "estimated_shape": float(estimated_shape[bin_index]),
                            "true_shape": float(true_density_shape[bin_index]),
                            "shape_error": float(
                                estimated_shape[bin_index]
                                - true_density_shape[bin_index]
                            ),
                        }
                    )
            estimate_values = np.concatenate(estimates)
            true_values = np.tile(true_density_shape, config.repetitions)
            finite = np.isfinite(estimate_values) & np.isfinite(true_values)
            finite_true = true_values[finite]
            finite_estimate = estimate_values[finite]
            shape_rmse = (
                rmse(finite_true, finite_estimate) if finite_true.size else np.nan
            )
            shape_mae = (
                mae(finite_true, finite_estimate) if finite_true.size else np.nan
            )
            shape_bias = (
                bias(finite_true, finite_estimate) if finite_true.size else np.nan
            )
            summary_rows.append(
                {
                    "profile": profile_label,
                    "N": float(n_events),
                    "repetitions": float(config.repetitions),
                    "num_bins": float(config.num_bins),
                    "shape_rmse": shape_rmse,
                    "shape_mae": shape_mae,
                    "shape_bias": shape_bias,
                }
            )
    return bin_rows, summary_rows


def _write_csv(rows: list[dict[str, float | str]], output_path: Path) -> Path:
    if not rows:
        raise RuntimeError("no rows to write")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def write_outputs(
    bin_rows: list[dict[str, float | str]],
    summary_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write local measure-profile outputs."""

    data_dir = output_dir / "data"
    return (
        _write_csv(bin_rows, data_dir / "local_measure_profile_bins.csv"),
        _write_csv(summary_rows, data_dir / "local_measure_profile_summary.csv"),
    )


def save_shape_plot(
    bin_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Save estimated versus true normalized shape for max N."""

    output_path = output_dir / "figures" / "local_measure_profile_shape.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    max_n = max(int(row["N"]) for row in bin_rows)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for profile in sorted({row["profile"] for row in bin_rows}):
        rows = [
            row
            for row in bin_rows
            if row["profile"] == profile and int(row["N"]) == max_n
        ]
        centers = sorted({float(row["bin_center"]) for row in rows})
        estimated = []
        truth = []
        for center in centers:
            center_rows = [row for row in rows if float(row["bin_center"]) == center]
            estimated.append(
                float(np.nanmean([row["estimated_shape"] for row in center_rows]))
            )
            truth.append(float(np.nanmean([row["true_shape"] for row in center_rows])))
        ax.plot(centers, truth, linestyle="--", label=f"{profile} true")
        ax.plot(centers, estimated, marker="o", label=f"{profile} estimated")
    ax.set_xlabel("Coordinate time bin center")
    ax.set_ylabel("Normalized density shape")
    ax.set_title("Local relative measure shape")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small", ncol=2)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_rmse_vs_n(
    summary_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Save shape RMSE versus N."""

    output_path = output_dir / "figures" / "local_measure_profile_rmse_vs_N.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    for profile in sorted({row["profile"] for row in summary_rows}):
        rows = [row for row in summary_rows if row["profile"] == profile]
        n_values = np.asarray([row["N"] for row in rows], dtype=float)
        values = np.asarray([row["shape_rmse"] for row in rows], dtype=float)
        ax.plot(n_values, values, marker="o", label=profile)
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Support event count N")
    ax.set_ylabel("Shape RMSE")
    ax.set_title("Local measure shape recovery")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_figures(
    bin_rows: list[dict[str, float | str]],
    summary_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Save local profile figures."""

    return (
        save_shape_plot(bin_rows, output_dir),
        save_rmse_vs_n(summary_rows, output_dir),
    )


def main() -> None:
    config = parse_args()
    bin_rows, summary_rows = run_experiment(config)
    bin_path, summary_path = write_outputs(bin_rows, summary_rows, config.output_dir)
    figure_paths = save_figures(bin_rows, summary_rows, config.output_dir)
    print(f"Wrote bin results: {bin_path}")
    print(f"Wrote summary: {summary_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
