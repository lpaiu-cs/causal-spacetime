"""Physical-measure sprinkling for conformal interval-volume reconstruction."""

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
    integrate_conformal_interval_volume_1p1,
    sinusoidal_time_profile,
)
from causal_spacetime_lab.estimators import bias, relative_rmse, rmse
from causal_spacetime_lab.measure_sprinkling import (
    conformal_full_diamond_volume_1p1,
    sprinkle_conformal_measure_1p1,
)
from causal_spacetime_lab.metrics import causal_diamond_volume_1p1
from causal_spacetime_lab.probes import sample_probe_timelike_pairs_1p1

DEFAULT_N_VALUES = (600, 1200, 2400)
DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for physical-measure sprinkling validation."""

    T: float = 2.0
    n_values: tuple[int, ...] = DEFAULT_N_VALUES
    repetitions: int = 5
    pairs_per_repetition: int = 200
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
        description="Validate count reconstruction when measure is in sprinkling."
    )
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument(
        "--n-values",
        nargs="+",
        default=[str(value) for value in DEFAULT_N_VALUES],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--pairs-per-repetition", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_values=_parse_int_values(args.n_values, "N"),
        repetitions=args.repetitions,
        pairs_per_repetition=args.pairs_per_repetition,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def default_profiles(T: float) -> tuple[ConformalProfile, ...]:
    """Return default profiles."""

    return flat_profile(), constant_profile(1.5), sinusoidal_time_profile(0.3, T)


def _profile_label(profile: ConformalProfile) -> str:
    if profile.name == "constant":
        return f"constant_{profile.parameters['scale']:g}"
    if profile.name == "sinusoidal_time":
        return f"sinusoidal_{profile.parameters['amplitude']:g}"
    return profile.name


def _interval_count(events: np.ndarray, p: np.ndarray, q: np.ndarray) -> int:
    dt_from_p = events[:, 0] - p[0]
    dx_from_p = events[:, 1] - p[1]
    dt_to_q = q[0] - events[:, 0]
    dx_to_q = q[1] - events[:, 1]
    inside = (
        (dt_from_p > 0.0)
        & (dt_from_p * dt_from_p >= dx_from_p * dx_from_p)
        & (dt_to_q > 0.0)
        & (dt_to_q * dt_to_q >= dx_to_q * dx_to_q)
    )
    return int(np.count_nonzero(inside))


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    """Run physical-measure sprinkling reconstruction."""

    if config.T <= 0.0:
        raise ValueError("T must be positive")
    if any(value <= 0 for value in config.n_values):
        raise ValueError("all n_values must be positive")
    if config.repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if config.pairs_per_repetition <= 0:
        raise ValueError("pairs_per_repetition must be positive")

    rows: list[dict[str, float | str]] = []
    coordinate_volume = causal_diamond_volume_1p1(config.T)
    for profile_index, profile in enumerate(default_profiles(config.T)):
        profile_label = _profile_label(profile)
        physical_volume = conformal_full_diamond_volume_1p1(config.T, profile)
        for n_index, n_events in enumerate(config.n_values):
            physical_density = n_events / physical_volume
            coordinate_density = n_events / coordinate_volume
            for repetition in range(config.repetitions):
                support_events = sprinkle_conformal_measure_1p1(
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
                probe_pairs = sample_probe_timelike_pairs_1p1(
                    config.pairs_per_repetition,
                    T=config.T,
                    seed=(
                        config.seed
                        + 2_000_000 * profile_index
                        + 10_000 * n_index
                        + repetition
                    ),
                    min_tau=0.10,
                )
                for pair_index, (p, q) in enumerate(probe_pairs):
                    count = _interval_count(support_events, p, q)
                    true_volume = integrate_conformal_interval_volume_1p1(
                        p,
                        q,
                        profile,
                    )
                    physical_estimate = count / physical_density
                    coordinate_estimate = count / coordinate_density
                    rows.append(
                        {
                            "profile": profile_label,
                            "N": float(n_events),
                            "repetition": float(repetition),
                            "pair_index": float(pair_index),
                            "count": float(count),
                            "true_volume": true_volume,
                            "physical_density_estimate": physical_estimate,
                            "coordinate_density_estimate": coordinate_estimate,
                            "physical_density_error": physical_estimate - true_volume,
                            "coordinate_density_error": coordinate_estimate
                            - true_volume,
                        }
                    )
    return rows, summarize_results(rows)


def summarize_results(
    rows: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    """Aggregate pair-level physical-measure rows."""

    summary: list[dict[str, float | str]] = []
    keys = sorted({(str(row["profile"]), int(row["N"])) for row in rows})
    for profile, n_events in keys:
        subset = [
            row
            for row in rows
            if row["profile"] == profile and int(row["N"]) == n_events
        ]
        true = np.asarray([row["true_volume"] for row in subset], dtype=float)
        physical = np.asarray(
            [row["physical_density_estimate"] for row in subset],
            dtype=float,
        )
        coordinate = np.asarray(
            [row["coordinate_density_estimate"] for row in subset],
            dtype=float,
        )
        summary.append(
            {
                "profile": profile,
                "N": float(n_events),
                "repetitions": float(len({int(row["repetition"]) for row in subset})),
                "pair_count": float(len(subset)),
                "true_volume_mean": float(np.mean(true)),
                "physical_density_rmse": rmse(true, physical),
                "physical_density_bias": bias(true, physical),
                "physical_density_relative_rmse": relative_rmse(true, physical),
                "coordinate_density_rmse": rmse(true, coordinate),
                "coordinate_density_bias": bias(true, coordinate),
                "coordinate_density_relative_rmse": relative_rmse(true, coordinate),
            }
        )
    return summary


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
    rows: list[dict[str, float | str]],
    summary: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write physical-measure sprinkling outputs."""

    data_dir = output_dir / "data"
    return (
        _write_csv(rows, data_dir / "physical_measure_sprinkling_pairs.csv"),
        _write_csv(summary, data_dir / "physical_measure_sprinkling_summary.csv"),
    )


def save_scatter(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Save true versus estimated physical volume scatter."""

    output_path = output_dir / "figures" / "physical_measure_volume_scatter.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    truth = np.asarray([row["true_volume"] for row in rows], dtype=float)
    physical = np.asarray(
        [row["physical_density_estimate"] for row in rows],
        dtype=float,
    )
    coordinate = np.asarray(
        [row["coordinate_density_estimate"] for row in rows],
        dtype=float,
    )
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    ax.scatter(truth, physical, s=6, alpha=0.25, label="physical density")
    ax.scatter(truth, coordinate, s=6, alpha=0.25, label="coordinate density")
    upper = float(max(np.max(truth), np.max(physical), np.max(coordinate)))
    ax.plot([0.0, upper], [0.0, upper], color="black", linewidth=1.0)
    ax.set_xlabel("True conformal physical volume")
    ax.set_ylabel("Estimated volume")
    ax.set_title("Physical-measure sprinkling volume reconstruction")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_rmse_vs_n(
    summary: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Save RMSE versus N by profile."""

    output_path = output_dir / "figures" / "physical_measure_rmse_vs_N.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for profile in sorted({row["profile"] for row in summary}):
        subset = [row for row in summary if row["profile"] == profile]
        n_values = np.asarray([row["N"] for row in subset], dtype=float)
        physical = np.asarray([row["physical_density_rmse"] for row in subset])
        coordinate = np.asarray([row["coordinate_density_rmse"] for row in subset])
        ax.plot(n_values, physical, marker="o", label=f"{profile} physical")
        ax.plot(n_values, coordinate, marker="s", label=f"{profile} coordinate")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Support event count N")
    ax.set_ylabel("Volume RMSE")
    ax.set_title("Physical-measure sprinkling RMSE")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small", ncol=2)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_bias_by_profile(
    summary: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Save bias by profile."""

    output_path = output_dir / "figures" / "physical_measure_bias_by_profile.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    profiles = sorted({str(row["profile"]) for row in summary})
    physical = []
    coordinate = []
    for profile in profiles:
        subset = [row for row in summary if row["profile"] == profile]
        physical.append(
            float(np.mean([row["physical_density_bias"] for row in subset]))
        )
        coordinate.append(
            float(np.mean([row["coordinate_density_bias"] for row in subset]))
        )
    x = np.arange(len(profiles))
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    ax.bar(x - 0.18, physical, width=0.36, label="physical density")
    ax.bar(x + 0.18, coordinate, width=0.36, label="coordinate density")
    ax.axhline(0.0, color="black", linewidth=0.9)
    ax.set_xticks(x, profiles, rotation=20)
    ax.set_ylabel("Volume bias")
    ax.set_title("Physical-measure sprinkling bias")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_figures(
    rows: list[dict[str, float | str]],
    summary: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path, Path]:
    """Save physical-measure sprinkling figures."""

    return (
        save_scatter(rows, output_dir),
        save_rmse_vs_n(summary, output_dir),
        save_bias_by_profile(summary, output_dir),
    )


def main() -> None:
    config = parse_args()
    rows, summary = run_experiment(config)
    pair_path, summary_path = write_outputs(rows, summary, config.output_dir)
    figure_paths = save_figures(rows, summary, config.output_dir)
    print(f"Wrote pair results: {pair_path}")
    print(f"Wrote summary: {summary_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
