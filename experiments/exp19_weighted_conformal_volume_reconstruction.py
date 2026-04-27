"""Weighted conformal interval-volume reconstruction experiment."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.conformal import (
    ConformalProfile,
    conformal_volume_density_1p1,
    constant_profile,
    flat_profile,
    integrate_conformal_interval_volume_1p1,
    sinusoidal_time_profile,
)
from causal_spacetime_lab.estimators import bias, relative_rmse, rmse
from causal_spacetime_lab.metrics import causal_diamond_volume_1p1
from causal_spacetime_lab.probes import sample_probe_timelike_pairs_1p1
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_N_VALUES = (600, 1200, 2400)
DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for weighted conformal volume reconstruction."""

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
        description="Compare unweighted and supplied-weight conformal volumes."
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
    """Return default conformal profiles for the experiment."""

    return (
        flat_profile(),
        constant_profile(1.5),
        sinusoidal_time_profile(0.3, T),
    )


def _profile_label(profile: ConformalProfile) -> str:
    if profile.name == "constant":
        return f"constant_{profile.parameters['scale']:g}"
    if profile.name == "sinusoidal_time":
        return f"sinusoidal_{profile.parameters['amplitude']:g}"
    return profile.name


def _interval_mask(events: np.ndarray, p: np.ndarray, q: np.ndarray) -> np.ndarray:
    dt_from_p = events[:, 0] - p[0]
    dx_from_p = events[:, 1] - p[1]
    dt_to_q = q[0] - events[:, 0]
    dx_to_q = q[1] - events[:, 1]
    return (
        (dt_from_p > 0.0)
        & (dt_from_p * dt_from_p >= dx_from_p * dx_from_p)
        & (dt_to_q > 0.0)
        & (dt_to_q * dt_to_q >= dx_to_q * dx_to_q)
    )


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    """Run weighted conformal volume reconstruction."""

    if config.T <= 0.0:
        raise ValueError("T must be positive")
    if any(value <= 0 for value in config.n_values):
        raise ValueError("all n_values must be positive")
    if config.repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if config.pairs_per_repetition <= 0:
        raise ValueError("pairs_per_repetition must be positive")

    profiles = default_profiles(config.T)
    pair_rows: list[dict[str, float | str]] = []
    global_volume = causal_diamond_volume_1p1(config.T)

    for n_index, n_events in enumerate(config.n_values):
        coordinate_density = n_events / global_volume
        for repetition in range(config.repetitions):
            support_events = sprinkle_1p1_causal_diamond(
                n_events,
                T=config.T,
                seed=config.seed + 10_000 * n_index + repetition,
            )
            probe_pairs = sample_probe_timelike_pairs_1p1(
                config.pairs_per_repetition,
                T=config.T,
                seed=config.seed + 100_000 * n_index + 1_000 * repetition,
                min_tau=0.10,
            )
            for profile in profiles:
                weights = conformal_volume_density_1p1(profile, support_events)
                for pair_index, (p, q) in enumerate(probe_pairs):
                    inside = _interval_mask(support_events, p, q)
                    true_volume = integrate_conformal_interval_volume_1p1(
                        p,
                        q,
                        profile,
                    )
                    unweighted = np.count_nonzero(inside) / coordinate_density
                    weighted = float(np.sum(weights[inside]) / coordinate_density)
                    pair_rows.append(
                        {
                            "profile": _profile_label(profile),
                            "N": float(n_events),
                            "repetition": float(repetition),
                            "pair_index": float(pair_index),
                            "true_volume": true_volume,
                            "unweighted_volume_estimate": float(unweighted),
                            "weighted_volume_estimate": weighted,
                            "unweighted_error": float(unweighted - true_volume),
                            "weighted_error": float(weighted - true_volume),
                        }
                    )

    return pair_rows, summarize_results(pair_rows)


def summarize_results(
    pair_rows: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    """Aggregate pair-level conformal volume results."""

    summary_rows: list[dict[str, float | str]] = []
    keys = sorted({(str(row["profile"]), int(row["N"])) for row in pair_rows})
    for profile, n_events in keys:
        rows = [
            row
            for row in pair_rows
            if row["profile"] == profile and int(row["N"]) == n_events
        ]
        true = np.asarray([row["true_volume"] for row in rows], dtype=float)
        unweighted = np.asarray(
            [row["unweighted_volume_estimate"] for row in rows],
            dtype=float,
        )
        weighted = np.asarray(
            [row["weighted_volume_estimate"] for row in rows],
            dtype=float,
        )
        summary_rows.append(
            {
                "profile": profile,
                "N": float(n_events),
                "repetitions": float(len({int(row["repetition"]) for row in rows})),
                "pair_count": float(len(rows)),
                "true_volume_mean": float(np.mean(true)),
                "unweighted_volume_rmse": rmse(true, unweighted),
                "weighted_volume_rmse": rmse(true, weighted),
                "unweighted_volume_bias": bias(true, unweighted),
                "weighted_volume_bias": bias(true, weighted),
                "unweighted_relative_rmse": relative_rmse(true, unweighted),
                "weighted_relative_rmse": relative_rmse(true, weighted),
            }
        )
    return summary_rows


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
    pair_rows: list[dict[str, float | str]],
    summary_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write weighted conformal volume outputs."""

    data_dir = output_dir / "data"
    pair_path = data_dir / "weighted_conformal_volume_pairs.csv"
    summary_path = data_dir / "weighted_conformal_volume_summary.csv"
    return _write_csv(pair_rows, pair_path), _write_csv(summary_rows, summary_path)


def save_scatter(pair_rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Save true versus estimated conformal volume scatter."""

    output_path = output_dir / "figures" / "weighted_conformal_volume_scatter.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    truth = np.asarray([row["true_volume"] for row in pair_rows], dtype=float)
    unweighted = np.asarray(
        [row["unweighted_volume_estimate"] for row in pair_rows],
        dtype=float,
    )
    weighted = np.asarray(
        [row["weighted_volume_estimate"] for row in pair_rows],
        dtype=float,
    )
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    ax.scatter(truth, unweighted, s=6, alpha=0.25, label="unweighted")
    ax.scatter(truth, weighted, s=6, alpha=0.25, label="weighted")
    lower = float(np.min(truth))
    upper = float(max(np.max(truth), np.max(unweighted), np.max(weighted)))
    ax.plot([lower, upper], [lower, upper], color="black", linewidth=1.0)
    ax.set_xlabel("True conformal interval volume")
    ax.set_ylabel("Estimated volume")
    ax.set_title("Weighted conformal volume reconstruction")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_rmse_vs_n(
    summary_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Save RMSE versus N by profile and estimator."""

    output_path = output_dir / "figures" / "weighted_conformal_volume_rmse_vs_N.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for profile in sorted({row["profile"] for row in summary_rows}):
        rows = [row for row in summary_rows if row["profile"] == profile]
        n_values = np.asarray([row["N"] for row in rows], dtype=float)
        unweighted = np.asarray(
            [row["unweighted_volume_rmse"] for row in rows],
            dtype=float,
        )
        weighted = np.asarray(
            [row["weighted_volume_rmse"] for row in rows],
            dtype=float,
        )
        ax.plot(n_values, unweighted, marker="o", label=f"{profile} unweighted")
        ax.plot(n_values, weighted, marker="s", label=f"{profile} weighted")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Support event count N")
    ax.set_ylabel("Volume RMSE")
    ax.set_title("Conformal volume RMSE")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small", ncol=2)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_bias_by_profile(
    summary_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Save estimator bias by profile."""

    output_path = (
        output_dir / "figures" / "weighted_conformal_volume_bias_by_profile.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    profiles = sorted({str(row["profile"]) for row in summary_rows})
    unweighted = []
    weighted = []
    for profile in profiles:
        rows = [row for row in summary_rows if row["profile"] == profile]
        unweighted.append(
            float(np.mean([row["unweighted_volume_bias"] for row in rows]))
        )
        weighted.append(float(np.mean([row["weighted_volume_bias"] for row in rows])))
    x = np.arange(len(profiles))
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    ax.bar(x - 0.18, unweighted, width=0.36, label="unweighted")
    ax.bar(x + 0.18, weighted, width=0.36, label="weighted")
    ax.axhline(0.0, color="black", linewidth=0.9)
    ax.set_xticks(x, profiles, rotation=20)
    ax.set_ylabel("Volume bias")
    ax.set_title("Conformal volume bias by supplied profile")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_figures(
    pair_rows: list[dict[str, float | str]],
    summary_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path, Path]:
    """Save weighted conformal volume figures."""

    return (
        save_scatter(pair_rows, output_dir),
        save_rmse_vs_n(summary_rows, output_dir),
        save_bias_by_profile(summary_rows, output_dir),
    )


def main() -> None:
    config = parse_args()
    pair_rows, summary_rows = run_experiment(config)
    pair_path, summary_path = write_outputs(pair_rows, summary_rows, config.output_dir)
    figure_paths = save_figures(pair_rows, summary_rows, config.output_dir)
    print(f"Wrote pair results: {pair_path}")
    print(f"Wrote summary: {summary_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
