"""Rindler horizon reconstruction-inaccessibility validation."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.discrete_radar import discrete_radar_coordinates_from_order
from causal_spacetime_lab.horizon_validation import (
    finite_coverage_mask_from_analytic_ticks,
    horizon_accessibility_confusion,
    radar_error_summary,
)
from causal_spacetime_lab.observer import observer_chain_indices
from causal_spacetime_lab.rindler import (
    analytic_rindler_radar_coordinates_1p1,
    analytic_rindler_radar_ticks_1p1,
    make_rindler_observer_chain_1p1,
    safe_tau_range_for_rindler_chain_1p1,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_ACCELERATIONS = (1.5, 2.0)
DEFAULT_N_VALUES = (600, 1200, 2400)
DEFAULT_TICK_VALUES = (32, 64, 128)
DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for Rindler horizon reconstruction validation."""

    T: float = 4.0
    acceleration_values: tuple[float, ...] = DEFAULT_ACCELERATIONS
    n_values: tuple[int, ...] = DEFAULT_N_VALUES
    tick_values: tuple[int, ...] = DEFAULT_TICK_VALUES
    repetitions: int = 5
    seed: int = 0
    direction: int = 1
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str], name: str) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    if not parsed:
        raise argparse.ArgumentTypeError(f"at least one {name} value is required")
    return tuple(parsed)


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
        description=(
            "Validate Rindler reconstruction-inaccessibility from causal order "
            "and a supplied accelerated observer clock protocol."
        )
    )
    parser.add_argument("--T", type=float, default=4.0)
    parser.add_argument(
        "--accelerations",
        nargs="+",
        default=[str(value) for value in DEFAULT_ACCELERATIONS],
    )
    parser.add_argument(
        "--n-values",
        nargs="+",
        default=[str(value) for value in DEFAULT_N_VALUES],
    )
    parser.add_argument(
        "--tick-values",
        nargs="+",
        default=[str(value) for value in DEFAULT_TICK_VALUES],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--direction", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        acceleration_values=_parse_float_values(args.accelerations, "acceleration"),
        n_values=_parse_int_values(args.n_values, "N"),
        tick_values=_parse_int_values(args.tick_values, "tick"),
        repetitions=args.repetitions,
        seed=args.seed,
        direction=args.direction,
        output_dir=args.output_dir,
    )


def _validate_config(config: ExperimentConfig) -> None:
    if config.T <= 0.0:
        raise ValueError("T must be positive")
    if any(accel <= 0.0 for accel in config.acceleration_values):
        raise ValueError("all accelerations must be positive")
    if any(n <= 0 for n in config.n_values):
        raise ValueError("all n_values must be positive")
    if any(ticks < 2 for ticks in config.tick_values):
        raise ValueError("all tick_values must be at least 2")
    if config.repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if config.direction not in {-1, 1}:
        raise ValueError("direction must be -1 or 1")


def _value_or_nan(value: float | None) -> float:
    return float("nan") if value is None else float(value)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float]], list[dict[str, float]]]:
    """Run Rindler horizon reconstruction validation."""

    _validate_config(config)
    event_rows: list[dict[str, float]] = []

    for accel_index, acceleration in enumerate(config.acceleration_values):
        tau_min, tau_max = safe_tau_range_for_rindler_chain_1p1(
            config.T,
            acceleration,
            direction=config.direction,
        )
        for n_index, n_events in enumerate(config.n_values):
            for repetition in range(config.repetitions):
                support_events = sprinkle_1p1_causal_diamond(
                    n_events,
                    T=config.T,
                    seed=(
                        config.seed
                        + 100_000 * accel_index
                        + 10_000 * n_index
                        + repetition
                    ),
                )
                tau_minus_true, tau_plus_true, wedge_accessible = (
                    analytic_rindler_radar_ticks_1p1(
                        support_events,
                        acceleration,
                        direction=config.direction,
                    )
                )
                radar_time_true, radar_distance_true, _ = (
                    analytic_rindler_radar_coordinates_1p1(
                        support_events,
                        acceleration,
                        direction=config.direction,
                    )
                )
                for tick_count in config.tick_values:
                    observer_events, clock_times = make_rindler_observer_chain_1p1(
                        acceleration,
                        tick_count,
                        tau_min,
                        tau_max,
                        direction=config.direction,
                    )
                    combined_events = np.vstack((support_events, observer_events))
                    causal_matrix = causal_matrix_1p1(combined_events)
                    target_indices = np.arange(n_events, dtype=int)
                    observer_indices = observer_chain_indices(n_events, tick_count)
                    reconstructions = discrete_radar_coordinates_from_order(
                        causal_matrix,
                        observer_indices,
                        target_indices,
                        clock_times,
                    )
                    finite_coverage = finite_coverage_mask_from_analytic_ticks(
                        tau_minus_true,
                        tau_plus_true,
                        float(clock_times[0]),
                        float(clock_times[-1]),
                    )

                    for target_index, result in enumerate(reconstructions):
                        reconstructed_accessible = bool(result.accessible)
                        time_hat = _value_or_nan(result.radar_time)
                        distance_hat = _value_or_nan(result.radar_distance)
                        event_rows.append(
                            {
                                "N": float(n_events),
                                "acceleration": float(acceleration),
                                "tick_count": float(tick_count),
                                "repetition": float(repetition),
                                "target_index": float(target_index),
                                "t": float(support_events[target_index, 0]),
                                "x": float(support_events[target_index, 1]),
                                "analytic_wedge_accessible": float(
                                    wedge_accessible[target_index]
                                ),
                                "finite_coverage_accessible": float(
                                    finite_coverage[target_index]
                                ),
                                "reconstructed_accessible": float(
                                    reconstructed_accessible
                                ),
                                "tau_minus_hat": _value_or_nan(result.tau_minus),
                                "tau_plus_hat": _value_or_nan(result.tau_plus),
                                "tau_minus_true": float(
                                    tau_minus_true[target_index]
                                ),
                                "tau_plus_true": float(tau_plus_true[target_index]),
                                "radar_time_hat": time_hat,
                                "radar_distance_hat": distance_hat,
                                "radar_time_true": float(
                                    radar_time_true[target_index]
                                ),
                                "radar_distance_true": float(
                                    radar_distance_true[target_index]
                                ),
                                "radar_time_error": time_hat
                                - radar_time_true[target_index]
                                if reconstructed_accessible
                                else float("nan"),
                                "radar_distance_error": distance_hat
                                - radar_distance_true[target_index]
                                if reconstructed_accessible
                                else float("nan"),
                            }
                        )

    return event_rows, summarize_results(event_rows)


def _rows_to_array(
    rows: list[dict[str, float]],
    key: str,
    dtype: type = float,
) -> np.ndarray:
    return np.asarray([row[key] for row in rows], dtype=dtype)


def summarize_results(event_rows: list[dict[str, float]]) -> list[dict[str, float]]:
    """Aggregate Rindler event-level results."""

    summary_rows: list[dict[str, float]] = []
    keys = sorted(
        {
            (
                row["acceleration"],
                int(row["N"]),
                int(row["tick_count"]),
            )
            for row in event_rows
        }
    )
    for acceleration, n_events, tick_count in keys:
        rows = [
            row
            for row in event_rows
            if row["acceleration"] == acceleration
            and int(row["N"]) == n_events
            and int(row["tick_count"]) == tick_count
        ]
        wedge = _rows_to_array(rows, "analytic_wedge_accessible", bool)
        finite = _rows_to_array(rows, "finite_coverage_accessible", bool)
        reconstructed = _rows_to_array(rows, "reconstructed_accessible", bool)
        confusion = horizon_accessibility_confusion(reconstructed, finite)
        errors = radar_error_summary(
            _rows_to_array(rows, "radar_time_hat"),
            _rows_to_array(rows, "radar_distance_hat"),
            _rows_to_array(rows, "radar_time_true"),
            _rows_to_array(rows, "radar_distance_true"),
            finite & reconstructed,
        )
        summary_rows.append(
            {
                "acceleration": float(acceleration),
                "N": float(n_events),
                "tick_count": float(tick_count),
                "repetitions": float(
                    len({int(row["repetition"]) for row in rows})
                ),
                "event_count": float(len(rows)),
                "wedge_accessible_fraction": float(np.mean(wedge)),
                "finite_coverage_accessible_fraction": float(np.mean(finite)),
                "reconstructed_accessible_fraction": float(np.mean(reconstructed)),
                "true_positive": confusion["true_positive"],
                "false_positive": confusion["false_positive"],
                "true_negative": confusion["true_negative"],
                "false_negative": confusion["false_negative"],
                "precision": confusion["precision"],
                "recall": confusion["recall"],
                "false_positive_rate": confusion["false_positive_rate"],
                "false_negative_rate": confusion["false_negative_rate"],
                "radar_time_rmse_finite_coverage": errors["radar_time_rmse"],
                "radar_distance_rmse_finite_coverage": errors[
                    "radar_distance_rmse"
                ],
                "radar_time_bias_finite_coverage": errors["radar_time_bias"],
                "radar_distance_bias_finite_coverage": errors[
                    "radar_distance_bias"
                ],
            }
        )
    return summary_rows


def _write_csv(rows: list[dict[str, float]], output_path: Path) -> Path:
    if not rows:
        raise RuntimeError("no rows to write")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def write_outputs(
    event_rows: list[dict[str, float]],
    summary_rows: list[dict[str, float]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write Rindler horizon CSV outputs."""

    data_dir = output_dir / "data"
    event_path = data_dir / "rindler_horizon_reconstruction_events.csv"
    summary_path = data_dir / "rindler_horizon_reconstruction_summary.csv"
    return _write_csv(event_rows, event_path), _write_csv(summary_rows, summary_path)


def _plot_identity(ax: plt.Axes, truth: np.ndarray, estimate: np.ndarray) -> None:
    if truth.size == 0:
        return
    lower = float(np.nanmin([np.nanmin(truth), np.nanmin(estimate)]))
    upper = float(np.nanmax([np.nanmax(truth), np.nanmax(estimate)]))
    ax.plot([lower, upper], [lower, upper], color="black", linewidth=1.0)


def save_accessibility_map(
    event_rows: list[dict[str, float]],
    config: ExperimentConfig,
) -> Path:
    """Save Rindler accessibility scatter plot."""

    output_path = config.output_dir / "figures" / "rindler_accessibility_map.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    acceleration = config.acceleration_values[0]
    n_events = max(config.n_values)
    tick_count = max(config.tick_values)
    rows = [
        row
        for row in event_rows
        if row["acceleration"] == acceleration
        and int(row["N"]) == n_events
        and int(row["tick_count"]) == tick_count
        and int(row["repetition"]) == 0
    ]
    t = _rows_to_array(rows, "t")
    x = _rows_to_array(rows, "x")
    reconstructed = _rows_to_array(rows, "reconstructed_accessible", bool)
    finite = _rows_to_array(rows, "finite_coverage_accessible", bool)
    tau_min, tau_max = safe_tau_range_for_rindler_chain_1p1(
        config.T,
        acceleration,
        direction=config.direction,
    )
    chain_events, _ = make_rindler_observer_chain_1p1(
        acceleration,
        200,
        tau_min,
        tau_max,
        direction=config.direction,
    )

    fig, ax = plt.subplots(figsize=(6.8, 6.0))
    colors = np.where(reconstructed, "tab:green", np.where(finite, "tab:orange", "0.7"))
    ax.scatter(t, x, c=colors, s=9, alpha=0.65)
    line = np.linspace(-config.T / 2.0, config.T / 2.0, 200)
    ax.plot(line, config.direction * line, color="black", linewidth=0.9)
    ax.plot(line, -config.direction * line, color="black", linewidth=0.9)
    ax.plot(chain_events[:, 0], chain_events[:, 1], color="tab:red", linewidth=1.4)
    ax.set_xlabel("t")
    ax.set_ylabel("x")
    ax.set_title("Rindler two-way radar accessibility")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_accessible_fraction_vs_ticks(
    summary_rows: list[dict[str, float]],
    output_dir: Path,
) -> Path:
    """Save wedge, finite-coverage, and reconstructed fractions."""

    output_path = (
        output_dir / "figures" / "rindler_accessible_fraction_vs_ticks.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ticks = sorted({int(row["tick_count"]) for row in summary_rows})
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    for key, label in [
        ("wedge_accessible_fraction", "infinite wedge"),
        ("finite_coverage_accessible_fraction", "finite chain coverage"),
        ("reconstructed_accessible_fraction", "reconstructed"),
    ]:
        values = [
            float(
                np.nanmean(
                    [
                        row[key]
                        for row in summary_rows
                        if int(row["tick_count"]) == tick
                    ]
                )
            )
            for tick in ticks
        ]
        ax.plot(ticks, values, marker="o", label=label)
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel("Accessible fraction")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("Rindler accessibility notions")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def _finite_reconstructed_arrays(
    rows: list[dict[str, float]],
    true_key: str,
    estimate_key: str,
) -> tuple[np.ndarray, np.ndarray]:
    selected = [
        row
        for row in rows
        if bool(row["finite_coverage_accessible"])
        and bool(row["reconstructed_accessible"])
        and np.isfinite(row[true_key])
        and np.isfinite(row[estimate_key])
    ]
    return (
        np.asarray([row[true_key] for row in selected], dtype=float),
        np.asarray([row[estimate_key] for row in selected], dtype=float),
    )


def save_radar_time_scatter(
    event_rows: list[dict[str, float]],
    output_dir: Path,
) -> Path:
    """Save analytic versus reconstructed Rindler radar time."""

    output_path = output_dir / "figures" / "rindler_radar_time_scatter.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    truth, estimate = _finite_reconstructed_arrays(
        event_rows,
        "radar_time_true",
        "radar_time_hat",
    )
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    ax.scatter(truth, estimate, s=7, alpha=0.35)
    _plot_identity(ax, truth, estimate)
    ax.set_xlabel("Analytic Rindler radar time")
    ax.set_ylabel("Reconstructed radar time")
    ax.set_title("Rindler radar time reconstruction")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_radar_distance_scatter(
    event_rows: list[dict[str, float]],
    output_dir: Path,
) -> Path:
    """Save analytic versus reconstructed Rindler radar distance."""

    output_path = output_dir / "figures" / "rindler_radar_distance_scatter.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    truth, estimate = _finite_reconstructed_arrays(
        event_rows,
        "radar_distance_true",
        "radar_distance_hat",
    )
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    ax.scatter(truth, estimate, s=7, alpha=0.35)
    _plot_identity(ax, truth, estimate)
    ax.set_xlabel("Analytic Rindler radar distance")
    ax.set_ylabel("Reconstructed radar distance")
    ax.set_title("Rindler radar distance reconstruction")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_error_vs_ticks(
    summary_rows: list[dict[str, float]],
    output_dir: Path,
) -> Path:
    """Save radar-coordinate RMSE versus tick count."""

    output_path = output_dir / "figures" / "rindler_error_vs_ticks.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ticks = sorted({int(row["tick_count"]) for row in summary_rows})
    time_values = []
    distance_values = []
    for tick in ticks:
        rows = [row for row in summary_rows if int(row["tick_count"]) == tick]
        time_values.append(
            float(np.nanmean([row["radar_time_rmse_finite_coverage"] for row in rows]))
        )
        distance_values.append(
            float(
                np.nanmean(
                    [row["radar_distance_rmse_finite_coverage"] for row in rows]
                )
            )
        )
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    ax.plot(ticks, time_values, marker="o", label="radar time RMSE")
    ax.plot(ticks, distance_values, marker="s", label="radar distance RMSE")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel("RMSE")
    ax.set_title("Rindler radar reconstruction error")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_false_rates_vs_ticks(
    summary_rows: list[dict[str, float]],
    output_dir: Path,
) -> Path:
    """Save false-positive and false-negative rates versus tick count."""

    output_path = (
        output_dir / "figures" / "rindler_false_positive_negative_vs_ticks.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ticks = sorted({int(row["tick_count"]) for row in summary_rows})
    fpr = []
    fnr = []
    for tick in ticks:
        rows = [row for row in summary_rows if int(row["tick_count"]) == tick]
        fpr.append(float(np.nanmean([row["false_positive_rate"] for row in rows])))
        fnr.append(float(np.nanmean([row["false_negative_rate"] for row in rows])))
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    ax.plot(ticks, fpr, marker="o", label="false positive rate")
    ax.plot(ticks, fnr, marker="s", label="false negative rate")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel("Rate relative to finite coverage")
    ax.set_title("Rindler accessibility confusion rates")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_figures(
    event_rows: list[dict[str, float]],
    summary_rows: list[dict[str, float]],
    config: ExperimentConfig,
) -> tuple[Path, Path, Path, Path, Path, Path]:
    """Save all Rindler horizon figures."""

    return (
        save_accessibility_map(event_rows, config),
        save_accessible_fraction_vs_ticks(summary_rows, config.output_dir),
        save_radar_time_scatter(event_rows, config.output_dir),
        save_radar_distance_scatter(event_rows, config.output_dir),
        save_error_vs_ticks(summary_rows, config.output_dir),
        save_false_rates_vs_ticks(summary_rows, config.output_dir),
    )


def main() -> None:
    config = parse_args()
    event_rows, summary_rows = run_experiment(config)
    event_path, summary_path = write_outputs(
        event_rows,
        summary_rows,
        config.output_dir,
    )
    figure_paths = save_figures(event_rows, summary_rows, config)
    print(f"Wrote event results: {event_path}")
    print(f"Wrote summary: {summary_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")
    print(
        "Rindler accessibility is a controlled flat-spacetime horizon analogue; "
        "finite-chain coverage is reported separately from wedge accessibility."
    )


if __name__ == "__main__":
    main()
