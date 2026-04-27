"""Oriented two-chain radar reconstruction and Lorentz-map recovery."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.discrete_radar import discrete_radar_coordinates_from_order
from causal_spacetime_lab.estimators import rmse
from causal_spacetime_lab.lorentz_maps import (
    fit_lorentz_beta_grid,
    lorentz_residual_rmse,
    true_lorentz_rest_coordinates_1p1,
)
from causal_spacetime_lab.observer import (
    make_inertial_observer_chain_1p1,
    observer_chain_indices,
    safe_tau_range_for_inertial_chain_1p1,
)
from causal_spacetime_lab.oriented_radar import (
    oriented_radar_coordinates_from_two_chains,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_N_VALUES = (300, 600, 1200)
DEFAULT_TICK_VALUES = (32, 64, 128)
DEFAULT_BETA_VALUES = (0.3, 0.6)
DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for oriented radar Lorentz-map recovery."""

    T: float = 2.0
    n_values: tuple[int, ...] = DEFAULT_N_VALUES
    tick_values: tuple[int, ...] = DEFAULT_TICK_VALUES
    beta_values: tuple[float, ...] = DEFAULT_BETA_VALUES
    beacon_separation: float = 0.15
    repetitions: int = 5
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


@dataclass(frozen=True)
class ChainIndices:
    """Observer-chain index block and clocks in the combined event array."""

    primary: np.ndarray
    beacon: np.ndarray
    clocks: np.ndarray


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
            "Recover oriented radar coordinates and Lorentz maps from two-chain "
            "observer protocols in known 1+1D Minkowski intervals."
        )
    )
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument(
        "--n-values",
        nargs="+",
        default=[str(n) for n in DEFAULT_N_VALUES],
    )
    parser.add_argument(
        "--tick-values",
        nargs="+",
        default=[str(tick) for tick in DEFAULT_TICK_VALUES],
    )
    parser.add_argument(
        "--beta-values",
        nargs="+",
        default=[str(beta) for beta in DEFAULT_BETA_VALUES],
    )
    parser.add_argument("--beacon-separation", type=float, default=0.15)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    return ExperimentConfig(
        T=args.T,
        n_values=_parse_int_values(args.n_values, "N"),
        tick_values=_parse_int_values(args.tick_values, "tick"),
        beta_values=_parse_float_values(args.beta_values, "beta"),
        beacon_separation=args.beacon_separation,
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _validate_config(config: ExperimentConfig) -> None:
    if config.T <= 0:
        raise ValueError("T must be positive")
    if any(n <= 0 for n in config.n_values):
        raise ValueError("all n_values must be positive")
    if any(tick < 2 for tick in config.tick_values):
        raise ValueError("all tick_values must be at least 2")
    if any(abs(beta) >= 1.0 for beta in config.beta_values):
        raise ValueError("all beta_values must satisfy abs(beta) < 1")
    if config.beacon_separation <= 0.0:
        raise ValueError("beacon_separation must be positive")
    if config.repetitions <= 0:
        raise ValueError("repetitions must be positive")


def _common_tau_range(
    T: float,
    beta: float,
    beacon_separation: float,
) -> tuple[float, float]:
    primary = safe_tau_range_for_inertial_chain_1p1(T, beta, x_prime=0.0)
    beacon = safe_tau_range_for_inertial_chain_1p1(
        T,
        beta,
        x_prime=beacon_separation,
    )
    tau_min = max(primary[0], beacon[0])
    tau_max = min(primary[1], beacon[1])
    if tau_max <= tau_min:
        raise ValueError("primary and beacon chains have no common safe tau range")
    return tau_min, tau_max


def _append_protocol_chains(
    chain_events: list[np.ndarray],
    beta: float,
    tick_count: int,
    tau_range: tuple[float, float],
    beacon_separation: float,
    start_index: int,
) -> tuple[ChainIndices, int]:
    primary_events, clocks = make_inertial_observer_chain_1p1(
        beta,
        tick_count,
        tau_range[0],
        tau_range[1],
        x_prime=0.0,
    )
    beacon_events, beacon_clocks = make_inertial_observer_chain_1p1(
        beta,
        tick_count,
        tau_range[0],
        tau_range[1],
        x_prime=beacon_separation,
    )
    if not np.allclose(clocks, beacon_clocks):
        raise RuntimeError("primary and beacon clocks are not synchronized")

    primary_indices = observer_chain_indices(start_index, tick_count)
    beacon_indices = observer_chain_indices(start_index + tick_count, tick_count)
    chain_events.extend([primary_events, beacon_events])
    return (
        ChainIndices(primary=primary_indices, beacon=beacon_indices, clocks=clocks),
        start_index + 2 * tick_count,
    )


def _oriented_coords_to_array(oriented_results) -> np.ndarray:  # type: ignore[no-untyped-def]
    return np.asarray(
        [
            [
                result.oriented_time if result.oriented_time is not None else np.nan,
                (
                    result.signed_position
                    if result.signed_position is not None
                    else np.nan
                ),
            ]
            for result in oriented_results
        ],
        dtype=float,
    )


def _fit_run_lorentz(
    lab_hat: np.ndarray,
    moving_hat: np.ndarray,
    beta: float,
) -> tuple[float, float, float]:
    if lab_hat.shape[0] == 0:
        return float("nan"), float("nan"), float("nan")
    known_residual = lorentz_residual_rmse(lab_hat, moving_hat, beta)
    fitted_beta, fitted_residual = fit_lorentz_beta_grid(lab_hat, moving_hat)
    return known_residual, fitted_beta, fitted_residual


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float]], list[dict[str, float]]]:
    """Run oriented radar reconstruction and return event and summary rows."""

    _validate_config(config)
    event_rows: list[dict[str, float]] = []
    run_rows: list[dict[str, float]] = []

    for n_index, n_events in enumerate(config.n_values):
        for beta in config.beta_values:
            for repetition in range(config.repetitions):
                support_events = sprinkle_1p1_causal_diamond(
                    n_events,
                    T=config.T,
                    seed=config.seed + 100_000 * n_index + 1_000 * repetition,
                )
                true_moving = true_lorentz_rest_coordinates_1p1(support_events, beta)
                for tick_count in config.tick_values:
                    stationary_range = _common_tau_range(
                        config.T,
                        0.0,
                        config.beacon_separation,
                    )
                    moving_range = _common_tau_range(
                        config.T,
                        beta,
                        config.beacon_separation,
                    )
                    chain_events: list[np.ndarray] = []
                    next_index = n_events
                    stationary_indices, next_index = _append_protocol_chains(
                        chain_events,
                        beta=0.0,
                        tick_count=tick_count,
                        tau_range=stationary_range,
                        beacon_separation=config.beacon_separation,
                        start_index=next_index,
                    )
                    moving_indices, next_index = _append_protocol_chains(
                        chain_events,
                        beta=beta,
                        tick_count=tick_count,
                        tau_range=moving_range,
                        beacon_separation=config.beacon_separation,
                        start_index=next_index,
                    )
                    combined_events = np.vstack((support_events, *chain_events))
                    causal_matrix = causal_matrix_1p1(combined_events)
                    target_indices = np.arange(n_events, dtype=int)

                    stationary_primary = discrete_radar_coordinates_from_order(
                        causal_matrix,
                        stationary_indices.primary,
                        target_indices,
                        stationary_indices.clocks,
                    )
                    stationary_beacon = discrete_radar_coordinates_from_order(
                        causal_matrix,
                        stationary_indices.beacon,
                        target_indices,
                        stationary_indices.clocks,
                    )
                    moving_primary = discrete_radar_coordinates_from_order(
                        causal_matrix,
                        moving_indices.primary,
                        target_indices,
                        moving_indices.clocks,
                    )
                    moving_beacon = discrete_radar_coordinates_from_order(
                        causal_matrix,
                        moving_indices.beacon,
                        target_indices,
                        moving_indices.clocks,
                    )
                    stationary_oriented = oriented_radar_coordinates_from_two_chains(
                        stationary_primary,
                        stationary_beacon,
                        config.beacon_separation,
                    )
                    moving_oriented = oriented_radar_coordinates_from_two_chains(
                        moving_primary,
                        moving_beacon,
                        config.beacon_separation,
                    )

                    lab_hat_all = _oriented_coords_to_array(stationary_oriented)
                    moving_hat_all = _oriented_coords_to_array(moving_oriented)
                    accessible = np.asarray(
                        [
                            lab.accessible and moving.accessible
                            for lab, moving in zip(
                                stationary_oriented,
                                moving_oriented,
                                strict=True,
                            )
                        ],
                        dtype=bool,
                    )
                    lab_hat = lab_hat_all[accessible]
                    moving_hat = moving_hat_all[accessible]
                    known_residual, fitted_beta, fitted_residual = _fit_run_lorentz(
                        lab_hat,
                        moving_hat,
                        beta,
                    )
                    run_rows.append(
                        {
                            "N": float(n_events),
                            "tick_count": float(tick_count),
                            "beta": float(beta),
                            "repetition": float(repetition),
                            "event_count": float(n_events),
                            "accessible_count": float(np.count_nonzero(accessible)),
                            "known_beta_lorentz_residual_rmse": known_residual,
                            "fitted_beta": fitted_beta,
                            "fitted_beta_error": fitted_beta - beta
                            if np.isfinite(fitted_beta)
                            else float("nan"),
                            "fitted_lorentz_residual_rmse": fitted_residual,
                        }
                    )

                    for target_index in range(n_events):
                        is_accessible = bool(accessible[target_index])
                        lab_hat_row = lab_hat_all[target_index]
                        moving_hat_row = moving_hat_all[target_index]
                        true_lab = support_events[target_index]
                        true_moving_row = true_moving[target_index]
                        event_rows.append(
                            {
                                "N": float(n_events),
                                "tick_count": float(tick_count),
                                "beta": float(beta),
                                "repetition": float(repetition),
                                "target_index": float(target_index),
                                "accessible": float(is_accessible),
                                "lab_time_hat": float(lab_hat_row[0]),
                                "lab_position_hat": float(lab_hat_row[1]),
                                "moving_time_hat": float(moving_hat_row[0]),
                                "moving_position_hat": float(moving_hat_row[1]),
                                "true_lab_time": float(true_lab[0]),
                                "true_lab_position": float(true_lab[1]),
                                "true_moving_time": float(true_moving_row[0]),
                                "true_moving_position": float(true_moving_row[1]),
                                "lab_time_error": float(lab_hat_row[0] - true_lab[0])
                                if is_accessible
                                else float("nan"),
                                "lab_position_error": float(
                                    lab_hat_row[1] - true_lab[1]
                                )
                                if is_accessible
                                else float("nan"),
                                "moving_time_error": float(
                                    moving_hat_row[0] - true_moving_row[0]
                                )
                                if is_accessible
                                else float("nan"),
                                "moving_position_error": float(
                                    moving_hat_row[1] - true_moving_row[1]
                                )
                                if is_accessible
                                else float("nan"),
                                "known_beta_lorentz_residual_rmse": known_residual,
                                "fitted_beta": fitted_beta,
                                "fitted_beta_error": fitted_beta - beta
                                if np.isfinite(fitted_beta)
                                else float("nan"),
                                "fitted_lorentz_residual_rmse": fitted_residual,
                            }
                        )

    return event_rows, summarize_results(event_rows, run_rows)


def _rmse_for_rows(rows: list[dict[str, float]], truth: str, estimate: str) -> float:
    accessible = [row for row in rows if bool(row["accessible"])]
    if not accessible:
        return float("nan")
    return rmse(
        [row[truth] for row in accessible],
        [row[estimate] for row in accessible],
    )


def summarize_results(
    event_rows: list[dict[str, float]],
    run_rows: list[dict[str, float]],
) -> list[dict[str, float]]:
    """Aggregate event-level and run-level oriented radar results."""

    summary_rows: list[dict[str, float]] = []
    keys = sorted(
        {
            (int(row["N"]), int(row["tick_count"]), float(row["beta"]))
            for row in event_rows
        }
    )
    for n_events, tick_count, beta in keys:
        events = [
            row
            for row in event_rows
            if int(row["N"]) == n_events
            and int(row["tick_count"]) == tick_count
            and row["beta"] == beta
        ]
        runs = [
            row
            for row in run_rows
            if int(row["N"]) == n_events
            and int(row["tick_count"]) == tick_count
            and row["beta"] == beta
        ]
        accessible_count = sum(1 for row in events if bool(row["accessible"]))
        fitted_betas = np.asarray([row["fitted_beta"] for row in runs], dtype=float)
        fitted_errors = np.asarray(
            [row["fitted_beta_error"] for row in runs],
            dtype=float,
        )
        finite_fit = np.isfinite(fitted_betas)
        summary_rows.append(
            {
                "N": float(n_events),
                "tick_count": float(tick_count),
                "beta": float(beta),
                "repetitions": float(len(runs)),
                "event_count": float(len(events)),
                "accessible_count": float(accessible_count),
                "accessible_fraction": accessible_count / len(events),
                "lab_time_rmse": _rmse_for_rows(
                    events,
                    "true_lab_time",
                    "lab_time_hat",
                ),
                "lab_position_rmse": _rmse_for_rows(
                    events,
                    "true_lab_position",
                    "lab_position_hat",
                ),
                "moving_time_rmse": _rmse_for_rows(
                    events,
                    "true_moving_time",
                    "moving_time_hat",
                ),
                "moving_position_rmse": _rmse_for_rows(
                    events,
                    "true_moving_position",
                    "moving_position_hat",
                ),
                "known_beta_lorentz_residual_rmse": float(
                    np.nanmean(
                        [row["known_beta_lorentz_residual_rmse"] for row in runs]
                    )
                ),
                "fitted_beta_mean": float(np.mean(fitted_betas[finite_fit]))
                if np.any(finite_fit)
                else float("nan"),
                "fitted_beta_std": float(np.std(fitted_betas[finite_fit], ddof=1))
                if np.count_nonzero(finite_fit) > 1
                else 0.0,
                "fitted_beta_bias": float(np.mean(fitted_errors[finite_fit]))
                if np.any(finite_fit)
                else float("nan"),
                "fitted_beta_rmse": float(
                    np.sqrt(np.mean(fitted_errors[finite_fit] ** 2))
                )
                if np.any(finite_fit)
                else float("nan"),
                "fitted_lorentz_residual_rmse": float(
                    np.nanmean([row["fitted_lorentz_residual_rmse"] for row in runs])
                ),
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
    """Write oriented radar Lorentz-map CSV outputs."""

    data_dir = output_dir / "data"
    event_path = data_dir / "oriented_radar_lorentz_events.csv"
    summary_path = data_dir / "oriented_radar_lorentz_summary.csv"
    return _write_csv(event_rows, event_path), _write_csv(summary_rows, summary_path)


def _accessible_arrays(
    rows: list[dict[str, float]],
    truth_key: str,
    estimate_key: str,
) -> tuple[np.ndarray, np.ndarray]:
    accessible = [row for row in rows if bool(row["accessible"])]
    truth = np.asarray([row[truth_key] for row in accessible], dtype=float)
    estimate = np.asarray([row[estimate_key] for row in accessible], dtype=float)
    return truth, estimate


def _plot_identity(ax: plt.Axes, truth: np.ndarray, estimate: np.ndarray) -> None:
    if truth.size == 0:
        return
    lower = float(min(np.min(truth), np.min(estimate)))
    upper = float(max(np.max(truth), np.max(estimate)))
    ax.plot([lower, upper], [lower, upper], color="black", linewidth=1.0)


def save_lab_position_scatter(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save true lab x versus oriented reconstructed lab x."""

    output_path = output_dir / "figures" / "oriented_radar_lab_position_scatter.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    truth, estimate = _accessible_arrays(rows, "true_lab_position", "lab_position_hat")
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    ax.scatter(truth, estimate, s=8, alpha=0.4)
    _plot_identity(ax, truth, estimate)
    ax.set_xlabel("True lab x")
    ax.set_ylabel("Reconstructed oriented lab x")
    ax.set_title("Two-chain oriented lab coordinate")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_moving_position_scatter(
    rows: list[dict[str, float]],
    output_dir: Path,
) -> Path:
    """Save true moving x' versus oriented reconstructed moving x'."""

    output_path = output_dir / "figures" / "oriented_radar_moving_position_scatter.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    truth, estimate = _accessible_arrays(
        rows,
        "true_moving_position",
        "moving_position_hat",
    )
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    ax.scatter(truth, estimate, s=8, alpha=0.4)
    _plot_identity(ax, truth, estimate)
    ax.set_xlabel("True moving-frame x'")
    ax.set_ylabel("Reconstructed oriented moving-frame x'")
    ax.set_title("Two-chain oriented moving coordinate")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_lorentz_residual_vs_ticks(
    summary_rows: list[dict[str, float]],
    output_dir: Path,
) -> Path:
    """Save Lorentz residual RMSE versus tick count."""

    output_path = (
        output_dir / "figures" / "oriented_radar_lorentz_residual_vs_ticks.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for beta in sorted({row["beta"] for row in summary_rows}):
        for n_events in sorted({int(row["N"]) for row in summary_rows}):
            rows = [
                row
                for row in summary_rows
                if row["beta"] == beta and int(row["N"]) == n_events
            ]
            ticks = np.asarray([row["tick_count"] for row in rows])
            residual = np.asarray(
                [row["known_beta_lorentz_residual_rmse"] for row in rows]
            )
            ax.plot(ticks, residual, marker="o", label=f"beta={beta}, N={n_events}")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel("Known-beta Lorentz residual RMSE")
    ax.set_title("Lorentz-map residual decreases with orientation resolution")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small", ncol=2)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_beta_fit_vs_ticks(
    summary_rows: list[dict[str, float]],
    output_dir: Path,
) -> Path:
    """Save fitted beta versus tick count."""

    output_path = output_dir / "figures" / "oriented_radar_beta_fit_vs_ticks.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for beta in sorted({row["beta"] for row in summary_rows}):
        ax.axhline(beta, color="0.7", linewidth=0.9)
        for n_events in sorted({int(row["N"]) for row in summary_rows}):
            rows = [
                row
                for row in summary_rows
                if row["beta"] == beta and int(row["N"]) == n_events
            ]
            ticks = np.asarray([row["tick_count"] for row in rows])
            fitted = np.asarray([row["fitted_beta_mean"] for row in rows])
            ax.plot(ticks, fitted, marker="o", label=f"beta={beta}, N={n_events}")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel("Fitted beta")
    ax.set_title("Lorentz beta fit from reconstructed coordinates")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small", ncol=2)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_accessible_fraction(
    summary_rows: list[dict[str, float]],
    output_dir: Path,
) -> Path:
    """Save accessible fraction versus tick count."""

    output_path = output_dir / "figures" / "oriented_radar_accessible_fraction.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    for beta in sorted({row["beta"] for row in summary_rows}):
        rows = [row for row in summary_rows if row["beta"] == beta]
        ticks = np.asarray([row["tick_count"] for row in rows])
        fractions = np.asarray([row["accessible_fraction"] for row in rows])
        ax.scatter(ticks, fractions, label=f"beta={beta}")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel("Accessible fraction")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("Four-chain accessibility")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_figures(
    event_rows: list[dict[str, float]],
    summary_rows: list[dict[str, float]],
    output_dir: Path,
) -> tuple[Path, Path, Path, Path, Path]:
    """Save all oriented radar Lorentz-map figures."""

    return (
        save_lab_position_scatter(event_rows, output_dir),
        save_moving_position_scatter(event_rows, output_dir),
        save_lorentz_residual_vs_ticks(summary_rows, output_dir),
        save_beta_fit_vs_ticks(summary_rows, output_dir),
        save_accessible_fraction(summary_rows, output_dir),
    )


def main() -> None:
    config = parse_args()
    event_rows, summary_rows = run_experiment(config)
    event_path, summary_path = write_outputs(
        event_rows,
        summary_rows,
        config.output_dir,
    )
    figure_paths = save_figures(event_rows, summary_rows, config.output_dir)
    print(f"Wrote event results: {event_path}")
    print(f"Wrote summary: {summary_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")
    print(
        "Signed coordinates require the supplied beacon separation; hidden "
        "coordinates are used only for validation."
    )


if __name__ == "__main__":
    main()
