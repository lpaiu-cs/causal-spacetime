"""Discrete observer-chain radar reconstruction from causal order."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.chains import longest_chain_indices
from causal_spacetime_lab.discrete_radar import discrete_radar_coordinates_from_order
from causal_spacetime_lab.estimators import bias, mae, rmse
from causal_spacetime_lab.observer import (
    make_stationary_observer_chain_1p1,
    observer_chain_indices,
)
from causal_spacetime_lab.radar_validation import (
    accessible_fraction,
    true_stationary_radar_coordinates_1p1,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_N_VALUES = (300, 600, 1200, 2400)
DEFAULT_TICK_VALUES = (16, 32, 64, 128)
DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for discrete observer-chain radar reconstruction."""

    T: float = 2.0
    n_values: tuple[int, ...] = DEFAULT_N_VALUES
    tick_values: tuple[int, ...] = DEFAULT_TICK_VALUES
    repetitions: int = 5
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR
    use_extracted_chain: bool = False


@dataclass(frozen=True)
class ProtocolData:
    """Combined-event protocol data for one radar reconstruction run."""

    combined_events: np.ndarray
    observer_indices: np.ndarray
    clock_times: np.ndarray
    target_indices: np.ndarray
    tick_count: int


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
        description=(
            "Validate causal-order-based radar reconstruction from a supplied "
            "observer chain and clock labels."
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
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--use-extracted-chain", action="store_true")
    args = parser.parse_args()

    return ExperimentConfig(
        T=args.T,
        n_values=_parse_int_values(args.n_values, "N"),
        tick_values=_parse_int_values(args.tick_values, "tick"),
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
        use_extracted_chain=args.use_extracted_chain,
    )


def _validate_config(config: ExperimentConfig) -> None:
    if config.T <= 0:
        raise ValueError("T must be positive")
    if any(n <= 0 for n in config.n_values):
        raise ValueError("all n_values must be positive")
    if any(ticks < 2 for ticks in config.tick_values):
        raise ValueError("all tick_values must be at least 2")
    if config.repetitions <= 0:
        raise ValueError("repetitions must be positive")


def extract_longest_observer_chain_between(
    causal_matrix: np.ndarray,
    start_index: int,
    end_index: int,
) -> np.ndarray:
    """Return one longest chain between endpoints as exploratory observer chain."""

    chain = longest_chain_indices(causal_matrix, start_index, end_index)
    if chain.size < 2:
        raise RuntimeError("could not extract an observer chain between endpoints")
    return chain


def _supplied_observer_protocol(
    support_events: np.ndarray,
    T: float,
    tick_count: int,
) -> ProtocolData:
    observer_events, clock_times = make_stationary_observer_chain_1p1(
        T,
        tick_count,
    )
    combined_events = np.vstack((support_events, observer_events))
    observer_indices = observer_chain_indices(support_events.shape[0], tick_count)
    target_indices = np.arange(support_events.shape[0], dtype=int)
    return ProtocolData(
        combined_events=combined_events,
        observer_indices=observer_indices,
        clock_times=clock_times,
        target_indices=target_indices,
        tick_count=tick_count,
    )


def _extracted_observer_protocol(
    support_events: np.ndarray,
    T: float,
) -> ProtocolData:
    p = np.array([-T / 2.0, 0.0])
    q = np.array([T / 2.0, 0.0])
    combined_events = np.vstack((support_events, p, q))
    start_index = support_events.shape[0]
    end_index = support_events.shape[0] + 1
    causal_matrix = causal_matrix_1p1(combined_events)
    observer_indices = extract_longest_observer_chain_between(
        causal_matrix,
        start_index,
        end_index,
    )
    clock_times = np.linspace(-T / 2.0, T / 2.0, observer_indices.size)
    support_indices = np.arange(support_events.shape[0], dtype=int)
    target_indices = np.setdiff1d(
        support_indices,
        observer_indices,
        assume_unique=False,
    )
    return ProtocolData(
        combined_events=combined_events,
        observer_indices=observer_indices,
        clock_times=clock_times,
        target_indices=target_indices,
        tick_count=int(observer_indices.size),
    )


def _radar_value(value: float | None) -> float:
    return float("nan") if value is None else float(value)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float]], list[dict[str, float]]]:
    """Run discrete radar reconstruction and return event and summary rows."""

    _validate_config(config)
    event_rows: list[dict[str, float]] = []

    for n_index, n_events in enumerate(config.n_values):
        for repetition in range(config.repetitions):
            support_events = sprinkle_1p1_causal_diamond(
                n_events,
                T=config.T,
                seed=config.seed + 10_000 * n_index + repetition,
            )
            tick_values = (0,) if config.use_extracted_chain else config.tick_values
            for requested_tick_count in tick_values:
                if config.use_extracted_chain:
                    protocol = _extracted_observer_protocol(support_events, config.T)
                    requested_tick = float("nan")
                else:
                    protocol = _supplied_observer_protocol(
                        support_events,
                        config.T,
                        requested_tick_count,
                    )
                    requested_tick = float(requested_tick_count)

                causal_matrix = causal_matrix_1p1(protocol.combined_events)
                reconstructions = discrete_radar_coordinates_from_order(
                    causal_matrix,
                    protocol.observer_indices,
                    protocol.target_indices,
                    protocol.clock_times,
                )
                true_time, true_distance = true_stationary_radar_coordinates_1p1(
                    protocol.combined_events[protocol.target_indices]
                )

                for local_number, (target_index, result) in enumerate(
                    zip(protocol.target_indices, reconstructions, strict=True)
                ):
                    time_hat = _radar_value(result.radar_time)
                    distance_hat = _radar_value(result.radar_distance)
                    time_error = (
                        time_hat - true_time[local_number]
                        if result.accessible
                        else float("nan")
                    )
                    distance_error = (
                        distance_hat - true_distance[local_number]
                        if result.accessible
                        else float("nan")
                    )
                    event_rows.append(
                        {
                            "N": float(n_events),
                            "tick_count": float(protocol.tick_count),
                            "requested_tick_count": requested_tick,
                            "repetition": float(repetition),
                            "target_index": float(target_index),
                            "accessible": float(result.accessible),
                            "tau_minus": _radar_value(result.tau_minus),
                            "tau_plus": _radar_value(result.tau_plus),
                            "radar_time_hat": time_hat,
                            "radar_distance_hat": distance_hat,
                            "true_radar_time": float(true_time[local_number]),
                            "true_radar_distance": float(true_distance[local_number]),
                            "radar_time_error": float(time_error),
                            "radar_distance_error": float(distance_error),
                            "absolute_radar_time_error": abs(float(time_error))
                            if result.accessible
                            else float("nan"),
                            "absolute_radar_distance_error": abs(
                                float(distance_error)
                            )
                            if result.accessible
                            else float("nan"),
                            "use_extracted_chain": float(config.use_extracted_chain),
                        }
                    )

    summary_rows = summarize_results(event_rows)
    return event_rows, summary_rows


def _summary_metrics(
    rows: list[dict[str, float]],
    estimate_key: str,
    truth_key: str,
) -> tuple[float, float, float]:
    accessible_rows = [row for row in rows if bool(row["accessible"])]
    if not accessible_rows:
        return float("nan"), float("nan"), float("nan")
    truth = np.asarray([row[truth_key] for row in accessible_rows], dtype=float)
    estimate = np.asarray([row[estimate_key] for row in accessible_rows], dtype=float)
    return mae(truth, estimate), rmse(truth, estimate), bias(truth, estimate)


def summarize_results(rows: list[dict[str, float]]) -> list[dict[str, float]]:
    """Aggregate radar reconstruction rows by event count and tick count."""

    summary_rows: list[dict[str, float]] = []
    keys = sorted({(int(row["N"]), int(row["tick_count"])) for row in rows})
    for n_events, tick_count in keys:
        group = [
            row
            for row in rows
            if int(row["N"]) == n_events and int(row["tick_count"]) == tick_count
        ]
        accessible_values = np.asarray([row["accessible"] for row in group], dtype=bool)
        time_metrics = _summary_metrics(group, "radar_time_hat", "true_radar_time")
        distance_metrics = _summary_metrics(
            group,
            "radar_distance_hat",
            "true_radar_distance",
        )
        summary_rows.append(
            {
                "N": float(n_events),
                "tick_count": float(tick_count),
                "repetitions": float(len({int(row["repetition"]) for row in group})),
                "event_count": float(len(group)),
                "accessible_count": float(np.count_nonzero(accessible_values)),
                "accessible_fraction": accessible_fraction(accessible_values),
                "radar_time_mae": time_metrics[0],
                "radar_time_rmse": time_metrics[1],
                "radar_time_bias": time_metrics[2],
                "radar_distance_mae": distance_metrics[0],
                "radar_distance_rmse": distance_metrics[1],
                "radar_distance_bias": distance_metrics[2],
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
    """Write event-level and aggregate CSV outputs."""

    data_dir = output_dir / "data"
    event_path = data_dir / "discrete_radar_reconstruction_events.csv"
    summary_path = data_dir / "discrete_radar_reconstruction_summary.csv"
    return _write_csv(event_rows, event_path), _write_csv(summary_rows, summary_path)


def _accessible_arrays(
    rows: list[dict[str, float]],
    truth_key: str,
    estimate_key: str,
) -> tuple[np.ndarray, np.ndarray]:
    accessible_rows = [row for row in rows if bool(row["accessible"])]
    truth = np.asarray([row[truth_key] for row in accessible_rows], dtype=float)
    estimate = np.asarray([row[estimate_key] for row in accessible_rows], dtype=float)
    return truth, estimate


def _plot_identity(ax: plt.Axes, truth: np.ndarray, estimate: np.ndarray) -> None:
    if truth.size == 0:
        return
    lower = float(min(np.min(truth), np.min(estimate)))
    upper = float(max(np.max(truth), np.max(estimate)))
    ax.plot([lower, upper], [lower, upper], color="black", linewidth=1.0)


def save_time_scatter(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save true radar time versus reconstructed radar time."""

    output_path = output_dir / "figures" / "discrete_radar_time_scatter.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    truth, estimate = _accessible_arrays(rows, "true_radar_time", "radar_time_hat")
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    ax.scatter(truth, estimate, s=8, alpha=0.4)
    _plot_identity(ax, truth, estimate)
    ax.set_xlabel("True radar time")
    ax.set_ylabel("Reconstructed radar time")
    ax.set_title("Causal-order-based radar time reconstruction")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_distance_scatter(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save true radar distance versus reconstructed radar distance."""

    output_path = output_dir / "figures" / "discrete_radar_distance_scatter.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    truth, estimate = _accessible_arrays(
        rows,
        "true_radar_distance",
        "radar_distance_hat",
    )
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    ax.scatter(truth, estimate, s=8, alpha=0.4)
    _plot_identity(ax, truth, estimate)
    ax.set_xlabel("True radar distance")
    ax.set_ylabel("Reconstructed radar distance")
    ax.set_title("Causal-order-based radar distance reconstruction")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_error_vs_ticks(summary_rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save radar reconstruction RMSE versus observer tick count."""

    output_path = output_dir / "figures" / "discrete_radar_error_vs_ticks.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for n_events in sorted({int(row["N"]) for row in summary_rows}):
        rows = [row for row in summary_rows if int(row["N"]) == n_events]
        ticks = np.asarray([row["tick_count"] for row in rows])
        time_rmse = np.asarray([row["radar_time_rmse"] for row in rows])
        distance_rmse = np.asarray([row["radar_distance_rmse"] for row in rows])
        ax.plot(ticks, time_rmse, marker="o", label=f"time N={n_events}")
        ax.plot(
            ticks,
            distance_rmse,
            marker="s",
            linestyle="--",
            label=f"distance N={n_events}",
        )
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel("RMSE")
    ax.set_title("Discrete radar error decreases with clock resolution")
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
    """Save accessible fraction versus observer tick count."""

    output_path = output_dir / "figures" / "discrete_radar_accessible_fraction.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    for n_events in sorted({int(row["N"]) for row in summary_rows}):
        rows = [row for row in summary_rows if int(row["N"]) == n_events]
        ticks = np.asarray([row["tick_count"] for row in rows])
        fractions = np.asarray([row["accessible_fraction"] for row in rows])
        ax.plot(ticks, fractions, marker="o", label=f"N={n_events}")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel("Accessible fraction")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("Radar accessibility from supplied observer chain")
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
) -> tuple[Path, Path, Path, Path]:
    """Save all discrete radar reconstruction figures."""

    return (
        save_time_scatter(event_rows, output_dir),
        save_distance_scatter(event_rows, output_dir),
        save_error_vs_ticks(summary_rows, output_dir),
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
        "Observer chain and clock labels are supplied protocol structure; "
        "hidden coordinates are used only for validation."
    )


if __name__ == "__main__":
    main()
