"""Non-tautological timelike pair reconstruction validation in 1+1D."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.chains import longest_chain_length
from causal_spacetime_lab.estimators import (
    bias,
    estimate_tau_from_interval_count_1p1,
    estimate_tau_from_longest_chain_1p1,
    global_density_1p1,
    mae,
    relative_rmse,
    rmse,
)
from causal_spacetime_lab.intervals import alexandrov_interval_indices
from causal_spacetime_lab.metrics import minkowski_tau_1p1
from causal_spacetime_lab.sampling import sample_timelike_pairs
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_N_VALUES = (300, 600, 1200, 2400)
DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for the timelike pair reconstruction experiment."""

    T: float = 2.0
    n_values: tuple[int, ...] = DEFAULT_N_VALUES
    repetitions: int = 5
    pairs_per_repetition: int = 100
    seed: int = 0
    max_pairs_for_chain: int = 25
    skip_chain: bool = False
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_n_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    if not parsed:
        raise argparse.ArgumentTypeError("at least one N value is required")
    if any(n <= 0 for n in parsed):
        raise argparse.ArgumentTypeError("all N values must be positive")
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Validate timelike proper-time reconstruction on internal event "
            "pairs in a known 1+1D Minkowski causal diamond."
        )
    )
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument(
        "--n-values",
        nargs="+",
        default=[str(n) for n in DEFAULT_N_VALUES],
        help="Event counts, supplied as space-separated values or comma lists.",
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--pairs-per-repetition", type=int, default=100)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-pairs-for-chain", type=int, default=25)
    parser.add_argument("--skip-chain", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    return ExperimentConfig(
        T=args.T,
        n_values=_parse_n_values(args.n_values),
        repetitions=args.repetitions,
        pairs_per_repetition=args.pairs_per_repetition,
        seed=args.seed,
        max_pairs_for_chain=args.max_pairs_for_chain,
        skip_chain=args.skip_chain,
        output_dir=args.output_dir,
    )


def _validate_config(config: ExperimentConfig) -> None:
    if config.T <= 0:
        raise ValueError("T must be positive")
    if any(n <= 0 for n in config.n_values):
        raise ValueError("all n_values must be positive")
    if config.repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if config.pairs_per_repetition <= 0:
        raise ValueError("pairs_per_repetition must be positive")
    if config.max_pairs_for_chain < 0:
        raise ValueError("max_pairs_for_chain must be non-negative")


def _chain_length_for_pair(
    causal_matrix: np.ndarray,
    events: np.ndarray,
    start: int,
    end: int,
) -> int:
    interval = alexandrov_interval_indices(causal_matrix, start, end)
    subset = np.concatenate(([start], interval, [end]))
    submatrix = causal_matrix[np.ix_(subset, subset)]
    return longest_chain_length(
        submatrix,
        start=0,
        end=subset.size - 1,
        event_times=events[subset, 0],
    )


def _nan_row_values() -> tuple[float, float, float, float]:
    return (float("nan"), float("nan"), float("nan"), float("nan"))


def _summary_metrics(
    rows: list[dict[str, float]],
    estimate_key: str,
) -> tuple[float, float, float, float]:
    true_values = np.asarray([row["true_tau"] for row in rows], dtype=float)
    estimates = np.asarray([row[estimate_key] for row in rows], dtype=float)
    keep = np.isfinite(estimates)
    if not np.any(keep):
        return _nan_row_values()
    return (
        mae(true_values[keep], estimates[keep]),
        rmse(true_values[keep], estimates[keep]),
        bias(true_values[keep], estimates[keep]),
        relative_rmse(true_values[keep], estimates[keep]),
    )


def _relative_error(estimate: float, truth: float) -> float:
    if truth == 0.0 or not np.isfinite(estimate):
        return float("nan")
    return (estimate - truth) / truth


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float]], list[dict[str, float]]]:
    """Run internal-pair reconstruction and return pair-level and summary rows."""

    _validate_config(config)
    pair_rows: list[dict[str, float]] = []

    for n_index, n_events in enumerate(config.n_values):
        rho = global_density_1p1(n_events, config.T)
        for repetition in range(config.repetitions):
            repetition_seed = config.seed + 10_000 * n_index + repetition
            pair_seed = config.seed + 20_000 * n_index + repetition
            events = sprinkle_1p1_causal_diamond(
                n_events,
                T=config.T,
                seed=repetition_seed,
            )
            causal_matrix = causal_matrix_1p1(events)
            pairs = sample_timelike_pairs(
                events,
                causal_matrix,
                num_pairs=config.pairs_per_repetition,
                seed=pair_seed,
            )

            chain_count = 0
            for pair_number, (i, j) in enumerate(pairs):
                true_tau = minkowski_tau_1p1(events[i], events[j])
                interval_count = int(
                    alexandrov_interval_indices(causal_matrix, i, j).size
                )
                volume_tau = estimate_tau_from_interval_count_1p1(
                    interval_count,
                    rho,
                )

                chain_length = float("nan")
                chain_tau = float("nan")
                if (
                    not config.skip_chain
                    and chain_count < config.max_pairs_for_chain
                ):
                    chain_length = float(
                        _chain_length_for_pair(causal_matrix, events, i, j)
                    )
                    chain_tau = estimate_tau_from_longest_chain_1p1(
                        int(chain_length),
                        rho,
                        chain_counts_endpoints=True,
                    )
                    chain_count += 1

                pair_rows.append(
                    {
                        "N": float(n_events),
                        "repetition": float(repetition),
                        "pair_number": float(pair_number),
                        "i": float(i),
                        "j": float(j),
                        "rho": rho,
                        "true_tau": true_tau,
                        "interval_count": float(interval_count),
                        "tau_volume_estimate": volume_tau,
                        "tau_volume_abs_error": abs(volume_tau - true_tau),
                        "tau_volume_relative_error": _relative_error(
                            volume_tau,
                            true_tau,
                        ),
                        "chain_length": chain_length,
                        "tau_chain_estimate": chain_tau,
                        "tau_chain_abs_error": abs(chain_tau - true_tau)
                        if np.isfinite(chain_tau)
                        else float("nan"),
                        "tau_chain_relative_error": _relative_error(
                            chain_tau,
                            true_tau,
                        ),
                    }
                )

    summary_rows = summarize_results(pair_rows)
    return pair_rows, summary_rows


def summarize_results(pair_rows: list[dict[str, float]]) -> list[dict[str, float]]:
    """Aggregate pair-level rows by event count."""

    summary_rows: list[dict[str, float]] = []
    n_values = sorted({int(row["N"]) for row in pair_rows})
    for n_events in n_values:
        rows = [row for row in pair_rows if int(row["N"]) == n_events]
        volume_metrics = _summary_metrics(rows, "tau_volume_estimate")
        chain_metrics = _summary_metrics(rows, "tau_chain_estimate")
        chain_pair_count = int(
            np.count_nonzero(
                np.isfinite([row["tau_chain_estimate"] for row in rows])
            )
        )

        summary_rows.append(
            {
                "N": float(n_events),
                "repetition_count": float(
                    len({int(row["repetition"]) for row in rows})
                ),
                "pair_count": float(len(rows)),
                "chain_pair_count": float(chain_pair_count),
                "tau_volume_mae": volume_metrics[0],
                "tau_volume_rmse": volume_metrics[1],
                "tau_volume_bias": volume_metrics[2],
                "tau_volume_relative_rmse": volume_metrics[3],
                "tau_chain_mae": chain_metrics[0],
                "tau_chain_rmse": chain_metrics[1],
                "tau_chain_bias": chain_metrics[2],
                "tau_chain_relative_rmse": chain_metrics[3],
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
    pair_rows: list[dict[str, float]],
    summary_rows: list[dict[str, float]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write pair-level and aggregate CSV outputs."""

    data_dir = output_dir / "data"
    pair_path = data_dir / "timelike_pair_reconstruction_pairs.csv"
    summary_path = data_dir / "timelike_pair_reconstruction_summary.csv"
    return _write_csv(pair_rows, pair_path), _write_csv(summary_rows, summary_path)


def _plot_identity(ax: plt.Axes, values: np.ndarray) -> None:
    if values.size == 0:
        return
    lower = 0.0
    upper = float(np.nanmax(values))
    ax.plot([lower, upper], [lower, upper], color="black", linewidth=1.0)


def save_scatter_plot(pair_rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save true-tau versus reconstructed-tau scatter plots."""

    figure_path = output_dir / "figures" / "timelike_pair_reconstruction_scatter.png"
    figure_path.parent.mkdir(parents=True, exist_ok=True)

    true_tau = np.asarray([row["true_tau"] for row in pair_rows], dtype=float)
    volume_tau = np.asarray(
        [row["tau_volume_estimate"] for row in pair_rows],
        dtype=float,
    )
    chain_tau = np.asarray(
        [row["tau_chain_estimate"] for row in pair_rows],
        dtype=float,
    )
    has_chain = np.any(np.isfinite(chain_tau))

    if has_chain:
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), sharex=True, sharey=True)
        ax_volume, ax_chain = axes
    else:
        fig, ax_volume = plt.subplots(figsize=(5.8, 4.8))
        ax_chain = None

    ax_volume.scatter(true_tau, volume_tau, s=10, alpha=0.55)
    _plot_identity(ax_volume, np.concatenate((true_tau, volume_tau)))
    ax_volume.set_title("Interval-count estimate")
    ax_volume.set_xlabel("True hidden proper time")
    ax_volume.set_ylabel("Estimated proper time")
    ax_volume.grid(True, alpha=0.25)

    if ax_chain is not None:
        keep = np.isfinite(chain_tau)
        ax_chain.scatter(true_tau[keep], chain_tau[keep], s=10, alpha=0.55)
        _plot_identity(ax_chain, np.concatenate((true_tau[keep], chain_tau[keep])))
        ax_chain.set_title("Longest-chain estimate")
        ax_chain.set_xlabel("True hidden proper time")
        ax_chain.grid(True, alpha=0.25)

    fig.suptitle("Internal timelike pair reconstruction in 1+1D Minkowski spacetime")
    fig.tight_layout()
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)
    return figure_path


def save_error_plot(summary_rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save reconstruction error versus event count."""

    figure_path = output_dir / "figures" / "timelike_pair_reconstruction_error_vs_N.png"
    figure_path.parent.mkdir(parents=True, exist_ok=True)

    n_values = np.asarray([row["N"] for row in summary_rows], dtype=float)
    volume_rmse = np.asarray([row["tau_volume_rmse"] for row in summary_rows])
    chain_rmse = np.asarray([row["tau_chain_rmse"] for row in summary_rows])

    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.plot(n_values, volume_rmse, marker="o", label="Interval cardinality")
    if np.any(np.isfinite(chain_rmse)):
        ax.plot(n_values, chain_rmse, marker="s", label="Longest chain")
    ax.set_xscale("log")
    ax.set_xlabel("Sprinkled internal events N")
    ax.set_ylabel("RMSE")
    ax.set_title("Timelike pair reconstruction error versus sprinkling density")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)
    return figure_path


def save_figures(
    pair_rows: list[dict[str, float]],
    summary_rows: list[dict[str, float]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Save experiment figures."""

    return save_scatter_plot(pair_rows, output_dir), save_error_plot(
        summary_rows,
        output_dir,
    )


def main() -> None:
    config = parse_args()
    pair_rows, summary_rows = run_experiment(config)
    pair_path, summary_path = write_outputs(pair_rows, summary_rows, config.output_dir)
    scatter_path, error_path = save_figures(pair_rows, summary_rows, config.output_dir)
    print(f"Wrote pair results: {pair_path}")
    print(f"Wrote summary: {summary_path}")
    print(f"Wrote scatter figure: {scatter_path}")
    print(f"Wrote error figure: {error_path}")
    print(
        "Density is supplied scale information; longest-chain normalization is "
        "asymptotic and finite-size sensitive."
    )


if __name__ == "__main__":
    main()
