"""Statistical calibration for independent probe-pair interval reconstruction."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.error_models import (
    bias,
    binned_summary,
    binomial_tau_abs_std_delta_1p1,
    binomial_tau_rel_std_delta_1p1,
    mae,
    poisson_tau_abs_std_delta_1p1,
    poisson_tau_rel_std_delta_1p1,
    relative_rmse,
    rmse,
)
from causal_spacetime_lab.estimators import (
    estimate_tau_from_interval_count_1p1,
    global_density_1p1,
)
from causal_spacetime_lab.metrics import minkowski_tau_1p1
from causal_spacetime_lab.probes import (
    count_interval_events_for_probe_pairs_1p1,
    sample_probe_timelike_pairs_1p1,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_N_VALUES = (300, 600, 1200, 2400)
DEFAULT_TAU_BINS = (0.10, 0.20, 0.35, 0.50, 0.75, 1.00, 1.50, 2.00)
DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for independent probe-pair statistical calibration."""

    T: float = 2.0
    n_values: tuple[int, ...] = DEFAULT_N_VALUES
    repetitions: int = 5
    pairs_per_repetition: int = 200
    seed: int = 0
    min_tau: float = 0.10
    max_tau: float | None = None
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
            "Calibrate interval-count timelike reconstruction using independent "
            "probe endpoint pairs in known 1+1D Minkowski spacetime."
        )
    )
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument(
        "--n-values",
        nargs="+",
        default=[str(n) for n in DEFAULT_N_VALUES],
        help="Support event counts, supplied as space-separated values or commas.",
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--pairs-per-repetition", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--min-tau", type=float, default=0.10)
    parser.add_argument("--max-tau", type=float, default=None)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    return ExperimentConfig(
        T=args.T,
        n_values=_parse_n_values(args.n_values),
        repetitions=args.repetitions,
        pairs_per_repetition=args.pairs_per_repetition,
        seed=args.seed,
        min_tau=args.min_tau,
        max_tau=args.max_tau,
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
    if config.min_tau < 0:
        raise ValueError("min_tau must be non-negative")
    if config.max_tau is not None and config.max_tau <= 0:
        raise ValueError("max_tau must be positive")
    if config.max_tau is not None and config.min_tau > config.max_tau:
        raise ValueError("min_tau must be less than or equal to max_tau")


def _relative_error(estimate: float, truth: float) -> float:
    if truth == 0.0:
        return float("nan")
    return (estimate - truth) / truth


def _safe_relative_std(tau: float, rho: float) -> float:
    if tau <= 0.0:
        return float("nan")
    return poisson_tau_rel_std_delta_1p1(tau, rho)


def _safe_binomial_relative_std(
    tau: float,
    T: float,
    n_support: int,
    rho: float,
) -> float:
    if tau <= 0.0:
        return float("nan")
    return binomial_tau_rel_std_delta_1p1(tau, T, n_support, rho)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[
    list[dict[str, float]],
    list[dict[str, float]],
    list[dict[str, float]],
]:
    """Run statistical calibration and return pair, summary, and binned rows."""

    _validate_config(config)
    pair_rows: list[dict[str, float]] = []

    for n_index, n_support in enumerate(config.n_values):
        rho = global_density_1p1(n_support, config.T)
        poisson_abs_std = poisson_tau_abs_std_delta_1p1(rho)
        for repetition in range(config.repetitions):
            support_seed = config.seed + 10_000 * n_index + repetition
            probe_seed = config.seed + 20_000 * n_index + repetition
            support_events = sprinkle_1p1_causal_diamond(
                n_support,
                T=config.T,
                seed=support_seed,
            )
            probe_pairs = sample_probe_timelike_pairs_1p1(
                config.pairs_per_repetition,
                T=config.T,
                seed=probe_seed,
                min_tau=config.min_tau,
                max_tau=config.max_tau,
            )
            interval_counts = count_interval_events_for_probe_pairs_1p1(
                support_events,
                probe_pairs,
            )

            for pair_number, (pair, interval_count) in enumerate(
                zip(probe_pairs, interval_counts, strict=True)
            ):
                true_tau = minkowski_tau_1p1(pair[0], pair[1])
                tau_hat = estimate_tau_from_interval_count_1p1(
                    int(interval_count),
                    rho,
                )
                binomial_abs_std = binomial_tau_abs_std_delta_1p1(
                    true_tau,
                    config.T,
                    n_support,
                    rho,
                )
                pair_rows.append(
                    {
                        "N": float(n_support),
                        "repetition": float(repetition),
                        "pair_number": float(pair_number),
                        "rho": rho,
                        "p_t": float(pair[0, 0]),
                        "p_x": float(pair[0, 1]),
                        "q_t": float(pair[1, 0]),
                        "q_x": float(pair[1, 1]),
                        "true_tau": true_tau,
                        "interval_count": float(interval_count),
                        "tau_hat": tau_hat,
                        "absolute_error": abs(tau_hat - true_tau),
                        "signed_error": tau_hat - true_tau,
                        "relative_error": _relative_error(tau_hat, true_tau),
                        "poisson_abs_std": poisson_abs_std,
                        "poisson_rel_std": _safe_relative_std(true_tau, rho),
                        "binomial_abs_std": binomial_abs_std,
                        "binomial_rel_std": _safe_binomial_relative_std(
                            true_tau,
                            config.T,
                            n_support,
                            rho,
                        ),
                    }
                )

    summary_rows = summarize_results(pair_rows)
    binned_rows = summarize_by_tau_bins(pair_rows, config.T, config.min_tau)
    return pair_rows, summary_rows, binned_rows


def summarize_results(pair_rows: list[dict[str, float]]) -> list[dict[str, float]]:
    """Aggregate probe-pair calibration rows by support event count."""

    summary_rows: list[dict[str, float]] = []
    n_values = sorted({int(row["N"]) for row in pair_rows})
    for n_support in n_values:
        rows = [row for row in pair_rows if int(row["N"]) == n_support]
        true_tau = np.asarray([row["true_tau"] for row in rows])
        tau_hat = np.asarray([row["tau_hat"] for row in rows])
        observed_rmse = rmse(true_tau, tau_hat)
        mean_poisson_abs = float(np.mean([row["poisson_abs_std"] for row in rows]))
        mean_binomial_abs = float(np.mean([row["binomial_abs_std"] for row in rows]))

        summary_rows.append(
            {
                "N": float(n_support),
                "repetitions": float(len({int(row["repetition"]) for row in rows})),
                "pair_count": float(len(rows)),
                "rho": float(rows[0]["rho"]),
                "observed_mae": mae(true_tau, tau_hat),
                "observed_rmse": observed_rmse,
                "observed_bias": bias(true_tau, tau_hat),
                "observed_relative_rmse": relative_rmse(true_tau, tau_hat),
                "mean_poisson_abs_std": mean_poisson_abs,
                "mean_poisson_rel_std": float(
                    np.mean([row["poisson_rel_std"] for row in rows])
                ),
                "mean_binomial_abs_std": mean_binomial_abs,
                "mean_binomial_rel_std": float(
                    np.mean([row["binomial_rel_std"] for row in rows])
                ),
                "rmse_to_poisson_abs_std_ratio": observed_rmse / mean_poisson_abs,
                "rmse_to_binomial_abs_std_ratio": observed_rmse / mean_binomial_abs,
            }
        )
    return summary_rows


def _tau_bins(T: float, min_tau: float) -> np.ndarray:
    clipped = [edge for edge in DEFAULT_TAU_BINS if min_tau <= edge <= T]
    edges = sorted({min_tau, *clipped, T})
    if len(edges) < 2:
        edges = [0.0, T]
    return np.asarray(edges, dtype=float)


def summarize_by_tau_bins(
    pair_rows: list[dict[str, float]],
    T: float,
    min_tau: float,
) -> list[dict[str, float]]:
    """Aggregate relative-error behavior by true-tau bins and support count."""

    bins = _tau_bins(T, min_tau)
    rows_out: list[dict[str, float]] = []
    for n_support in sorted({int(row["N"]) for row in pair_rows}):
        rows = [row for row in pair_rows if int(row["N"]) == n_support]
        tau = np.asarray([row["true_tau"] for row in rows])
        rel_error = np.asarray([row["relative_error"] for row in rows])
        abs_error = np.asarray([row["signed_error"] for row in rows])
        binned_rel = binned_summary(tau, rel_error, bins)
        binned_abs = binned_summary(tau, abs_error, bins)
        for rel_row, abs_row in zip(binned_rel, binned_abs, strict=True):
            in_bin = (tau >= rel_row["bin_left"]) & (tau < rel_row["bin_right"])
            if rel_row["bin_right"] == bins[-1]:
                in_bin = (tau >= rel_row["bin_left"]) & (tau <= rel_row["bin_right"])
            bin_rows = [row for row, keep in zip(rows, in_bin, strict=True) if keep]
            rows_out.append(
                {
                    "N": float(n_support),
                    "bin_left": rel_row["bin_left"],
                    "bin_right": rel_row["bin_right"],
                    "bin_center": rel_row["bin_center"],
                    "pair_count": rel_row["count"],
                    "relative_error_bias": rel_row["mean"],
                    "relative_error_rmse": rel_row["rmse"],
                    "absolute_error_bias": abs_row["mean"],
                    "absolute_error_rmse": abs_row["rmse"],
                    "mean_poisson_rel_std": float(
                        np.mean([row["poisson_rel_std"] for row in bin_rows])
                    )
                    if bin_rows
                    else float("nan"),
                    "mean_binomial_rel_std": float(
                        np.mean([row["binomial_rel_std"] for row in bin_rows])
                    )
                    if bin_rows
                    else float("nan"),
                }
            )
    return rows_out


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
    binned_rows: list[dict[str, float]],
    output_dir: Path,
) -> tuple[Path, Path, Path]:
    """Write pair-level, aggregate, and binned CSV outputs."""

    data_dir = output_dir / "data"
    pair_path = data_dir / "probe_pair_statistical_calibration_pairs.csv"
    summary_path = data_dir / "probe_pair_statistical_calibration_summary.csv"
    binned_path = data_dir / "probe_pair_statistical_calibration_binned_by_tau.csv"
    return (
        _write_csv(pair_rows, pair_path),
        _write_csv(summary_rows, summary_path),
        _write_csv(binned_rows, binned_path),
    )


def save_tau_scatter(pair_rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save true tau versus interval-count tau estimate."""

    figure_path = output_dir / "figures" / "probe_pair_tau_scatter.png"
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    true_tau = np.asarray([row["true_tau"] for row in pair_rows])
    tau_hat = np.asarray([row["tau_hat"] for row in pair_rows])

    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    ax.scatter(true_tau, tau_hat, s=9, alpha=0.45)
    upper = float(max(np.max(true_tau), np.max(tau_hat)))
    ax.plot([0.0, upper], [0.0, upper], color="black", linewidth=1.0)
    ax.set_xlabel("True hidden proper time")
    ax.set_ylabel("Interval-count estimate")
    ax.set_title("Independent probe-pair tau reconstruction")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)
    return figure_path


def save_error_vs_tau(pair_rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save signed reconstruction error versus true tau."""

    figure_path = output_dir / "figures" / "probe_pair_error_vs_tau.png"
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    true_tau = np.asarray([row["true_tau"] for row in pair_rows])
    errors = np.asarray([row["signed_error"] for row in pair_rows])

    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.scatter(true_tau, errors, s=9, alpha=0.45)
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.set_xlabel("True hidden proper time")
    ax.set_ylabel("tau_hat - true_tau")
    ax.set_title("Probe-pair reconstruction error versus interval size")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)
    return figure_path


def save_rmse_vs_n(summary_rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save observed RMSE versus predicted finite-sampling standard deviations."""

    figure_path = output_dir / "figures" / "probe_pair_rmse_vs_N.png"
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    n_values = np.asarray([row["N"] for row in summary_rows])
    observed = np.asarray([row["observed_rmse"] for row in summary_rows])
    poisson = np.asarray([row["mean_poisson_abs_std"] for row in summary_rows])
    binomial = np.asarray([row["mean_binomial_abs_std"] for row in summary_rows])

    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.plot(n_values, observed, marker="o", label="Observed RMSE")
    ax.plot(n_values, poisson, marker="s", label="Poisson delta prediction")
    ax.plot(n_values, binomial, marker="^", label="Fixed-N binomial prediction")
    ax.set_xscale("log")
    ax.set_xlabel("Support event count N")
    ax.set_ylabel("Absolute tau error scale")
    ax.set_title("Probe-pair RMSE compared with sampling-noise predictions")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)
    return figure_path


def save_relative_error_by_tau_bin(
    binned_rows: list[dict[str, float]],
    output_dir: Path,
) -> Path:
    """Save relative RMSE by true-tau bin."""

    figure_path = output_dir / "figures" / "probe_pair_relative_error_by_tau_bin.png"
    figure_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for n_support in sorted({int(row["N"]) for row in binned_rows}):
        rows = [
            row
            for row in binned_rows
            if int(row["N"]) == n_support and row["pair_count"] > 0
        ]
        centers = np.asarray([row["bin_center"] for row in rows])
        rel_rmse_values = np.asarray([row["relative_error_rmse"] for row in rows])
        ax.plot(centers, rel_rmse_values, marker="o", label=f"N={n_support}")
    ax.set_xlabel("True-tau bin center")
    ax.set_ylabel("Relative RMSE")
    ax.set_title("Relative error is largest for small intervals")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)
    return figure_path


def save_figures(
    pair_rows: list[dict[str, float]],
    summary_rows: list[dict[str, float]],
    binned_rows: list[dict[str, float]],
    output_dir: Path,
) -> tuple[Path, Path, Path, Path]:
    """Save all statistical calibration figures."""

    return (
        save_tau_scatter(pair_rows, output_dir),
        save_error_vs_tau(pair_rows, output_dir),
        save_rmse_vs_n(summary_rows, output_dir),
        save_relative_error_by_tau_bin(binned_rows, output_dir),
    )


def main() -> None:
    config = parse_args()
    pair_rows, summary_rows, binned_rows = run_experiment(config)
    pair_path, summary_path, binned_path = write_outputs(
        pair_rows,
        summary_rows,
        binned_rows,
        config.output_dir,
    )
    figure_paths = save_figures(
        pair_rows,
        summary_rows,
        binned_rows,
        config.output_dir,
    )
    print(f"Wrote pair results: {pair_path}")
    print(f"Wrote summary: {summary_path}")
    print(f"Wrote binned summary: {binned_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")
    print(
        "Probe endpoints are independent of support events; density remains "
        "supplied scale information."
    )


if __name__ == "__main__":
    main()
