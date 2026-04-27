"""Lightweight finite-size calibration for longest chains in probe intervals."""

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
    estimate_tau_from_longest_chain_1p1,
    global_density_1p1,
)
from causal_spacetime_lab.metrics import minkowski_tau_1p1
from causal_spacetime_lab.probes import sample_probe_timelike_pairs_1p1
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for the optional longest-chain calibration."""

    T: float = 2.0
    n_support: int = 600
    repetitions: int = 2
    pairs_per_repetition: int = 30
    seed: int = 0
    min_tau: float = 0.20
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Calibrate finite-size longest-chain behavior inside independent "
            "probe intervals."
        )
    )
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument("--n-support", type=int, default=600)
    parser.add_argument("--repetitions", type=int, default=2)
    parser.add_argument("--pairs-per-repetition", type=int, default=30)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--min-tau", type=float, default=0.20)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_support=args.n_support,
        repetitions=args.repetitions,
        pairs_per_repetition=args.pairs_per_repetition,
        seed=args.seed,
        min_tau=args.min_tau,
        output_dir=args.output_dir,
    )


def _validate_config(config: ExperimentConfig) -> None:
    if config.T <= 0:
        raise ValueError("T must be positive")
    if config.n_support <= 0:
        raise ValueError("n_support must be positive")
    if config.repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if config.pairs_per_repetition <= 0:
        raise ValueError("pairs_per_repetition must be positive")
    if config.min_tau < 0:
        raise ValueError("min_tau must be non-negative")


def _events_inside_probe_interval(
    support_events: np.ndarray,
    p: np.ndarray,
    q: np.ndarray,
) -> np.ndarray:
    dt_from_p = support_events[:, 0] - p[0]
    dx_from_p = support_events[:, 1] - p[1]
    dt_to_q = q[0] - support_events[:, 0]
    dx_to_q = q[1] - support_events[:, 1]
    inside = (
        (dt_from_p > 0.0)
        & (dt_from_p * dt_from_p >= dx_from_p * dx_from_p)
        & (dt_to_q > 0.0)
        & (dt_to_q * dt_to_q >= dx_to_q * dx_to_q)
    )
    return support_events[inside]


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run longest-chain calibration using independent probe intervals."""

    _validate_config(config)
    rows: list[dict[str, float]] = []
    rho = global_density_1p1(config.n_support, config.T)

    for repetition in range(config.repetitions):
        support_events = sprinkle_1p1_causal_diamond(
            config.n_support,
            T=config.T,
            seed=config.seed + repetition,
        )
        probe_pairs = sample_probe_timelike_pairs_1p1(
            config.pairs_per_repetition,
            T=config.T,
            seed=config.seed + 10_000 + repetition,
            min_tau=config.min_tau,
        )
        for pair_number, pair in enumerate(probe_pairs):
            p = pair[0]
            q = pair[1]
            interval_events = _events_inside_probe_interval(support_events, p, q)
            local_events = np.vstack((p, interval_events, q))
            local_causal_matrix = causal_matrix_1p1(local_events)
            chain_length = longest_chain_length(
                local_causal_matrix,
                start=0,
                end=local_events.shape[0] - 1,
                event_times=local_events[:, 0],
            )
            true_tau = minkowski_tau_1p1(p, q)
            predicted_chain_length = np.sqrt(2.0 * rho) * true_tau
            tau_chain = estimate_tau_from_longest_chain_1p1(
                chain_length,
                rho,
                chain_counts_endpoints=True,
            )
            rows.append(
                {
                    "N": float(config.n_support),
                    "repetition": float(repetition),
                    "pair_number": float(pair_number),
                    "rho": rho,
                    "true_tau": true_tau,
                    "interval_count": float(interval_events.shape[0]),
                    "chain_length_including_endpoints": float(chain_length),
                    "effective_chain_length_minus_endpoints": float(
                        max(chain_length - 2, 0)
                    ),
                    "asymptotic_chain_length": float(predicted_chain_length),
                    "tau_chain_estimate": tau_chain,
                    "tau_chain_error": tau_chain - true_tau,
                }
            )
    return rows


def write_summary(
    rows: list[dict[str, float]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write longest-chain calibration rows."""

    if not rows:
        raise RuntimeError("no rows to write")
    output_path = output_dir / "data" / "longest_chain_calibration_summary.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_plot(
    rows: list[dict[str, float]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Save longest-chain finite-size calibration plot."""

    output_path = output_dir / "figures" / "longest_chain_calibration.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    true_tau = np.asarray([row["true_tau"] for row in rows])
    chain_length = np.asarray(
        [row["chain_length_including_endpoints"] for row in rows]
    )
    effective_length = np.asarray(
        [row["effective_chain_length_minus_endpoints"] for row in rows]
    )
    predicted = np.asarray([row["asymptotic_chain_length"] for row in rows])

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    ax.scatter(true_tau, chain_length, s=14, alpha=0.6, label="Chain incl. endpoints")
    ax.scatter(
        true_tau,
        effective_length,
        s=14,
        alpha=0.6,
        label="Chain minus endpoints",
    )
    order = np.argsort(true_tau)
    ax.plot(true_tau[order], predicted[order], color="black", label="Asymptotic scale")
    ax.set_xlabel("True hidden proper time")
    ax.set_ylabel("Longest-chain length")
    ax.set_title("Finite-size longest-chain calibration in probe intervals")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    summary_path = write_summary(rows, config.output_dir)
    figure_path = save_plot(rows, config.output_dir)
    print(f"Wrote summary: {summary_path}")
    print(f"Wrote figure: {figure_path}")
    print(
        "Longest-chain behavior is finite-size sensitive; this calibration "
        "reports endpoint-convention effects rather than tuning them away."
    )


if __name__ == "__main__":
    main()

