"""Random thinning and density-rescaling stability experiment."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.coarse_graining import (
    rescaled_density_after_thinning,
    thin_events,
)
from causal_spacetime_lab.dimension import (
    estimate_myrheim_meyer_dimension,
    relation_fraction,
)
from causal_spacetime_lab.estimators import bias, rmse
from causal_spacetime_lab.metrics import causal_diamond_volume_1p1
from causal_spacetime_lab.probes import sample_probe_timelike_pairs_1p1
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_KEEP_PROBABILITIES = (1.0, 0.75, 0.5, 0.25)
DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for thinning stability validation."""

    T: float = 2.0
    n_events: int = 4800
    repetitions: int = 5
    keep_probabilities: tuple[float, ...] = DEFAULT_KEEP_PROBABILITIES
    pairs_per_repetition: int = 200
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


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
        description="Test reconstruction stability under random thinning."
    )
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument("--N", type=int, default=4800)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument(
        "--keep-probabilities",
        nargs="+",
        default=[str(value) for value in DEFAULT_KEEP_PROBABILITIES],
    )
    parser.add_argument("--pairs-per-repetition", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_events=args.N,
        repetitions=args.repetitions,
        keep_probabilities=_parse_float_values(args.keep_probabilities, "keep"),
        pairs_per_repetition=args.pairs_per_repetition,
        seed=args.seed,
        output_dir=args.output_dir,
    )


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


def _flat_probe_volume(p: np.ndarray, q: np.ndarray) -> float:
    dt = q[0] - p[0]
    dx = q[1] - p[1]
    return float(0.5 * (dt * dt - dx * dx))


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float]], list[dict[str, float]]]:
    """Run thinning stability experiment."""

    if config.T <= 0.0:
        raise ValueError("T must be positive")
    if config.n_events <= 0:
        raise ValueError("N must be positive")
    if config.repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if config.pairs_per_repetition <= 0:
        raise ValueError("pairs_per_repetition must be positive")
    if any(not 0.0 <= value <= 1.0 for value in config.keep_probabilities):
        raise ValueError("keep probabilities must satisfy 0 <= p <= 1")

    global_volume = causal_diamond_volume_1p1(config.T)
    original_density = config.n_events / global_volume
    pair_rows: list[dict[str, float]] = []
    run_rows: list[dict[str, float]] = []
    for repetition in range(config.repetitions):
        support_events = sprinkle_1p1_causal_diamond(
            config.n_events,
            T=config.T,
            seed=config.seed + repetition,
        )
        probe_pairs = sample_probe_timelike_pairs_1p1(
            config.pairs_per_repetition,
            T=config.T,
            seed=config.seed + 10_000 + repetition,
            min_tau=0.10,
        )
        true_volumes = np.asarray([_flat_probe_volume(p, q) for p, q in probe_pairs])
        for keep_probability in config.keep_probabilities:
            thinned = thin_events(
                support_events,
                keep_probability,
                seed=config.seed + 100_000 * repetition + int(1000 * keep_probability),
            )
            corrected_density = rescaled_density_after_thinning(
                original_density,
                keep_probability,
            )
            realized_density = thinned.shape[0] / global_volume
            dimension_estimate = float("nan")
            if thinned.shape[0] >= 100:
                fraction = relation_fraction(causal_matrix_1p1(thinned))
                dimension_estimate = estimate_myrheim_meyer_dimension(fraction)
            run_rows.append(
                {
                    "keep_probability": float(keep_probability),
                    "repetition": float(repetition),
                    "thinned_count": float(thinned.shape[0]),
                    "dimension_estimate": dimension_estimate,
                }
            )
            for pair_index, ((p, q), true_volume) in enumerate(
                zip(probe_pairs, true_volumes, strict=True)
            ):
                count = _interval_count(thinned, p, q)
                corrected = (
                    count / corrected_density if corrected_density > 0 else np.nan
                )
                realized = count / realized_density if realized_density > 0 else np.nan
                uncorrected = count / original_density
                pair_rows.append(
                    {
                        "keep_probability": float(keep_probability),
                        "repetition": float(repetition),
                        "pair_index": float(pair_index),
                        "thinned_count": float(thinned.shape[0]),
                        "interval_count": float(count),
                        "true_volume": float(true_volume),
                        "corrected_volume": float(corrected),
                        "realized_density_volume": float(realized),
                        "uncorrected_volume": float(uncorrected),
                    }
                )
    return pair_rows, summarize_results(pair_rows, run_rows, config)


def _finite_pair_arrays(
    rows: list[dict[str, float]],
    estimate_key: str,
) -> tuple[np.ndarray, np.ndarray]:
    true = np.asarray([row["true_volume"] for row in rows], dtype=float)
    estimate = np.asarray([row[estimate_key] for row in rows], dtype=float)
    finite = np.isfinite(true) & np.isfinite(estimate)
    return true[finite], estimate[finite]


def summarize_results(
    pair_rows: list[dict[str, float]],
    run_rows: list[dict[str, float]],
    config: ExperimentConfig,
) -> list[dict[str, float]]:
    """Aggregate thinning stability results."""

    summary: list[dict[str, float]] = []
    for keep_probability in sorted({row["keep_probability"] for row in pair_rows}):
        pairs = [
            row for row in pair_rows if row["keep_probability"] == keep_probability
        ]
        runs = [row for row in run_rows if row["keep_probability"] == keep_probability]
        true_c, corrected = _finite_pair_arrays(pairs, "corrected_volume")
        true_r, realized = _finite_pair_arrays(pairs, "realized_density_volume")
        true_u, uncorrected = _finite_pair_arrays(pairs, "uncorrected_volume")
        dimensions = np.asarray(
            [row["dimension_estimate"] for row in runs],
            dtype=float,
        )
        dimensions = dimensions[np.isfinite(dimensions)]
        summary.append(
            {
                "keep_probability": float(keep_probability),
                "repetitions": float(config.repetitions),
                "mean_thinned_count": float(
                    np.mean([row["thinned_count"] for row in runs])
                ),
                "expected_thinned_count": float(keep_probability * config.n_events),
                "corrected_volume_rmse": rmse(true_c, corrected),
                "corrected_volume_bias": bias(true_c, corrected),
                "realized_density_volume_rmse": rmse(true_r, realized),
                "realized_density_volume_bias": bias(true_r, realized),
                "uncorrected_volume_rmse": rmse(true_u, uncorrected),
                "uncorrected_volume_bias": bias(true_u, uncorrected),
                "dimension_mean": float(np.mean(dimensions))
                if dimensions.size
                else float("nan"),
                "dimension_std": float(np.std(dimensions, ddof=1))
                if dimensions.size > 1
                else 0.0,
                "dimension_bias": float(np.mean(dimensions) - 2.0)
                if dimensions.size
                else float("nan"),
            }
        )
    return summary


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
    """Write thinning CSV outputs."""

    data_dir = output_dir / "data"
    return (
        _write_csv(pair_rows, data_dir / "thinning_coarse_graining_pairs.csv"),
        _write_csv(summary_rows, data_dir / "thinning_coarse_graining_summary.csv"),
    )


def save_volume_error_vs_keep(
    summary_rows: list[dict[str, float]],
    output_dir: Path,
) -> Path:
    """Save volume RMSE versus keep probability."""

    output_path = (
        output_dir / "figures" / "thinning_volume_error_vs_keep_probability.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    keep = np.asarray([row["keep_probability"] for row in summary_rows], dtype=float)
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    for key, label in [
        ("corrected_volume_rmse", "expected-density corrected"),
        ("realized_density_volume_rmse", "realized-density corrected"),
        ("uncorrected_volume_rmse", "uncorrected original density"),
    ]:
        values = np.asarray([row[key] for row in summary_rows], dtype=float)
        ax.plot(keep, values, marker="o", label=label)
    ax.set_xlabel("Keep probability")
    ax.set_ylabel("Volume RMSE")
    ax.set_title("Thinning volume reconstruction")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_dimension_vs_keep(
    summary_rows: list[dict[str, float]],
    output_dir: Path,
) -> Path:
    """Save Myrheim-Meyer dimension versus keep probability."""

    output_path = output_dir / "figures" / "thinning_dimension_vs_keep_probability.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    keep = np.asarray([row["keep_probability"] for row in summary_rows], dtype=float)
    dims = np.asarray([row["dimension_mean"] for row in summary_rows], dtype=float)
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    ax.plot(keep, dims, marker="o")
    ax.axhline(2.0, color="black", linewidth=0.9)
    ax.set_xlabel("Keep probability")
    ax.set_ylabel("Estimated dimension")
    ax.set_title("Dimension under thinning")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_count_ratio(
    summary_rows: list[dict[str, float]],
    config: ExperimentConfig,
) -> Path:
    """Save realized count ratio versus expected keep probability."""

    output_path = config.output_dir / "figures" / "thinning_count_ratio.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    keep = np.asarray([row["keep_probability"] for row in summary_rows], dtype=float)
    ratios = np.asarray(
        [row["mean_thinned_count"] / config.n_events for row in summary_rows],
        dtype=float,
    )
    fig, ax = plt.subplots(figsize=(6.5, 4.6))
    ax.plot(keep, ratios, marker="o", label="realized")
    ax.plot(keep, keep, color="black", linewidth=1.0, label="expected")
    ax.set_xlabel("Keep probability")
    ax.set_ylabel("Count ratio")
    ax.set_title("Random thinning count ratio")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_figures(
    summary_rows: list[dict[str, float]],
    config: ExperimentConfig,
) -> tuple[Path, Path, Path]:
    """Save thinning figures."""

    return (
        save_volume_error_vs_keep(summary_rows, config.output_dir),
        save_dimension_vs_keep(summary_rows, config.output_dir),
        save_count_ratio(summary_rows, config),
    )


def main() -> None:
    config = parse_args()
    pair_rows, summary_rows = run_experiment(config)
    pair_path, summary_path = write_outputs(pair_rows, summary_rows, config.output_dir)
    figure_paths = save_figures(summary_rows, config)
    print(f"Wrote pair results: {pair_path}")
    print(f"Wrote summary: {summary_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
