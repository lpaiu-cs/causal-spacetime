"""Motif density and delay-rank resolution diagnostics."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change_echo_diagnostics import (
    echo_order_resolution_summary,
)
from causal_spacetime_lab.state_change_echo_motif_validation import (
    motif_order_recovery_rate,
    motif_recovery_report,
)
from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    insert_multiple_echo_motifs,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)

DEFAULT_OUTPUT_DIR = Path("outputs")
DELAY_RANK_SETS = {
    "small": (2, 3),
    "medium": (2, 3, 5, 8),
    "large": (2, 3, 5, 8, 13, 21),
}


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for motif density and rank resolution."""

    reference_length: int = 64
    motif_counts: tuple[int, ...] = (10, 30, 100, 200)
    delay_rank_sets: tuple[str, ...] = ("small", "medium", "large")
    repetitions: int = 5
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Echo motif density resolution.")
    parser.add_argument("--reference-length", type=int, default=64)
    parser.add_argument("--motif-counts", nargs="+", default=["10", "30", "100", "200"])
    parser.add_argument(
        "--delay-rank-sets",
        nargs="+",
        default=["small", "medium", "large"],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        motif_counts=_parse_int_values(args.motif_counts),
        delay_rank_sets=tuple(args.delay_rank_sets),
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _random_specs(
    reference_length: int,
    motif_count: int,
    delay_ranks: tuple[int, ...],
    seed: int,
) -> list[EchoMotifSpec]:
    rng = np.random.default_rng(seed)
    valid_delays = np.asarray(
        [delay for delay in delay_ranks if delay < reference_length]
    )
    specs: list[EchoMotifSpec] = []
    for _ in range(motif_count):
        delay = int(rng.choice(valid_delays))
        emission = int(rng.integers(0, reference_length - delay))
        specs.append(EchoMotifSpec(emission, delay))
    return specs


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run motif density and rank-resolution diagnostics."""

    rows: list[dict[str, float | str]] = []
    for motif_count in config.motif_counts:
        for delay_set_name in config.delay_rank_sets:
            if delay_set_name not in DELAY_RANK_SETS:
                raise ValueError(f"unknown delay-rank set: {delay_set_name}")
            delay_ranks = DELAY_RANK_SETS[delay_set_name]
            for repetition in range(config.repetitions):
                seed = config.seed + 1000 * motif_count + repetition
                network, reference = generate_reference_backbone_network(
                    config.reference_length
                )
                specs = _random_specs(
                    config.reference_length,
                    motif_count,
                    delay_ranks,
                    seed,
                )
                network, motifs = insert_multiple_echo_motifs(network, reference, specs)
                closure = transitive_closure_dag(immediate_trigger_adjacency(network))
                motif_rows = [
                    motif_recovery_report(closure, reference, motif)
                    for motif in motifs
                ]
                recovered = np.asarray(
                    [float(row["recovered_delay_rank"]) for row in motif_rows],
                    dtype=float,
                )
                reachable = np.isfinite(recovered)
                resolution = echo_order_resolution_summary(
                    np.where(reachable, recovered, -1).astype(int),
                    reachable,
                )
                order_summary = motif_order_recovery_rate(motif_rows)
                rows.append(
                    {
                        "reference_length": float(config.reference_length),
                        "motif_count": float(motif_count),
                        "delay_rank_set": delay_set_name,
                        "available_delay_rank_count": float(len(delay_ranks)),
                        "repetition": float(repetition),
                        "distinct_recovered_delay_ranks": float(
                            len(set(recovered[reachable].astype(int).tolist()))
                        )
                        if np.any(reachable)
                        else 0.0,
                        "tied_target_fraction": resolution["tied_pair_fraction"],
                        "strict_order_pair_fraction": resolution[
                            "strict_order_pair_fraction"
                        ],
                        "order_inversion_rate": order_summary[
                            "order_inversion_rate"
                        ],
                        "order_agreement_rate": order_summary["order_agreement_rate"],
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write motif density-resolution rows."""

    output_path = output_dir / "data" / "echo_motif_density_resolution.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _mean_by_set(
    rows: list[dict[str, float | str]],
    key: str,
) -> dict[str, tuple[list[float], list[float]]]:
    result: dict[str, tuple[list[float], list[float]]] = {}
    for delay_set in sorted({str(row["delay_rank_set"]) for row in rows}):
        counts = sorted(
            {
                float(row["motif_count"])
                for row in rows
                if row["delay_rank_set"] == delay_set
            }
        )
        means: list[float] = []
        for count in counts:
            values = [
                float(row[key])
                for row in rows
                if row["delay_rank_set"] == delay_set
                and float(row["motif_count"]) == count
                and np.isfinite(float(row[key]))
            ]
            means.append(float(np.mean(values)) if values else float("nan"))
        result[delay_set] = (counts, means)
    return result


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Save motif density-resolution figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    tie_path = figure_dir / "echo_motif_tie_fraction.png"
    strict_path = figure_dir / "echo_motif_strict_order_fraction.png"
    for path, key, ylabel in (
        (tie_path, "tied_target_fraction", "Tied target-pair fraction"),
        (strict_path, "strict_order_pair_fraction", "Strictly ordered pair fraction"),
    ):
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        for label, (xs, ys) in _mean_by_set(rows, key).items():
            ax.plot(xs, ys, marker="o", label=label)
        ax.set_xlabel("Planted motif count")
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, alpha=0.3)
        ax.legend(title="Delay-rank set")
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
    return tie_path, strict_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote echo motif density resolution: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
