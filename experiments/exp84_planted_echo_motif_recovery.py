"""Recovery of planted echo-delay ranks in clean controlled motif networks."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change_echo_motif_validation import (
    motif_order_recovery_rate,
    motif_recovery_report,
    summarize_motif_recovery,
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


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for planted echo motif recovery."""

    reference_lengths: tuple[int, ...] = (16, 32, 64)
    motif_counts: tuple[int, ...] = (10, 30, 60)
    delay_ranks: tuple[int, ...] = (2, 3, 5, 8)
    repetitions: int = 5
    seed: int = 0
    outward_steps: int = 1
    return_steps: int = 1
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Planted echo motif recovery.")
    parser.add_argument("--reference-lengths", nargs="+", default=["16", "32", "64"])
    parser.add_argument("--motif-counts", nargs="+", default=["10", "30", "60"])
    parser.add_argument("--delay-ranks", nargs="+", default=["2", "3", "5", "8"])
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--outward-steps", type=int, default=1)
    parser.add_argument("--return-steps", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_lengths=_parse_int_values(args.reference_lengths),
        motif_counts=_parse_int_values(args.motif_counts),
        delay_ranks=_parse_int_values(args.delay_ranks),
        repetitions=args.repetitions,
        seed=args.seed,
        outward_steps=args.outward_steps,
        return_steps=args.return_steps,
        output_dir=args.output_dir,
    )


def _random_specs(
    reference_length: int,
    motif_count: int,
    delay_ranks: tuple[int, ...],
    outward_steps: int,
    return_steps: int,
    seed: int,
) -> list[EchoMotifSpec]:
    rng = np.random.default_rng(seed)
    valid_delays = np.asarray(
        [delay for delay in delay_ranks if delay < reference_length]
    )
    if valid_delays.size == 0:
        return []
    specs: list[EchoMotifSpec] = []
    for _ in range(motif_count):
        delay = int(rng.choice(valid_delays))
        emission = int(rng.integers(0, reference_length - delay))
        specs.append(
            EchoMotifSpec(
                emission_position=emission,
                planted_delay_rank=delay,
                outward_steps=outward_steps,
                return_steps=return_steps,
            )
        )
    return specs


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    """Run clean planted echo motif recovery diagnostics."""

    motif_rows: list[dict[str, float | str]] = []
    summary_rows: list[dict[str, float | str]] = []
    for reference_length in config.reference_lengths:
        for motif_count in config.motif_counts:
            for repetition in range(config.repetitions):
                seed = (
                    config.seed
                    + 1000 * reference_length
                    + 10 * motif_count
                    + repetition
                )
                network, reference = generate_reference_backbone_network(
                    reference_length
                )
                specs = _random_specs(
                    reference_length,
                    motif_count,
                    config.delay_ranks,
                    config.outward_steps,
                    config.return_steps,
                    seed,
                )
                network, motifs = insert_multiple_echo_motifs(network, reference, specs)
                closure = transitive_closure_dag(immediate_trigger_adjacency(network))
                rows = [
                    motif_recovery_report(closure, reference, motif)
                    for motif in motifs
                ]
                for index, row in enumerate(rows):
                    row.update(
                        {
                            "reference_length": float(reference_length),
                            "motif_count_setting": float(motif_count),
                            "repetition": float(repetition),
                            "motif_index": float(index),
                        }
                    )
                    motif_rows.append(row)
                summary = summarize_motif_recovery(rows)
                order_summary = motif_order_recovery_rate(rows)
                summary_rows.append(
                    {
                        "reference_length": float(reference_length),
                        "motif_count_setting": float(motif_count),
                        "repetition": float(repetition),
                        **summary,
                        **order_summary,
                    }
                )
    return motif_rows, summary_rows


def write_outputs(
    motif_rows: list[dict[str, float | str]],
    summary_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write motif-level and summary rows."""

    data_dir = output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    motif_path = data_dir / "planted_echo_motif_recovery.csv"
    summary_path = data_dir / "planted_echo_motif_recovery_summary.csv"
    with motif_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(motif_rows[0]))
        writer.writeheader()
        writer.writerows(motif_rows)
    with summary_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    return motif_path, summary_path


def _mean_by_reference_length(
    rows: list[dict[str, float | str]],
    key: str,
) -> tuple[list[float], list[float]]:
    lengths = sorted({float(row["reference_length"]) for row in rows})
    means: list[float] = []
    for length in lengths:
        values = [
            float(row[key])
            for row in rows
            if float(row["reference_length"]) == length
            and np.isfinite(float(row[key]))
        ]
        means.append(float(np.mean(values)) if values else float("nan"))
    return lengths, means


def save_figures(
    summary_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Save planted motif recovery figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    exact_path = figure_dir / "planted_echo_exact_recovery_fraction.png"
    order_path = figure_dir / "planted_echo_order_recovery.png"

    for path, key, ylabel in (
        (exact_path, "exact_recovery_fraction", "Exact recovery fraction"),
        (order_path, "order_agreement_rate", "Delay-order agreement rate"),
    ):
        xs, ys = _mean_by_reference_length(summary_rows, key)
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        ax.plot(xs, ys, marker="o")
        ax.set_xlabel("Reference-chain length")
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
    return exact_path, order_path


def main() -> None:
    config = parse_args()
    motif_rows, summary_rows = run_experiment(config)
    motif_path, summary_path = write_outputs(
        motif_rows,
        summary_rows,
        config.output_dir,
    )
    figure_paths = save_figures(summary_rows, config.output_dir)
    print(f"Wrote planted echo motif recovery: {motif_path}")
    print(f"Wrote planted echo motif recovery summary: {summary_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
