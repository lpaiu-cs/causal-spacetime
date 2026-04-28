"""Reference-chain subsampling and echo-rank resolution loss."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.ordinal import order_inversion_rate
from causal_spacetime_lab.state_change_echo import echo_delay_rank_for_emission
from causal_spacetime_lab.state_change_echo_coarse_graining import (
    coarse_emission_position_for_motif,
    expected_coarse_delay_rank_for_motif,
    subsample_reference_chain_positions,
)
from causal_spacetime_lab.state_change_echo_diagnostics import (
    echo_order_resolution_summary,
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
    """Configuration for reference-chain subsampling diagnostics."""

    reference_length: int = 128
    motif_count: int = 100
    delay_ranks: tuple[int, ...] = (2, 3, 5, 8, 13, 21)
    strides: tuple[int, ...] = (1, 2, 4, 8)
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

    parser = argparse.ArgumentParser(description="Echo reference subsampling.")
    parser.add_argument("--reference-length", type=int, default=128)
    parser.add_argument("--motif-count", type=int, default=100)
    parser.add_argument(
        "--delay-ranks",
        nargs="+",
        default=["2", "3", "5", "8", "13", "21"],
    )
    parser.add_argument("--strides", nargs="+", default=["1", "2", "4", "8"])
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        motif_count=args.motif_count,
        delay_ranks=_parse_int_values(args.delay_ranks),
        strides=_parse_int_values(args.strides),
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _random_specs(config: ExperimentConfig, seed: int) -> list[EchoMotifSpec]:
    rng = np.random.default_rng(seed)
    valid_delays = np.asarray(
        [delay for delay in config.delay_ranks if delay < config.reference_length],
        dtype=int,
    )
    specs: list[EchoMotifSpec] = []
    for _ in range(config.motif_count):
        delay = int(rng.choice(valid_delays))
        emission = int(rng.integers(0, config.reference_length - delay))
        specs.append(EchoMotifSpec(emission, delay))
    return specs


def _safe_inversion(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 2:
        return float("nan")
    return order_inversion_rate(a, b, ignore_ties=True)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run reference-chain subsampling resolution diagnostics."""

    rows: list[dict[str, float | str]] = []
    for stride in config.strides:
        for repetition in range(config.repetitions):
            seed = config.seed + 1000 * repetition + stride
            network, reference = generate_reference_backbone_network(
                config.reference_length
            )
            network, motifs = insert_multiple_echo_motifs(
                network,
                reference,
                _random_specs(config, seed),
            )
            closure = transitive_closure_dag(immediate_trigger_adjacency(network))
            emissions = np.asarray([motif.emission_position for motif in motifs])
            subsampling = subsample_reference_chain_positions(
                reference,
                stride,
                protected_positions=emissions,
            )
            planted: list[float] = []
            expected: list[float] = []
            recovered: list[float] = []
            for motif in motifs:
                coarse_emission = coarse_emission_position_for_motif(
                    motif,
                    subsampling,
                )
                expected_delay = expected_coarse_delay_rank_for_motif(
                    motif,
                    subsampling,
                )
                if coarse_emission is None or expected_delay is None:
                    continue
                actual_delay = echo_delay_rank_for_emission(
                    closure,
                    subsampling.subsampled_reference_chain,
                    motif.target_event_id,
                    coarse_emission,
                )
                if actual_delay is None:
                    continue
                planted.append(float(motif.planted_delay_rank))
                expected.append(float(expected_delay))
                recovered.append(float(actual_delay))
            recovered_array = np.asarray(recovered, dtype=float)
            expected_array = np.asarray(expected, dtype=float)
            planted_array = np.asarray(planted, dtype=float)
            reachable = np.ones(recovered_array.size, dtype=bool)
            resolution = echo_order_resolution_summary(
                recovered_array.astype(int),
                reachable,
            )
            rows.append(
                {
                    "reference_length": float(config.reference_length),
                    "motif_count_setting": float(config.motif_count),
                    "stride": float(stride),
                    "repetition": float(repetition),
                    "retained_reference_fraction": float(
                        subsampling.subsampled_reference_chain.size / reference.size
                    ),
                    "exact_expected_coarse_recovery_fraction": float(
                        np.mean(recovered_array == expected_array)
                    )
                    if recovered_array.size
                    else float("nan"),
                    "tied_pair_fraction": resolution["tied_pair_fraction"],
                    "strict_order_pair_fraction": resolution[
                        "strict_order_pair_fraction"
                    ],
                    "order_inversion_vs_planted": _safe_inversion(
                        planted_array,
                        recovered_array,
                    ),
                    "order_inversion_vs_expected": _safe_inversion(
                        expected_array,
                        recovered_array,
                    ),
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write reference subsampling rows."""

    output_path = output_dir / "data" / "echo_reference_subsampling_resolution.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _mean_by_stride(
    rows: list[dict[str, float | str]],
    key: str,
) -> tuple[list[float], list[float]]:
    strides = sorted({float(row["stride"]) for row in rows})
    means: list[float] = []
    for stride in strides:
        values = [
            float(row[key])
            for row in rows
            if float(row["stride"]) == stride and np.isfinite(float(row[key]))
        ]
        means.append(float(np.mean(values)) if values else float("nan"))
    return strides, means


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Save reference subsampling figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    ties_path = figure_dir / "echo_reference_subsampling_ties.png"
    order_path = figure_dir / "echo_reference_subsampling_order_loss.png"
    for path, key, ylabel in (
        (ties_path, "tied_pair_fraction", "Tied target-pair fraction"),
        (
            order_path,
            "order_inversion_vs_planted",
            "Order inversion versus planted ranks",
        ),
    ):
        xs, ys = _mean_by_stride(rows, key)
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        ax.plot(xs, ys, marker="o")
        ax.set_xlabel("Reference-chain subsampling stride")
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
    return ties_path, order_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote echo reference subsampling resolution: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
