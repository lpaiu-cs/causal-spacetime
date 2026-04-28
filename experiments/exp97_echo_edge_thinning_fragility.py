"""Immediate-edge thinning fragility for echo motif recovery."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change_echo_motif_validation import (
    motif_order_recovery_rate,
)
from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    insert_multiple_echo_motifs,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_edge_coarse_graining import (
    compare_recovery_before_after_edge_thinning,
    protected_edge_keys_for_motifs,
    protected_source_target_pairs_for_reference_chain,
    thin_immediate_trigger_edges,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for immediate-edge thinning diagnostics."""

    reference_length: int = 64
    motif_count: int = 60
    delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    removal_probabilities: tuple[float, ...] = (0.0, 0.05, 0.15, 0.30, 0.60)
    repetitions: int = 5
    seed: int = 0
    preserve_motif_edges_options: tuple[bool, ...] = (True, False)
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_float_values(values: list[str]) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    return tuple(parsed)


def _parse_int_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Echo edge thinning fragility.")
    parser.add_argument("--reference-length", type=int, default=64)
    parser.add_argument("--motif-count", type=int, default=60)
    parser.add_argument("--delay-ranks", nargs="+", default=["3", "5", "8", "13"])
    parser.add_argument(
        "--removal-probabilities",
        nargs="+",
        default=["0.0", "0.05", "0.15", "0.30", "0.60"],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--preserve-motif-edges-options",
        nargs="+",
        default=["true", "false"],
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    preserve_options = tuple(
        value.lower() in {"1", "true", "yes"}
        for value in args.preserve_motif_edges_options
    )
    return ExperimentConfig(
        reference_length=args.reference_length,
        motif_count=args.motif_count,
        delay_ranks=_parse_int_values(args.delay_ranks),
        removal_probabilities=_parse_float_values(args.removal_probabilities),
        repetitions=args.repetitions,
        seed=args.seed,
        preserve_motif_edges_options=preserve_options,
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


def _effect_fraction(rows: list[dict[str, float | str]], effect: str) -> float:
    if not rows:
        return float("nan")
    return float(np.mean([row["effect"] == effect for row in rows]))


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run immediate-edge thinning fragility diagnostics."""

    rows: list[dict[str, float | str]] = []
    for preserve in config.preserve_motif_edges_options:
        for probability in config.removal_probabilities:
            for repetition in range(config.repetitions):
                seed = config.seed + 1000 * repetition + int(probability * 1000)
                network, reference = generate_reference_backbone_network(
                    config.reference_length
                )
                network, motifs = insert_multiple_echo_motifs(
                    network,
                    reference,
                    _random_specs(config, seed),
                )
                before = transitive_closure_dag(immediate_trigger_adjacency(network))
                protected_pairs: set[tuple[int, int]] = set()
                if preserve:
                    reference_pairs = (
                        protected_source_target_pairs_for_reference_chain(reference)
                    )
                    protected_pairs |= reference_pairs
                    protected_pairs |= protected_edge_keys_for_motifs(motifs)
                thinned, removed = thin_immediate_trigger_edges(
                    network,
                    probability,
                    protected_source_target_pairs=protected_pairs,
                    seed=seed + 500,
                )
                after = transitive_closure_dag(immediate_trigger_adjacency(thinned))
                comparison = compare_recovery_before_after_edge_thinning(
                    before,
                    after,
                    reference,
                    motifs,
                )
                order_rows = [
                    {
                        "planted_delay_rank": row["planted_delay_rank"],
                        "recovered_delay_rank": row["after_delay_rank"],
                    }
                    for row in comparison
                ]
                order_summary = motif_order_recovery_rate(order_rows)
                exact_after = float(
                    np.mean(
                        [
                            np.isfinite(float(row["after_delay_rank"]))
                            and float(row["after_delay_rank"])
                            == float(row["planted_delay_rank"])
                            for row in comparison
                        ]
                    )
                )
                rows.append(
                    {
                        "reference_length": float(config.reference_length),
                        "motif_count_setting": float(config.motif_count),
                        "removal_probability": float(probability),
                        "preserve_motif_edges": float(preserve),
                        "repetition": float(repetition),
                        "removed_edge_count": float(removed),
                        "missing_fraction": _effect_fraction(
                            comparison,
                            "became_missing",
                        ),
                        "delayed_fraction": _effect_fraction(comparison, "delayed"),
                        "exact_recovery_fraction": float(exact_after),
                        "order_inversion_rate": order_summary[
                            "order_inversion_rate"
                        ],
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write edge thinning rows."""

    output_path = output_dir / "data" / "echo_edge_thinning_fragility.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _mean_by_probability(
    rows: list[dict[str, float | str]],
    preserve: float,
    key: str,
) -> tuple[list[float], list[float]]:
    probabilities = sorted({float(row["removal_probability"]) for row in rows})
    means: list[float] = []
    for probability in probabilities:
        values = [
            float(row[key])
            for row in rows
            if float(row["removal_probability"]) == probability
            and float(row["preserve_motif_edges"]) == preserve
            and np.isfinite(float(row[key]))
        ]
        means.append(float(np.mean(values)) if values else float("nan"))
    return probabilities, means


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Save edge thinning figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    missing_path = figure_dir / "echo_edge_thinning_missing_fraction.png"
    exact_path = figure_dir / "echo_edge_thinning_exact_recovery.png"
    for path, key, ylabel in (
        (missing_path, "missing_fraction", "Missing recovery fraction"),
        (exact_path, "exact_recovery_fraction", "Exact recovery fraction"),
    ):
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        for preserve in sorted({float(row["preserve_motif_edges"]) for row in rows}):
            xs, ys = _mean_by_probability(rows, preserve, key)
            ax.plot(xs, ys, marker="o", label=f"preserve motif edges={bool(preserve)}")
        ax.set_xlabel("Immediate-edge removal probability")
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
    return missing_path, exact_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote echo edge thinning fragility: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
