"""Controlled shortcut-injection sweep for echo-response motifs."""

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
from causal_spacetime_lab.state_change_echo_interference import (
    return_spectrum_report_for_motif,
    summarize_return_spectrum_reports,
)
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
from causal_spacetime_lab.state_change_echo_shortcuts import (
    ShortcutInjectionSpec,
    inject_shortcut_returns,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for controlled shortcut-injection sweeps."""

    reference_length: int = 64
    motif_count: int = 60
    delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    shortcut_probabilities: tuple[float, ...] = (0.0, 0.1, 0.3, 0.6, 1.0)
    shortcut_modes: tuple[str, ...] = (
        "target_to_early_reference",
        "decoy_path_to_early_reference",
    )
    repetitions: int = 5
    seed: int = 0
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

    parser = argparse.ArgumentParser(description="Echo shortcut injection sweep.")
    parser.add_argument("--reference-length", type=int, default=64)
    parser.add_argument("--motif-count", type=int, default=60)
    parser.add_argument("--delay-ranks", nargs="+", default=["3", "5", "8", "13"])
    parser.add_argument(
        "--shortcut-probabilities",
        nargs="+",
        default=["0.0", "0.1", "0.3", "0.6", "1.0"],
    )
    parser.add_argument(
        "--shortcut-modes",
        nargs="+",
        default=["target_to_early_reference", "decoy_path_to_early_reference"],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        motif_count=args.motif_count,
        delay_ranks=_parse_int_values(args.delay_ranks),
        shortcut_probabilities=_parse_float_values(args.shortcut_probabilities),
        shortcut_modes=tuple(args.shortcut_modes),
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


def _resolution_from_reports(rows: list[dict[str, float | str]]) -> dict[str, float]:
    recovered = np.asarray([float(row["recovered_delay_rank"]) for row in rows])
    reachable = np.isfinite(recovered)
    return echo_order_resolution_summary(
        np.where(reachable, recovered, -1).astype(int),
        reachable,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run shortcut injection sweep."""

    rows: list[dict[str, float | str]] = []
    for mode in config.shortcut_modes:
        for probability in config.shortcut_probabilities:
            for repetition in range(config.repetitions):
                seed = config.seed + 1000 * repetition + int(1000 * probability)
                network, reference = generate_reference_backbone_network(
                    config.reference_length
                )
                network, motifs = insert_multiple_echo_motifs(
                    network,
                    reference,
                    _random_specs(config, seed),
                )
                network, shortcut_records = inject_shortcut_returns(
                    network,
                    reference,
                    motifs,
                    ShortcutInjectionSpec(probability=probability, mode=mode),
                    seed=seed + 500,
                )
                closure = transitive_closure_dag(immediate_trigger_adjacency(network))
                reports = [
                    return_spectrum_report_for_motif(closure, reference, motif)
                    for motif in motifs
                ]
                summary = summarize_return_spectrum_reports(reports)
                order_summary = motif_order_recovery_rate(reports)
                resolution = _resolution_from_reports(reports)
                rows.append(
                    {
                        "reference_length": float(config.reference_length),
                        "motif_count_setting": float(config.motif_count),
                        "shortcut_probability": float(probability),
                        "shortcut_mode": mode,
                        "repetition": float(repetition),
                        "shortcut_injection_count": float(len(shortcut_records)),
                        **summary,
                        "order_inversion_rate": order_summary[
                            "order_inversion_rate"
                        ],
                        "order_agreement_rate": order_summary["order_agreement_rate"],
                        "tied_pair_fraction": resolution["tied_pair_fraction"],
                        "strict_order_pair_fraction": resolution[
                            "strict_order_pair_fraction"
                        ],
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write shortcut-injection sweep rows."""

    output_path = output_dir / "data" / "echo_shortcut_injection_sweep.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _mean_by_probability(
    rows: list[dict[str, float | str]],
    mode: str,
    key: str,
) -> tuple[list[float], list[float]]:
    probabilities = sorted({float(row["shortcut_probability"]) for row in rows})
    means: list[float] = []
    for probability in probabilities:
        values = [
            float(row[key])
            for row in rows
            if row["shortcut_mode"] == mode
            and float(row["shortcut_probability"]) == probability
            and np.isfinite(float(row[key]))
        ]
        means.append(float(np.mean(values)) if values else float("nan"))
    return probabilities, means


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path, Path]:
    """Save shortcut-injection sweep figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    specs = (
        (
            "echo_shortcut_fraction_vs_probability.png",
            "shortcut_fraction",
            "Shortcut fraction",
        ),
        (
            "echo_shortcut_depth_vs_probability.png",
            "mean_shortcut_depth",
            "Mean shortcut depth",
        ),
        (
            "echo_shortcut_order_inversion_vs_probability.png",
            "order_inversion_rate",
            "Order inversion against planted ranks",
        ),
    )
    paths: list[Path] = []
    modes = sorted({str(row["shortcut_mode"]) for row in rows})
    for filename, key, ylabel in specs:
        path = figure_dir / filename
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        for mode in modes:
            xs, ys = _mean_by_probability(rows, mode, key)
            ax.plot(xs, ys, marker="o", label=mode)
        ax.set_xlabel("Shortcut injection probability")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend(title="Injection mode")
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
        paths.append(path)
    return tuple(paths)  # type: ignore[return-value]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote echo shortcut injection sweep: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
