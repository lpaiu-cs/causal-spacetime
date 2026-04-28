"""Motif path-length robustness for shortcut-return stress tests."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change_echo_interference import (
    return_spectrum_report_for_motif,
    summarize_return_spectrum_reports,
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
SHORTCUT_MODES = ("target_to_early_reference", "return_path_to_early_reference")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for motif path-length robustness diagnostics."""

    reference_length: int = 64
    motif_count: int = 40
    delay_ranks: tuple[int, ...] = (5, 8, 13)
    outward_steps_values: tuple[int, ...] = (0, 1, 3)
    return_steps_values: tuple[int, ...] = (0, 1, 3)
    shortcut_probability: float = 0.3
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

    parser = argparse.ArgumentParser(description="Echo motif path-length robustness.")
    parser.add_argument("--reference-length", type=int, default=64)
    parser.add_argument("--motif-count", type=int, default=40)
    parser.add_argument("--delay-ranks", nargs="+", default=["5", "8", "13"])
    parser.add_argument(
        "--outward-steps-values",
        nargs="+",
        default=["0", "1", "3"],
    )
    parser.add_argument(
        "--return-steps-values",
        nargs="+",
        default=["0", "1", "3"],
    )
    parser.add_argument("--shortcut-probability", type=float, default=0.3)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        motif_count=args.motif_count,
        delay_ranks=_parse_int_values(args.delay_ranks),
        outward_steps_values=_parse_int_values(args.outward_steps_values),
        return_steps_values=_parse_int_values(args.return_steps_values),
        shortcut_probability=args.shortcut_probability,
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _random_specs(
    config: ExperimentConfig,
    outward_steps: int,
    return_steps: int,
    seed: int,
) -> list[EchoMotifSpec]:
    rng = np.random.default_rng(seed)
    valid_delays = np.asarray(
        [delay for delay in config.delay_ranks if delay < config.reference_length],
        dtype=int,
    )
    specs: list[EchoMotifSpec] = []
    for _ in range(config.motif_count):
        delay = int(rng.choice(valid_delays))
        emission = int(rng.integers(0, config.reference_length - delay))
        specs.append(
            EchoMotifSpec(
                emission_position=emission,
                planted_delay_rank=delay,
                outward_steps=outward_steps,
                return_steps=return_steps,
            )
        )
    return specs


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run path-length robustness diagnostics."""

    rows: list[dict[str, float | str]] = []
    for outward_steps in config.outward_steps_values:
        for return_steps in config.return_steps_values:
            for mode in SHORTCUT_MODES:
                for repetition in range(config.repetitions):
                    seed = (
                        config.seed
                        + 10000 * outward_steps
                        + 1000 * return_steps
                        + 100 * repetition
                    )
                    network, reference = generate_reference_backbone_network(
                        config.reference_length
                    )
                    specs = _random_specs(
                        config,
                        outward_steps,
                        return_steps,
                        seed,
                    )
                    network, motifs = insert_multiple_echo_motifs(
                        network,
                        reference,
                        specs,
                    )
                    network, shortcut_records = inject_shortcut_returns(
                        network,
                        reference,
                        motifs,
                        ShortcutInjectionSpec(
                            probability=config.shortcut_probability,
                            mode=mode,
                        ),
                        seed=seed + 500,
                    )
                    closure = transitive_closure_dag(
                        immediate_trigger_adjacency(network)
                    )
                    reports = [
                        return_spectrum_report_for_motif(closure, reference, motif)
                        for motif in motifs
                    ]
                    summary = summarize_return_spectrum_reports(reports)
                    rows.append(
                        {
                            "reference_length": float(config.reference_length),
                            "motif_count_setting": float(config.motif_count),
                            "outward_steps": float(outward_steps),
                            "return_steps": float(return_steps),
                            "shortcut_mode": mode,
                            "shortcut_probability": float(
                                config.shortcut_probability
                            ),
                            "repetition": float(repetition),
                            "shortcut_injection_count": float(
                                len(shortcut_records)
                            ),
                            **summary,
                        }
                    )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write path-length robustness rows."""

    output_path = output_dir / "data" / "echo_motif_path_length_robustness.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _mean_by_path(
    rows: list[dict[str, float | str]],
    key: str,
) -> tuple[list[str], list[float]]:
    labels = sorted(
        {
            f"o{int(float(row['outward_steps']))}/r{int(float(row['return_steps']))}"
            for row in rows
        }
    )
    means: list[float] = []
    for label in labels:
        values = []
        for row in rows:
            row_label = (
                f"o{int(float(row['outward_steps']))}/"
                f"r{int(float(row['return_steps']))}"
            )
            value = float(row[key])
            if row_label == label and np.isfinite(value):
                values.append(value)
        means.append(float(np.mean(values)) if values else float("nan"))
    return labels, means


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Save path-length robustness figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    shortcut_path = figure_dir / "echo_path_length_shortcut_fraction.png"
    spectrum_path = figure_dir / "echo_path_length_spectrum_size.png"
    for path, key, ylabel in (
        (shortcut_path, "shortcut_fraction", "Shortcut fraction"),
        (spectrum_path, "mean_spectrum_size", "Mean return spectrum size"),
    ):
        labels, means = _mean_by_path(rows, key)
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        ax.bar(labels, means)
        ax.set_xlabel("Outward/return intermediate steps")
        ax.set_ylabel(ylabel)
        ax.grid(True, axis="y", alpha=0.3)
        fig.autofmt_xdate(rotation=25)
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
    return shortcut_path, spectrum_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote echo motif path-length robustness: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
