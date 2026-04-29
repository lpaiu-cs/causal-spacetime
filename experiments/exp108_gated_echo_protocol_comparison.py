"""Compare earliest-return and predeclared gated echo protocols."""

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
from causal_spacetime_lab.state_change_gated_echo import (
    gated_response_signature_from_motifs,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)
from causal_spacetime_lab.state_change_response_layers import (
    planted_response_signature_from_motifs,
)
from causal_spacetime_lab.state_change_response_signature import (
    echo_response_signature_from_motifs,
    signature_reachable_fraction,
    signature_tie_fraction,
)
from causal_spacetime_lab.state_change_response_signature_comparison import (
    compare_response_signatures,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for gated echo protocol comparison."""

    reference_length: int = 64
    motif_count: int = 60
    delay_ranks: tuple[int, ...] = (5, 8, 13, 21)
    shortcut_probability_values: tuple[float, ...] = (0.0, 0.2, 0.5)
    gate_delay_ranks: tuple[int, ...] = (1, 3, 5, 8)
    repetitions: int = 5
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def _parse_float_values(values: list[str]) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Gated echo protocol comparison.")
    parser.add_argument("--reference-length", type=int, default=64)
    parser.add_argument("--motif-count", type=int, default=60)
    parser.add_argument("--delay-ranks", nargs="+", default=["5", "8", "13", "21"])
    parser.add_argument(
        "--shortcut-probability-values",
        nargs="+",
        default=["0.0", "0.2", "0.5"],
    )
    parser.add_argument("--gate-delay-ranks", nargs="+", default=["1", "3", "5", "8"])
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        motif_count=args.motif_count,
        delay_ranks=_parse_int_values(args.delay_ranks),
        shortcut_probability_values=_parse_float_values(
            args.shortcut_probability_values
        ),
        gate_delay_ranks=_parse_int_values(args.gate_delay_ranks),
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


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run gated versus earliest-return comparisons."""

    rows: list[dict[str, float | str]] = []
    for probability in config.shortcut_probability_values:
        for repetition in range(config.repetitions):
            seed = config.seed + repetition + int(probability * 1000)
            network, reference = generate_reference_backbone_network(
                config.reference_length
            )
            network, motifs = insert_multiple_echo_motifs(
                network,
                reference,
                _random_specs(config, seed),
            )
            network, _ = inject_shortcut_returns(
                network,
                reference,
                motifs,
                ShortcutInjectionSpec(probability=probability),
                seed=seed + 1000,
            )
            closure = transitive_closure_dag(immediate_trigger_adjacency(network))
            planted = planted_response_signature_from_motifs(motifs)
            earliest = echo_response_signature_from_motifs(
                closure,
                reference,
                motifs,
                label="earliest",
            )
            reports = [
                return_spectrum_report_for_motif(closure, reference, motif)
                for motif in motifs
            ]
            shortcut_summary = summarize_return_spectrum_reports(reports)
            earliest_comparison = compare_response_signatures(planted, earliest)
            rows.append(
                {
                    "shortcut_probability": float(probability),
                    "gate_delay_rank": float("nan"),
                    "protocol": "earliest",
                    "repetition": float(repetition),
                    "reachable_fraction": signature_reachable_fraction(earliest),
                    "pair_agreement_fraction": earliest_comparison[
                        "pair_agreement_fraction"
                    ],
                    "tie_fraction": signature_tie_fraction(earliest),
                    "shortcut_fraction": shortcut_summary["shortcut_fraction"],
                }
            )
            for gate in config.gate_delay_ranks:
                gated = gated_response_signature_from_motifs(
                    closure,
                    reference,
                    motifs,
                    gate,
                    label=f"gated_{gate}",
                )
                comparison = compare_response_signatures(planted, gated)
                rows.append(
                    {
                        "shortcut_probability": float(probability),
                        "gate_delay_rank": float(gate),
                        "protocol": "gated",
                        "repetition": float(repetition),
                        "reachable_fraction": signature_reachable_fraction(gated),
                        "pair_agreement_fraction": comparison[
                            "pair_agreement_fraction"
                        ],
                        "tie_fraction": signature_tie_fraction(gated),
                        "shortcut_fraction": shortcut_summary["shortcut_fraction"],
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write experiment CSV."""

    path = output_dir / "data" / "gated_echo_protocol_comparison.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _group_key(row: dict[str, float | str]) -> str:
    if row["protocol"] == "earliest":
        return "earliest"
    return f"gate_{int(float(row['gate_delay_rank']))}"


def _means_by_probability(
    rows: list[dict[str, float | str]],
    key: str,
    protocol_key: str,
) -> tuple[list[float], list[float]]:
    probabilities = sorted({float(row["shortcut_probability"]) for row in rows})
    means: list[float] = []
    for probability in probabilities:
        values = [
            float(row[key])
            for row in rows
            if _group_key(row) == protocol_key
            and float(row["shortcut_probability"]) == probability
            and np.isfinite(float(row[key]))
        ]
        means.append(float(np.mean(values)) if values else float("nan"))
    return probabilities, means


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save summary figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    protocol_keys = sorted({_group_key(row) for row in rows})
    paths: list[Path] = []
    for value, filename, ylabel in [
        ("pair_agreement_fraction", "gated_echo_pair_agreement.png", "pair agreement"),
        ("reachable_fraction", "gated_echo_reachability.png", "reachable fraction"),
    ]:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for protocol_key in protocol_keys:
            xs, ys = _means_by_probability(rows, value, protocol_key)
            ax.plot(xs, ys, marker="o", label=protocol_key)
        ax.set_xlabel("shortcut probability")
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, alpha=0.3)
        ax.legend()
        path = figure_dir / filename
        fig.tight_layout()
        fig.savefig(path, dpi=160)
        plt.close(fig)
        paths.append(path)
    return paths


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote gated echo protocol comparison: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()

