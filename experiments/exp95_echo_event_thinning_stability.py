"""Closure-preserving event thinning stability for echo return spectra."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change_echo_coarse_graining import (
    protected_indices_for_reference_and_motifs,
    remap_echo_motif_record_for_event_thinning,
    remap_reference_chain,
    restrict_transitive_order_to_retained_events,
    return_spectrum_stability_report,
    sample_retained_indices,
)
from causal_spacetime_lab.state_change_echo_interference import (
    return_delay_spectrum,
    return_spectrum_report_for_motif,
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
    """Configuration for closure-preserving event thinning."""

    reference_length: int = 64
    motif_count: int = 60
    delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    keep_probabilities: tuple[float, ...] = (1.0, 0.7, 0.4, 0.2)
    repetitions: int = 5
    seed: int = 0
    include_motif_paths_options: tuple[bool, ...] = (False, True)
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

    parser = argparse.ArgumentParser(description="Echo event thinning stability.")
    parser.add_argument("--reference-length", type=int, default=64)
    parser.add_argument("--motif-count", type=int, default=60)
    parser.add_argument("--delay-ranks", nargs="+", default=["3", "5", "8", "13"])
    parser.add_argument(
        "--keep-probabilities",
        nargs="+",
        default=["1.0", "0.7", "0.4", "0.2"],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--include-motif-paths-options",
        nargs="+",
        default=["false", "true"],
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    include_options = tuple(
        value.lower() in {"1", "true", "yes"}
        for value in args.include_motif_paths_options
    )
    return ExperimentConfig(
        reference_length=args.reference_length,
        motif_count=args.motif_count,
        delay_ranks=_parse_int_values(args.delay_ranks),
        keep_probabilities=_parse_float_values(args.keep_probabilities),
        repetitions=args.repetitions,
        seed=args.seed,
        include_motif_paths_options=include_options,
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


def _classification(row: dict[str, float | str]) -> tuple[float, float, float, float]:
    return (
        float(row["exact_recovery"]),
        float(row["early_shortcut"]),
        float(row["missing_return"]),
        float(row["late_recovery"]),
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run closure-preserving event thinning stability diagnostics."""

    rows: list[dict[str, float | str]] = []
    for include_paths in config.include_motif_paths_options:
        for keep_probability in config.keep_probabilities:
            for repetition in range(config.repetitions):
                seed = config.seed + 1000 * repetition + int(keep_probability * 1000)
                network, reference = generate_reference_backbone_network(
                    config.reference_length
                )
                network, motifs = insert_multiple_echo_motifs(
                    network,
                    reference,
                    _random_specs(config, seed),
                )
                closure = transitive_closure_dag(immediate_trigger_adjacency(network))
                protected = protected_indices_for_reference_and_motifs(
                    reference,
                    motifs,
                    include_motif_paths=include_paths,
                )
                retained = sample_retained_indices(
                    len(network.events),
                    keep_probability,
                    protected,
                    seed=seed + 500,
                )
                thinning = restrict_transitive_order_to_retained_events(
                    closure,
                    retained,
                )
                coarse_reference = remap_reference_chain(reference, thinning.old_to_new)
                jaccards: list[float] = []
                earliest_preserved: list[float] = []
                shifts: list[float] = []
                agreements: list[float] = []
                for motif in motifs:
                    coarse_motif = remap_echo_motif_record_for_event_thinning(
                        motif,
                        thinning.old_to_new,
                    )
                    if coarse_motif is None:
                        continue
                    baseline_spectrum = return_delay_spectrum(
                        closure,
                        reference,
                        motif.target_event_id,
                        motif.emission_position,
                    )
                    coarse_spectrum = return_delay_spectrum(
                        thinning.restricted_order_matrix,
                        coarse_reference,
                        coarse_motif.target_event_id,
                        coarse_motif.emission_position,
                    )
                    stability = return_spectrum_stability_report(
                        baseline_spectrum,
                        coarse_spectrum,
                    )
                    baseline_report = return_spectrum_report_for_motif(
                        closure,
                        reference,
                        motif,
                    )
                    coarse_report = return_spectrum_report_for_motif(
                        thinning.restricted_order_matrix,
                        coarse_reference,
                        coarse_motif,
                    )
                    jaccards.append(stability["spectrum_jaccard"])
                    earliest_preserved.append(stability["exact_earliest_preserved"])
                    if np.isfinite(stability["earliest_delay_shift"]):
                        shifts.append(stability["earliest_delay_shift"])
                    agreements.append(
                        float(
                            _classification(baseline_report)
                            == _classification(coarse_report)
                        )
                    )
                rows.append(
                    {
                        "reference_length": float(config.reference_length),
                        "motif_count_setting": float(config.motif_count),
                        "keep_probability": float(keep_probability),
                        "include_motif_paths": float(include_paths),
                        "repetition": float(repetition),
                        "retained_event_fraction": float(
                            retained.size / len(network.events)
                        ),
                        "mean_spectrum_jaccard": float(np.mean(jaccards)),
                        "earliest_preserved_fraction": float(
                            np.mean(earliest_preserved)
                        ),
                        "mean_earliest_delay_shift": float(np.mean(shifts))
                        if shifts
                        else float("nan"),
                        "classification_agreement_fraction": float(
                            np.mean(agreements)
                        ),
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write event thinning rows."""

    output_path = output_dir / "data" / "echo_event_thinning_stability.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _mean_by_keep(
    rows: list[dict[str, float | str]],
    include_paths: float,
    key: str,
) -> tuple[list[float], list[float]]:
    probabilities = sorted({float(row["keep_probability"]) for row in rows})
    means: list[float] = []
    for probability in probabilities:
        values = [
            float(row[key])
            for row in rows
            if float(row["keep_probability"]) == probability
            and float(row["include_motif_paths"]) == include_paths
            and np.isfinite(float(row[key]))
        ]
        means.append(float(np.mean(values)) if values else float("nan"))
    return probabilities, means


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Save event thinning figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    jaccard_path = figure_dir / "echo_event_thinning_spectrum_jaccard.png"
    earliest_path = figure_dir / "echo_event_thinning_earliest_preservation.png"
    for path, key, ylabel in (
        (jaccard_path, "mean_spectrum_jaccard", "Mean spectrum Jaccard"),
        (
            earliest_path,
            "earliest_preserved_fraction",
            "Earliest-return preserved fraction",
        ),
    ):
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        include_values = sorted(
            {float(row["include_motif_paths"]) for row in rows}
        )
        for include_paths in include_values:
            xs, ys = _mean_by_keep(rows, include_paths, key)
            label = f"include paths={bool(include_paths)}"
            ax.plot(xs, ys, marker="o", label=label)
        ax.set_xlabel("Event keep probability")
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
    return jaccard_path, earliest_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote echo event thinning stability: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
