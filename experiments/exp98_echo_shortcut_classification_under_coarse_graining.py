"""Shortcut-classification stability under echo coarse-graining choices."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change_echo_coarse_graining import (
    coarse_emission_position_for_motif,
    expected_coarse_delay_rank_for_motif,
    protected_indices_for_reference_and_motifs,
    remap_echo_motif_record_for_event_thinning,
    remap_reference_chain,
    restrict_transitive_order_to_retained_events,
    sample_retained_indices,
    subsample_reference_chain_positions,
)
from causal_spacetime_lab.state_change_echo_interference import (
    return_spectrum_report_for_motif,
)
from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifRecord,
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
from causal_spacetime_lab.state_change_edge_coarse_graining import (
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
    """Configuration for shortcut classification coarse-graining."""

    reference_length: int = 64
    motif_count: int = 60
    delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    shortcut_probability: float = 0.3
    event_keep_probability: float = 0.4
    edge_removal_probability: float = 0.15
    reference_strides: tuple[int, ...] = (1, 2, 4)
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

    parser = argparse.ArgumentParser(
        description="Echo shortcut classification under coarse-graining."
    )
    parser.add_argument("--reference-length", type=int, default=64)
    parser.add_argument("--motif-count", type=int, default=60)
    parser.add_argument("--delay-ranks", nargs="+", default=["3", "5", "8", "13"])
    parser.add_argument("--shortcut-probability", type=float, default=0.3)
    parser.add_argument("--event-keep-probability", type=float, default=0.4)
    parser.add_argument("--edge-removal-probability", type=float, default=0.15)
    parser.add_argument("--reference-strides", nargs="+", default=["1", "2", "4"])
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        motif_count=args.motif_count,
        delay_ranks=_parse_int_values(args.delay_ranks),
        shortcut_probability=args.shortcut_probability,
        event_keep_probability=args.event_keep_probability,
        edge_removal_probability=args.edge_removal_probability,
        reference_strides=_parse_int_values(args.reference_strides),
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


def _classification(row: dict[str, float | str]) -> str:
    if float(row["exact_recovery"]) == 1.0:
        return "exact"
    if float(row["early_shortcut"]) == 1.0:
        return "shortcut"
    if float(row["missing_return"]) == 1.0:
        return "missing"
    if float(row["late_recovery"]) == 1.0:
        return "late"
    return "other"


def _report_counts(
    baseline_labels: list[str],
    coarse_labels: list[str],
) -> dict[str, float]:
    total = len(coarse_labels)
    if total == 0:
        return {
            "classification_agreement_fraction": float("nan"),
            "exact_fraction": float("nan"),
            "shortcut_fraction": float("nan"),
            "missing_fraction": float("nan"),
            "late_fraction": float("nan"),
        }
    agreements = [
        a == b for a, b in zip(baseline_labels, coarse_labels, strict=True)
    ]
    return {
        "classification_agreement_fraction": float(np.mean(agreements)),
        "exact_fraction": float(np.mean([label == "exact" for label in coarse_labels])),
        "shortcut_fraction": float(
            np.mean([label == "shortcut" for label in coarse_labels])
        ),
        "missing_fraction": float(
            np.mean([label == "missing" for label in coarse_labels])
        ),
        "late_fraction": float(np.mean([label == "late" for label in coarse_labels])),
    }


def _reports_for_event_thinning(
    closure: np.ndarray,
    reference: np.ndarray,
    motifs: list[EchoMotifRecord],
    event_keep_probability: float,
    seed: int,
) -> list[dict[str, float | str]]:
    protected = protected_indices_for_reference_and_motifs(reference, motifs)
    retained = sample_retained_indices(
        closure.shape[0],
        event_keep_probability,
        protected,
        seed=seed,
    )
    thinning = restrict_transitive_order_to_retained_events(closure, retained)
    coarse_reference = remap_reference_chain(reference, thinning.old_to_new)
    reports: list[dict[str, float | str]] = []
    for motif in motifs:
        coarse_motif = remap_echo_motif_record_for_event_thinning(
            motif,
            thinning.old_to_new,
        )
        if coarse_motif is None:
            continue
        reports.append(
            return_spectrum_report_for_motif(
                thinning.restricted_order_matrix,
                coarse_reference,
                coarse_motif,
            )
        )
    return reports


def _reports_for_reference_subsampling(
    closure: np.ndarray,
    reference: np.ndarray,
    motifs: list[EchoMotifRecord],
    stride: int,
) -> list[dict[str, float | str]]:
    emissions = np.asarray([motif.emission_position for motif in motifs])
    subsampling = subsample_reference_chain_positions(
        reference,
        stride,
        protected_positions=emissions,
    )
    reports: list[dict[str, float | str]] = []
    for motif in motifs:
        coarse_emission = coarse_emission_position_for_motif(motif, subsampling)
        expected_delay = expected_coarse_delay_rank_for_motif(motif, subsampling)
        if coarse_emission is None or expected_delay is None:
            continue
        coarse_motif = replace(
            motif,
            emission_position=coarse_emission,
            planted_delay_rank=expected_delay,
        )
        reports.append(
            return_spectrum_report_for_motif(
                closure,
                subsampling.subsampled_reference_chain,
                coarse_motif,
            )
        )
    return reports


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run shortcut classification stability diagnostics."""

    rows: list[dict[str, float | str]] = []
    for repetition in range(config.repetitions):
        seed = config.seed + repetition
        network, reference = generate_reference_backbone_network(
            config.reference_length
        )
        network, motifs = insert_multiple_echo_motifs(
            network,
            reference,
            _random_specs(config, seed + 1000),
        )
        network, _ = inject_shortcut_returns(
            network,
            reference,
            motifs,
            ShortcutInjectionSpec(probability=config.shortcut_probability),
            seed=seed + 2000,
        )
        closure = transitive_closure_dag(immediate_trigger_adjacency(network))
        baseline_reports = [
            return_spectrum_report_for_motif(closure, reference, motif)
            for motif in motifs
        ]
        baseline_labels = [_classification(row) for row in baseline_reports]

        event_reports = _reports_for_event_thinning(
            closure,
            reference,
            motifs,
            config.event_keep_probability,
            seed + 3000,
        )
        event_labels = [_classification(row) for row in event_reports]
        rows.append(
            {
                "coarse_graining_type": "closure_preserving_event_thinning",
                "reference_stride": float("nan"),
                "repetition": float(repetition),
                **_report_counts(baseline_labels, event_labels),
            }
        )

        protected_pairs = protected_source_target_pairs_for_reference_chain(reference)
        edge_network, _ = thin_immediate_trigger_edges(
            network,
            config.edge_removal_probability,
            protected_source_target_pairs=protected_pairs,
            seed=seed + 4000,
        )
        edge_closure = transitive_closure_dag(immediate_trigger_adjacency(edge_network))
        edge_reports = [
            return_spectrum_report_for_motif(edge_closure, reference, motif)
            for motif in motifs
        ]
        edge_labels = [_classification(row) for row in edge_reports]
        rows.append(
            {
                "coarse_graining_type": "immediate_edge_thinning",
                "reference_stride": float("nan"),
                "repetition": float(repetition),
                **_report_counts(baseline_labels, edge_labels),
            }
        )

        for stride in config.reference_strides:
            reference_reports = _reports_for_reference_subsampling(
                closure,
                reference,
                motifs,
                stride,
            )
            reference_labels = [_classification(row) for row in reference_reports]
            rows.append(
                {
                    "coarse_graining_type": "reference_chain_subsampling",
                    "reference_stride": float(stride),
                    "repetition": float(repetition),
                    **_report_counts(baseline_labels, reference_labels),
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write shortcut classification coarse-graining rows."""

    output_path = (
        output_dir
        / "data"
        / "echo_shortcut_classification_under_coarse_graining.csv"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _mean_by_type(
    rows: list[dict[str, float | str]],
    key: str,
) -> tuple[list[str], list[float]]:
    labels = sorted(
        {
            str(row["coarse_graining_type"])
            if not np.isfinite(float(row["reference_stride"]))
            else f"reference_stride_{int(float(row['reference_stride']))}"
            for row in rows
        }
    )
    means: list[float] = []
    for label in labels:
        values = []
        for row in rows:
            row_label = (
                str(row["coarse_graining_type"])
                if not np.isfinite(float(row["reference_stride"]))
                else f"reference_stride_{int(float(row['reference_stride']))}"
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
    """Save shortcut classification figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    agreement_path = figure_dir / "echo_shortcut_classification_agreement.png"
    changes_path = figure_dir / "echo_shortcut_classification_changes.png"
    for path, key, ylabel in (
        (
            agreement_path,
            "classification_agreement_fraction",
            "Classification agreement fraction",
        ),
        (changes_path, "shortcut_fraction", "Shortcut classification fraction"),
    ):
        labels, means = _mean_by_type(rows, key)
        fig, ax = plt.subplots(figsize=(8.5, 4.8))
        ax.bar(labels, means)
        ax.set_xlabel("Coarse-graining choice")
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, axis="y", alpha=0.3)
        fig.autofmt_xdate(rotation=25)
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
    return agreement_path, changes_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote echo shortcut classification under coarse-graining: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
