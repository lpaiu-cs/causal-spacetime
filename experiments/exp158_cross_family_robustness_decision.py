"""Apply cross-family robustness criteria to Milestone 33 outputs."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from cross_family_robustness_experiment_helpers import (
    decision_to_row,
    load_metrics_and_decisions,
)
from manifest_family_experiment_helpers import save_bar_figure, write_csv

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for robustness decision experiment."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Cross-family robustness decisions.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    """Aggregate metrics and apply fixed carry-forward decisions."""

    metrics, decisions, missing = load_metrics_and_decisions(config.output_dir)
    decision_rows = [
        decision_to_row(decision, missing_inputs=missing)
        for decision in decisions
    ]
    if missing and not decision_rows:
        decision_rows.append(
            {
                "family_name": "__missing_inputs__",
                "family_kind": "report_only",
                "decision": "blocked",
                "passed": 0.0,
                "failed_reasons": "missing_family_outputs",
                "warning_reasons": "",
                "missing_inputs": ";".join(missing),
            }
        )
    return metrics, decision_rows


def write_outputs(
    metrics: list[dict[str, float | str]],
    decisions: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write robustness metrics and decisions."""

    data_dir = output_dir / "data"
    metrics_path = write_csv(
        metrics,
        data_dir / "cross_family_robustness_metrics.csv",
        [
            "family_name",
            "family_kind",
            "manifest_count",
            "fitted_count",
            "no_fit_count",
            "fitted_fraction",
            "no_fit_fraction",
            "mean_heldout_violation",
            "mean_generalization_gap",
            "stricter_threshold_pass_fraction",
            "destructive_null_gap",
            "symmetry_control_gap",
            "target_coverage_fraction",
            "pair_node_coverage_fraction",
            "no_retuning_audit_pass",
            "failed_accounting_present",
        ],
    )
    decisions_path = write_csv(
        decisions,
        data_dir / "cross_family_robustness_decisions.csv",
        [
            "family_name",
            "family_kind",
            "decision",
            "passed",
            "failed_reasons",
            "warning_reasons",
            "missing_inputs",
        ],
    )
    return metrics_path, decisions_path


def save_figures(
    decisions: list[dict[str, float | str]],
    output_dir: Path,
) -> list[Path]:
    """Save robustness decision summary figures."""

    counts: dict[str, float] = {}
    heldout_by_decision: dict[str, list[float]] = {}
    for row in decisions:
        decision = str(row["decision"])
        counts[decision] = counts.get(decision, 0.0) + 1.0
        try:
            heldout = float(row["mean_heldout_violation"])
        except (KeyError, TypeError, ValueError):
            continue
        if heldout == heldout:
            heldout_by_decision.setdefault(decision, []).append(heldout)
    labels = sorted(counts)
    heldout_labels = sorted(heldout_by_decision)
    return [
        save_bar_figure(
            labels,
            [counts[label] for label in labels],
            output_dir / "figures" / "cross_family_decision_counts.png",
            ylabel="family count",
        ),
        save_bar_figure(
            heldout_labels,
            [
                sum(heldout_by_decision[label]) / len(heldout_by_decision[label])
                for label in heldout_labels
            ],
            output_dir / "figures" / "cross_family_heldout_vs_decision.png",
            ylabel="mean held-out violation",
        ),
    ]


def main() -> None:
    config = parse_args()
    metrics, decisions = run_experiment(config)
    paths = write_outputs(metrics, decisions, config.output_dir)
    figure_paths = save_figures(decisions, config.output_dir)
    print(f"Wrote cross-family robustness decisions: {', '.join(map(str, paths))}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
