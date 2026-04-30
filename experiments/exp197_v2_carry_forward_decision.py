"""Apply fixed M34 carry-forward criteria to v2 metrics."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import save_metric_bar

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v2_carry_forward import (
    decide_v2_family_robustness,
    v2_decision_to_row,
    v2_diagnostic_completeness_by_family,
    v2_metrics_rows_from_bundle,
)
from causal_spacetime_lab.state_change_manifest_v2_outputs import (
    load_v2_output_bundle,
    missing_v2_bundle_inputs,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 carry-forward decisions."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 carry-forward decision.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    """Apply fixed criteria to v2 metric rows."""

    bundle = load_v2_output_bundle(config.output_dir)
    missing = missing_v2_bundle_inputs(bundle)
    completeness = v2_diagnostic_completeness_by_family(bundle)
    metrics = v2_metrics_rows_from_bundle(bundle)
    decisions = decide_v2_family_robustness(
        metrics,
        default_cross_family_robustness_criteria(),
    )
    decision_rows = [
        v2_decision_to_row(
            decision,
            missing_inputs=missing,
            diagnostic_complete=completeness.get(decision.family_name, False),
        )
        for decision in decisions
    ]
    return metrics, decision_rows


def write_outputs(
    metrics: list[dict[str, float | str]],
    decisions: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write v2 decision metrics and decisions."""

    return (
        write_csv(
            metrics,
            output_dir / "data" / "v2_cross_family_robustness_decision_metrics.csv",
            [
                "family_name",
                "family_kind",
                "manifest_count",
                "mean_heldout_violation",
            ],
        ),
        write_csv(
            decisions,
            output_dir / "data" / "v2_cross_family_robustness_decisions.csv",
            [
                "family_name",
                "family_kind",
                "decision",
                "passed",
                "failed_reasons",
                "warning_reasons",
                "missing_inputs",
                "diagnostic_complete",
            ],
        ),
    )


def save_figures(
    decisions: list[dict[str, float | str]],
    output_dir: Path,
) -> list[Path]:
    """Save v2 carry-forward decision figures."""

    counts: dict[str, float] = {}
    heldout_rows: list[dict[str, float | str]] = []
    for row in decisions:
        decision = str(row["decision"])
        counts[decision] = counts.get(decision, 0.0) + 1.0
        heldout_rows.append(
            {
                "family_name": f"{row['family_name']}:{decision}",
                "mean_heldout_violation": row.get("mean_heldout_violation", "nan"),
            }
        )
    count_rows = [
        {"decision": decision, "count": count}
        for decision, count in sorted(counts.items())
    ]
    return [
        save_metric_bar(
            count_rows,
            label_key="decision",
            value_key="count",
            path=output_dir / "figures" / "v2_cross_family_decision_counts.png",
            ylabel="family count",
        ),
        save_metric_bar(
            heldout_rows,
            label_key="family_name",
            value_key="mean_heldout_violation",
            path=output_dir
            / "figures"
            / "v2_cross_family_heldout_vs_decision.png",
            ylabel="mean held-out violation",
        ),
    ]


def main() -> None:
    config = parse_args()
    metrics, decisions = run_experiment(config)
    paths = write_outputs(metrics, decisions, config.output_dir)
    figure_paths = save_figures(decisions, config.output_dir)
    print("Wrote v2 carry-forward decisions: " + ", ".join(str(path) for path in paths))
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
