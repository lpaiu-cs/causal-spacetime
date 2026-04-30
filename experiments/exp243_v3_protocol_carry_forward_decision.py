"""Apply fixed M34 criteria to the M41 v3 protocol bundle."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import save_metric_bar
from v3_protocol_carry_forward_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_carry_forward import (
    decide_v3_protocol_family_robustness,
    v3_protocol_decision_to_row,
    v3_protocol_diagnostic_completeness_by_family,
    v3_protocol_metrics_rows_from_bundle,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_outputs import (
    load_v3_protocol_output_bundle,
    missing_v3_protocol_bundle_inputs,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_preconditions import (
    evaluate_v3_protocol_preconditions,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")
    manifest_dir: Path = Path("outputs/manifests_v3")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 protocol carry-forward decision.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument(
        "--manifest-dir",
        type=Path,
        default=Path("outputs/manifests_v3"),
    )
    args = parser.parse_args()
    return ExperimentConfig(args.output_dir, args.manifest_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    bundle = load_v3_protocol_output_bundle(config.output_dir)
    missing = missing_v3_protocol_bundle_inputs(bundle)
    metrics = v3_protocol_metrics_rows_from_bundle(bundle)
    metrics_by_family = {str(row.get("family_name", "")): row for row in metrics}
    completeness = v3_protocol_diagnostic_completeness_by_family(bundle)
    preconditions = evaluate_v3_protocol_preconditions(config.manifest_dir, bundle)
    precondition_by_family = {item.family_name: item for item in preconditions}
    decisions = decide_v3_protocol_family_robustness(
        metrics,
        default_cross_family_robustness_criteria(),
        preconditions,
    )
    decision_rows = []
    for decision in decisions:
        precondition = precondition_by_family.get(decision.family_name)
        row = v3_protocol_decision_to_row(
            decision,
            diagnostic_complete=completeness.get(decision.family_name, False),
            preconditions_passed=(
                precondition.preconditions_passed if precondition else False
            ),
            failed_preconditions=(
                precondition.failed_preconditions if precondition else []
            ),
            missing_inputs=missing,
        )
        metric_row = metrics_by_family.get(decision.family_name, {})
        for key in [
            "dominant_handoff_provenance_type",
            "top_down_template_count",
            "hybrid_manifest_count",
            "bottom_up_manifest_count",
            "report_only_manifest_count",
            "all_manifests_have_provenance",
        ]:
            if key in metric_row:
                row[key] = metric_row[key]
        decision_rows.append(row)
    return metrics, decision_rows


def write_outputs(
    metrics: list[dict[str, float | str]],
    decisions: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    return (
        write_csv(
            metrics,
            data_path(
                output_dir,
                "v3_protocol_cross_family_robustness_decision_metrics.csv",
            ),
            ["family_name", "family_kind", "manifest_count"],
        ),
        write_csv(
            decisions,
            data_path(output_dir, "v3_protocol_cross_family_robustness_decisions.csv"),
            [
                "family_name",
                "family_kind",
                "decision",
                "passed",
                "failed_reasons",
                "warning_reasons",
                "missing_inputs",
                "diagnostic_complete",
                "preconditions_passed",
                "failed_preconditions",
            ],
        ),
    )


def save_figures(
    decisions: list[dict[str, float | str]],
    output_dir: Path,
) -> list[Path]:
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
            path=output_dir / "figures" / "v3_protocol_decision_counts.png",
            ylabel="family count",
        ),
        save_metric_bar(
            heldout_rows,
            label_key="family_name",
            value_key="mean_heldout_violation",
            path=output_dir / "figures" / "v3_protocol_heldout_vs_decision.png",
            ylabel="mean held-out violation",
        ),
    ]


def main() -> None:
    config = parse_args()
    metrics, decisions = run_experiment(config)
    paths = write_outputs(metrics, decisions, config.output_dir)
    figures = save_figures(decisions, config.output_dir)
    print("Wrote v3 protocol decisions: " + ", ".join(str(path) for path in paths))
    for figure in figures:
        print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
