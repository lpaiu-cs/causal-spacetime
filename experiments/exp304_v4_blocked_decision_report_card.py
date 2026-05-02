"""Aggregate M46 v4 blocked-decision audit outputs."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v4_blocked_v5_experiment_helpers import data_path, read_csv_rows


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V4 blocked report card.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def _dominant_root(summary_rows: list[dict[str, str]], family: str) -> str:
    counts: dict[str, float] = {}
    for row in summary_rows:
        if row.get("family_name") == family and row.get("status") == "fail":
            root = row.get("root_cause_category", "")
            counts[root] = counts.get(root, 0.0) + float(row.get("count", "0") or 0)
    if not counts:
        return "none"
    return max(counts, key=counts.get)


def _linked_v5(root: str) -> str:
    mapping = {
        "heldout_failure": "rank_gap_earliest_full_manifest_transfer_v5",
        "latent_order_instability": "rank_gap_latent_stability_replicated_v5",
        "null_separation_failure": "combined_null_separation_v5",
        "stricter_pass_failure": "rank_gap_low_tie_reference_diverse_v5",
        "coverage_failure": "rank_gap_high_coverage_strict_pair_v5",
    }
    return mapping.get(root, "v5_protocol_family_design")


def _delta_summary(delta_rows: list[dict[str, str]], family: str) -> str:
    matches = [row for row in delta_rows if row.get("v4_family_name") == family]
    if not matches:
        return ""
    row = matches[0]
    return (
        f"improved={row.get('improved_metric_count', '')};"
        f"worsened={row.get('worsened_metric_count', '')};"
        f"dominant_improvement={row.get('dominant_improvement', '')};"
        f"dominant_regression={row.get('dominant_regression', '')}"
    )


def _risk_warning(risk_rows: list[dict[str, str]], linked_design: str) -> str:
    for row in risk_rows:
        if row.get("proposed_action") == linked_design:
            return f"{row.get('risk_level', '')}:{row.get('risk_category', '')}"
    return ""


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    decisions = read_csv_rows(
        data_path(
            config.output_dir,
            "v4_protocol_cross_family_robustness_decisions.csv",
        )
    )
    margins = {
        row.get("family_name", ""): row
        for row in read_csv_rows(
            data_path(config.output_dir, "v4_criterion_margin_summary.csv")
        )
    }
    summary_rows = read_csv_rows(
        data_path(config.output_dir, "v4_blocked_root_cause_summary.csv")
    )
    counterfactual_rows = read_csv_rows(
        data_path(config.output_dir, "v4_single_fix_counterfactual.csv")
    )
    delta_rows = read_csv_rows(
        data_path(config.output_dir, "v3_to_v4_metric_delta_summary.csv")
    )
    risk_rows = read_csv_rows(
        data_path(config.output_dir, "v5_remediation_iteration_risk_audit.csv")
    )
    rows: list[dict[str, float | str]] = []
    for decision in decisions:
        family = decision.get("family_name", "")
        root = _dominant_root(summary_rows, family)
        margin = margins.get(family, {})
        remaining = [
            row
            for row in counterfactual_rows
            if row.get("family_name") == family
            and row.get("ignored_root_cause") == root
        ]
        linked = _linked_v5(root)
        rows.append(
            {
                "family_name": family,
                "decision": decision.get("decision", ""),
                "dominant_root_cause": root,
                "worst_criterion": margin.get("worst_criterion", ""),
                "worst_margin": margin.get("worst_margin", ""),
                "measured_blocker_count": margin.get("measured_failure_count", ""),
                "counterfactual_remaining_blocked": (
                    remaining[0].get("would_remain_blocked", "") if remaining else ""
                ),
                "v3_to_v4_improvement_summary": _delta_summary(delta_rows, family),
                "linked_v5_design": linked,
                "iteration_risk_warning": _risk_warning(risk_rows, linked),
                "interpretation_warning": (
                    "blocked-decision audit is not metric geometry"
                ),
            }
        )
    return rows


def main() -> None:
    config = parse_args()
    path = write_csv(
        run_experiment(config),
        data_path(config.output_dir, "v4_blocked_decision_report_card.csv"),
        ["family_name", "decision", "dominant_root_cause"],
    )
    print(f"Wrote v4 blocked-decision report card: {path}")


if __name__ == "__main__":
    main()
