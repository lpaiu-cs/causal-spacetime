"""Aggregate M43 v3 protocol blocked-decision audit outputs."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v3_protocol_blocked_v4_experiment_helpers import data_path, read_csv_rows


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 blocked report card.")
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


def _linked_v4(root: str) -> str:
    mapping = {
        "heldout_failure": "rank_gap_earliest_full_stability_v4",
        "latent_order_instability": "rank_gap_gated_full_stability_v4",
        "null_separation_failure": "combined_earliest_full_stability_v4",
        "coverage_failure": "rank_gap_high_resolution_reference_set_v4",
        "restart_instability": "rank_gap_earliest_retained_resolution_v4",
    }
    return mapping.get(root, "v4_protocol_family_design")


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    decisions = read_csv_rows(
        data_path(
            config.output_dir,
            "v3_protocol_cross_family_robustness_decisions.csv",
        )
    )
    margins = {
        row.get("family_name", ""): row
        for row in read_csv_rows(
            data_path(config.output_dir, "v3_protocol_criterion_margin_summary.csv")
        )
    }
    summary_rows = read_csv_rows(
        data_path(config.output_dir, "v3_protocol_blocked_root_cause_summary.csv")
    )
    counterfactual_rows = read_csv_rows(
        data_path(config.output_dir, "v3_protocol_single_fix_counterfactual.csv")
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
                "linked_v4_design": _linked_v4(root),
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
        data_path(config.output_dir, "v3_protocol_blocked_decision_report_card.csv"),
        ["family_name", "decision", "dominant_root_cause"],
    )
    print(f"Wrote v3 protocol blocked-decision report card: {path}")


if __name__ == "__main__":
    main()
