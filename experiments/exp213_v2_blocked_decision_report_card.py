"""Aggregate M39 v2 blocked-decision audit into one report card."""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path, read_csv_rows


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for blocked-decision report card."""

    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V2 blocked-decision report card.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def _counts_by_family(
    rows: list[dict[str, str]],
    *,
    blocking_type: str,
) -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    for row in rows:
        if row.get("blocking_type") == blocking_type and row.get("status") == "fail":
            counts[row["family_name"]] += int(float(row.get("count", "0")))
    return counts


def _summary_lookup(path: Path) -> dict[str, dict[str, str]]:
    return {row.get("family_name", ""): row for row in read_csv_rows(path)}


def _linked_v3_design(family_name: str) -> str:
    if family_name.startswith("rank_gap_more_protocol_columns"):
        return "rank_gap_more_protocol_columns_v3"
    if family_name.startswith("rank_gap_rank_resolution"):
        return "rank_gap_rank_resolution_enriched_v3"
    if family_name.startswith("rank_gap_coverage"):
        return "rank_gap_coverage_enriched_v3"
    if family_name.startswith("combined"):
        return "combined_diagnostic_complete_v3"
    if family_name.startswith("failed"):
        return "failed_controls_v3"
    return "rank_gap_multi_manifest_v3"


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    output_dir = config.output_dir
    decision_rows = read_csv_rows(
        data_path(output_dir, "v2_cross_family_robustness_decisions.csv")
    )
    blocking_summary = read_csv_rows(
        data_path(output_dir, "v2_blocked_root_cause_summary.csv")
    )
    margin_lookup = _summary_lookup(
        data_path(output_dir, "v2_criterion_margin_summary.csv")
    )
    remaining_lookup = _summary_lookup(
        data_path(
            output_dir,
            "v2_remaining_measured_blockers_after_structural_ignore.csv",
        )
    )
    structural_counts = _counts_by_family(
        blocking_summary,
        blocking_type="structural_blocking",
    )
    measured_counts = _counts_by_family(
        blocking_summary,
        blocking_type="measured_blocking",
    )
    diagnostic_counts = _counts_by_family(
        blocking_summary,
        blocking_type="diagnostic_blocking",
    )
    rows: list[dict[str, float | str]] = []
    for decision in decision_rows:
        family = decision.get("family_name", "")
        margin = margin_lookup.get(family, {})
        remaining = remaining_lookup.get(family, {})
        structural = structural_counts.get(family, 0)
        measured = measured_counts.get(family, 0)
        structural_alone = structural > 0 and measured == 0
        rows.append(
            {
                "family_name": family,
                "decision": decision.get("decision", ""),
                "structural_blockers": float(structural),
                "measured_blockers": float(measured),
                "diagnostic_blockers": float(diagnostic_counts.get(family, 0)),
                "worst_criterion": margin.get("worst_criterion", ""),
                "structural_count_alone_explains_block": float(structural_alone),
                "measured_blockers_remain": remaining.get(
                    "would_remain_blocked",
                    float(measured > 0),
                ),
                "linked_v3_design_proposal": _linked_v3_design(family),
                "interpretation_warning": (
                    "v3 design is preregistered but not executed"
                ),
            }
        )
    return rows


def main() -> None:
    config = parse_args()
    path = write_csv(
        run_experiment(config),
        data_path(config.output_dir, "v2_blocked_decision_report_card.csv"),
        [
            "family_name",
            "decision",
            "structural_blockers",
            "measured_blockers",
            "diagnostic_blockers",
            "worst_criterion",
            "structural_count_alone_explains_block",
            "measured_blockers_remain",
            "linked_v3_design_proposal",
            "interpretation_warning",
        ],
    )
    print(f"Wrote v2 blocked-decision report card: {path}")


if __name__ == "__main__":
    main()

