"""Audit existing v2 manifests and planned v3 design for protocol metadata."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path, read_csv_rows

from causal_spacetime_lab.state_change_profile_protocol_invariance_audit import (
    ProfileProtocolInvarianceAuditRow,
    audit_manifest_directory_protocol_invariance,
    audit_row_to_dict,
    audit_v3_design_protocol_invariance,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(
        description="Current response-profile protocol-invariance audit."
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> list[dict[str, float | str]]:
    rows = audit_manifest_directory_protocol_invariance(
        config.output_dir / "manifests_v2",
        source="v2_manifest",
    )
    design_rows = read_csv_rows(
        data_path(config.output_dir, "v3_manifest_family_design.csv")
    )
    rows.extend(audit_v3_design_protocol_invariance(design_rows))
    _write_markdown_report(rows)
    return [audit_row_to_dict(row) for row in rows]


def _write_markdown_report(rows: list[ProfileProtocolInvarianceAuditRow]) -> Path:
    path = Path("docs/validation/response_profile_protocol_invariance_audit.md")
    path.parent.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}
    for row in rows:
        counts[row.profile_invariance_status] = (
            counts.get(row.profile_invariance_status, 0) + 1
        )
    lines = [
        "# Response Profile Protocol Invariance Audit",
        "",
        "Milestone 40 audits response-profile protocol invariance before v3 execution.",
        "This is a pre-execution design correction, not threshold retuning.",
        "It does not change M34, M38, or M39 decisions.",
        "It does not generate v3 manifests, run stress tests, or fit new "
        "representation models.",
        "",
        "## Status Counts",
        "",
    ]
    for status in sorted(counts):
        lines.append(f"- {status}: {counts[status]}")
    lines.extend(["", "## Audit Rows", ""])
    for row in rows:
        lines.append(
            f"- {row.source}: {row.family_id} -> "
            f"{row.profile_invariance_status} "
            f"({row.reason_if_blocked or 'admissible'})"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_csv(
        rows,
        data_path(config.output_dir, "response_profile_protocol_invariance_audit.csv"),
        [
            "family_id",
            "manifest_id",
            "measurement_protocol_ids_used",
            "reference_set_id",
            "number_of_protocols_mixed",
            "number_of_references",
            "profile_invariance_status",
            "admissible_for_pairwise_dissimilarity",
            "reason_if_blocked",
            "source",
        ],
    )
    # Normalize field order for empty-output compatibility.
    if rows:
        with path.open(newline="", encoding="utf-8") as csv_file:
            list(csv.DictReader(csv_file))
    print(f"Wrote current profile protocol-invariance audit: {path}")


if __name__ == "__main__":
    main()
