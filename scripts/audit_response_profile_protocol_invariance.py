"""Audit response-profile protocol invariance without executing v3."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from causal_spacetime_lab.state_change_profile_protocol_invariance_audit import (
    ProfileProtocolInvarianceAuditRow,
    audit_manifest_directory_protocol_invariance,
    audit_row_to_dict,
    audit_v3_design_protocol_invariance,
)


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    """Read CSV rows, returning no rows when absent."""

    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as csv_file:
        return list(csv.DictReader(csv_file))


def write_csv_rows(rows: list[dict[str, float | str]], path: Path) -> Path:
    """Write audit rows to CSV."""

    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
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
    ]
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


def audit_current_outputs(output_dir: Path) -> list[ProfileProtocolInvarianceAuditRow]:
    """Audit v2 manifests, planned v3 design, and optional patched M40 spec."""

    rows = audit_manifest_directory_protocol_invariance(
        output_dir / "manifests_v2",
        source="v2_manifest",
    )
    design_rows = read_csv_rows(output_dir / "data" / "v3_manifest_family_design.csv")
    rows.extend(audit_v3_design_protocol_invariance(design_rows))
    patch_path = (
        output_dir
        / "remediation"
        / "v3_protocol_patched_preregistration_m40.json"
    )
    if patch_path.exists():
        rows.extend(_audit_patched_preregistration(patch_path))
    return rows


def write_markdown_report(
    rows: list[ProfileProtocolInvarianceAuditRow],
    path: Path,
) -> Path:
    """Write a compact validation report."""

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
        "This report does not modify v2 decisions, generate v3 manifests, "
        "run fit diagnostics, or run carry-forward evaluation.",
        "",
        "## Status Counts",
        "",
    ]
    for status in sorted(counts):
        lines.append(f"- {status}: {counts[status]}")
    lines.extend(["", "## Rows", ""])
    for row in rows:
        lines.append(
            "- "
            f"{row.source} / {row.family_id}: "
            f"{row.profile_invariance_status}; "
            f"admissible={row.admissible_for_pairwise_dissimilarity}; "
            f"reason={row.reason_if_blocked}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def run_audit(output_dir: Path) -> list[dict[str, float | str]]:
    """Run the current-output audit and write validation artifacts."""

    rows = audit_current_outputs(output_dir)
    row_dicts = [audit_row_to_dict(row) for row in rows]
    write_csv_rows(
        row_dicts,
        output_dir / "data" / "response_profile_protocol_invariance_audit.csv",
    )
    write_markdown_report(
        rows,
        Path("docs/validation/response_profile_protocol_invariance_audit.md"),
    )
    return row_dicts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit response-profile protocol invariance."
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = run_audit(args.output_dir)
    print(
        "Wrote response-profile protocol invariance audit: "
        f"{args.output_dir / 'data' / 'response_profile_protocol_invariance_audit.csv'}"
    )
    print(f"Audit rows: {len(rows)}")


def _audit_patched_preregistration(
    path: Path,
) -> list[ProfileProtocolInvarianceAuditRow]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    patched = payload.get("patched_families", [])
    rows: list[ProfileProtocolInvarianceAuditRow] = []
    if not isinstance(patched, list):
        return rows
    for item in patched:
        if not isinstance(item, dict):
            continue
        family_id = str(item.get("patched_family_name", ""))
        admissible = bool(item.get("admissible_for_pairwise_dissimilarity", False))
        family_kind = str(item.get("family_kind", ""))
        status = "protocol_invariant"
        reason = ""
        if family_kind in {"failed_control", "report_only"} or not admissible:
            reason = "explicit report-only or failed-control family"
        rows.append(
            ProfileProtocolInvarianceAuditRow(
                family_id=family_id,
                manifest_id="",
                measurement_protocol_ids_used=[
                    str(item.get("measurement_protocol_id", ""))
                ],
                reference_set_id="planned_v3_reference_set",
                number_of_protocols_mixed=1,
                number_of_references=0,
                profile_invariance_status=status,
                admissible_for_pairwise_dissimilarity=admissible,
                reason_if_blocked=reason,
                source="v3_protocol_patch",
            )
        )
    return rows


if __name__ == "__main__":
    main()
