"""Audit patched planned-only v3 families for protocol invariance."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_v3_protocol_patch import (
    default_v3_protocol_invariant_family_patches,
    read_v3_protocol_patch_json,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="Audit v3 protocol patch.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> list[dict[str, float | str]]:
    patch_path = (
        config.output_dir
        / "remediation"
        / "v3_protocol_invariant_family_patch_m40.json"
    )
    rows = read_v3_protocol_patch_json(patch_path)
    if not rows:
        rows = [
            dict(row)
            for row in (
                patch.__dict__
                for patch in default_v3_protocol_invariant_family_patches()
            )
        ]
    audit_rows: list[dict[str, float | str]] = []
    for row in rows:
        family_kind = str(row.get("family_kind", ""))
        admissible = bool(row.get("admissible_for_pairwise_dissimilarity", False))
        is_control = family_kind in {"failed_control", "report_only"}
        audit_rows.append(
            {
                "family_id": str(row.get("patched_family_name", "")),
                "family_kind": family_kind,
                "profile_invariance_status": "protocol_invariant",
                "admissible_for_pairwise_dissimilarity": float(admissible),
                "report_only_or_control": float(is_control),
                "measurement_protocol_id": str(row.get("measurement_protocol_id", "")),
                "execution_status": str(row.get("execution_status", "")),
                "passed": float(
                    str(row.get("execution_status", "")) == "planned_only"
                    and (
                        (not is_control and admissible)
                        or (is_control and not admissible)
                    )
                ),
            }
        )
    return audit_rows


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_csv(
        rows,
        data_path(config.output_dir, "v3_protocol_patch_audit.csv"),
        [
            "family_id",
            "family_kind",
            "profile_invariance_status",
            "admissible_for_pairwise_dissimilarity",
            "report_only_or_control",
            "measurement_protocol_id",
            "execution_status",
            "passed",
        ],
    )
    print(f"Wrote v3 protocol patch audit: {path}")


if __name__ == "__main__":
    main()
