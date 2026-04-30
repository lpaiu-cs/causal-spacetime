"""Report diagnostic completeness for the patched v3 output bundle."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path, read_csv_rows


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 protocol bundle report.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return ExperimentConfig(args.output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    rows = read_csv_rows(
        data_path(config.output_dir, "v3_protocol_diagnostic_completeness_check.csv")
    )
    report: list[dict[str, float | str]] = []
    for row in rows:
        report.append(
            {
                "family_name": row.get("family_name", ""),
                "diagnostic_complete": float(row.get("diagnostic_complete", 0.0)),
                "missing_metrics": row.get("missing_metrics", ""),
                "not_carry_forward_evaluated": 1.0,
                "profile_invariance_status": row.get(
                    "profile_invariance_status",
                    "",
                ),
                "admissible_for_pairwise_dissimilarity": float(
                    row.get("admissible_for_pairwise_dissimilarity", 0.0)
                ),
                "dominant_handoff_provenance_type": row.get(
                    "dominant_handoff_provenance_type",
                    "",
                ),
                "all_manifests_have_provenance": float(
                    row.get("all_manifests_have_provenance", 0.0)
                ),
                "interpretation_warning": (
                    "diagnostic-complete v3 output bundle is not "
                    "carry-forward evaluated"
                ),
            }
        )
    return report


def main() -> None:
    config = parse_args()
    path = write_csv(
        run_experiment(config),
        data_path(
            config.output_dir, "v3_protocol_diagnostic_complete_bundle_report.csv"
        ),
        ["family_name", "diagnostic_complete", "missing_metrics"],
    )
    print(f"Wrote v3 protocol diagnostic-complete bundle report: {path}")


if __name__ == "__main__":
    main()
