"""Export planned-only v3 preregistration JSON."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path, read_csv_rows

from causal_spacetime_lab.state_change_manifest_v3_design import (
    V3ManifestFamilyDesign,
    default_v3_manifest_family_designs,
)
from causal_spacetime_lab.state_change_manifest_v3_preregistration import (
    build_v3_preregistration_spec,
    write_v3_preregistration_spec,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 preregistration export.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def _bool(value: str) -> bool:
    return value.lower() in {"1", "1.0", "true", "yes"}


def _designs_from_csv(output_dir: Path) -> list[V3ManifestFamilyDesign]:
    rows = read_csv_rows(data_path(output_dir, "v3_manifest_family_design.csv"))
    if not rows:
        return default_v3_manifest_family_designs()
    return [
        V3ManifestFamilyDesign(
            family_name=row["family_name"],
            family_kind=row["family_kind"],
            planned_manifest_count=int(float(row["planned_manifest_count"])),
            profile_column_policy=row["profile_column_policy"],
            target_inclusion_policy=row["target_inclusion_policy"],
            rank_resolution_policy=row["rank_resolution_policy"],
            coverage_policy=row["coverage_policy"],
            null_taxonomy_policy=row["null_taxonomy_policy"],
            restart_stability_policy=row["restart_stability_policy"],
            latent_order_stability_policy=row["latent_order_stability_policy"],
            comparison_method=row["comparison_method"],
            min_margin=float(row["min_margin"]),
            execution_status=row["execution_status"],
        )
        for row in rows
        if row
    ]


def run_experiment(
    config: ExperimentConfig,
) -> tuple[object, list[dict[str, float | str]]]:
    designs = _designs_from_csv(config.output_dir)
    spec = build_v3_preregistration_spec(designs)
    summary = [
        {
            "spec_id": spec.spec_id,
            "created_by_milestone": spec.created_by_milestone,
            "planned_family_count": float(len(spec.planned_families)),
            "execution_allowed_in_current_milestone": float(
                spec.execution_allowed_in_current_milestone
            ),
        }
    ]
    return spec, summary


def main() -> None:
    config = parse_args()
    spec, summary = run_experiment(config)
    json_path = write_v3_preregistration_spec(
        spec,
        config.output_dir / "remediation" / "v3_preregistration_spec_m39.json",
    )
    csv_path = write_csv(
        summary,
        data_path(config.output_dir, "v3_preregistration_export.csv"),
        [
            "spec_id",
            "created_by_milestone",
            "planned_family_count",
            "execution_allowed_in_current_milestone",
        ],
    )
    print(f"Wrote v3 preregistration export: {csv_path}")
    print(f"Wrote v3 preregistration JSON: {json_path}")


if __name__ == "__main__":
    main()
