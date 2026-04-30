"""Export planned-only protocol-patched v3 preregistration."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_v3_protocol_patch import (
    default_v3_protocol_invariant_family_patches,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_preregistration import (
    build_v3_protocol_patched_preregistration,
    write_v3_protocol_patched_preregistration,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(
        description="Export protocol-patched v3 preregistration."
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[object, list[dict[str, float | str]]]:
    base_spec_path = (
        config.output_dir / "remediation" / "v3_preregistration_spec_m39.json"
    )
    patches = default_v3_protocol_invariant_family_patches()
    spec = build_v3_protocol_patched_preregistration(base_spec_path, patches)
    summary = [
        {
            "spec_id": spec.spec_id,
            "created_by_milestone": spec.created_by_milestone,
            "base_spec_present": float(base_spec_path.exists()),
            "patched_family_count": float(len(spec.patched_families)),
            "execution_allowed_in_current_milestone": float(
                spec.execution_allowed_in_current_milestone
            ),
        }
    ]
    return spec, summary


def main() -> None:
    config = parse_args()
    spec, summary = run_experiment(config)
    json_path = write_v3_protocol_patched_preregistration(
        spec,
        config.output_dir
        / "remediation"
        / "v3_protocol_patched_preregistration_m40.json",
    )
    csv_path = write_csv(
        summary,
        data_path(config.output_dir, "v3_protocol_patched_preregistration_export.csv"),
        [
            "spec_id",
            "created_by_milestone",
            "base_spec_present",
            "patched_family_count",
            "execution_allowed_in_current_milestone",
        ],
    )
    print(f"Wrote v3 protocol-patched preregistration export: {csv_path}")
    print(f"Wrote v3 protocol-patched preregistration JSON: {json_path}")


if __name__ == "__main__":
    main()
