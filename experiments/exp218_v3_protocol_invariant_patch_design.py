"""Export planned-only protocol-invariant v3 family patches."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_v3_protocol_patch import (
    default_v3_protocol_invariant_family_patches,
    v3_protocol_patch_table,
    write_v3_protocol_patch_json,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 protocol patch design.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> list[dict[str, float | str]]:
    patches = default_v3_protocol_invariant_family_patches()
    write_v3_protocol_patch_json(
        patches,
        config.output_dir
        / "remediation"
        / "v3_protocol_invariant_family_patch_m40.json",
    )
    return v3_protocol_patch_table(patches)


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_csv(
        rows,
        data_path(config.output_dir, "v3_protocol_invariant_family_patch.csv"),
        [
            "original_family_name",
            "patched_family_name",
            "family_kind",
            "measurement_protocol_id",
            "measurement_protocol_hash",
            "echo_rule",
            "emission_rule",
            "gate_rule",
            "reference_subsampling_rule",
            "spectrum_type",
            "normalization_rule",
            "missing_policy",
            "tie_policy",
            "margin_policy",
            "planned_manifest_count",
            "profile_family_semantics",
            "admissible_for_pairwise_dissimilarity",
            "execution_status",
        ],
    )
    print(f"Wrote v3 protocol-invariant family patch: {path}")


if __name__ == "__main__":
    main()
