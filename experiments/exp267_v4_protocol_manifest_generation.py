"""Generate preregistered v4 protocol handoff manifests."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_v4_execution_spec import (
    load_v4_protocol_execution_specs,
    mark_v4_protocol_execution_specs_executed,
    v4_protocol_execution_spec_table,
)
from causal_spacetime_lab.state_change_manifest_v4_generation import (
    V4ProtocolManifestGenerationConfig,
    build_v4_protocol_handoff_manifests,
    write_v4_protocol_handoff_manifests,
)
from causal_spacetime_lab.state_change_manifest_v4_profiles import (
    default_v4_protocol_profile_configs,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")
    seed: int = 0
    max_constraints: int = 1000
    bootstrap_count: int = 5
    null_repetitions: int = 3


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="Generate v4 protocol manifests.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-constraints", type=int, default=1000)
    parser.add_argument("--bootstrap-count", type=int, default=5)
    parser.add_argument("--null-repetitions", type=int, default=3)
    args = parser.parse_args()
    return ExperimentConfig(
        output_dir=args.output_dir,
        seed=args.seed,
        max_constraints=args.max_constraints,
        bootstrap_count=args.bootstrap_count,
        null_repetitions=args.null_repetitions,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    spec_path = (
        config.output_dir
        / "remediation"
        / "v4_protocol_preregistration_spec_m43.json"
    )
    specs = load_v4_protocol_execution_specs(spec_path)
    if not specs:
        return [
            {
                "status": "missing_input",
                "missing_input": str(spec_path),
                "manifest_count": 0.0,
            }
        ]
    profile_configs = default_v4_protocol_profile_configs(specs, seed=config.seed)
    manifests = build_v4_protocol_handoff_manifests(
        specs,
        profile_configs,
        V4ProtocolManifestGenerationConfig(
            max_constraints=config.max_constraints,
            bootstrap_count=config.bootstrap_count,
            null_repetitions=config.null_repetitions,
            constraint_seed=config.seed,
        ),
    )
    paths = write_v4_protocol_handoff_manifests(
        manifests,
        config.output_dir / "manifests_v4",
        overwrite=True,
    )
    rows = v4_protocol_execution_spec_table(
        mark_v4_protocol_execution_specs_executed(specs)
    )
    for row in rows:
        family = str(row["family_name"])
        row["status"] = "generated"
        row["manifest_count"] = float(
            sum(manifest.profile_label.startswith(family) for manifest in manifests)
        )
        row["written_path_count"] = float(len(paths))
        row["not_carry_forward_evaluated"] = 1.0
    return rows


def main() -> None:
    config = parse_args()
    path = write_csv(
        run_experiment(config),
        data_path(config.output_dir, "v4_protocol_manifest_generation.csv"),
        ["status", "missing_input", "manifest_count"],
    )
    print(f"Wrote v4 protocol manifest generation: {path}")


if __name__ == "__main__":
    main()
