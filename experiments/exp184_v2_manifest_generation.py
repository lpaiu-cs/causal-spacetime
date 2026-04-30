"""Generate production v2 handoff manifests from the M36 plan."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import (
    DEFAULT_OUTPUT_DIR,
    data_path,
    load_v2_specs_with_fallback,
    missing_input_rows,
    read_json_object,
    remediation_path,
)

from causal_spacetime_lab.state_change_manifest_v2_generation import (
    V2ManifestGenerationConfig,
    build_v2_handoff_manifests,
    write_v2_handoff_manifests,
)
from causal_spacetime_lab.state_change_manifest_v2_profiles import (
    default_v2_profile_configs,
)
from causal_spacetime_lab.state_change_manifest_v2_spec import (
    mark_v2_family_specs_executed,
)


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 manifest generation."""

    output_dir: Path = DEFAULT_OUTPUT_DIR
    seed: int = 0
    max_constraints: int = 5000
    min_margin: float = 0.05
    bootstrap_count: int = 20
    null_repetitions: int = 5


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="Generate v2 handoff manifests.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-constraints", type=int, default=5000)
    parser.add_argument("--min-margin", type=float, default=0.05)
    parser.add_argument("--bootstrap-count", type=int, default=20)
    parser.add_argument("--null-repetitions", type=int, default=5)
    args = parser.parse_args()
    return ExperimentConfig(
        output_dir=args.output_dir,
        seed=args.seed,
        max_constraints=args.max_constraints,
        min_margin=args.min_margin,
        bootstrap_count=args.bootstrap_count,
        null_repetitions=args.null_repetitions,
    )


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[Path], list[dict[str, float | str]]]:
    """Build and write v2 handoff manifests."""

    plan_path = remediation_path(config.output_dir, "remediation_plan_m36.json")
    spec_path = remediation_path(config.output_dir, "future_manifest_run_spec_m36.json")
    missing = missing_input_rows([plan_path, spec_path])
    if missing:
        rows = [
            {
                "row_type": "missing_input",
                "family_name": "",
                "manifest_id": "",
                "eligible": 0.0,
                "failed_reasons": str(row["input_path"]),
                "constraint_count": 0.0,
                "target_coverage_fraction": float("nan"),
                "pair_node_coverage_fraction": float("nan"),
                "manifest_path": "",
                "execution_status": "not_executed",
            }
            for row in missing
        ]
        return [], rows

    future_spec = read_json_object(spec_path)
    specs = mark_v2_family_specs_executed(
        load_v2_specs_with_fallback(config.output_dir)
    )
    profile_configs = default_v2_profile_configs(specs, seed=config.seed)
    generation_config = V2ManifestGenerationConfig(
        max_constraints=config.max_constraints,
        min_margin=config.min_margin,
        bootstrap_count=config.bootstrap_count,
        null_repetitions=config.null_repetitions,
        constraint_seed=config.seed,
    )
    manifests = build_v2_handoff_manifests(specs, profile_configs, generation_config)
    paths = write_v2_handoff_manifests(manifests, config.output_dir / "manifests_v2")
    rows: list[dict[str, float | str]] = []
    for manifest, path in zip(manifests, paths, strict=True):
        summary = manifest.validation_summary
        rows.append(
            {
                "row_type": "manifest",
                "family_name": manifest.profile_label,
                "manifest_id": manifest.manifest_id,
                "eligible": float(manifest.handoff_decision.eligible),
                "failed_reasons": ";".join(manifest.handoff_decision.failed_reasons),
                "constraint_count": float(summary.constraint_count),
                "target_coverage_fraction": float(summary.target_coverage_fraction),
                "pair_node_coverage_fraction": float(
                    summary.pair_node_coverage_fraction
                ),
                "manifest_path": str(path),
                "execution_status": "executed",
                "future_run_status": str(
                    future_spec.get("run_status", "requires_new_preregistration")
                ),
            }
        )
    return paths, rows


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write generation summary."""

    return write_csv(
        rows,
        data_path(output_dir, "v2_manifest_generation.csv"),
        [
            "row_type",
            "family_name",
            "manifest_id",
            "eligible",
            "failed_reasons",
            "constraint_count",
            "target_coverage_fraction",
            "pair_node_coverage_fraction",
            "manifest_path",
            "execution_status",
            "future_run_status",
        ],
    )


def main() -> None:
    config = parse_args()
    paths, rows = run_experiment(config)
    output_path = write_outputs(rows, config.output_dir)
    print(f"Wrote v2 manifest generation summary: {output_path}")
    print(f"Wrote {len(paths)} v2 manifests")


if __name__ == "__main__":
    main()
