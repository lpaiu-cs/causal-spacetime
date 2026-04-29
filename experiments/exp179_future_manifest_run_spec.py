"""Export a future manifest-run spec derived from the remediation plan."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from remediation_plan_experiment_helpers import (
    read_or_build_remediation_plan,
    write_table,
)

from causal_spacetime_lab.state_change_manifest_future_run_spec import (
    build_future_manifest_run_spec_from_plan,
    future_run_spec_to_jsonable,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for future manifest-run spec export."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Future manifest-run spec export.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[Path, list[dict[str, float | str]]]:
    """Build and export a future-run spec."""

    plan = read_or_build_remediation_plan(config.output_dir)
    spec = build_future_manifest_run_spec_from_plan(plan)
    path = config.output_dir / "remediation" / "future_manifest_run_spec_m36.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(future_run_spec_to_jsonable(spec), indent=2, sort_keys=True)
    )
    row = asdict(spec)
    row["allowed_to_execute_now"] = float(spec.allowed_to_execute_now)
    row["planned_manifest_families"] = ";".join(spec.planned_manifest_families)
    row["required_output_files"] = ";".join(spec.required_output_files)
    row["required_metrics"] = ";".join(spec.required_metrics)
    row["forbidden_actions"] = ";".join(spec.forbidden_actions)
    row["json_path"] = str(path)
    return path, [row]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write future-run spec summary."""

    return write_table(
        rows,
        output_dir,
        "future_manifest_run_spec.csv",
        [
            "run_name",
            "run_status",
            "allowed_to_execute_now",
            "prerequisite_milestone",
            "planned_manifest_families",
            "required_output_files",
            "required_metrics",
            "fixed_threshold_source",
            "forbidden_actions",
            "json_path",
        ],
    )


def main() -> None:
    config = parse_args()
    spec_path, rows = run_experiment(config)
    summary_path = write_outputs(rows, config.output_dir)
    print(f"Wrote future manifest-run spec: {summary_path} and {spec_path}")


if __name__ == "__main__":
    main()
