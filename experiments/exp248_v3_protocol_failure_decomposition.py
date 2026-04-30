"""Decompose v3 protocol carry-forward failures."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import save_metric_bar
from v3_protocol_carry_forward_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_carry_forward import (
    v3_protocol_metrics_rows_from_bundle,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_failure_decomposition import (  # noqa: E501
    decompose_v3_protocol_family_failures,
    summarize_v3_protocol_failure_records,
    v3_protocol_failure_record_rows,
    v3_protocol_missing_or_failed_precondition_rows,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_outputs import (
    load_v3_protocol_output_bundle,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_preconditions import (
    evaluate_v3_protocol_preconditions,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")
    manifest_dir: Path = Path("outputs/manifests_v3")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 protocol failure decomposition.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument(
        "--manifest-dir", type=Path, default=Path("outputs/manifests_v3")
    )
    args = parser.parse_args()
    return ExperimentConfig(args.output_dir, args.manifest_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[
    list[dict[str, float | str]],
    list[dict[str, float | str]],
    list[dict[str, float | str]],
]:
    bundle = load_v3_protocol_output_bundle(config.output_dir)
    preconditions = evaluate_v3_protocol_preconditions(config.manifest_dir, bundle)
    records = decompose_v3_protocol_family_failures(
        v3_protocol_metrics_rows_from_bundle(bundle),
        default_cross_family_robustness_criteria(),
        preconditions,
    )
    return (
        v3_protocol_failure_record_rows(records),
        summarize_v3_protocol_failure_records(records),
        v3_protocol_missing_or_failed_precondition_rows(preconditions),
    )


def main() -> None:
    config = parse_args()
    rows, summary, precondition_rows = run_experiment(config)
    paths = [
        write_csv(
            rows,
            data_path(
                config.output_dir, "v3_protocol_carry_forward_failure_decomposition.csv"
            ),
            ["family_name", "criterion_name", "status"],
        ),
        write_csv(
            summary,
            data_path(
                config.output_dir, "v3_protocol_carry_forward_failure_summary.csv"
            ),
            ["family_name", "root_cause_category", "status", "count"],
        ),
        write_csv(
            precondition_rows,
            data_path(config.output_dir, "v3_protocol_precondition_failure_rows.csv"),
            ["family_name", "precondition", "status"],
        ),
    ]
    root_rows: dict[str, float] = {}
    for row in summary:
        if row["status"] == "measured_failure":
            root = str(row["root_cause_category"])
            root_rows[root] = root_rows.get(root, 0.0) + float(row["count"])
    figure = save_metric_bar(
        [
            {"root_cause_category": root, "count": count}
            for root, count in sorted(root_rows.items())
        ],
        label_key="root_cause_category",
        value_key="count",
        path=config.output_dir / "figures" / "v3_protocol_failure_root_causes.png",
        ylabel="failure count",
    )
    print("Wrote v3 protocol failure decomposition: " + ", ".join(map(str, paths)))
    print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
