"""Report failed manifests and verify no-fit behavior by default."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

from manifest_representation_experiment_helpers import build_exact_manifest, write_csv

from causal_spacetime_lab.state_change_manifest_dataset import (
    load_manifest_dataset,
    load_manifest_datasets,
)
from causal_spacetime_lab.state_change_manifest_representation import (
    ManifestRepresentationConfig,
    fit_manifest_ordinal_representation,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for failed-manifest no-fit controls."""

    manifest_dir: Path = Path("outputs/manifests")
    generate_failed_controls: bool = True
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Failed manifest no-fit controls.")
    parser.add_argument("--manifest-dir", type=Path, default=Path("outputs/manifests"))
    parser.add_argument(
        "--generate-failed-controls",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        manifest_dir=args.manifest_dir,
        generate_failed_controls=args.generate_failed_controls,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _write_failed_control(output_dir: Path) -> Path:
    base_path = build_exact_manifest(
        output_dir,
        "response_handoff_failed_control_base.json",
    )
    manifest = json.loads(base_path.read_text(encoding="utf-8"))
    manifest["manifest_id"] = "failed_control_manifest"
    manifest["handoff_decision"]["eligible"] = False
    manifest["handoff_decision"]["failed_reasons"] = [
        "synthetic_failed_control",
    ]
    manifest_path = output_dir / "manifests" / "response_handoff_failed_control.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest_path


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Record explicit no-fit rows for failed manifests."""

    datasets = load_manifest_datasets(config.manifest_dir, include_ineligible=True)
    if (
        config.generate_failed_controls
        and not any(not dataset.eligible for dataset in datasets)
    ):
        failed_path = _write_failed_control(config.output_dir)
        datasets.append(load_manifest_dataset(failed_path))

    rows: list[dict[str, float | str]] = []
    fit_config = ManifestRepresentationConfig(embedding_dim=1, steps=100, seed=0)
    for dataset in datasets:
        fit = fit_manifest_ordinal_representation(dataset, fit_config)
        rows.append(
            {
                "manifest_id": dataset.manifest_id,
                "eligible": float(dataset.eligible),
                "failed_reasons": ";".join(dataset.failed_reasons),
                "fit_attempted": float(fit.fitted),
                "reason_not_fit": fit.reason_not_fit,
                "constraint_count": float(dataset.constraints.shape[0]),
            }
        )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write failed-manifest control CSV."""

    return write_csv(
        rows,
        output_dir / "data" / "failed_manifest_no_fit_controls.csv",
        [
            "manifest_id",
            "eligible",
            "failed_reasons",
            "fit_attempted",
            "reason_not_fit",
            "constraint_count",
        ],
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    print(f"Wrote failed manifest no-fit controls: {data_path}")


if __name__ == "__main__":
    main()
