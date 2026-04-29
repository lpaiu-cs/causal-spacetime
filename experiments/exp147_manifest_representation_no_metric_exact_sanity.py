"""Exact no-metric checks for manifest representation diagnostics."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from manifest_representation_experiment_helpers import build_exact_manifest

from causal_spacetime_lab.state_change_manifest_dataset import (
    load_manifest_dataset,
    manifest_integrity_report,
)
from causal_spacetime_lab.state_change_manifest_representation import (
    ManifestRepresentationConfig,
    fit_manifest_ordinal_representation,
    representation_fit_to_row,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def _failed_manifest_from(path: Path, output_dir: Path) -> Path:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest["manifest_id"] = "manifest_representation_ineligible_exact"
    manifest["handoff_decision"]["eligible"] = False
    manifest["handoff_decision"]["failed_reasons"] = ["exact_ineligible_control"]
    failed_path = output_dir / "manifests" / "response_handoff_no_metric_failed.json"
    failed_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return failed_path


def run_experiment(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, float | str]]:
    """Run exact no-metric representation checks."""

    manifest_path = build_exact_manifest(
        output_dir,
        "response_handoff_no_metric_exact.json",
    )
    dataset = load_manifest_dataset(manifest_path)
    integrity = manifest_integrity_report(dataset)
    fit = fit_manifest_ordinal_representation(
        dataset,
        ManifestRepresentationConfig(embedding_dim=1, steps=100, restarts=1),
    )
    row = representation_fit_to_row(fit)
    failed_path = _failed_manifest_from(manifest_path, output_dir)
    failed_dataset = load_manifest_dataset(failed_path)
    failed_fit = fit_manifest_ordinal_representation(
        failed_dataset,
        ManifestRepresentationConfig(embedding_dim=1, steps=100, restarts=1),
    )
    train_indices = set(dataset.train_constraint_indices.tolist())
    heldout_indices = set(dataset.heldout_constraint_indices.tolist())
    checks = [
        (
            "loaded_manifest_has_no_metric_fields",
            integrity["has_metric_fields"] == 0.0,
            integrity["has_metric_fields"],
        ),
        (
            "fit_row_excludes_embedding_by_default",
            "embedding_values" not in row,
            sorted(row),
        ),
        (
            "embedding_internal_only",
            fit.embedding is not None and row["representation_kind"]
            == "latent_ordinal_representation",
            fit.embedding.shape if fit.embedding is not None else None,
        ),
        (
            "forbidden_interpretations_include_metric_terms",
            any("metric" in item for item in dataset.forbidden_interpretations),
            dataset.forbidden_interpretations,
        ),
        (
            "ineligible_skipped_by_default",
            not failed_fit.fitted
            and failed_fit.reason_not_fit == "manifest_ineligible",
            failed_fit.reason_not_fit,
        ),
        (
            "train_heldout_split_respected",
            not train_indices & heldout_indices
            and fit.train_constraint_count == dataset.train_constraints.shape[0],
            f"{fit.train_constraint_count}/{fit.heldout_constraint_count}",
        ),
    ]
    return [
        {"check": name, "passed": float(passed), "value": str(value)}
        for name, passed, value in checks
    ]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write exact no-metric sanity CSV."""

    path = output_dir / "data" / "manifest_representation_no_metric_exact_sanity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote manifest representation no-metric exact sanity: {output_path}")


if __name__ == "__main__":
    main()
