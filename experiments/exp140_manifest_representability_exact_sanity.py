"""Exact sanity checks for frozen-manifest representation diagnostics."""

from __future__ import annotations

import csv
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


def run_experiment(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, float | str]]:
    """Run exact frozen-manifest representation checks."""

    manifest_path = build_exact_manifest(
        output_dir,
        "response_handoff_manifest_representation_exact.json",
    )
    dataset = load_manifest_dataset(manifest_path)
    integrity = manifest_integrity_report(dataset)
    fit = fit_manifest_ordinal_representation(
        dataset,
        ManifestRepresentationConfig(embedding_dim=1, steps=200, restarts=1, seed=0),
    )
    row = representation_fit_to_row(fit)
    checks = [
        (
            "integrity_passed",
            integrity["passed_integrity"] == 1.0,
            integrity["passed_integrity"],
        ),
        (
            "train_constraints_used",
            fit.train_constraint_count == dataset.train_constraints.shape[0]
            and fit.train_constraint_count > 0,
            fit.train_constraint_count,
        ),
        (
            "heldout_constraints_evaluated",
            fit.heldout_constraint_count == dataset.heldout_constraints.shape[0]
            and fit.heldout_constraint_count > 0,
            fit.heldout_constraint_count,
        ),
        (
            "manifest_has_no_embedding_fields",
            integrity["has_embedding_fields"] == 0.0,
            integrity["has_embedding_fields"],
        ),
        (
            "fit_is_latent_representation",
            row["representation_kind"] == "latent_ordinal_representation",
            row["representation_kind"],
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
    """Write exact sanity CSV."""

    path = output_dir / "data" / "manifest_representability_exact_sanity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote manifest representability exact sanity: {output_path}")


if __name__ == "__main__":
    main()
