"""Exact sanity checks for response-comparison handoff manifests."""

from __future__ import annotations

import csv
from pathlib import Path

from response_handoff_experiment_helpers import deterministic_handoff_profile

from causal_spacetime_lab.state_change_response_constraint_validation import (
    ConstraintValidationGate,
)
from causal_spacetime_lab.state_change_response_handoff import (
    manifest_digest,
    manifest_to_jsonable,
    read_handoff_manifest,
    write_handoff_manifest,
)
from causal_spacetime_lab.state_change_response_handoff_pipeline import (
    build_candidate_handoff_manifest,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def run_experiment(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, float | str]]:
    """Run exact handoff manifest sanity checks."""

    manifest = build_candidate_handoff_manifest(
        deterministic_handoff_profile(),
        PairwiseResponseComparisonProtocol("gap", "rank_gap_mean"),
        ConstraintValidationGate("exact", min_constraint_count=1),
        max_constraints=40,
        min_margin=0.0,
        train_fraction=0.6,
        constraint_seed=0,
        bootstrap_count=3,
        null_repetitions=2,
        source_label="exact_sanity",
    )
    jsonable = manifest_to_jsonable(manifest)
    digest = manifest_digest(jsonable)
    manifest_path = output_dir / "manifests" / "response_handoff_exact_sanity.json"
    write_handoff_manifest(manifest, manifest_path)
    loaded = read_handoff_manifest(manifest_path)
    train = set(int(value) for value in manifest.train_constraint_indices)
    heldout = set(int(value) for value in manifest.heldout_constraint_indices)
    checks = [
        (
            "manifest_has_constraints",
            manifest.constraints.shape[0] > 0,
            manifest.constraints.shape[0],
        ),
        ("digest_stable", digest == manifest.manifest_id, digest),
        (
            "json_serialization",
            loaded["manifest_id"] == manifest.manifest_id,
            manifest_path,
        ),
        ("split_disjoint", not train & heldout, f"{len(train)}/{len(heldout)}"),
        (
            "forbidden_interpretations_nonempty",
            bool(manifest.forbidden_interpretations),
            manifest.forbidden_interpretations,
        ),
        (
            "decision_fields_exist",
            isinstance(manifest.handoff_decision.eligible, bool),
            manifest.handoff_decision,
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

    path = output_dir / "data" / "response_handoff_exact_sanity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote response handoff exact sanity: {output_path}")


if __name__ == "__main__":
    main()
