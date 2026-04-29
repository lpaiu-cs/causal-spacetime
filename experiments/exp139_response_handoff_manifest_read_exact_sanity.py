"""Exact read/write sanity checks for response handoff manifests."""

from __future__ import annotations

import csv
from pathlib import Path

from response_handoff_experiment_helpers import deterministic_handoff_profile

from causal_spacetime_lab.state_change_response_constraint_validation import (
    ConstraintValidationGate,
)
from causal_spacetime_lab.state_change_response_handoff import (
    manifest_digest,
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
    """Run manifest read/write exact sanity checks."""

    manifest = build_candidate_handoff_manifest(
        deterministic_handoff_profile(),
        PairwiseResponseComparisonProtocol("gap", "rank_gap_mean"),
        ConstraintValidationGate("read", min_constraint_count=1),
        max_constraints=40,
        bootstrap_count=3,
        null_repetitions=2,
        source_label="read_sanity",
    )
    path = output_dir / "manifests" / "response_handoff_read_sanity.json"
    write_handoff_manifest(manifest, path)
    loaded = read_handoff_manifest(path)
    digest = manifest_digest(loaded)
    forbidden = [str(value) for value in loaded.get("forbidden_interpretations", [])]
    keys = set(str(key) for key in loaded)
    checks = [
        ("json_written", path.exists(), path),
        (
            "json_read",
            loaded["manifest_id"] == manifest.manifest_id,
            loaded["manifest_id"],
        ),
        ("digest_stable", digest == manifest.manifest_id, digest),
        ("forbidden_interpretations_present", bool(forbidden), forbidden),
        (
            "no_embedding_coordinates",
            "embedding_coordinates" not in keys and "coordinates" not in keys,
            sorted(keys),
        ),
        (
            "no_metric_interpretation_fields",
            "metric_reconstruction" not in keys and "spatial_distance" not in keys,
            sorted(keys),
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
    """Write exact manifest-read sanity CSV."""

    path = output_dir / "data" / "response_handoff_manifest_read_exact_sanity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote response handoff manifest read exact sanity: {output_path}")


if __name__ == "__main__":
    main()
