"""Aggregate frozen-manifest representation diagnostics."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from manifest_representation_experiment_helpers import write_csv

DEFAULT_OUTPUT_DIR = Path("outputs")


def _read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as csv_file:
        return list(csv.DictReader(csv_file))


def _best_fit_rows(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    by_manifest: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        if row.get("fitted") not in {"1.0", "1", "True", "true"}:
            continue
        by_manifest.setdefault(row["manifest_id"], []).append(row)
    best: dict[str, dict[str, str]] = {}
    for manifest_id, manifest_rows in by_manifest.items():
        best[manifest_id] = min(
            manifest_rows,
            key=lambda row: float(row["heldout_violation_rate"]),
        )
    return best


def run_experiment(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, float | str]]:
    """Aggregate available Milestone 32 outputs by manifest."""

    data_dir = output_dir / "data"
    fit_rows = _read_rows(data_dir / "frozen_manifest_ordinal_representation.csv")
    null_rows = _read_rows(data_dir / "frozen_manifest_representation_nulls.csv")
    stability_rows = _read_rows(data_dir / "frozen_manifest_fit_stability.csv")
    no_fit_rows = _read_rows(data_dir / "failed_manifest_no_fit_controls.csv")

    manifest_ids = set()
    for rows in (fit_rows, null_rows, stability_rows, no_fit_rows):
        manifest_ids.update(
            row["manifest_id"] for row in rows if row.get("manifest_id")
        )

    best_fit = _best_fit_rows(fit_rows)
    stability_by_manifest = {
        row["manifest_id"]: row for row in stability_rows if row.get("manifest_id")
    }
    no_fit_by_manifest = {
        row["manifest_id"]: row for row in no_fit_rows if row.get("manifest_id")
    }
    null_by_manifest: dict[str, list[dict[str, str]]] = {}
    for row in null_rows:
        null_by_manifest.setdefault(row["manifest_id"], []).append(row)

    rows: list[dict[str, float | str]] = []
    for manifest_id in sorted(manifest_ids):
        fit = best_fit.get(manifest_id, {})
        nulls = null_by_manifest.get(manifest_id, [])
        stability = stability_by_manifest.get(manifest_id, {})
        no_fit = no_fit_by_manifest.get(manifest_id, {})
        null_delta_values = np.asarray(
            [
                float(row["structured_minus_null_mean"])
                for row in nulls
                if row.get("structured_minus_null_mean")
            ],
            dtype=float,
        )
        rows.append(
            {
                "manifest_id": manifest_id,
                "eligible": float(no_fit.get("eligible", fit.get("eligible", "nan"))),
                "best_heldout_violation": float(
                    fit.get("heldout_violation_rate", "nan")
                ),
                "best_dimension": float(fit.get("embedding_dim", "nan")),
                "null_summary_count": float(len(nulls)),
                "mean_structured_minus_null": (
                    float(np.nanmean(null_delta_values))
                    if null_delta_values.size
                    else float("nan")
                ),
                "restart_mean_heldout_violation": float(
                    stability.get("mean_heldout_violation_rate", "nan")
                ),
                "latent_order_disagreement": float(
                    stability.get("mean_pair_order_disagreement", "nan")
                ),
                "no_fit_status": no_fit.get("reason_not_fit", ""),
                "failed_reasons": no_fit.get("failed_reasons", ""),
            }
        )
    return rows


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write aggregate summary CSV."""

    return write_csv(
        rows,
        output_dir / "data" / "frozen_manifest_representation_summary.csv",
        [
            "manifest_id",
            "eligible",
            "best_heldout_violation",
            "best_dimension",
            "null_summary_count",
            "mean_structured_minus_null",
            "restart_mean_heldout_violation",
            "latent_order_disagreement",
            "no_fit_status",
            "failed_reasons",
        ],
    )


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote frozen manifest representation summary: {output_path}")


if __name__ == "__main__":
    main()
