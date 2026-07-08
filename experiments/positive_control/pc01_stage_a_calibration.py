"""PC-V1 Stage A: instrument calibration on geometric scenes plus controls.

Exploratory by design: no gates are applied here. The output is a metric
distribution table and a MECHANICAL threshold proposal (preregistration
Section 10). Freezing the proposal is a separate human step (Section 12);
this script never writes into docs/prereg/frozen/.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pc_common import (
    DEFAULT_OUTPUT_DIR,
    cohens_d,
    git_describe,
    heldout_at_gate_dim,
    parse_seed_spec,
    primary_scene_config,
    run_seed_conditions,
    smoke_fit_policy,
    smoke_scene_config,
    structured_row_at_gate_dim,
    write_rows_csv,
)

from causal_spacetime_lab.positive_control.gates import (
    HARD_HELDOUT_CEILING,
    RepresentabilityFitPolicy,
    propose_thresholds,
    save_thresholds,
)

STAGE_A_CONDITIONS = [
    "structured",
    "column_shuffled",
    "relabel_symmetry",
    "random_order",
]
EFFECT_SIZE_FLOOR = 2.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", default="0-9", help="preregistered: 0-9")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="small engineering run; results carry no preregistered meaning",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    seeds = parse_seed_spec(args.seeds)
    policy = smoke_fit_policy() if args.smoke else RepresentabilityFitPolicy()
    code_version = git_describe()
    stage = "A-smoke" if args.smoke else "A"

    all_rows: list[dict[str, float | str]] = []
    for seed in seeds:
        config = (
            smoke_scene_config(seed) if args.smoke else primary_scene_config(seed)
        )
        rows = run_seed_conditions(
            config, policy, STAGE_A_CONDITIONS, stage, code_version
        )
        all_rows.extend(rows)
        print(f"seed {seed}: {sum(r['status'] == 'ok' for r in rows)} ok rows")

    output_dir = args.output_dir
    suffix = "_smoke" if args.smoke else ""
    csv_path = output_dir / f"stage_a_calibration{suffix}.csv"
    write_rows_csv(csv_path, all_rows)

    gate_dim = policy.gate_dim
    structured_heldout: list[float] = []
    null_heldout: list[float] = []
    null_gaps: list[float] = []
    structured_rows: list[dict[str, float]] = []
    for seed in seeds:
        structured = structured_row_at_gate_dim(
            all_rows, "structured", seed, gate_dim
        )
        null_value = heldout_at_gate_dim(
            all_rows, "column_shuffled", seed, gate_dim
        )
        if structured is None or null_value is None:
            continue
        structured_rows.append(structured)
        structured_heldout.append(structured["heldout_violation"])
        null_heldout.append(null_value)
        null_gaps.append(null_value - structured["heldout_violation"])

    if not structured_rows:
        print("no valid seeds; calibration failed structurally")
        return

    median_heldout = sorted(structured_heldout)[len(structured_heldout) // 2]
    effect_size = cohens_d(null_heldout, structured_heldout)
    f1_ok = median_heldout <= HARD_HELDOUT_CEILING
    f2_ok = effect_size >= EFFECT_SIZE_FLOOR

    summary = {
        "stage": stage,
        "valid_seed_count": len(structured_rows),
        "median_structured_heldout": median_heldout,
        "median_null_gap": sorted(null_gaps)[len(null_gaps) // 2],
        "effect_size_cohens_d": effect_size,
        "hard_floor_f1_pass": f1_ok,
        "hard_floor_f2_pass": f2_ok,
        "code_version": code_version,
    }

    proposal_path = output_dir / f"stage_a_threshold_proposal{suffix}.json"
    if f1_ok:
        thresholds = propose_thresholds(structured_rows, null_gaps, gate_dim)
        save_thresholds(thresholds, proposal_path)
        summary["threshold_proposal_path"] = str(proposal_path)
    else:
        summary["threshold_proposal_path"] = "WITHHELD (F1 hard floor failed)"

    summary_path = output_dir / f"stage_a_summary{suffix}.json"
    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    print(json.dumps(summary, indent=2))
    if args.smoke:
        print("SMOKE RUN: engineering check only; no preregistered meaning.")
    elif f1_ok and f2_ok:
        print(
            "Stage A complete. To freeze (Section 12): review the proposal, "
            "copy it to docs/prereg/frozen/pc_v1_thresholds.json with the "
            "freeze commit hash, flip the prereg status to FROZEN, commit."
        )
    else:
        print(
            "Hard floor failed (F1 or F2). Do not freeze; follow the "
            "escalation/stop rules in preregistration Section 11."
        )


if __name__ == "__main__":
    main()
