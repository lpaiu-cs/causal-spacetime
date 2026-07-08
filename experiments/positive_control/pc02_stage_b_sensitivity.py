"""PC-V1 Stage B: confirmatory sensitivity run on fresh seeds.

Refuses to run unless the frozen thresholds file exists (preregistration
Section 12). Applies frozen gates G1-G4 at the gate dimension to seeds
100-119. Decision rule: H-SENS supported iff >= 16 of 20 valid seeds pass.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from pc_common import (
    DEFAULT_OUTPUT_DIR,
    FROZEN_THRESHOLDS_PATH,
    git_describe,
    heldout_at_gate_dim,
    parse_seed_spec,
    primary_scene_config,
    run_seed_conditions,
    structured_row_at_gate_dim,
    write_rows_csv,
)

from causal_spacetime_lab.positive_control.gates import (
    RepresentabilityFitPolicy,
    evaluate_gates,
    load_frozen_thresholds,
)

STAGE_B_CONDITIONS = ["structured", "column_shuffled"]
PASS_RULE_NUMERATOR = 16
PASS_RULE_DENOMINATOR = 20


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", default="100-119", help="preregistered: 100-119")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        thresholds = load_frozen_thresholds(FROZEN_THRESHOLDS_PATH)
    except FileNotFoundError as error:
        print(error)
        sys.exit(2)

    seeds = parse_seed_spec(args.seeds)
    policy = RepresentabilityFitPolicy()
    code_version = git_describe()

    all_rows: list[dict[str, float | str]] = []
    decisions: list[dict[str, float | str]] = []
    for seed in seeds:
        rows = run_seed_conditions(
            primary_scene_config(seed), policy, STAGE_B_CONDITIONS, "B", code_version
        )
        all_rows.extend(rows)
        structured = structured_row_at_gate_dim(
            rows, "structured", seed, thresholds.gate_dim
        )
        null_value = heldout_at_gate_dim(
            rows, "column_shuffled", seed, thresholds.gate_dim
        )
        if structured is None or null_value is None:
            decisions.append(
                {"seed": float(seed), "status": "invalid", "all_gates_pass": 0.0}
            )
            continue
        gate_row = evaluate_gates(
            structured,
            null_value - structured["heldout_violation"],
            thresholds,
            apply_truth_gate=True,
        )
        decisions.append({"seed": float(seed), "status": "ok", **gate_row})

    passes = sum(int(d.get("all_gates_pass", 0.0)) for d in decisions)
    valid = sum(1 for d in decisions if d["status"] == "ok")
    supported = passes >= PASS_RULE_NUMERATOR and len(seeds) == PASS_RULE_DENOMINATOR

    registry = {
        "stage": "B",
        "frozen_commit": thresholds.frozen_commit,
        "seed_count": len(seeds),
        "valid_seed_count": valid,
        "gate_pass_count": passes,
        "pass_rule": f">={PASS_RULE_NUMERATOR}/{PASS_RULE_DENOMINATOR}",
        "h_sens_supported": supported,
        "code_version": code_version,
        "decisions": decisions,
    }

    write_rows_csv(args.output_dir / "stage_b_sensitivity.csv", all_rows)
    registry_path = args.output_dir / "stage_b_decision_registry.json"
    registry_path.parent.mkdir(parents=True, exist_ok=True)
    registry_path.write_text(json.dumps(registry, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in registry.items() if k != "decisions"}, indent=2))
    print(f"Commit {registry_path} under docs/prereg/frozen/ per Section 12.")


if __name__ == "__main__":
    main()
