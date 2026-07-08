"""PC-V1 Stage C: confirmatory specificity run on geometry-free controls.

Refuses to run unless the frozen thresholds file exists. For each control
family, H-SPEC is supported iff >= 8 of 10 valid seeds FAIL gate G1 or G2
(the truth gate G4 is not applicable to controls). A structural failure to
build a valid constraint pool counts as a block (gate G5).
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

CONTROL_FAMILIES = {
    "random_order": ("random_order", "random_order_shuffled"),
    "column_shuffled_geometric": ("column_shuffled", "column_shuffled_repeat"),
}
BLOCK_RULE_NUMERATOR = 8
BLOCK_RULE_DENOMINATOR = 10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", default="200-209", help="preregistered: 200-209")
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
    conditions = sorted({c for pair in CONTROL_FAMILIES.values() for c in pair})

    all_rows: list[dict[str, float | str]] = []
    for seed in seeds:
        all_rows.extend(
            run_seed_conditions(
                primary_scene_config(seed), policy, conditions, "C", code_version
            )
        )

    families: dict[str, dict[str, float | str | list]] = {}
    for family, (main_condition, null_condition) in CONTROL_FAMILIES.items():
        decisions: list[dict[str, float | str]] = []
        for seed in seeds:
            main_row = structured_row_at_gate_dim(
                all_rows, main_condition, seed, thresholds.gate_dim
            )
            null_value = heldout_at_gate_dim(
                all_rows, null_condition, seed, thresholds.gate_dim
            )
            if main_row is None or null_value is None:
                decisions.append(
                    {"seed": float(seed), "status": "structural_block", "blocked": 1.0}
                )
                continue
            gate_row = evaluate_gates(
                main_row,
                null_value - main_row["heldout_violation"],
                thresholds,
                apply_truth_gate=False,
            )
            blocked = not (
                bool(gate_row["g1_heldout_pass"]) and bool(gate_row["g2_null_gap_pass"])
            )
            decisions.append(
                {
                    "seed": float(seed),
                    "status": "ok",
                    "blocked": float(blocked),
                    **gate_row,
                }
            )
        block_count = sum(int(d["blocked"]) for d in decisions)
        families[family] = {
            "block_count": block_count,
            "block_rule": f">={BLOCK_RULE_NUMERATOR}/{BLOCK_RULE_DENOMINATOR}",
            "h_spec_supported": block_count >= BLOCK_RULE_NUMERATOR
            and len(seeds) == BLOCK_RULE_DENOMINATOR,
            "decisions": decisions,
        }

    registry = {
        "stage": "C",
        "frozen_commit": thresholds.frozen_commit,
        "seed_count": len(seeds),
        "families": families,
        "code_version": code_version,
    }

    write_rows_csv(args.output_dir / "stage_c_specificity.csv", all_rows)
    registry_path = args.output_dir / "stage_c_decision_registry.json"
    registry_path.parent.mkdir(parents=True, exist_ok=True)
    registry_path.write_text(json.dumps(registry, indent=2) + "\n", encoding="utf-8")
    for family, result in families.items():
        print(
            f"{family}: blocked {result['block_count']}/{len(seeds)} "
            f"h_spec_supported={result['h_spec_supported']}"
        )
    print(f"Commit {registry_path} under docs/prereg/frozen/ per Section 12.")


if __name__ == "__main__":
    main()
