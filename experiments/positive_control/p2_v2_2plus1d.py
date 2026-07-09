"""P2-v2: remediated 2+1D robustness (see docs/prereg/p2_v2_2plus1d.md).

Two upstream changes vs frozen P2: N = 3600 (scene-invalidity -> ~0) and a
coverage-aware confirmatory denominator (P1's rule). Reuses P2's measurement/
fit code (`sweep_seed`) and the midpoint gate rule unchanged; only the scene N
and the decision denominator differ. P2v2-A calibration seeds 0-9;
P2v2-B confirmatory fresh seeds 500-519.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from p2_2plus1d import sweep_seed  # reuse the P2 measurement/fit sweep
from pc_common import DEFAULT_OUTPUT_DIR, git_describe, parse_seed_spec, write_rows_csv

from causal_spacetime_lab.positive_control.gates import RepresentabilityFitPolicy
from causal_spacetime_lab.positive_control.scene_2d import Scene2DConfig

N_EVENTS_V2 = 3600
COVERAGE_FLOOR = 18
PASS_FRACTION = 0.80
FROZEN_CONSTANTS_PATH = Path("docs/prereg/frozen/p2_v2_test_constants.json")


def _median(values):
    vals = sorted(v for v in values if not (isinstance(v, float) and math.isnan(v)))
    return vals[len(vals) // 2] if vals else None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["a", "b"], default="a")
    parser.add_argument("--seeds", default=None, help="default 0-9 (a) / 500-519 (b)")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    stage = args.stage
    seeds = parse_seed_spec(args.seeds or ("0-9" if stage == "a" else "500-519"))
    policy = RepresentabilityFitPolicy()
    code_version = git_describe()

    if stage == "b" and not FROZEN_CONSTANTS_PATH.exists():
        raise SystemExit(
            f"frozen P2-v2 constants not found at {FROZEN_CONSTANTS_PATH}; P2v2-B "
            "may only run after the P2-v2 freeze (preregistration Section 9)"
        )

    rows = [
        sweep_seed(Scene2DConfig(n_events=N_EVENTS_V2, seed=s), policy,
                   code_version, f"P2v2-{stage.upper()}")
        for s in seeds
    ]
    ok = [r for r in rows if r["status"] == "ok"]
    invalid = [int(r["seed"]) for r in rows if r["status"] != "ok"]
    for r in rows:
        if r["status"] != "ok":
            print(f"seed {int(r['seed'])}: {r['status']}")
            continue
        print(
            f"seed {int(r['seed'])}: d1t={r['d1_truth']:.3f} d2t={r['d2_truth']:.3f} "
            f"d3t={r['d3_truth']:.3f} d2h={r['d2_heldout']:.3f} "
            f"ctrl_d2h={r['control_d2_heldout']:.3f}({r['control_status'][:10]})"
        )
    write_rows_csv(args.output_dir / f"p2_v2_2plus1d_{stage}.csv", rows)

    summary = {
        "stage": f"P2v2-{stage.upper()}", "code_version": code_version,
        "seed_count": len(seeds), "valid_seed_count": len(ok),
        "scene_invalid_seeds": invalid,
        "median_d1_truth": _median([r["d1_truth"] for r in ok]),
        "median_d2_truth": _median([r["d2_truth"] for r in ok]),
        "median_d3_truth": _median([r["d3_truth"] for r in ok]),
        "median_d2_heldout": _median([r["d2_heldout"] for r in ok]),
    }
    if stage == "a" and ok:
        def rnd(v):
            return float(round(v / 0.05) * 0.05)

        d2_truth_max = max(r["d2_truth"] for r in ok)
        d1_truth_min = min(r["d1_truth"] for r in ok)
        d2_heldout_max = max(r["d2_heldout"] for r in ok)
        control_vals = [
            r["control_d2_heldout"] for r in ok
            if r["control_status"] == "ok" and not math.isnan(r["control_d2_heldout"])
        ]
        control_min = min(control_vals) if control_vals else 0.15
        summary["proposed_gate_truth"] = min(
            0.20, rnd(0.5 * (d2_truth_max + d1_truth_min))
        )
        summary["proposed_gate_heldout"] = min(
            0.10, rnd(0.5 * (d2_heldout_max + control_min))
        )
        summary["d1_truth_min"] = d1_truth_min
        summary["d2_truth_max"] = d2_truth_max

    (args.output_dir / f"p2_v2_summary_{stage}.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))

    if stage == "b":
        _decide_confirmatory(ok, invalid, len(seeds), args.output_dir, code_version)


def _decide_confirmatory(ok, invalid, seed_count, output_dir, code_version):
    """Coverage-aware confirmatory decision (preregistration Section 7)."""

    const = json.loads(FROZEN_CONSTANTS_PATH.read_text(encoding="utf-8"))
    gt, gh = float(const["gate_truth"]), float(const["gate_heldout"])
    valid = len(ok)
    need = math.ceil(PASS_FRACTION * valid)
    covered = valid >= COVERAGE_FLOOR

    sens = sum(1 for r in ok if r["d2_heldout"] <= gh and r["d2_truth"] <= gt)
    hdim = sum(1 for r in ok if r["d1_truth"] > gt and r["d2_truth"] <= gt)
    saturate = sum(1 for r in ok if r["d2_truth"] <= r["d3_truth"] + 0.05)
    spec = sum(
        1 for r in ok
        if r["control_status"] != "ok" or r["control_d2_heldout"] > gh
    )
    registry = {
        "stage": "P2v2-B", "frozen_commit": const.get("frozen_commit"),
        "code_version": code_version, "seed_count": seed_count,
        "valid_seed_count": valid, "scene_invalid_seeds": invalid,
        "coverage_floor": COVERAGE_FLOOR, "coverage_ok": covered,
        "pass_needed": need, "gate_truth": gt, "gate_heldout": gh,
        "h_sens_2d_pass_count": sens,
        "h_sens_2d_supported": covered and sens >= need,
        "h_dim_pass_count": hdim, "h_dim_supported": covered and hdim >= need,
        "h_dim_saturation_count": saturate,
        "h_spec_2d_block_count": spec,
        "h_spec_2d_supported": covered and spec >= need,
    }
    path = output_dir / "p2_v2_stage_b_decision_registry.json"
    path.write_text(json.dumps(registry, indent=2) + "\n", encoding="utf-8")
    print(
        f"\nP2v2-B: valid {valid}/{seed_count} (coverage_ok={covered}, need {need}) "
        f"| H-SENS-2D={registry['h_sens_2d_supported']} ({sens}) "
        f"H-DIM={registry['h_dim_supported']} ({hdim}, sat {saturate}) "
        f"H-SPEC-2D={registry['h_spec_2d_supported']} ({spec})"
    )
    print(f"Commit {path} under docs/prereg/frozen/ per Section 9.")


if __name__ == "__main__":
    main()
