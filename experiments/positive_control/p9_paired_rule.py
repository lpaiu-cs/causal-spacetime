"""P9: 3+1D dimension selection under a within-seed decision rule.

See docs/prereg/p9_paired_dimension_rule.md. P9 changes ONLY the decision
rule relative to P8: the scene, the measurement, and every pipeline
primitive are P8's, reused by import. There is no calibration stage and
no gate -- every constant in the rule is a count or is inherited from
P2's frozen decision rule -- so this script has exactly one confirmatory
mode: seeds 500-519, run once.

The three legs, per seed with truth errors t_d, gate-dimension held-out
h3 and control held-out c3:

    underfit        t3 < t2            support requires >= 16 of 20
    no-improvement  t4 < t3 - 0.05     support requires <= 4 of 20
    paired spec     h3 < c3            support requires >= 16 of 20

H-KNEE is the conjunction of the first two; P9 overall is the
conjunction of all three. H-CEILING (the fraction with h3 <= 0.10) is
reported as a diagnostic and gates nothing.

The no-improvement leg is an equivalence claim and its pass is a
NON-REJECTION, not support; the registry records that distinction in its
own field rather than leaving it to prose.

Usage:
    python p9_paired_rule.py            # confirmatory, seeds 500-519
    python p9_paired_rule.py --seeds N  # explicitly non-confirmatory
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from p8_3plus1d import sweep_seed
from pc_common import DEFAULT_OUTPUT_DIR, git_describe, parse_seed_spec, write_rows_csv

from causal_spacetime_lab.positive_control.gates import RepresentabilityFitPolicy
from causal_spacetime_lab.positive_control.scene_3d import Scene3DConfig

CONFIRMATORY_SEEDS = "500-519"
PASS_MIN = 16          # P2's pass_min, frozen since P2
DENOMINATOR = 20       # P2's denominator
SATURATION_TOL = 0.05  # P2's saturation tolerance, frozen since P2
CEILING = 0.10         # programme-wide held-out ceiling (diagnostic only)


def decide(ok: list[dict], seed_count: int) -> dict:
    """Apply the preregistered paired rule to the valid rows.

    Counts only; no fitted quantity anywhere. The rule is evaluated on
    the valid-seed denominator as the prereg specifies (invalid seeds
    are excluded and never topped up), while the pass thresholds stay
    at their absolute 16-of-20 / 4-of-20 values -- so lost seeds make
    the rule strictly harder, never easier.
    """

    underfit = sum(1 for r in ok if r["d3_truth"] < r["d2_truth"])
    improve = sum(
        1 for r in ok if r["d4_truth"] < r["d3_truth"] - SATURATION_TOL
    )
    paired_spec = sum(
        1
        for r in ok
        if r["control_status"] != "ok"
        or (
            not math.isnan(r["control_d3_heldout"])
            and r["d3_heldout"] < r["control_d3_heldout"]
        )
    )
    under_ceiling = sum(1 for r in ok if r["d3_heldout"] <= CEILING)

    h_knee_underfit = underfit >= PASS_MIN
    h_knee_no_improvement = improve <= (DENOMINATOR - PASS_MIN)
    h_knee = h_knee_underfit and h_knee_no_improvement
    h_spec = paired_spec >= PASS_MIN

    return {
        "seed_count": seed_count,
        "valid_seed_count": len(ok),
        "pass_min": PASS_MIN,
        "denominator": DENOMINATOR,
        "saturation_tolerance": SATURATION_TOL,
        "underfit_count": underfit,
        "improvement_count": improve,
        "paired_spec_count": paired_spec,
        "h_knee_underfit_supported": h_knee_underfit,
        "h_knee_no_improvement_not_rejected": h_knee_no_improvement,
        "h_knee_no_improvement_is_a_non_rejection": True,
        "h_knee_supported": h_knee,
        "h_spec_paired_supported": h_spec,
        "p9_supported": h_knee and h_spec,
        "h_ceiling_fraction": (
            under_ceiling / len(ok) if ok else float("nan")
        ),
        "h_ceiling_note": (
            "descriptive diagnostic; gates nothing (prereg Section 2)"
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--seeds", default=None,
        help=f"default {CONFIRMATORY_SEEDS} (confirmatory, run once); "
             "anything else is labelled non-confirmatory",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    seed_spec = args.seeds or CONFIRMATORY_SEEDS
    confirmatory = seed_spec == CONFIRMATORY_SEEDS
    seeds = parse_seed_spec(seed_spec)
    policy = RepresentabilityFitPolicy()
    code_version = git_describe()
    stage = "P9" if confirmatory else "P9-smoke"

    rows = [
        sweep_seed(
            Scene3DConfig(seed=s, chain_count=12, min_bracketing_chains=12),
            policy, code_version, stage,
        )
        for s in seeds
    ]
    for r in rows:
        if r["status"] != "ok":
            print(f"seed {int(r['seed'])}: {r['status']}")
            continue
        print(
            f"seed {int(r['seed'])}: "
            f"d2t={r['d2_truth']:.3f} d3t={r['d3_truth']:.3f} "
            f"d4t={r['d4_truth']:.3f} | d3h={r['d3_heldout']:.3f} "
            f"ctrl={r['control_d3_heldout']:.3f} | "
            f"underfit={'Y' if r['d3_truth'] < r['d2_truth'] else 'n'} "
            f"improve={'Y' if r['d4_truth'] < r['d3_truth'] - SATURATION_TOL else 'n'} "
            f"spec={'Y' if r['d3_heldout'] < r['control_d3_heldout'] else 'n'}"
        )

    suffix = "" if confirmatory else "_smoke"
    write_rows_csv(args.output_dir / f"p9_paired_rule{suffix}.csv", rows)

    ok = [r for r in rows if r["status"] == "ok"]
    registry = {
        "stage": stage,
        "confirmatory": confirmatory,
        "code_version": code_version,
        "seed_spec": seed_spec,
        **decide(ok, len(seeds)),
    }
    path = args.output_dir / f"p9_decision_registry{suffix}.json"
    path.write_text(json.dumps(registry, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(registry, indent=2))


if __name__ == "__main__":
    main()
