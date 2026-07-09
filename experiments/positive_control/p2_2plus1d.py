"""P2: 2+1D robustness and dimension selection (see docs/prereg/p2_2plus1d.md).

The novel 2+1D test: with two spatial dimensions, the effective embedding
dimension must be 2. We fit d = 1, 2, 3 on measured 2+1D geometric order and on
a density-matched geometry-free control, and test H-SENS-2D (geometric passes
and recovers true 2D order at d = 2), H-DIM (d = 1 underfits while d = 2
suffices), and H-SPEC-2D (geometry-free blocks). Reuses the frozen PC-V1
pipeline primitives unchanged.

P2-A calibration seeds 0-9 propose thresholds mechanically; P2-B confirmatory
seeds 400-419 apply the frozen thresholds.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from pc_common import DEFAULT_OUTPUT_DIR, git_describe, parse_seed_spec, write_rows_csv

from causal_spacetime_lab.ordinal_embedding import (
    embedding_distance_order_error,
    fit_ordinal_embedding_gradient_descent,
    quadruplet_violation_rate,
)
from causal_spacetime_lab.positive_control.dissimilarity import (
    build_constraint_split,
    margin_from_probe_quantile,
    profile_dissimilarity_matrix,
)
from causal_spacetime_lab.positive_control.echo_profiles import (
    measure_bracket_echo_profiles,
)
from causal_spacetime_lab.positive_control.epsilon_sweep import build_epsilon_scene
from causal_spacetime_lab.positive_control.gates import RepresentabilityFitPolicy
from causal_spacetime_lab.positive_control.rewire import (
    geometric_post_closure_density,
)
from causal_spacetime_lab.positive_control.scene import SceneValidityError
from causal_spacetime_lab.positive_control.scene_2d import (
    Scene2DConfig,
    build_scene_2plus1d,
    target_positions_2d,
)

DIMS = (1, 2, 3)
GATE_DIM = 2
FROZEN_CONSTANTS_PATH = Path("docs/prereg/frozen/p2_test_constants.json")


def _fit_dims(profiles, truth_xy, policy, seed):
    """Return {dim: (heldout, truth_order_error)} for the 2+1D profiles."""

    dissimilarity = profile_dissimilarity_matrix(profiles, policy.min_common_columns)
    margin = margin_from_probe_quantile(
        dissimilarity, quantile=policy.margin_quantile, seed=seed + 3
    )
    split = build_constraint_split(
        dissimilarity,
        policy.train_constraints,
        policy.heldout_constraints,
        margin,
        train_fraction=policy.pair_train_fraction,
        seed=seed + 5,
    )
    out: dict[int, tuple[float, float]] = {}
    for dim in DIMS:
        coords, _ = fit_ordinal_embedding_gradient_descent(
            profiles.target_count, dim, split.train,
            steps=policy.steps, learning_rate=policy.learning_rate,
            seed=seed + 100 * dim, restarts=policy.restarts,
        )
        heldout = quadruplet_violation_rate(coords, split.heldout)
        truth = float("nan")
        if truth_xy is not None:
            truth = embedding_distance_order_error(
                coords, truth_xy, num_pair_comparisons=policy.truth_comparisons,
                seed=seed + 9,
            )
        out[dim] = (float(heldout), float(truth))
    return out


def sweep_seed(config, policy, code_version, stage) -> dict:
    base = {"stage": stage, "seed": float(config.seed), "code_version": code_version}
    try:
        scene = build_scene_2plus1d(config)
    except SceneValidityError as error:
        return {**base, "status": f"scene_invalid: {error}"}

    profiles = measure_bracket_echo_profiles(scene)
    truth_xy = target_positions_2d(scene)
    geo = _fit_dims(profiles, truth_xy, policy, config.seed)

    # geometry-free control: density-matched random order (epsilon = 1)
    target_density, _, _ = geometric_post_closure_density(scene)
    control_status = "ok"
    control_d2_heldout = float("nan")
    try:
        control, _ = build_epsilon_scene(scene, 1.0, config.seed + 41, target_density)
        control_profiles = measure_bracket_echo_profiles(control)
        control_d2_heldout = _fit_dims(control_profiles, None, policy, config.seed)[
            GATE_DIM
        ][0]
    except (SceneValidityError, ValueError) as error:
        control_status = f"structural_block: {str(error)[:60]}"

    row = {
        **base,
        **scene.provenance_row(),
        "status": "ok",
        "control_status": control_status,
        "control_d2_heldout": control_d2_heldout,
    }
    for dim in DIMS:
        row[f"d{dim}_heldout"] = geo[dim][0]
        row[f"d{dim}_truth"] = geo[dim][1]
    return row


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["a", "b"], default="a")
    parser.add_argument("--seeds", default=None, help="default 0-9 (a) / 400-419 (b)")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def _q(values, q):
    vals = sorted(values)
    return vals[min(len(vals) - 1, int(q * (len(vals) - 1)))] if vals else float("nan")


def main() -> None:
    args = parse_args()
    stage = args.stage
    seeds = parse_seed_spec(args.seeds or ("0-9" if stage == "a" else "400-419"))
    policy = RepresentabilityFitPolicy()
    code_version = git_describe()

    if stage == "b" and not FROZEN_CONSTANTS_PATH.exists():
        raise SystemExit(
            f"frozen P2 constants not found at {FROZEN_CONSTANTS_PATH}; P2-B may "
            "only run after the P2 freeze (preregistration Section 9)"
        )

    stage_label = f"P2-{stage.upper()}"
    rows = [
        sweep_seed(Scene2DConfig(seed=s), policy, code_version, stage_label)
        for s in seeds
    ]
    for r in rows:
        if r["status"] != "ok":
            print(f"seed {int(r['seed'])}: {r['status']}")
            continue
        print(
            f"seed {int(r['seed'])}: "
            f"d1t={r['d1_truth']:.3f} d2t={r['d2_truth']:.3f} d3t={r['d3_truth']:.3f} "
            f"d2h={r['d2_heldout']:.3f} "
            f"ctrl_d2h={r['control_d2_heldout']:.3f}({r['control_status'][:12]})"
        )
    write_rows_csv(args.output_dir / f"p2_2plus1d_{stage}.csv", rows)

    ok = [r for r in rows if r["status"] == "ok"]
    summary = {
        "stage": f"P2-{stage.upper()}", "code_version": code_version,
        "seed_count": len(seeds), "valid_seed_count": len(ok),
        "median_d1_truth": _median([r["d1_truth"] for r in ok]),
        "median_d2_truth": _median([r["d2_truth"] for r in ok]),
        "median_d3_truth": _median([r["d3_truth"] for r in ok]),
        "median_d2_heldout": _median([r["d2_heldout"] for r in ok]),
        "median_control_d2_heldout": _median(
            [
                r["control_d2_heldout"]
                for r in ok
                if not math.isnan(r["control_d2_heldout"])
            ]
        ),
        "control_block_count": sum(
            1
            for r in ok
            if r["control_status"] != "ok" or r["control_d2_heldout"] > 0.05
        ),
    }
    if stage == "a":
        def rnd(v):
            return float(round(v / 0.05) * 0.05)

        # Gates are placed at the MIDPOINT between the pass cluster (d=2 / the
        # geometric fit) and the fail cluster (d=1 truth for the truth gate; the
        # geometry-free control for the held-out gate), so both sides keep
        # margin. Anchoring to the pass cluster's p90 (the first rule) hugged
        # the d=2 cluster and left no margin; resolution does not move the
        # d=2-truth floor, so margin must come from gate placement (deviation
        # D1).
        d2_truth_max = max(r["d2_truth"] for r in ok)
        d1_truth_min = min(r["d1_truth"] for r in ok)
        d2_heldout_max = max(r["d2_heldout"] for r in ok)
        control_vals = [
            r["control_d2_heldout"]
            for r in ok
            if r["control_status"] == "ok"
            and not math.isnan(r["control_d2_heldout"])
        ]
        control_heldout_min = min(control_vals) if control_vals else 0.15
        summary["proposed_gate_truth"] = min(
            0.20, rnd(0.5 * (d2_truth_max + d1_truth_min))
        )
        summary["proposed_gate_heldout"] = min(
            0.10, rnd(0.5 * (d2_heldout_max + control_heldout_min))
        )
        summary["d1_truth_min"] = d1_truth_min
        summary["d2_truth_max"] = d2_truth_max
        summary["d2_heldout_max"] = d2_heldout_max
        summary["control_heldout_min"] = control_heldout_min
    (args.output_dir / f"p2_summary_{stage}.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))

    if stage == "b":
        _decide_confirmatory(ok, len(seeds), args.output_dir, code_version)


def _decide_confirmatory(ok, seed_count, output_dir, code_version):
    """Apply frozen P2 gates: H-SENS-2D, H-DIM, H-SPEC-2D (Section 7)."""

    const = json.loads(FROZEN_CONSTANTS_PATH.read_text(encoding="utf-8"))
    gt, gh = float(const["gate_truth"]), float(const["gate_heldout"])
    pass_min = int(const["pass_min"])

    sens = sum(1 for r in ok if r["d2_heldout"] <= gh and r["d2_truth"] <= gt)
    hdim = sum(1 for r in ok if r["d1_truth"] > gt and r["d2_truth"] <= gt)
    saturate = sum(1 for r in ok if r["d2_truth"] <= r["d3_truth"] + 0.05)
    spec = sum(
        1
        for r in ok
        if r["control_status"] != "ok" or r["control_d2_heldout"] > gh
    )
    registry = {
        "stage": "P2-B",
        "frozen_commit": const.get("frozen_commit"),
        "code_version": code_version,
        "seed_count": seed_count,
        "valid_seed_count": len(ok),
        "gate_truth": gt,
        "gate_heldout": gh,
        "pass_rule": f">={pass_min}/{seed_count}",
        "h_sens_2d_pass_count": sens,
        "h_sens_2d_supported": sens >= pass_min,
        "h_dim_pass_count": hdim,
        "h_dim_supported": hdim >= pass_min,
        "h_dim_saturation_count": saturate,
        "h_spec_2d_block_count": spec,
        "h_spec_2d_supported": spec >= pass_min,
    }
    path = output_dir / "p2_stage_b_decision_registry.json"
    path.write_text(json.dumps(registry, indent=2) + "\n", encoding="utf-8")
    print(
        f"\nP2-B: H-SENS-2D={registry['h_sens_2d_supported']} ({sens}/{len(ok)}) "
        f"H-DIM={registry['h_dim_supported']} ({hdim}/{len(ok)}, "
        f"saturate {saturate}) H-SPEC-2D={registry['h_spec_2d_supported']} "
        f"({spec}/{len(ok)})"
    )
    print(f"Commit {path} under docs/prereg/frozen/ per Section 9.")


def _median(values):
    vals = sorted(v for v in values if not (isinstance(v, float) and math.isnan(v)))
    return vals[len(vals) // 2] if vals else None


if __name__ == "__main__":
    main()
