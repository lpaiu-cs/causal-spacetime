"""P8: 3+1D robustness and dimension selection (see docs/prereg/p8_3plus1d.md).

The P2 question one dimension up. With three spatial dimensions the effective
embedding dimension must be 3. We fit d = 2, 3, 4 on measured 3+1D geometric
order and on a density-matched geometry-free control, and test H-SENS-3D
(geometric passes and recovers true 3D order at d = 3), H-DIM-3D (d = 2
underfits while d = 3 suffices), and H-SPEC-3D (geometry-free blocks). Reuses
the frozen PC-V1 pipeline primitives unchanged.

P8-A calibration seeds 100-119 propose thresholds mechanically; P8-B
confirmatory seeds 400-419 apply the frozen thresholds. Stage A deliberately
avoids seeds 0-3, which informed the scene choice.

Two numbers here differ from P2 and both follow from one measured fact: the
pass and fail clusters sit about four times closer together in 3+1D than in
2+1D (preregistration Section 4).

* Twelve observer chains rather than eight, which is what restores the
  separation to a workable size.
* A gate rounding grid of 0.01 rather than P2's 0.05. P2's grid is comparable
  to P8's ENTIRE cluster gap, so rounding could carry the gate out of the gap
  it is meant to sit in. A grid has to be small against the interval it
  discretises.
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
from causal_spacetime_lab.positive_control.scene_3d import (
    Scene3DConfig,
    build_scene_3plus1d,
    target_positions_3d,
)

DIMS = (2, 3, 4)
GATE_DIM = 3
GATE_GRID = 0.01
FROZEN_CONSTANTS_PATH = Path("docs/prereg/frozen/p8_test_constants.json")


def _median(values):
    vals = sorted(values)
    if not vals:
        return float("nan")
    mid = len(vals) // 2
    return vals[mid] if len(vals) % 2 else 0.5 * (vals[mid - 1] + vals[mid])


def _fit_dims(profiles, truth_xyz, policy, seed):
    """Return {dim: (heldout, truth_order_error)} for the 3+1D profiles."""

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
        if truth_xyz is not None:
            truth = embedding_distance_order_error(
                coords, truth_xyz, num_pair_comparisons=policy.truth_comparisons,
                seed=seed + 9,
            )
        out[dim] = (float(heldout), float(truth))
    return out


def sweep_seed(config, policy, code_version, stage) -> dict:
    base = {
        "stage": stage,
        "seed": float(config.seed),
        "chain_count": float(config.chain_count),
        "code_version": code_version,
    }
    try:
        scene = build_scene_3plus1d(config)
    except SceneValidityError as error:
        return {**base, "status": f"scene_invalid: {error}"}

    profiles = measure_bracket_echo_profiles(scene)
    truth_xyz = target_positions_3d(scene)
    geo = _fit_dims(profiles, truth_xyz, policy, config.seed)

    # geometry-free control: density-matched random order (epsilon = 1)
    target_density, _, _ = geometric_post_closure_density(scene)
    control_status = "ok"
    control_gate_heldout = float("nan")
    try:
        control, _ = build_epsilon_scene(scene, 1.0, config.seed + 41, target_density)
        control_profiles = measure_bracket_echo_profiles(control)
        control_gate_heldout = _fit_dims(
            control_profiles, None, policy, config.seed
        )[GATE_DIM][0]
    except SceneValidityError as error:
        control_status = f"structural_block: {str(error)[:60]}"

    row = {
        **base,
        **scene.provenance_row(),
        "status": "ok",
        "control_status": control_status,
        "control_d3_heldout": control_gate_heldout,
    }
    for dim in DIMS:
        row[f"d{dim}_heldout"] = geo[dim][0]
        row[f"d{dim}_truth"] = geo[dim][1]
    return row


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["a", "b"], default="a")
    parser.add_argument(
        "--seeds", default=None, help="default 100-119 (a) / 400-419 (b)"
    )
    parser.add_argument(
        "--chains", type=int, default=12,
        help="observer chains; 12 is the preregistered primary, 8 the "
             "descriptive H-OBS contrast",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    stage = args.stage
    seeds = parse_seed_spec(
        args.seeds or ("100-119" if stage == "a" else "400-419")
    )
    policy = RepresentabilityFitPolicy()
    code_version = git_describe()

    if stage == "b" and not FROZEN_CONSTANTS_PATH.exists():
        raise SystemExit(
            f"frozen P8 constants not found at {FROZEN_CONSTANTS_PATH}; P8-B may "
            "only run after the P8 freeze (preregistration Section 9)"
        )

    suffix = "" if args.chains == 12 else f"_r{args.chains}"
    stage_label = f"P8-{stage.upper()}{suffix}"
    rows = [
        sweep_seed(
            Scene3DConfig(
                seed=s, chain_count=args.chains,
                min_bracketing_chains=args.chains,
            ),
            policy, code_version, stage_label,
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
            f"d4t={r['d4_truth']:.3f} d3h={r['d3_heldout']:.3f} "
            f"ctrl_d3h={r['control_d3_heldout']:.3f}({r['control_status'][:12]})"
        )
    write_rows_csv(args.output_dir / f"p8_3plus1d_{stage}{suffix}.csv", rows)

    ok = [r for r in rows if r["status"] == "ok"]
    summary = {
        "stage": stage_label, "code_version": code_version,
        "chain_count": args.chains,
        "seed_count": len(seeds), "valid_seed_count": len(ok),
        "median_d2_truth": _median([r["d2_truth"] for r in ok]),
        "median_d3_truth": _median([r["d3_truth"] for r in ok]),
        "median_d4_truth": _median([r["d4_truth"] for r in ok]),
        "median_d3_heldout": _median([r["d3_heldout"] for r in ok]),
        "median_control_d3_heldout": _median(
            [
                r["control_d3_heldout"]
                for r in ok
                if not math.isnan(r["control_d3_heldout"])
            ]
        ),
        "control_block_count": sum(
            1
            for r in ok
            if r["control_status"] != "ok" or r["control_d3_heldout"] > 0.05
        ),
    }
    if stage == "a" and ok:
        def rnd(value):
            return float(round(value / GATE_GRID) * GATE_GRID)

        # Same midpoint rule as P2: the gate sits between the pass cluster
        # (the geometric fit at the gate dimension) and the fail cluster
        # (d = 2 truth for the truth gate; the geometry-free control for the
        # held-out gate), so both sides keep margin. The grid is finer than
        # P2's because the gap it discretises is about four times narrower.
        d3_truth_max = max(r["d3_truth"] for r in ok)
        d2_truth_min = min(r["d2_truth"] for r in ok)
        d3_heldout_max = max(r["d3_heldout"] for r in ok)
        control_vals = [
            r["control_d3_heldout"]
            for r in ok
            if r["control_status"] == "ok"
            and not math.isnan(r["control_d3_heldout"])
        ]
        control_heldout_min = min(control_vals) if control_vals else 0.15
        clusters_overlap = d3_truth_max >= d2_truth_min
        summary["clusters_overlap"] = clusters_overlap
        summary["d2_truth_min"] = d2_truth_min
        summary["d3_truth_max"] = d3_truth_max
        summary["d3_heldout_max"] = d3_heldout_max
        summary["control_heldout_min"] = control_heldout_min
        summary["cluster_gap_truth"] = d2_truth_min - d3_truth_max
        if clusters_overlap:
            # Preregistration Section 6: no gate is placed, and that is a
            # reportable outcome rather than grounds for changing the scene.
            summary["proposed_gate_truth"] = None
            summary["proposed_gate_heldout"] = None
            summary["note"] = (
                "pass and fail truth clusters overlap; no gate placed "
                "(prereg Section 6)"
            )
        else:
            summary["proposed_gate_truth"] = min(
                0.20, rnd(0.5 * (d3_truth_max + d2_truth_min))
            )
            summary["proposed_gate_heldout"] = min(
                0.10, rnd(0.5 * (d3_heldout_max + control_heldout_min))
            )
    (args.output_dir / f"p8_summary_{stage}{suffix}.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))

    if stage == "b":
        _decide_confirmatory(ok, len(seeds), args.output_dir, code_version, suffix)


def _decide_confirmatory(ok, seed_count, output_dir, code_version, suffix):
    """Apply frozen P8 gates: H-SENS-3D, H-DIM-3D, H-SPEC-3D (Section 7)."""

    const = json.loads(FROZEN_CONSTANTS_PATH.read_text(encoding="utf-8"))
    gt, gh = float(const["gate_truth"]), float(const["gate_heldout"])
    pass_min = int(const["pass_min"])

    sens = sum(1 for r in ok if r["d3_heldout"] <= gh and r["d3_truth"] <= gt)
    hdim = sum(1 for r in ok if r["d2_truth"] > gt and r["d3_truth"] <= gt)
    saturate = sum(1 for r in ok if r["d3_truth"] <= r["d4_truth"] + 0.05)
    spec = sum(
        1
        for r in ok
        if r["control_status"] != "ok" or r["control_d3_heldout"] > gh
    )
    registry = {
        "stage": f"P8-B{suffix}",
        "frozen_commit": const.get("frozen_commit"),
        "code_version": code_version,
        "seed_count": seed_count,
        "valid_seed_count": len(ok),
        "gate_truth": gt,
        "gate_heldout": gh,
        "pass_rule": f">={pass_min}/{seed_count}",
        "h_sens_3d_pass_count": sens,
        "h_sens_3d_supported": sens >= pass_min,
        "h_dim_3d_pass_count": hdim,
        "h_dim_3d_supported": hdim >= pass_min,
        "h_dim_saturation_count": saturate,
        "h_spec_3d_block_count": spec,
        "h_spec_3d_supported": spec >= pass_min,
    }
    path = output_dir / f"p8_stage_b_decision_registry{suffix}.json"
    path.write_text(json.dumps(registry, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(registry, indent=2))


if __name__ == "__main__":
    main()
