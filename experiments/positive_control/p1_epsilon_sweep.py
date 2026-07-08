"""P1 geometry-dilution epsilon sweep (see docs/prereg/p1_epsilon_sweep.md).

P1-A (calibration, seeds 0-9): measure the truth-recovery dose-response curve
across the epsilon grid and propose the H-MONO test constant. No frozen test
is applied. P1-B (confirmatory, seeds 100-119): load the frozen constants,
apply H-MONO / H-THRESH / H-LAG and endpoint reproduction.

The frozen PC-V1 instrument is imported unchanged; P1 adds only the
epsilon-diluted order generator.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from pc_common import (
    DEFAULT_OUTPUT_DIR,
    git_describe,
    parse_seed_spec,
    primary_scene_config,
    smoke_scene_config,
    write_rows_csv,
)

from causal_spacetime_lab.positive_control.echo_profiles import (
    measure_bracket_echo_profiles,
)
from causal_spacetime_lab.positive_control.epsilon_sweep import build_epsilon_scene
from causal_spacetime_lab.positive_control.gates import (
    RepresentabilityFitPolicy,
    analyze_profiles,
    load_frozen_thresholds,
)
from causal_spacetime_lab.positive_control.rewire import (
    geometric_post_closure_density,
)
from causal_spacetime_lab.positive_control.scene import (
    SceneValidityError,
    build_positive_control_scene,
)

EPSILON_GRID = (0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.0)
SMOKE_GRID = (0.0, 0.5, 1.0)
MIN_TEST_CELLS = 6
FROZEN_CONSTANTS_PATH = Path("docs/prereg/frozen/p1_test_constants.json")
PC_V1_THRESHOLDS_PATH = Path("docs/prereg/frozen/pc_v1_thresholds.json")
DENSITY_TOLERANCE = 0.02


def _rank(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    for position, index in enumerate(order):
        ranks[index] = float(position)
    return ranks


def spearman(xs: list[float], ys: list[float]) -> float:
    if len(xs) < 3:
        return float("nan")
    rx, ry = _rank(xs), _rank(ys)
    mx = sum(rx) / len(rx)
    my = sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    if dx == 0.0 or dy == 0.0:
        return float("nan")
    return num / (dx * dy)


def crossing_epsilon(
    epsilons: list[float], values: list[float], level: float
) -> float | None:
    """First epsilon where ``values`` crosses ``level`` (linear interp)."""

    for i in range(1, len(epsilons)):
        lo, hi = values[i - 1], values[i]
        if (lo < level <= hi) or (lo >= level > hi):
            if hi == lo:
                return epsilons[i]
            frac = (level - lo) / (hi - lo)
            return epsilons[i - 1] + frac * (epsilons[i] - epsilons[i - 1])
    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["a", "b"], default="a")
    parser.add_argument("--seeds", default=None, help="default 0-9 (a) / 300-319 (b)")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def sweep_seed(config, grid, policy, code_version, stage) -> list[dict]:
    """Return one row per (epsilon, valid) at the gate dim for one seed."""

    seed = config.seed
    scene = build_positive_control_scene(config)
    target_density, _, _ = geometric_post_closure_density(scene)
    rows: list[dict] = []
    for epsilon in grid:
        base = {
            "stage": stage,
            "seed": float(seed),
            "epsilon": float(epsilon),
            "code_version": code_version,
        }
        try:
            eps_scene, density = build_epsilon_scene(
                scene, epsilon, seed + 41, target_density
            )
            profiles = measure_bracket_echo_profiles(eps_scene)
            metric_rows = analyze_profiles(
                profiles, eps_scene.target_coords[:, 1], policy, seed
            )
        except (SceneValidityError, ValueError) as error:
            rows.append({**base, "status": f"invalid: {error}"})
            continue
        gate_row = next(r for r in metric_rows if int(r["embedding_dim"]) == 1)
        rows.append(
            {
                **base,
                **eps_scene.provenance_row(),
                "status": "ok",
                "achieved_density": float(density),
                "target_density": float(target_density),
                "density_held": float(
                    abs(density - target_density) <= DENSITY_TOLERANCE
                ),
                "truth_order_error": gate_row["truth_distance_order_error"],
                "heldout_violation": gate_row["heldout_violation"],
                "restart_order_disagreement": gate_row["restart_order_disagreement"],
            }
        )
    return rows


def per_seed_curve(
    rows: list[dict], seed: int
) -> tuple[list[float], list[float], list[float]]:
    cells = sorted(
        (
            r
            for r in rows
            if r["status"] == "ok"
            and float(r["seed"]) == seed
            and float(r["density_held"]) == 1.0
        ),
        key=lambda r: r["epsilon"],
    )
    eps = [r["epsilon"] for r in cells]
    truth = [r["truth_order_error"] for r in cells]
    heldout = [r["heldout_violation"] for r in cells]
    return eps, truth, heldout


def main() -> None:
    args = parse_args()
    stage = args.stage
    default_seeds = "0-9" if stage == "a" else "300-319"
    seeds = parse_seed_spec(args.seeds or default_seeds)
    grid = SMOKE_GRID if args.smoke else EPSILON_GRID
    policy = RepresentabilityFitPolicy()
    code_version = git_describe()
    suffix = "_smoke" if args.smoke else ""

    if stage == "b" and not args.smoke and not FROZEN_CONSTANTS_PATH.exists():
        raise SystemExit(
            f"frozen P1 constants not found at {FROZEN_CONSTANTS_PATH}; P1-B may "
            "only run after the P1 freeze (preregistration Section 9)"
        )

    all_rows: list[dict] = []
    for seed in seeds:
        config = smoke_scene_config(seed) if args.smoke else primary_scene_config(seed)
        rows = sweep_seed(config, grid, policy, code_version, f"P1-{stage.upper()}")
        all_rows.extend(rows)
        held = sum(r["status"] == "ok" and r.get("density_held") == 1.0 for r in rows)
        print(f"seed {seed}: {held}/{len(grid)} density-held cells")

    write_rows_csv(args.output_dir / f"p1_epsilon_sweep_{stage}{suffix}.csv", all_rows)

    # A seed enters the monotonicity test only if its density-held curve has
    # enough coverage to span the grid (>= MIN_TEST_CELLS cells, both endpoints
    # present). Under-covered seeds are recorded and excluded from the test
    # denominator -- a 3-point curve cannot test monotonicity across the grid.
    min_cells = min(MIN_TEST_CELLS, len(grid))
    covered: list[dict] = []
    insufficient: list[int] = []
    for seed in seeds:
        eps, truth, heldout = per_seed_curve(all_rows, seed)
        if len(eps) < min_cells or eps[0] != 0.0 or eps[-1] != max(grid):
            insufficient.append(seed)
            continue
        covered.append(
            {
                "seed": seed,
                "rho": spearman(eps, truth),
                "cross_truth": crossing_epsilon(eps, truth, 0.15),
                "cross_heldout": crossing_epsilon(eps, heldout, 0.05),
                "eps0_truth": truth[0],
                "eps0_heldout": heldout[0],
                "eps1_truth": truth[-1],
                "eps1_heldout": heldout[-1],
            }
        )

    valid_rho = [c["rho"] for c in covered if c["rho"] == c["rho"]]
    lag = [
        c["cross_heldout"] - c["cross_truth"]
        for c in covered
        if c["cross_heldout"] is not None and c["cross_truth"] is not None
    ]
    summary = {
        "stage": f"P1-{stage.upper()}",
        "code_version": code_version,
        "seed_count": len(seeds),
        "covered_seed_count": len(covered),
        "insufficient_coverage_seeds": insufficient,
        "epsilon_grid": list(grid),
        "median_monotonicity_rho": _median(valid_rho),
        "min_monotonicity_rho": min(valid_rho) if valid_rho else None,
        "median_truth_crossing_eps": _median(
            [c["cross_truth"] for c in covered if c["cross_truth"] is not None]
        ),
        "median_heldout_crossing_eps": _median(
            [c["cross_heldout"] for c in covered if c["cross_heldout"] is not None]
        ),
        "median_heldout_minus_truth_crossing": _median(lag),
    }

    if stage == "a":
        p10 = _quantile(valid_rho, 0.10) if valid_rho else float("nan")
        rho_min = min(0.85, math.floor(p10 / 0.05) * 0.05) if valid_rho else None
        summary["proposed_rho_min"] = rho_min
        summary["hard_floor_hf1_note"] = "check epsilon=0 truth<=0.15 & heldout<=0.05"
        summary["hard_floor_hf2_note"] = "check epsilon=1 truth>=0.40"

    (args.output_dir / f"p1_summary_{stage}{suffix}.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))

    if stage == "b" and not args.smoke:
        _decide_confirmatory(covered, summary, args.output_dir, code_version)
    if args.smoke:
        print("SMOKE RUN: engineering check only; no preregistered meaning.")


def _decide_confirmatory(
    covered: list[dict], summary: dict, output_dir: Path, code_version: str
) -> None:
    """Apply frozen P1 constants + frozen PC-V1 gates (preregistration Section 7)."""

    constants = _load_p1_constants(FROZEN_CONSTANTS_PATH)
    pcv1 = load_frozen_thresholds(PC_V1_THRESHOLDS_PATH)
    rho_min = float(constants["rho_min"])
    n_covered = len(covered)

    mono_pass = sum(1 for c in covered if c["rho"] >= rho_min)
    h_mono = n_covered > 0 and mono_pass >= constants["mono_pass_fraction"] * n_covered

    lag = [
        c["cross_heldout"] - c["cross_truth"]
        for c in covered
        if c["cross_heldout"] is not None and c["cross_truth"] is not None
    ]
    h_lag = bool(lag) and (_median(lag) or 0.0) > 0.0

    thresh_estimable = sum(1 for c in covered if c["cross_truth"] is not None)
    h_thresh = thresh_estimable >= constants["thresh_estimable_min"]

    def gate_pass(truth: float, heldout: float) -> bool:
        return truth <= pcv1.truth_error_max and heldout <= pcv1.heldout_max

    endpoints_ok = sum(
        1
        for c in covered
        if gate_pass(c["eps0_truth"], c["eps0_heldout"])
        and not gate_pass(c["eps1_truth"], c["eps1_heldout"])
    )

    registry = {
        "stage": "P1-B",
        "frozen_commit": constants.get("frozen_commit"),
        "pc_v1_frozen_commit": pcv1.frozen_commit,
        "code_version": code_version,
        "covered_seed_count": n_covered,
        "rho_min": rho_min,
        "h_mono_pass_count": mono_pass,
        "h_mono_supported": bool(h_mono),
        "h_thresh_estimable_count": thresh_estimable,
        "h_thresh_supported": bool(h_thresh),
        "h_lag_median_gap": _median(lag),
        "h_lag_supported": bool(h_lag),
        "endpoint_reproduction_count": endpoints_ok,
        "endpoint_reproduction_supported": bool(
            endpoints_ok >= constants["endpoint_reproduction_min"]
        ),
        "median_truth_crossing_eps": summary["median_truth_crossing_eps"],
        "per_seed": covered,
    }
    path = output_dir / "p1_stage_b_decision_registry.json"
    path.write_text(json.dumps(registry, indent=2) + "\n", encoding="utf-8")
    print(
        f"\nP1-B: H-MONO={h_mono} ({mono_pass}/{n_covered}>=rho_min {rho_min}) "
        f"H-THRESH={h_thresh} H-LAG={h_lag} "
        f"endpoints={endpoints_ok}/{n_covered} reproduce"
    )
    print(f"Commit {path} under docs/prereg/frozen/ per Section 9.")


def _load_p1_constants(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(
            f"frozen P1 constants not found at {path}; P1-B may only run after "
            "the P1 freeze (preregistration Section 9)"
        )
    return json.loads(path.read_text(encoding="utf-8"))


def _median(values: list[float]) -> float | None:
    vals = sorted(v for v in values if v is not None)
    if not vals:
        return None
    return vals[len(vals) // 2]


def _quantile(values: list[float], q: float) -> float:
    vals = sorted(values)
    if not vals:
        return float("nan")
    idx = int(q * (len(vals) - 1))
    return vals[idx]


if __name__ == "__main__":
    main()
