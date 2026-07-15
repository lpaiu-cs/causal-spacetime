"""T1 numerical verification harness (docs/theory/t1_parallax_identifiability.md
Section 7).

Verifies the [PROVED] Model-D statements of the T1 document against the code
that the document claims to describe. This harness is analysis-only and
carries no gate: it verifies theory, it does not certify the instrument. But
because every claim it exercises is tagged [PROVED], any violation here is a
bug -- in the theory, in the instrument, or in this harness -- and is
reported as a hard failure, not a statistic.

The five checks, in the document's numbering:

1. Slope/quantization (Lemma 2, Model D): every measured bracket width obeys
   W = 2|dx|/delta + 1 + theta with theta in [-1, 1), on a controlled chain
   and on the full PC-V1 scene pipeline.
2. Fold (Theorem 1 sharpness): targets mirrored about an observer's position
   produce *identical* widths on that observer -- one observer cannot tell
   left from right -- and a second, offset observer separates them.
3. Resolution scaling (Section 5, Model D): position error is bounded by
   delta/2 and its RMSE falls like 1/K; varying the bulk sprinkling density
   at fixed K changes nothing (the direct falsifier for any accidental
   density dependence -- Model D's clock does not know about rho).
4. Centered residue (Lemma 3b): a pure time translation of a target moves
   the centered profile by O(1) ranks (the review's counterexample,
   reproduced exactly), and the residue never exceeds the proved bound.

Usage:
    python t1_verification.py            # run all checks, write JSON summary
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.discrete_radar import find_radar_ticks_from_order
from causal_spacetime_lab.observer import make_stationary_observer_chain_1p1
from causal_spacetime_lab.positive_control.echo_profiles import (
    measure_bracket_echo_profiles,
)
from causal_spacetime_lab.positive_control.scene import (
    PositiveControlSceneConfig,
    build_positive_control_scene,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "theory"

# Float fuzz allowance on the closed side of the proved band. The band itself
# is exact integer arithmetic; the tolerance only covers tick times built by
# np.linspace and the atol inside the causal matrix.
BAND_TOL = 1e-9


def bracket_widths_ranks(
    targets: np.ndarray,
    x0: float,
    span: float,
    ticks: int,
    extra_events: np.ndarray | None = None,
) -> np.ndarray:
    """Measured W (rank units) of each target against one stationary chain.

    Targets come first in the combined event array and the chain after them,
    so chain indices are stable when ``extra_events`` (bulk sprinkling that
    must NOT matter) are appended at the end. Unreachable targets get NaN.
    """

    chain_events, _ = make_stationary_observer_chain_1p1(span, ticks, x=x0)
    blocks = [np.asarray(targets, dtype=float), chain_events]
    if extra_events is not None:
        blocks.append(np.asarray(extra_events, dtype=float))
    events = np.vstack(blocks)
    causal = causal_matrix_1p1(events)
    chain_indices = np.arange(len(targets), len(targets) + ticks)
    ranks = np.arange(ticks, dtype=np.float64)

    widths = np.full(len(targets), np.nan)
    for j in range(len(targets)):
        bracket = find_radar_ticks_from_order(causal, chain_indices, j, ranks)
        if bracket is not None:
            rank_minus, rank_plus = bracket
            widths[j] = rank_plus - rank_minus
    return widths


def predicted_center(abs_dx: np.ndarray, delta: float) -> np.ndarray:
    """Lemma 2 (Model D) center: W = 2|dx|/delta + 1 + theta, theta in [-1, 1)."""

    return 2.0 * np.asarray(abs_dx) / delta + 1.0


def _sample_reachable_targets(
    rng: np.random.Generator,
    n: int,
    x0: float,
    span: float,
    margin: float = 0.05,
) -> np.ndarray:
    """Targets whose radar brackets fit inside the chain's tick range."""

    half = span / 2.0
    abs_dx = rng.uniform(0.01, half * 0.45, size=n)
    side = rng.choice([-1.0, 1.0], size=n)
    x = x0 + side * abs_dx
    t_room = half - abs_dx - margin
    t = rng.uniform(-t_room, t_room)
    return np.column_stack([t, x])


def check_quantization_band(
    seed: int = 0,
    n_targets: int = 400,
    span: float = 1.4,
    ticks: int = 96,
    x0: float = 0.0,
) -> dict:
    """Check 1 (controlled): residual W - (2|dx|/delta + 1) in [-1, 1)."""

    rng = np.random.default_rng(seed)
    delta = span / (ticks - 1)
    targets = _sample_reachable_targets(rng, n_targets, x0, span)
    widths = bracket_widths_ranks(targets, x0, span, ticks)
    reachable = ~np.isnan(widths)

    abs_dx = np.abs(targets[reachable, 1] - x0)
    residuals = widths[reachable] - predicted_center(abs_dx, delta)
    in_band = (residuals >= -1.0 - BAND_TOL) & (residuals < 1.0 + BAND_TOL)

    return {
        "n_reachable": int(reachable.sum()),
        "residual_min": float(residuals.min()),
        "residual_max": float(residuals.max()),
        "violations": int((~in_band).sum()),
        "passed": bool(in_band.all() and reachable.sum() > 0),
    }


def check_fold(span: float = 1.4, ticks: int = 96) -> dict:
    """Check 2: one observer folds mirrored targets; a second separates them.

    The pair (t, x0+d) / (t, x0-d) has identical |dx| to the observer at x0,
    hence identical predecessor/successor tick classes, hence *exactly* equal
    W -- integer equality, not approximate. An observer offset by s sees
    |d-s| vs |d+s| and separates the pair by ~4*min(d,s)/delta ranks.
    """

    x0, d, s, t = 0.0, 0.10, 0.15, 0.02
    delta = span / (ticks - 1)
    pair = np.array([[t, x0 + d], [t, x0 - d]])

    on_axis = bracket_widths_ranks(pair, x0, span, ticks)
    offset = bracket_widths_ranks(pair, x0 + s, span, ticks)

    separation_predicted = 2.0 * (abs(d + s) - abs(d - s)) / delta
    separation_measured = abs(offset[0] - offset[1])

    return {
        "on_axis_widths": [float(v) for v in on_axis],
        "fold_exact": bool(on_axis[0] == on_axis[1]),
        "offset_widths": [float(v) for v in offset],
        "separation_measured": float(separation_measured),
        "separation_predicted": float(separation_predicted),
        # The proved band allows each width +/-1 rank of quantization, so the
        # separation may deviate from its center by at most 2.
        "passed": bool(
            on_axis[0] == on_axis[1]
            and abs(separation_measured - separation_predicted) <= 2.0
        ),
    }


def check_density_invariance(
    seed: int = 1,
    n_targets: int = 60,
    bulk_sizes: tuple[int, ...] = (0, 300, 3000),
    span: float = 1.4,
    ticks: int = 96,
) -> dict:
    """Check 3b (falsifier): bulk sprinkling density must not move any width.

    Model D's clock is a deterministic grid appended to the scene; the bulk
    events are causally related to targets and ticks but never enter the
    bracket computation. If growing the bulk by 10x changed a single width,
    the document's central claim about the instrument would be false.
    """

    rng = np.random.default_rng(seed)
    targets = _sample_reachable_targets(rng, n_targets, 0.0, span)

    widths_by_bulk = []
    for n_bulk in bulk_sizes:
        extra = (
            sprinkle_1p1_causal_diamond(n_bulk, T=2.0, seed=seed + n_bulk)
            if n_bulk
            else None
        )
        widths_by_bulk.append(
            bracket_widths_ranks(targets, 0.0, span, ticks, extra_events=extra)
        )

    baseline = widths_by_bulk[0]
    identical = all(
        np.array_equal(baseline, other, equal_nan=True)
        for other in widths_by_bulk[1:]
    )
    return {
        "bulk_sizes": list(bulk_sizes),
        "identical_across_density": bool(identical),
        "passed": bool(identical),
    }


def check_resolution_scaling(
    seed: int = 2,
    n_targets: int = 200,
    tick_ladder: tuple[int, ...] = (12, 24, 48, 96, 192, 384),
    span: float = 1.4,
) -> dict:
    """Check 3a: position error bounded by delta/2, RMSE falling like 1/K.

    The position estimate is |dx|_hat = (W - 1) * delta / 2; the proved band
    |W - (2|dx|/delta + 1)| <= 1 gives |dx|_hat error <= delta/2 pointwise.
    The log-log RMSE slope against (K - 1) should sit near -1: quantization,
    not sampling noise.
    """

    rng = np.random.default_rng(seed)
    targets = _sample_reachable_targets(rng, n_targets, 0.0, span, margin=0.1)
    abs_dx = np.abs(targets[:, 1])

    rows = []
    for ticks in tick_ladder:
        delta = span / (ticks - 1)
        widths = bracket_widths_ranks(targets, 0.0, span, ticks)
        reachable = ~np.isnan(widths)
        estimate = (widths[reachable] - 1.0) * delta / 2.0
        errors = np.abs(estimate - abs_dx[reachable])
        rows.append({
            "ticks": ticks,
            "delta": delta,
            "n_reachable": int(reachable.sum()),
            "max_error": float(errors.max()),
            "rmse": float(np.sqrt(np.mean(errors**2))),
            "bound_delta_over_2": delta / 2.0,
            "within_bound": bool((errors <= delta / 2.0 + BAND_TOL).all()),
        })

    log_k = np.log([row["ticks"] - 1 for row in rows])
    log_rmse = np.log([row["rmse"] for row in rows])
    slope = float(np.polyfit(log_k, log_rmse, 1)[0])

    return {
        "ladder": rows,
        "rmse_loglog_slope": slope,
        "passed": bool(
            all(row["within_bound"] for row in rows) and -1.3 < slope < -0.7
        ),
    }


def pc_v1_centered_profile(t: float, x: float) -> np.ndarray:
    """Centered six-chain profile P[.,r] of one target on the PC-V1 grid."""

    config = PositiveControlSceneConfig()
    target = np.array([[t, x]])
    widths = np.array([
        bracket_widths_ranks(
            target, x0, config.chain_span, config.ticks_per_chain
        )[0]
        for x0 in config.chain_positions
    ])
    return widths - widths.mean()


def check_centered_residue(x: float = 0.10, t_shift: float = 0.003) -> dict:
    """Check 4: the residue centering does NOT remove, on the real grid.

    Reproduces the review's counterexample: translating a target purely in
    time moves the centered profile by O(1) ranks (so "no temporal quantity
    survives" holds only up to quantization), while the residue against the
    ideal profile 2*lambda*(|dx_r| - mean |dx|) stays within the proved
    bound of 2 ranks componentwise.
    """

    config = PositiveControlSceneConfig()
    delta = config.chain_span / (config.ticks_per_chain - 1)

    profile_a = pc_v1_centered_profile(0.0, x)
    profile_b = pc_v1_centered_profile(t_shift, x)

    abs_dx = np.abs(x - np.asarray(config.chain_positions))
    ideal = (2.0 * abs_dx / delta) - (2.0 * abs_dx / delta).mean()
    residue_a = profile_a - ideal
    residue_b = profile_b - ideal

    max_residue = float(max(np.abs(residue_a).max(), np.abs(residue_b).max()))
    time_shift_moved_profile = bool(np.any(profile_a != profile_b))

    return {
        "profile_t0": [float(v) for v in profile_a],
        "profile_shifted": [float(v) for v in profile_b],
        "time_shift_moved_profile": time_shift_moved_profile,
        "max_abs_residue": max_residue,
        "residue_bound": 2.0,
        "passed": bool(time_shift_moved_profile and max_residue < 2.0 + BAND_TOL),
    }


def check_pipeline_band(seed: int = 0) -> dict:
    """Check 1 (pipeline): the band holds through the full PC-V1 scene path.

    Same [PROVED] statement, but measured via build_positive_control_scene()
    and measure_bracket_echo_profiles() -- the instrument end to end -- so a
    divergence between the controlled harness and the real pipeline cannot
    hide in glue code.
    """

    config = PositiveControlSceneConfig(seed=seed)
    scene = build_positive_control_scene(config)
    profiles = measure_bracket_echo_profiles(scene)

    delta = config.chain_span / (config.ticks_per_chain - 1)
    coords = scene.target_coords
    residuals = []
    for row in range(profiles.target_count):
        for column, x0 in enumerate(config.chain_positions):
            if not profiles.reachable[row, column]:
                continue
            abs_dx = abs(coords[row, 1] - x0)
            width = profiles.delay_ranks[row, column]
            residuals.append(width - float(predicted_center(abs_dx, delta)))

    residuals = np.array(residuals)
    in_band = (residuals >= -1.0 - BAND_TOL) & (residuals < 1.0 + BAND_TOL)
    return {
        "scene_seed": seed,
        "n_measurements": int(residuals.size),
        "residual_min": float(residuals.min()),
        "residual_max": float(residuals.max()),
        "violations": int((~in_band).sum()),
        "passed": bool(in_band.all() and residuals.size > 0),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    checks = {
        "1_quantization_band_controlled": check_quantization_band(),
        "1_quantization_band_pipeline": check_pipeline_band(),
        "2_fold": check_fold(),
        "3a_resolution_scaling": check_resolution_scaling(),
        "3b_density_invariance": check_density_invariance(),
        "4_centered_residue": check_centered_residue(),
    }

    all_passed = all(result["passed"] for result in checks.values())
    for name, result in checks.items():
        status = "PASS" if result["passed"] else "FAIL"
        print(f"[{status}] {name}")
        for key, value in result.items():
            if key in ("passed", "ladder"):
                continue
            print(f"    {key} = {value}")
        if "ladder" in result:
            for row in result["ladder"]:
                print(
                    f"    K={row['ticks']:4d}  rmse={row['rmse']:.6f}  "
                    f"max={row['max_error']:.6f}  "
                    f"bound={row['bound_delta_over_2']:.6f}  "
                    f"within={row['within_bound']}"
                )

    summary = {
        "all_passed": all_passed,
        "checks": checks,
        "status": "analysis-only harness; no gate; verifies "
        "docs/theory/t1_parallax_identifiability.md Section 7",
    }
    (OUT / "t1_verification_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(f"\nall_passed = {all_passed}")
    print(f"wrote {OUT / 't1_verification_summary.json'}")

    if not all_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
