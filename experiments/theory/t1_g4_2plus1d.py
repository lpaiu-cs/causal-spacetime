"""G4a: which T1 statements survive in 2+1D, verified on the frozen scene.

T1's results were stated and checked in 1+1D. This harness separates the
ones whose proofs never used the dimension from the ones that did, and
verifies the former against the *existing* 2+1D instrument
(`build_scene_2plus1d`, frozen for P2/P2-v2: eight stationary chains on
a non-collinear ring, 96-tick deterministic grids).

Dimension enters T1's Lemmas 1-3 and Theorem 2 in exactly one place:
`|dx|`, the observer-to-target spatial separation, becomes a Euclidean
norm instead of an absolute value. Every step downstream -- the
partition of ticks into predecessor / spacelike / successor, the rank-gap
identity `W = N + 1`, the Model-D quantization band, Campbell's formula,
the shared-region cancellation and the same-slice concentricity -- is a
statement about times along one worldline and is blind to how many
spatial dimensions produced that separation. Checks 1-4 and 6-7 below
are those statements, exercised in 2+1D.

What genuinely changes is the *inverse problem*, and it changes twice.

- The fold is bigger (check 5). In 1+1D one observer confuses a mirror
  pair; in 2+1D it confuses an entire circle of equal radius, and two
  observers still confuse the pair mirrored across the line joining
  them. Three non-collinear observers separate them -- the geometric
  content of the roadmap's ">= 3 non-collinear observers" phrasing.
- Labeled identifiability becomes multilateration (check 8), with a
  bound that has no 1+1D analogue: the observer configuration enters
  through the conditioning of a linear system, so "non-collinear" is
  quantitative, not merely a genericity clause. Writing
  `d_r = |x - p_r|` for the target's distance to observer `r`,

      d_r^2 - d_1^2 = -2 x . (p_r - p_1) + |p_r|^2 - |p_1|^2,

  so the squared-distance differences are LINEAR in the unknown `x`:
  `A x = b` with rows `A_r = 2 (p_r - p_1)^T` and
  `b_r = |p_r|^2 - |p_1|^2 - (d_r^2 - d_1^2)`. `A` has rank 2 exactly
  when the observers are non-collinear, and then `x` is determined.
  Feeding the measured `d_hat_r = (W_r - 1) delta / 2`, which obeys
  `|d_hat_r - d_r| <= delta / 2` by Lemma 2's band, gives

      |b_hat_r - b_r| <= delta (d_r + d_1) + delta^2 / 2,
      |x_hat - x|     <= ||A^+||_2 ||b_hat - b||_2,

  which is `O(delta)`: positional recovery in 2+1D is
  resolution-limited exactly as in 1+1D, with a constant set by the
  observer geometry.

NOT in scope (this is G4a, not G4b): the unlabeled statements. Lemma 4
recovers the spatial *order* from the parallax dissimilarity `D` alone
via a strict Robinson (seriation) structure -- and "order" is a 1D
notion with no 2D analogue, while the piecewise-linearity that drove
the proof becomes conic in the plane. The 2+1D counterpart would be
"`D` determines the configuration up to similarity", a different
theorem needing distance-geometry / rigidity arguments. It is left
open, deliberately and explicitly.

Theory-track: no frozen gate consumes this, and nothing here is
preregistered. The scene builder it runs against IS frozen, and is used
unmodified. Deterministic seeds; the table is written to the tracked
path docs/theory/t1_g4_2plus1d_results.json.

Usage:
    python t1_g4_2plus1d.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from causal_spacetime_lab.causal import causal_matrix_minkowski
from causal_spacetime_lab.discrete_radar import find_radar_ticks_from_order
from causal_spacetime_lab.positive_control.scene_2d import (
    Scene2DConfig,
    build_scene_2plus1d,
)

ROOT = Path(__file__).resolve().parents[2]
RESULTS_PATH = ROOT / "docs" / "theory" / "t1_g4_2plus1d_results.json"

BAND_TOL = 1e-9
DIAMOND_T = 2.0
CHAIN_SPAN = 1.4


def residuals_in_band(residuals: np.ndarray) -> np.ndarray:
    """The proved half-open band ``[-1, 1)``; the upper edge stays open
    for the same reason as in 1+1D (a tick spacetime-coincident with the
    target realizes residual exactly +1 and must NOT pass)."""

    residuals = np.asarray(residuals)
    return (residuals >= -1.0 - BAND_TOL) & (residuals < 1.0 - BAND_TOL)


def chain_events(position: np.ndarray, ticks: int, span: float = CHAIN_SPAN):
    """A stationary 2+1D worldline: uniform tick grid at fixed (x, y)."""

    clock = np.linspace(-span / 2.0, span / 2.0, ticks)
    px, py = float(position[0]), float(position[1])
    return np.column_stack(
        (clock, np.full(clock.size, px), np.full(clock.size, py))
    )


def measure_widths(
    targets: np.ndarray,
    positions: np.ndarray,
    ticks: int,
    span: float = CHAIN_SPAN,
    extra_events: np.ndarray | None = None,
) -> np.ndarray:
    """Bracket widths ``W[j, r]`` of 2+1D targets against stationary chains.

    Targets come first in the event array and chains after them, so
    target and chain indices are stable when ``extra_events`` (bulk
    sprinkling that must not matter) is appended at the end.
    """

    targets = np.asarray(targets, dtype=float)
    blocks = [targets]
    chain_slices = []
    start = targets.shape[0]
    for position in positions:
        block = chain_events(np.asarray(position, dtype=float), ticks, span)
        blocks.append(block)
        chain_slices.append(np.arange(start, start + block.shape[0], dtype=int))
        start += block.shape[0]
    if extra_events is not None and len(extra_events):
        blocks.append(np.asarray(extra_events, dtype=float))

    events = np.vstack(blocks)
    causal = causal_matrix_minkowski(events)
    ranks = np.arange(ticks, dtype=np.float64)

    widths = np.full((targets.shape[0], len(chain_slices)), np.nan)
    for j in range(targets.shape[0]):
        for r, chain in enumerate(chain_slices):
            bracket = find_radar_ticks_from_order(causal, chain, j, ranks)
            if bracket is not None:
                widths[j, r] = bracket[1] - bracket[0]
    return widths


def radial_distances(targets: np.ndarray, positions: np.ndarray) -> np.ndarray:
    """Euclidean spatial separations ``|x_j - p_r|`` -- the only place the
    spatial dimension enters Lemmas 1-3."""

    targets = np.asarray(targets, dtype=float)
    positions = np.asarray(positions, dtype=float)
    delta = targets[:, None, 1:3] - positions[None, :, :]
    return np.sqrt(np.sum(delta * delta, axis=2))


def _ring(count: int, radius: float) -> np.ndarray:
    angles = np.linspace(0.0, 2.0 * np.pi, count, endpoint=False)
    return radius * np.column_stack((np.cos(angles), np.sin(angles)))


def _sample_targets(
    rng: np.random.Generator, n: int, band_t: float, band_r: float
) -> np.ndarray:
    """Targets in a small central band, well inside every chain's reach."""

    angle = rng.uniform(0.0, 2.0 * np.pi, size=n)
    radius = band_r * np.sqrt(rng.uniform(0.0, 1.0, size=n))
    t = rng.uniform(-band_t, band_t, size=n)
    return np.column_stack((t, radius * np.cos(angle), radius * np.sin(angle)))


def check_band(seed: int = 0, n_targets: int = 400, ticks: int = 96) -> dict:
    """Check 1: Lemma 2's band holds verbatim in 2+1D.

    The identity and band never used the dimension; if this fails, the
    claim that they transfer is simply false.
    """

    rng = np.random.default_rng(seed)
    delta = CHAIN_SPAN / (ticks - 1)
    positions = _ring(8, 0.25)
    targets = _sample_targets(rng, n_targets, 0.10, 0.22)
    widths = measure_widths(targets, positions, ticks)
    distances = radial_distances(targets, positions)

    reachable = ~np.isnan(widths)
    residuals = widths[reachable] - (2.0 * distances[reachable] / delta + 1.0)
    in_band = residuals_in_band(residuals)
    return {
        "n_measurements": int(reachable.sum()),
        "n_unreachable": int((~reachable).sum()),
        "residual_min": float(residuals.min()),
        "residual_max": float(residuals.max()),
        "violations": int((~in_band).sum()),
        "passed": bool(in_band.all() and reachable.all()),
    }


def check_resolution_scaling(
    seed: int = 2,
    n_targets: int = 200,
    tick_ladder: tuple[int, ...] = (12, 24, 48, 96, 192, 384),
) -> dict:
    """Check 2: the radial-distance estimator is resolution-limited in
    2+1D too -- pointwise error under ``delta / 2``, RMSE falling like
    ``1/K``."""

    rng = np.random.default_rng(seed)
    positions = _ring(3, 0.25)
    targets = _sample_targets(rng, n_targets, 0.08, 0.18)
    distances = radial_distances(targets, positions)

    rows = []
    for ticks in tick_ladder:
        delta = CHAIN_SPAN / (ticks - 1)
        widths = measure_widths(targets, positions, ticks)
        reachable = ~np.isnan(widths)
        estimate = (widths[reachable] - 1.0) * delta / 2.0
        errors = np.abs(estimate - distances[reachable])
        rows.append({
            "ticks": ticks,
            "delta": delta,
            "n_measurements": int(reachable.sum()),
            "max_error": float(errors.max()),
            "rmse": float(np.sqrt(np.mean(errors**2))),
            "bound_delta_over_2": delta / 2.0,
            "within_bound": bool((errors <= delta / 2.0 + BAND_TOL).all()),
        })

    slope = float(np.polyfit(
        np.log([row["ticks"] - 1 for row in rows]),
        np.log([row["rmse"] for row in rows]),
        1,
    )[0])
    return {
        "ladder": rows,
        "rmse_loglog_slope": slope,
        "passed": bool(
            all(row["within_bound"] for row in rows) and -1.3 < slope < -0.7
        ),
    }


def check_density_invariance(seed: int = 5, ticks: int = 96) -> dict:
    """Check 3: the Model-D clock does not know the sprinkling density,
    in 2+1D as in 1+1D -- widths must be bit-identical when bulk events
    are added at a fixed tick grid."""

    rng = np.random.default_rng(seed)
    positions = _ring(4, 0.25)
    targets = _sample_targets(rng, 40, 0.08, 0.18)

    reference = measure_widths(targets, positions, ticks)
    identical = True
    for n_bulk in (300, 900, 2700):
        bulk_rng = np.random.default_rng(1000 + n_bulk)
        angle = bulk_rng.uniform(0.0, 2.0 * np.pi, size=n_bulk)
        radius = 0.4 * np.sqrt(bulk_rng.uniform(0.0, 1.0, size=n_bulk))
        t = bulk_rng.uniform(-0.5, 0.5, size=n_bulk)
        bulk = np.column_stack(
            (t, radius * np.cos(angle), radius * np.sin(angle))
        )
        widths = measure_widths(targets, positions, ticks, extra_events=bulk)
        identical = identical and bool(
            np.array_equal(
                np.nan_to_num(widths, nan=-1.0),
                np.nan_to_num(reference, nan=-1.0),
            )
        )
    return {"bulk_sizes": [300, 900, 2700], "bit_identical": identical,
            "passed": identical}


def check_builder_density_invariance(ticks_seed: int = 0) -> dict:
    """Check 4: the same falsifier through the FROZEN 2+1D builder --
    scenes at different ``n_events`` must place bit-identical chain
    worldlines (the grid is independent of the sprinkling)."""

    # sizes chosen so every scene satisfies the builder's own validity
    # preconditions at this seed (below ~2600 the eligible-target floor
    # bites -- the scene-invalidity that motivated P2-v2's N = 3600)
    grid = (2600, 3000, 3600)
    blocks = []
    for n_events in grid:
        config = Scene2DConfig(n_events=n_events, seed=ticks_seed)
        scene = build_scene_2plus1d(config)
        chains = np.vstack([
            scene.events[chain] for chain in scene.chain_index_arrays
        ])
        blocks.append(chains)
    identical = all(
        np.array_equal(blocks[0], other) for other in blocks[1:]
    )
    return {
        "n_events_grid": list(grid),
        "chain_worldlines_bit_identical": identical,
        "passed": identical,
    }


def check_fold_structure(ticks: int = 96) -> dict:
    """Check 5: the fold is bigger in 2+1D, and three non-collinear
    observers close it.

    One observer cannot distinguish ANY two targets at equal radius --
    a whole circle, not just a mirror pair. Two observers still cannot
    distinguish the pair mirrored across the line joining them. A third
    observer off that line does, by more than the quantization spread.
    """

    t0 = 0.0
    ring_radius = 0.25
    # two observers on the x-axis: the mirror line is the x-axis itself
    p1 = np.array([-ring_radius, 0.0])
    p2 = np.array([ring_radius, 0.0])
    p3 = np.array([0.0, ring_radius])  # off the mirror line

    # a circle of targets at equal radius from p1 (single-observer fold)
    circle_radius = 0.12
    angles = np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False)
    circle = np.column_stack((
        np.full(angles.size, t0),
        p1[0] + circle_radius * np.cos(angles),
        p1[1] + circle_radius * np.sin(angles),
    ))
    circle_widths = measure_widths(circle, [p1], ticks)[:, 0]
    circle_identical = bool(
        np.all(~np.isnan(circle_widths))
        and np.all(circle_widths == circle_widths[0])
    )

    # a mirror pair across the p1-p2 line (the y -> -y reflection)
    mirror = np.array([[t0, 0.05, 0.09], [t0, 0.05, -0.09]])
    two_widths = measure_widths(mirror, [p1, p2], ticks)
    two_identical = bool(
        np.all(~np.isnan(two_widths))
        and np.array_equal(two_widths[0], two_widths[1])
    )
    three_widths = measure_widths(mirror, [p1, p2, p3], ticks)
    third_gap = float(abs(three_widths[0, 2] - three_widths[1, 2]))

    return {
        "circle_targets": int(angles.size),
        "single_observer_folds_the_circle": circle_identical,
        "two_observers_fold_the_mirror_pair": two_identical,
        "third_observer_width_gap": third_gap,
        "passed": bool(
            circle_identical and two_identical and third_gap > 2.0
        ),
    }


def multilaterate(
    widths: np.ndarray, positions: np.ndarray, delta: float
) -> tuple[np.ndarray, float]:
    """Solve ``A x = b`` from measured widths; also return ``||A^+||_2``.

    The linear system of the module docstring, using the anchor
    observer ``r = 1``. Least squares handles ``R > 3``.
    """

    positions = np.asarray(positions, dtype=float)
    estimates = (np.asarray(widths, dtype=float) - 1.0) * delta / 2.0
    anchor, rest = positions[0], positions[1:]
    matrix = 2.0 * (rest - anchor)
    rhs = (
        np.sum(rest * rest, axis=1)
        - float(anchor @ anchor)
        - (estimates[1:] ** 2 - estimates[0] ** 2)
    )
    solution, *_ = np.linalg.lstsq(matrix, rhs, rcond=None)
    pinv_norm = float(np.linalg.norm(np.linalg.pinv(matrix), 2))
    return solution, pinv_norm


def check_multilateration(
    seed: int = 7, n_targets: int = 200, ticks: int = 96
) -> dict:
    """Check 6 (Theorem 1, labeled, 2+1D): three non-collinear observers
    determine the target's 2D position, with the proved O(delta) bound.

    Also measured: the same solve on three NEARLY collinear observers,
    where the bound's ``||A^+||`` blows up -- the quantitative content of
    "non-collinear".
    """

    rng = np.random.default_rng(seed)
    delta = CHAIN_SPAN / (ticks - 1)
    targets = _sample_targets(rng, n_targets, 0.08, 0.18)
    truth = targets[:, 1:3]

    layouts = {
        "non_collinear": _ring(3, 0.25),
        "near_collinear": np.array([
            [-0.25, 0.0], [0.0, 0.004], [0.25, 0.0]
        ]),
    }
    report: dict = {}
    for name, positions in layouts.items():
        widths = measure_widths(targets, positions, ticks)
        distances = radial_distances(targets, positions)
        errors, bounds = [], []
        for j in range(targets.shape[0]):
            if np.isnan(widths[j]).any():
                continue
            estimate, pinv_norm = multilaterate(widths[j], positions, delta)
            errors.append(float(np.linalg.norm(estimate - truth[j])))
            # ||b_hat - b||_2 bound from the docstring, per component
            per_component = (
                delta * (distances[j, 1:] + distances[j, 0]) + delta**2 / 2.0
            )
            bounds.append(
                float(pinv_norm * np.linalg.norm(per_component))
            )
        errors_arr = np.asarray(errors)
        bounds_arr = np.asarray(bounds)
        report[name] = {
            "n_solved": int(errors_arr.size),
            "pinv_norm": float(
                np.linalg.norm(
                    np.linalg.pinv(2.0 * (positions[1:] - positions[0])), 2
                )
            ),
            "max_position_error": float(errors_arr.max()),
            "median_position_error": float(np.median(errors_arr)),
            "max_proved_bound": float(bounds_arr.max()),
            "all_within_proved_bound": bool(
                np.all(errors_arr <= bounds_arr + BAND_TOL)
            ),
        }
    report["passed"] = bool(
        report["non_collinear"]["all_within_proved_bound"]
        and report["near_collinear"]["all_within_proved_bound"]
        # the layout that is nearly degenerate must be measurably worse
        and report["near_collinear"]["median_position_error"]
        > 10.0 * report["non_collinear"]["median_position_error"]
    )
    return report


def _poisson_ticks(rng: np.random.Generator, lam: float, span: float):
    count = rng.poisson(lam * span)
    return np.sort(rng.uniform(-span / 2.0, span / 2.0, size=count))


def check_model_p_same_slice(
    seed: int = 11, lam: float = 40.0, gap: float = 0.02, trials: int = 4000
) -> dict:
    """Check 7 (Theorem 2 in 2+1D): same-slice pairs stay pathwise
    monotone, with the exact tie probability the 1+1D proof predicts.

    "Same slice" now means equal time and different *radial* distance;
    the brackets are still concentric, so the argument is unchanged --
    ties are possible, strict inversions are not.

    The tie constant counts observers, not dimensions. One observer's
    two annular regions have total measure ``2g``, so a tie needs an
    empty set of measure ``2g``: ``exp(-2 lambda g)``. The 1+1D
    document quotes ``exp(-4 lambda g)`` for the *flanking* estimator,
    which is two observers, each contributing ``2g`` -- verified here
    with a symmetric 2+1D layout that makes both radial gaps exactly
    ``g`` (observers at ``(-a, 0)`` and ``(a, 0)``, targets on the
    ``y`` axis, so each target is equidistant from both).
    """

    rng = np.random.default_rng(seed)
    span = 4.0

    # --- one observer: radial distances near / near + gap -------------
    near, far = 0.10, 0.10 + gap
    inversions = ties = 0
    for _ in range(trials):
        ticks = _poisson_ticks(rng, lam, span)
        # widths from the rank-gap identity: W = N + 1 with N the count
        # of ticks strictly inside the open radar interval
        w_near = int(np.sum((ticks > -near) & (ticks < near))) + 1
        w_far = int(np.sum((ticks > -far) & (ticks < far))) + 1
        if w_far < w_near:
            inversions += 1
        elif w_far == w_near:
            ties += 1
    single = {
        "strict_inversions": inversions,
        "tie_rate": ties / trials,
        "tie_rate_predicted": float(np.exp(-2.0 * lam * gap)),
    }

    # --- two observers, symmetric layout: both radial gaps equal g ----
    a = 0.25
    y_near = 0.10
    d_near = float(np.hypot(a, y_near))
    d_far = d_near + gap
    y_far = float(np.sqrt(d_far**2 - a**2))
    inversions2 = ties2 = 0
    for _ in range(trials):
        ticks_1 = _poisson_ticks(rng, lam, span)
        ticks_2 = _poisson_ticks(rng, lam, span)
        total_near = total_far = 0
        for ticks in (ticks_1, ticks_2):
            total_near += int(np.sum((ticks > -d_near) & (ticks < d_near))) + 1
            total_far += int(np.sum((ticks > -d_far) & (ticks < d_far))) + 1
        if total_far < total_near:
            inversions2 += 1
        elif total_far == total_near:
            ties2 += 1
    flanking = {
        "observer_offset": a,
        "radial_gap": float(d_far - d_near),
        "y_near": y_near,
        "y_far": y_far,
        "strict_inversions": inversions2,
        "tie_rate": ties2 / trials,
        "tie_rate_predicted": float(np.exp(-4.0 * lam * gap)),
    }

    def _agrees(block: dict) -> bool:
        # binomial standard error at these rates is ~0.006; 3 se
        return abs(block["tie_rate"] - block["tie_rate_predicted"]) < 0.02

    return {
        "trials": trials,
        "lam": lam,
        "gap": gap,
        "single_observer": single,
        "two_observer_flanking": flanking,
        "passed": bool(
            single["strict_inversions"] == 0
            and flanking["strict_inversions"] == 0
            and _agrees(single)
            and _agrees(flanking)
        ),
    }


def check_pipeline_band(seed: int = 0) -> dict:
    """Check 8: the band end to end through the frozen 2+1D builder --
    every target, every chain, real sprinkled scene."""

    scene = build_scene_2plus1d(Scene2DConfig(seed=seed))
    delta = CHAIN_SPAN / (Scene2DConfig().ticks_per_chain - 1)
    positions = np.array([
        scene.events[chain[0]][1:3] for chain in scene.chain_index_arrays
    ])
    targets = scene.events[scene.target_indices]
    distances = radial_distances(targets, positions)

    residuals = []
    for j, target in enumerate(scene.target_indices):
        for r, chain in enumerate(scene.chain_index_arrays):
            bracket = find_radar_ticks_from_order(
                scene.causal, chain, int(target), scene.tick_ranks
            )
            if bracket is None:
                continue
            width = bracket[1] - bracket[0]
            residuals.append(width - (2.0 * distances[j, r] / delta + 1.0))
    residuals_arr = np.asarray(residuals)
    in_band = residuals_in_band(residuals_arr)
    return {
        "n_targets": int(scene.target_indices.size),
        "n_measurements": int(residuals_arr.size),
        "residual_min": float(residuals_arr.min()),
        "residual_max": float(residuals_arr.max()),
        "violations": int((~in_band).sum()),
        "passed": bool(in_band.all()),
    }


CHECKS = (
    ("band", check_band),
    ("resolution_scaling", check_resolution_scaling),
    ("density_invariance", check_density_invariance),
    ("builder_density_invariance", check_builder_density_invariance),
    ("fold_structure", check_fold_structure),
    ("multilateration", check_multilateration),
    ("model_p_same_slice", check_model_p_same_slice),
    ("pipeline_band", check_pipeline_band),
)


def main() -> None:
    results: dict = {
        "scope": (
            "G4a: T1 statements whose proofs never used the spatial "
            "dimension, verified in 2+1D on the frozen scene builder, "
            "plus labeled multilateration. G4b (unlabeled recovery of "
            "the configuration from D alone) is NOT addressed: Lemma 4's "
            "Robinson/seriation engine has no 2D analogue."
        ),
        "checks": {},
    }
    for name, check in CHECKS:
        outcome = check()
        results["checks"][name] = outcome
        flag = "PASS" if outcome["passed"] else "FAIL"
        print(f"[{flag}] {name}: {json.dumps(outcome)[:160]}")

    all_passed = all(row["passed"] for row in results["checks"].values())
    results["all_passed"] = bool(all_passed)
    results["status"] = (
        "theory-track verification; no gate; no frozen artifact modified "
        "(the 2+1D scene builder is used unchanged)"
    )
    RESULTS_PATH.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"\nall_passed = {all_passed}")
    print(f"wrote {RESULTS_PATH}")
    if not all_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
