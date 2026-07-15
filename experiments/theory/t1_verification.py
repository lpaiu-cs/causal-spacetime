"""T1 numerical verification harness (docs/theory/t1_parallax_identifiability.md
Section 7).

Verifies the [PROVED] Model-D statements of the T1 document against the code
that the document claims to describe. This harness is analysis-only and
carries no gate: it verifies theory, it does not certify the instrument. But
because every claim it exercises is tagged [PROVED], any violation here is a
bug -- in the theory, in the instrument, or in this harness -- and is
reported as a hard failure, not a statistic.

The checks, in the document's numbering:

1. Slope/quantization (Lemma 2, Model D): every measured bracket width obeys
   W = 2|dx|/delta + 1 + theta with theta in [-1, 1), on a controlled chain
   and on the full PC-V1 scene pipeline.
2. Fold (Theorem 1 sharpness): targets mirrored about an observer's position
   produce *identical* widths on that observer -- one observer cannot tell
   left from right -- and a second, offset observer separates them.
3. Resolution scaling (Section 5, Model D): position error is bounded by
   delta/2 and its RMSE falls like 1/K; varying the bulk sprinkling density
   at fixed K changes nothing (the falsifier for any accidental density
   dependence -- Model D's clock does not know about rho), asserted both on
   a hand-built order (unit invariant) and through
   build_positive_control_scene() itself (instrument regression).
4. Centered residue (Lemma 3b): a pure time translation of a target moves
   the centered profile by O(1) ranks (the review's counterexample,
   reproduced exactly), and the residue never exceeds the proved bound.
5. Edge and convention pinning: the two deterministic cases where a wrong
   instrument would otherwise agree with every sampled check -- the
   tick-coincident orphan (the one true band violation; the predicate must
   reject it) and Lemma 2's null-aligned example (the one configuration
   where the null-inclusive and strict orders give different widths; pins
   the convention the document declares load-bearing).
6. Unlabeled decoding (Lemma 4): the exact-model dissimilarity D is
   strictly Robinson (gap-direction inner products 4 a_k b_l / R, strict
   quadrance superadditivity) and the anchor decoder recovers the spatial
   order up to reversal from D alone, including at R = 2; on the pipeline,
   |D_measured - D_exact| < 4 (under the asserted shared-centering
   hypothesis) and every anchor comparison with exact margin above 8 is
   ordered correctly -- the Model-D claim at its proved strength,
   deliberately NOT asserting full-order recovery on sprinkled targets.
   Sharpness (Lemma 4f): the review's counterexample -- two
   configurations with identical D but affinely inequivalent targets --
   is reproduced to float precision, pinning that D carries order and
   nothing metric.
7. Theorem 2 (Model P): direct seeded simulation of the stated
   stochastic model -- Poisson tick chains, widths via the rank-gap
   identity, cross-checked exactly against find_radar_ticks_from_order.
   Same-slice pairs must NEVER strictly invert (pathwise claim: one
   inversion falsifies the theorem); tie rate / mean / variance match
   exp(-4 lam g) / 4 lam g / 4 lam g; arbitrary-time error rates and
   full-order failure rates are dominated by the Claim 1 / Claim 3
   bounds at parameters where those bounds are nontrivial. Verifies the
   THEOREM by simulation; no instrument realizes Model P (G2), and no
   rho-scaling claim is tested.

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
from causal_spacetime_lab.positive_control.dissimilarity import (
    profile_dissimilarity_matrix,
)
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

# Float fuzz allowance around the proved band [-1, 1). The band itself is
# exact integer arithmetic; the tolerance only covers tick times built by
# np.linspace and the atol inside the causal matrix.
BAND_TOL = 1e-9


def residuals_in_band(residuals: np.ndarray) -> np.ndarray:
    """Membership test for the proved half-open band ``[-1, 1)``.

    The asymmetry is deliberate. The lower edge is closed, so it gets a +tol
    allowance against float fuzz. The upper edge is OPEN, and there exists a
    true residual of exactly +1: a target spacetime-coincident with a tick
    (|dx| = 0, t on the grid) leaves that tick unrelated under the dt > 0
    convention -- neither predecessor nor successor -- so W = 2 against a
    predicted center of 1. A symmetric ``< 1 + tol`` check would wave that
    genuine violation through as fuzz; shaving the open edge by tol instead
    means a true +1 can never pass, at the cost of (measure-zero) legitimate
    theta within 1e-9 of 1.
    """

    residuals = np.asarray(residuals)
    return (residuals >= -1.0 - BAND_TOL) & (residuals < 1.0 - BAND_TOL)


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
    in_band = residuals_in_band(residuals)

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
    """Check 3b, direct-order form: bulk density must not move any width.

    Model D's clock is a deterministic grid appended to the scene; the bulk
    events are causally related to targets and ticks but never enter the
    bracket computation. This is a *unit invariant* of the order machinery:
    it hand-builds the chain, so it holds by the pairwise nature of
    causal_matrix_1p1() no matter what the scene builder does -- which is
    why it cannot, by itself, catch the builder ever deriving its clock
    from n_events. check_builder_density_invariance() covers that.
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


def check_builder_density_invariance(
    n_targets: int = 20,
    n_events_grid: tuple[int, ...] = (300, 900, 2700),
    scene_seed: int = 0,
    target_seed: int = 3,
) -> dict:
    """Check 3b, builder form: the scene BUILDER must not couple its clock
    to the sprinkling density.

    Scenes come from build_positive_control_scene() at a 9x density range.
    If the builder ever started deriving ticks_per_chain, chain span, or
    tick placement from n_events -- the instrument regression a density
    falsifier exists to catch -- the chains below would differ and this
    check would fail, where the direct-order form stays green by
    construction. Asserted: (i) the builder's chain worldlines are
    bit-identical across densities, and (ii) one fixed set of controlled
    targets appended to each built scene measures bit-identical widths on
    every chain.
    """

    rng = np.random.default_rng(target_seed)
    fixed_targets = np.column_stack([
        rng.uniform(-0.10, 0.10, size=n_targets),
        rng.uniform(-0.25, 0.25, size=n_targets),
    ])

    chain_coords_by_density = []
    widths_by_density = []
    for n_events in n_events_grid:
        # min_targets=0: the scene's own targets are drawn from the bulk and
        # legitimately vary with density; only the builder's chains and these
        # fixed appended targets are under test.
        config = PositiveControlSceneConfig(
            n_events=n_events, seed=scene_seed, min_targets=0
        )
        scene = build_positive_control_scene(config)
        chain_coords_by_density.append(
            np.vstack([scene.events[idx] for idx in scene.chain_index_arrays])
        )
        combined = np.vstack([scene.events, fixed_targets])
        causal = causal_matrix_1p1(combined)
        widths = np.full((n_targets, len(scene.chain_index_arrays)), np.nan)
        for row in range(n_targets):
            for column, chain in enumerate(scene.chain_index_arrays):
                bracket = find_radar_ticks_from_order(
                    causal, chain, len(scene.events) + row, scene.tick_ranks
                )
                if bracket is not None:
                    widths[row, column] = bracket[1] - bracket[0]
        widths_by_density.append(widths)

    chains_identical = all(
        np.array_equal(chain_coords_by_density[0], other)
        for other in chain_coords_by_density[1:]
    )
    widths_identical = all(
        np.array_equal(widths_by_density[0], other, equal_nan=True)
        for other in widths_by_density[1:]
    )
    n_reachable = int(np.sum(~np.isnan(widths_by_density[0])))
    return {
        "n_events_grid": list(n_events_grid),
        "chains_identical": bool(chains_identical),
        "widths_identical": bool(widths_identical),
        "n_reachable": n_reachable,
        "passed": bool(chains_identical and widths_identical and n_reachable > 0),
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


def check_coincident_tick_orphan(span: float = 1.4, ticks: int = 96) -> dict:
    """The one true band violation, pinned: a target spacetime-coincident
    with a tick.

    With |dx| = 0 and the target's time exactly on the grid, the coincident
    tick has dt = 0 and is unrelated under the dt > 0 convention -- neither
    predecessor, spacelike, nor successor -- so the partition behind
    W = N + 1 fails: W = 2 against a predicted center of 1, residual exactly
    +1, outside the half-open [-1, 1). Lemma 2 excludes this by hypothesis
    (no tick coincident with the target; measure zero for sprinkled targets).
    This check exists so the band predicate demonstrably REJECTS the case
    rather than waving it through as float fuzz.
    """

    delta = span / (ticks - 1)
    tick_times = np.linspace(-span / 2.0, span / 2.0, ticks)
    target = np.array([[float(tick_times[ticks // 3]), 0.0]])
    width = bracket_widths_ranks(target, 0.0, span, ticks)[0]
    residual = float(width - predicted_center(np.array([0.0]), delta)[0])
    rejected = not bool(residuals_in_band(np.array([residual]))[0])

    return {
        "width": float(width),
        "residual": residual,
        "band_rejects_it": rejected,
        "passed": bool(width == 2.0 and residual == 1.0 and rejected),
    }


def check_null_aligned_tick() -> dict:
    """The load-bearing convention, pinned: Lemma 2's null-aligned example.

    The identity W = N + 1 holds under the null-inclusive order
    (dt > 0 and dt^2 >= dx^2) and fails under the strict chronological
    order -- and a target whose light cone passes exactly through ticks is
    the one deterministic configuration where the two conventions measure
    different widths. Every other check samples general-position targets,
    where the conventions agree; the second review round demonstrated that
    swapping causal_matrix_1p1 for the strict relation at runtime left the
    whole suite green. This check closes that hole with the document's own
    example: target (0, 0.5) against ticks -0.75, -0.5, ..., 0.75
    (delta = 0.25, every coordinate exactly representable). Inclusive: the
    null-aligned ticks at t = -/+0.5 are the radar endpoints, W = 4 = N + 1,
    theta = -1 exactly (both interval endpoints grid-aligned -- the closed
    edge of the band). Strict: both null ticks are orphaned, W = 6,
    residual +1, outside the band. The strict width is also computed here
    as an explicit counterfactual, so the discrimination gap is recorded,
    not assumed.
    """

    span, ticks = 1.5, 7
    delta = span / (ticks - 1)
    abs_dx = 2.0 * delta
    target = np.array([[0.0, abs_dx]])

    width = bracket_widths_ranks(target, 0.0, span, ticks)[0]
    residual = float(width - predicted_center(np.array([abs_dx]), delta)[0])
    in_band = bool(residuals_in_band(np.array([residual]))[0])

    chain_events, _ = make_stationary_observer_chain_1p1(span, ticks, x=0.0)
    events = np.vstack([target, chain_events])
    dt = events[None, :, 0] - events[:, None, 0]
    dx = events[None, :, 1] - events[:, None, 1]
    strict_causal = (dt > 0.0) & (dt * dt - dx * dx > 1e-12)
    bracket = find_radar_ticks_from_order(
        strict_causal,
        np.arange(1, 1 + ticks),
        0,
        np.arange(ticks, dtype=np.float64),
    )
    strict_width = float(bracket[1] - bracket[0]) if bracket else float("nan")

    return {
        "width_inclusive": float(width),
        "residual": residual,
        "in_band": in_band,
        "width_strict_counterfactual": strict_width,
        "passed": bool(
            width == 4.0
            and residual == -1.0
            and in_band
            and strict_width == 6.0
        ),
    }


def exact_centered_profiles(
    positions: np.ndarray,
    chain_positions: np.ndarray,
    lam: float,
) -> np.ndarray:
    """Lemma 2 expectation profiles, parallax-centered: phi(x) of Lemma 4a."""

    widths = (
        2.0 * lam * np.abs(positions[:, None] - chain_positions[None, :]) + 1.0
    )
    return widths - widths.mean(axis=1, keepdims=True)


def exact_dissimilarity(
    positions: np.ndarray,
    chain_positions: np.ndarray,
    lam: float,
) -> np.ndarray:
    """Exact-model D(i,j): RMS over columns of centered-profile differences."""

    profiles = exact_centered_profiles(positions, chain_positions, lam)
    diffs = profiles[:, None, :] - profiles[None, :, :]
    return np.sqrt(np.mean(diffs * diffs, axis=2))


def anchor_decode(dissimilarity: np.ndarray) -> tuple[np.ndarray, tuple[int, int]]:
    """Lemma 4c decoder: argmax pair anchors, one anchor row sorts.

    Returns (order, anchor_pair) where order[k] is the index of the k-th
    target counted from the first anchor. Uses D alone -- no labels, no
    coordinates, no profile access.
    """

    matrix = np.asarray(dissimilarity, dtype=float)
    i, j = np.unravel_index(np.nanargmax(matrix), matrix.shape)
    return np.argsort(matrix[i], kind="stable"), (int(i), int(j))


def check_unlabeled_decoding_exact(
    seed: int = 5,
    n_targets: int = 18,
) -> dict:
    """Check 6 (exact model): Lemma 4 end to end, at R = 2, 3 and 6.

    Three assertions per observer count: (a) the gap-direction inner
    products match 4 a_k b_l / R exactly (finite differences across gap
    midpoints); (b) D is strictly Robinson -- quadrance superadditivity
    D(x,z)^2 > D(x,y)^2 + D(y,z)^2 on every ordered triple; (c) the anchor
    decoder returns the true spatial order up to global reversal, with the
    argmax pair equal to the extreme pair. R = 2 is deliberately included:
    v0.2 expected R >= 3 to be the clean hypothesis, and the proof showed
    R >= 2 suffices -- a claim worth pinning, not just stating.
    """

    rng = np.random.default_rng(seed)
    observer_sets = {
        2: np.array([-0.25, 0.25]),
        3: np.array([-0.25, 0.0, 0.25]),
        6: np.array([-0.25, -0.15, -0.05, 0.05, 0.15, 0.25]),
    }
    lam = 95.0 / 1.4  # the PC-V1 rate; any positive value proves the same
    results = {}
    all_ok = True
    for count, observers in sorted(observer_sets.items()):
        positions = np.sort(
            rng.uniform(observers[0] + 0.01, observers[-1] - 0.01, size=n_targets)
        )

        # (a) gap-direction inner products, by finite differences.
        eps = 1e-7
        gap_mids = (observers[:-1] + observers[1:]) / 2.0
        directions = []
        for mid in gap_mids:
            pair = exact_centered_profiles(
                np.array([mid - eps, mid + eps]), observers, lam
            )
            directions.append((pair[1] - pair[0]) / (2.0 * eps * 2.0 * lam))
        formula_ok = True
        for k in range(len(directions)):
            for line in range(k, len(directions)):
                a_k, b_l = k + 1, count - (line + 1)
                got = float(np.dot(directions[k], directions[line]))
                if abs(got - 4.0 * a_k * b_l / count) > 1e-5:
                    formula_ok = False

        # (b) strict Robinson on every ordered triple.
        matrix = exact_dissimilarity(positions, observers, lam)
        robinson_ok = True
        for i in range(n_targets):
            for j in range(i + 1, n_targets):
                for k in range(j + 1, n_targets):
                    lhs = matrix[i, k] ** 2
                    rhs = matrix[i, j] ** 2 + matrix[j, k] ** 2
                    if lhs <= rhs + 1e-9:
                        robinson_ok = False

        # (c) decoder: anchors are the extremes, order exact up to reversal.
        order, anchor_pair = anchor_decode(matrix)
        anchors_ok = sorted(anchor_pair) == [0, n_targets - 1]
        forward = np.array_equal(order, np.arange(n_targets))
        reverse = np.array_equal(order, np.arange(n_targets)[::-1])

        ok = formula_ok and robinson_ok and anchors_ok and (forward or reverse)
        all_ok = all_ok and ok
        results[f"R={count}"] = {
            "inner_product_formula": formula_ok,
            "strictly_robinson": robinson_ok,
            "anchors_are_extremes": anchors_ok,
            "order_recovered_up_to_reversal": bool(forward or reverse),
        }

    return {**results, "passed": bool(all_ok)}


def check_unlabeled_decoding_pipeline(seed: int = 0) -> dict:
    """Check 6 (pipeline): the Model-D decoding claim at its proved strength.

    Measured D comes from the real instrument path
    (build_positive_control_scene -> measure_bracket_echo_profiles ->
    profile_dissimilarity_matrix); exact D from the targets' true
    coordinates and the Lemma 2 expectations. Asserted, per Lemma 4e:
    (0) the shared-centering hypothesis -- every chain reaches every
    target, so the measured per-target centering uses the same column
    set as the exact comparison. PC-V1 guarantees this via scene
    validity (min_bracketing_chains = R), but the bound below is FALSE
    without it (unequal reachable sets shift row means at coordinate
    scale, not rank scale), so the check fails rather than assumes;
    (i) |D_measured - D_exact| < 4 on every defined pair; (ii) every
    anchor-row comparison whose exact-model margin exceeds 8 is ordered
    correctly by the measured decoder. Full-order recovery is NOT
    asserted: sprinkled targets may sit closer than the decodable
    separation, and on this scene a few sub-delta pairs do invert --
    that is the resolution limit, not a failure of the lemma.
    """

    config = PositiveControlSceneConfig(seed=seed)
    scene = build_positive_control_scene(config)
    profiles = measure_bracket_echo_profiles(scene)
    measured = profile_dissimilarity_matrix(profiles)

    # Lemma 4e's shared-centering hypothesis, asserted rather than assumed.
    all_reachable = bool(np.asarray(profiles.reachable).all())

    lam = (config.ticks_per_chain - 1) / config.chain_span
    positions = scene.target_coords[:, 1]
    exact = exact_dissimilarity(
        positions, np.asarray(config.chain_positions), lam
    )

    defined = np.isfinite(measured)
    perturbation = float(np.max(np.abs(measured[defined] - exact[defined])))

    order, anchor_pair = anchor_decode(measured)
    anchor = anchor_pair[0]
    rank = np.empty(len(positions), dtype=int)
    rank[order] = np.arange(len(positions))

    margin_comparisons = 0
    violations = 0
    for i in range(len(positions)):
        for j in range(len(positions)):
            if i == j:
                continue
            if exact[anchor, j] - exact[anchor, i] > 8.0:
                margin_comparisons += 1
                if rank[j] <= rank[i]:
                    violations += 1

    return {
        "n_targets": int(len(positions)),
        "shared_centering_hypothesis": all_reachable,
        "max_perturbation": perturbation,
        "perturbation_bound": 4.0,
        "margin_comparisons": margin_comparisons,
        "margin_violations": violations,
        "anchors_are_true_extremes": bool(
            sorted(anchor_pair)
            == sorted([int(np.argmin(positions)), int(np.argmax(positions))])
        ),
        "passed": bool(
            all_reachable
            and perturbation < 4.0
            and margin_comparisons > 0
            and violations == 0
        ),
    }


def check_d_spacing_sharpness() -> dict:
    """Check 6 (sharpness, Lemma 4f): D determines order, NOT spacings.

    The review's counterexample, reproduced to float precision: two
    hypothesis-satisfying configurations (R = 6, lambda = 1) with the
    SAME dissimilarity matrix whose target sets are affinely
    inequivalent -- the affine invariant (x2 - x1)/(x3 - x1) differs by
    about 0.06. The per-gap D-speeds 4 a_k b_l / R depend on the unknown
    observer layout, which absorbs spacing information while preserving
    every pairwise D. Pins that Theorem 1's positive-affine clause
    belongs to the labeled flanking decoder only: any future "decoder"
    reading D alone and emitting spacings must fail this check's premise,
    because the spacings are demonstrably not in D.
    """

    observers_a = np.array([0.0, 0.1, 0.21, 0.34, 0.44, 0.65])
    targets_a = np.array([0.02, 0.28, 0.53])
    observers_b = np.array([
        0.0, 0.02928039325682, 0.075237960253012,
        0.300251599099648, 0.311938592558539, 0.65,
    ])
    targets_b = np.array([0.02, 0.241172809035204, 0.511495131189338])

    d_a = exact_dissimilarity(targets_a, observers_a, lam=1.0)
    d_b = exact_dissimilarity(targets_b, observers_b, lam=1.0)
    max_difference = float(np.max(np.abs(d_a - d_b)))

    invariant_a = float(
        (targets_a[1] - targets_a[0]) / (targets_a[2] - targets_a[0])
    )
    invariant_b = float(
        (targets_b[1] - targets_b[0]) / (targets_b[2] - targets_b[0])
    )

    hypotheses_ok = bool(
        np.all(np.diff(observers_a) > 0)
        and np.all(np.diff(observers_b) > 0)
        and np.all((targets_a > observers_a[0]) & (targets_a < observers_a[-1]))
        and np.all((targets_b > observers_b[0]) & (targets_b < observers_b[-1]))
    )

    return {
        "max_d_difference": max_difference,
        "affine_invariant_a": invariant_a,
        "affine_invariant_b": invariant_b,
        "hypotheses_satisfied": hypotheses_ok,
        "passed": bool(
            hypotheses_ok
            and max_difference < 1e-12
            and abs(invariant_a - invariant_b) > 0.05
        ),
    }


def _poisson_chain(rng: np.random.Generator, lam: float, span: float) -> np.ndarray:
    """One realization of a rate-lam Poisson tick process on [-span/2, span/2]."""

    count = rng.poisson(lam * span)
    return np.sort(rng.uniform(-span / 2.0, span / 2.0, size=count))


def _width_by_identity(ticks: np.ndarray, t: float, d: float) -> tuple[int, bool]:
    """Rank-gap identity: W = N + 1, N = ticks in the open radar interval.

    Returns (width, reachable). Poisson ticks a.s. avoid the interval
    endpoints, so the open/closed distinction never binds in simulation.
    """

    lo = int(np.searchsorted(ticks, t - d, side="right"))
    hi = int(np.searchsorted(ticks, t + d, side="left"))
    reachable = lo > 0 and hi < ticks.size
    return hi - lo + 1, reachable


def check_model_p_theorem2(
    seed: int = 7,
    trials: int = 4000,
    order_trials: int = 1500,
) -> dict:
    """Check 7: Theorem 2 by direct simulation of the stated model.

    This verifies the THEOREM, not the instrument: Poisson chains are
    simulated here precisely because no code constructs them as an
    instrument (G2), and nothing rho-dependent is claimed. Four parts:

    - identity cross-check: widths from the rank-gap identity equal
      widths from causal_matrix_1p1 + find_radar_ticks_from_order on the
      same Poisson draws, exactly (ties Theorem 2's W to the machinery
      the rest of the harness uses);
    - same-slice pair (Claim 2): strict inversions must be ZERO -- the
      pathwise-monotonicity claim is falsified by a single one -- and
      the tie rate / mean / variance must match exp(-4 lam g) /
      4 lam g / 4 lam g within generous MC tolerances;
    - arbitrary-time pair (Claim 1): the empirical error rate must be
      dominated by exp(-2 lam g^2 / (L + g/3)), and the mean/variance
      must match the shared-region-cancellation bookkeeping (Step 2's
      symmetric-difference formula);
    - full order (Claim 3): the empirical failure rate must be dominated
      by the (n-1)-term union bound, at parameters where that bound is
      nontrivial (< 1).

    Every simulated width's reachability flag is collected and asserted
    zero-failures per section: _width_by_identity still returns the
    synthetic value N + 1 when the observable is undefined (an empty
    chain gives (1, False)), so silently discarding the flag could let
    invalid draws into the statistics under a different seed or
    parameter set. A nonzero count fails the check loudly instead.
    """

    rng = np.random.default_rng(seed)
    lam, span = 40.0, 3.0
    x0_left, x0_right = -0.3, 0.3

    # --- identity cross-check against the instrument machinery ------------
    identity_matches = True
    for _ in range(20):
        ticks = _poisson_chain(rng, lam, span)
        if ticks.size < 2:
            continue
        target_t = float(rng.uniform(-0.4, 0.4))
        target_x = float(rng.uniform(-0.2, 0.2))
        x0 = 0.25
        events = np.vstack([
            np.array([[target_t, target_x]]),
            np.column_stack([ticks, np.full(ticks.size, x0)]),
        ])
        causal = causal_matrix_1p1(events)
        bracket = find_radar_ticks_from_order(
            causal,
            np.arange(1, 1 + ticks.size),
            0,
            np.arange(ticks.size, dtype=np.float64),
        )
        width_id, reachable_id = _width_by_identity(
            ticks, target_t, abs(target_x - x0)
        )
        if bracket is None:
            if reachable_id:
                identity_matches = False
        elif not reachable_id or int(bracket[1] - bracket[0]) != width_id:
            identity_matches = False

    # --- Claim 2: same-slice pair ------------------------------------------
    xi, xj, t_common = 0.04, 0.02, 0.0
    gap = xi - xj
    deltas = np.empty(trials)
    unreachable_same = 0
    for k in range(trials):
        left = _poisson_chain(rng, lam, span)
        right = _poisson_chain(rng, lam, span)
        wi_l, ok_a = _width_by_identity(left, t_common, xi - x0_left)
        wj_l, ok_b = _width_by_identity(left, t_common, xj - x0_left)
        wi_r, ok_c = _width_by_identity(right, t_common, x0_right - xi)
        wj_r, ok_d = _width_by_identity(right, t_common, x0_right - xj)
        if not (ok_a and ok_b and ok_c and ok_d):
            unreachable_same += 1
        deltas[k] = (wi_l - wi_r) - (wj_l - wj_r)

    strict_inversions = int(np.sum(deltas < 0))
    tie_rate = float(np.mean(deltas == 0))
    tie_exact = float(np.exp(-4.0 * lam * gap))
    tie_sigma = float(np.sqrt(tie_exact * (1 - tie_exact) / trials))
    mean_ok = abs(float(deltas.mean()) - 4 * lam * gap) < 5 * np.sqrt(
        4 * lam * gap / trials
    )
    var_ok = abs(float(deltas.var()) - 4 * lam * gap) < 5 * np.sqrt(
        2.0 / trials
    ) * 4 * lam * gap + 0.5

    same_slice = {
        "strict_inversions": strict_inversions,
        "unreachable_trials": unreachable_same,
        "tie_rate": tie_rate,
        "tie_exact": tie_exact,
        "mean_matches": bool(mean_ok),
        "variance_matches": bool(var_ok),
        "tie_rate_matches": bool(abs(tie_rate - tie_exact) < 5 * tie_sigma),
    }

    # --- Claim 1: arbitrary-time pair (partial overlap on one chain) -------
    ti, tj = 0.31, -0.27
    xi2, xj2 = 0.06, 0.01
    gap2 = xi2 - xj2
    d_il, d_jl = xi2 - x0_left, xj2 - x0_left
    d_ir, d_jr = x0_right - xi2, x0_right - xj2
    length_bound = 2 * max(d_il, d_jl, d_ir, d_jr)

    def _symdiff(t1, d1, t2, d2):
        overlap = max(0.0, min(t1 + d1, t2 + d2) - max(t1 - d1, t2 - d2))
        return 2 * d1 + 2 * d2 - 2 * overlap

    predicted_var = lam * (
        _symdiff(ti, d_il, tj, d_jl) + _symdiff(ti, d_ir, tj, d_jr)
    )
    deltas2 = np.empty(trials)
    unreachable_general = 0
    for k in range(trials):
        left = _poisson_chain(rng, lam, span)
        right = _poisson_chain(rng, lam, span)
        wi_l, ok_a = _width_by_identity(left, ti, d_il)
        wj_l, ok_b = _width_by_identity(left, tj, d_jl)
        wi_r, ok_c = _width_by_identity(right, ti, d_ir)
        wj_r, ok_d = _width_by_identity(right, tj, d_jr)
        if not (ok_a and ok_b and ok_c and ok_d):
            unreachable_general += 1
        deltas2[k] = (wi_l - wi_r) - (wj_l - wj_r)

    claim1_bound = float(
        np.exp(-2 * lam * gap2**2 / (length_bound + gap2 / 3.0))
    )
    error_rate = float(np.mean(deltas2 <= 0))
    general_pair = {
        "error_rate": error_rate,
        "unreachable_trials": unreachable_general,
        "claim1_bound": claim1_bound,
        "bound_dominates": bool(error_rate <= claim1_bound),
        "mean_matches": bool(
            abs(float(deltas2.mean()) - 4 * lam * gap2)
            < 5 * np.sqrt(predicted_var / trials)
        ),
        "variance_matches": bool(
            abs(float(deltas2.var()) - predicted_var)
            < 5 * np.sqrt(2.0 / trials) * predicted_var
        ),
    }

    # --- Claim 3: full order, at parameters where the bound is < 1 ---------
    lam_order = 300.0
    positions = np.linspace(-0.15, 0.15, 6)
    times = rng.uniform(-0.25, 0.25, size=positions.size)
    g_min = float(np.min(np.diff(positions)))
    length_all = 2 * max(
        float(np.max(positions - x0_left)), float(np.max(x0_right - positions))
    )
    per_pair = np.exp(-2 * lam_order * g_min**2 / (length_all + g_min / 3.0))
    union_bound = float((positions.size - 1) * per_pair)

    failures = 0
    unreachable_order = 0
    for _ in range(order_trials):
        left = _poisson_chain(rng, lam_order, span)
        right = _poisson_chain(rng, lam_order, span)
        flanking = np.empty(positions.size)
        trial_reachable = True
        for m in range(positions.size):
            w_l, ok_l = _width_by_identity(
                left, times[m], positions[m] - x0_left
            )
            w_r, ok_r = _width_by_identity(
                right, times[m], x0_right - positions[m]
            )
            trial_reachable = trial_reachable and ok_l and ok_r
            flanking[m] = w_l - w_r
        if not trial_reachable:
            unreachable_order += 1
        if np.any(np.diff(flanking) <= 0):
            failures += 1

    full_order = {
        "failure_rate": failures / order_trials,
        "unreachable_trials": unreachable_order,
        "union_bound": union_bound,
        "bound_nontrivial": bool(union_bound < 1.0),
        "bound_dominates": bool(failures / order_trials <= union_bound),
    }

    passed = bool(
        identity_matches
        and unreachable_same == 0
        and unreachable_general == 0
        and unreachable_order == 0
        and strict_inversions == 0
        and same_slice["tie_rate_matches"]
        and same_slice["mean_matches"]
        and same_slice["variance_matches"]
        and general_pair["bound_dominates"]
        and general_pair["mean_matches"]
        and general_pair["variance_matches"]
        and full_order["bound_nontrivial"]
        and full_order["bound_dominates"]
    )
    return {
        "identity_cross_check": bool(identity_matches),
        "same_slice": same_slice,
        "general_pair": general_pair,
        "full_order": full_order,
        "passed": passed,
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
    in_band = residuals_in_band(residuals)
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
        "3b_density_invariance_direct": check_density_invariance(),
        "3b_density_invariance_builder": check_builder_density_invariance(),
        "4_centered_residue": check_centered_residue(),
        "5_coincident_tick_orphan": check_coincident_tick_orphan(),
        "5_null_aligned_tick": check_null_aligned_tick(),
        "6_unlabeled_decoding_exact": check_unlabeled_decoding_exact(),
        "6_unlabeled_decoding_pipeline": check_unlabeled_decoding_pipeline(),
        "6_d_spacing_sharpness": check_d_spacing_sharpness(),
        "7_model_p_theorem2": check_model_p_theorem2(),
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
