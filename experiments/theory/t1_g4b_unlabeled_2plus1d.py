"""G4b: what the parallax dissimilarity determines in 2+1D.

G4a settled the statements whose proofs never used the spatial
dimension, and recorded G4b -- the *unlabeled* clause -- as open,
because Lemma 4's engine is seriation and "spatial order" is a
one-dimensional notion. That diagnosis of the obstruction was right and
is unchanged: no Robinson structure survives into the plane. The
conclusion drawn from it was wrong. The 2+1D counterpart is not a
weaker statement about order, it is a STRONGER metric one, and this
harness establishes its shape.

Setup (the exact model, as in Lemma 4). Observers at `p_1..p_R` in
`R^d`, targets at `x_1..x_n`. The profile is `Phi(x)_r = |x - p_r|`,
centered as `Phi~ = Phi - mean_r Phi` (PC-V1's parallax step), and the
pipeline consumes only

    D(j,k) = || Phi~(x_j) - Phi~(x_k) || / sqrt(R),

an `n x n` matrix carrying NO observer labels and NO observer
positions. The G4b question: what does `D` alone determine?

Three findings, in the order they force themselves on you.

1. The centering fold is NOT the obstruction. One might expect
   centering to cost an observer: `Phi~(x') = Phi~(x)` iff
   `|x'-p_r| - |x-p_r| = c` for every `r`, which is the classical TDOA
   condition, and TDOA localization from three receivers is famously
   two-valued. It is not, for targets in the hull: a dense scan finds
   no twin, and `d(Phi~)/dx` has full rank `d` throughout. The
   *labeled* centered profile already determines the target.

2. `R = 3` is nevertheless never enough, and the failure is not a
   fold but a continuous flex. With `R = 3` the centered profiles lie
   in the 2-dimensional mean-zero subspace of `R^3`, so the profile
   surface fills its ambient space; `D` is then just the distance
   matrix of `n` points in a plane, invariant under any isometry of
   that plane, and the freedom is realizable by moving the observers.
   Measured: exactly 6 flexes beyond the rigid-motion gauge, at every
   `n` from 6 to 34 -- adding targets never helps. Flowing along one
   flex produces an explicit scene with the SAME `D` to machine
   precision and a demonstrably different shape.

3. `R >= 4` with enough targets is rigid. The profile surface is then
   a genuinely curved 2-surface in a `>= 3`-dimensional subspace, and
   curvature rigidifies: the nullity drops to exactly `d(d+1)/2` = 3,
   the rigid-motion gauge. So `D` determines targets AND observers up
   to Euclidean congruence -- including the absolute scale, since `D`
   is homogeneous of degree 1 in the scene. Measured thresholds:
   `R = 4` needs `n >= 11`, `R = 5` and `R = 8` need `n >= 9`,
   `R = 6` needs `n >= 8`. The frozen 2+1D instrument (8 chains,
   34-44 targets) sits well inside the rigid regime.

The 1+1D control is what makes this readable: there the harness finds
extra flexes at EVERY `(n, R)` tested, which is Lemma 4f seen through
the same instrument -- `Phi~` is piecewise linear in 1D, so the profile
"surface" is a polyline whose per-cell slopes trade off against
spacings, and no amount of data removes that. The dimensional story is
therefore the opposite of the one G4a's text anticipated: the unlabeled
observable is strictly MORE informative in 2+1D than in 1+1D.

Scope, stated plainly. These are computations in the exact model, not
written proofs, and the rigidity they establish is INFINITESIMAL --
local uniqueness generically, not global. Rank is lower semicontinuous,
so full rank at a sampled configuration does imply full rank on a dense
open set, which is why the thresholds are reported as generic; but a
global-uniqueness claim would need a separate argument, and the
measured-data counterpart (Lemma 4e's role for Lemma 4) is not
attempted here. Nothing is frozen and no gate consumes any of it.

Usage:
    python t1_g4b_unlabeled_2plus1d.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RESULTS_PATH = ROOT / "docs" / "theory" / "t1_g4b_unlabeled_results.json"

RING_RADIUS = 0.25
TARGET_RADIUS = 0.16


# --------------------------------------------------------------------
# the observable
# --------------------------------------------------------------------

def dissimilarity(theta: np.ndarray, n: int, R: int, d: int = 2):
    """``D`` over target pairs, from the flat scene vector ``theta``.

    Written with an explicit ``sqrt(sum of squares)`` rather than
    ``np.linalg.norm`` so it stays analytic under a complex step, which
    is what makes the Jacobian below exact to machine precision.
    """

    X = theta[: n * d].reshape(n, d)
    P = theta[n * d:].reshape(R, d)
    sep = X[:, None, :] - P[None, :, :]
    phi = np.sqrt(np.sum(sep * sep, axis=2))
    tilde = phi - phi.sum(axis=1, keepdims=True) / R
    rows, cols = np.triu_indices(n, k=1)
    diff = tilde[rows] - tilde[cols]
    return np.sqrt(np.sum(diff * diff, axis=1)) / np.sqrt(R)


def jacobian(theta: np.ndarray, n: int, R: int, d: int = 2, h: float = 1e-30):
    """d(D)/d(scene) by complex-step differentiation.

    A real finite difference would leave the rank decision resting on a
    spectrum blurred by subtractive cancellation; the complex step has
    no cancellation at all, so the gap between the true zeros and the
    smallest true singular value comes out around eleven orders of
    magnitude and the nullity is unambiguous.
    """

    width = dissimilarity(theta, n, R, d).size
    out = np.zeros((width, theta.size))
    for k in range(theta.size):
        probe = theta.astype(complex)
        probe[k] += 1j * h
        out[:, k] = np.imag(dissimilarity(probe, n, R, d)) / h
    return out


def rigid_motion_gauge(theta: np.ndarray, d: int) -> np.ndarray:
    """Basis of the exact invariances of ``D``: translations and (in the
    plane) rotation, applied to the WHOLE scene.

    Scale is deliberately absent. ``D`` is homogeneous of degree 1, so a
    global scaling multiplies it rather than preserving it -- which is
    why the absolute scale is recoverable and the gauge is ``d(d+1)/2``,
    not one more.
    """

    points = theta.reshape(-1, d)
    columns = []
    for axis in range(d):
        direction = np.zeros_like(points)
        direction[:, axis] = 1.0
        columns.append(direction.ravel())
    if d == 2:
        columns.append(
            np.column_stack((-points[:, 1], points[:, 0])).ravel()
        )
    return np.column_stack(columns)


def gauge_dimension(d: int) -> int:
    return d * (d + 1) // 2


def nullity(theta: np.ndarray, n: int, R: int, d: int = 2) -> tuple[int, np.ndarray]:
    """Dimension of the flex space, counted against the PARAMETER count.

    ``np.linalg.svd`` returns only ``min(m, k)`` singular values, so a
    Jacobian with more parameters than constraints hides its automatic
    null directions; the nullity has to be ``parameters - rank``.
    """

    matrix = jacobian(theta, n, R, d)
    spectrum = np.linalg.svd(matrix, compute_uv=False)
    if spectrum.size == 0 or spectrum[0] == 0.0:
        return theta.size, spectrum
    rank = int(np.sum(spectrum > spectrum[0] * 1e-9))
    return theta.size - rank, spectrum


# --------------------------------------------------------------------
# scenes
# --------------------------------------------------------------------

def observer_ring(R: int, radius: float = RING_RADIUS) -> np.ndarray:
    angles = np.linspace(0.0, 2.0 * np.pi, R, endpoint=False) + 0.11
    return radius * np.column_stack((np.cos(angles), np.sin(angles)))


def scene(n: int, R: int, seed: int, d: int = 2):
    rng = np.random.default_rng(seed)
    if d == 1:
        P = np.sort(rng.uniform(-1.0, 1.0, size=(R, 1)), axis=0)
        X = np.sort(
            rng.uniform(P.min() + 0.05, P.max() - 0.05, size=(n, 1)), axis=0
        )
        return X, P
    P = observer_ring(R)
    radius = TARGET_RADIUS * np.sqrt(rng.uniform(0.0, 1.0, size=n))
    angle = rng.uniform(0.0, 2.0 * np.pi, size=n)
    X = np.column_stack((radius * np.cos(angle), radius * np.sin(angle)))
    return X, P


def flatten(X: np.ndarray, P: np.ndarray) -> np.ndarray:
    return np.concatenate([X.ravel(), P.ravel()])


def scene_edm(theta: np.ndarray, d: int = 2) -> np.ndarray:
    """Labelled Euclidean distance matrix of the whole scene. Two scenes
    share it exactly when they are congruent (reflections included), so
    it separates a genuine ambiguity from gauge motion."""

    points = theta.reshape(-1, d)
    sep = points[:, None, :] - points[None, :, :]
    return np.sqrt(np.sum(sep * sep, axis=2))


def shape_gap(a: np.ndarray, b: np.ndarray, d: int = 2) -> float:
    return float(np.max(np.abs(scene_edm(a, d) - scene_edm(b, d))))


# --------------------------------------------------------------------
# checks
# --------------------------------------------------------------------

def check_machinery(seed: int = 5) -> dict:
    """Check 0: validate the instrument before trusting any verdict.

    Rigid motions must be exact null directions of the Jacobian, and a
    global scaling must NOT be -- if either fails, every nullity below
    is meaningless.
    """

    X, P = scene(10, 5, seed)
    theta = flatten(X, P)
    matrix = jacobian(theta, 10, 5, 2)
    base = float(np.linalg.norm(dissimilarity(theta, 10, 5, 2)))

    gauge = rigid_motion_gauge(theta, 2)
    residuals = [
        float(np.linalg.norm(matrix @ gauge[:, k]) / base)
        for k in range(gauge.shape[1])
    ]
    scale_response = float(np.linalg.norm(matrix @ theta) / base)
    return {
        "gauge_residuals": residuals,
        "gauge_directions_are_null": bool(max(residuals) < 1e-10),
        "scale_response": scale_response,
        "scale_is_not_a_gauge": bool(scale_response > 0.5),
        "passed": bool(max(residuals) < 1e-10 and scale_response > 0.5),
    }


def check_labeled_centered_map_is_injective(
    R_values: tuple[int, ...] = (3, 4, 8), grid: int = 700
) -> dict:
    """Check 1: the centering fold does not exist on the hull.

    ``Phi~(x') = Phi~(x)`` is the TDOA condition, and three-receiver
    TDOA is two-valued in general -- so one might expect centering to
    cost an observer. A dense scan of ``||Phi~(x') - Phi~(x)||`` over a
    wide box finds no second zero for targets in the hull, and the
    differential has full rank, so the labeled centered profile already
    pins the target down. Whatever defeats R = 3 below, it is not this.
    """

    rows = []
    ok = True
    for R in R_values:
        P = observer_ring(R)
        axis = np.linspace(-1.2, 1.2, grid)
        gx, gy = np.meshgrid(axis, axis, indexing="ij")
        candidates = np.column_stack((gx.ravel(), gy.ravel()))
        sep = candidates[:, None, :] - P[None, :, :]
        field = np.sqrt(np.sum(sep * sep, axis=2))
        field = field - field.mean(axis=1, keepdims=True)

        worst_rank = 2
        twins = 0
        for target in ([0.07, -0.05], [0.0, 0.0], [0.11, 0.06]):
            x = np.array(target, dtype=float)
            ref = np.linalg.norm(x[None, :] - P, axis=1)
            ref = ref - ref.mean()
            gaps = np.linalg.norm(field - ref[None, :], axis=1)
            far = np.linalg.norm(candidates - x[None, :], axis=1) > 0.05
            twins += int(np.sum(gaps[far] < 1e-4))

            step = 1e-6
            jac = np.zeros((R, 2))
            for k in range(2):
                bump = np.zeros(2)
                bump[k] = step
                plus = np.linalg.norm((x + bump)[None, :] - P, axis=1)
                minus = np.linalg.norm((x - bump)[None, :] - P, axis=1)
                jac[:, k] = (
                    (plus - plus.mean()) - (minus - minus.mean())
                ) / (2 * step)
            singular = np.linalg.svd(jac, compute_uv=False)
            worst_rank = min(worst_rank, int(np.sum(singular > 1e-8)))

        rows.append({
            "R": R,
            "twins_found": twins,
            "min_differential_rank": worst_rank,
        })
        ok = ok and twins == 0 and worst_rank == 2
    return {"rows": rows, "passed": bool(ok)}


def check_1p1d_control(
    cells: tuple[tuple[int, int], ...] = ((6, 3), (10, 4), (16, 6), (24, 8)),
) -> dict:
    """Check 2 (control): 1+1D is never rigid.

    Lemma 4f says D carries order and nothing metric, so extra flexes
    MUST appear here. If they do not, the harness is not measuring what
    it claims and nothing else in this file can be believed.
    """

    rows = []
    ok = True
    for n, R in cells:
        X, P = scene(n, R, seed=100 + n + 7 * R, d=1)
        theta = flatten(X, P)
        null, _ = nullity(theta, n, R, 1)
        extra = null - gauge_dimension(1)
        rows.append({
            "n": n, "R": R, "nullity": null,
            "gauge": gauge_dimension(1), "extra_flexes": extra,
        })
        ok = ok and extra > 0
    return {
        "rows": rows,
        "always_has_extra_flexes": bool(ok),
        "passed": bool(ok),
    }


def check_r3_is_never_rigid(
    n_values: tuple[int, ...] = (6, 10, 14, 20, 34),
) -> dict:
    """Check 3: with three observers the flex count is 6, whatever n.

    The mechanism is dimensional: three centered profiles span the
    2-dimensional mean-zero subspace of R^3, so the profile surface has
    no room to curve inside its ambient space, and D degenerates to the
    distance matrix of n coplanar points.
    """

    rows = []
    ok = True
    for n in n_values:
        X, P = scene(n, 3, seed=200 + n)
        theta = flatten(X, P)
        null, _ = nullity(theta, n, 3, 2)
        extra = null - gauge_dimension(2)
        rows.append({"n": n, "nullity": null, "extra_flexes": extra})
        ok = ok and extra == 6
    return {
        "rows": rows,
        "extra_flexes_constant_at_6": bool(ok),
        "passed": bool(ok),
    }


def _null_direction_off_gauge(theta, n, R, d=2):
    matrix = jacobian(theta, n, R, d)
    _, spectrum, vt = np.linalg.svd(matrix)
    rank = int(np.sum(spectrum > spectrum[0] * 1e-9))
    null_space = vt[rank:].T
    if null_space.shape[1] == 0:
        return None
    basis, _ = np.linalg.qr(rigid_motion_gauge(theta, d))
    cleaned = null_space - basis @ (basis.T @ null_space)
    norms = np.linalg.norm(cleaned, axis=0)
    if not np.any(norms > 1e-8):
        return None
    column = cleaned[:, int(np.argmax(norms))]
    return column / np.linalg.norm(column)


def _restore_to_level_set(theta, target, n, R, d=2, iters=60):
    for _ in range(iters):
        residual = dissimilarity(theta, n, R, d) - target
        if np.linalg.norm(residual) < 1e-14:
            break
        step, *_ = np.linalg.lstsq(jacobian(theta, n, R, d), -residual,
                                   rcond=None)
        theta = theta + step
    return theta


def check_r3_explicit_counterexample(
    n: int = 12, steps: int = 200, step_size: float = 0.01
) -> dict:
    """Check 4: build the counterexample, do not merely count it.

    Flow along a non-gauge flex, projecting back onto the exact D level
    set at every step. What comes out is a scene with the same
    dissimilarity to machine precision and a different shape -- the
    2+1D analogue of Lemma 4f, and the decisive form of "R = 3 does not
    suffice".
    """

    X, P = scene(n, 3, seed=3)
    start = flatten(X, P)
    target = dissimilarity(start, n, 3, 2)

    theta = start.copy()
    travelled = 0.0
    for _ in range(steps):
        direction = _null_direction_off_gauge(theta, n, 3, 2)
        if direction is None:
            break
        theta = _restore_to_level_set(
            theta + step_size * direction, target, n, 3, 2
        )
        travelled += step_size

    drift = float(np.max(np.abs(dissimilarity(theta, n, 3, 2) - target)))
    gap = shape_gap(theta, start)
    spread = float(np.max(scene_edm(start)[:n, :n]))
    triangle_before = scene_edm(start)[n:, n:][np.triu_indices(3, 1)]
    triangle_after = scene_edm(theta)[n:, n:][np.triu_indices(3, 1)]
    return {
        "n": n,
        "arc_length": travelled,
        "max_D_drift": drift,
        "scene_shape_gap": gap,
        "target_spread": spread,
        "observer_triangle_before": triangle_before.tolist(),
        "observer_triangle_after": triangle_after.tolist(),
        # both scenes are returned so the claim can be checked from
        # outside this function -- a counterexample nobody can recompute
        # is an assertion, not evidence
        "scene_before": start.tolist(),
        "scene_after": theta.tolist(),
        "passed": bool(drift < 1e-12 and gap > 0.02 * spread),
    }


def check_rigidity_thresholds(
    R_values: tuple[int, ...] = (4, 5, 6, 8),
    n_range: tuple[int, ...] = tuple(range(6, 16)),
    seeds: int = 4,
) -> dict:
    """Check 5: locate, per observer count, the target count from which
    the nullity equals the rigid-motion gauge on every seed."""

    rows = []
    ok = True
    for R in R_values:
        first = None
        per_n = []
        for n in n_range:
            hits = 0
            for seed in range(seeds):
                X, P = scene(n, R, seed=500 + 13 * seed + n)
                theta = flatten(X, P)
                null, _ = nullity(theta, n, R, 2)
                hits += int(null == gauge_dimension(2))
            per_n.append({"n": n, "rigid_seeds": hits, "of": seeds})
            if hits == seeds and first is None:
                first = n
        rows.append({"R": R, "first_rigid_n": first, "detail": per_n})
        ok = ok and first is not None
    return {"rows": rows, "all_reached_rigidity": bool(ok), "passed": bool(ok)}


def check_frozen_scene_is_in_the_rigid_regime() -> dict:
    """Check 6: the configuration the instrument actually uses.

    The frozen 2+1D builder places 8 chains on a non-collinear ring and
    keeps 34-44 targets, so it should sit inside the rigid regime with
    room to spare. Verified on the builder's own geometry rather than
    on a generic ring.
    """

    from causal_spacetime_lab.positive_control.scene_2d import (
        Scene2DConfig,
        build_scene_2plus1d,
        target_positions_2d,
    )

    config = Scene2DConfig(seed=0)
    built = build_scene_2plus1d(config)
    P = np.array([
        built.events[chain[0]][1:3] for chain in built.chain_index_arrays
    ])
    X = target_positions_2d(built)
    theta = flatten(X, P)
    n, R = X.shape[0], P.shape[0]
    null, spectrum = nullity(theta, n, R, 2)
    gauge = gauge_dimension(2)
    return {
        "n_targets": int(n),
        "n_observers": int(R),
        "nullity": null,
        "gauge": gauge,
        "extra_flexes": null - gauge,
        "smallest_nonzero_singular_value": float(
            spectrum[spectrum > spectrum[0] * 1e-9].min()
        ),
        "passed": bool(null == gauge),
    }


CHECKS = (
    ("machinery", check_machinery),
    ("labeled_centered_map_injective", check_labeled_centered_map_is_injective),
    ("control_1p1d_never_rigid", check_1p1d_control),
    ("r3_never_rigid", check_r3_is_never_rigid),
    ("r3_explicit_counterexample", check_r3_explicit_counterexample),
    ("rigidity_thresholds", check_rigidity_thresholds),
    ("frozen_scene_rigid", check_frozen_scene_is_in_the_rigid_regime),
)


def main() -> None:
    results: dict = {
        "scope": (
            "G4b: what D determines in 2+1D. Exact model, infinitesimal "
            "rigidity, numerical -- not written proofs, and not global "
            "uniqueness. Nothing frozen; no gate consumes this."
        ),
        "checks": {},
    }
    for name, check in CHECKS:
        outcome = check()
        results["checks"][name] = outcome
        flag = "PASS" if outcome["passed"] else "FAIL"
        print(f"[{flag}] {name}")

    all_passed = all(row["passed"] for row in results["checks"].values())
    results["all_passed"] = bool(all_passed)
    results["headline"] = {
        "r3_not_identifiable": True,
        "r3_extra_flexes": 6,
        "first_rigid_n_by_R": {
            str(row["R"]): row["first_rigid_n"]
            for row in results["checks"]["rigidity_thresholds"]["rows"]
        },
        "recovered_up_to": (
            "Euclidean congruence of the whole scene, absolute scale "
            "included (D is homogeneous of degree 1)"
        ),
        "contrast_with_1p1d": (
            "1+1D keeps extra flexes at every (n, R) tested -- Lemma 4f "
            "through the same instrument -- so the unlabeled observable "
            "is strictly more informative in 2+1D than in 1+1D"
        ),
    }
    RESULTS_PATH.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"\nall_passed = {all_passed}")
    print(f"wrote {RESULTS_PATH}")
    if not all_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
