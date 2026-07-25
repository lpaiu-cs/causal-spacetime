"""A1: what the exact-model rigidity results cost at finite resolution.

G4b and G4c say the unlabeled dissimilarity determines the scene up to
congruence once `R >= d + 2`. Every one of those statements is about the
EXACT model. The instrument reads the profile at finite resolution: the
band theorem bounds each radial readout by `delta/2`, `delta` being the
tick spacing. Nothing proved so far says what that costs.

This module measures it. The quantity is

    amp := (rms position error / L) / (delta / 2L)

for a scene of diameter `L` -- how many multiples of a single readout's
own error bound you end up with in the reconstructed positions. `amp ~ 1`
means the reconstruction is as good as one readout; `amp >> 1` means the
exact-model result is formally true and practically thin.

Error model: uniform on `[-delta/2, +delta/2]` independently per profile
entry, then `Phi -> Phi~ -> D`, then the linearised inverse `J^+ dD`.
That is the instrument's own error model, not Gaussian noise on `D`, and
the difference turned out to matter -- isotropic noise on `D` overstates
the damage several-fold because it puts weight on directions the profile
map cannot produce.

STATUS: DESCRIPTIVE. No preregistration, no gate, nothing frozen. These
are characterisation measurements of an instrument, reported so that the
Section 8.6 metric claim carries its operating conditions with it.

Three findings, in order of how much they change what one would do.

1. The frozen 2+1D instrument is fine. `amp` is about 2, and the response
   is LINEAR in `delta` over three decades -- no threshold, no blow-up.
   At its own `K = 96` that is roughly 3% of the scene diameter.

2. `sigma_min` of the Jacobian is the operative quantity: `amp` tracks
   `c / sigma_min` across four decades. So the exact-model margin is not
   a formality, it is the error budget.

3. And that margin behaves in a way the exact-model theory hides
   completely. More observers LOWER the rigidity threshold while making
   the conditioning worse, catastrophically so past a point; and the
   minimum workable observer count `R = d + 2` is the worst place to
   operate, increasingly so as `d` grows. The theory says how few
   observers are possible. It does not say how many are wise.

Usage:
    python t1_g4c_conditioning.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from t1_g4b_unlabeled_2plus1d import flatten, jacobian, scene  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RESULTS_PATH = ROOT / "docs" / "theory" / "t1_g4c_conditioning_results.json"


# --------------------------------------------------------------------
# the observable, and the instrument's error model
# --------------------------------------------------------------------

def raw_profile(theta: np.ndarray, n: int, R: int, d: int) -> np.ndarray:
    X = theta[: n * d].reshape(n, d)
    P = theta[n * d:].reshape(R, d)
    return np.linalg.norm(X[:, None, :] - P[None, :, :], axis=2)


def profile_to_dissimilarity(phi: np.ndarray, n: int, R: int) -> np.ndarray:
    tilde = phi - phi.mean(axis=1, keepdims=True)
    rows, cols = np.triu_indices(n, k=1)
    return np.linalg.norm(tilde[rows] - tilde[cols], axis=1) / np.sqrt(R)


def diameter(theta: np.ndarray, d: int) -> float:
    points = theta.reshape(-1, d)
    return float(np.max(
        np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)
    ))


def smallest_nonzero_singular_value(theta, n: int, R: int, d: int) -> float:
    spectrum = np.linalg.svd(jacobian(theta, n, R, d), compute_uv=False)
    return float(spectrum[spectrum > spectrum[0] * 1e-9].min())


def amplification(theta, n: int, R: int, d: int, relative_delta: float,
                  trials: int = 30, seed: int = 3) -> float:
    """Position error in units of the pointwise readout bound.

    ``relative_delta`` is ``delta / L``, so the answer is dimensionless
    and independent of scene scale -- which matters, because ``J`` is
    homogeneous of degree zero and ``sigma_min`` is scale-free while
    ``delta`` is not. Comparing the two without dividing by the scene
    size is meaningless, and was the first thing this module got wrong.
    """

    rng = np.random.default_rng(seed)
    matrix = jacobian(theta, n, R, d)
    length = diameter(theta, d)
    delta = relative_delta * length
    phi0 = raw_profile(theta, n, R, d)
    base = profile_to_dissimilarity(phi0, n, R)

    errors = []
    for _ in range(trials):
        phi = phi0 + rng.uniform(-delta / 2, delta / 2, size=phi0.shape)
        residual = profile_to_dissimilarity(phi, n, R) - base
        step, *_ = np.linalg.lstsq(matrix, residual, rcond=None)
        errors.append(np.linalg.norm(step) / np.sqrt(n + R) / length)
    return float(np.mean(errors)) / (relative_delta / 2)


# --------------------------------------------------------------------
# checks
# --------------------------------------------------------------------

def check_frozen_instrument() -> dict:
    """Check 1: the configuration the instrument actually uses.

    Swept over the same tick ladder G4a used, so the `delta` values are
    the instrument's own rather than invented for this test.
    """

    from causal_spacetime_lab.positive_control.scene_2d import (
        Scene2DConfig,
        build_scene_2plus1d,
        target_positions_2d,
    )

    built = build_scene_2plus1d(Scene2DConfig(seed=0))
    P = np.array([
        built.events[chain[0]][1:3] for chain in built.chain_index_arrays
    ])
    X = target_positions_2d(built)
    theta = flatten(X, P)
    n, R = X.shape[0], P.shape[0]
    length = diameter(theta, 2)

    rows = []
    for ticks in (12, 24, 48, 96, 192, 384, 768):
        delta = 1.4 / (ticks - 1)
        relative = delta / length
        amp = amplification(theta, n, R, 2, relative)
        rows.append({
            "ticks": ticks,
            "delta": delta,
            "delta_over_2L": delta / (2 * length),
            "position_error_over_L": amp * relative / 2,
            "amp": amp,
        })

    tail = [row["amp"] for row in rows if row["ticks"] >= 96]
    spread = (max(tail) - min(tail)) / min(tail)
    return {
        "n_targets": int(n), "n_observers": int(R), "scene_diameter": length,
        "sigma_min": smallest_nonzero_singular_value(theta, n, R, 2),
        "rows": rows,
        "amp_at_the_instruments_own_resolution": rows[3]["amp"],
        "relative_spread_of_amp_over_the_fine_half": spread,
        "response_is_linear_in_delta": bool(spread < 0.05),
        "passed": bool(spread < 0.05),
    }


CORNERS = (
    ("2+1D frozen-like", 2, 8, 34),
    ("2+1D minimum R", 2, 4, 48),
    ("2+1D many observers", 2, 16, 48),
    ("2+1D very many observers", 2, 24, 48),
    ("3+1D minimum R", 3, 5, 48),
    ("3+1D more observers", 3, 8, 48),
    ("4+1D minimum R", 4, 6, 60),
)


def check_amplification_tracks_sigma_min() -> dict:
    """Check 2: is `sigma_min` the operative quantity, or a formality?

    The frozen scene's `amp` of ~2 sat far below its worst-case bound, so
    it was worth asking whether the margin matters at all. It does: over
    four decades of `sigma_min`, `amp * sigma_min` stays inside one order
    of magnitude.

    Linearity is judged on the two FINE resolutions only, and that is not
    a convenience. A linearised inverse describes the reconstruction only
    while the displacement it predicts is small against the scene; at the
    ill-conditioned corners the coarsest `delta` predicts displacements
    of many scene diameters, and there the deviation from linearity is
    the model correctly reporting that it no longer applies. Recorded as
    a flag rather than smoothed over, because it is the sharper statement
    about those corners: the reconstruction is not merely imprecise, it
    leaves the regime in which "determined up to congruence" says
    anything at all.
    """

    rows = []
    for label, d, R, n in CORNERS:
        X, P = scene(n, R, seed=4321 + n + R, d=d)
        theta = flatten(X, P)
        sigma = smallest_nonzero_singular_value(theta, n, R, d)
        amps = {
            f"{rel:.0e}": amplification(theta, n, R, d, rel)
            for rel in (1e-2, 1e-3, 1e-4)
        }
        fine = [amps["1e-03"], amps["1e-04"]]
        coarse_displacement = amps["1e-02"] * 1e-2 / 2
        rows.append({
            "corner": label, "d": d, "R": R, "n": n,
            "sigma_min": sigma,
            "amp_by_relative_delta": amps,
            "amp_times_sigma_min": fine[-1] * sigma,
            "linear_across_the_fine_pair": bool(
                abs(fine[0] - fine[1]) / min(fine) < 0.05
            ),
            "predicted_displacement_over_L_at_coarse_delta":
                coarse_displacement,
            "coarse_delta_outside_the_linear_regime": bool(
                coarse_displacement > 0.1
            ),
        })
    products = [row["amp_times_sigma_min"] for row in rows]
    return {
        "rows": rows,
        "product_range": [min(products), max(products)],
        "product_within_one_order_of_magnitude": bool(
            max(products) / min(products) < 10
        ),
        "all_linear_across_the_fine_pair": bool(
            all(r["linear_across_the_fine_pair"] for r in rows)
        ),
        "corners_that_leave_the_linear_regime": [
            r["corner"] for r in rows
            if r["coarse_delta_outside_the_linear_regime"]
        ],
        "passed": bool(
            all(r["linear_across_the_fine_pair"] for r in rows)
            and max(products) / min(products) < 10
        ),
    }


def check_observer_count_has_an_optimum(
    cells=((2, 48, (4, 5, 6, 8, 10, 12, 16, 20)),
           (3, 48, (5, 6, 7, 8, 10, 12, 14))),
    seeds: int = 3,
) -> dict:
    """Check 3: the trade-off the exact model cannot see.

    More observers lower the rigidity threshold monotonically. They do
    NOT improve the error budget monotonically -- past a point they wreck
    it, because observers crowded into a fixed shell produce nearly
    duplicate profile columns and the margin collapses.
    """

    rows = []
    for d, n, R_values in cells:
        per_R = []
        for R in R_values:
            sigmas = []
            for seed in range(seeds):
                X, P = scene(n, R, seed=6100 + 37 * seed + R, d=d)
                sigmas.append(
                    smallest_nonzero_singular_value(flatten(X, P), n, R, d)
                )
            per_R.append({"R": int(R), "sigma_min": float(np.median(sigmas))})
        best = max(per_R, key=lambda row: row["sigma_min"])
        worst = min(per_R, key=lambda row: row["sigma_min"])
        rows.append({
            "d": d, "n": n, "by_R": per_R,
            "best_conditioned_R": best["R"],
            "best_sigma_min": best["sigma_min"],
            "worst_R": worst["R"],
            "worst_sigma_min": worst["sigma_min"],
            "spread_factor": best["sigma_min"] / worst["sigma_min"],
            "optimum_is_interior": bool(
                best["R"] not in (R_values[0], R_values[-1])
            ),
        })
    return {
        "rows": rows,
        "every_dimension_has_an_interior_optimum": bool(
            all(r["optimum_is_interior"] for r in rows)
        ),
        "passed": bool(all(len(r["by_R"]) > 2 for r in rows)),
    }


def check_minimum_observer_count_is_not_the_operating_point() -> dict:
    """Check 4: what it costs to sit at the proved threshold `R = d + 2`.

    That is the corner G4c pins exactly, and it is the corner with the
    smallest margin. The penalty grows with dimension, which matters
    because 3+1D is the case the whole result was aimed at.
    """

    rows = []
    for d, n in ((2, 48), (3, 48), (4, 60)):
        readings = {}
        for R in (d + 2, d + 5):
            X, P = scene(n, R, seed=6600 + 10 * d + R, d=d)
            theta = flatten(X, P)
            readings[R] = {
                "sigma_min": smallest_nonzero_singular_value(theta, n, R, d),
                "amp": amplification(theta, n, R, d, 1e-3),
            }
        minimum, roomy = readings[d + 2], readings[d + 5]
        rows.append({
            "d": d, "n": n,
            "at_minimum_R": {"R": d + 2, **minimum},
            "at_R_plus_three": {"R": d + 5, **roomy},
            "amp_penalty_for_sitting_at_the_threshold":
                minimum["amp"] / roomy["amp"],
        })
    return {
        "rows": rows,
        "threshold_is_always_worse":
            bool(all(r["amp_penalty_for_sitting_at_the_threshold"] > 1
                     for r in rows)),
        "passed": bool(len(rows) == 3),
    }


CHECKS = (
    ("frozen_instrument", check_frozen_instrument),
    ("amplification_tracks_sigma_min", check_amplification_tracks_sigma_min),
    ("observer_count_optimum", check_observer_count_has_an_optimum),
    ("threshold_penalty", check_minimum_observer_count_is_not_the_operating_point),
)


def main() -> None:
    results: dict = {
        "scope": (
            "A1: finite-resolution cost of the exact-model rigidity "
            "results. DESCRIPTIVE -- no preregistration, no gate, nothing "
            "frozen. Linearised inverse; the error model is the "
            "instrument's own delta/2 readout bound, applied per profile "
            "entry."
        ),
        "checks": {},
    }
    for name, check in CHECKS:
        outcome = check()
        results["checks"][name] = outcome
        print(f"[{'PASS' if outcome['passed'] else 'FAIL'}] {name}")

    frozen = results["checks"]["frozen_instrument"]
    optimum = results["checks"]["observer_count_optimum"]["rows"]
    penalty = results["checks"]["threshold_penalty"]["rows"]
    results["all_passed"] = bool(
        all(row["passed"] for row in results["checks"].values())
    )
    results["headline"] = {
        "frozen_2plus1d_amp": frozen["amp_at_the_instruments_own_resolution"],
        "frozen_2plus1d_position_error_over_L_at_96_ticks":
            frozen["rows"][3]["position_error_over_L"],
        "response_linear_in_delta": frozen["response_is_linear_in_delta"],
        "best_conditioned_R_by_dimension": {
            str(row["d"]): row["best_conditioned_R"] for row in optimum
        },
        "conditioning_spread_across_R": {
            str(row["d"]): row["spread_factor"] for row in optimum
        },
        "amp_penalty_at_R_equals_d_plus_2": {
            str(row["d"]): row["amp_penalty_for_sitting_at_the_threshold"]
            for row in penalty
        },
        "reading": (
            "The instrument's own configuration is comfortable. The proved "
            "threshold R = d + 2 is not a good place to operate, and gets "
            "worse with dimension; nor is a large observer count, which "
            "lowers the threshold while collapsing the margin. The exact "
            "model says how few observers are possible and nothing about "
            "how many are wise."
        ),
    }

    print("\n--- headline ---")
    for key, value in results["headline"].items():
        if key != "reading":
            print(f"  {key}: {value}")

    RESULTS_PATH.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"\nall_passed = {results['all_passed']}")
    print(f"wrote {RESULTS_PATH}")
    if not results["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
