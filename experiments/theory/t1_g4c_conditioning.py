"""A1: what the exact-model rigidity results cost at finite resolution.

G4b and G4c say the unlabeled dissimilarity determines the scene up to
congruence once `R >= d + 2`. Every one of those statements is about the
EXACT model. The instrument reads the profile at finite resolution: the
band theorem bounds each radial readout by `delta/2`, `delta` being the
tick spacing. Nothing proved so far says what that costs.

This module measures it. The quantity is

    amp := (rms position error / L) / (delta / 2L)

for a scene of diameter `L` -- how many multiples of a single readout's
own error bound you end up with in the reconstructed positions.

STATUS: DESCRIPTIVE. No preregistration, no gate, nothing frozen. `passed`
means the measurement ran and the instrument validated at that
configuration; it deliberately does NOT encode any threshold on the
measured values, because every such threshold here would have been chosen
after seeing the data. The numbers are pinned by the regression tests
instead.

What is assumed, stated up front because two of the three are choices
rather than consequences:

* **Support** -- `|readout error| <= delta/2` per profile entry. PROVED
  (the band theorem), and the only part that is.
* **Distribution** -- uniform on that interval. ASSUMED.
* **Independence** across `(target, observer)`. ASSUMED, and it is
  load-bearing: `Phi~ = M Phi` centers across observers, so any error
  component common to all observers of a target is annihilated outright.
  Check 3 measures the swing -- `amp` runs from ~2.1 under independence
  to ~0 under perfect correlation. Independence is the conservative
  choice among those tested, so the headline is pessimistic in that
  respect, but it is a choice and not a derivation.

And one more, on the inverse itself: the reconstruction here is the
LINEARISED, MINIMUM-NORM solution `J^+ dD`. Minimum-norm is optimistic --
a real nonlinear solver has no guarantee of landing on the smallest
displacement consistent with the data -- so every `amp` below is a best
case, not a typical one.

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

from t1_g4b_unlabeled_2plus1d import (  # noqa: E402
    RANK_TOLERANCE,
    flatten,
    gauge_dimension,
    jacobian,
    profile_to_dissimilarity,
    radial_profile,
    scene,
)

ROOT = Path(__file__).resolve().parents[2]
RESULTS_PATH = ROOT / "docs" / "theory" / "t1_g4c_conditioning_results.json"
G4A_RESULTS_PATH = ROOT / "docs" / "theory" / "t1_g4_2plus1d_results.json"

#: How far above ``RANK_TOLERANCE`` a margin must sit before it can be
#: read as a margin at all. See :func:`margin`.
MIN_HEADROOM = 100.0

#: The frozen instrument's own tick count, used for the headline number.
INSTRUMENT_TICKS = 96


# --------------------------------------------------------------------
# margin, with the trap it sits next to
# --------------------------------------------------------------------

def margin(theta, n: int, R: int, d: int) -> dict:
    """Smallest surviving singular value of ``J``, WITH the context that
    makes it interpretable.

    Reporting `sigma_min` alone is unsafe in precisely the regime this
    module studies. `sigma_min` is "smallest above `RANK_TOLERANCE`", not
    "smallest nonzero"; once a configuration's true margin falls under the
    cutoff it loses a dimension of rank and the function then returns the
    NEXT singular value up. Measured at `d = 2, R = 32, n = 48`: moving
    the cutoff one notch reclassifies one direction as null and the
    reported margin jumps 31-fold. The worse configuration reports the
    better number, silently.

    So: nullity comes back too, and so does how many multiples of the
    cutoff the margin is sitting on.
    """

    spectrum = np.linalg.svd(jacobian(theta, n, R, d), compute_uv=False)
    largest = float(spectrum[0])
    keep = spectrum[spectrum > largest * RANK_TOLERANCE]
    smallest = float(keep.min())
    nullity = int(theta.size - keep.size)
    headroom = (smallest / largest) / RANK_TOLERANCE
    return {
        "sigma_min": smallest,
        "sigma_max": largest,
        "ratio": smallest / largest,
        "nullity": nullity,
        "gauge": gauge_dimension(d),
        "rigid": bool(nullity == gauge_dimension(d)),
        "headroom_over_rank_tolerance": headroom,
        "margin_is_readable": bool(
            nullity == gauge_dimension(d) and headroom > MIN_HEADROOM
        ),
    }


# --------------------------------------------------------------------
# the instrument's readout error
# --------------------------------------------------------------------

def diameter(theta: np.ndarray, d: int) -> float:
    points = theta.reshape(-1, d)
    return float(np.max(
        np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)
    ))


def independent_error(rng, shape, delta):
    """The shipped model: uniform, independent per (target, observer)."""
    return rng.uniform(-delta / 2, delta / 2, size=shape)


def common_mode_error(rng, shape, delta):
    """One draw per target, shared by every observer. Annihilated by the
    centering, so this is the model under which the instrument is blind
    to its own readout error."""
    n, R = shape
    return np.repeat(rng.uniform(-delta / 2, delta / 2, size=(n, 1)), R, axis=1)


def half_common_error(rng, shape, delta):
    """Equal parts of the two, as a midpoint rather than an extreme."""
    return 0.5 * common_mode_error(rng, shape, delta) \
        + 0.5 * independent_error(rng, shape, delta)


ERROR_MODELS = {
    "independent": independent_error,
    "common_mode": common_mode_error,
    "half_common": half_common_error,
}


def amplification(theta, n: int, R: int, d: int, relative_delta: float,
                  trials: int = 30, seed: int = 3,
                  model: str = "independent") -> dict:
    """Position error in units of the pointwise readout bound.

    ``relative_delta`` is ``delta / L``, so the answer is dimensionless
    and independent of scene scale -- which matters, because ``J`` is
    homogeneous of degree zero and ``sigma_min`` is scale-free while
    ``delta`` is not. Comparing the two without dividing by the scene
    size is meaningless, and was the first thing this module got wrong.

    The ``seed`` is deliberately shared across calls at different
    ``relative_delta``: the perturbation DIRECTIONS are then identical
    and only the amplitude changes, so a difference in the answer is
    nonlinearity of ``Phi -> D`` and nothing else. Varying the seed per
    resolution would fold sampling scatter into the linearity check and
    blunt it.

    Returns the spread as well as the mean. A bare mean would sit below
    the standard this track holds its scaling exponents to.
    """

    rng = np.random.default_rng(seed)
    matrix = jacobian(theta, n, R, d)
    length = diameter(theta, d)
    delta = relative_delta * length
    phi0 = radial_profile(theta, n, R, d)
    base = profile_to_dissimilarity(phi0, R)
    draw = ERROR_MODELS[model]

    errors = []
    for _ in range(trials):
        phi = phi0 + draw(rng, phi0.shape, delta)
        residual = profile_to_dissimilarity(phi, R) - base
        step, *_ = np.linalg.lstsq(matrix, residual, rcond=None)
        errors.append(np.linalg.norm(step) / np.sqrt(n + R) / length)
    scale = relative_delta / 2
    values = np.array(errors) / scale
    return {
        "relative_delta": relative_delta,
        "model": model,
        "trials": trials,
        "amp": float(values.mean()),
        "amp_sd": float(values.std(ddof=1)),
        "amp_min": float(values.min()),
        "amp_max": float(values.max()),
        "position_error_over_L": float(values.mean()) * scale,
    }


# --------------------------------------------------------------------
# checks
# --------------------------------------------------------------------

def instrument_tick_ladder() -> dict:
    """The tick ladder, taken from G4a's tracked table rather than
    restated.

    G4a's ladder is the authority on what `delta` the instrument has at a
    given tick count. Recomputing it here from a hard-coded span would
    drift silently the first time that configuration changed, so the
    shared rungs are checked against it and the two extra rungs are
    flagged as extensions of the same rule.
    """

    g4a = json.loads(G4A_RESULTS_PATH.read_text(encoding="utf-8"))
    authoritative = {
        int(row["ticks"]): float(row["delta"])
        for row in g4a["checks"]["resolution_scaling"]["ladder"]
    }
    span = 1.4
    rungs, mismatches = [], []
    for ticks in sorted(set(authoritative) | {768, 1536}):
        delta = span / (ticks - 1)
        from_g4a = authoritative.get(ticks)
        if from_g4a is not None and abs(from_g4a - delta) > 1e-12:
            mismatches.append(ticks)
        rungs.append({
            "ticks": ticks,
            "delta": from_g4a if from_g4a is not None else delta,
            "source": "g4a table" if from_g4a is not None
                      else "extension of the same rule",
        })
    return {
        "rungs": rungs,
        "span_used_for_extensions": span,
        "shared_rungs_agree_with_g4a": not mismatches,
        "mismatched_ticks": mismatches,
    }


def frozen_instrument_scene():
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
    return flatten(X, P), X.shape[0], P.shape[0]


def check_frozen_instrument() -> dict:
    """Check 1: the configuration the instrument actually uses."""

    theta, n, R = frozen_instrument_scene()
    ladder = instrument_tick_ladder()
    length = diameter(theta, 2)
    margin_here = margin(theta, n, R, 2)

    rows = []
    for rung in ladder["rungs"]:
        relative = rung["delta"] / length
        reading = amplification(theta, n, R, 2, relative)
        rows.append({
            "ticks": rung["ticks"],
            "delta": rung["delta"],
            "delta_over_2L": rung["delta"] / (2 * length),
            **reading,
        })

    by_ticks = {row["ticks"]: row for row in rows}
    headline = by_ticks[INSTRUMENT_TICKS]
    fine = [row["amp"] for row in rows if row["ticks"] >= INSTRUMENT_TICKS]
    spread = (max(fine) - min(fine)) / min(fine)
    return {
        "n_targets": int(n), "n_observers": int(R), "scene_diameter": length,
        "margin": margin_here,
        "tick_ladder": ladder,
        "rows": rows,
        "instrument_ticks": INSTRUMENT_TICKS,
        "amp_at_the_instruments_own_resolution": headline["amp"],
        "amp_sd_there": headline["amp_sd"],
        "position_error_over_L_there": headline["position_error_over_L"],
        "relative_spread_of_amp_over_the_fine_half": spread,
        "passed": bool(
            margin_here["margin_is_readable"]
            and ladder["shared_rungs_agree_with_g4a"]
        ),
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

RESOLUTIONS = (1e-2, 1e-3, 1e-4)


def check_amplification_tracks_sigma_min() -> dict:
    """Check 2: is `sigma_min` the operative quantity, or a formality?

    Linearity is judged on the two FINE resolutions only, and that is not
    a convenience. A linearised inverse describes the reconstruction only
    while the displacement it predicts is small against the scene; at the
    ill-conditioned corners the coarsest `delta` predicts displacements of
    many scene diameters, and there the departure from linearity is the
    model correctly reporting that it no longer applies. Recorded as a
    flag rather than smoothed over, because it is the sharper statement
    about those corners: the reconstruction is not merely imprecise, it
    leaves the regime in which "determined up to congruence" says
    anything at all.
    """

    rows = []
    for label, d, R, n in CORNERS:
        X, P = scene(n, R, seed=4321 + n + R, d=d)
        theta = flatten(X, P)
        margin_here = margin(theta, n, R, d)
        readings = [
            amplification(theta, n, R, d, rel) for rel in RESOLUTIONS
        ]
        coarse, *_, finest = readings
        fine_pair = readings[-2:]
        coarse_displacement = coarse["position_error_over_L"]
        rows.append({
            "corner": label, "d": d, "R": R, "n": n,
            "margin": margin_here,
            "readings": readings,
            "amp_times_sigma_min": finest["amp"] * margin_here["sigma_min"],
            "linear_across_the_fine_pair": bool(
                abs(fine_pair[0]["amp"] - fine_pair[1]["amp"])
                / min(r["amp"] for r in fine_pair) < 0.05
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
        "product_spread_factor": max(products) / min(products),
        "all_linear_across_the_fine_pair": bool(
            all(r["linear_across_the_fine_pair"] for r in rows)
        ),
        "corners_that_leave_the_linear_regime": [
            r["corner"] for r in rows
            if r["coarse_delta_outside_the_linear_regime"]
        ],
        "passed": bool(
            all(r["margin"]["margin_is_readable"] for r in rows)
        ),
    }


def check_error_model_sensitivity() -> dict:
    """Check 3: how much of the answer is the independence assumption?

    `Phi~ = M Phi` subtracts the observer mean, so an error common to
    every observer of a target vanishes before `D` is formed. Only the
    SUPPORT of the readout error is proved; its correlation structure is
    a modelling choice, and this measures what that choice is worth.
    """

    theta, n, R = frozen_instrument_scene()
    rows = []
    for name in ERROR_MODELS:
        reading = amplification(theta, n, R, 2, 1e-3, model=name)
        rows.append({"model": name, **reading})
    by_model = {row["model"]: row["amp"] for row in rows}
    return {
        "rows": rows,
        "amp_by_model": by_model,
        "independent_is_the_conservative_choice": bool(
            by_model["independent"] >= max(
                by_model["common_mode"], by_model["half_common"]
            )
        ),
        "common_mode_is_annihilated_by_centering": bool(
            by_model["common_mode"] < 1e-6
        ),
        "passed": True,
    }


def check_observer_count_has_an_optimum(
    cells=((2, 48, (4, 5, 6, 8, 10, 12, 16, 20)),
           (3, 48, (5, 6, 7, 8, 10, 12, 14))),
    seeds: int = 3,
) -> dict:
    """Check 4: the trade-off the exact model cannot see.

    Every margin here carries its nullity, so a configuration that has
    quietly stopped being rigid cannot masquerade as a well-conditioned
    one -- which is the failure mode :func:`margin` exists to prevent.
    """

    rows = []
    for d, n, R_values in cells:
        per_R = []
        for R in R_values:
            readings = []
            for seed in range(seeds):
                X, P = scene(n, R, seed=6100 + 37 * seed + R, d=d)
                readings.append(margin(flatten(X, P), n, R, d))
            per_R.append({
                "R": int(R),
                "sigma_min": float(np.median(
                    [row["sigma_min"] for row in readings]
                )),
                "all_rigid": bool(all(row["rigid"] for row in readings)),
                "min_headroom": float(min(
                    row["headroom_over_rank_tolerance"] for row in readings
                )),
                "all_readable": bool(
                    all(row["margin_is_readable"] for row in readings)
                ),
            })
        usable = [row for row in per_R if row["all_readable"]]
        best = max(usable, key=lambda row: row["sigma_min"])
        worst = min(usable, key=lambda row: row["sigma_min"])
        rows.append({
            "d": d, "n": n, "by_R": per_R,
            "R_values_with_an_unreadable_margin": [
                row["R"] for row in per_R if not row["all_readable"]
            ],
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
        "passed": bool(all(r["by_R"] for r in rows)),
    }


def check_minimum_observer_count_is_not_the_operating_point() -> dict:
    """Check 5: what it costs to sit at the proved threshold `R = d + 2`."""

    rows = []
    for d, n in ((2, 48), (3, 48), (4, 60)):
        readings = {}
        for R in (d + 2, d + 5):
            X, P = scene(n, R, seed=6600 + 10 * d + R, d=d)
            theta = flatten(X, P)
            readings[R] = {
                "margin": margin(theta, n, R, d),
                **amplification(theta, n, R, d, 1e-3),
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
        "threshold_is_always_worse": bool(all(
            r["amp_penalty_for_sitting_at_the_threshold"] > 1 for r in rows
        )),
        "passed": bool(all(
            side["margin"]["margin_is_readable"]
            for r in rows for side in (r["at_minimum_R"], r["at_R_plus_three"])
        )),
    }


CHECKS = (
    ("frozen_instrument", check_frozen_instrument),
    ("amplification_tracks_sigma_min", check_amplification_tracks_sigma_min),
    ("error_model_sensitivity", check_error_model_sensitivity),
    ("observer_count_optimum", check_observer_count_has_an_optimum),
    ("threshold_penalty", check_minimum_observer_count_is_not_the_operating_point),
)


def main() -> None:
    results: dict = {
        "scope": (
            "A1: finite-resolution cost of the exact-model rigidity "
            "results. DESCRIPTIVE -- no preregistration, no gate, nothing "
            "frozen. 'passed' means the measurement ran and the "
            "instrument validated; it encodes no threshold on any "
            "measured value, since every such threshold would have been "
            "chosen after seeing the data. The regression tests pin the "
            "numbers."
        ),
        "assumptions": {
            "readout_support": "|error| <= delta/2 per profile entry (PROVED)",
            "readout_distribution": "uniform on that interval (ASSUMED)",
            "readout_independence": (
                "independent across (target, observer) (ASSUMED, and "
                "load-bearing -- see check error_model_sensitivity)"
            ),
            "inverse": (
                "linearised, minimum-norm J^+ dD -- optimistic, since a "
                "nonlinear solver need not find the smallest consistent "
                "displacement"
            ),
        },
        "checks": {},
    }
    for name, check in CHECKS:
        outcome = check()
        results["checks"][name] = outcome
        print(f"[{'PASS' if outcome['passed'] else 'FAIL'}] {name}")

    frozen = results["checks"]["frozen_instrument"]
    optimum = results["checks"]["observer_count_optimum"]["rows"]
    penalty = results["checks"]["threshold_penalty"]["rows"]
    sensitivity = results["checks"]["error_model_sensitivity"]
    results["all_passed"] = bool(
        all(row["passed"] for row in results["checks"].values())
    )
    results["headline"] = {
        "frozen_2plus1d_amp": frozen["amp_at_the_instruments_own_resolution"],
        "frozen_2plus1d_amp_sd": frozen["amp_sd_there"],
        "frozen_2plus1d_position_error_over_L":
            frozen["position_error_over_L_there"],
        "frozen_2plus1d_ticks": frozen["instrument_ticks"],
        "amp_spread_over_the_fine_half":
            frozen["relative_spread_of_amp_over_the_fine_half"],
        "amp_by_error_model": sensitivity["amp_by_model"],
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
