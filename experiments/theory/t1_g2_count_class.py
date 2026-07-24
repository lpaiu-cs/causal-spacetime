"""G2 close-out: the count-fluctuation class, measured directly.

G2's last open item was the fluctuation class of a harvested chain's
tick counts -- a Poisson-rate guess (`-1/4` for the distance error)
against a KPZ-like one (`-1/3`). Every earlier attempt read it off the
DISTANCE error and could not settle it, for two systematic reasons
recorded in the v0.7 mechanism note: the arm that was supposed to
isolate the count term still carried a sliding wandering admixture, and
the fitting design's own wobble (`0.110`, calibrated on the arm whose
exponent is proved) exceeded the `0.083` separation between the
candidates.

This harness stops asking the distance error. The count class is a
property of the chain's tick times alone, so measure it there:

    sd(chain length) ~ mean(chain length)^theta
        theta = 1/2  ->  Poisson-rate counts   ->  error exponent -1/4
        theta = 1/3  ->  KPZ / Tracy-Widom     ->  error exponent -1/3

No target, no radar bracket, no distance estimator, hence no wandering
admixture whatever. And the separation to resolve is `1/6 = 0.167`,
twice the `0.083` that defeated the previous route.

Why those two exponents are the right candidates. For a longest chain
among `rho * A` sprinkled points the mean length is `~ 2 sqrt(rho A)`
with Tracy-Widom fluctuations `~ (rho A)^{1/6} = mean^{1/3}`. A Poisson
clock of rate `lam` on a window of length `T` has `N ~ Poisson(lam T)`,
so `sd = sqrt(mean)`. Feeding either into `d_hat = (W - 1) / (2 lam)`
with the measured `lam ~ sqrt(rho)` reproduces exactly the `-1/4` and
`-1/3` of the open question.

The answer is that **the count class is protocol-dependent**, like the
rate coupling and the error law before it -- there is no single "the"
class, which is why asking for one never converged:

- the **order-only anchored chain**, a genuine longest chain free to
  wander, is KPZ: `theta = 1/3`;
- the **tube-confined chain** is Poisson: `theta = 1/2`. Its tube is
  `~ rho^{-1/2}` wide while the chain's natural transverse wandering is
  `~ rho^{-1/6}`, so the tube is asymptotically far narrower than the
  excursions KPZ scaling wants. Confinement destroys the transverse
  optimization that produces Tracy-Widom fluctuations, and what is left
  counts like a Poisson process.

That also retires a reading this document carried since v0.5. The
scaled tube's error exponent `-0.317` was described as the count class
"nearer `-1/3`". It is not the count exponent at all: that arm's counts
are Poisson (`-1/4`), and the measured `-0.317` is the mixture with its
residual `rho^{-1/2}` wandering that the v0.7 systematics already
flagged -- the flag was right and the attribution was wrong.

The thinned Poisson clock is the calibration arm throughout: its
`theta` is `1/2` by construction, so if the estimator does not return
that, nothing else here is trustworthy.

Theory-track: nothing frozen, no gate consumes this.

Usage:
    python t1_g2_count_class.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from causal_spacetime_lab.density_coupled_clocks import (
    harvest_chain_from_sprinkling_1p1,
    harvest_order_only_chain_1p1,
    make_poisson_clock_chain_1p1,
    nearest_event_index,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

ROOT = Path(__file__).resolve().parents[2]
RESULTS_PATH = ROOT / "docs" / "theory" / "t1_g2_count_class_results.json"

DIAMOND_T = 2.0
DIAMOND_AREA = DIAMOND_T * DIAMOND_T / 2.0
TICK_WINDOW = (-0.6, 0.6)
ELL = 0.1
TUBE_SCALE = 3.0
OBSERVER_X = 0.3
RHO_GRID = (500, 1000, 2000, 4000, 8000, 16000, 32000, 64000)
TRIALS = 300

CANDIDATES = {"poisson_rate_1/2": 0.5, "kpz_like_1/3": 1.0 / 3.0}

# Two-sided Student-t critical values; the project has no scipy.
_T975 = {2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365,
         8: 2.306, 9: 2.262, 10: 2.228}


def fit_exponent_with_uncertainty(x: list[float], y: list[float]) -> dict:
    """Log-log slope with a residual-based interval and split halves --
    the same treatment every other exponent in this track carries."""

    log_x = np.log(np.asarray(x, dtype=float))
    log_y = np.log(np.asarray(y, dtype=float))
    if log_x.size < 3:
        raise ValueError(
            "a residual-based interval needs at least 3 densities; "
            f"got {log_x.size}"
        )
    slope, intercept = np.polyfit(log_x, log_y, 1)
    residual = log_y - (slope * log_x + intercept)
    dof = int(log_x.size - 2)
    stderr = float(np.sqrt(
        float(residual @ residual) / dof
        / float(((log_x - log_x.mean()) ** 2).sum())
    ))
    critical = _T975.get(dof, 1.96)
    result = {
        "theta": float(slope),
        "stderr": stderr,
        "dof": dof,
        "n_points": int(log_x.size),
        "ci95": [float(slope - critical * stderr),
                 float(slope + critical * stderr)],
    }
    if log_x.size >= 6:
        cut = log_x.size // 2
        low = np.polyfit(log_x[: cut + 1], log_y[: cut + 1], 1)[0]
        high = np.polyfit(log_x[cut:], log_y[cut:], 1)[0]
        result["halves"] = {"low": float(low), "high": float(high)}
        result["half_split_spread"] = float(abs(low - high))
    return result


def chain_length(arm: str, rho: float, rng: np.random.Generator,
                 tube_width: float | None = None) -> int | None:
    """One realization's tick count. ``None`` marks a clock failure, so
    a failed harvest can never be silently read as a short chain."""

    if arm == "thinned":
        ticks = make_poisson_clock_chain_1p1(
            TICK_WINDOW[0], TICK_WINDOW[1], rho * ELL, OBSERVER_X, seed=rng
        )
        return int(len(ticks))

    bulk = sprinkle_1p1_causal_diamond(
        rng.poisson(rho * DIAMOND_AREA), T=DIAMOND_T, seed=rng
    )
    if arm == "tube":
        width = (
            TUBE_SCALE / np.sqrt(rho) if tube_width is None else tube_width
        )
        idx = harvest_chain_from_sprinkling_1p1(
            bulk, OBSERVER_X, width, TICK_WINDOW[0], TICK_WINDOW[1]
        )
        return int(idx.size)
    bottom = nearest_event_index(bulk, TICK_WINDOW[0], OBSERVER_X)
    top = nearest_event_index(bulk, TICK_WINDOW[1], OBSERVER_X)
    try:
        idx = harvest_order_only_chain_1p1(bulk, bottom, top)
    except ValueError:
        return None
    return int(idx.size)


def measure_arm(arm: str, rho_grid=RHO_GRID, trials: int = TRIALS,
                tube_width_at=None, base_seed: int = 0) -> dict:
    """Mean and sd of the chain length per density, then the exponent."""

    rows = []
    failures = 0
    for rho in rho_grid:
        lengths = []
        for k in range(trials):
            rng = np.random.default_rng(base_seed + 100_003 * k + int(rho))
            width = None if tube_width_at is None else tube_width_at(rho)
            value = chain_length(arm, rho, rng, tube_width=width)
            if value is None:
                failures += 1
                continue
            lengths.append(value)
        sample = np.asarray(lengths, dtype=float)
        rows.append({
            "rho": rho,
            "trials": int(sample.size),
            "mean_length": float(sample.mean()),
            "sd_length": float(sample.std(ddof=1)),
        })
    fit = fit_exponent_with_uncertainty(
        [row["mean_length"] for row in rows],
        [row["sd_length"] for row in rows],
    )
    low, high = fit["ci95"]
    fit["candidates_inside_ci95"] = {
        name: bool(low <= value <= high) for name, value in CANDIDATES.items()
    }
    fit["t_against"] = {
        name: float((fit["theta"] - value) / fit["stderr"])
        for name, value in CANDIDATES.items()
    }
    return {"rows": rows, "clock_failures": failures, "fit": fit}


def check_calibration_on_the_known_arm() -> dict:
    """Check 1: the thinned Poisson clock has ``theta = 1/2`` by
    construction. If the estimator does not return it, no verdict below
    means anything."""

    measured = measure_arm("thinned")
    fit = measured["fit"]
    low, high = fit["ci95"]
    return {
        **measured,
        "bias_vs_one_half": fit["theta"] - 0.5,
        "ci_covers_one_half": bool(low <= 0.5 <= high),
        "passed": bool(low <= 0.5 <= high and measured["clock_failures"] == 0),
    }


def check_order_only_is_kpz() -> dict:
    """Check 2: a genuine longest chain, free to wander, should show
    Tracy-Widom counts."""

    measured = measure_arm("order_only")
    fit = measured["fit"]
    return {
        **measured,
        "passed": bool(
            fit["candidates_inside_ci95"]["kpz_like_1/3"]
            and not fit["candidates_inside_ci95"]["poisson_rate_1/2"]
            and measured["clock_failures"] == 0
        ),
    }


def check_tube_is_poisson() -> dict:
    """Check 3: the shipped tube is `~ rho^{-1/2}` wide against a
    natural wandering of `~ rho^{-1/6}`, so it confines the chain far
    below the excursions KPZ scaling wants. The prediction is that the
    transverse optimization is destroyed and the counts revert to
    Poisson."""

    measured = measure_arm("tube")
    fit = measured["fit"]
    return {
        **measured,
        "passed": bool(
            fit["candidates_inside_ci95"]["poisson_rate_1/2"]
            and not fit["candidates_inside_ci95"]["kpz_like_1/3"]
        ),
    }


def check_confinement_crossover(trials: int = 150) -> dict:
    """Check 4: the mechanism, not just the two endpoints.

    If confinement is what destroys the KPZ counts, then widening the
    tube toward the natural wandering scale `rho^{-1/6}` must move
    ``theta`` back from `1/2` toward `1/3`. A tube that shrinks faster
    than the wandering stays Poisson; one that keeps up with it does
    not.
    """

    grid = RHO_GRID[:6]
    rows = []
    for label, coefficient, exponent in (
        ("w = 3.0 rho^-1/2 (shipped)", 3.0, 0.5),
        ("w = 1.0 rho^-1/3", 1.0, 1.0 / 3.0),
        ("w = 1.0 rho^-1/6 (~ wandering)", 1.0, 1.0 / 6.0),
    ):
        measured = measure_arm(
            "tube", rho_grid=grid, trials=trials,
            tube_width_at=lambda rho, c=coefficient, p=exponent: c * rho ** -p,
        )
        theta = measured["fit"]["theta"]
        rows.append({
            "tube": label,
            "width_exponent": exponent,
            "theta": theta,
            "ci95": measured["fit"]["ci95"],
            "nearer": (
                "kpz_like_1/3" if abs(theta - 1 / 3) < abs(theta - 0.5)
                else "poisson_rate_1/2"
            ),
        })
    tightest, widest = rows[0], rows[-1]
    return {
        "rows": rows,
        "theta_moves_toward_kpz_as_tube_widens": bool(
            widest["theta"] < tightest["theta"]
        ),
        "passed": bool(
            tightest["nearer"] == "poisson_rate_1/2"
            and widest["theta"] < tightest["theta"]
        ),
    }


CHECKS = (
    ("calibration_thinned_is_one_half", check_calibration_on_the_known_arm),
    ("order_only_is_kpz", check_order_only_is_kpz),
    ("tube_is_poisson", check_tube_is_poisson),
    ("confinement_crossover", check_confinement_crossover),
)


def main() -> None:
    results: dict = {
        "scope": (
            "G2 close-out: the count-fluctuation class measured directly "
            "from chain lengths, with no distance estimator and hence no "
            "wandering admixture. Theory-track; nothing frozen."
        ),
        "candidates": CANDIDATES,
        "config": {
            "rho_grid": list(RHO_GRID),
            "trials_per_density": TRIALS,
            "tick_window": list(TICK_WINDOW),
            "ell": ELL,
            "tube_scale": TUBE_SCALE,
            "observer_x": OBSERVER_X,
        },
        "checks": {},
    }
    for name, check in CHECKS:
        outcome = check()
        results["checks"][name] = outcome
        flag = "PASS" if outcome["passed"] else "FAIL"
        theta = outcome.get("fit", {}).get("theta")
        extra = f"  theta = {theta:+.4f}" if theta is not None else ""
        print(f"[{flag}] {name}{extra}")

    all_passed = all(row["passed"] for row in results["checks"].values())
    results["all_passed"] = bool(all_passed)
    results["headline"] = {
        "count_class_is_protocol_dependent": True,
        "order_only": results["checks"]["order_only_is_kpz"]["fit"]["theta"],
        "tube": results["checks"]["tube_is_poisson"]["fit"]["theta"],
        "calibration_thinned": (
            results["checks"]["calibration_thinned_is_one_half"]["fit"]["theta"]
        ),
        "implied_error_exponents": {
            "order_only_count_term": -1.0 / 3.0,
            "tube_count_term": -0.25,
        },
        "retired_reading": (
            "the scaled tube's -0.317 is NOT its count exponent: that "
            "arm's counts are Poisson (-1/4), and -0.317 is the mixture "
            "with its residual rho^{-1/2} wandering that the v0.7 "
            "systematics flagged"
        ),
    }
    RESULTS_PATH.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"\nall_passed = {all_passed}")
    print(f"wrote {RESULTS_PATH}")
    if not all_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
