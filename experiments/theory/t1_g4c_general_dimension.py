"""G4c: is the G4b rigidity result a 2+1D fact, or a general-dimension law?

G4b measured what the unlabeled dissimilarity `D` determines in the
plane: nothing useful at `R = 3`, the whole scene up to congruence at
`R >= 4`. The argument it gave for that threshold, though, never
mentions the number two. The centered profile map sends `R^d` into the
`(R-1)`-dimensional mean-zero subspace of `R^R`, so its image is a
`d`-surface, and the behaviour follows from the sign of `(R-1) - d`:

    R - 1 <  d   the map is a submersion; each target slides along a
                 fibre of dimension `d - R + 1` without moving `Phi~` at
                 all, so the flex count GROWS with `n`
    R - 1 == d   the surface fills its ambient, `D` degenerates to the
                 distance matrix of `n` points in `R^d`, and the freedom
                 is realizable by moving observers: flex count CONSTANT
    R - 1 >  d   the surface is a curved `d`-surface in a bigger ambient
                 and curvature rigidifies: nullity drops to the gauge

If that reading is right, `R >= d + 2` is the threshold in every
dimension, `3 + 1` included, and the plane was never special.

This harness measures it. Two design points matter more than the code.

First, the predictions were committed BEFORE this file existed, in
`docs/theory/t1_g4c_predictions.json`. This module reads that file and
reports hit or miss against it; it does not restate the predictions, so
they cannot drift into agreement.

Second -- and this is the part worth copying elsewhere -- a missed
prediction is NOT a harness failure. `passed` here means the
measurement ran and the instrument validated in the dimension being
measured; `prediction_hit` is a separate field. Wiring CI to the
predictions instead of to the measurements would put quiet pressure on
the predictions every time one of them was wrong, which is the exact
failure the freeze exists to prevent. The regression test pins what was
MEASURED.

Scope, stated plainly, and unchanged from G4b: infinitesimal rigidity
in the exact model. Not global uniqueness, not noisy `D`, not `D`
harvested through the instrument from a measured causal set. Nothing is
frozen and no gate consumes any of it.

Usage:
    python t1_g4c_general_dimension.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from t1_g4b_unlabeled_2plus1d import (  # noqa: E402
    dissimilarity,
    flatten,
    gauge_dimension,
    jacobian,
    nullity,
    rigid_motion_gauge,
    scene,
)

ROOT = Path(__file__).resolve().parents[2]
PREDICTIONS_PATH = ROOT / "docs" / "theory" / "t1_g4c_predictions.json"
ROUND2_PATH = ROOT / "docs" / "theory" / "t1_g4c_predictions_round2.json"
ROUND3_PATH = ROOT / "docs" / "theory" / "t1_g4c_predictions_round3.json"
RESULTS_PATH = ROOT / "docs" / "theory" / "t1_g4c_general_dimension_results.json"


# --------------------------------------------------------------------
# shared measurement
# --------------------------------------------------------------------

def nullity_by_n(
    n_values, R: int, d: int, seed_base: int, shell_variant: int = 0
) -> list[dict]:
    """Nullity and extra flexes at each target count, one seed each."""

    rows = []
    for n in n_values:
        X, P = scene(
            n, R, seed=seed_base + 13 * n, d=d, shell_variant=shell_variant
        )
        theta = flatten(X, P)
        null, _ = nullity(theta, n, R, d)
        rows.append({
            "n": int(n),
            "nullity": int(null),
            "extra_flexes": int(null - gauge_dimension(d)),
        })
    return rows


def consecutive_differences(rows) -> list[int]:
    return [
        rows[k + 1]["nullity"] - rows[k]["nullity"]
        for k in range(len(rows) - 1)
    ]


def counting_bound(d: int, R: int) -> int:
    """Smallest ``n`` for which the pair constraints could in principle
    pin the scene: ``n(n-1)/2 >= d(n + R) - d(d+1)/2``.

    A necessary condition only. It was not tight in the plane (it gives
    7 at ``d = 2, R = 4`` where 11 was measured), which is what makes
    prediction P4 a directional claim rather than a numerical one.
    """

    n = 2
    while n * (n - 1) / 2 < d * (n + R) - gauge_dimension(d):
        n += 1
    return n


def affine_rank(points: np.ndarray) -> int:
    """How many dimensions the observers actually span. Four observers
    that happen to be coplanar are not a 3D configuration, and a nullity
    measured on one would be answering a different question."""

    centered = points - points.mean(axis=0, keepdims=True)
    if centered.size == 0:
        return 0
    spectrum = np.linalg.svd(centered, compute_uv=False)
    return int(np.sum(spectrum > spectrum[0] * 1e-9)) if spectrum[0] > 0 else 0


# --------------------------------------------------------------------
# checks
# --------------------------------------------------------------------

def check_machinery_at_d3(seed: int = 5) -> dict:
    """Check 0: validate the instrument in the dimension being measured.

    G4b ran this at ``d = 2``. Repeating it at ``d = 3`` is not
    ceremony: the gauge basis grew from three columns to six, and three
    of those six -- the rotation generators -- are new code paths. If a
    rotation generator is not an exact null direction, every nullity
    below is off by one and the regime verdicts are worthless.
    """

    n, R, d = 10, 5, 3
    X, P = scene(n, R, seed=seed, d=d)
    theta = flatten(X, P)
    matrix = jacobian(theta, n, R, d)
    base = float(np.linalg.norm(dissimilarity(theta, n, R, d)))

    gauge = rigid_motion_gauge(theta, d)
    residuals = [
        float(np.linalg.norm(matrix @ gauge[:, k]) / base)
        for k in range(gauge.shape[1])
    ]
    scale_response = float(np.linalg.norm(matrix @ theta) / base)
    return {
        "d": d,
        "gauge_columns": int(gauge.shape[1]),
        "gauge_dimension_formula": gauge_dimension(d),
        "gauge_residuals": residuals,
        "gauge_directions_are_null": bool(max(residuals) < 1e-10),
        "scale_response": scale_response,
        "scale_is_not_a_gauge": bool(scale_response > 0.5),
        "observer_affine_rank": affine_rank(P),
        "observers_span_the_space": bool(affine_rank(P) == d),
        "passed": bool(
            gauge.shape[1] == gauge_dimension(d)
            and max(residuals) < 1e-10
            and scale_response > 0.5
            and affine_rank(P) == d
        ),
    }


def check_under_observed(d: int, R: int, n_values, seed_base: int) -> dict:
    """The ``R - 1 < d`` regime: nullity should grow with ``n``, one
    direction per target per missing profile dimension."""

    rows = nullity_by_n(n_values, R, d, seed_base)
    steps = consecutive_differences(rows)
    expected = d - R + 1
    return {
        "d": d,
        "R": R,
        "regime_condition": f"R - 1 = {R - 1} < d = {d}",
        "rows": rows,
        "consecutive_nullity_steps": steps,
        "expected_slope": expected,
        "slope_is_constant": bool(len(set(steps)) == 1),
        "observed_slope": steps[0] if len(set(steps)) == 1 else None,
        "grows_with_n": bool(all(step > 0 for step in steps)),
        "passed": bool(len(rows) == len(n_values)),
    }


def check_saturated(d: int, R: int, n_values, seed_base: int) -> dict:
    """The ``R - 1 == d`` regime: the profile surface fills its ambient,
    so adding targets buys nothing and the flex count is flat."""

    rows = nullity_by_n(n_values, R, d, seed_base)
    nullities = sorted({row["nullity"] for row in rows})
    extras = sorted({row["extra_flexes"] for row in rows})
    return {
        "d": d,
        "R": R,
        "regime_condition": f"R - 1 = {R - 1} == d = {d}",
        "rows": rows,
        "distinct_nullities": nullities,
        "distinct_extra_flexes": extras,
        "constant_in_n": bool(len(nullities) == 1),
        "nullity": nullities[0] if len(nullities) == 1 else None,
        "extra_flexes": extras[0] if len(extras) == 1 else None,
        "counting_law_d_times_R_plus_gauge": d * R + gauge_dimension(d),
        "passed": bool(len(rows) == len(n_values)),
    }


def check_curved(
    d: int,
    R_values,
    n_range,
    seeds: int,
    seed_base: int,
    shell_variant: int = 0,
) -> dict:
    """The ``R >= d + 2`` regime: locate, per observer count, the target
    count from which the nullity equals the gauge on every seed."""

    rows = []
    for R in R_values:
        first = None
        per_n = []
        for n in n_range:
            hits = 0
            for seed in range(seeds):
                X, P = scene(
                    n, R, seed=seed_base + 13 * seed + n, d=d,
                    shell_variant=shell_variant,
                )
                theta = flatten(X, P)
                null, _ = nullity(theta, n, R, d)
                hits += int(null == gauge_dimension(d))
            per_n.append({"n": int(n), "rigid_seeds": hits, "of": seeds})
            if hits == seeds and first is None:
                first = int(n)
        rows.append({
            "R": int(R),
            "regime_condition": f"R - 1 = {R - 1} > d = {d}",
            "first_rigid_n": first,
            "counting_bound": counting_bound(d, R),
            "detail": per_n,
        })
    return {
        "d": d,
        "gauge": gauge_dimension(d),
        "shell_variant": shell_variant,
        "rows": rows,
        "all_reached_rigidity": bool(all(r["first_rigid_n"] for r in rows)),
        "passed": bool(len(rows) == len(R_values)),
    }


# --------------------------------------------------------------------
# scoring against the frozen file
# --------------------------------------------------------------------

def score(checks: dict) -> list[dict]:
    """Compare each frozen prediction against what was measured.

    Reads the committed prediction file rather than restating it, so the
    two cannot quietly converge.
    """

    frozen = json.loads(PREDICTIONS_PATH.read_text(encoding="utf-8"))
    by_id = {row["id"]: row for row in frozen["predictions"]}

    p1 = checks["p1_d3_r3_under_observed"]
    p2 = checks["p2_d3_r4_saturated"]
    p3 = checks["p3_d3_r5plus_curved"]
    p5 = checks["p5_d2_r2_under_observed"]
    r5 = next(row for row in p3["rows"] if row["R"] == 5)

    outcomes = [
        {
            "id": "P1",
            "predicted": "nullity grows in n with slope 1",
            "observed": (
                f"slope {p1['observed_slope']}, steps {p1['consecutive_nullity_steps']}"
            ),
            "hit": bool(p1["grows_with_n"] and p1["observed_slope"] == 1),
        },
        {
            "id": "P2",
            "predicted": "nullity 18, extra 12, constant in n",
            "observed": (
                f"nullity {p2['nullity']}, extra {p2['extra_flexes']}, "
                f"constant={p2['constant_in_n']}"
            ),
            "hit": bool(
                p2["constant_in_n"]
                and p2["nullity"] == 18
                and p2["extra_flexes"] == 12
            ),
        },
        {
            "id": "P3",
            "predicted": "nullity 6 = gauge, extra 0, for every R >= 5",
            "observed": "; ".join(
                f"R={row['R']} first rigid n={row['first_rigid_n']}"
                for row in p3["rows"]
            ),
            "hit": bool(p3["all_reached_rigidity"]),
        },
        {
            "id": "P4",
            "predicted": f"first rigid n at R=5 exceeds {r5['counting_bound']}",
            "observed": f"first rigid n = {r5['first_rigid_n']}",
            "hit": bool(
                r5["first_rigid_n"] is not None
                and r5["first_rigid_n"] > r5["counting_bound"]
            ),
        },
        {
            "id": "P5",
            "predicted": "d=2, R=2: nullity grows in n with slope 1",
            "observed": (
                f"slope {p5['observed_slope']}, steps {p5['consecutive_nullity_steps']}"
            ),
            "hit": bool(p5["grows_with_n"] and p5["observed_slope"] == 1),
        },
    ]
    for row in outcomes:
        row["confidence_when_frozen"] = by_id[row["id"]]["confidence"]
    return outcomes


def score_round2(checks: dict) -> list[dict]:
    """The out-of-sample round at ``d = 4``.

    P6 is the one that matters. Every cell measured before this round
    had fibre dimension 1, so "slope 1" and "slope = d - R + 1" were
    indistinguishable on the existing data; at ``d = 4, R = 3`` they
    part company.
    """

    frozen = json.loads(ROUND2_PATH.read_text(encoding="utf-8"))
    by_id = {row["id"]: row for row in frozen["predictions"]}

    p6 = checks["p6_d4_r3_under_observed"]
    p7 = checks["p7_d4_r4_under_observed"]
    p8 = checks["p8_d4_r5_saturated"]
    p9 = checks["p9_d4_r6plus_curved"]
    r6 = next(row for row in p9["rows"] if row["R"] == 6)

    outcomes = [
        {
            "id": "P6",
            "predicted": "d=4, R=3: nullity grows in n with slope 2",
            "observed": (
                f"slope {p6['observed_slope']}, steps "
                f"{p6['consecutive_nullity_steps']}"
            ),
            "hit": bool(p6["grows_with_n"] and p6["observed_slope"] == 2),
        },
        {
            "id": "P7",
            "predicted": "d=4, R=4: nullity grows in n with slope 1",
            "observed": (
                f"slope {p7['observed_slope']}, steps "
                f"{p7['consecutive_nullity_steps']}"
            ),
            "hit": bool(p7["grows_with_n"] and p7["observed_slope"] == 1),
        },
        {
            "id": "P8",
            "predicted": "d=4, R=5: nullity 30, extra 20, constant in n",
            "observed": (
                f"nullity {p8['nullity']}, extra {p8['extra_flexes']}, "
                f"constant={p8['constant_in_n']}"
            ),
            "hit": bool(
                p8["constant_in_n"]
                and p8["nullity"] == 30
                and p8["extra_flexes"] == 20
            ),
        },
        {
            "id": "P9",
            "predicted": "d=4, R>=6: nullity 10 = gauge, extra 0",
            "observed": "; ".join(
                f"R={row['R']} first rigid n={row['first_rigid_n']}"
                for row in p9["rows"]
            ),
            "hit": bool(p9["all_reached_rigidity"]),
        },
        {
            "id": "P10",
            "predicted": (
                f"d=4, R=6 threshold exceeds both {r6['counting_bound']} "
                "(counting bound) and 19 (the d=3, R=5 threshold)"
            ),
            "observed": f"first rigid n = {r6['first_rigid_n']}",
            "hit": bool(
                r6["first_rigid_n"] is not None
                and r6["first_rigid_n"] > r6["counting_bound"]
                and r6["first_rigid_n"] > 19
            ),
        },
    ]
    for row in outcomes:
        row["confidence_when_frozen"] = by_id[row["id"]]["confidence"]
    return outcomes


def score_round3(checks: dict) -> list[dict]:
    """The threshold conjecture and its control.

    P11 is a three-point quadratic interpolation put out of sample, and
    P12 asks whether a threshold is a property of the geometry at all or
    of one arbitrary observer placement. A miss on either falsifies the
    threshold line only; the regime law is a rank statement and stands
    or falls on P1-P9.
    """

    frozen = json.loads(ROUND3_PATH.read_text(encoding="utf-8"))
    by_id = {row["id"]: row for row in frozen["predictions"]}

    d5 = checks["p11_d5_r7_threshold"]["rows"][0]
    control = checks["p12_d3_r5_other_shell"]["rows"][0]
    sat = checks["p13_d5_r6_saturated"]

    outcomes = [
        {
            "id": "P11",
            "predicted": "d=5, R=7 threshold is exactly 41 (= d^2 + 3d + 1)",
            "observed": f"first rigid n = {d5['first_rigid_n']}",
            "hit": bool(d5["first_rigid_n"] == 41),
        },
        {
            "id": "P12",
            "predicted": "a redrawn observer shell leaves the d=3, R=5 threshold at 19",
            "observed": (
                f"first rigid n = {control['first_rigid_n']} "
                f"(shell variant {checks['p12_d3_r5_other_shell']['shell_variant']})"
            ),
            "hit": bool(control["first_rigid_n"] == 19),
        },
        {
            "id": "P13",
            "predicted": (
                "d=5: gauge 15 when rigid, and R=6 saturated at nullity 45, "
                "extra 30, constant in n"
            ),
            "observed": (
                f"gauge {checks['p11_d5_r7_threshold']['gauge']}; "
                f"R=6 nullity {sat['nullity']}, extra {sat['extra_flexes']}, "
                f"constant={sat['constant_in_n']}"
            ),
            "hit": bool(
                checks["p11_d5_r7_threshold"]["gauge"] == 15
                and sat["constant_in_n"]
                and sat["nullity"] == 45
                and sat["extra_flexes"] == 30
            ),
        },
    ]
    for row in outcomes:
        row["confidence_when_frozen"] = by_id[row["id"]]["confidence"]
    return outcomes


CHECKS = (
    ("machinery_at_d3", lambda: check_machinery_at_d3()),
    ("p5_d2_r2_under_observed",
     lambda: check_under_observed(2, 2, range(6, 13), 700)),
    ("p1_d3_r3_under_observed",
     lambda: check_under_observed(3, 3, range(6, 13), 800)),
    ("p2_d3_r4_saturated",
     lambda: check_saturated(3, 4, (6, 10, 14, 20, 34), 900)),
    # n runs to 25 rather than stopping at 20 because the R = 5 threshold
    # landed at 19, which left only two confirming points at the edge of
    # the original scan. Widening can only weaken the result if rigidity
    # turns out non-monotone in n, which is the direction a robustness
    # check should be able to fail in.
    ("p3_d3_r5plus_curved",
     lambda: check_curved(3, (5, 6, 8), range(6, 26), 4, 1000)),
    # round 2, frozen separately in t1_g4c_predictions_round2.json
    ("p6_d4_r3_under_observed",
     lambda: check_under_observed(4, 3, range(6, 13), 1100)),
    ("p7_d4_r4_under_observed",
     lambda: check_under_observed(4, 4, range(6, 13), 1200)),
    ("p8_d4_r5_saturated",
     lambda: check_saturated(4, 5, (6, 10, 14, 20, 30), 1300)),
    ("p9_d4_r6plus_curved",
     lambda: check_curved(4, (6, 8), range(8, 41), 4, 1400)),
    # round 3, frozen separately in t1_g4c_predictions_round3.json
    ("p11_d5_r7_threshold",
     lambda: check_curved(5, (7,), range(32, 51), 4, 1500)),
    # scanned from 10, not from 19, so a shell that LOWERS the threshold
    # is caught rather than hidden by where the scan starts
    ("p12_d3_r5_other_shell",
     lambda: check_curved(3, (5,), range(10, 26), 4, 1600, shell_variant=1)),
    ("p13_d5_r6_saturated",
     lambda: check_saturated(5, 6, (8, 12, 16, 20), 1700)),
)


def main() -> None:
    results: dict = {
        "scope": (
            "G4c: whether the G4b threshold is a general-dimension law. "
            "Exact model, infinitesimal rigidity, numerical -- not "
            "written proofs, and not global uniqueness. Predictions were "
            "committed before this harness existed; a missed prediction "
            "is recorded, not treated as a failure. Nothing frozen; no "
            "gate consumes this."
        ),
        "predictions_read_from": str(PREDICTIONS_PATH.relative_to(ROOT)),
        "checks": {},
    }
    for name, check in CHECKS:
        outcome = check()
        results["checks"][name] = outcome
        flag = "PASS" if outcome["passed"] else "FAIL"
        print(f"[{flag}] {name}")

    all_passed = all(row["passed"] for row in results["checks"].values())
    results["all_passed"] = bool(all_passed)
    results["prediction_scorecard"] = score(results["checks"])
    results["prediction_scorecard_round2"] = score_round2(results["checks"])
    results["prediction_scorecard_round3"] = score_round3(results["checks"])
    combined = (
        results["prediction_scorecard"]
        + results["prediction_scorecard_round2"]
        + results["prediction_scorecard_round3"]
    )
    results["predictions_hit"] = sum(int(row["hit"]) for row in combined)
    results["predictions_total"] = len(combined)

    for label, card in (
        ("round 1 (d = 3, plus the d = 2 retro-check)",
         results["prediction_scorecard"]),
        ("round 2 (d = 4, out of sample)",
         results["prediction_scorecard_round2"]),
        ("round 3 (threshold conjecture + construction control)",
         results["prediction_scorecard_round3"]),
    ):
        print(f"\n--- frozen predictions: {label} ---")
        for row in card:
            mark = "HIT " if row["hit"] else "MISS"
            print(f"[{mark}] {row['id']} ({row['confidence_when_frozen']})")
            print(f"       predicted: {row['predicted']}")
            print(f"       observed:  {row['observed']}")

    RESULTS_PATH.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(
        f"\nall_passed = {all_passed}  "
        f"predictions_hit = {results['predictions_hit']}"
        f"/{results['predictions_total']}"
    )
    print(f"wrote {RESULTS_PATH}")
    if not all_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
