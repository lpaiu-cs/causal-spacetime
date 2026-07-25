"""G4c: verify the proof program's LEMMAS, and run the round-4 tests.

`t1_g4c_general_dimension.py` measures the law. This module checks the
derivation in `docs/theory/t1_g4c_proof.md` step by step, which is a
different job: a conclusion can come out right while the argument behind
it is wrong, and only the second kind of error is cheap to fix later.

The derivation rests on one observation. Since

    R * D(j,k)^2 = || w_j - w_k ||^2,   w_j = M Phi(x_j),  M = I - (1/R) 1 1^T

the flex condition is EXACTLY that the profile cloud {w_j} moves as an
infinitesimal isometry of V = ker(1^T), a space of dimension m = R - 1.
Writing L for the derivative of scene -> cloud and F for the cloud's flex
space, that gives

    nullity = dim ker L + dim(Im L  cap  F),      dim F = n(m-q) + q(q+1)/2

with q the affine span of the cloud, and both closed forms drop out.

Checks 1-6 verify the lemmas as lemmas. Check 7 is the d = 1
counterexample, which matters more than it looks: it is a case where the
necessary count is satisfied and rigidity still fails, so it is standing
proof that the sufficiency step (Theorem 2b, the one open item) needs
curvature and cannot be closed by counting alone.

Checks 8-12 are round 4: out-of-sample cells for the derived closed
forms, frozen in `t1_g4c_predictions_round4.json` before this file
existed. As in the earlier rounds a missed prediction is recorded rather
than failing the run -- CI pins measurements, never predictions.

Usage:
    python t1_g4c_proof_check.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from t1_g4b_unlabeled_2plus1d import (  # noqa: E402
    flatten,
    gauge_dimension,
    jacobian,
    nullity,
    observer_shell,
    scene,
)
from t1_g4c_general_dimension import (  # noqa: E402
    check_curved,
    check_saturated,
    check_under_observed,
)

ROOT = Path(__file__).resolve().parents[2]
ROUND4_PATH = ROOT / "docs" / "theory" / "t1_g4c_predictions_round4.json"
RESULTS_PATH = ROOT / "docs" / "theory" / "t1_g4c_proof_check_results.json"

TOL = 1e-9


# --------------------------------------------------------------------
# the objects the proof talks about
# --------------------------------------------------------------------

def rank_of(matrix: np.ndarray) -> int:
    if matrix.size == 0:
        return 0
    spectrum = np.linalg.svd(matrix, compute_uv=False)
    if spectrum[0] == 0.0:
        return 0
    return int(np.sum(spectrum > spectrum[0] * TOL))


def null_basis(matrix: np.ndarray) -> np.ndarray:
    _, spectrum, vt = np.linalg.svd(matrix)
    rank = 0 if spectrum.size == 0 or spectrum[0] == 0 else int(
        np.sum(spectrum > spectrum[0] * TOL)
    )
    return vt[rank:].T


def mean_zero_basis(R: int) -> np.ndarray:
    """Orthonormal basis of V = ker(1^T) in R^R."""
    basis, _, _ = np.linalg.svd(np.eye(R) - np.ones((R, R)) / R)
    return basis[:, : R - 1]


def profile_cloud(theta: np.ndarray, n: int, R: int, d: int) -> np.ndarray:
    X = theta[: n * d].reshape(n, d)
    P = theta[n * d:].reshape(R, d)
    phi = np.linalg.norm(X[:, None, :] - P[None, :, :], axis=2)
    return phi - phi.mean(axis=1, keepdims=True)


def operator_L(theta: np.ndarray, n: int, R: int, d: int,
               h: float = 1e-30) -> np.ndarray:
    """Derivative of scene -> profile cloud, as (n*m) x d(n+R) in the
    mean-zero coordinates. Complex step, so the rank decision is not
    resting on a spectrum blurred by cancellation."""

    basis = mean_zero_basis(R)
    out = np.zeros((n * (R - 1), theta.size))
    for k in range(theta.size):
        probe = theta.astype(complex)
        probe[k] += 1j * h
        X = probe[: n * d].reshape(n, d)
        P = probe[n * d:].reshape(R, d)
        sep = X[:, None, :] - P[None, :, :]
        phi = np.sqrt(np.sum(sep * sep, axis=2))
        w = phi - phi.sum(axis=1, keepdims=True) / R
        out[:, k] = (np.imag(w) / h @ basis).ravel()
    return out


def cloud_flex_space(theta: np.ndarray, n: int, R: int, d: int):
    """Basis of F = {dw : <w_j - w_k, dw_j - dw_k> = 0}, in mean-zero
    coordinates, plus the measured affine span q."""

    coords = profile_cloud(theta, n, R, d) @ mean_zero_basis(R)
    m = coords.shape[1]
    rows = []
    for j in range(n):
        for k in range(j + 1, n):
            row = np.zeros(n * m)
            diff = coords[j] - coords[k]
            row[j * m:(j + 1) * m] = diff
            row[k * m:(k + 1) * m] = -diff
            rows.append(row)
    flex = null_basis(np.array(rows))
    centered = coords - coords.mean(axis=0, keepdims=True)
    return flex, rank_of(centered), m


def intersection_dim(first: np.ndarray, second: np.ndarray) -> int:
    """dim(span(first) cap span(second)) by the rank identity."""
    if first.size == 0 or second.size == 0:
        return 0
    return (rank_of(first) + rank_of(second)
            - rank_of(np.hstack([first, second])))


def observer_edm(P: np.ndarray) -> np.ndarray:
    """The NON-squared observer distance matrix. Squared EDMs have rank
    at most d + 2; this one is generically full rank, which is the whole
    input to Lemma E."""
    return np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)


def singular_points(P: np.ndarray) -> np.ndarray:
    """w(p_r) = M E_r -- the conical singularities of Sigma_P.

    Phi is not differentiable at an observer, so the profile surface has
    a genuine cone point over each one. A CONTINUOUS group of ambient
    isometries cannot permute a finite set, so it must fix all R of
    them, and that is what Lemma E turns into a contradiction.
    """
    E = observer_edm(P)
    return E - E.mean(axis=1, keepdims=True)


def singular_gap_jacobian(P: np.ndarray, h: float = 1e-30) -> np.ndarray:
    """d/dP of the pairwise distances between singular points.

    An ambient isometry fixes the singular points, so any flex that
    moves the observers must preserve || M(E_r - E_s) ||. Hypothesis
    (G3) is that this map is an immersion modulo congruence.
    """
    R, d = P.shape
    flat = P.ravel()
    rows, cols = np.triu_indices(R, k=1)
    out = np.zeros((rows.size, flat.size))
    for k in range(flat.size):
        probe = flat.astype(complex)
        probe[k] += 1j * h
        Q = probe.reshape(R, d)
        E = np.sqrt(np.sum((Q[:, None, :] - Q[None, :, :]) ** 2, axis=2))
        W = E - E.sum(axis=1, keepdims=True) / R
        diff = W[rows] - W[cols]
        out[:, k] = np.imag(np.sqrt(np.sum(diff * diff, axis=1))) / h
    return out


def sufficiency_upper_bound(d: int, R: int) -> int:
    """Lemma F: dim Lambda - dim K(infinity), the number of strict drops
    a nested chain of subspaces can make. Generic targets each force at
    least one, so this many targets always suffice."""
    m = R - 1
    return d * R + m * (m + 1) // 2 - d * (d + 1) // 2


def threshold_count(d: int, R: int) -> int | None:
    """N(d,R) of Theorem 2a: the proved lower bound on target count."""
    m = R - 1
    if m <= d:
        return None
    return math.ceil(
        (d * R + m * (m + 1) / 2 - d * (d + 1) / 2) / (m - d)
    )


def closed_form_nullity(d: int, R: int, n: int) -> int | None:
    """Theorem 1, valid for R <= d + 1."""
    if R - 1 > d:
        return None
    return d * R + n * (d - R + 1) + R * (R - 1) // 2


# --------------------------------------------------------------------
# lemma checks
# --------------------------------------------------------------------

LEMMA_CELLS = (
    (2, 2, 6), (2, 2, 9), (2, 3, 6), (2, 3, 10), (2, 4, 8), (2, 4, 11),
    (3, 3, 6), (3, 3, 9), (3, 4, 6), (3, 4, 12), (3, 5, 12), (3, 5, 19),
    (4, 3, 6), (4, 4, 6), (4, 5, 6),
)


def _cell(d: int, R: int, n: int):
    X, P = scene(n, R, seed=4242 + n, d=d)
    return flatten(X, P)


def check_lemma_a_reduction() -> dict:
    """Check 1: every flex of D induces a cloud motion that is an
    infinitesimal isometry, and the two spaces agree in dimension.

    Containment is tested vector by vector rather than by comparing
    dimensions alone, because two subspaces of equal dimension need not
    be the same subspace.
    """

    rows = []
    ok = True
    for d, R, n in LEMMA_CELLS:
        theta = _cell(d, R, n)
        flexes = null_basis(jacobian(theta, n, R, d))
        L = operator_L(theta, n, R, d)
        flex_cloud, _, m = cloud_flex_space(theta, n, R, d)

        worst = 0.0
        for k in range(flexes.shape[1]):
            image = L @ flexes[:, k]
            if np.linalg.norm(image) < TOL:
                continue
            residual = image - flex_cloud @ (flex_cloud.T @ image)
            worst = max(worst, float(
                np.linalg.norm(residual) / np.linalg.norm(image)
            ))
        rows.append({
            "d": d, "R": R, "n": n, "m": m,
            "flex_dim": int(flexes.shape[1]),
            "max_relative_residual": worst,
            "every_flex_lands_in_F": bool(worst < 1e-8),
        })
        ok = ok and worst < 1e-8
    return {"rows": rows, "passed": bool(ok)}


def check_lemma_b_flex_dimension() -> dict:
    """Check 2: dim F = n(m - q) + q(q+1)/2 with q the measured span."""

    rows = []
    ok = True
    for d, R, n in LEMMA_CELLS:
        theta = _cell(d, R, n)
        flex, q, m = cloud_flex_space(theta, n, R, d)
        predicted = n * (m - q) + q * (q + 1) // 2
        rows.append({
            "d": d, "R": R, "n": n, "m": m, "q": q,
            "dim_F": int(flex.shape[1]), "formula": int(predicted),
            "agrees": bool(flex.shape[1] == predicted),
        })
        ok = ok and flex.shape[1] == predicted
    return {"rows": rows, "passed": bool(ok)}


def check_lemma_c_decomposition() -> dict:
    """Check 3: nullity = dim ker L + dim(Im L cap F), on rigid cells
    as well as flexible ones."""

    rows = []
    ok = True
    for d, R, n in LEMMA_CELLS:
        theta = _cell(d, R, n)
        measured, _ = nullity(theta, n, R, d)
        L = operator_L(theta, n, R, d)
        flex, _, _ = cloud_flex_space(theta, n, R, d)
        ker_L = theta.size - rank_of(L)
        overlap = intersection_dim(L, flex)
        rows.append({
            "d": d, "R": R, "n": n,
            "measured_nullity": int(measured),
            "dim_ker_L": int(ker_L),
            "dim_ImL_cap_F": int(overlap),
            "sum": int(ker_L + overlap),
            "agrees": bool(measured == ker_L + overlap),
        })
        ok = ok and measured == ker_L + overlap
    return {"rows": rows, "passed": bool(ok)}


def check_theorem_1_closed_form() -> dict:
    """Check 4: the R <= d + 1 closed form against measured nullity."""

    rows = []
    ok = True
    for d, R, n in LEMMA_CELLS:
        if R - 1 > d:
            continue
        theta = _cell(d, R, n)
        measured, _ = nullity(theta, n, R, d)
        formula = closed_form_nullity(d, R, n)
        rows.append({
            "d": d, "R": R, "n": n,
            "measured": int(measured), "formula": int(formula),
            "agrees": bool(measured == formula),
        })
        ok = ok and measured == formula
    return {"rows": rows, "passed": bool(ok)}


#: Every (d, R) threshold ever measured. Rounds 1-3 supplied the first
#: ten; round 4 added the five marked below, each committed in advance.
MEASURED_THRESHOLDS = {
    (2, 4): 11, (2, 5): 9, (2, 6): 8, (2, 7): 8, (2, 8): 9,
    (3, 5): 19, (3, 6): 14, (3, 7): 12, (3, 8): 12,
    (4, 6): 29, (4, 7): 20, (4, 8): 17,
    (5, 7): 41, (5, 8): 27,
    (6, 8): 55,
}
ROUND4_THRESHOLD_CELLS = ((2, 7), (3, 7), (4, 7), (5, 8), (6, 8))


def check_theorem_2a_necessity() -> dict:
    """Check 5: no measured rigid cell sits below the proved bound.

    Also records whether the bound is ATTAINED, which is the separate,
    unproved half (Theorem 2b).
    """

    rows = []
    respected = True
    attained = 0
    for (d, R), measured in sorted(MEASURED_THRESHOLDS.items()):
        bound = threshold_count(d, R)
        rows.append({
            "d": d, "R": R, "m_minus_d": R - 1 - d,
            "bound_N": bound, "measured_threshold": measured,
            "bound_respected": bool(measured >= bound),
            "bound_attained": bool(measured == bound),
        })
        respected = respected and measured >= bound
        attained += int(measured == bound)
    return {
        "rows": rows,
        "bound_respected_everywhere": bool(respected),
        "cells_where_bound_is_attained": attained,
        "cells_total": len(rows),
        "passed": bool(respected),
    }


def check_d1_defeats_counting_alone() -> dict:
    """Check 6: in one dimension the necessary count is satisfied and
    rigidity still fails.

    This is the reason Theorem 2b is [CONJECTURED] rather than
    [PROVABLE]. If counting sufficed, 1+1D would be rigid, and Lemma 4f
    says it is not.
    """

    rows = []
    ok = True
    for n, R in ((6, 3), (10, 4), (16, 6), (24, 8)):
        X, P = scene(n, R, seed=100 + n + 7 * R, d=1)
        theta = flatten(X, P)
        null, _ = nullity(theta, n, R, 1)
        flex, q, m = cloud_flex_space(theta, n, R, 1)
        slack = 1 * (n + R) + flex.shape[1] - m * n
        rows.append({
            "n": n, "R": R, "m": m, "q": q,
            "nullity": int(null), "gauge": gauge_dimension(1),
            "extra_flexes": int(null - gauge_dimension(1)),
            "counting_lower_bound_on_nullity": int(slack),
            "count_permits_rigidity": bool(slack <= gauge_dimension(1)),
            "actually_rigid": bool(null == gauge_dimension(1)),
        })
        ok = ok and null > gauge_dimension(1)
    return {
        "rows": rows,
        "counting_permits_but_reality_refuses": bool(
            ok and any(r["count_permits_rigidity"] for r in rows)
        ),
        "passed": bool(ok),
    }


# --------------------------------------------------------------------
# round 4
# --------------------------------------------------------------------

def score_round4(checks: dict) -> list[dict]:
    frozen = json.loads(ROUND4_PATH.read_text(encoding="utf-8"))
    by_id = {row["id"]: row for row in frozen["predictions"]}

    p14 = checks["p14_d6_r4_slope_three"]
    p15 = checks["p15_d6_r7_saturated"]
    thresholds = checks["p16_to_p20_thresholds"]["rows"]
    found = {(row["d"], row["R"]): row["first_rigid_n"] for row in thresholds}

    outcomes = [
        {
            "id": "P14",
            "predicted": "d=6, R=4: nullity = 3n + 30 (slope 3, intercept 30)",
            "observed": (
                f"slope {p14['observed_slope']}, nullities "
                f"{[row['nullity'] for row in p14['rows']]}"
            ),
            "hit": bool(
                p14["observed_slope"] == 3
                and all(row["nullity"] == 3 * row["n"] + 30
                        for row in p14["rows"])
            ),
        },
        {
            "id": "P15",
            "predicted": "d=6, R=7: nullity 63, constant in n",
            "observed": (
                f"nullity {p15['nullity']}, constant={p15['constant_in_n']}"
            ),
            "hit": bool(p15["constant_in_n"] and p15["nullity"] == 63),
        },
    ]
    for pid, (d, R), claim in (
        ("P16", (2, 7), 8), ("P17", (3, 7), 12), ("P18", (4, 7), 20),
        ("P19", (5, 8), 27), ("P20", (6, 8), 55),
    ):
        outcomes.append({
            "id": pid,
            "predicted": f"d={d}, R={R}: first rigid n = {claim}",
            "observed": f"first rigid n = {found.get((d, R))}",
            "hit": bool(found.get((d, R)) == claim),
        })
    for row in outcomes:
        row["confidence_when_frozen"] = by_id[row["id"]]["confidence"]
    return outcomes


def check_round4_thresholds() -> dict:
    """Checks 8-12: the five out-of-sample threshold cells.

    Each window is centred on the derived value but reaches well past it
    on both sides, so a wrong prediction reports the number it actually
    found instead of reporting 'not reached'.
    """

    windows = (
        (2, 7, range(5, 15)),
        (3, 7, range(8, 19)),
        (4, 7, range(14, 27)),
        (5, 8, range(20, 34)),
        (6, 8, range(48, 63)),
    )
    rows = []
    for d, R, window in windows:
        outcome = check_curved(d, (R,), window, 4, 2000 + 100 * d + R)
        row = outcome["rows"][0]
        row["d"] = d
        row["scanned"] = [int(window[0]), int(window[-1])]
        row["derived_bound"] = threshold_count(d, R)
        rows.append(row)
    return {
        "rows": rows,
        "all_reached_rigidity": bool(all(r["first_rigid_n"] for r in rows)),
        "passed": bool(len(rows) == len(windows)),
    }


SUFFICIENCY_CELLS = tuple(
    (d, R) for d in range(2, 8) for R in range(d + 2, d + 5)
)


def check_lemma_e_inputs() -> dict:
    """Check 7: the two facts Lemma E rests on.

    (i) The non-squared observer distance matrix has full rank `R`.
    (ii) Hence the singular points `M E_r` span `V`, so their AFFINE span
         is at least `m - 1`.

    A nonzero infinitesimal isometry of `R^m` has zero set of dimension
    at most `m - 2`, because a nonzero skew matrix has even rank >= 2.
    So (ii) leaves no room: no nonzero ambient isometry can fix all the
    cone points, and `Sigma_P` has no continuous symmetry.
    """

    rows = []
    ok = True
    for d, R in SUFFICIENCY_CELLS:
        P = observer_shell(R, d)
        rank_E = rank_of(observer_edm(P))
        span = rank_of(singular_points(P))
        rows.append({
            "d": d, "R": R,
            "rank_of_observer_edm": rank_E, "needs": R,
            "linear_span_of_cone_points": span, "m": R - 1,
            "affine_span_at_least": span - 1,
            "zero_set_of_a_nonzero_isometry_at_most": R - 1 - 2,
            "no_room": bool(span - 1 > R - 3),
        })
        ok = ok and rank_E == R and span == R - 1
    return {"rows": rows, "passed": bool(ok)}


def check_g3_observer_gap_immersion() -> dict:
    """Check 8: hypothesis (G3).

    An ambient isometry fixes the cone points, so a flex that moves the
    observers must preserve every `|| M(E_r - E_s) ||`. (G3) says that
    map is an immersion modulo congruence -- rank `dR - d(d+1)/2` -- and
    then the only surviving observer motions are rigid ones.

    This is the single hypothesis the sufficiency proof still leans on
    that is not itself proved here. It is explicit, finite dimensional,
    and about one concrete map, which is a long way from where the
    conjecture stood.
    """

    rows = []
    ok = True
    for d, R in SUFFICIENCY_CELLS:
        P = observer_shell(R, d)
        measured = rank_of(singular_gap_jacobian(P))
        target = d * R - d * (d + 1) // 2
        rows.append({
            "d": d, "R": R, "jacobian_rank": measured,
            "needed_for_immersion_mod_congruence": target,
            "n_gap_constraints": R * (R - 1) // 2,
            "agrees": bool(measured == target),
        })
        ok = ok and measured == target
    return {"rows": rows, "passed": bool(ok)}


def check_sufficiency_brackets_the_threshold() -> dict:
    """Check 9: the proved lower and upper bounds, and where they meet.

    Theorem 2a gives `n >= N(d,R)`. Lemma F gives rigidity once
    `n >= dim Lambda - dim K(inf) = dR + m(m+1)/2 - d(d+1)/2`, because a
    nested chain of subspaces admits at most that many strict drops.

    At `R = d + 2` the two coincide -- `m - d = 1` makes the ceiling in
    `N` vacuous -- so the threshold is PINNED, not bracketed. That is the
    case that fixes 3+1D at 19.
    """

    rows = []
    ok = True
    for d, R in SUFFICIENCY_CELLS:
        lower = threshold_count(d, R)
        upper = sufficiency_upper_bound(d, R)
        measured = MEASURED_THRESHOLDS.get((d, R))
        pinned = lower == upper
        rows.append({
            "d": d, "R": R, "m_minus_d": R - 1 - d,
            "lower_bound_theorem_2a": lower,
            "upper_bound_lemma_f": upper,
            "bounds_coincide": bool(pinned),
            "measured": measured,
            "measured_inside_bounds": (
                None if measured is None else bool(lower <= measured <= upper)
            ),
        })
        ok = ok and lower <= upper
        if pinned:
            ok = ok and (R == d + 2)
        if measured is not None:
            ok = ok and lower <= measured <= upper
    return {
        "rows": rows,
        "pinned_cells": [
            [r["d"], r["R"]] for r in rows if r["bounds_coincide"]
        ],
        "pinned_exactly_when_R_is_d_plus_2": bool(
            all((r["R"] == d + 2) == r["bounds_coincide"]
                for r in rows for d in [r["d"]])
        ),
        "passed": bool(ok),
    }


CHECKS = (
    ("lemma_a_reduction", check_lemma_a_reduction),
    ("lemma_b_flex_dimension", check_lemma_b_flex_dimension),
    ("lemma_c_decomposition", check_lemma_c_decomposition),
    ("theorem_1_closed_form", check_theorem_1_closed_form),
    ("theorem_2a_necessity", check_theorem_2a_necessity),
    ("d1_defeats_counting_alone", check_d1_defeats_counting_alone),
    ("lemma_e_inputs", check_lemma_e_inputs),
    ("g3_observer_gap_immersion", check_g3_observer_gap_immersion),
    ("sufficiency_brackets", check_sufficiency_brackets_the_threshold),
    ("p14_d6_r4_slope_three",
     lambda: check_under_observed(6, 4, range(6, 13), 1800)),
    ("p15_d6_r7_saturated",
     lambda: check_saturated(6, 7, (8, 12, 16, 20), 1900)),
    ("p16_to_p20_thresholds", check_round4_thresholds),
)


def main() -> None:
    results: dict = {
        "scope": (
            "G4c proof program: the lemmas of docs/theory/t1_g4c_proof.md "
            "verified as lemmas, plus the round-4 out-of-sample tests of "
            "the derived closed forms. Exact model, infinitesimal "
            "rigidity. Nothing frozen; no gate consumes this."
        ),
        "predictions_read_from": str(ROUND4_PATH.relative_to(ROOT)),
        "checks": {},
    }
    for name, check in CHECKS:
        outcome = check()
        results["checks"][name] = outcome
        print(f"[{'PASS' if outcome['passed'] else 'FAIL'}] {name}")

    all_passed = all(row["passed"] for row in results["checks"].values())
    results["all_passed"] = bool(all_passed)
    results["prediction_scorecard_round4"] = score_round4(results["checks"])
    results["predictions_hit"] = sum(
        int(row["hit"]) for row in results["prediction_scorecard_round4"]
    )
    results["predictions_total"] = len(results["prediction_scorecard_round4"])

    print("\n--- frozen predictions: round 4 (derived formulas, out of sample) ---")
    for row in results["prediction_scorecard_round4"]:
        mark = "HIT " if row["hit"] else "MISS"
        print(f"[{mark}] {row['id']} ({row['confidence_when_frozen']})")
        print(f"       predicted: {row['predicted']}")
        print(f"       observed:  {row['observed']}")

    RESULTS_PATH.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(
        f"\nall_passed = {all_passed}  predictions_hit = "
        f"{results['predictions_hit']}/{results['predictions_total']}"
    )
    print(f"wrote {RESULTS_PATH}")
    if not all_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
