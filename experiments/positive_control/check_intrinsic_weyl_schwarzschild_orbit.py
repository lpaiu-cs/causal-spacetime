"""Independent post-completion audit of the circular-orbit Schwarzschild test.

Nothing here calls the runner's statistics, gates or evaluation code.

- Statistics are rebuilt from the saved pure fits in component space: the 3x3
  electric and magnetic tensors of the four-index basis, Cholesky solves, and a
  covariance mapped to the 18 tensor entries. Quantiles come from closed-form
  survival functions by bisection, not from scipy.stats.
- The truth comes from the Riemann tensor in areal Schwarzschild coordinates,
  by finite differences of closed-form Christoffel symbols, contracted with a
  tetrad built from the circular orbit's conserved energy and angular momentum.
  Its frame in each fit comes from an O(3) alignment of all landmarks, not from
  the two-landmark rule.
- A sample of the archived order bits is recomputed with the Clairaut solver of
  s1_schwarzschild_cost in areal coordinates, which shares no code with the
  generator; flat cases use the Minkowski closed form.
- The first measurement events are tested against the proper-volume measure,
  with interval bounds from the Clairaut solver and the areal volume element at
  the generator's quadrature nodes, and the boundary of those nodes is checked
  with the same solver.
- Refits run in fresh single-threaded processes and are compared after the
  axis-sign gauge that the inverse leaves free is aligned (the calibration audit
  showed that raw fields can differ by such signs alone).

Every `posthoc` entry was added after the outcomes were known and decides nothing.
"""

import argparse
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from check_intrinsic_weyl_bombelli_benchmark import (
    ELECTRIC,
    GRAM2,
    MAGNETIC,
    MIRROR,
    close,
    distance,
    framed,
    norm,
    statistic,
    tensors,
)
from check_intrinsic_weyl_weak_accuracy import quantile10, read, sha
from replay_intrinsic_weyl_threads import THREAD_VARIABLES, differences, write_json

ROOT = Path(__file__).resolve().parents[2]
LEVELS = (1048576, 4194304, 16777216, 67108864)
KINDS = ("orbit", "orbit", "rest", "rest", "flat", "flat")
SEED, FALSE_POSITIVE = 2026112200, 1e-3
# Fixed in the protocol (Sect. 9) before the run.
FROZEN_REFITS = ((0, 4194304), (1, 16777216))
MASS, APEX, PROPER_TIME, SPEED = 1.0, 6.0, 3.0, 0.5
EPS3 = np.zeros((3, 3, 3))
for _i, _j, _k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
    EPS3[_i, _j, _k], EPS3[_i, _k, _j] = 1.0, -1.0


# -- quantiles without scipy.stats -----------------------------------------------------
def survival5(x):
    """Chi-square survival with five degrees of freedom, closed form."""
    return math.erfc(math.sqrt(x / 2)) + math.sqrt(2 * x / math.pi) * math.exp(
        -x / 2
    ) * (1 + x / 3)


def quantile5(probability):
    lo, hi = 0.0, 128.0
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if survival5(mid) > 1 - probability else (lo, mid)
    return (lo + hi) / 2


def normal_quantile(probability):
    lo, hi = 0.0, 10.0
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (
            (mid, hi)
            if 0.5 * math.erfc(mid / math.sqrt(2)) > 1 - probability
            else (lo, mid)
        )
    return (lo + hi) / 2


# -- truth in areal coordinates --------------------------------------------------------
def areal_christoffel(x, mass):
    """Gamma^a_bc of Schwarzschild in (t, r, theta, phi), closed form."""
    _, r, theta, _ = x
    f = 1 - 2 * mass / r
    g = np.zeros((4, 4, 4))
    g[0, 0, 1] = g[0, 1, 0] = mass / (r * r * f)
    g[1, 0, 0] = mass * f / (r * r)
    g[1, 1, 1] = -mass / (r * r * f)
    g[1, 2, 2] = -r * f
    g[1, 3, 3] = -r * f * math.sin(theta) ** 2
    g[2, 1, 2] = g[2, 2, 1] = 1 / r
    g[2, 3, 3] = -math.sin(theta) * math.cos(theta)
    g[3, 1, 3] = g[3, 3, 1] = 1 / r
    g[3, 2, 3] = g[3, 3, 2] = math.cos(theta) / math.sin(theta)
    return g


def areal_metric(x, mass):
    _, r, theta, _ = x
    f = 1 - 2 * mass / r
    return np.diag([-f, 1 / f, r * r, (r * math.sin(theta)) ** 2])


def areal_riemann(x, mass, step=1e-3):
    """R_abcd from central differences (fourth order) of the closed-form Gammas."""
    x = np.asarray(x, float)
    gamma = areal_christoffel(x, mass)
    d = np.zeros((4, 4, 4, 4))  # d[c, a, b, e] = d_c Gamma^a_be
    for c in (1, 2):  # the metric depends on r and theta only
        e = np.eye(4)[c] * step
        d[c] = (
            -areal_christoffel(x + 2 * e, mass)
            + 8 * areal_christoffel(x + e, mass)
            - 8 * areal_christoffel(x - e, mass)
            + areal_christoffel(x - 2 * e, mass)
        ) / (12 * step)
    up = np.zeros((4, 4, 4, 4))
    for a, b, c, e in np.ndindex(4, 4, 4, 4):
        up[a, b, c, e] = (
            d[c, a, b, e]
            - d[e, a, b, c]
            + gamma[a, c, :] @ gamma[:, b, e]
            - gamma[a, e, :] @ gamma[:, b, c]
        )
    return np.einsum("af,fbcd->abcd", areal_metric(x, mass), up)


def areal_frame(mass, radius, moving):
    """Rows (u, r-hat, motion, z-hat) at theta = pi/2, phi = 0, z-hat = -theta-hat."""
    f = 1 - 2 * mass / radius
    if moving:
        # Circular geodesic: E = f u^t, L = r^2 u^phi with the standard values.
        energy = f / math.sqrt(1 - 3 * mass / radius)
        angular = math.sqrt(mass * radius) / math.sqrt(1 - 3 * mass / radius)
        u = np.array([energy / f, 0.0, 0.0, angular / radius**2])
    else:
        u = np.array([1 / math.sqrt(f), 0.0, 0.0, 0.0])
    g = areal_metric([0.0, radius, math.pi / 2, 0.0], mass)
    # The unit vector in the (t, phi) plane orthogonal to u, along +phi.
    candidate = np.array([0.0, 0.0, 0.0, 1.0])
    candidate = candidate + (candidate @ g @ u) * u
    candidate /= math.sqrt(candidate @ g @ candidate)
    return np.array(
        [
            u,
            [0.0, math.sqrt(f), 0.0, 0.0],
            candidate,
            [0.0, 0.0, -1 / radius, 0.0],
        ]
    )


def areal_truth(mass, radius, proper_time, moving):
    """T^2 E and T^2 B in the frame (radial, motion, z) from areal coordinates."""
    frame = areal_frame(mass, radius, moving)
    g = areal_metric([0.0, radius, math.pi / 2, 0.0], mass)
    orthonormal = frame @ g @ frame.T
    tensor = proper_time**2 * np.einsum(
        "ae,bf,cg,dh,efgh->abcd",
        frame,
        frame,
        frame,
        frame,
        areal_riemann([0.0, radius, math.pi / 2, 0.0], mass),
    )
    electric = tensor[0, 1:, 0, 1:]
    magnetic = 0.5 * np.einsum("ikl,klj->ij", EPS3, tensor[1:, 1:, 0, 1:])
    return electric, magnetic, float(abs(orthonormal - np.diag([-1, 1, 1, 1])).max())


def coefficients_of(electric, magnetic):
    """Least squares on the 3x3 parts against the tensor-basis parts."""
    design = np.column_stack(
        [np.r_[ELECTRIC[k].ravel(), MAGNETIC[k].ravel()] for k in range(10)]
    )
    target = np.r_[electric.ravel(), magnetic.ravel()]
    value, *_ = np.linalg.lstsq(design, target, rcond=None)
    if np.max(abs(design @ value - target)) > 1e-12 * max(1, np.max(abs(target))):
        raise ArithmeticError("Not an algebraic Weyl pair")
    return value


def aligned_truth(fit, normal, electric, magnetic):
    """Truth in the fitted frame by O(3) alignment of all interior landmarks."""
    points = np.array(fit["coordinates_inferred_from_counts"])
    points[:, 1] *= np.sign(points[1, 1])
    points[:, 2] *= np.sign(points[2, 2])
    left, _, right = np.linalg.svd(normal[1:-1, 1:].T @ points[1:-1, 1:])
    rotation = left @ right
    residual = np.c_[normal[:, :1], normal[:, 1:] @ rotation] - points
    # B is a pseudo-tensor: an improper alignment flips its sign.
    turned = (
        rotation.T @ electric @ rotation,
        np.linalg.det(rotation) * rotation.T @ magnetic @ rotation,
    )
    return coefficients_of(*turned), float(np.sqrt((residual**2).sum(axis=1).mean()))


# -- statistics in component space -----------------------------------------------------
def component_statistics(c, s):
    """Q_B, |E|^2 - |B|^2 debiased with its plug-in Gaussian sigma, E.B likewise."""
    lower = np.linalg.cholesky(s[5:, 5:])
    q_magnetic = float(np.sum(np.linalg.solve(lower, c[5:]) ** 2))
    jacobian = np.vstack(
        [ELECTRIC.reshape(10, 9).T, MAGNETIC.reshape(10, 9).T]
    )  # 18 x 10
    mean = jacobian @ c
    cov = jacobian @ s @ jacobian.T
    out = {"q_magnetic": q_magnetic}
    for name, form in (
        ("invariant", np.diag(np.r_[np.ones(9), -np.ones(9)])),
        (
            "electric_dot_magnetic",
            0.5
            * np.block([[np.zeros((9, 9)), np.eye(9)], [np.eye(9), np.zeros((9, 9))]]),
        ),
    ):
        value = float(mean @ form @ mean - np.trace(form @ cov))
        variance = float(
            4 * mean @ form @ cov @ form @ mean + 2 * np.trace(form @ cov @ form @ cov)
        )
        out[name] = value
        out[name + "_se"] = math.sqrt(variance)
    return out


def verdicts(rows):
    """The three frozen gates and the primary claim, restated, bad news first.

    rows are dicts in case order; a failed row has valid_fit False. A valid
    control exceedance (or a valid reversed interval) is reported even when
    another row failed; only then does a failure make a gate INCOMPLETE.
    """
    q10, q5 = quantile10(1 - FALSE_POSITIVE), quantile5(1 - FALSE_POSITIVE)
    if [r["kind"] for r in rows] != list(KINDS):
        raise ValueError("All six assigned cases are required, in order")
    good = [r for r in rows if r["valid_fit"]]
    result = {}
    used = [r for r in rows if r["kind"] in ("orbit", "flat")]
    if any(r["q"] > q10 for r in good if r["kind"] == "flat"):
        result["field"] = "CONTROL_FAILED"
    elif not all(r["valid_fit"] for r in used):
        result["field"] = "INCOMPLETE"
    else:
        hit = [r["q"] > q10 for r in used if r["kind"] == "orbit"]
        result["field"] = (
            "DETECTED" if all(hit) else "MIXED" if any(hit) else "NOT_DETECTED"
        )
    if any(r["q_magnetic"] > q5 for r in good if r["kind"] != "orbit"):
        result["magnetic"] = "CONTROL_FAILED"
    elif len(good) < len(rows):
        result["magnetic"] = "INCOMPLETE"
    else:
        hit = [r["q_magnetic"] > q5 for r in rows if r["kind"] == "orbit"]
        result["magnetic"] = (
            "MAGNETIC_DETECTED" if all(hit) else "MIXED" if any(hit) else "NOT_DETECTED"
        )
    orbit = [r for r in rows if r["kind"] == "orbit"]
    if any(r["invariant_interval_95"][1] < 0 for r in orbit if r["valid_fit"]):
        result["unequal_norms"] = "REVERSED"
    elif not all(r["valid_fit"] for r in orbit):
        result["unequal_norms"] = "INCOMPLETE"
    else:
        hit = [r["invariant_interval_95"][0] > 0 for r in orbit]
        result["unequal_norms"] = (
            "UNEQUAL" if all(hit) else "MIXED" if any(hit) else "NOT_RESOLVED"
        )
    return result


def primary(level):
    """The frozen headline: all three gates pass at the top level."""
    top = level[str(LEVELS[-1])]
    if "INCOMPLETE" in top.values():
        return "INCOMPLETE"
    return (
        "SUPPORTED"
        if (top["field"], top["magnetic"], top["unequal_norms"])
        == ("DETECTED", "MAGNETIC_DETECTED", "UNEQUAL")
        else "NOT_SUPPORTED"
    )


def audit_row(folder, case, n, recorded, normal, electric, magnetic, z95):
    pure = read(folder / f"case_{case}_n{n}_pure.json")
    full = read(folder / f"case_{case}_n{n}.json")
    archive = folder / f"case_{case}_n{n}.npz"
    fit = pure["fit"]
    c, s, valid = framed(fit)
    q = statistic(fit)
    extra = component_statistics(c, s)
    interval = [
        extra["invariant"] - z95 * extra["invariant_se"],
        extra["invariant"] + z95 * extra["invariant_se"],
    ]
    largest = float(np.linalg.eigvals(s @ GRAM2).real.max())
    radius = (math.sqrt(quantile10(0.95) * largest) + 0.02 * norm(c)) / 0.98
    if normal is not None:
        independent, residual = aligned_truth(fit, normal, electric, magnetic)
    else:
        independent, residual = None, None
    detection, evaluation = recorded["detection"], recorded.get("evaluation")
    # A failed truth step keeps the gate statistics; compare with our own truth.
    truth = (
        np.asarray(evaluation["truth_coefficients_xy_positive"])
        if evaluation
        else independent
    )
    size = float(math.sqrt(np.sum(electric**2) + np.sum(magnetic**2)))
    e_hat, b_hat = tensors(c)
    row = {
        "case": case,
        "kind": KINDS[case],
        "measurement_n": n,
        "valid_fit": valid,
        "q": q,
        "q_magnetic": extra["q_magnetic"],
        "invariant": extra["invariant"],
        "invariant_se": extra["invariant_se"],
        "invariant_interval_95": interval,
        "electric_dot_magnetic": extra["electric_dot_magnetic"],
        "electric_dot_magnetic_se": extra["electric_dot_magnetic_se"],
        "variance_times_N": n * largest,
        "radius_95": radius,
        "error_mod_reflection": distance(c, truth) if truth is not None else None,
        "error_with_independent_truth": (
            distance(c, independent) if independent is not None else None
        ),
        "truth_vectors_differ_by": (
            distance(truth, independent)
            if truth is not None and independent is not None
            else None
        ),
        "landmark_coordinate_rms_residual": residual,
        "electric_eigenvalues": np.linalg.eigvalsh(e_hat).tolist(),
        "magnetic_eigenvalues": np.linalg.eigvalsh(b_hat).tolist(),
        "electric_norm": float(np.linalg.norm(e_hat)),
        "magnetic_norm": float(np.linalg.norm(b_hat)),
        "blind_boundary_kept": pure["oracle_modules_before_return"] == []
        and [Path(p) for p in pure["opened_archives"]] == [archive.resolve()]
        and fit == full["full"]["fit"],
        "archive_hash_matches": sha(archive)
        == recorded["archive_sha256"]
        == pure["archive_sha256"]
        == fit["order_archive_sha256"],
    }
    agree = [
        close(q, detection["q"]),
        close(extra["q_magnetic"], detection["q_magnetic"]),
        close(extra["invariant"], detection["invariant"]),
        close(extra["invariant_se"], detection["invariant_se"]),
        all(
            close(a, b)
            for a, b in zip(interval, detection["invariant_interval_95"], strict=True)
        ),
        close(extra["electric_dot_magnetic"], detection["electric_dot_magnetic"]),
        close(n * largest, detection["variance_times_N"]),
        np.allclose(row["electric_eigenvalues"], detection["electric_eigenvalues"]),
        np.allclose(row["magnetic_eigenvalues"], detection["magnetic_eigenvalues"]),
    ]
    row["recorded_statistics_reproduced"] = all(bool(v) for v in agree)
    row["evaluation_recorded"] = evaluation is not None
    if evaluation is None:
        agree = []
    agree += (
        [
            close(radius, evaluation["confidence_radius_95"]),
            close(
                row["error_mod_reflection"],
                evaluation["full_tidal_error_mod_reflection"],
            ),
            bool(row["error_mod_reflection"] <= radius)
            == evaluation["nominal_ball_contains_truth"],
            close(norm(truth), size) and close(evaluation["truth_norm"], size, 1e-6),
            # A wrong frame or magnetic parity would move the truth by its own norm.
            size == 0
            or (
                row["truth_vectors_differ_by"] is not None
                and row["truth_vectors_differ_by"] <= 0.1 * size
            ),
        ]
        if evaluation is not None
        else []
    )
    row["recorded_values_reproduced"] = row["recorded_statistics_reproduced"] and all(
        bool(v) for v in agree
    )
    # Post hoc, descriptive. least_squares.cost is half the whitened chi-square.
    freedom = fit["measured_intervals"] - fit["nuisance_parameters"] - 10
    row["posthoc"] = {
        "fit_chi2": 2 * fit["curved_cost"],
        "fit_degrees_of_freedom": freedom,
        "fit_chi2_z": (2 * fit["curved_cost"] - freedom) / math.sqrt(2 * freedom),
    }
    if truth is not None:
        candidates = (truth, MIRROR @ truth)
        row["posthoc"].update(
            noncentrality_of_the_truth=min(
                float(t @ np.linalg.solve(s, t)) for t in candidates
            ),
            magnetic_noncentrality_of_the_truth=min(
                float(t[5:] @ np.linalg.solve(s[5:, 5:], t[5:])) for t in candidates
            ),
            q_about_the_truth=min(
                float((c - t) @ np.linalg.solve(s, c - t)) for t in candidates
            ),
        )
        if size:
            row["relative_error"] = row["error_mod_reflection"] / size
            row["relative_radius"] = radius / size
    return row


# -- order bits and measure ------------------------------------------------------------
def areal(rho, mass):
    return rho * (1 + mass / (2 * rho)) ** 2


def independent_time(first, second, mass):
    """Clairaut solver of s1_schwarzschild_cost in areal coordinates: (time, bound)."""
    from s1_schwarzschild_cost import flight_time

    angle = math.atan2(np.linalg.norm(np.cross(first, second)), first @ second)
    return flight_time(
        areal(np.linalg.norm(first), mass),
        areal(np.linalg.norm(second), mass),
        angle,
        mass,
        1e-12,
    )


def measure_check(geometry, mass, events, grid=(12, 16, 12)):
    """Moments against the areal volume element and Clairaut interval bounds."""
    (p, wp, low), (q, wq, high) = geometry.quadrature(*grid)
    tips = geometry.tips

    def solve(first, second):
        """(time, error bound): Clairaut for curved cases, closed form for flat."""
        if mass:
            return independent_time(first, second, mass)
        return float(np.linalg.norm(second - first)), 0.0

    # Disagreements are judged against the independent solver's own error
    # bound: the Clairaut solver stops at a tolerance in impact parameter and
    # reports the time error that leaves (up to 4e-7 here, plumbing run).
    below = [solve(tips[0, 1:], x) for x in p]
    above = [solve(x, tips[1, 1:]) for x in q]
    own_low = tips[0, 0] + np.array([t for t, _ in below])
    own_high = tips[1, 0] - np.array([t for t, _ in above])
    bound = np.r_[[b for _, b in below], [b for _, b in above]]
    difference = np.r_[abs(own_low - low), abs(own_high - high)]
    # The ray ends: the independent solver must put them on the null boundary.
    omega, _ = geometry.directions(grid[0], grid[1])
    edge_gap, edge_excess = 0.0, 0.0
    for focus in geometry.foci:
        edge = geometry.ray_edges(focus, omega)
        for x in focus + edge[:, None] * omega:
            (top, b_top), (bottom, b_bottom) = (
                solve(x, tips[1, 1:]),
                solve(tips[0, 1:], x),
            )
            height = (tips[1, 0] - top) - (tips[0, 0] + bottom)
            edge_gap = max(edge_gap, abs(height))
            edge_excess = max(edge_excess, abs(height) - b_top - b_bottom)

    def density(x):
        rho = np.linalg.norm(x, axis=1)
        # sqrt(-g) = r^2 sin(theta) in areal coordinates, pulled back to isotropic.
        return (areal(rho, mass) / rho) ** 2 * (1 - (mass / (2 * rho)) ** 2)

    def integral(spatial, antiderivative):
        return float(
            np.sum(wq * density(q) * spatial(q) * antiderivative(own_high))
            - np.sum(wp * density(p) * spatial(p) * antiderivative(own_low))
        )

    c = geometry.base
    tests = (
        (lambda x: x[:, 0] - c[0], lambda t: t, events[:, 1] - c[0]),
        (lambda x: x[:, 1] - c[1], lambda t: t, events[:, 2] - c[1]),
        (
            lambda x: np.sum((x - c) ** 2, axis=1),
            lambda t: t,
            np.sum((events[:, 1:] - c) ** 2, axis=1),
        ),
        (
            lambda x: np.linalg.norm(x, axis=1),
            lambda t: t,
            np.linalg.norm(events[:, 1:], axis=1),
        ),
        (lambda x: np.ones(len(x)), lambda t: t**3 / 3, events[:, 0] ** 2),
        (
            lambda x: x[:, 1] - c[1],
            lambda t: t**2 / 2,
            events[:, 0] * (events[:, 2] - c[1]),
        ),
    )
    volume = integral(lambda x: np.ones(len(x)), lambda t: t)
    return {
        "events": len(events),
        "largest_bound_difference_from_generator": float(difference.max()),
        "largest_difference_beyond_independent_error_bound": float(
            (difference - bound).max()
        ),
        "largest_boundary_height_by_independent_solver": edge_gap,
        "largest_boundary_height_beyond_independent_error_bound": edge_excess,
        "proper_volume_over_flat": volume / (math.pi * geometry.proper_time**4 / 24),
        "moment_z_scores": [
            float(
                (sample.mean() - integral(f, g) / volume)
                / (sample.std() / math.sqrt(len(sample)))
            )
            for f, g, sample in tests
        ],
    }


def audit_case(folder, case, normal, pairs=1000):
    """Nested prefixes, a CPU replay, independent order bits and the measure."""
    from schwarzschild_orbit import MovingSchwarzschildDiamond

    kind = KINDS[case]
    mass = 0.0 if kind == "flat" else MASS
    speed = 0.0 if kind == "rest" else SPEED
    with np.load(folder / f"case_{case}_n{LEVELS[-1]}.npz", allow_pickle=False) as file:
        top = {key: file[key] for key in file.files}
    nested = int(top["measurement_n"]) == LEVELS[-1]
    for n in LEVELS[:-1]:
        with np.load(folder / f"case_{case}_n{n}.npz", allow_pickle=False) as file:
            nested &= set(file.files) == set(top) and int(file["measurement_n"]) == n
            for key in ("after", "before"):
                nested &= np.array_equal(file[key], top[key][:, : n // 8])
            for key in ("selection_order", "selection_design", "selected_ids"):
                nested &= np.array_equal(file[key], top[key])
    with np.load(folder / f"truth_{case}.npz", allow_pickle=False) as file:
        points = file["points"]
        labelled = str(file["kind"].item()) == kind and float(
            file["proper_time"]
        ) == float(PROPER_TIME)
    geometry = MovingSchwarzschildDiamond(mass, APEX, PROPER_TIME, speed)
    stream = np.random.default_rng([SEED + case, 1])
    sample = np.vstack([geometry.sample(4096, stream) for _ in range(256)])
    events = sample[:4096]
    stored = [
        np.unpackbits(top[key][:, :512], axis=1).astype(bool)
        for key in ("after", "before")
    ]
    replayed = geometry.relations(points[:, None], events[None])
    dt = events[None, :, 0] - points[:, None, 0]
    margin = abs(dt) - geometry.travel_time(points[:, None, 1:], events[None, :, 1:])
    report = {
        "case": case,
        "kind": kind,
        "landmarks_including_tips": len(points),
        "truth_file_labelled": bool(labelled),
        "tips_at_half_unit_time_in_normal_coordinates": float(
            abs(normal[[0, -1]] - [[-0.5, 0, 0, 0], [0.5, 0, 0, 0]]).max()
        ),
        "nested_prefixes_bit_exact": bool(nested),
        "tips_see_every_event": bool(
            (top["after"][0] == 255).all() and (top["before"][-1] == 255).all()
        ),
        "cpu_replay_of_first_4096_events_bit_exact": all(
            np.array_equal(a, b) for a, b in zip(replayed, stored, strict=True)
        ),
        "measure": measure_check(geometry, mass, sample),
    }
    if kind == "flat":
        flat = abs(dt) - np.linalg.norm(
            events[None, :, 1:] - points[:, None, 1:], axis=2
        )
        clear = abs(flat) > 1e-12
        related = (flat > 0) & clear
        report["closed_form_pairs_checked"] = int(clear.sum())
        report["closed_form_bits_match"] = bool(
            np.array_equal(stored[0][clear], (related & (dt > 0))[clear])
            and np.array_equal(stored[1][clear], (related & (dt < 0))[clear])
        )
        return report
    # The pairs nearest the light cone are the ones a solver error could flip.
    chosen = np.unique(
        np.r_[
            np.argsort(abs(margin), axis=None)[:pairs],
            np.random.default_rng([SEED + case, 2]).choice(margin.size, pairs, False),
        ]
    )
    matched, undecided, largest = 0, 0, 0.0
    for landmark, event in zip(*np.unravel_index(chosen, margin.shape), strict=True):
        duration, bound = independent_time(
            points[landmark, 1:], events[event, 1:], mass
        )
        gap = abs(dt[landmark, event]) - duration
        largest = max(largest, abs(gap - margin[landmark, event]))
        if abs(gap) <= bound:
            undecided += 1
            continue
        ahead = dt[landmark, event] > 0
        bits = (stored[0][landmark, event], stored[1][landmark, event])
        matched += bits == (gap > 0 and ahead, gap > 0 and not ahead)
    report.update(
        independent_pairs_checked=len(chosen),
        independent_pairs_inside_the_error_bound=undecided,
        independent_bits_match=int(matched),
        largest_margin_difference_between_solvers=largest,
        smallest_margin_checked=float(abs(margin).min()),
    )
    return report


# -- refits ----------------------------------------------------------------------------
def gauge_aligned(old, new):
    """Compare two fits after aligning the free signs of the three spatial axes."""
    a = np.asarray(old["coordinates_inferred_from_counts"])
    b = np.asarray(new["coordinates_inferred_from_counts"])
    signs = np.ones(4)
    for axis, row in ((1, 1), (2, 2), (3, None)):
        if row is not None:
            signs[axis] = np.sign(a[row, axis]) * np.sign(b[row, axis])
        else:
            # z: the sign that best maps the new cloud onto the old.
            signs[axis] = 1.0 if np.sum(a[:, 3] * b[:, 3]) >= 0 else -1.0
    turn = np.asarray(_reflection_matrix(signs[1:]))
    c_old = np.asarray(old["coefficients"])
    c_new = turn @ np.asarray(new["coefficients"])
    s_old = np.asarray(old["nominal_component_covariance"])
    s_new = turn @ np.asarray(new["nominal_component_covariance"]) @ turn.T
    scale_c = max(1e-12, float(np.max(abs(c_old))))
    scale_s = max(1e-30, float(np.max(abs(s_old))))
    return {
        "axis_signs_applied": signs[1:].tolist(),
        "coordinates_max_difference": float(np.max(abs(a - b * signs))),
        "coefficients_relative_difference": float(np.max(abs(c_old - c_new)) / scale_c),
        "covariance_relative_difference": float(np.max(abs(s_old - s_new)) / scale_s),
        # Lists that do not depend on the axis signs.
        "invariant_lists_differing": sorted(
            key
            for key in ("electric_eigenvalues", "flat_costs", "selected_ids")
            if not np.allclose(old[key], new[key], rtol=1e-8, atol=1e-12)
        ),
        "scalar_fields_differing": sorted(
            key
            for key in old
            if not isinstance(old[key], list)
            and key in new
            and old[key] != new[key]
            and not (isinstance(old[key], float) and close(old[key], new[key], 1e-8))
        ),
    }


def _reflection_matrix(signs):
    from check_intrinsic_weyl_bombelli_benchmark import reflection

    return reflection(list(signs))


def refit(folder, case, n):
    """Re-estimate one row in a fresh single-threaded process."""
    target = folder.with_name(folder.name + "_audit")
    target.mkdir(exist_ok=True)
    output = target / f"case_{case}_n{n}_{int(time.time())}.json"
    env = {k: v for k, v in os.environ.items() if k not in THREAD_VARIABLES}
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    subprocess.run(
        [
            sys.executable,
            str(Path(__file__).with_name("intrinsic_weyl_schwarzschild_orbit.py")),
            "--worker",
            str(folder / f"case_{case}_n{n}.npz"),
            "--truth",
            str(folder / f"truth_{case}.npz"),
            "--output",
            str(output),
        ],
        check=True,
        env=env,
    )
    new = read(output.with_name(output.stem + "_pure.json"))["fit"]
    old = read(folder / f"case_{case}_n{n}_pure.json")["fit"]
    aligned = gauge_aligned(old, new)
    return {
        "case": case,
        "measurement_n": n,
        "q_old": statistic(old),
        "q_new": statistic(new),
        "raw_fields_identical": new == old,
        "raw_differences_beyond_1e-8": len(differences(new, old)),
        "gauge_aligned": aligned,
        "agrees_after_gauge_alignment": bool(
            close(statistic(old), statistic(new))
            and aligned["coordinates_max_difference"] < 1e-8
            and aligned["coefficients_relative_difference"] < 1e-8
            and aligned["covariance_relative_difference"] < 1e-8
            and not aligned["invariant_lists_differing"]
            and not aligned["scalar_fields_differing"]
        ),
    }


def case_passed(c):
    """A skipped case passes only if the run itself recorded its failure."""
    if "skipped" in c:
        return bool(c["run_recorded_a_failure"])
    return bool(
        c["truth_file_labelled"]
        and c["tips_at_half_unit_time_in_normal_coordinates"] < 1e-9
        and c["nested_prefixes_bit_exact"]
        and c["tips_see_every_event"]
        and c["cpu_replay_of_first_4096_events_bit_exact"]
        and c.get("closed_form_bits_match", True)
        and c.get("independent_bits_match", 0)
        == c.get("independent_pairs_checked", 0)
        - c.get("independent_pairs_inside_the_error_bound", 0)
        and c["measure"]["largest_difference_beyond_independent_error_bound"] < 1e-9
        and c["measure"]["largest_boundary_height_beyond_independent_error_bound"]
        < 1e-9
        and max(map(abs, c["measure"]["moment_z_scores"])) < 4.5
    )


def audit(folder, receipt, refits):
    from schwarzschild_orbit import MovingSchwarzschildDiamond
    from scipy.optimize import brentq
    from scipy.stats import chi2, ncx2
    from scipy.stats import norm as gaussian

    q10, q5 = quantile10(1 - FALSE_POSITIVE), quantile5(1 - FALSE_POSITIVE)
    z95 = normal_quantile(0.975)
    # The independent quantiles agree with the library ones the runner used.
    assert abs(q10 - chi2.ppf(1 - FALSE_POSITIVE, 10)) < 1e-9
    assert abs(q5 - chi2.ppf(1 - FALSE_POSITIVE, 5)) < 1e-9
    assert abs(z95 - gaussian.ppf(0.975)) < 1e-9
    summary = read(folder / "summary.json")
    stamp = summary["code_version"]

    def git(*arguments):
        return subprocess.run(
            ["git", *arguments], cwd=ROOT, capture_output=True, text=True
        )

    truth = {}
    tensors_by_kind = {}
    for kind, moving, mass in (
        ("orbit", True, MASS),
        ("rest", False, MASS),
        ("flat", True, 0.0),
    ):
        if mass:
            electric, magnetic, frame_error = areal_truth(
                mass, APEX, PROPER_TIME, moving
            )
        else:
            electric, magnetic, frame_error = np.zeros((3, 3)), np.zeros((3, 3)), 0.0
        recorded = summary["geometry"][kind]
        tensors_by_kind[kind] = (electric, magnetic)
        truth[kind] = {
            "areal_electric": electric.tolist(),
            "areal_magnetic": magnetic.tolist(),
            "orthonormality_error_of_areal_frame": frame_error,
            "difference_from_generator_closed_form": float(
                max(
                    abs(electric - np.asarray(recorded["true_electric"])).max(),
                    abs(magnetic - np.asarray(recorded["true_magnetic"])).max(),
                )
            ),
            "invariant": float(np.sum(electric**2) - np.sum(magnetic**2)),
            "electric_dot_magnetic": float(np.sum(electric * magnetic)),
        }
    normal, missing = {}, {}
    for case, kind in enumerate(KINDS):
        mass = 0.0 if kind == "flat" else MASS
        speed = 0.0 if kind == "rest" else SPEED
        path = folder / f"truth_{case}.npz"
        normal[case] = None
        if not path.exists():
            missing[case] = "no truth file: the selection did not complete"
            continue
        try:
            with np.load(path, allow_pickle=False) as file:
                normal[case] = MovingSchwarzschildDiamond(
                    mass, APEX, PROPER_TIME, speed
                ).normal_coordinates(file["points"])
        except (ValueError, ArithmeticError) as error:
            missing[case] = f"{type(error).__name__}: {error}"
    recorded = {(r["case"], r["measurement_n"]): r for r in summary["results"]}
    rows = []
    for case in range(len(KINDS)):
        for n in LEVELS:
            entry = recorded[case, n]
            if "failure" in entry or "detection" not in entry:
                rows.append(
                    {
                        "case": case,
                        "kind": KINDS[case],
                        "measurement_n": n,
                        "valid_fit": False,
                        "failure": entry.get("failure"),
                    }
                )
                continue
            rows.append(
                audit_row(
                    folder,
                    case,
                    n,
                    entry,
                    normal[case],
                    *tensors_by_kind[KINDS[case]],
                    z95,
                )
            )
    level = {}
    for n in LEVELS:
        level[str(n)] = verdicts([r for r in rows if r["measurement_n"] == n])
    recorded_verdicts = {
        str(n): {
            gate: summary["level_verdicts"][gate][str(n)]
            for gate in summary["level_verdicts"]
        }
        for n in LEVELS
    }
    cases = []
    for case in range(len(KINDS)):
        top = folder / f"case_{case}_n{LEVELS[-1]}.npz"
        if normal[case] is None or not top.exists():
            # Allowed only for a failure the run itself recorded (checked below).
            cases.append(
                {
                    "case": case,
                    "kind": KINDS[case],
                    "skipped": missing.get(case, "top-level archive missing"),
                    "run_recorded_a_failure": any(
                        "failure" in r for r in summary["results"] if r["case"] == case
                    ),
                }
            )
            continue
        cases.append(audit_case(folder, case, normal[case]))
    power = {}
    for case in (0, 1):
        top = next(
            r for r in rows if r["case"] == case and r["measurement_n"] == LEVELS[-1]
        )
        if not top.get("valid_fit"):
            continue
        rates = {
            "field": (
                top["posthoc"]["noncentrality_of_the_truth"] / LEVELS[-1],
                10,
                q10,
            ),
            "magnetic": (
                top["posthoc"]["magnetic_noncentrality_of_the_truth"] / LEVELS[-1],
                5,
                q5,
            ),
        }
        power[case] = {
            name: {
                "noncentrality_per_event": rate,
                **{
                    f"events_for_power_{int(100 * p)}": brentq(
                        lambda n, p=p, rate=rate, df=df, limit=limit: (
                            ncx2.sf(limit, df, rate * n) - p
                        ),
                        1e3,
                        1e11,
                    )
                    for p in (0.5, 0.9)
                },
            }
            for name, (rate, df, limit) in rates.items()
        }
    result = {
        "audited_folder": str(folder),
        "audit_code_version": git("rev-parse", "HEAD").stdout.strip(),
        "audit_worktree_clean": not git("status", "--porcelain").stdout.strip(),
        "run_code_version": stamp,
        "run_stamp_is_ancestor_of_head": git(
            "merge-base", "--is-ancestor", stamp, "HEAD"
        ).returncode
        == 0,
        "summary_sha256": sha(folder / "summary.json"),
        "manifest_sha256": sha(folder / "manifest.json"),
        "thresholds": {"field_chi2_10": q10, "magnetic_chi2_5": q5, "normal_975": z95},
        "truth": truth,
        "rows": rows,
        "level_verdicts": level,
        "level_verdicts_match_run": level == recorded_verdicts,
        "primary_claim": primary(level),
        "primary_claim_matches_run": primary(level) == summary["primary_claim"],
        # Descriptive label of the protocol, restated from our own containment.
        "consistency_label": {
            str(n): (
                "UNAVAILABLE"
                if any(
                    r.get("error_mod_reflection") is None
                    for r in rows
                    if r["measurement_n"] == n and r["kind"] == "orbit"
                )
                else "CONSISTENT"
                if all(
                    r["error_mod_reflection"] <= r["radius_95"]
                    for r in rows
                    if r["measurement_n"] == n and r["kind"] == "orbit"
                )
                else "INCONSISTENT"
            )
            for n in LEVELS
        },
        "cases": cases,
        # Unbiased-estimate assumption: the model bias is not measured in this run.
        "posthoc_power_if_unbiased": power,
        "refits": [
            refit(folder, case, n)
            if (folder / f"case_{case}_n{n}.npz").exists()
            else {"case": case, "measurement_n": n, "not_applicable": "no archive"}
            for case, n in dict.fromkeys((*FROZEN_REFITS, *refits))
        ],
    }
    result["passed"] = bool(
        result["run_stamp_is_ancestor_of_head"]
        and result["level_verdicts_match_run"]
        and result["primary_claim_matches_run"]
        and result["consistency_label"] == summary["consistency_label"]
        and all(
            t["difference_from_generator_closed_form"] < 1e-7 for t in truth.values()
        )
        and all(
            t["orthonormality_error_of_areal_frame"] < 1e-12 for t in truth.values()
        )
        and all(
            r.get("failure") is not None
            or (
                r["valid_fit"]
                and r["recorded_values_reproduced"]
                and r["blind_boundary_kept"]
                and r["archive_hash_matches"]
            )
            for r in rows
        )
        and all(case_passed(c) for c in cases)
        and all(
            r.get("agrees_after_gauge_alignment", "not_applicable" in r)
            for r in result["refits"]
        )
    )
    write_json(receipt, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", type=Path)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--refit", action="append", default=[], metavar="CASE:N")
    args = parser.parse_args()
    outcome = audit(
        args.folder.resolve(),
        args.receipt or args.folder.with_name(args.folder.name + ".audit.json"),
        [tuple(int(part) for part in item.split(":")) for item in args.refit],
    )
    for line in outcome["rows"]:
        print(json.dumps(line))
    for line in outcome["cases"]:
        print(json.dumps(line))
    print(
        json.dumps(
            {
                key: outcome[key]
                for key in (
                    "truth",
                    "level_verdicts",
                    "level_verdicts_match_run",
                    "primary_claim",
                    "primary_claim_matches_run",
                    "posthoc_power_if_unbiased",
                    "refits",
                    "passed",
                )
            }
        )
    )
    sys.exit(0 if outcome["passed"] else 1)
