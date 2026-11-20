"""Pre-registered test of the unchanged order-only Weyl inverse on an interval
whose tips move on a circular Schwarzschild geodesic.

Two observers pass through one event, the apex at areal radius 6M. The one on
the circular geodesic (speed 1/2 relative to the static observer) sees
T^2 E = diag(-1/8, 1/24, 1/12) and a magnetic part T^2 B_xz = 1/12, so that
|E| != |B| != 0; the one at rest sees T^2 E = diag(-1/12, 1/24, 1/24) and B = 0.
Orbit intervals are the test, rest intervals the magnetic-null control, and
flat intervals with the same moving tips the flat control. Only the generator
is new; landmark selection, clocks, covariance, the fifth-degree inverse with
curvature-square terms and the confidence construction are the frozen ones.
See docs/paper/paper_a/intrinsic_weyl/schwarzschild_orbit_protocol.md.
"""

import argparse
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

from intrinsic_weyl_bombelli_benchmark import lf_sha, shutil_free
from replay_intrinsic_weyl_threads import (
    ROOT,
    THREAD_VARIABLES,
    differences,
    sha,
    write_json,
)

PROTOCOL = "docs/paper/paper_a/intrinsic_weyl/schwarzschild_orbit_protocol.md"
MASS, APEX, PROPER_TIME, SPEED = 1.0, 6.0, 3.0, 0.5
SELECTION_N, CHUNK, BATCH = 4096, 4096, 65536
LEVELS = (1048576, 4194304, 16777216, 67108864)
SEED = 2026112200
FALSE_POSITIVE, COVERAGE = 1e-3, 0.95
GATES = ("field", "magnetic", "unequal_norms")
# kind, seed offset. All six are assigned at every level.
CASES = (
    ("orbit", 0),
    ("orbit", 1),
    ("rest", 2),
    ("rest", 3),
    ("flat", 4),
    ("flat", 5),
)
# Engineering population: a stronger field on a larger interval, disjoint seed.
ENGINEERING_LEVELS = (1048576, 4194304, 16777216)
ENGINEERING_CASES = (("orbit", 90),)
ENGINEERING_PROPER_TIME = 4.0
ORACLES = (
    "schwarzschild_orbit",
    "schwarzschild_diamond",
    "profile_plane_wave",
    "rotating_plane_wave",
    "rotating_plane_wave_gpu",
    "window_response",
    "p14_plane_wave",
    "intrinsic_weyl_rotating_pilot",
    "intrinsic_weyl_joint_audit",
    "intrinsic_weyl_joint_pilot",
)


def build(kind, proper_time):
    from schwarzschild_orbit import MovingSchwarzschildDiamond, circular_velocity

    if abs(circular_velocity(MASS, APEX) - SPEED) > 1e-15:
        raise ValueError("The assigned speed is not the circular speed")
    if kind == "orbit":
        return MovingSchwarzschildDiamond(MASS, APEX, proper_time, SPEED)
    if kind == "rest":
        return MovingSchwarzschildDiamond(MASS, APEX, proper_time, 0.0)
    if kind == "flat":
        return MovingSchwarzschildDiamond(0.0, APEX, proper_time, SPEED)
    raise ValueError("Unknown test geometry")


def thresholds():
    from scipy.stats import chi2, norm

    return {
        "field_chi2_10": float(chi2.ppf(1 - FALSE_POSITIVE, 10)),
        "magnetic_chi2_5": float(chi2.ppf(1 - FALSE_POSITIVE, 5)),
        "interval_normal_quantile": float(norm.ppf(1 - (1 - COVERAGE) / 2)),
    }


# -- frozen gates ----------------------------------------------------------------------
def _assigned(rows):
    if len(rows) != len(CASES) or [r["kind"] for r in rows] != [k for k, _ in CASES]:
        raise ValueError("All six assigned cases are required, in order")


def _failed(row):
    return "failure" in row or "detection" not in row


def field_verdict(rows):
    """G1: the field departs from flat in both orbit cases and in no flat control.

    Bad news first: a valid flat exceedance is CONTROL_FAILED even when another
    used row failed; only then does a failure make the level INCOMPLETE.
    """
    _assigned(rows)
    used = [r for r in rows if r["kind"] in ("orbit", "flat")]
    limit = thresholds()["field_chi2_10"]
    valid = [r for r in used if not _failed(r)]
    if any(r["detection"]["q"] > limit for r in valid if r["kind"] == "flat"):
        return "CONTROL_FAILED"
    if len(valid) < len(used):
        return "INCOMPLETE"
    orbit = [r["detection"]["q"] > limit for r in used if r["kind"] == "orbit"]
    if all(orbit):
        return "DETECTED"
    return "MIXED" if any(orbit) else "NOT_DETECTED"


def magnetic_verdict(rows):
    """G2: a magnetic part in both orbit cases and in none of the four controls."""
    _assigned(rows)
    limit = thresholds()["magnetic_chi2_5"]
    valid = [r for r in rows if not _failed(r)]
    if any(r["detection"]["q_magnetic"] > limit for r in valid if r["kind"] != "orbit"):
        return "CONTROL_FAILED"
    if len(valid) < len(rows):
        return "INCOMPLETE"
    orbit = [r["detection"]["q_magnetic"] > limit for r in rows if r["kind"] == "orbit"]
    if all(orbit):
        return "MAGNETIC_DETECTED"
    return "MIXED" if any(orbit) else "NOT_DETECTED"


def unequal_verdict(rows):
    """G3: |E|^2 - |B|^2 above zero on its 95% interval in both orbit cases."""
    _assigned(rows)
    used = [r for r in rows if r["kind"] == "orbit"]
    valid = [r for r in used if not _failed(r)]
    if any(r["detection"]["invariant_interval_95"][1] < 0 for r in valid):
        return "REVERSED"
    if len(valid) < len(used):
        return "INCOMPLETE"
    above = [r["detection"]["invariant_interval_95"][0] > 0 for r in used]
    if all(above):
        return "UNEQUAL"
    return "MIXED" if any(above) else "NOT_RESOLVED"


def consistency(rows):
    """Descriptive label, fixed in advance: the truth inside both orbit balls."""
    _assigned(rows)
    used = [r for r in rows if r["kind"] == "orbit"]
    if any(_failed(r) or "evaluation" not in r for r in used):
        return "UNAVAILABLE"
    inside = all(r["evaluation"]["nominal_ball_contains_truth"] for r in used)
    return "CONSISTENT" if inside else "INCONSISTENT"


def primary_claim(level_verdicts):
    """The headline, fixed in advance: all three gates pass at the top level."""
    top = str(LEVELS[-1])
    outcome = {gate: level_verdicts[gate][top] for gate in GATES}
    if any(v == "INCOMPLETE" for v in outcome.values()):
        return "INCOMPLETE"
    passed = (
        outcome["field"] == "DETECTED"
        and outcome["magnetic"] == "MAGNETIC_DETECTED"
        and outcome["unequal_norms"] == "UNEQUAL"
    )
    return "SUPPORTED" if passed else "NOT_SUPPORTED"


def selfcheck():
    """The gates on synthetic rows: every branch, and the failure rules."""
    limits = thresholds()
    q10, q5 = limits["field_chi2_10"], limits["magnetic_chi2_5"]
    kinds = [k for k, _ in CASES]

    def rows(q=None, qb=None, bounds=None, failed=()):
        out = []
        for i, kind in enumerate(kinds):
            if i in failed:
                out.append({"kind": kind, "failure": "kept"})
                continue
            out.append(
                {
                    "kind": kind,
                    "detection": {
                        "q": (q or [1.0] * 6)[i],
                        "q_magnetic": (qb or [1.0] * 6)[i],
                        "invariant_interval_95": (bounds or [(-1.0, 1.0)] * 6)[i],
                    },
                }
            )
        return out

    hi10, lo10, hi5, lo5 = q10 + 1, q10 - 1, q5 + 1, q5 - 1
    assert field_verdict(rows([hi10, hi10, hi10, hi10, lo10, lo10])) == "DETECTED"
    assert field_verdict(rows([hi10, lo10, lo10, lo10, lo10, lo10])) == "MIXED"
    assert field_verdict(rows([q10, hi10, 0, 0, 0, 0])) == "MIXED"
    assert field_verdict(rows([lo10] * 6)) == "NOT_DETECTED"
    assert field_verdict(rows([hi10, hi10, 0, 0, hi10, 0])) == "CONTROL_FAILED"
    # The rest cases are not part of G1.
    assert field_verdict(rows([hi10, hi10, 0, 0, 0, 0], failed=(2,))) == "DETECTED"
    assert field_verdict(rows([hi10] * 4 + [0, 0], failed=(4,))) == "INCOMPLETE"
    # Bad news is not hidden by a failure elsewhere.
    assert field_verdict(rows([hi10] * 6, failed=(4,))) == "CONTROL_FAILED"
    assert field_verdict(rows([0, 0, 0, 0, 0, hi10], failed=(0,))) == "CONTROL_FAILED"
    assert magnetic_verdict(rows(qb=[hi5, hi5, lo5, lo5, lo5, lo5])) == (
        "MAGNETIC_DETECTED"
    )
    assert magnetic_verdict(rows(qb=[hi5, lo5, lo5, lo5, lo5, lo5])) == "MIXED"
    assert magnetic_verdict(rows(qb=[q5, q5, 0, 0, 0, 0])) == "NOT_DETECTED"
    assert magnetic_verdict(rows(qb=[hi5, hi5, hi5, 0, 0, 0])) == "CONTROL_FAILED"
    assert magnetic_verdict(rows(qb=[hi5, hi5, 0, 0, 0, hi5])) == "CONTROL_FAILED"
    assert magnetic_verdict(rows(qb=[hi5] * 2 + [0] * 4, failed=(3,))) == ("INCOMPLETE")
    assert magnetic_verdict(rows(qb=[0, 0, hi5, 0, 0, 0], failed=(0,))) == (
        "CONTROL_FAILED"
    )
    up, down, mid = (0.1, 0.2), (-0.2, -0.1), (-0.1, 0.1)
    assert unequal_verdict(rows(bounds=[up, up] + [mid] * 4)) == "UNEQUAL"
    assert unequal_verdict(rows(bounds=[up, mid] + [mid] * 4)) == "MIXED"
    assert unequal_verdict(rows(bounds=[mid, mid] + [up] * 4)) == "NOT_RESOLVED"
    assert unequal_verdict(rows(bounds=[up, down] + [mid] * 4)) == "REVERSED"
    assert unequal_verdict(rows(bounds=[(0.0, 0.2), up] + [mid] * 4)) == "MIXED"
    assert unequal_verdict(rows(bounds=[up, up] + [mid] * 4, failed=(5,))) == (
        "UNEQUAL"
    )
    assert unequal_verdict(rows(bounds=[up, up] + [mid] * 4, failed=(1,))) == (
        "INCOMPLETE"
    )
    assert unequal_verdict(rows(bounds=[down, up] + [mid] * 4, failed=(1,))) == (
        "REVERSED"
    )
    every = {g: {str(n): "x" for n in LEVELS} for g in GATES}
    top = str(LEVELS[-1])
    every["field"][top], every["magnetic"][top] = "DETECTED", "MAGNETIC_DETECTED"
    every["unequal_norms"][top] = "UNEQUAL"
    assert primary_claim(every) == "SUPPORTED"
    every["unequal_norms"][top] = "MIXED"
    assert primary_claim(every) == "NOT_SUPPORTED"
    every["magnetic"][top] = "INCOMPLETE"
    assert primary_claim(every) == "INCOMPLETE"
    try:
        field_verdict(rows()[:5])
    except ValueError:
        pass
    else:
        raise AssertionError("A missing assignment was accepted")
    return {
        **limits,
        "strict_exceedance_required": True,
        "control_exceedance_invalidates": True,
        "bad_news_precedes_incomplete": True,
        "primary_level": LEVELS[-1],
        "failure_cannot_disappear": True,
        "all_six_assigned": True,
    }


# -- statistics of one fit -------------------------------------------------------------
def statistics(inference):
    """Q, Q_B, the invariant |E|^2 - |B|^2 with its interval, from one inference."""
    import numpy as np
    from intrinsic_weyl import WEYL

    c = np.asarray(inference["coefficients_xy_positive"])
    s = np.asarray(inference["nominal_covariance_xy_positive"])
    basis = WEYL[:5, 0, 1:, 0, 1:]
    gram = np.einsum("kij,lij->kl", basis, basis)
    zero = np.zeros((5, 5))
    form = np.block([[gram, zero], [zero, -gram]])
    cross = np.block([[zero, gram / 2], [gram / 2, zero]])
    limits = thresholds()

    def quadratic(matrix):
        """Debiased c^T A c and the variance of c^T A c for Gaussian c, plug-in."""
        value = float(c @ matrix @ c - np.trace(matrix @ s))
        variance = float(
            4 * c @ matrix @ s @ matrix @ c + 2 * np.trace(matrix @ s @ matrix @ s)
        )
        return value, float(np.sqrt(variance))

    invariant, invariant_se = quadratic(form)
    product, product_se = quadratic(cross)
    z = limits["interval_normal_quantile"]
    electric = np.einsum("k,kij->ij", c[:5], basis)
    magnetic = np.einsum("k,kij->ij", c[5:], basis)
    return {
        "q": float(c @ np.linalg.solve(s, c)),
        "q_magnetic": float(c[5:] @ np.linalg.solve(s[5:, 5:], c[5:])),
        "q_electric": float(c[:5] @ np.linalg.solve(s[:5, :5], c[:5])),
        "invariant": invariant,
        "invariant_se": invariant_se,
        "invariant_interval_95": [
            invariant - z * invariant_se,
            invariant + z * invariant_se,
        ],
        "invariant_raw": float(c @ form @ c),
        "electric_dot_magnetic": product,
        "electric_dot_magnetic_se": product_se,
        "electric_eigenvalues": np.linalg.eigvalsh(electric).tolist(),
        "magnetic_eigenvalues": np.linalg.eigvalsh(magnetic).tolist(),
        "electric_norm": float(np.linalg.norm(electric)),
        "magnetic_norm": float(np.linalg.norm(magnetic)),
    }


def truth_coefficients(geometry, points):
    """True E and B in the internal gauge of the landmarks' normal coordinates."""
    import numpy as np
    from intrinsic_weyl import WEYL

    # The inverse fixes its frame by the first two landmarks after the past tip.
    normal = geometry.normal_coordinates(points[1:3])
    first = normal[0, 1:] / np.linalg.norm(normal[0, 1:])
    second = normal[1, 1:] - (normal[1, 1:] @ first) * first
    second /= np.linalg.norm(second)
    frame = np.column_stack((first, second, np.cross(first, second)))
    basis = WEYL[:5, 0, 1:, 0, 1:]
    gram = np.einsum("kij,lij->kl", basis, basis)

    def project(tensor):
        rotated = frame.T @ tensor @ frame
        return np.linalg.solve(gram, np.einsum("kij,ij->k", basis, rotated))

    # A proper rotation: the magnetic part transforms as a tensor here.
    return np.r_[project(geometry.tidal()), project(geometry.magnetic())]


def readings(estimate, covariance, target):
    """S^-1 weighted amplitude of the estimate along a target, with its sigma."""
    import numpy as np

    weight = np.linalg.solve(covariance, target)
    information = float(target @ weight)
    return float(weight @ estimate) / information, 1 / float(np.sqrt(information))


def detect(fit, n):
    """Gate statistics from the pure fit alone; no truth and no generator."""
    from intrinsic_weyl_joint import confidence_ball
    from scipy.stats import chi2

    inference = confidence_ball(fit, coverage=COVERAGE)
    detection = statistics(inference)
    detection.update(
        thresholds(),
        nominal_p_value=float(chi2.sf(detection["q"], 10)),
        nominal_p_value_magnetic=float(chi2.sf(detection["q_magnetic"], 5)),
        tidal_norm_estimate=inference["tidal_norm"],
        largest_tidal_variance=inference["largest_tidal_variance"],
        variance_times_N=n * inference["largest_tidal_variance"],
    )
    return inference, detection


def evaluate_fit(inference, truth_path, n):
    """Oracle evaluation, after the pure fit and its gate statistics are saved."""
    import numpy as np
    from intrinsic_weyl_joint import TIDAL_ROOT, axis_signs, tidal_distance

    estimate = np.asarray(inference["coefficients_xy_positive"])
    covariance = np.asarray(inference["nominal_covariance_xy_positive"])
    with np.load(truth_path, allow_pickle=False) as saved:
        points, kind = saved["points"], str(saved["kind"].item())
        proper_time = float(saved["proper_time"])
    geometry = build(kind, proper_time)
    expected = truth_coefficients(geometry, points)
    electric, magnetic = geometry.tidal(), geometry.magnetic()
    norm = float(np.sqrt(np.sum(electric**2) + np.sum(magnetic**2)))
    error = tidal_distance(estimate, expected)
    evaluation = {
        "truth_coefficients_xy_positive": expected.tolist(),
        "truth_norm": norm,
        "truth_electric_norm": float(np.linalg.norm(electric)),
        "truth_magnetic_norm": float(np.linalg.norm(magnetic)),
        "truth_electric_eigenvalues": np.linalg.eigvalsh(electric).tolist(),
        "truth_magnetic_eigenvalues": np.linalg.eigvalsh(magnetic).tolist(),
        "truth_invariant": float(np.sum(electric**2) - np.sum(magnetic**2)),
        "full_tidal_error_mod_reflection": error,
        "relative_full_tidal_error": error / norm if norm else None,
        "confidence_radius_95": inference["radius"],
        "nominal_ball_contains_truth": bool(error <= inference["radius"]),
        "continuum_bias": None,
        "continuum_bias_note": "Exact interval volumes between landmarks are not "
        "implemented for this geometry",
    }
    candidates = (expected, axis_signs([1, 1, -1]) * expected)
    # The reflection candidate nearer to the estimate, in the tidal metric.
    matched = min(
        candidates, key=lambda item: np.linalg.norm(TIDAL_ROOT @ (estimate - item))
    )
    evaluation["q_about_truth"] = min(
        float((estimate - t) @ np.linalg.solve(covariance, estimate - t))
        for t in candidates
    )
    if norm:
        rows = {}
        for name, part in (
            ("full", slice(0, 10)),
            ("electric", slice(0, 5)),
            ("magnetic", slice(5, 10)),
        ):
            target = matched[part]
            if not np.any(target):
                continue
            block = covariance[part, part]
            value, sigma = readings(estimate[part], block, target)
            snr = min(
                float(np.sqrt(t[part] @ np.linalg.solve(block, t[part])))
                for t in candidates
            )
            rows[name] = {
                "reading_matched_reflection": value,
                "reading_sigma": sigma,
                "expected_snr_less_favorable": snr,
                "events_for_expected_snr_one": n / snr**2,
                # TIDAL_ROOT is block diagonal, so a slice is that part's error.
                "part_error_matched_reflection": float(
                    np.linalg.norm((TIDAL_ROOT @ (estimate - matched))[part])
                ),
            }
        evaluation["truth_pattern"] = rows
    return evaluation


# -- GPU checks before data ------------------------------------------------------------
def relation_selfcheck(proper_time):
    """GPU relation bits against the direct solver, for the three classes."""
    import cupy as cp
    import numpy as np

    report = {}
    for kind in ("orbit", "rest", "flat"):
        geometry = build(kind, proper_time)
        rng = np.random.default_rng(20261128)
        landmarks = np.vstack(
            (geometry.tips[0], geometry.sample(62, rng), geometry.tips[1])
        )
        events = geometry.sample(4096, rng)
        packed = geometry.packed_relations(landmarks, events, cp)
        direct = geometry.relations(landmarks[:, None], events[None])
        assert all(
            np.array_equal(np.unpackbits(p, axis=1).astype(bool), d)
            for p, d in zip(packed, direct, strict=True)
        )
        report[kind] = {
            "pairs": int(direct[0].size),
            "related": int(direct[0].sum() + direct[1].sum()),
            "resolved_by_the_direct_solver": geometry.resolved_margins,
        }
    cp.get_default_memory_pool().free_all_blocks()
    return report


# -- worker, generation and the run ----------------------------------------------------
def worker(archive, truth, output):
    """Pure order-only fit first; oracle modules load only after it is saved."""
    import intrinsic_weyl_joint as joint
    import numpy as np
    from intrinsic_weyl_gpu_counts import clock_second_moment
    from intrinsic_weyl_high_count import peak_working_set

    joint._clock_second_moment = clock_second_moment
    wall, cpu = time.perf_counter(), time.process_time()
    record = {
        "archive_sha256": sha(archive),
        "thread_environment": {key: os.environ.get(key) for key in THREAD_VARIABLES},
    }
    original_load, opened = np.load, []

    def guard(selected, *args, **options):
        if Path(selected).resolve() != archive.resolve():
            raise ValueError("Only the assigned order archive may be opened")
        opened.append(str(Path(selected).resolve()))
        return original_load(selected, *args, **options)

    try:
        assert not any(name in sys.modules for name in ORACLES)
        np.load = guard
        try:
            fit = joint.fit_order_archive(
                archive, model_order=5, square_correction=True
            )
        finally:
            np.load = original_load
        imported = [name for name in ORACLES if name in sys.modules]
        if imported:
            raise RuntimeError(f"Oracle imported before pure-fit boundary: {imported}")
        record["finite_fit_wall_seconds"] = time.perf_counter() - wall
        write_json(
            output.with_name(output.stem + "_pure.json"),
            {
                "fit": fit,
                "opened_archives": opened,
                "oracle_modules_before_return": imported,
                "inverse_module": joint.__file__,
                "archive_sha256": sha(archive),
            },
        )
        inference, detection = detect(fit, fit["measurement_n"])
        record["full"] = {"fit": fit, "inference": inference, "detection": detection}
        # The gate statistics are kept even if the truth step below fails.
        write_json(output.with_name(output.stem + "_detection.json"), record)
        try:
            record["full"]["evaluation"] = evaluate_fit(
                inference, truth, fit["measurement_n"]
            )
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
            record["evaluation_failure"] = f"{type(error).__name__}: {error}"
    except (ValueError, ArithmeticError, MemoryError, np.linalg.LinAlgError) as error:
        record["failure"] = f"{type(error).__name__}: {error}"
    record.update(
        wall_seconds=time.perf_counter() - wall,
        process_cpu_seconds=time.process_time() - cpu,
        peak_working_set_bytes=peak_working_set(),
    )
    write_json(output, record)


def generate(folder, index, kind, proper_time, seed, levels, status):
    """Selection order, landmarks, and nested measurement prefixes for one case."""
    import cupy as cp
    import numpy as np
    from intrinsic_weyl_high_count import peak_working_set
    from intrinsic_weyl_joint import selection

    geometry = build(kind, proper_time)
    wall, cpu = time.perf_counter(), time.process_time()
    generation = {
        "selection_events_drawn": 0,
        "sampled_events": 0,
        "classified_events": 0,
    }
    truth_path = folder / f"truth_{index}.npz"
    try:
        status(state="RUNNING", phase="selection", case=index)
        rng = np.random.default_rng([seed, 0])
        points = np.vstack(
            (geometry.tips[0], geometry.sample(SELECTION_N, rng), geometry.tips[1])
        )
        generation["selection_events_drawn"] = SELECTION_N
        # Labels are randomized independently of coordinates; only the two
        # supplied tips have distinguished IDs.
        permutation = np.r_[
            0, rng.permutation(np.arange(1, len(points) - 1)), len(points) - 1
        ]
        points = points[permutation]
        order = geometry.order(points)
        try:
            ids, _ = selection(order, design="nested")
        except ValueError:
            np.savez(
                folder / f"case_{index}_selection_failed.npz",
                selection_order=order,
                selection_design="nested",
            )
            raise
        selected = {
            "selection_order": order,
            "selection_design": np.array("nested"),
            "selected_ids": ids,
        }
        points = points[ids]
        np.savez(
            truth_path,
            points=points,
            kind=np.array(kind),
            proper_time=np.array(float(proper_time)),
        )
        generation["landmarks_including_tips"] = len(points)
        top = levels[-1]
        after = np.empty((len(points), top // 8), dtype=np.uint8)
        before = np.empty_like(after)
        rng = np.random.default_rng([seed, 1])
        for start in range(0, top, BATCH):
            end = min(start + BATCH, top)
            events = np.vstack(
                [geometry.sample(CHUNK, rng) for _ in range((end - start) // CHUNK)]
            )
            generation["sampled_events"] += len(events)
            future, past = geometry.packed_relations(points, events, cp)
            after[:, start // 8 : end // 8], before[:, start // 8 : end // 8] = (
                future,
                past,
            )
            generation["classified_events"] = end
            if end in levels:
                with (folder / f"case_{index}_n{end}.npz").open("xb") as stream:
                    np.savez(
                        stream,
                        **selected,
                        after=after[:, : end // 8],
                        before=before[:, : end // 8],
                        measurement_n=end,
                    )
            if end % (64 * BATCH) == 0 or end == top:
                status(
                    state="RUNNING",
                    phase="generate",
                    case=index,
                    events=end,
                    target=top,
                )
    except (ValueError, ArithmeticError, MemoryError) as error:
        generation["failure"] = f"{type(error).__name__}: {error}"
    # Return pooled GPU memory before the worker process needs the device.
    cp.get_default_memory_pool().free_all_blocks()
    generation.update(
        margins_resolved_by_the_direct_solver=geometry.resolved_margins,
        wall_seconds=time.perf_counter() - wall,
        process_cpu_seconds=time.process_time() - cpu,
        peak_working_set_bytes=peak_working_set(),
    )
    if "landmarks_including_tips" in generation:
        generation["directional_predicates_classified"] = (
            2 * generation["landmarks_including_tips"] * generation["classified_events"]
        )
    write_json(folder / f"case_{index}_generation.json", generation)
    return generation, truth_path


def geometry_record(proper_time):
    import numpy as np

    record = {}
    for kind in ("orbit", "rest", "flat"):
        g = build(kind, proper_time)
        e, b = g.tidal(), g.magnetic()
        record[kind] = {
            "mass": g.mass,
            "apex_areal_radius": g.apex_areal_radius,
            "speed_relative_to_static_observer": g.velocity,
            "tip_proper_time": g.proper_time,
            "apex_isotropic_radius": g.apex,
            "tips": g.tips.tolist(),
            "bounding_ball_centre": g.base.tolist(),
            "bounding_ball_radius": g.radius,
            "true_electric": e.tolist(),
            "true_magnetic": b.tolist(),
            "true_electric_eigenvalues": np.linalg.eigvalsh(e).tolist(),
            "true_magnetic_eigenvalues": np.linalg.eigvalsh(b).tolist(),
            "true_full_tidal_norm": float(np.sqrt(np.sum(e * e) + np.sum(b * b))),
            "true_invariant": float(np.sum(e * e) - np.sum(b * b)),
        }
    return record


def run(destination, engineering=False):
    import cupy as cp
    import numpy as np
    import scipy
    from intrinsic_weyl_gpu_counts import selfcheck as gpu_selfcheck
    from intrinsic_weyl_joint import selfcheck_jet
    from schwarzschild_orbit import selfcheck as geometry_selfcheck

    if os.name != "nt":
        raise RuntimeError("The resource accounting of this runner uses Windows")
    if destination.exists():
        raise FileExistsError("Preserve every existing run")
    dirty = subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip()
    if dirty and not engineering:
        raise RuntimeError("Commit the protocol and runner before test data")
    stamp = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    here = Path(__file__).resolve().parent
    sources = [
        here / name
        for name in (
            Path(__file__).name,
            "schwarzschild_orbit.py",
            "schwarzschild_diamond.py",
            "intrinsic_weyl_bombelli_benchmark.py",
            "replay_intrinsic_weyl_threads.py",
            "intrinsic_weyl_joint.py",
            "intrinsic_weyl.py",
            "intrinsic_weyl_square.py",
            "intrinsic_weyl_gpu_counts.py",
            "intrinsic_weyl_high_count.py",
        )
    ]
    if not engineering:
        sources.append(ROOT / PROTOCOL)
    source_hashes = {
        str(p.relative_to(ROOT)).replace("\\", "/"): lf_sha(p) for p in sources
    }
    if engineering:
        cases, levels, proper_time = (
            ENGINEERING_CASES,
            ENGINEERING_LEVELS,
            ENGINEERING_PROPER_TIME,
        )
    else:
        cases, levels, proper_time = CASES, LEVELS, PROPER_TIME
    checked = {
        "generator": geometry_selfcheck(),
        "gpu_relation_bits": relation_selfcheck(proper_time),
        "jet_inverse": selfcheck_jet(),
        "gpu_second_moment": gpu_selfcheck(),
        "verdict": selfcheck(),
    }
    if shutil_free(destination.parent) < 40 * 2**30:
        raise RuntimeError("At least 40 GiB free storage is required")
    destination.mkdir(parents=True)
    started = time.time()

    def status(**fields):
        temporary = destination / ".status.tmp"
        temporary.write_text(
            json.dumps(
                {"pid": os.getpid(), "started_unix": started, **fields}, indent=2
            )
            + "\n"
        )
        os.replace(temporary, destination / "status.json")

    worker_env = os.environ.copy()
    for key in THREAD_VARIABLES:
        worker_env.pop(key, None)
    worker_env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    manifest = {
        "engineering_only": bool(engineering),
        "code_version": stamp,
        "worktree_clean_at_start": not dirty,
        "source_sha256_lf": source_hashes,
        "selfcheck": checked,
        "geometry": geometry_record(proper_time),
        "cases": [
            {"case": i, "kind": k, "seed": SEED + offset}
            for i, (k, offset) in enumerate(cases)
        ],
        "measurement_levels": levels,
        "selection_n": SELECTION_N,
        "sampling_chunk": CHUNK,
        "classification_batch": BATCH,
        "false_positive_level_per_test": FALSE_POSITIVE,
        "interval_coverage": COVERAGE,
        "thresholds": thresholds(),
        "worker_thread_environment": {
            key: worker_env.get(key) for key in THREAD_VARIABLES
        },
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "cupy": cp.__version__,
            "device": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
        },
        "scope": (
            "Engineering replay on a stronger field; no test verdict"
            if engineering
            else "Existence gates for the field and its magnetic part and an "
            "ordering gate for |E| > |B| on one circular-orbit interval; no "
            "accuracy, generality or superiority claim"
        ),
    }
    write_json(destination / "manifest.json", manifest)

    def invoke(archive, truth, output):
        with output.with_suffix(".log").open("xb") as log:
            child = subprocess.Popen(
                [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    "--worker",
                    str(archive),
                    "--truth",
                    str(truth),
                    "--output",
                    str(output),
                ],
                env=worker_env,
                stdout=log,
                stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW,
            )
            status(
                state="RUNNING", phase="fit", output=output.name, child_pid=child.pid
            )
            code = child.wait()
        if code:
            return {
                "failure": f"Worker exit {code}; see {output.with_suffix('.log').name}"
            }
        return json.loads(output.read_text())

    rows = []
    try:
        for index, (kind, offset) in enumerate(cases):
            seed = SEED + offset
            generation, truth = generate(
                destination, index, kind, proper_time, seed, levels, status
            )
            for n in levels:
                row = {"case": index, "kind": kind, "seed": seed, "measurement_n": n}
                archive = destination / f"case_{index}_n{n}.npz"
                if not archive.exists():
                    row["failure"] = generation.get(
                        "failure", "Measurement prefix not completed"
                    )
                else:
                    result = invoke(archive, truth, archive.with_suffix(".json"))
                    row.update({k: v for k, v in result.items() if k != "full"})
                    if "full" in result:
                        row["detection"] = result["full"]["detection"]
                        if "evaluation" in result["full"]:
                            row["evaluation"] = result["full"]["evaluation"]
                    row.update(
                        archive_sha256=sha(archive),
                        archive_bytes=archive.stat().st_size,
                    )
                write_json(destination / f"case_{index}_n{n}_result.json", row)
                rows.append(row)
                detection = row.get("detection", {})
                print(
                    json.dumps(
                        {
                            "case": index,
                            "kind": kind,
                            "n": n,
                            "q": detection.get("q"),
                            "q_magnetic": detection.get("q_magnetic"),
                            "invariant_interval": detection.get(
                                "invariant_interval_95"
                            ),
                            "relative_error": row.get("evaluation", {}).get(
                                "relative_full_tidal_error"
                            ),
                            "failure": row.get("failure"),
                            "seconds": row.get("wall_seconds"),
                            "peak_gib": row.get("peak_working_set_bytes", 0) / 2**30,
                        }
                    ),
                    flush=True,
                )
        summary = {
            **manifest,
            "results": rows,
            "elapsed_seconds": time.time() - started,
        }
        if not engineering:
            summary["level_verdicts"] = {
                gate: {
                    str(n): judge([r for r in rows if r["measurement_n"] == n])
                    for n in levels
                }
                for gate, judge in zip(
                    GATES,
                    (field_verdict, magnetic_verdict, unequal_verdict),
                    strict=True,
                )
            }
            summary["primary_claim"] = primary_claim(summary["level_verdicts"])
            summary["consistency_label"] = {
                str(n): consistency([r for r in rows if r["measurement_n"] == n])
                for n in levels
            }
        for name, digest in source_hashes.items():
            if lf_sha(ROOT / name) != digest:
                raise RuntimeError(f"A bound source changed during execution: {name}")
        write_json(destination / "summary.json", summary)
        status(
            state="COMPLETED",
            assignments=len(rows),
            level_verdicts=summary.get("level_verdicts"),
            primary_claim=summary.get("primary_claim"),
        )
    except BaseException as error:
        status(state="STOPPED", error=f"{type(error).__name__}: {error}")
        raise


def compare(first, second):
    """Engineering replays must agree: exact order bits, 1e-8 on every fit field."""
    import numpy as np

    report = {"first": str(first), "second": str(second), "levels": {}}
    for n in ENGINEERING_LEVELS:
        name = f"case_0_n{n}"
        with (
            np.load(first / f"{name}.npz", allow_pickle=False) as a,
            np.load(second / f"{name}.npz", allow_pickle=False) as b,
        ):
            bits = set(a.files) == set(b.files) and all(
                np.array_equal(a[key], b[key]) for key in a.files
            )
        old = json.loads((first / f"{name}.json").read_text())["full"]
        new = json.loads((second / f"{name}.json").read_text())["full"]
        # The zip container stores timestamps; the order bits are compared above.
        for record in (old, new):
            record["fit"].pop("order_archive_sha256")
        report["levels"][str(n)] = {
            "order_archives_bit_exact": bool(bits),
            "fit_and_evaluation_differences": differences(new, old),
        }
    report["reproduced"] = all(
        level["order_archives_bit_exact"]
        and not level["fit_and_evaluation_differences"]
        for level in report["levels"].values()
    )
    write_json(second / "replay_comparison.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--engineering", action="store_true")
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--truth", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--compare", type=Path, nargs=2)
    args = parser.parse_args()
    if args.compare:
        outcome = compare(*args.compare)
        print(json.dumps({"reproduced": outcome["reproduced"]}))
        sys.exit(0 if outcome["reproduced"] else 1)
    elif args.check:
        print(json.dumps({"verdict": selfcheck()}))
    elif args.worker:
        if args.truth is None or args.output is None:
            parser.error("truth and output are required for a worker")
        worker(args.worker, args.truth, args.output)
    else:
        if args.destination is None:
            parser.error("destination is required")
        run(args.destination, engineering=args.engineering)
