"""Order-only frame link between worldline windows (deterministic).

Development diagnostic: no pre-specified criterion, no sampling noise. The
multi-window study (window_study_results.md, code d4ebd51) removed the
fifth-degree bias with windows along the observer's worldline, but aligned the
windows' frames with the truth; misaligning the side windows by one degree
already dominated the corrected bias.

Here every window off the centre shares landmarks with its neighbour toward
the centre: the neighbour's benchmark landmarks that lie inside it, with both
of its tip clocks above the selection floor, join its own landmark set. The
rotation between two windows' gauge frames is the proper Procrustes rotation
of the shared landmarks' inferred spatial coordinates; windows two steps out
are chained through the middle ones. Fits, alias map, principal null
direction and frame links then use only inferred coordinates, clock counts
and clock covariances: the corrected estimator is order-only end to end. The
truth enters only the evaluation.

Measured, per benchmark case and wavelength:
1. Link angle errors against the oracle relative rotation, and the part that
   curvature alone produces (Procrustes on the true normal coordinates).
2. Corrected biases with linked windows, with links from the true normal
   coordinates, and with oracle alignment.
3. The statistical link error at N = 2^24 and its effect on the corrected 95%
   radius, by sampling each window's linearized covariance of its Weyl
   coefficients and linked-landmark coordinates (independent windows).
4. The link errors when only the first 4, 8 or 16 shared landmarks of each
   link (in selection order) enter the Procrustes rotation.

Setting (idealized, as the landmark study's direct-placement mode): a landmark
pair is a clock when its flat expected selection count 4096 tau^4 is at least
two, and every landmark keeps its two tip clocks; each window's clock
covariance is the flat expected covariance at N = 2^24 from a scrambled Sobol
sample of 2^20 events (on the GPU when CuPy is present; checked against the
CPU harness). Masks, clock covariances and linearized covariances are those of
the half-cycle windows for every wavelength. The oracle-aligned baseline is
recomputed in this setting.
"""

import argparse
import json
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
REFERENCE_N = 2**24
SELECTION_N = 4096
TAU_MIN = (2 / SELECTION_N) ** 0.25
SHIFTS = (-1.0, -0.5, 0.0, 0.5, 1.0)
WAVELENGTHS = (2.0, 3.0, 4.0, 6.0)
QMC_LOG2 = 20
QMC_SEED = 2026112101
DRAWS = 4000
SUBSETS = (4, 8, 16)
PROBED = (
    ("single_window", 0.0, "none"),
    ("three_point_first", 0.5, "plane_wave_estimated_k"),
    ("three_point_first", 0.5, "minimum_norm_time_time"),
    ("five_point_second", 0.5, "plane_wave_estimated_k"),
    ("five_point_second", 0.5, "minimum_norm_time_time"),
)


# ---------------------------------------------------------------------------
# Layouts, clocks and flat covariances
# ---------------------------------------------------------------------------


def inner_neighbour(shift):
    """The window one step toward the centre, or None for the centre."""
    return None if shift == 0 else shift - 0.5 * float(np.sign(shift))


def build_layouts(points, amplitude):
    """Landmark sets of the five windows, each in its own Brinkmann chart.

    Every window keeps the benchmark landmarks at the same chart coordinates.
    A window off the centre also takes the benchmark landmarks of its inner
    neighbour that lie strictly inside it and whose flat proper times to both
    of its tips clear 1.02 times the selection floor; they are appended before
    the last tip. A benchmark landmark has the same index in every window.
    """
    from intrinsic_weyl_window_study import build_window

    tips = points[[0, -1]]
    layouts = {}
    for shift in SHIFTS:
        neighbour = inner_neighbour(shift)
        layout = {"neighbour": neighbour, "shared_neighbour": [], "shared_own": []}
        if neighbour is None:
            layout["points"] = points.copy()
            layouts[shift] = layout
            continue
        wave = build_window("sinusoid", amplitude, 2.0, shift)
        move = np.array([neighbour - shift, shift - neighbour, 0.0, 0.0])
        candidates = points[1:-1] + move
        index = np.flatnonzero(
            (candidates[:, 0] > 1e-6) & (candidates[:, 0] < 1 - 1e-6)
        )
        later, _ = wave.relations(tips[:1], candidates[index])
        earlier, _ = wave.relations(candidates[index], tips[1:])
        index = index[np.asarray(later) & np.asarray(earlier)]
        normal = wave.normal_coordinates(
            np.vstack((tips[:1], candidates[index], tips[1:]))
        )
        inner = normal[1:-1]
        radius2 = (inner[:, 1:] ** 2).sum(axis=1)
        floor = (1.02 * TAU_MIN) ** 2
        keep = ((inner[:, 0] + 0.5) ** 2 - radius2 >= floor) & (
            (0.5 - inner[:, 0]) ** 2 - radius2 >= floor
        )
        index = index[keep]
        layout["points"] = np.vstack((points[:-1], candidates[index], points[-1:]))
        layout["shared_neighbour"] = (index + 1).tolist()
        layout["shared_own"] = list(
            range(len(points) - 1, len(points) - 1 + len(index))
        )
        layouts[shift] = layout
    return layouts


def clock_mask(normal):
    """Selection-floor clocks from flat proper times of normal coordinates.

    A landmark pair is a clock when 4096 tau^4 >= 2; every landmark keeps its
    two tip clocks; the tip-to-tip pair is not a clock.
    """
    dt = normal[None, :, 0] - normal[:, None, 0]
    dx2 = ((normal[None, :, 1:] - normal[:, None, 1:]) ** 2).sum(axis=-1)
    tau2 = np.where(dt > 0, dt * dt - dx2, -1.0)
    tau2 = np.maximum(tau2, tau2.T)
    mask = tau2 >= TAU_MIN**2
    tip = np.zeros_like(mask)
    tip[0, 1:-1] = tip[1:-1, -1] = True
    tip |= tip.T
    mask |= tip & (tau2 > 0)
    np.fill_diagonal(mask, False)
    mask[0, -1] = mask[-1, 0] = False
    if not np.all(mask[0, 1:-1] & mask[1:-1, -1]):
        raise ValueError("A landmark lacks a tip clock")
    return mask


def clock_ends(normal, mask):
    """Earlier and later point of every clock, in fit_joint's pair order."""
    upper = np.triu(mask, 1)
    upper[0, -1] = False
    pairs = np.column_stack(np.nonzero(upper))
    forward = normal[pairs[:, 1], 0] > normal[pairs[:, 0], 0]
    return (
        np.where(forward, pairs[:, 0], pairs[:, 1]),
        np.where(forward, pairs[:, 1], pairs[:, 0]),
    )


def gpu_membership_moments(
    points, a, b, log2_size=QMC_LOG2, seed=QMC_SEED, chunk=32768
):
    """The landmark harness's flat_membership_moments on the GPU (same events).

    Products of 0/1 memberships summed over at most 2^15 events are exact in
    float32 (TF32 stays off), so the moments are exact integers as on the CPU.
    """
    os.environ["CUPY_TF32"] = "0"
    import cupy as cp
    from intrinsic_weyl_landmark_design import diamond_sobol

    events = diamond_sobol(log2_size, seed)
    p, ia, ib = cp.asarray(points), cp.asarray(a), cp.asarray(b)
    first = cp.zeros(len(a))
    second = cp.zeros((len(a), len(a)))
    for start in range(0, len(events), chunk):
        e = cp.asarray(events[start : start + chunk])
        dt = e[None, :, 0] - p[:, None, 0]
        dx2 = ((e[None, :, 1:] - p[:, None, 1:]) ** 2).sum(axis=-1)
        light = dt * dt > dx2
        member = ((light & (dt > 0))[ia] & (light & (dt < 0))[ib]).astype(cp.float32)
        first += member.sum(axis=1)
        second += member @ member.T
    return cp.asnumpy(first), cp.asnumpy(second), len(events)


def membership_moments(points, a, b):
    try:
        return gpu_membership_moments(points, a, b), "gpu"
    except ImportError:
        from intrinsic_weyl_landmark_design import flat_membership_moments

        return flat_membership_moments(points, a, b, QMC_LOG2, QMC_SEED), "cpu"


def prepare_case(case, archive, cache, amplitude):
    """Layouts, clock masks and flat clock covariances of one case's windows."""
    from intrinsic_weyl_landmark_design import flat_clock_covariance
    from intrinsic_weyl_window_study import build_window

    with np.load(Path(archive) / f"truth_{case}.npz") as saved:
        points = saved["points"]
    layouts = build_layouts(points, amplitude)
    report = {}
    for shift, layout in layouts.items():
        path = Path(cache) / f"layout_case{case}_shift{shift:+.1f}.npz"
        wave = build_window("sinusoid", amplitude, 2.0, shift)
        normal = wave.normal_coordinates(layout["points"])
        mask = clock_mask(normal)
        a, b = clock_ends(normal, mask)
        begin = time.time()
        (first, second, size), device = membership_moments(normal, a, b)
        covariance, p = flat_clock_covariance(first, second, size, REFERENCE_N)
        np.savez(
            path,
            points=layout["points"],
            mask=mask,
            covariance=covariance,
            shared_neighbour=np.array(layout["shared_neighbour"], dtype=int),
            shared_own=np.array(layout["shared_own"], dtype=int),
        )
        exact = (
            (normal[b, 0] - normal[a, 0]) ** 2
            - ((normal[b, 1:] - normal[a, 1:]) ** 2).sum(1)
        ) ** 2
        report[f"{shift:+.1f}"] = {
            "landmarks": len(layout["points"]) - 2,
            "shared_with_inner_neighbour": len(layout["shared_own"]),
            "clocks": len(a),
            "qmc_device": device,
            "qmc_seconds": round(time.time() - begin, 1),
            "clock_probability_range": [float(p.min()), float(p.max())],
            "qmc_vs_exact_flat_probability_max_relative": float(
                np.max(abs(p / exact - 1))
            ),
        }
        print(case, shift, report[f"{shift:+.1f}"], flush=True)
    report["chart_move_volume_check"] = chart_move_check(layouts, amplitude)
    return report


def chart_move_check(layouts, amplitude):
    """Shared landmarks are the same events in both windows' charts.

    The exact squared interval of every shared pair with distinct u, from each
    window's own wave and chart, must agree (largest difference relative to
    the largest squared interval). The tolerance 1e-10 sits above the
    generator's float error budget (1e-13 per stored fundamental value) and
    far below any curvature effect.
    """
    from intrinsic_weyl_second_continuum import exact_squared_interval, oriented_pairs
    from intrinsic_weyl_window_study import build_window

    worst = 0.0
    for shift in (0.5, -0.5, 1.0, -1.0):
        layout = layouts[shift]
        neighbour = layout["neighbour"]
        own = layout["points"][layout["shared_own"]]
        theirs = layouts[neighbour]["points"][layout["shared_neighbour"]]
        pairs = oriented_pairs(own)
        mine = exact_squared_interval(
            build_window("sinusoid", amplitude, 2.0, shift), own, pairs
        )
        other = exact_squared_interval(
            build_window("sinusoid", amplitude, 2.0, neighbour), theirs, pairs
        )
        worst = max(worst, float(np.max(abs(mine - other)) / np.max(abs(mine))))
    if worst > 1e-10:
        raise ArithmeticError("Shared landmarks are not the same events in both charts")
    return worst


# ---------------------------------------------------------------------------
# Window fits and linearized covariances
# ---------------------------------------------------------------------------


def fit_job(job):
    """Continuum fifth-degree fit of one window at one wavelength."""
    from intrinsic_weyl_joint import gauge
    from intrinsic_weyl_landmark_design import continuum_fit
    from intrinsic_weyl_window_study import build_window, gauge_frame

    begin = time.time()
    with np.load(job["layout"]) as saved:
        points, mask, covariance = saved["points"], saved["mask"], saved["covariance"]
    wave = build_window("sinusoid", job["amplitude"], job["wavelength"], job["shift"])
    normal = wave.normal_coordinates(points)
    volume = wave.volumes(points)
    upper = np.triu(mask, 1)
    upper[0, -1] = False
    pairs = np.column_stack(np.nonzero(upper))
    if not np.all(volume[pairs[:, 0], pairs[:, 1]] > 0):
        raise ArithmeticError("A clock has no exact wave volume")
    start = gauge(normal)
    fit = continuum_fit(
        REFERENCE_N * volume / volume[0, -1], mask, covariance, 5, start
    )
    return {
        "case": job["case"],
        "wavelength": job["wavelength"],
        "shift": job["shift"],
        "coefficients": fit["coefficients"],
        "covariance": fit["nominal_component_covariance"],
        "coordinates": fit["coordinates_inferred_from_counts"],
        "normal": normal.tolist(),
        "frame": gauge_frame(normal).tolist(),
        "optimizer_success": fit["optimizer_success"],
        "function_evaluations": fit["function_evaluations"],
        "projected_curvature_rank": fit["projected_curvature_rank"],
        "seconds": round(time.time() - begin, 1),
    }


def linked_landmarks(layouts_shared, shift):
    """Landmarks of window ``shift`` that enter any frame link."""
    own = set(layouts_shared[shift]["shared_own"])
    for other, shared in layouts_shared.items():
        if inner_neighbour(other) == shift:
            own |= set(shared["shared_neighbour"])
    return sorted(own)


def covariance_job(job):
    """Linearized covariance of Weyl and linked-landmark coordinates; alias map.

    At each wavelength's inferred coordinates of the window, zero curvature,
    with the window's whitening. Gauge-fixed coordinates have no row (-1).
    """
    from intrinsic_weyl_window_study import fifth_degree_design, whitening

    begin = time.time()
    with np.load(job["layout"]) as saved:
        mask, covariance = saved["mask"], saved["covariance"]
    whiten = whitening(covariance)
    result = {"case": job["case"], "shift": job["shift"], "blocks": {}, "alias": {}}
    for wavelength, coordinates in job["coordinates"].items():
        design = fifth_degree_design(np.asarray(coordinates), mask, whiten)
        inverse = np.linalg.pinv(design["j5"])
        dimension = design["dimension"]
        rows = list(range(dimension, dimension + 10))
        slots = []
        for landmark in job["linked"]:
            for component in (1, 2, 3):
                column = int(design["index"][landmark, component])
                slots.append(len(rows) if column >= 0 else -1)
                if column >= 0:
                    rows.append(column)
        influence = inverse[rows]
        result["blocks"][wavelength] = {
            "covariance": influence @ influence.T,
            "slots": slots,
        }
        if job["shift"] == 0.0:
            result["alias"][wavelength] = (inverse @ design["jd"])[
                dimension : dimension + 10
            ]
    result["seconds"] = round(time.time() - begin, 1)
    return result


# ---------------------------------------------------------------------------
# Links and corrections
# ---------------------------------------------------------------------------


def procrustes(target, source):
    """Proper rotation Q minimizing |target - source @ Q| (rows are points)."""
    u, _, vt = np.linalg.svd(source.T @ target)
    sign = np.sign(np.linalg.det(u @ vt))
    return u @ np.diag([1.0, 1.0, sign]) @ vt


def rotation_angle(rotation):
    return float(np.degrees(np.arccos(np.clip((np.trace(rotation) - 1) / 2, -1, 1))))


def rotation_vector(rotation):
    """Small-angle rotation vector of a rotation matrix (radians)."""
    return (
        np.array(
            [
                rotation[2, 1] - rotation[1, 2],
                rotation[0, 2] - rotation[2, 0],
                rotation[1, 0] - rotation[0, 1],
            ]
        )
        / 2
    )


def links(coordinates, shared):
    """Rotations Q[shift] with x_centre ~ x_shift @ Q[shift], chained outward.

    ``coordinates[shift]`` are a window's points (spatial part used);
    ``shared[shift]`` holds the shared index lists of the window.
    """
    result = {0.0: np.eye(3)}
    for shift in (0.5, -0.5, 1.0, -1.0):
        neighbour = inner_neighbour(shift)
        target = coordinates[neighbour][shared[shift]["shared_neighbour"], 1:]
        source = coordinates[shift][shared[shift]["shared_own"], 1:]
        result[shift] = procrustes(target, source) @ result[neighbour]
    return result


def first_shared(shared, count):
    """The first ``count`` shared landmarks of every link, in selection order."""
    return {
        shift: {key: value[key][:count] for key in ("shared_neighbour", "shared_own")}
        for shift, value in shared.items()
    }


def link_residuals(coordinates, shared, rotations):
    """RMS Procrustes residual of each step, in units of the shared radius."""
    result = {}
    for shift in (0.5, -0.5, 1.0, -1.0):
        neighbour = inner_neighbour(shift)
        target = coordinates[neighbour][shared[shift]["shared_neighbour"], 1:]
        source = coordinates[shift][shared[shift]["shared_own"], 1:]
        step = rotations[shift] @ rotations[neighbour].T
        residual = target - source @ step
        result[shift] = float(
            np.sqrt((residual**2).sum(1).mean() / (target**2).sum(1).mean())
        )
    return result


def analyse(case, wavelength, fits, blocks, alias, shared, amplitude, rng):
    """Links, corrections and the statistical link error at one wavelength."""
    from intrinsic_weyl import second_derivative_basis
    from intrinsic_weyl_window_study import (
        correction_specs,
        covector,
        evaluate_specs,
        minimum_norm_operator,
        principal_null_direction,
        radius_95,
        tidal_norm,
        unit_truth,
        weyl_rotation,
    )

    frames = {s: np.asarray(fits[s]["frame"]) for s in SHIFTS}
    fitted = {s: np.asarray(fits[s]["coordinates"]) for s in SHIFTS}
    normal = {s: np.asarray(fits[s]["normal"]) for s in SHIFTS}
    estimate = {s: np.asarray(fits[s]["coefficients"]) for s in SHIFTS}
    covariance = {s: np.asarray(fits[s]["covariance"]) for s in SHIFTS}
    frame0 = frames[0.0]
    truth = weyl_rotation(frame0) @ (amplitude * unit_truth())
    oracle = {s: frames[s].T @ frame0 for s in SHIFTS}
    linked = links(fitted, shared)
    # Links from the true normal coordinates, in the windows' truth gauge frames.
    gauged = {s: np.c_[normal[s][:, :1], normal[s][:, 1:] @ frames[s]] for s in SHIFTS}
    geometric = links(gauged, shared)
    result = {"case": case, "wavelength": wavelength, "links": {}}
    residuals = link_residuals(fitted, shared, linked)
    for shift in (0.5, -0.5, 1.0, -1.0):
        result["links"][f"{shift:+.1f}"] = {
            "shared_landmarks": len(shared[shift]["shared_own"]),
            "geometric_error_deg": rotation_angle(geometric[shift] @ oracle[shift].T),
            "total_error_deg": rotation_angle(linked[shift] @ oracle[shift].T),
            "procrustes_relative_residual": residuals[shift],
        }
    # Fewer shared landmarks: the first m of each list (selection order) only.
    subsets = {m: first_shared(shared, m) for m in SUBSETS}
    result["subsets"] = {}
    for m, chosen in subsets.items():
        rotations = links(fitted, chosen)
        result["subsets"][str(m)] = {
            "total_error_deg": {
                f"{shift:+.1f}": rotation_angle(rotations[shift] @ oracle[shift].T)
                for shift in (0.5, -0.5, 1.0, -1.0)
            }
        }
    # Corrections in the centre window's fitted gauge frame. The wave covector
    # (1, 0, 0, 1) of the generator frame and the propagation direction
    # (0, 0, -1) are carried into that frame for the oracle and the check.
    k_oracle = np.r_[1.0, np.array([0.0, 0.0, 1.0]) @ frame0]
    n_hat, _ = principal_null_direction(estimate[0.0])
    k_hat = covector(n_hat)
    result["principal_null_direction_error_deg"] = float(
        np.degrees(np.arccos(np.clip(-n_hat @ k_oracle[1:], -1, 1)))
    )
    minimum_norm, _ = minimum_norm_operator(alias @ second_derivative_basis())
    specs = correction_specs(list(SHIFTS), alias, k_oracle, k_hat, minimum_norm)

    def window_set(rotations):
        mapped = {}
        for s in SHIFTS:
            turn = weyl_rotation(rotations[s])
            mapped[s] = (turn @ estimate[s], turn @ covariance[s] @ turn.T)
        return mapped

    kept = (
        "correction",
        "h",
        "operator",
        "windows",
        "relative_bias",
        "relative_radius_95",
        "radius_ratio_to_single",
        "events_ratio_to_single",
    )
    result["corrections"] = {}
    for label, rotations in (
        ("oracle", oracle),
        ("true_coordinate_link", geometric),
        ("order_only_link", linked),
    ):
        rows = evaluate_specs(window_set(rotations), specs, truth)
        result["corrections"][label] = [{k: r[k] for k in kept} for r in rows]
    # Statistical link error: sample each window's linearized covariance of its
    # Weyl coefficients and linked-landmark coordinates; windows independent.
    probed = [s for s in specs if (s[0], s[1], s[2]) in PROBED]
    roots, targets = {}, {}
    for s in SHIFTS:
        values, vectors = np.linalg.eigh(blocks[s]["covariance"])
        roots[s] = vectors * np.sqrt(np.clip(values, 0, None))
        landmarks = np.repeat(linked_landmarks(shared, s), 3)
        components = np.tile([1, 2, 3], len(landmarks) // 3)
        slot = np.asarray(blocks[s]["slots"])
        free = slot >= 0
        targets[s] = (landmarks[free], components[free], slot[free])
    linked_draws = [[] for _ in probed]
    oracle_draws = [[] for _ in probed]
    angle_draws = {shift: [] for shift in (0.5, -0.5, 1.0, -1.0)}
    oracle_turns = {s: weyl_rotation(oracle[s]) for s in SHIFTS}
    subset_links = {m: links(fitted, chosen) for m, chosen in subsets.items()}
    subset_draws = {m: {shift: [] for shift in (0.5, -0.5, 1.0, -1.0)} for m in subsets}
    for _ in range(DRAWS):
        perturbed, coefficients = {}, {}
        for s in SHIFTS:
            delta = roots[s] @ rng.normal(size=roots[s].shape[1])
            coefficients[s] = estimate[s] + delta[:10]
            points = fitted[s].copy()
            landmarks, components, slot = targets[s]
            points[landmarks, components] += delta[slot]
            perturbed[s] = points
        rotations = links(perturbed, shared)
        for shift, draws in angle_draws.items():
            draws.append(rotation_vector(rotations[shift] @ linked[shift].T))
        for m, chosen in subsets.items():
            fewer = links(perturbed, chosen)
            for shift, draws in subset_draws[m].items():
                draws.append(rotation_vector(fewer[shift] @ subset_links[m][shift].T))
        mapped = {s: weyl_rotation(rotations[s]) @ coefficients[s] for s in SHIFTS}
        fixed = {s: oracle_turns[s] @ coefficients[s] for s in SHIFTS}
        for i, (_, _, _, chosen, weights) in enumerate(probed):
            pairs = list(zip(chosen, weights, strict=True))
            linked_draws[i].append(sum(w @ mapped[s] for s, w in pairs))
            oracle_draws[i].append(sum(w @ fixed[s] for s, w in pairs))
    truth_norm = tidal_norm(truth)
    result["statistical"] = {
        "draws": DRAWS,
        "link_rms_angle_deg_at_reference_n": {
            f"{shift:+.1f}": float(
                np.degrees(np.sqrt((np.array(v) ** 2).sum(1).mean()))
            )
            for shift, v in angle_draws.items()
        },
        "subset_link_rms_angle_deg_at_reference_n": {
            str(m): {
                f"{shift:+.1f}": float(
                    np.degrees(np.sqrt((np.array(v) ** 2).sum(1).mean()))
                )
                for shift, v in draws.items()
            }
            for m, draws in subset_draws.items()
        },
        "corrections": [],
    }
    for i, (name, h, operator, _, _) in enumerate(probed):
        with_link = radius_95(np.cov(np.array(linked_draws[i]).T))
        without = radius_95(np.cov(np.array(oracle_draws[i]).T))
        result["statistical"]["corrections"].append(
            {
                "correction": name,
                "h": h,
                "operator": operator,
                "relative_radius_95_oracle": without / truth_norm,
                "relative_radius_95_linked": with_link / truth_norm,
                "radius_ratio_linked_to_oracle": with_link / without,
            }
        )
    return result


# ---------------------------------------------------------------------------
# Selfcheck and driver
# ---------------------------------------------------------------------------


def selfcheck():
    from intrinsic_weyl_landmark_design import diamond_sobol, flat_membership_moments
    from intrinsic_weyl_window_study import axis_rotation, weyl_rotation

    rng = np.random.default_rng(2026112105)
    report = {}
    # Procrustes recovers a proper rotation exactly.
    rotation = axis_rotation(rng.normal(size=3), 0.7)
    source = rng.normal(size=(40, 3))
    recovered = procrustes(source @ rotation, source)
    report["procrustes_error"] = float(np.max(abs(recovered - rotation)))
    assert report["procrustes_error"] < 1e-12 and np.linalg.det(recovered) > 0
    # Weyl rotations compose as coordinate maps: first B, then A.
    first, second = (
        axis_rotation(rng.normal(size=3), 0.3),
        axis_rotation(rng.normal(size=3), 1.1),
    )
    report["weyl_rotation_composition_error"] = float(
        np.max(
            abs(
                weyl_rotation(second) @ weyl_rotation(first)
                - weyl_rotation(first @ second)
            )
        )
    )
    assert report["weyl_rotation_composition_error"] < 1e-12
    # Chained links recover the rotations; index i is the same event in every
    # window, as for benchmark landmarks.
    shared = {
        s: {"shared_neighbour": list(range(20 * k, 20 * k + 20))}
        for k, s in enumerate((0.5, -0.5, 1.0, -1.0))
    }
    for value in shared.values():
        value["shared_own"] = value["shared_neighbour"]
    truth = {s: axis_rotation(rng.normal(size=3), rng.uniform(0, 3)) for s in SHIFTS}
    truth[0.0] = np.eye(3)
    base = rng.normal(size=(80, 4))
    # x_s = x_centre @ truth[s]^T, so that x_centre = x_s @ truth[s].
    coordinates = {s: np.c_[base[:, :1], base[:, 1:] @ truth[s].T] for s in SHIFTS}
    found = links(coordinates, shared)
    report["chained_link_error"] = float(
        max(np.max(abs(found[s] - truth[s])) for s in SHIFTS)
    )
    assert report["chained_link_error"] < 1e-12
    # GPU moments equal the CPU harness exactly.
    try:
        import cupy  # noqa: F401

        events = diamond_sobol(12, 3)
        points = np.vstack(([-0.5, 0, 0, 0], events[:20], [0.5, 0, 0, 0]))
        dt = points[None, :, 0] - points[:, None, 0]
        dx2 = ((points[None, :, 1:] - points[:, None, 1:]) ** 2).sum(-1)
        mask = (dt * dt > dx2) & (dt > 0)
        a, b = np.nonzero(mask)
        cpu = flat_membership_moments(points, a, b, 14, 11)
        gpu = gpu_membership_moments(points, a, b, 14, 11)
        report["gpu_moments_equal_cpu"] = bool(
            np.array_equal(cpu[0], gpu[0]) and np.array_equal(cpu[1], gpu[1])
        )
        assert report["gpu_moments_equal_cpu"]
    except ImportError:
        report["gpu_moments_equal_cpu"] = None
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cases", type=int, nargs="+", default=[0, 1])
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--threads", type=int, default=3)
    parser.add_argument("--selfcheck-only", action="store_true")
    args = parser.parse_args()
    from intrinsic_weyl_bombelli_benchmark import benchmark_amplitude
    from intrinsic_weyl_landmark_design import code_state

    amplitude = benchmark_amplitude()
    commit, clean = code_state()
    report = selfcheck()
    print(json.dumps(report, indent=1), flush=True)
    if args.selfcheck_only:
        return
    args.cache.mkdir(parents=True, exist_ok=True)
    begin = time.time()
    preparation = {
        case: prepare_case(case, args.archive, args.cache, amplitude)
        for case in args.cases
    }
    preparation_seconds = round(time.time() - begin, 1)
    shared = {}
    for case in args.cases:
        shared[case] = {}
        for s in SHIFTS:
            with np.load(args.cache / f"layout_case{case}_shift{s:+.1f}.npz") as saved:
                shared[case][s] = {
                    "shared_neighbour": saved["shared_neighbour"].tolist(),
                    "shared_own": saved["shared_own"].tolist(),
                }
    for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ[name] = str(args.threads)
    jobs = [
        {
            "case": case,
            "wavelength": wavelength,
            "shift": s,
            "amplitude": str(amplitude),
            "layout": str(args.cache / f"layout_case{case}_shift{s:+.1f}.npz"),
        }
        for case in args.cases
        for wavelength in WAVELENGTHS
        for s in SHIFTS
    ]
    begin = time.time()
    fits = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for fit in pool.map(fit_job, jobs):
            fits[fit["case"], fit["wavelength"], fit["shift"]] = fit
            print(
                fit["case"],
                fit["wavelength"],
                fit["shift"],
                f"{fit['seconds']}s",
                fit["optimizer_success"],
                flush=True,
            )
        covariance_jobs = [
            {
                "case": case,
                "shift": s,
                "layout": str(args.cache / f"layout_case{case}_shift{s:+.1f}.npz"),
                "linked": linked_landmarks(shared[case], s),
                "coordinates": {
                    w: fits[case, w, s]["coordinates"] for w in WAVELENGTHS
                },
            }
            for case in args.cases
            for s in SHIFTS
        ]
        blocks = {}
        for item in pool.map(covariance_job, covariance_jobs):
            blocks[item["case"], item["shift"]] = item
            print(
                "covariance",
                item["case"],
                item["shift"],
                f"{item['seconds']}s",
                flush=True,
            )
    fit_seconds = round(time.time() - begin, 1)
    rng = np.random.default_rng(2026112106)
    rows = []
    for case in args.cases:
        for wavelength in WAVELENGTHS:
            rows.append(
                analyse(
                    case,
                    wavelength,
                    {s: fits[case, wavelength, s] for s in SHIFTS},
                    {s: blocks[case, s]["blocks"][wavelength] for s in SHIFTS},
                    blocks[case, 0.0]["alias"][wavelength],
                    shared[case],
                    float(amplitude),
                    rng,
                )
            )
            print(json.dumps(rows[-1]["links"]), flush=True)
    summary = {
        "diagnostic": (
            "order-only frame link between worldline windows; deterministic continuum; "
            "idealized clocks and flat covariances; no criterion"
        ),
        "code_version": commit,
        "worktree_clean": clean,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "benchmark_amplitude": str(amplitude),
        "reference_n": REFERENCE_N,
        "selfcheck": report,
        "preparation": {str(k): v for k, v in preparation.items()},
        "preparation_seconds": preparation_seconds,
        "fit_seconds": fit_seconds,
        "fits": [
            {
                k: fit[k]
                for k in (
                    "case",
                    "wavelength",
                    "shift",
                    "optimizer_success",
                    "function_evaluations",
                    "projected_curvature_rank",
                    "seconds",
                )
            }
            for fit in fits.values()
        ],
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
