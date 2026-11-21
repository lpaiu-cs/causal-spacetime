"""Finite-N test of the order-only five-window correction of the Weyl inverse.

The geometry is the Ibrahim-Bombelli half-cycle plane wave of the benchmark
(h0 = 1e-2). Five windows lie along the observer's geodesic x = y = 0 of the
Brinkmann chart, centred at proper-time offsets -1, -1/2, 0, 1/2, 1 of the tip
separation. The window at offset delta is the standard window of the profile
a(u + delta), and an event's chart coordinates in windows delta1 and delta2
differ by (delta1 - delta2, delta2 - delta1, 0, 0).

Each window is an independent sprinkle. Its landmarks are the unchanged nested
selection on its own 4096 selection events, plus, off the centre, the
landmarks that its inner neighbour's own selection chose, kept when at least
two of the window's own selection events lie between each tip and the
landmark (order-only). Clocks are the pairs with at least two of the window's
own selection events between them. Each window is fitted by the unchanged
fifth-degree inverse with curvature-square terms behind the blind boundary.

The correction uses only the five pure fits and the centre window's clock
covariance: O(3) Procrustes links of the shared landmarks' inferred spatial
coordinates (the order does not fix an embedding's handedness), the alias map
at the centre window's inferred coordinates, the minimum-norm operator and the
five-point second-order stencil at h = 1/2. It is saved before any generator
or truth module is imported in the analysing process.

See docs/paper/paper_a/intrinsic_weyl/window_test_protocol.md.
"""

import argparse
import hashlib
import itertools
import json
import math
import os
import platform
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from replay_intrinsic_weyl_threads import (
    ROOT,
    THREAD_VARIABLES,
    differences,
    sha,
    write_json,
)

# clock_second_moment relies on exact FP32 products: cuBLAS must not use TF32.
os.environ["CUPY_TF32"] = "0"
PROTOCOL = "docs/paper/paper_a/intrinsic_weyl/window_test_protocol.md"
SELECTION_N, CHUNK, BATCH = 4096, 4096, 65536
# Window index w -> proper-time offset of its centre, in tip separations.
WINDOWS = (-1.0, -0.5, 0.0, 0.5, 1.0)
CENTRE = 2
GENERATION_ORDER = (2, 1, 3, 0, 4)
INNER = {0: 1, 1: 2, 3: 2, 4: 3}
TIPS = np.array([[0.0, 0.5, 0.0, 0.0], [1.0, -0.5, 0.0, 0.0]])
WAVELENGTH = 2.0  # diamond heights: the crest-centred half cycle
LEVELS = (33554432,)
REPLICATES = 64
SEED = 2026112600
KEPT_ARCHIVES = (0, 31)
TOLERANCE = 0.05
COVERAGE = 0.95
MINIMUM_COMPLETE = 58  # ceil(0.9 * REPLICATES)
MINIMUM_SHARED = 8
ENGINEERING_LEVELS = (4194304, 16777216)
ENGINEERING_SEEDS = (2026112700, 2026112701)
# Generated windows whose fits have not finished; bounds disk and memory.
OUTSTANDING = 8
SINGLE = ("single_window", 0.0, "none")
PRIMARY = ("five_point_second", 0.5, "minimum_norm_time_time")
CORRECTIONS = {
    "single": SINGLE,
    "primary": PRIMARY,
    "five_point_plane_wave": ("five_point_second", 0.5, "plane_wave_estimated_k"),
    "three_point_plane_wave": ("three_point_first", 0.5, "plane_wave_estimated_k"),
    "three_point_minimum_norm": ("three_point_first", 0.5, "minimum_norm_time_time"),
}
ELECTRIC = np.r_[np.ones(5), np.zeros(5)]
MAGNETIC = np.r_[np.zeros(5), np.ones(5)]
BOOTSTRAP, BOOTSTRAP_SEED = 10000, 2026112601
ARCHIVE_FIELDS = {
    "selection_order",
    "selection_design",
    "selection_own",
    "candidate_source",
    "selected_ids",
    "after",
    "before",
    "measurement_n",
}
ORACLES = (
    "profile_plane_wave",
    "intrinsic_weyl_joint_audit",
    "intrinsic_weyl_bombelli_benchmark",
    "intrinsic_weyl_second_continuum",
    "intrinsic_weyl_joint_pilot",
    "intrinsic_weyl_rotating_pilot",
    "rotating_plane_wave",
    "window_response",
    "p14_plane_wave",
)


# ---------------------------------------------------------------------------
# Statistics and verdicts
# ---------------------------------------------------------------------------


def t_interval(values, coverage=COVERAGE):
    """Mean, sample sd and two-sided Student-t interval; None below two values."""
    from scipy.stats import t as student

    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return None
    mean, sd = float(values.mean()), float(values.std(ddof=1))
    quantile = float(student.ppf(0.5 + coverage / 2, len(values) - 1))
    half = quantile * sd / math.sqrt(len(values))
    return {
        "replicates": int(len(values)),
        "mean": mean,
        "sd": sd,
        "coverage": coverage,
        "t_quantile": quantile,
        "interval": [mean - half, mean + half],
    }


def accuracy_verdict(interval, tolerance=TOLERANCE):
    """G1: strict containment in [1 - t, 1 + t]; strict exclusion is INACCURATE."""
    low, high = interval
    if low > 1 - tolerance and high < 1 + tolerance:
        return "ACCURATE"
    if high < 1 - tolerance or low > 1 + tolerance:
        return "INACCURATE"
    return "UNRESOLVED"


def control_verdict(interval, tolerance=TOLERANCE):
    """G2: the single window's interval lies strictly below 1 - t."""
    return "BIASED_LOW" if interval[1] < 1 - tolerance else "NOT_BIASED"


def primary_claim(complete, accuracy, control, minimum=MINIMUM_COMPLETE):
    if complete < minimum:
        return "INCOMPLETE"
    if accuracy == "ACCURATE" and control == "BIASED_LOW":
        return "SUPPORTED"
    return "NOT_SUPPORTED"


def judge(rows, tolerance=TOLERANCE, minimum=MINIMUM_COMPLETE):
    """Both gates on the completed replicates of one level; failures are kept."""
    complete = [row for row in rows if "failure" not in row]
    corrected = t_interval([row["readings"]["primary"]["value"] for row in complete])
    single = t_interval([row["readings"]["single"]["value"] for row in complete])
    accuracy = accuracy_verdict(corrected["interval"], tolerance) if corrected else None
    control = control_verdict(single["interval"], tolerance) if single else None
    return {
        "assigned": len(rows),
        "complete": len(complete),
        "minimum_complete": minimum,
        "tolerance": tolerance,
        "corrected_reading": corrected,
        "single_window_reading": single,
        "g1_accuracy": accuracy,
        "g2_control": control,
        "primary_claim": primary_claim(len(complete), accuracy, control, minimum),
    }


def selfcheck():
    """The gates as frozen, including strictness and the completion rule."""

    def row(primary, single, failed=False):
        if failed:
            return {"failure": "kept"}
        return {
            "readings": {"primary": {"value": primary}, "single": {"value": single}}
        }

    assert MINIMUM_COMPLETE == math.ceil(0.9 * REPLICATES)
    pair = t_interval([0.0, 1.0])
    assert abs(pair["interval"][1] - (0.5 + 12.706204736174698 * 0.5)) < 1e-9
    assert accuracy_verdict([0.951, 1.049]) == "ACCURATE"
    assert accuracy_verdict([0.95, 1.0]) == "UNRESOLVED"
    assert accuracy_verdict([0.96, 1.05]) == "UNRESOLVED"
    assert accuracy_verdict([0.90, 0.9499]) == "INACCURATE"
    assert accuracy_verdict([0.90, 0.95]) == "UNRESOLVED"
    assert accuracy_verdict([1.0501, 1.1]) == "INACCURATE"
    assert accuracy_verdict([0.94, 1.06]) == "UNRESOLVED"
    assert control_verdict([0.80, 0.9499]) == "BIASED_LOW"
    assert control_verdict([0.80, 0.95]) == "NOT_BIASED"
    assert primary_claim(58, "ACCURATE", "BIASED_LOW") == "SUPPORTED"
    assert primary_claim(57, "ACCURATE", "BIASED_LOW") == "INCOMPLETE"
    assert primary_claim(64, "ACCURATE", "NOT_BIASED") == "NOT_SUPPORTED"
    assert primary_claim(64, "UNRESOLVED", "BIASED_LOW") == "NOT_SUPPORTED"
    rng = np.random.default_rng(2026112605)
    good = [row(1 + 0.01 * rng.normal(), 0.88 + 0.01 * rng.normal()) for _ in range(64)]
    assert judge(good)["primary_claim"] == "SUPPORTED"
    failed = good[:57] + [row(0, 0, failed=True)] * 7
    verdict = judge(failed)
    assert verdict["complete"] == 57 and verdict["primary_claim"] == "INCOMPLETE"
    assert verdict["g1_accuracy"] == "ACCURATE"
    return {
        "tolerance": TOLERANCE,
        "coverage_two_sided": COVERAGE,
        "minimum_complete": MINIMUM_COMPLETE,
        "strict_bounds": True,
        "failures_kept_and_counted": True,
    }


# ---------------------------------------------------------------------------
# Order-only landmarks and clocks of one window
# ---------------------------------------------------------------------------


def window_selection(order, own):
    """Landmarks, clock mask and kept candidates from an extended selection order.

    ``order`` covers [tip, own selection events, candidates, tip]. The own
    landmarks are the unchanged nested selection on the own events alone; a
    candidate is kept when at least two own events lie between each tip and it.
    Counts between landmarks use own events only. Kept candidates follow the
    own landmarks, before the last tip.
    """
    from intrinsic_weyl_joint import selection

    order = np.asarray(order, dtype=bool)
    total = len(order)
    if order.ndim != 2 or order.shape != (total, total) or total < own + 2:
        raise ValueError("Extended selection order has the wrong shape")
    if np.any(np.diag(order)) or np.any(order & order.T):
        raise ValueError("Strict selection order required")
    if not np.all(order[0, 1:]) or not np.all(order[:-1, -1]):
        raise ValueError("Selection events and candidates must lie in the window")
    own_index = np.r_[: own + 1, total - 1]
    ids, _ = selection(order[np.ix_(own_index, own_index)], design="nested")
    ids = own_index[ids]
    inner = np.arange(1, own + 1)

    def between(rows, columns):
        return order[np.ix_(rows, inner)].astype(float) @ order[
            np.ix_(inner, columns)
        ].astype(float)

    candidates = np.arange(own + 1, total - 1)
    low = between([0], candidates)[0]
    high = between(candidates, [total - 1])[:, 0]
    kept = candidates[(low >= 2) & (high >= 2)]
    combined = np.r_[ids[:-1], kept, ids[-1]]
    forward = between(combined, combined)
    counts = forward + forward.T
    clock_mask = counts >= 2
    clock_mask[0, -1] = clock_mask[-1, 0] = False
    if not np.all(clock_mask[0, 1:-1] & clock_mask[1:-1, -1]):
        raise ValueError("A landmark lacks a tip clock")
    return combined, clock_mask, kept


def clock_second_moment(after, before, a, b, chunk=1024):
    """The inverse's exact integer second moment, bounded in memory.

    The grouping of repeated membership patterns is that of
    intrinsic_weyl_joint._clock_second_moment. Memberships are formed on the
    GPU from one chunk of packed patterns at a time, and every weight is split
    into two parts of at most 2**13, so each FP32 product and partial sum of a
    1024-pattern chunk is an integer of at most 2**23 (exact; TF32 off); the
    chunks accumulate in float64. The result is therefore the same integers as
    the inverse's float64 paths, which moment_selfcheck verifies, while host
    transfers stay at one chunk. (Converting whole float64 membership blocks
    on the host, as the GPU helper of the inverse does, grew pinned staging
    memory to 35-65 GB per side window in the first pilot.)
    """
    import cupy as cp

    landmarks, bytes_per_landmark = after.shape
    n = 8 * bytes_per_landmark
    if n > 2**26:
        raise ValueError("The 13-bit weight split requires at most 2**26 events")
    signatures = np.empty((n, (2 * landmarks + 7) // 8), dtype=np.uint8)
    for start in range(0, bytes_per_landmark, 4096):
        end = min(start + 4096, bytes_per_landmark)
        bits = np.unpackbits(
            np.vstack((after[:, start:end], before[:, start:end])), axis=1
        )
        signatures[8 * start : 8 * end] = np.packbits(bits, axis=0).T
    width = signatures.shape[1]
    unique, weights = np.unique(
        signatures.view(np.dtype(f"V{width}")).ravel(), return_counts=True
    )
    del signatures
    unique = unique.view(np.uint8).reshape(-1, width)
    # np.packbits is big-endian: pattern row r sits in byte r // 8, bit 7 - r % 8.
    rows_a, rows_b = np.asarray(a), np.asarray(b) + landmarks
    byte_a, byte_b = cp.asarray(rows_a // 8), cp.asarray(rows_b // 8)
    shift_a = cp.asarray(7 - rows_a % 8, dtype=cp.uint8)
    shift_b = cp.asarray(7 - rows_b % 8, dtype=cp.uint8)
    totals = [cp.zeros((len(a), len(a)), dtype=cp.float64) for _ in range(2)]
    for start in range(0, len(unique), chunk):
        packed = cp.asarray(unique[start : start + chunk])
        member = (
            (packed[:, byte_a] >> shift_a) & (packed[:, byte_b] >> shift_b) & 1
        ).astype(cp.float32)
        weight = weights[start : start + chunk]
        for total, part in zip(totals, (weight >> 13, weight & 8191), strict=True):
            scaled = member * cp.asarray(part, dtype=cp.float32)[:, None]
            total += (scaled.T @ member).astype(cp.float64)
    second = cp.asnumpy(totals[0] * 8192 + totals[1])
    del totals, packed, member, scaled
    cp.get_default_memory_pool().free_all_blocks()
    cp.get_default_pinned_memory_pool().free_all_blocks()
    return second


def moment_selfcheck():
    """clock_second_moment equals both float64 paths of the inverse, bit for
    bit, on random patterns: few and many distinct patterns, weights beyond
    2**13 (the high half), and more patterns than one chunk."""
    from intrinsic_weyl_gpu_counts import clock_second_moment as gpu_float64
    from intrinsic_weyl_joint import _clock_second_moment as cpu_float64

    rng = np.random.default_rng(2026112606)
    cases = []
    # (landmarks, events, distinct patterns, share of events on pattern 0)
    for landmarks, n, distinct, dominant in (
        (13, 8192, 7, 0.0),
        (9, 512, 512, 0.0),
        (11, 256, 1, 0.0),
        (17, 2**21, 5, 0.0),
        (21, 2**16, 3000, 0.0),
        (12, 2**22, 2500, 0.6),
    ):
        patterns = rng.integers(0, 2, (2 * landmarks, distinct), dtype=np.uint8)
        index = rng.integers(0, distinct, n)
        index[rng.random(n) < dominant] = 0
        index[:distinct] = np.arange(distinct)
        bits = patterns[:, index]
        after = np.packbits(bits[:landmarks], axis=1)
        before = np.packbits(bits[landmarks:], axis=1)
        pairs = [(i, j) for i in range(landmarks) for j in range(landmarks) if i != j]
        a, b = (np.array(side) for side in zip(*pairs, strict=True))
        expected = cpu_float64(after, before, a, b)
        assert np.array_equal(clock_second_moment(after, before, a, b), expected)
        assert np.array_equal(gpu_float64(after, before, a, b), expected)
        cases.append(
            {
                "landmarks": landmarks,
                "events": n,
                "distinct_patterns": distinct,
                "largest_weight": int(np.bincount(index).max()),
            }
        )
    return {"bit_exact_against_cpu_and_gpu_float64": True, "cases": cases}


def read_window_clocks(path):
    """Clocks, their covariance and the shared rows, from the saved order bits.

    The measurement part is read_order_clocks of the inverse, with the
    bit-identical clock_second_moment.
    """
    with np.load(path, allow_pickle=False) as source:
        if set(source.files) != ARCHIVE_FIELDS:
            raise ValueError("Archive contains unexpected fields")
        order = source["selection_order"]
        own = int(source["selection_own"])
        origins = source["candidate_source"]
        ids = source["selected_ids"]
        after, before = source["after"], source["before"]
        n = int(source["measurement_n"])
        design = str(source["selection_design"].item())
    if design != "nested":
        raise ValueError("The nested selection design is required")
    if len(origins) != len(order) - own - 2:
        raise ValueError("Every candidate needs its source landmark")
    expected, clock_mask, _ = window_selection(order, own)
    if not np.array_equal(ids, expected):
        raise ValueError("Saved landmarks differ from the order-only selection")
    if after.shape != (len(ids), n // 8) or before.shape != after.shape or n % 8:
        raise ValueError("Packed order dimensions are inconsistent")
    if not (
        np.all(after[0] == 255)
        and np.all(before[-1] == 255)
        and np.all(before[0] == 0)
        and np.all(after[-1] == 0)
    ):
        raise ValueError("Measurement order is not contained in the outer interval")
    counts = np.zeros((len(ids), len(ids)))
    restricted = order[np.ix_(ids, ids)]
    for i, j in itertools.combinations(range(len(ids)), 2):
        if restricted[i, j]:
            a, b = i, j
        elif restricted[j, i]:
            a, b = j, i
        else:
            continue
        counts[i, j] = counts[j, i] = np.bitwise_count(after[a] & before[b]).sum()
    pairs = np.column_stack(np.nonzero(np.triu(clock_mask, 1)))
    forward = restricted[pairs[:, 0], pairs[:, 1]]
    a = np.where(forward, pairs[:, 0], pairs[:, 1])
    b = np.where(forward, pairs[:, 1], pairs[:, 0])
    p = counts[pairs[:, 0], pairs[:, 1]] / n
    if np.any(p <= 0):
        raise ValueError("A selected clock has zero measurement count")
    second = clock_second_moment(after, before, a, b)
    indicator_covariance = (second - n * np.outer(p, p)) / (n - 1)
    derivative = -0.5 / np.sqrt(p)
    clock_covariance = indicator_covariance * np.outer(derivative, derivative) / n
    rows = np.flatnonzero((ids > own) & (ids < len(order) - 1))
    return (
        counts,
        clock_mask,
        clock_covariance,
        ids,
        n,
        rows,
        origins[ids[rows] - own - 1],
    )


def flat_diamond(count, rng):
    """Uniform events in the flat unit-height diamond, coordinates (t, x, y, z)."""
    accepted, total = [], 0
    while total < count:
        t = rng.uniform(-0.5, 0.5, 4 * count)
        x = rng.uniform(-0.5, 0.5, (4 * count, 3))
        keep = np.linalg.norm(x, axis=1) < 0.5 - abs(t)
        accepted.append(np.column_stack((t, x))[keep])
        total += int(keep.sum())
    return np.vstack(accepted)[:count]


def flat_order(points, block=512):
    order = np.zeros((len(points), len(points)), dtype=bool)
    for start in range(0, len(points), block):
        p = points[start : start + block]
        dt = points[None, :, 0] - p[:, None, 0]
        dx2 = ((points[None, :, 1:] - p[:, None, 1:]) ** 2).sum(-1)
        order[start : start + block] = (dt > 0) & (dt * dt > dx2)
    return order


def selection_selfcheck():
    """window_selection equals the unchanged selection without candidates, and
    applies the shared rule and the own-event counts with them."""
    from intrinsic_weyl_joint import selection

    rng = np.random.default_rng(2026112602)
    own = np.vstack(([-0.5, 0, 0, 0], flat_diamond(SELECTION_N, rng), [0.5, 0, 0, 0]))
    extra = flat_diamond(48, rng)
    extended = np.vstack((own[:-1], extra, own[-1:]))
    order = flat_order(extended)
    total = len(order)
    own_index = np.r_[: SELECTION_N + 1, total - 1]
    base = order[np.ix_(own_index, own_index)]
    ids0, mask0 = selection(base, design="nested")
    ids, mask, kept = window_selection(base, SELECTION_N)
    assert np.array_equal(ids, ids0) and np.array_equal(mask, mask0) and not len(kept)
    ids, mask, kept = window_selection(order, SELECTION_N)
    inner = np.arange(1, SELECTION_N + 1)
    candidates = np.arange(SELECTION_N + 1, total - 1)
    low = order[np.ix_(inner, candidates)].sum(axis=0)
    high = order[np.ix_(candidates, inner)].sum(axis=1)
    assert np.array_equal(kept, candidates[(low >= 2) & (high >= 2)])
    assert np.array_equal(ids, np.r_[own_index[ids0[:-1]], kept, total - 1])
    direct = np.zeros((len(ids), len(ids)), dtype=int)
    for i, j in itertools.combinations(range(len(ids)), 2):
        direct[i, j] = direct[j, i] = int(
            np.sum(order[ids[i], inner] & order[inner, ids[j]])
            + np.sum(order[ids[j], inner] & order[inner, ids[i]])
        )
    expected = direct >= 2
    expected[0, -1] = expected[-1, 0] = False
    assert np.array_equal(mask, expected)
    return {
        "own_selection_unchanged": True,
        "candidates": len(candidates),
        "kept_by_count_rule": len(kept),
        "clock_mask_from_own_events": True,
    }


# ---------------------------------------------------------------------------
# Links and the order-only correction
# ---------------------------------------------------------------------------


def procrustes_o3(target, source):
    """Orthogonal Q, proper or improper, minimizing |target - source @ Q|."""
    u, _, vt = np.linalg.svd(source.T @ target)
    return u @ vt


def rotation_angle(rotation):
    return float(np.degrees(np.arccos(np.clip((np.trace(rotation) - 1) / 2, -1, 1))))


def link_windows(coordinates, shared, minimum=MINIMUM_SHARED):
    """O(3) maps Q[shift] with x_centre ~ x_shift @ Q[shift], chained outward.

    ``coordinates[shift]`` are a window's inferred points; ``shared[shift]``
    holds (rows in the window, rows of the same events in its inner neighbour).
    """
    rotations, record = {0.0: np.eye(3)}, {}
    for w in (1, 3, 0, 4):
        shift, neighbour = WINDOWS[w], WINDOWS[INNER[w]]
        rows, sources = (np.asarray(v, dtype=int) for v in shared[shift])
        if len(rows) < minimum:
            raise ValueError(
                f"Link {shift:+.1f} has {len(rows)} shared landmarks, "
                f"{minimum} required"
            )
        target = coordinates[neighbour][sources, 1:]
        source = coordinates[shift][rows, 1:]
        step = procrustes_o3(target, source)
        rotations[shift] = step @ rotations[neighbour]
        residual = target - source @ step
        record[f"{shift:+.1f}"] = {
            "shared_landmarks": int(len(rows)),
            "step_parity": int(round(np.linalg.det(step))),
            "parity_to_centre": int(round(np.linalg.det(rotations[shift]))),
            "procrustes_relative_residual": float(
                np.sqrt((residual**2).sum(1).mean() / (target**2).sum(1).mean())
            ),
        }
    return rotations, record


def correct(fits, mask, covariance):
    """The order-only corrections from five pure fits; no truth enters."""
    from intrinsic_weyl import second_derivative_basis
    from intrinsic_weyl_window_study import (
        alias_map,
        combine,
        correction_specs,
        covector,
        minimum_norm_operator,
        principal_null_direction,
        weyl_rotation,
        whitening,
    )

    coordinates = {
        s: np.asarray(fits[s]["coordinates_inferred_from_counts"]) for s in WINDOWS
    }
    shared = {
        s: (fits[s]["shared_rows"], fits[s]["shared_source_rows"])
        for s in WINDOWS
        if s != 0.0
    }
    rotations, links = link_windows(coordinates, shared)
    mapped = {}
    for s in WINDOWS:
        turn = weyl_rotation(rotations[s])
        estimate = np.asarray(fits[s]["coefficients"])
        nominal = np.asarray(fits[s]["nominal_component_covariance"])
        mapped[s] = (turn @ estimate, turn @ nominal @ turn.T)
    alias, consistency, _ = alias_map(coordinates[0.0], mask, whitening(covariance))
    minimum_norm, rank = minimum_norm_operator(alias @ second_derivative_basis())
    n_hat, residual = principal_null_direction(np.asarray(fits[0.0]["coefficients"]))
    specs = correction_specs(
        list(WINDOWS), alias, covector(n_hat), covector(n_hat), minimum_norm
    )
    corrections = {}
    for label, key in CORRECTIONS.items():
        _, _, _, shifts, weights = next(s for s in specs if s[:3] == key)
        estimate, nominal = combine([mapped[s] for s in shifts], weights)
        corrections[label] = {
            "key": list(key),
            "shifts": list(shifts),
            "weights": [np.asarray(w).tolist() for w in weights],
            "estimate": estimate.tolist(),
            "covariance": nominal.tolist(),
        }
    return {
        "links": {
            "steps": links,
            "rotations": {f"{s:+.1f}": rotations[s].tolist() for s in WINDOWS},
        },
        "alias_column_consistency": consistency,
        "minimum_norm_slice_rank": rank,
        "principal_null_direction": n_hat.tolist(),
        "principal_null_residual": residual,
        "corrections": corrections,
    }


def link_selfcheck():
    """O(3) Procrustes, chained links with mixed parity, Weyl maps of O(3)."""
    from intrinsic_weyl_joint import axis_signs
    from intrinsic_weyl_window_study import axis_rotation, weyl_rotation

    rng = np.random.default_rng(2026112604)
    mirror = np.diag([1.0, 1.0, -1.0])
    report = {}
    for label, matrix in (
        ("proper", axis_rotation(rng.normal(size=3), 0.7)),
        ("improper", axis_rotation(rng.normal(size=3), 1.3) @ mirror),
    ):
        source = rng.normal(size=(40, 3))
        report[f"{label}_procrustes_error"] = float(
            np.max(abs(procrustes_o3(source @ matrix, source) - matrix))
        )
        assert report[f"{label}_procrustes_error"] < 1e-12
    first = axis_rotation(rng.normal(size=3), 0.4) @ mirror
    second = axis_rotation(rng.normal(size=3), 2.1)
    report["weyl_map_composition_error"] = float(
        np.max(
            abs(
                weyl_rotation(second) @ weyl_rotation(first)
                - weyl_rotation(first @ second)
            )
        )
    )
    assert report["weyl_map_composition_error"] < 1e-12
    report["mirror_map_equals_axis_signs"] = float(
        np.max(abs(weyl_rotation(mirror) - np.diag(axis_signs([1, 1, -1]))))
    )
    assert report["mirror_map_equals_axis_signs"] < 1e-12
    truth = {0.0: np.eye(3)}
    for s in (0.5, -0.5, 1.0, -1.0):
        turn = axis_rotation(rng.normal(size=3), rng.uniform(0, 3))
        truth[s] = turn @ mirror if s in (0.5, -1.0) else turn
    events = {
        0.0: np.arange(0, 60),
        0.5: np.arange(30, 90),
        -0.5: np.r_[0:20, 90:130],
        1.0: np.arange(60, 110),
        -1.0: np.r_[95:130, 0:5],
    }
    base = rng.normal(size=(130, 4))
    coordinates = {
        s: np.c_[base[events[s], :1], base[events[s], 1:] @ truth[s].T] for s in WINDOWS
    }
    shared = {}
    for w in (1, 3, 0, 4):
        s, neighbour = WINDOWS[w], WINDOWS[INNER[w]]
        _, rows, sources = np.intersect1d(
            events[s], events[neighbour], return_indices=True
        )
        shared[s] = (rows, sources)
    found, steps = link_windows(coordinates, shared)
    report["chained_link_error"] = float(
        max(np.max(abs(found[s] - truth[s])) for s in WINDOWS)
    )
    assert report["chained_link_error"] < 1e-12
    assert {k: v["parity_to_centre"] for k, v in steps.items()} == {
        "-0.5": 1,
        "+0.5": -1,
        "-1.0": -1,
        "+1.0": 1,
    }
    shared[1.0] = (shared[1.0][0][:7], shared[1.0][1][:7])
    try:
        link_windows(coordinates, shared)
    except ValueError:
        report["link_below_eight_shared_fails"] = True
    else:
        raise AssertionError("A link with seven shared landmarks was accepted")
    return report


# ---------------------------------------------------------------------------
# Workers: pure fits and the correction precede every truth import
# ---------------------------------------------------------------------------


def guarded_load(allowed, opened):
    original = np.load

    def guard(selected, *args, **options):
        if Path(selected).resolve() not in allowed:
            raise ValueError("Only the assigned order-only input may be opened")
        opened.append(str(Path(selected).resolve()))
        return original(selected, *args, **options)

    return original, guard


def window_worker(archive, output, clocks=None):
    """Pure fit of one window; no generator or truth module is ever imported."""
    import intrinsic_weyl_joint as joint
    from intrinsic_weyl_high_count import peak_working_set

    wall, cpu = time.perf_counter(), time.process_time()
    record = {
        "archive_sha256": sha(archive),
        "thread_environment": {key: os.environ.get(key) for key in THREAD_VARIABLES},
    }
    opened = []
    original, guard = guarded_load({archive.resolve()}, opened)
    try:
        assert not any(name in sys.modules for name in ORACLES)
        np.load = guard
        try:
            counts, mask, covariance, ids, n, rows, origins = read_window_clocks(
                archive
            )
        finally:
            np.load = original
        record["clock_seconds"] = time.perf_counter() - wall
        fit = joint.fit_joint(
            counts,
            clock_mask=mask,
            clock_covariance=covariance,
            model_order=5,
            square_correction=True,
        )
        fit.update(
            selected_ids=ids.tolist(),
            measurement_n=n,
            smallest_measurement_clock=int(counts[mask].min()),
            order_archive_sha256=record["archive_sha256"],
            clocks=int(np.triu(mask, 1).sum()),
            shared_rows=rows.tolist(),
            shared_source_rows=origins.tolist(),
        )
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
                "archive_sha256": record["archive_sha256"],
            },
        )
        if clocks is not None:
            with clocks.open("xb") as stream:
                np.savez(stream, mask=mask, covariance=covariance)
        # The unchanged validity rule: success, ranks, positive covariance.
        joint.confidence_ball(fit)
        record["fit_summary"] = {
            key: fit[key]
            for key in (
                "landmarks",
                "measured_intervals",
                "smallest_measurement_clock",
                "optimizer_success",
                "function_evaluations",
                "coordinate_rank",
                "coordinate_parameters",
                "projected_curvature_rank",
                "boundary_excess",
            )
        }
        record["fit_summary"]["shared_landmarks"] = len(rows)
    except (ValueError, ArithmeticError, MemoryError, np.linalg.LinAlgError) as error:
        record["failure"] = f"{type(error).__name__}: {error}"
    record.update(
        wall_seconds=time.perf_counter() - wall,
        process_cpu_seconds=time.process_time() - cpu,
        peak_working_set_bytes=peak_working_set(),
    )
    write_json(output, record)


def combine_worker(folder, prefix, n, output):
    """The order-only correction of one replicate, saved before the truth step."""
    from intrinsic_weyl_high_count import peak_working_set

    wall = time.perf_counter()
    record = {}
    clocks = folder / f"{prefix}_win{CENTRE}_n{n}_clocks.npz"
    opened = []
    original, guard = guarded_load({clocks.resolve()}, opened)
    try:
        assert not any(name in sys.modules for name in ORACLES)
        fits, digests = {}, {}
        for w, shift in enumerate(WINDOWS):
            path = folder / f"{prefix}_win{w}_n{n}_pure.json"
            fits[shift] = json.loads(path.read_text(encoding="utf-8"))["fit"]
            digests[path.name] = sha(path)
        np.load = guard
        try:
            with np.load(clocks, allow_pickle=False) as saved:
                mask, covariance = saved["mask"], saved["covariance"]
        finally:
            np.load = original
        correction = correct(fits, mask, covariance)
        imported = [name for name in ORACLES if name in sys.modules]
        if imported:
            raise RuntimeError(f"Oracle imported before the correction: {imported}")
        write_json(
            output.with_name(output.stem + "_pure.json"),
            {
                "correction": correction,
                "pure_fit_sha256": digests,
                "opened_inputs": opened,
                "oracle_modules_before_return": imported,
            },
        )
        record["correction_seconds"] = time.perf_counter() - wall
        record["links"] = correction["links"]["steps"]
        record["evaluation"] = evaluate(folder, prefix, fits, correction)
    except (ValueError, ArithmeticError, MemoryError, np.linalg.LinAlgError) as error:
        record["failure"] = f"{type(error).__name__}: {error}"
    record.update(
        wall_seconds=time.perf_counter() - wall,
        peak_working_set_bytes=peak_working_set(),
    )
    write_json(output, record)


# ---------------------------------------------------------------------------
# Truth-side evaluation (after the correction has been saved)
# ---------------------------------------------------------------------------


def oracle_alignment(normal, inferred):
    """O(3) R, inferred ~ normal @ R on interior landmarks (aligned_tensor_error)."""
    left, _, right = np.linalg.svd(normal[1:-1, 1:].T @ inferred[1:-1, 1:])
    return left @ right


def evaluate(folder, prefix, fits, correction):
    """Readings against the truth in the centre window's inferred frame."""
    from intrinsic_weyl_bombelli_benchmark import benchmark_amplitude
    from intrinsic_weyl_joint import TIDAL_GRAM
    from intrinsic_weyl_joint_audit import aligned_tensor_error
    from intrinsic_weyl_window_study import (
        build_window,
        tidal_norm,
        unit_truth,
        weyl_rotation,
    )

    amplitude = benchmark_amplitude()
    value = float(amplitude)
    inferred = {
        s: np.asarray(fits[s]["coordinates_inferred_from_counts"]) for s in WINDOWS
    }

    def normal(w):
        with np.load(
            folder / f"truth_{prefix}_win{w}.npz", allow_pickle=False
        ) as saved:
            points = saved["points"]
        return build_window(
            "sinusoid", amplitude, WAVELENGTH, WINDOWS[w]
        ).normal_coordinates(points)

    centre = normal(CENTRE)
    aligned = aligned_tensor_error(
        {
            "coefficients": np.zeros(10),
            "coordinates_inferred_from_counts": inferred[0.0],
        },
        centre,
        value,
    )
    truth = np.asarray(aligned["aligned_truth_coefficients"])
    rotation0 = oracle_alignment(centre, inferred[0.0])
    generator_truth = value * unit_truth()
    check = float(np.max(abs(weyl_rotation(rotation0) @ generator_truth - truth)))
    if check > 1e-12:
        raise ArithmeticError("Aligned truth and the alignment rotation disagree")
    back = weyl_rotation(rotation0.T)

    def reading(estimate, covariance, part=None):
        mask = np.ones(10) if part is None else part
        target = truth * mask
        norm2 = float(target @ TIDAL_GRAM @ target)
        projected = TIDAL_GRAM @ target
        return {
            "value": float(projected @ (estimate * mask)) / norm2,
            "nominal_sd": float(
                np.sqrt(projected @ (mask[:, None] * covariance * mask) @ projected)
            )
            / norm2,
        }

    readings = {}
    for label, item in correction["corrections"].items():
        estimate = np.asarray(item["estimate"])
        covariance = np.asarray(item["covariance"])
        readings[label] = {
            **reading(estimate, covariance),
            "electric": reading(estimate, covariance, ELECTRIC),
            "magnetic": reading(estimate, covariance, MAGNETIC),
            "relative_full_tidal_error": tidal_norm(estimate - truth)
            / tidal_norm(truth),
            "generator_frame_estimate": (back @ estimate).tolist(),
        }
    result = {
        "truth_coefficients_centre_frame": truth.tolist(),
        "centre_alignment_parity": aligned["alignment_parity"],
        "truth_frame_check": check,
        "readings": readings,
    }
    # Descriptive: links against the oracle relative frames, and the primary
    # correction with oracle links. The normal frames of all windows are
    # parallel transports of one another, so the oracle map is R_s^T R_0.
    try:
        rotations = {CENTRE: rotation0}
        for w in (0, 1, 3, 4):
            rotations[w] = oracle_alignment(normal(w), inferred[WINDOWS[w]])
        oracle = {WINDOWS[w]: rotations[w].T @ rotation0 for w in rotations}
        linked = {
            s: np.asarray(correction["links"]["rotations"][f"{s:+.1f}"])
            for s in WINDOWS
        }
        links = {}
        for s in (0.5, -0.5, 1.0, -1.0):
            difference = linked[s] @ oracle[s].T
            agrees = bool(np.linalg.det(difference) > 0)
            links[f"{s:+.1f}"] = {
                "parity_agrees_with_oracle": agrees,
                "angle_error_deg": rotation_angle(difference) if agrees else None,
            }
        primary = correction["corrections"]["primary"]
        estimate = sum(
            np.asarray(weight)
            @ weyl_rotation(oracle[s])
            @ np.asarray(fits[s]["coefficients"])
            for s, weight in zip(primary["shifts"], primary["weights"], strict=True)
        )
        result["oracle_links"] = {
            "links": links,
            "window_parity_vs_truth": {
                f"{WINDOWS[w]:+.1f}": int(round(np.linalg.det(r)))
                for w, r in rotations.items()
            },
            "primary_reading_with_oracle_links": reading(estimate, np.zeros((10, 10)))[
                "value"
            ],
        }
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
        result["oracle_links_failure"] = f"{type(error).__name__}: {error}"
    return result


# ---------------------------------------------------------------------------
# Generation
# ---------------------------------------------------------------------------


def build_waves(amplitude):
    from intrinsic_weyl_window_study import build_window

    return {s: build_window("sinusoid", amplitude, WAVELENGTH, s) for s in WINDOWS}


def generate_window(folder, prefix, w, seed, wave, levels, status, neighbour):
    """Selection order with candidates, landmarks, measurement prefixes."""
    import cupy as cp
    from intrinsic_weyl_high_count import peak_working_set

    shift = WINDOWS[w]
    wall, cpu, resolved = (
        time.perf_counter(),
        time.process_time(),
        wave.resolved_margins,
    )
    generation = {
        "window": w,
        "shift": shift,
        "selection_events_drawn": 0,
        "sampled_events": 0,
        "classified_events": 0,
    }
    info = None
    try:
        status(state="RUNNING", phase="selection", replicate=prefix, window=w)
        rng = np.random.default_rng([seed, w, 0])
        points = np.vstack((TIPS[:1], wave.sample(SELECTION_N, rng), TIPS[1:]))
        generation["selection_events_drawn"] = SELECTION_N
        # Labels are randomized independently of coordinates; only the two
        # supplied tips have distinguished IDs.
        permutation = np.r_[
            0, rng.permutation(np.arange(1, len(points) - 1)), len(points) - 1
        ]
        points = points[permutation]
        candidates, origins = np.empty((0, 4)), np.empty(0, dtype=np.int64)
        if neighbour is not None:
            inner_points, own_rows, inner_shift = neighbour
            rows = np.arange(1, own_rows)
            moved = inner_points[rows] + np.array(
                [inner_shift - shift, shift - inner_shift, 0.0, 0.0]
            )
            charted = np.flatnonzero((moved[:, 0] > 0) & (moved[:, 0] < 1))
            later, _ = wave.relations(TIPS[:1], moved[charted])
            earlier, _ = wave.relations(moved[charted], TIPS[1:])
            inside = charted[np.asarray(later) & np.asarray(earlier)]
            candidates, origins = moved[inside], rows[inside]
            generation.update(
                candidates_offered=len(rows),
                candidates_outside_chart=int(len(rows) - len(charted)),
                candidates_outside_interval=int(len(charted) - len(inside)),
                candidates_in_interval=len(inside),
            )
        extended = np.vstack((points[:-1], candidates, points[-1:]))
        order = wave.order(extended)
        try:
            ids, _, kept = window_selection(order, SELECTION_N)
        except ValueError:
            np.savez(
                folder / f"{prefix}_win{w}_selection_failed.npz",
                selection_order=order,
                selection_own=np.array(SELECTION_N),
                candidate_source=origins,
            )
            raise
        selected = {
            "selection_order": order,
            "selection_design": np.array("nested"),
            "selection_own": np.array(SELECTION_N),
            "candidate_source": origins,
            "selected_ids": ids,
        }
        np.savez(
            folder / f"{prefix}_win{w}_selection.npz",
            selection_order_packed=np.packbits(order, axis=1),
            size=np.array(len(order)),
            selection_own=np.array(SELECTION_N),
            candidate_source=origins,
            selected_ids=ids,
        )
        landmarks = extended[ids]
        own_rows = int(np.sum(ids <= SELECTION_N))
        np.savez(
            folder / f"truth_{prefix}_win{w}.npz",
            points=landmarks,
            shift=np.array(shift),
            own_rows=np.array(own_rows),
            shared_source=origins[kept - SELECTION_N - 1],
        )
        info = (landmarks, own_rows)
        generation.update(
            landmarks_including_tips=len(landmarks),
            own_landmarks_including_first_tip=own_rows,
            shared_landmarks=int(len(kept)),
        )
        top = levels[-1]
        after = np.empty((len(landmarks), top // 8), dtype=np.uint8)
        before = np.empty_like(after)
        rng = np.random.default_rng([seed, w, 1])
        for start in range(0, top, BATCH):
            end = min(start + BATCH, top)
            events = np.vstack(
                [wave.sample(CHUNK, rng) for _ in range((end - start) // CHUNK)]
            )
            generation["sampled_events"] += len(events)
            future, past = wave.packed_relations(landmarks, events, cp)
            after[:, start // 8 : end // 8], before[:, start // 8 : end // 8] = (
                future,
                past,
            )
            generation["classified_events"] = end
            if end in levels:
                with (folder / f"{prefix}_win{w}_n{end}.npz").open("xb") as stream:
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
                    replicate=prefix,
                    window=w,
                    events=end,
                    target=top,
                )
    except (ValueError, ArithmeticError, MemoryError) as error:
        generation["failure"] = f"{type(error).__name__}: {error}"
    cp.get_default_memory_pool().free_all_blocks()
    generation.update(
        margins_resolved_in_decimal=wave.resolved_margins - resolved,
        wall_seconds=time.perf_counter() - wall,
        process_cpu_seconds=time.process_time() - cpu,
        peak_working_set_bytes=peak_working_set(),
    )
    write_json(folder / f"{prefix}_win{w}_generation.json", generation)
    return generation, info


def relation_selfcheck(waves):
    """GPU relation bits against the direct CPU path for three windows."""
    import cupy as cp

    report = {}
    for shift in (0.0, 0.5, -1.0):
        wave = waves[shift]
        rng = np.random.default_rng(20261126)
        landmarks = np.vstack((TIPS[:1], wave.sample(62, rng), TIPS[1:]))
        events = wave.sample(4096, rng)
        packed = wave.packed_relations(landmarks, events, cp)
        direct = wave.relations(landmarks[:, None], events[None])
        assert all(
            np.array_equal(np.unpackbits(p, axis=1).astype(bool), d)
            for p, d in zip(packed, direct, strict=True)
        )
        report[f"{shift:+.1f}"] = {
            "pairs": int(direct[0].size),
            "related": int(direct[0].sum() + direct[1].sum()),
        }
    cp.get_default_memory_pool().free_all_blocks()
    return report


def chart_move_selfcheck(waves):
    """The same events in neighbouring windows' charts: equal order, and equal
    exact squared intervals to 1e-10 of the largest."""
    from intrinsic_weyl_second_continuum import exact_squared_interval, oriented_pairs

    rng = np.random.default_rng(2026112603)
    report = {}
    for outer, inner in ((0.5, 0.0), (-0.5, 0.0), (1.0, 0.5), (-1.0, -0.5)):
        events = waves[outer].sample(400, rng)
        moved = events + np.array([outer - inner, inner - outer, 0.0, 0.0])
        both = np.flatnonzero((moved[:, 0] > 0) & (moved[:, 0] < 1))[:160]
        first, second = events[both], moved[both]
        equal = bool(
            np.array_equal(waves[outer].order(first), waves[inner].order(second))
        )
        pairs = oriented_pairs(first)
        mine = exact_squared_interval(waves[outer], first, pairs)
        other = exact_squared_interval(waves[inner], second, pairs)
        difference = float(np.max(abs(mine - other)) / np.max(abs(mine)))
        report[f"{outer:+.1f}->{inner:+.1f}"] = {
            "events": len(both),
            "orders_equal": equal,
            "squared_interval_relative_difference": difference,
        }
        assert equal and difference < 1e-10
    return report


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------


def lf_sha(path):
    """Source digest under the repository's LF storage-boundary rule."""
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def free_bytes(path):
    import shutil

    return shutil.disk_usage(path).free


def bootstrap_norm(errors, truth, rng, draws=BOOTSTRAP):
    """Full-tidal norm of the mean error vector and its percentile interval."""
    from intrinsic_weyl_joint import TIDAL_ROOT

    errors = np.asarray(errors)
    scale = float(np.linalg.norm(TIDAL_ROOT @ truth))
    index = rng.integers(0, len(errors), (draws, len(errors)))
    norms = np.linalg.norm(errors[index].mean(axis=1) @ TIDAL_ROOT.T, axis=1) / scale
    return {
        "relative_norm": float(
            np.linalg.norm(TIDAL_ROOT @ errors.mean(axis=0)) / scale
        ),
        "bootstrap_95": np.quantile(norms, [0.025, 0.975]).tolist(),
        "draws": draws,
    }


def describe(rows):
    """Descriptive statistics of one level; no verdict."""
    from intrinsic_weyl_bombelli_benchmark import benchmark_amplitude
    from intrinsic_weyl_window_study import unit_truth

    complete = [row for row in rows if "failure" not in row]
    result = {"complete": len(complete), "readings": {}}
    for label in CORRECTIONS:
        result["readings"][label] = {
            part: t_interval(
                [
                    (
                        row["readings"][label]
                        if part == "full"
                        else row["readings"][label][part]
                    )["value"]
                    for row in complete
                ]
            )
            for part in ("full", "electric", "magnetic")
        }
        result["readings"][label]["nominal_sd_median"] = (
            float(np.median([row["readings"][label]["nominal_sd"] for row in complete]))
            if complete
            else None
        )
    if len(complete) >= 2:
        truth = float(benchmark_amplitude()) * unit_truth()
        rng = np.random.default_rng(BOOTSTRAP_SEED)
        result["mean_error_vector"] = {
            label: bootstrap_norm(
                [
                    np.asarray(row["readings"][label]["generator_frame_estimate"])
                    - truth
                    for row in complete
                ],
                truth,
                rng,
            )
            for label in ("primary", "single")
        }
    steps = [row["links"] for row in complete if "links" in row]
    if steps:
        result["links"] = {
            shift: {
                "shared_min_median_max": [
                    int(np.min([s[shift]["shared_landmarks"] for s in steps])),
                    float(np.median([s[shift]["shared_landmarks"] for s in steps])),
                    int(np.max([s[shift]["shared_landmarks"] for s in steps])),
                ],
                "improper_to_centre": int(
                    sum(s[shift]["parity_to_centre"] < 0 for s in steps)
                ),
                "residual_median": float(
                    np.median([s[shift]["procrustes_relative_residual"] for s in steps])
                ),
            }
            for shift in steps[0]
        }
    oracle = [row["oracle_links"] for row in complete if "oracle_links" in row]
    if oracle:
        result["oracle_links"] = {
            "parity_disagreements": int(
                sum(
                    not link["parity_agrees_with_oracle"]
                    for o in oracle
                    for link in o["links"].values()
                )
            ),
            "angle_error_deg_max": max(
                (
                    link["angle_error_deg"]
                    for o in oracle
                    for link in o["links"].values()
                    if link["angle_error_deg"] is not None
                ),
                default=None,
            ),
            "primary_reading_with_oracle_links": t_interval(
                [o["primary_reading_with_oracle_links"] for o in oracle]
            ),
        }
    return result


def replicate_row(prefix, replicate, seed, n, generations, outcome):
    """One replicate at one level; a failure of any window fails the replicate."""
    row = {
        "replicate": replicate,
        "prefix": prefix,
        "seed": seed,
        "measurement_n": n,
        "generation": generations,
        "windows": outcome.get("windows", {}),
    }
    if "failure" in outcome:
        row["failure"] = outcome["failure"]
        return row
    row["links"] = outcome["links"]
    evaluation = outcome["evaluation"]
    row["readings"] = evaluation["readings"]
    row["centre_alignment_parity"] = evaluation["centre_alignment_parity"]
    if "oracle_links" in evaluation:
        row["oracle_links"] = evaluation["oracle_links"]
    else:
        row["oracle_links_failure"] = evaluation["oracle_links_failure"]
    row["combine_wall_seconds"] = outcome["wall_seconds"]
    return row


def run(destination, engineering=False, parallel=3):
    import cupy as cp
    import scipy
    from intrinsic_weyl_bombelli_benchmark import benchmark_amplitude
    from intrinsic_weyl_gpu_counts import selfcheck as gpu_selfcheck
    from intrinsic_weyl_joint import selfcheck_jet
    from profile_plane_wave import selfcheck as wave_selfcheck

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
            "intrinsic_weyl_window_study.py",
            "profile_plane_wave.py",
            "intrinsic_weyl_bombelli_benchmark.py",
            "replay_intrinsic_weyl_threads.py",
            "intrinsic_weyl_joint.py",
            "intrinsic_weyl.py",
            "intrinsic_weyl_square.py",
            "intrinsic_weyl_gpu_counts.py",
            "intrinsic_weyl_high_count.py",
            "intrinsic_weyl_joint_audit.py",
            "intrinsic_weyl_second_continuum.py",
        )
    ]
    if not engineering:
        sources.append(ROOT / PROTOCOL)
    source_hashes = {
        str(p.relative_to(ROOT)).replace("\\", "/"): lf_sha(p) for p in sources
    }
    if engineering:
        seeds, levels = ENGINEERING_SEEDS, ENGINEERING_LEVELS
        kept = set(range(len(seeds)))
    else:
        seeds = tuple(SEED + r for r in range(REPLICATES))
        levels, kept = LEVELS, set(KEPT_ARCHIVES)
    amplitude = benchmark_amplitude()
    waves = build_waves(amplitude)
    checked = {
        "generator": wave_selfcheck(),
        "gpu_relation_bits": relation_selfcheck(waves),
        "chart_move": chart_move_selfcheck(waves),
        "jet_inverse": selfcheck_jet(),
        "gpu_second_moment": gpu_selfcheck(),
        "bounded_second_moment": moment_selfcheck(),
        "window_selection": selection_selfcheck(),
        "links": link_selfcheck(),
        "verdict": selfcheck(),
    }
    if free_bytes(destination.parent) < 80 * 2**30:
        raise RuntimeError("At least 80 GiB free storage is required")
    destination.mkdir(parents=True)
    started = time.time()
    lock = threading.Lock()

    def status(**fields):
        with lock:
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
    worker_env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", CUPY_TF32="0")
    manifest = {
        "engineering_only": bool(engineering),
        "code_version": stamp,
        "worktree_clean_at_start": not dirty,
        "source_sha256_lf": source_hashes,
        "selfcheck": checked,
        "benchmark_amplitude": str(amplitude),
        "windows": list(WINDOWS),
        "generation_order": list(GENERATION_ORDER),
        "replicates": [
            {"replicate": r, "prefix": f"rep{r:02d}", "seed": seed}
            for r, seed in enumerate(seeds)
        ],
        "rng_streams": "selection default_rng([seed, w, 0]); measurement [seed, w, 1]",
        "measurement_levels": levels,
        "selection_n": SELECTION_N,
        "sampling_chunk": CHUNK,
        "classification_batch": BATCH,
        "archives_kept_for_replicates": sorted(kept),
        "tolerance": TOLERANCE,
        "coverage_two_sided": COVERAGE,
        "minimum_complete": MINIMUM_COMPLETE if not engineering else None,
        "minimum_shared_per_link": MINIMUM_SHARED,
        "primary_correction": list(PRIMARY),
        "parallel_fit_workers": parallel,
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
            "Engineering pilot; no test verdict"
            if engineering
            else "Accuracy gate on the amplitude along the true tensor for one plane "
            "wave, amplitude and interval size, with the single window as control; "
            "no full-tensor, generality or single-replicate claim"
        ),
    }
    write_json(destination / "manifest.json", manifest)

    def invoke(arguments, output):
        with output.with_suffix(".log").open("xb") as log:
            child = subprocess.Popen(
                [sys.executable, str(Path(__file__).resolve()), *map(str, arguments)],
                env=worker_env,
                stdout=log,
                stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW,
            )
            code = child.wait()
        if code:
            return {
                "failure": f"Worker exit {code}; see {output.with_suffix('.log').name}"
            }
        return json.loads(output.read_text(encoding="utf-8"))

    slots = threading.BoundedSemaphore(OUTSTANDING * len(levels))

    def fit_task(prefix, w, n, keep):
        archive = destination / f"{prefix}_win{w}_n{n}.npz"
        output = archive.with_suffix(".json")
        arguments = ["--worker", archive, "--output", output]
        if w == CENTRE:
            arguments += ["--clocks", destination / f"{prefix}_win{w}_n{n}_clocks.npz"]
        try:
            record = invoke(arguments, output)
            record["archive_bytes"] = archive.stat().st_size
            if not keep:
                archive.unlink()
                record["archive_deleted_after_fit"] = True
            return record
        finally:
            slots.release()

    def combine_task(prefix, replicate, seed, n, generations, futures, keep):
        records = {str(w): futures[w].result() for w in sorted(futures)}
        missing = [w for w in range(len(WINDOWS)) if w not in futures]
        failed = {w: r["failure"] for w, r in records.items() if "failure" in r}
        if missing or failed:
            reasons = {
                str(w): generations.get(w, {}).get("failure")
                or generations.get(w, {}).get("skipped")
                for w in missing
            }
            outcome = {
                "windows": records,
                "failure": f"Window failures: missing {reasons}, fits {failed}",
            }
        else:
            output = destination / f"{prefix}_n{n}_combined.json"
            outcome = invoke(
                [
                    "--combine",
                    destination,
                    "--prefix",
                    prefix,
                    "--level",
                    n,
                    "--output",
                    output,
                ],
                output,
            )
            outcome["windows"] = records
            clocks = destination / f"{prefix}_win{CENTRE}_n{n}_clocks.npz"
            if not keep and clocks.exists():
                outcome["centre_clocks_sha256"] = sha(clocks)
                clocks.unlink()
        row = replicate_row(prefix, replicate, seed, n, generations, outcome)
        write_json(destination / f"{prefix}_n{n}_result.json", row)
        primary = row.get("readings", {}).get("primary", {})
        single = row.get("readings", {}).get("single", {})
        print(
            json.dumps(
                {
                    "replicate": replicate,
                    "n": n,
                    "primary": primary.get("value"),
                    "primary_nominal_sd": primary.get("nominal_sd"),
                    "single": single.get("value"),
                    "failure": row.get("failure"),
                }
            ),
            flush=True,
        )
        return row

    fit_pool = ThreadPoolExecutor(max_workers=parallel)
    combine_pool = ThreadPoolExecutor(max_workers=2)
    pending = []
    try:
        for replicate, seed in enumerate(seeds):
            prefix = f"rep{replicate:02d}"
            keep = replicate in kept
            landmarks, generations = {}, {}
            futures = {n: {} for n in levels}
            stopped = None
            for w in GENERATION_ORDER:
                if stopped is not None:
                    generations[w] = {"skipped": stopped}
                    continue
                for _ in levels:
                    slots.acquire()
                neighbour = None
                if w in INNER:
                    points, own_rows = landmarks[INNER[w]]
                    neighbour = (points, own_rows, WINDOWS[INNER[w]])
                generation, info = generate_window(
                    destination,
                    prefix,
                    w,
                    seed,
                    waves[WINDOWS[w]],
                    levels,
                    status,
                    neighbour,
                )
                generations[w] = generation
                if info is not None:
                    landmarks[w] = info
                produced = 0
                for n in levels:
                    if (destination / f"{prefix}_win{w}_n{n}.npz").exists():
                        futures[n][w] = fit_pool.submit(fit_task, prefix, w, n, keep)
                        produced += 1
                for _ in range(len(levels) - produced):
                    slots.release()
                if "failure" in generation:
                    stopped = f"window {w} failed in generation"
            for n in levels:
                pending.append(
                    combine_pool.submit(
                        combine_task,
                        prefix,
                        replicate,
                        seed,
                        n,
                        dict(generations),
                        futures[n],
                        keep,
                    )
                )
        status(state="RUNNING", phase="draining fits")
        rows = [future.result() for future in pending]
        summary = {
            **manifest,
            "results": rows,
            "elapsed_seconds": time.time() - started,
        }
        summary["levels"] = {}
        for n in levels:
            level_rows = [row for row in rows if row["measurement_n"] == n]
            entry = {"descriptive": describe(level_rows)}
            if not engineering:
                entry["verdict"] = judge(level_rows)
            summary["levels"][str(n)] = entry
        if not engineering:
            summary["primary_claim"] = summary["levels"][str(LEVELS[0])]["verdict"][
                "primary_claim"
            ]
        for name, digest in source_hashes.items():
            if lf_sha(ROOT / name) != digest:
                raise RuntimeError(f"A bound source changed during execution: {name}")
        write_json(destination / "summary.json", summary)
        status(
            state="COMPLETED",
            replicates=len(seeds),
            primary_claim=summary.get("primary_claim"),
        )
    except BaseException as error:
        status(state="STOPPED", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        fit_pool.shutdown(wait=True, cancel_futures=True)
        combine_pool.shutdown(wait=True, cancel_futures=True)


def compare(first, second):
    """Engineering replays must agree: exact order bits, 1e-8 on every fit field."""
    report = {"first": str(first), "second": str(second), "items": {}}
    for replicate in range(len(ENGINEERING_SEEDS)):
        prefix = f"rep{replicate:02d}"
        for w in range(len(WINDOWS)):
            name = f"{prefix}_win{w}_selection.npz"
            with (
                np.load(first / name, allow_pickle=False) as a,
                np.load(second / name, allow_pickle=False) as b,
            ):
                report["items"][name] = bool(
                    set(a.files) == set(b.files)
                    and all(np.array_equal(a[k], b[k]) for k in a.files)
                )
            for n in ENGINEERING_LEVELS:
                name = f"{prefix}_win{w}_n{n}"
                with (
                    np.load(first / f"{name}.npz", allow_pickle=False) as a,
                    np.load(second / f"{name}.npz", allow_pickle=False) as b,
                ):
                    bits = set(a.files) == set(b.files) and all(
                        np.array_equal(a[k], b[k]) for k in a.files
                    )
                old, new = (
                    json.loads((f / f"{name}_pure.json").read_text())["fit"]
                    for f in (first, second)
                )
                for record in (old, new):
                    record.pop("order_archive_sha256")
                report["items"][name] = {
                    "order_archive_bit_exact": bool(bits),
                    "fit_differences": differences(new, old),
                }
        for n in ENGINEERING_LEVELS:
            name = f"{prefix}_n{n}_combined"
            old, new = (
                json.loads((f / f"{name}_pure.json").read_text())["correction"]
                for f in (first, second)
            )
            old_eval, new_eval = (
                json.loads((f / f"{name}.json").read_text()).get("evaluation")
                for f in (first, second)
            )
            report["items"][name] = {
                "correction_differences": differences(new, old),
                "evaluation_differences": differences(new_eval, old_eval),
            }
    report["reproduced"] = all(
        item
        if isinstance(item, bool)
        else all(not v if isinstance(v, dict) else v for v in item.values())
        for item in report["items"].values()
    )
    write_json(second / "replay_comparison.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--engineering", action="store_true")
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--clocks", type=Path)
    parser.add_argument("--combine", type=Path)
    parser.add_argument("--prefix")
    parser.add_argument("--level", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--parallel", type=int, default=3)
    parser.add_argument("--compare", type=Path, nargs=2)
    args = parser.parse_args()
    if args.compare:
        outcome = compare(*args.compare)
        print(json.dumps({"reproduced": outcome["reproduced"]}))
        sys.exit(0 if outcome["reproduced"] else 1)
    elif args.check:
        print(
            json.dumps(
                {
                    "verdict": selfcheck(),
                    "window_selection": selection_selfcheck(),
                    "links": link_selfcheck(),
                    "bounded_second_moment": moment_selfcheck(),
                },
                indent=1,
            )
        )
    elif args.worker:
        if args.output is None:
            parser.error("output is required for a worker")
        window_worker(args.worker, args.output, args.clocks)
    elif args.combine:
        if args.prefix is None or args.level is None or args.output is None:
            parser.error("prefix, level and output are required to combine")
        combine_worker(args.combine, args.prefix, args.level, args.output)
    else:
        if args.destination is None:
            parser.error("destination is required")
        run(args.destination, engineering=args.engineering, parallel=args.parallel)
