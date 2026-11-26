"""Independent post-completion audit of the five-window link test run.

A reviewed copy of the frozen checker ``check_intrinsic_weyl_window_test.py``
(commit 1ca7af2) for ``intrinsic_weyl_window_link_test.py``: windows at
offsets 0, +-1/3, +-2/3, the five-point stencil at h = 1/3, and the
falsification readings (centre term alone; identity, handedness-flipped and
Haar-random side links; each side window's linked reading) with their gates.

Nothing here calls the runner's gates, links, correction or evaluation code.
From the saved pure fits and pure corrections it rebuilds the frame links (a
polar decomposition, not the runner's SVD path), the Weyl maps of O(3) through
the 3x3 electric and magnetic tensors (B is a pseudotensor), the structure of
the stencil weights, every corrected estimate, every reading and counterfactual
reading (the Haar draws from their seed), the intervals and the verdicts.
Every window's landmarks and shared rows are re-derived from its saved
selection order with the checker's own count rule; only the nested selection
of the window's own events is the frozen inverse's.

For the replicates whose archives are kept it also re-derives the clock mask,
replays the first 4096 measurement events of every window on the CPU, checks a
sample of order bits with a general-purpose ODE integrator that shares no code
with the generator, checks that excluded candidates lie outside the window,
checks that shared landmarks are the same events in both charts, and recomputes
the minimum-norm operator. Two windows are refitted in fresh processes. The
truth vectors are reused from the run after their invariants and orientation
are checked. Every `posthoc` entry decides nothing.

It also restates every interval the record quotes (the corrected, single,
centre-term, identity, flipped and side-share readings, the Haar summary and
the four link-sign intervals), recomputes each side link's parity and angle
against the oracle relative frames by polar decomposition of the truth
landmarks (G5), checks the correction keys and shifts, the run's constants
and seeds against the protocol, the run's source digests against the stamped
commit, and that the bound sources at this checkout are the stamped ones (so
that the refits run the stamped code).

Shared with the runner, so not independent: the frozen inverse's nested
selection and fit (refits), the generator's window geometry (truth side), and
the alias map, whitening and minimum-norm operator of
``intrinsic_weyl_window_study``, which recompute the operator for the two
replicates whose clocks are kept; for the others only the five-point
structure of the weights is checked.
"""

import argparse
import hashlib
import itertools
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from check_intrinsic_weyl_weak_accuracy import read, sha
from intrinsic_weyl import WEYL
from replay_intrinsic_weyl_threads import THREAD_VARIABLES, differences, write_json

ROOT = Path(__file__).resolve().parents[2]
RUNNER = Path(__file__).with_name("intrinsic_weyl_window_link_test.py")
SPACING = 1 / 3
WINDOWS = tuple(k * SPACING for k in (-2, -1, 0, 1, 2))
INNER = {0: 1, 1: 2, 3: 2, 4: 3}
CENTRE, SELECTION_N, H0 = 2, 4096, 0.01
SIDES = (0, 1, 3, 4)
TOLERANCE, COVERAGE, MINIMUM_SHARED = 0.05, 0.95, 8
HAAR_SEED, HAAR_DRAWS = 2026113050, 64
SIDE_SHARE_BAND = (-0.22, -0.17)
# The protocol's run constants (Sects. 6, 8 and 10).
PROTOCOL = "docs/paper/paper_a/intrinsic_weyl/window_link_test_protocol.md"
REPLICATES, SEED, MINIMUM_COMPLETE = 96, 2026112800, 87
LEVELS, KEPT_ARCHIVES = (33554432,), (0, 48)
ENGINEERING_SEEDS, ENGINEERING_LEVELS = (2026113000, 2026113001), (4194304, 16777216)
CORRECTIONS = {
    "single": ["single_window", 0.0, "none"],
    "primary": ["five_point_second", SPACING, "minimum_norm_time_time"],
    "five_point_plane_wave": ["five_point_second", SPACING, "plane_wave_estimated_k"],
    "three_point_plane_wave": ["three_point_first", SPACING, "plane_wave_estimated_k"],
    "three_point_minimum_norm": [
        "three_point_first",
        SPACING,
        "minimum_norm_time_time",
    ],
}
SHIFTS = {
    "single_window": [WINDOWS[CENTRE]],
    "three_point_first": [WINDOWS[1], WINDOWS[2], WINDOWS[3]],
    "five_point_second": list(WINDOWS),
}
EPSILON = np.zeros((3, 3, 3))
for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
    EPSILON[i, j, k], EPSILON[i, k, j] = 1, -1
# E_ij = C_0i0j and B_ij = eps_ikl C_kl0j / 2, taken from the full tensor basis.
ELECTRIC = WEYL[:, 0, 1:, 0, 1:]
MAGNETIC = 0.5 * np.einsum("ikl,nklj->nij", EPSILON, WEYL[:, 1:, 1:, 0, 1:])
GRAM2 = sum(np.einsum("nij,mij->nm", part, part) for part in (ELECTRIC, MAGNETIC))
# Rows: the nine E and nine B components; columns: the ten coefficients.
BASIS = np.vstack((ELECTRIC.reshape(10, 9).T, MAGNETIC.reshape(10, 9).T))


def tensors(c):
    """Electric and magnetic 3x3 parts of a ten-coefficient vector."""
    return tuple(np.einsum("n,nij->ij", c, part) for part in (ELECTRIC, MAGNETIC))


def coefficients(electric, magnetic):
    target = np.r_[electric.ravel(), magnetic.ravel()]
    c, *_ = np.linalg.lstsq(BASIS, target, rcond=None)
    if np.max(abs(BASIS @ c - target)) > 1e-12 * max(1.0, np.max(abs(target))):
        raise ArithmeticError("Tensor pair outside the Weyl basis")
    return c


def weyl_map(rotation, c):
    """Coefficients in coordinates x' = x @ R; B changes sign with det R."""
    e, b = tensors(c)
    parity = np.linalg.det(rotation)
    return coefficients(rotation.T @ e @ rotation, parity * rotation.T @ b @ rotation)


def weyl_matrix(rotation):
    return np.column_stack([weyl_map(rotation, unit) for unit in np.eye(10)])


def polar(target, source):
    """Orthogonal Q minimizing |target - source @ Q| by polar decomposition."""
    m = source.T @ target
    values, vectors = np.linalg.eigh(m.T @ m)
    if values[0] <= 1e-12 * values[-1]:
        raise ArithmeticError("Degenerate shared landmark cloud")
    return m @ vectors @ np.diag(values**-0.5) @ vectors.T


def norm(c):
    return math.sqrt(float(c @ GRAM2 @ c))


def close(left, right, tolerance=1e-9):
    left, right = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
    scale = max(1.0, float(np.max(abs(left))), float(np.max(abs(right))))
    return bool(np.max(abs(left - right)) <= tolerance * scale)


def interval(values):
    from scipy.stats import t as student

    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return None
    half = student.ppf(0.5 + COVERAGE / 2, len(values) - 1) * values.std(ddof=1)
    half /= math.sqrt(len(values))
    return [float(values.mean() - half), float(values.mean() + half)]


def label(shift):
    """Key of a window in the run's records."""
    return f"{shift:+.4f}"


def same(ours, theirs, tolerance=1e-9):
    """Equal to the tolerance, where None (fewer than two values) equals None."""
    if ours is None or theirs is None:
        return ours is None and theirs is None
    return close(ours, theirs, tolerance)


def bounds_of(reading):
    """The interval of one of the run's t-interval records, or None."""
    return None if reading is None else reading["interval"]


def gates(primary, single, minimum, counterfactuals, oracle):
    """The gates of the protocol, restated (Sect. 5), with every interval the
    record quotes. ``primary``, ``single`` and ``counterfactuals`` run over the
    completed replicates in one order; ``oracle`` holds for each of them the
    side links' parity agreements with the oracle, or None if unavailable."""
    corrected, control = interval(primary), interval(single)
    g1 = g2 = f1 = f2 = g4 = None
    if corrected:
        if corrected[0] > 1 - TOLERANCE and corrected[1] < 1 + TOLERANCE:
            g1 = "ACCURATE"
        elif corrected[1] < 1 - TOLERANCE or corrected[0] > 1 + TOLERANCE:
            g1 = "INACCURATE"
        else:
            g1 = "UNRESOLVED"
    if control:
        g2 = "BIASED_LOW" if control[1] < 1 - TOLERANCE else "NOT_BIASED"

    def above(bounds):
        return bounds is not None and bounds[0] > 1 + TOLERANCE

    centre = interval([c["centre_term_only"] for c in counterfactuals])
    identity = interval([c["identity_links"] for c in counterfactuals])
    flipped = interval([c["flipped_links"] for c in counterfactuals])
    haar = [
        interval([c["haar_links"][d] for c in counterfactuals])
        for d in range(HAAR_DRAWS)
    ]
    if centre:
        f1 = "CENTRE_TERM_MISSES" if above(centre) else "NOT_SHOWN"
        wrong = [identity, flipped, *haar]
        f2 = "WRONG_LINKS_MISS" if all(above(v) for v in wrong) else "NOT_SHOWN"
    signs, sign_intervals = {}, {}
    for w in SIDES:
        key = label(WINDOWS[w])
        bounds = interval([c["linked_side_readings"][key] for c in counterfactuals])
        sign_intervals[key] = bounds
        if bounds is None:
            signs[key] = None
        elif abs(w - CENTRE) == 1:
            signs[key] = "LINKED" if bounds[0] > 0 else "NOT_LINKED"
        else:
            signs[key] = "LINKED" if bounds[1] < 0 else "NOT_LINKED"
    share = interval(
        [
            p - c["centre_term_only"]
            for p, c in zip(primary, counterfactuals, strict=True)
        ]
    )
    if share:
        low, high = SIDE_SHARE_BAND
        if share[0] > low and share[1] < high:
            g4 = "AS_PREDICTED"
        elif share[1] < low or share[0] > high:
            g4 = "NOT_AS_PREDICTED"
        else:
            g4 = "UNRESOLVED"
    available = [flags for flags in oracle if flags is not None]
    disagreeing = int(sum(not flag for flags in available for flag in flags))
    missing = len(oracle) - len(available)
    if disagreeing:
        g5 = "DISAGREE"
    elif missing or not oracle:
        g5 = "UNAVAILABLE"
    else:
        g5 = "ALL_AGREE"
    if len(primary) < minimum:
        claim = "INCOMPLETE"
    elif (g1, g2, f1, f2) == (
        "ACCURATE",
        "BIASED_LOW",
        "CENTRE_TERM_MISSES",
        "WRONG_LINKS_MISS",
    ) and all(signs[label(WINDOWS[w])] == "LINKED" for w in SIDES):
        claim = "SUPPORTED"
    else:
        claim = "NOT_SUPPORTED"
    middles = [(h[0] + h[1]) / 2 for h in haar if h]
    return {
        "complete": len(primary),
        "minimum_complete": minimum,
        "interval": corrected,
        "control_interval": control,
        "centre_term_only_interval": centre,
        "identity_links_interval": identity,
        "flipped_links_interval": flipped,
        "haar_draws_missing_band": int(sum(above(h) for h in haar)),
        "haar_lowest_lower_bound": min((h[0] for h in haar if h), default=None),
        "haar_mean_range": [min(middles), max(middles)] if middles else [None, None],
        "link_sign_intervals": sign_intervals,
        "side_share_interval": share,
        "g1": g1,
        "g2": g2,
        "f1": f1,
        "f2": f2,
        "g3": signs,
        "g4": g4,
        "g5": g5,
        "g5_links_compared": int(sum(len(flags) for flags in available)),
        "g5_links_disagreeing": disagreeing,
        "g5_replicates_without_oracle": missing,
        "claim": claim,
    }


def matches_verdict(restated, verdict):
    """Item-by-item comparison of the restated gates with the run's verdict."""
    signs = verdict["g3_link_signs"]
    share = verdict["g4_side_share_secondary"]
    handedness = verdict["g5_oracle_handedness_secondary"]
    haar = verdict["haar_links"]
    checks = {
        "complete": restated["complete"] == verdict["complete"],
        "minimum_complete": restated["minimum_complete"] == verdict["minimum_complete"],
        "g1": restated["g1"] == verdict["g1_accuracy"],
        "g2": restated["g2"] == verdict["g2_control"],
        "f1": restated["f1"] == verdict["f1_centre_term"],
        "f2": restated["f2"] == verdict["f2_wrong_links"],
        "g3": restated["g3"] == {k: v["verdict"] for k, v in signs.items()},
        "g4": restated["g4"] == share["verdict"]
        and list(SIDE_SHARE_BAND) == share["band"],
        "g5": restated["g5"] == handedness["verdict"]
        and restated["g5_links_compared"] == handedness["links_compared"]
        and restated["g5_links_disagreeing"] == handedness["links_disagreeing"]
        and restated["g5_replicates_without_oracle"]
        == handedness["replicates_without_oracle_links"],
        "claim": restated["claim"] == verdict["primary_claim"],
        "corrected_interval": same(
            restated["interval"], bounds_of(verdict["corrected_reading"])
        ),
        "single_interval": same(
            restated["control_interval"], bounds_of(verdict["single_window_reading"])
        ),
        "centre_term_interval": same(
            restated["centre_term_only_interval"],
            bounds_of(verdict["centre_term_only_reading"]),
        ),
        "identity_links_interval": same(
            restated["identity_links_interval"],
            bounds_of(verdict["identity_links_reading"]),
        ),
        "flipped_links_interval": same(
            restated["flipped_links_interval"],
            bounds_of(verdict["flipped_links_reading"]),
        ),
        "haar_summary": haar["draws"] == HAAR_DRAWS
        and restated["haar_draws_missing_band"] == haar["draws_missing_band"]
        and same(restated["haar_lowest_lower_bound"], haar["lowest_lower_bound"])
        and all(
            same(a, b)
            for a, b in zip(
                restated["haar_mean_range"], haar["mean_range"], strict=True
            )
        ),
        "link_sign_intervals": set(restated["link_sign_intervals"]) == set(signs)
        and all(
            same(bounds, bounds_of(signs[key]["linked_reading"]))
            for key, bounds in restated["link_sign_intervals"].items()
        ),
        "side_share_interval": same(
            restated["side_share_interval"], bounds_of(share["reading"])
        ),
    }
    return all(checks.values()), checks


def stencil(weights, shifts):
    """B from the weights of c0 - B D2 + B^2 D4 at h = 1/3, and the residual of
    the structure. With d2 = (-1, 16, -30, 16, -1)/12 and d4 = (1, -4, 6, -4, 1),
    W(+-h) = -(4/3) B / h^2 - 4 B^2 / h^4 and W(+-2h) = B / (12 h^2) + B^2 / h^4,
    so B = -h^2 (W(h) + 4 W(2h)); at h = 1/3: W(+-2h) = 0.75 B + 81 B^2,
    W(+-h) = -12 B - 324 B^2, W(0) = I + 22.5 B + 486 B^2."""
    w = {
        int(round(s / SPACING)): np.asarray(m)
        for s, m in zip(shifts, weights, strict=True)
    }
    b = -(SPACING**2) * (w[1] + 4 * w[2])
    b2 = b @ b
    expected = {
        -2: 0.75 * b + 81 * b2,
        -1: -12 * b - 324 * b2,
        0: np.eye(10) + 22.5 * b + 486 * b2,
        1: -12 * b - 324 * b2,
        2: 0.75 * b + 81 * b2,
    }
    scale = max(1.0, float(np.max(abs(b))))
    return b, max(float(np.max(abs(w[k] - expected[k]))) for k in expected) / scale


def haar_o3(rng):
    """The protocol's Haar draw on O(3): QR of a standard normal matrix with the
    sign fix (the definition of the draw, restated here)."""
    q, r = np.linalg.qr(rng.standard_normal((3, 3)))
    return q * np.sign(np.diag(r))


def audit_replicate(folder, row):
    """Links, maps, stencil, estimates, readings and the blind boundary."""
    prefix, n = row["prefix"], row["measurement_n"]
    pure = {w: read(folder / f"{prefix}_win{w}_n{n}_pure.json") for w in range(5)}
    fits = {WINDOWS[w]: pure[w]["fit"] for w in range(5)}
    combined = read(folder / f"{prefix}_n{n}_combined_pure.json")
    correction = combined["correction"]
    evaluation = read(folder / f"{prefix}_n{n}_combined.json")["evaluation"]
    report = {"replicate": row["replicate"], "measurement_n": n}
    report["blind_boundary_kept"] = bool(
        all(p["oracle_modules_before_return"] == [] for p in pure.values())
        and all(
            [Path(x) for x in p["opened_archives"]]
            == [(folder / f"{prefix}_win{w}_n{n}.npz").resolve()]
            for w, p in pure.items()
        )
        and combined["oracle_modules_before_return"] == []
        and [Path(x) for x in combined["opened_inputs"]]
        == [(folder / f"{prefix}_win{CENTRE}_n{n}_clocks.npz").resolve()]
        and all(
            combined["pure_fit_sha256"][f"{prefix}_win{w}_n{n}_pure.json"]
            == sha(folder / f"{prefix}_win{w}_n{n}_pure.json")
            for w in range(5)
        )
    )
    # Links by polar decomposition.
    inferred = {
        s: np.asarray(f["coordinates_inferred_from_counts"]) for s, f in fits.items()
    }
    centre = WINDOWS[CENTRE]
    sides = [WINDOWS[w] for w in SIDES]
    rotations = {centre: np.eye(3)}
    shared_ok = True
    for w in (1, 3, 0, 4):
        s, neighbour = WINDOWS[w], WINDOWS[INNER[w]]
        rows = np.asarray(fits[s]["shared_rows"], dtype=int)
        sources = np.asarray(fits[s]["shared_source_rows"], dtype=int)
        shared_ok &= len(rows) >= MINIMUM_SHARED
        step = polar(inferred[neighbour][sources, 1:], inferred[s][rows, 1:])
        rotations[s] = step @ rotations[neighbour]
    recorded = {
        s: np.asarray(correction["links"]["rotations"][label(s)]) for s in WINDOWS
    }
    report["links_reproduced"] = all(
        close(rotations[s], recorded[s], 1e-9) for s in WINDOWS
    )
    report["every_link_has_eight_shared"] = bool(shared_ok)
    report["link_parities_to_centre"] = {
        label(s): int(round(np.linalg.det(rotations[s]))) for s in sides
    }
    # Estimates through the checker's own Weyl maps; stencil structure.
    own = {s: np.asarray(fits[s]["coefficients"]) for s in WINDOWS}
    mapped = {s: weyl_matrix(rotations[s]) @ own[s] for s in WINDOWS}
    estimates_ok, structure = True, {}
    for name, item in correction["corrections"].items():
        estimate = sum(
            np.asarray(weight) @ mapped[s]
            for s, weight in zip(item["shifts"], item["weights"], strict=True)
        )
        estimates_ok &= close(estimate, item["estimate"], 1e-9)
        if item["key"][0] == "five_point_second":
            b, residual = stencil(item["weights"], item["shifts"])
            structure[name] = residual
            if name == "primary":
                report["primary_operator"] = b
    report["estimates_reproduced"] = bool(estimates_ok)
    report["five_point_structure_residual"] = structure
    report["corrections_as_protocol"] = bool(
        set(correction["corrections"]) == set(CORRECTIONS)
        and all(
            item["key"] == CORRECTIONS[name]
            and item["shifts"] == SHIFTS[item["key"][0]]
            for name, item in correction["corrections"].items()
        )
    )
    # Truth: invariants and orientation, then the readings.
    truth = np.asarray(evaluation["truth_coefficients_centre_frame"])
    e, b = tensors(truth)
    amplitude = math.pi**2 * H0 / 2
    target = amplitude * np.array([-1.0, 0.0, 1.0])
    report["truth_invariant_error"] = float(
        max(
            abs(np.linalg.eigvalsh(e) - target).max(),
            abs(np.linalg.eigvalsh(b) - target).max(),
            np.linalg.norm(e @ b + b @ e),
        )
    )
    readings, agree = {}, []
    for name, item in correction["corrections"].items():
        estimate = np.asarray(item["estimate"])
        value = float(truth @ GRAM2 @ estimate) / float(truth @ GRAM2 @ truth)
        readings[name] = value
        agree.append(close(value, evaluation["readings"][name]["value"], 1e-9))
        agree.append(close(value, row["readings"][name]["value"], 1e-9))
    report["readings"] = readings
    report["readings_reproduced"] = all(agree)
    # Counterfactual readings with the checker's own maps and Haar draws.
    primary = correction["corrections"]["primary"]
    weights = {
        s: np.asarray(m)
        for s, m in zip(primary["shifts"], primary["weights"], strict=True)
    }

    def along(vector):
        return float(truth @ GRAM2 @ vector) / float(truth @ GRAM2 @ truth)

    def with_side_maps(maps):
        return along(
            weights[centre] @ own[centre]
            + sum(weights[s] @ maps[s] @ own[s] for s in sides)
        )

    rng = np.random.default_rng([HAAR_SEED, int(row["replicate"])])
    haar = [
        with_side_maps({s: weyl_matrix(haar_o3(rng)) for s in sides})
        for _ in range(HAAR_DRAWS)
    ]
    ours = {
        "centre_term_only": along(weights[centre] @ own[centre]),
        "identity_links": with_side_maps({s: np.eye(10) for s in sides}),
        "flipped_links": with_side_maps({s: weyl_matrix(-rotations[s]) for s in sides}),
        "haar_links": haar,
        "linked_side_readings": {label(s): along(mapped[s]) for s in sides},
    }
    saved = evaluation["counterfactuals"]
    report["counterfactuals"] = ours
    report["counterfactuals_reproduced"] = bool(
        all(
            close(ours[k], saved[k], 1e-9) and close(ours[k], row["counterfactuals"][k])
            for k in (
                "centre_term_only",
                "identity_links",
                "flipped_links",
                "haar_links",
            )
        )
        and set(saved["linked_side_readings"]) == set(ours["linked_side_readings"])
        and set(row["counterfactuals"]["linked_side_readings"])
        == set(ours["linked_side_readings"])
        and all(
            close(ours["linked_side_readings"][k], saved["linked_side_readings"][k])
            and close(
                ours["linked_side_readings"][k],
                row["counterfactuals"]["linked_side_readings"][k],
            )
            for k in ours["linked_side_readings"]
        )
        and saved["haar_rng"] == [HAAR_SEED, int(row["replicate"])]
    )
    # Oracle relative frames (truth side), by polar decomposition of the truth
    # landmarks: each side link's parity and angle against the oracle (G5), and
    # the primary reading with oracle links.
    recorded = row.get("oracle_links")
    report["oracle_parity_flags"] = None
    try:
        frames = {
            shift: polar(
                inferred[shift][1:-1, 1:],
                normal_coordinates(folder, prefix, w)[1:-1, 1:],
            )
            for w, shift in enumerate(WINDOWS)
        }
    except ArithmeticError as error:
        report["oracle_failure"] = f"{type(error).__name__}: {error}"
        report["oracle_reproduced"] = recorded is None
        return report, truth
    oracle = {s: frames[s].T @ frames[centre] for s in WINDOWS}
    flags, angles = {}, {}
    for s in sides:
        difference = rotations[s] @ oracle[s].T
        flags[label(s)] = bool(np.linalg.det(difference) > 0)
        cosine = float(np.clip((np.trace(difference) - 1) / 2, -1, 1))
        angles[label(s)] = math.degrees(math.acos(cosine)) if flags[label(s)] else None
    with_oracle = along(
        sum(weights[s] @ weyl_matrix(oracle[s]) @ own[s] for s in WINDOWS)
    )
    report["oracle_angles_deg"] = angles
    if recorded is None:
        # The run recorded no oracle links for this replicate: G5 counts it as
        # unavailable, so the checker does too.
        report["oracle_reproduced"] = None
        return report, truth
    links = recorded["links"]
    report["oracle_parity_flags"] = [flags[label(s)] for s in sides]
    report["oracle_reproduced"] = bool(
        set(links) == set(flags)
        and all(links[k]["parity_agrees_with_oracle"] == flags[k] for k in flags)
        and all(
            (links[k]["angle_error_deg"] is None) == (angles[k] is None)
            and (
                angles[k] is None or abs(links[k]["angle_error_deg"] - angles[k]) < 1e-5
            )
            for k in flags
        )
        and close(with_oracle, recorded["primary_reading_with_oracle_links"], 1e-9)
    )
    return report, truth


def normal_coordinates(folder, prefix, w):
    """True normal coordinates of one window's landmarks (truth side)."""
    from intrinsic_weyl_bombelli_benchmark import benchmark_amplitude
    from intrinsic_weyl_window_study import build_window

    with np.load(folder / f"truth_{prefix}_win{w}.npz", allow_pickle=False) as saved:
        points = saved["points"]
    return build_window(
        "sinusoid", benchmark_amplitude(), 2.0, WINDOWS[w]
    ).normal_coordinates(points)


def truth_orientation(folder, row, truth):
    """The centre-frame truth maps back to a plus-polarized wave along z in the
    generator frame through the alignment of the inferred to the true
    coordinates (computed here by polar decomposition)."""
    from intrinsic_weyl_window_study import build_window

    prefix, n = row["prefix"], row["measurement_n"]
    with np.load(
        folder / f"truth_{prefix}_win{CENTRE}.npz", allow_pickle=False
    ) as saved:
        points = saved["points"]
    from intrinsic_weyl_bombelli_benchmark import benchmark_amplitude

    normal = build_window(
        "sinusoid", benchmark_amplitude(), 2.0, 0.0
    ).normal_coordinates(points)
    fit = read(folder / f"{prefix}_win{CENTRE}_n{n}_pure.json")["fit"]
    inferred = np.asarray(fit["coordinates_inferred_from_counts"])
    rotation = polar(inferred[1:-1, 1:], normal[1:-1, 1:])  # inferred ~ normal @ R
    generator = weyl_map(rotation.T, truth)
    e, b = tensors(generator)
    a = math.pi**2 * H0 / 2
    return float(
        max(
            abs(e[2]).max(),
            abs(e[:, 2]).max(),
            abs(e[0, 1]),
            abs(abs(e[0, 0]) - a),
            abs(e[0, 0] + e[1, 1]),
            abs(b[2]).max(),
            abs(b[:, 2]).max(),
            abs(b[0, 0]),
            abs(b[1, 1]),
            abs(abs(b[0, 1]) - a),
        )
    )


def own_selection(order, own, with_mask=True):
    """Landmarks and clock mask: the frozen nested selection on the own events,
    then the shared rule and the clock rule by direct boolean counts."""
    from intrinsic_weyl_joint import selection

    total = len(order)
    own_index = np.r_[: own + 1, total - 1]
    ids, _ = selection(order[np.ix_(own_index, own_index)], design="nested")
    ids = own_index[ids]
    inner = np.arange(1, own + 1)
    candidates = np.arange(own + 1, total - 1)
    below = order[np.ix_(inner, candidates)].sum(axis=0)
    above = order[np.ix_(candidates, inner)].sum(axis=1)
    kept = candidates[(below >= 2) & (above >= 2)]
    combined = np.r_[ids[:-1], kept, ids[-1]]
    if not with_mask:
        return combined, None
    mask = np.zeros((len(combined), len(combined)), dtype=bool)
    for i, j in itertools.combinations(range(len(combined)), 2):
        count = np.sum(order[combined[i], inner] & order[inner, combined[j]]) + np.sum(
            order[combined[j], inner] & order[inner, combined[i]]
        )
        mask[i, j] = mask[j, i] = count >= 2
    mask[0, -1] = mask[-1, 0] = False
    return combined, mask


def audit_selection(folder, prefix, w, fit_rows=None):
    """Re-derive one window's landmarks and shared rows from its selection record."""
    with np.load(
        folder / f"{prefix}_win{w}_selection.npz", allow_pickle=False
    ) as saved:
        size = int(saved["size"])
        order = np.unpackbits(saved["selection_order_packed"], axis=1)[:, :size].astype(
            bool
        )
        own = int(saved["selection_own"])
        origins = saved["candidate_source"]
        ids = saved["selected_ids"]
    combined, _ = own_selection(order, own, with_mask=False)
    rows = np.flatnonzero((combined > own) & (combined < size - 1))
    result = {
        "landmarks_reproduced": bool(np.array_equal(combined, ids)),
        "shared": int(len(rows)),
        "candidates": int(size - own - 2),
        "shared_rows": rows.tolist(),
        "shared_sources": origins[combined[rows] - own - 1].tolist(),
    }
    if fit_rows is not None:
        result["fit_shared_rows_reproduced"] = (
            fit_rows[0] == result["shared_rows"]
            and fit_rows[1] == result["shared_sources"]
        )
    return result


def independent_margin(first, last, amplitude, shift):
    """Causal margin from a direct ODE solution of both transverse Jacobi
    equations for the profile a(u) = A0 sin(pi (u + shift))."""
    from scipy.integrate import solve_ivp

    def rhs(u, y):
        a = amplitude * math.sin(math.pi * (u + shift))
        return [y[1], a * y[0], y[3], a * y[2], y[5], -a * y[4], y[7], -a * y[6]]

    end = solve_ivp(
        rhs,
        (first[0], last[0]),
        [1, 0, 0, 1, 1, 0, 0, 1],
        method="DOP853",
        rtol=1e-13,
        atol=1e-15,
    ).y[:, -1]
    cost = 0.0
    for axis, (c, dc, s, ds) in enumerate((end[:4], end[4:])):
        start, stop = first[2 + axis], last[2 + axis]
        slope = (stop - start * c) / s
        cost += 0.5 * (stop * (start * dc + slope * ds) - start * slope)
    return -(last[1] - first[1]) - cost


def audit_archives(folder, prefix, seed, levels, pairs=300):
    """Kept archives: clocks, CPU replay, ODE bits, exclusions, chart moves."""
    from intrinsic_weyl_bombelli_benchmark import benchmark_amplitude
    from intrinsic_weyl_second_continuum import exact_squared_interval, oriented_pairs
    from intrinsic_weyl_window_study import build_window

    amplitude = benchmark_amplitude()
    value = float(amplitude)
    truth, report = {}, {"prefix": prefix, "windows": {}}
    for w in range(5):
        with np.load(
            folder / f"truth_{prefix}_win{w}.npz", allow_pickle=False
        ) as saved:
            truth[w] = {key: saved[key] for key in saved.files}
    for w in range(5):
        shift = WINDOWS[w]
        wave = build_window("sinusoid", amplitude, 2.0, shift)
        top = levels[-1]
        with np.load(
            folder / f"{prefix}_win{w}_n{top}.npz", allow_pickle=False
        ) as saved:
            archive = {key: saved[key] for key in saved.files}
        nested = True
        for n in levels[:-1]:
            with np.load(
                folder / f"{prefix}_win{w}_n{n}.npz", allow_pickle=False
            ) as saved:
                for key in ("after", "before"):
                    nested &= np.array_equal(saved[key], archive[key][:, : n // 8])
                for key in ("selection_order", "selected_ids", "candidate_source"):
                    nested &= np.array_equal(saved[key], archive[key])
        combined, mask = own_selection(
            archive["selection_order"], int(archive["selection_own"])
        )
        points = truth[w]["points"]
        events = wave.sample(4096, np.random.default_rng([seed, w, 1]))
        stored = [
            np.unpackbits(archive[k][:, :512], axis=1).astype(bool)
            for k in ("after", "before")
        ]
        replayed = wave.relations(points[:, None], events[None])
        # ODE bits for the pairs nearest the flat cone and a random sample.
        left, right = points[:, None], events[None]
        du = right[..., 0] - left[..., 0]
        span = np.where(du == 0, 1.0, abs(du))
        flat = -np.sign(du) * (right[..., 1] - left[..., 1]) - (
            (right[..., 2:] - left[..., 2:]) ** 2
        ).sum(axis=-1) / (2 * span)
        nearest = np.argsort(abs(flat), axis=None)[:pairs]
        chosen = np.unique(
            np.r_[
                nearest,
                np.random.default_rng([seed, w, 2]).choice(flat.size, pairs, False),
            ]
        )
        matched = 0
        for landmark, event in zip(*np.unravel_index(chosen, flat.shape), strict=True):
            ahead = du[landmark, event] > 0
            first, last = (
                (points[landmark], events[event])
                if ahead
                else (events[event], points[landmark])
            )
            margin = independent_margin(first, last, value, shift)
            bits = (stored[0][landmark, event], stored[1][landmark, event])
            matched += bits == ((margin > 0) and ahead, (margin > 0) and not ahead)
        fit = read(folder / f"{prefix}_win{w}_n{top}_pure.json")["fit"]
        entry = {
            "nested_prefixes_bit_exact": bool(nested),
            "landmarks_reproduced": bool(
                np.array_equal(combined, archive["selected_ids"])
            ),
            "clock_count_reproduced": int(np.triu(mask, 1).sum()) == fit["clocks"],
            "truth_landmark_count_matches": bool(len(points) == len(combined)),
            "cpu_replay_first_4096_bit_exact": all(
                np.array_equal(a, b) for a, b in zip(replayed, stored, strict=True)
            ),
            "ode_pairs_checked": int(len(chosen)),
            "ode_bits_match": int(matched),
        }
        if w in INNER:
            inner = INNER[w]
            own_rows = int(truth[inner]["own_rows"])
            moved = truth[inner]["points"][1:own_rows] + np.array(
                [WINDOWS[inner] - shift, shift - WINDOWS[inner], 0.0, 0.0]
            )
            included = set(archive["candidate_source"].tolist())
            excluded_inside = 0
            for row, point in enumerate(moved, start=1):
                if row in included or not 0 < point[0] < 1:
                    continue
                lower = independent_margin(
                    np.array([0.0, 0.5, 0, 0]), point, value, shift
                )
                upper = independent_margin(
                    point, np.array([1.0, -0.5, 0, 0]), value, shift
                )
                excluded_inside += (lower > 0) and (upper > 0)
            own = int(truth[w]["own_rows"])
            sources = truth[w]["shared_source"]
            same = np.array_equal(
                points[own:-1],
                truth[inner]["points"][sources]
                + np.array([WINDOWS[inner] - shift, shift - WINDOWS[inner], 0.0, 0.0]),
            )
            shared_points = points[own:-1]
            theirs = truth[inner]["points"][sources]
            pairs_shared = oriented_pairs(shared_points)
            mine = exact_squared_interval(wave, shared_points, pairs_shared)
            other = exact_squared_interval(
                build_window("sinusoid", amplitude, 2.0, WINDOWS[inner]),
                theirs,
                pairs_shared,
            )
            entry.update(
                excluded_candidates_inside_window=int(excluded_inside),
                shared_points_are_moved_neighbour_landmarks=bool(same),
                shared_squared_interval_relative_difference=float(
                    np.max(abs(mine - other)) / np.max(abs(mine))
                ),
            )
        report["windows"][str(w)] = entry
    return report


def recompute_covariance(folder, prefix, n):
    """The centre window's clock mask and covariance from its archive with the
    inverse's CPU float64 second moment; the run used the bounded GPU variant
    and saved both, so they must agree bit for bit."""
    from intrinsic_weyl_joint import _clock_second_moment

    with np.load(
        folder / f"{prefix}_win{CENTRE}_n{n}.npz", allow_pickle=False
    ) as saved:
        order, own = saved["selection_order"], int(saved["selection_own"])
        after, before = saved["after"], saved["before"]
    ids, mask = own_selection(order, own)
    restricted = order[np.ix_(ids, ids)]
    pairs = np.column_stack(np.nonzero(np.triu(mask, 1)))
    forward = restricted[pairs[:, 0], pairs[:, 1]]
    a = np.where(forward, pairs[:, 0], pairs[:, 1])
    b = np.where(forward, pairs[:, 1], pairs[:, 0])
    p = (
        np.array(
            [
                np.bitwise_count(after[i] & before[j]).sum()
                for i, j in zip(a, b, strict=True)
            ]
        )
        / n
    )
    second = _clock_second_moment(after, before, a, b)
    # The operation order of read_order_clocks, so that equality can be exact.
    derivative = -0.5 / np.sqrt(p)
    covariance = (
        (second - n * np.outer(p, p)) / (n - 1) * np.outer(derivative, derivative) / n
    )
    with np.load(
        folder / f"{prefix}_win{CENTRE}_n{n}_clocks.npz", allow_pickle=False
    ) as saved:
        return {
            "clock_mask_equal": bool(np.array_equal(mask, saved["mask"])),
            "covariance_bit_identical": bool(
                np.array_equal(covariance, saved["covariance"])
            ),
            "clocks": int(len(a)),
        }


def recompute_operator(folder, prefix, n, recorded):
    """Minimum-norm operator from the centre fit and clocks (kept replicates)."""
    from intrinsic_weyl import second_derivative_basis
    from intrinsic_weyl_window_study import alias_map, minimum_norm_operator, whitening

    fit = read(folder / f"{prefix}_win{CENTRE}_n{n}_pure.json")["fit"]
    with np.load(
        folder / f"{prefix}_win{CENTRE}_n{n}_clocks.npz", allow_pickle=False
    ) as saved:
        mask, covariance = saved["mask"], saved["covariance"]
    alias, _, _ = alias_map(
        np.asarray(fit["coordinates_inferred_from_counts"]), mask, whitening(covariance)
    )
    operator, _ = minimum_norm_operator(alias @ second_derivative_basis())
    return float(np.max(abs(operator - recorded)) / np.max(abs(operator)))


def refit(folder, prefix, w, n):
    """Re-estimate one window in a fresh single-threaded process."""
    target = folder.with_name(folder.name + "_audit")
    target.mkdir(exist_ok=True)
    output = target / f"{prefix}_win{w}_n{n}_{int(time.time())}.json"
    env = {k: v for k, v in os.environ.items() if k not in THREAD_VARIABLES}
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    subprocess.run(
        [
            sys.executable,
            str(RUNNER),
            "--worker",
            str(folder / f"{prefix}_win{w}_n{n}.npz"),
            "--output",
            str(output),
        ],
        check=True,
        env=env,
    )
    new = read(output.with_name(output.stem + "_pure.json"))["fit"]
    old = read(folder / f"{prefix}_win{w}_n{n}_pure.json")["fit"]
    return {
        "prefix": prefix,
        "window": w,
        "measurement_n": n,
        "fields_identical": new == old,
        "differences_beyond_1e-8": differences(new, old),
    }


def lf_digest(data):
    """SHA-256 under the repository's LF storage-boundary rule."""
    return hashlib.sha256(data.replace(b"\r\n", b"\n")).hexdigest()


def sources_as_stamped(stamp, recorded):
    """The run's source digests against the stamped commit's blobs, and the
    bound sources of this checkout against the run's digests."""
    at_stamp, at_checkout = {}, {}
    for name, digest in recorded.items():
        blob = subprocess.run(
            ["git", "show", f"{stamp}:{name}"], cwd=ROOT, capture_output=True
        )
        at_stamp[name] = blob.returncode == 0 and lf_digest(blob.stdout) == digest
        path = ROOT / name
        at_checkout[name] = path.exists() and lf_digest(path.read_bytes()) == digest
    return at_stamp, at_checkout


def constants_as_protocol(summary, engineering):
    """Seeds, levels, windows and gate constants of the run against the
    protocol, not taken from the run."""
    if engineering:
        seeds, levels = list(ENGINEERING_SEEDS), list(ENGINEERING_LEVELS)
        kept = list(range(len(seeds)))
    else:
        seeds = [SEED + r for r in range(REPLICATES)]
        levels, kept = list(LEVELS), list(KEPT_ARCHIVES)
    haar = (
        f"default_rng([{HAAR_SEED}, replicate]), {HAAR_DRAWS} draws of the four "
        "side windows in the order of their offsets"
    )
    checks = {
        "seeds": [r["seed"] for r in summary["replicates"]] == seeds,
        "prefixes": [r["prefix"] for r in summary["replicates"]]
        == [f"rep{r:02d}" for r in range(len(seeds))],
        "levels": summary["measurement_levels"] == levels,
        "windows": summary["windows"] == list(WINDOWS),
        "spacing": summary["spacing_in_tip_separations"] == SPACING,
        "tolerance": summary["tolerance"] == TOLERANCE,
        "coverage": summary["coverage_two_sided"] == COVERAGE,
        "minimum_complete": summary["minimum_complete"]
        == (None if engineering else MINIMUM_COMPLETE),
        "minimum_shared": summary["minimum_shared_per_link"] == MINIMUM_SHARED,
        "kept_archives": summary["archives_kept_for_replicates"] == kept,
        "primary_correction": summary["primary_correction"] == CORRECTIONS["primary"],
        "haar_rng": summary["falsification"]["haar_rng"] == haar,
        "rows": sorted((r["replicate"], r["measurement_n"]) for r in summary["results"])
        == sorted((i, n) for i in range(len(seeds)) for n in levels),
        "clean_worktree_at_start": summary["worktree_clean_at_start"] is True,
        "protocol_bound": engineering or PROTOCOL in summary["source_sha256_lf"],
    }
    return all(checks.values()), checks


def audit(folder, receipt, refits):
    from intrinsic_weyl_joint import TIDAL_GRAM

    assert np.allclose(GRAM2, TIDAL_GRAM, atol=1e-14)
    summary = read(folder / "summary.json")
    stamp = summary["code_version"]
    engineering = summary["engineering_only"]
    levels = summary["measurement_levels"]
    kept = set(summary["archives_kept_for_replicates"])

    def git(*arguments):
        return subprocess.run(
            ["git", *arguments], cwd=ROOT, capture_output=True, text=True
        )

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
        "replicates": [],
        "selections": {},
        "levels": {},
    }
    result["constants_as_protocol"], result["constants"] = constants_as_protocol(
        summary, engineering
    )
    at_stamp, at_checkout = sources_as_stamped(stamp, summary["source_sha256_lf"])
    result["sources_match_stamp"] = all(at_stamp.values())
    result["sources_at_checkout_match_run"] = all(at_checkout.values())
    result["sources"] = {"at_stamp": at_stamp, "at_checkout": at_checkout}
    operators_ok = []
    for row in summary["results"]:
        prefix, n = row["prefix"], row["measurement_n"]
        for w in range(5):
            name = f"{prefix}_win{w}"
            if (
                name in result["selections"]
                or not (folder / f"{name}_selection.npz").exists()
            ):
                continue
            pure = folder / f"{name}_n{n}_pure.json"
            rows = None
            if pure.exists():
                fit = read(pure)["fit"]
                rows = (fit["shared_rows"], fit["shared_source_rows"])
            result["selections"][name] = audit_selection(folder, prefix, w, rows)
        if "failure" in row:
            result["replicates"].append(
                {
                    "replicate": row["replicate"],
                    "measurement_n": n,
                    "recorded_failure": row["failure"],
                }
            )
            continue
        report, truth = audit_replicate(folder, row)
        report["truth_orientation_error"] = truth_orientation(folder, row, truth)
        if row["replicate"] in kept:
            report["operator_relative_difference"] = recompute_operator(
                folder, prefix, n, report["primary_operator"]
            )
            report["centre_clocks"] = recompute_covariance(folder, prefix, n)
            operators_ok.append(
                report["operator_relative_difference"] < 1e-9
                and report["centre_clocks"]["clock_mask_equal"]
                and report["centre_clocks"]["covariance_bit_identical"]
            )
        report.pop("primary_operator", None)
        result["replicates"].append(report)
    for n in levels:
        audited = [r for r in result["replicates"] if r["measurement_n"] == n]
        done = [r for r in audited if "readings" in r]
        minimum = math.ceil(0.9 * len(audited)) if engineering else MINIMUM_COMPLETE
        restated = gates(
            [r["readings"]["primary"] for r in done],
            [r["readings"]["single"] for r in done],
            minimum,
            [r["counterfactuals"] for r in done],
            [r.get("oracle_parity_flags") for r in done],
        )
        entry = {"assigned": len(audited), "restated": restated}
        if not engineering:
            verdict = summary["levels"][str(n)]["verdict"]
            matches, entry["matches"] = matches_verdict(restated, verdict)
            entry["matches_run"] = bool(
                matches and len(audited) == REPLICATES == verdict["assigned"]
            )
        result["levels"][str(n)] = entry
    result["archives"] = [
        audit_archives(folder, f"rep{r:02d}", summary["replicates"][r]["seed"], levels)
        for r in sorted(kept)
    ]
    result["refits"] = [refit(folder, prefix, w, n) for prefix, w, n in refits]
    ok = [r for r in result["replicates"] if "readings" in r]
    result["passed"] = bool(
        result["run_stamp_is_ancestor_of_head"]
        and result["constants_as_protocol"]
        and result["sources_match_stamp"]
        and result["sources_at_checkout_match_run"]
        and all(
            s["landmarks_reproduced"] and s.get("fit_shared_rows_reproduced", True)
            for s in result["selections"].values()
        )
        and all(
            r["blind_boundary_kept"]
            and r["links_reproduced"]
            and r["every_link_has_eight_shared"]
            and r["estimates_reproduced"]
            and r["readings_reproduced"]
            and r["counterfactuals_reproduced"]
            and r["corrections_as_protocol"]
            and r["oracle_reproduced"] is not False
            and max(r["five_point_structure_residual"].values()) < 1e-9
            and r["truth_invariant_error"] < 1e-12
            and r["truth_orientation_error"] < 1e-9
            for r in ok
        )
        and all(operators_ok)
        and all(e.get("matches_run", True) for e in result["levels"].values())
        and all(
            entry["nested_prefixes_bit_exact"]
            and entry["landmarks_reproduced"]
            and entry["clock_count_reproduced"]
            and entry["truth_landmark_count_matches"]
            and entry["cpu_replay_first_4096_bit_exact"]
            and entry["ode_bits_match"] == entry["ode_pairs_checked"]
            and entry.get("excluded_candidates_inside_window", 0) == 0
            and entry.get("shared_points_are_moved_neighbour_landmarks", True)
            and entry.get("shared_squared_interval_relative_difference", 0) < 1e-10
            for a in result["archives"]
            for entry in a["windows"].values()
        )
        and all(r["fields_identical"] for r in result["refits"])
    )
    write_json(receipt, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", type=Path)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--refit", action="append", default=[], metavar="PREFIX:W:N")
    args = parser.parse_args()
    outcome = audit(
        args.folder.resolve(),
        args.receipt or args.folder.with_name(args.folder.name + ".audit.json"),
        [(p, int(w), int(n)) for p, w, n in (item.split(":") for item in args.refit)],
    )
    print(
        json.dumps({"levels": outcome["levels"], "passed": outcome["passed"]}, indent=1)
    )
    sys.exit(0 if outcome["passed"] else 1)
