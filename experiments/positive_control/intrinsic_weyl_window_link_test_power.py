"""Geometry, power and cost of a finite-N test in which the side windows and the
order-only links of the five-window correction carry the result (design aid).

Design aid for ``docs/paper/paper_a/intrinsic_weyl/window_link_test_protocol.md``
(the draft until the freeze). No measurement event is drawn and nothing here is
a test result.

Why. The frozen test (``window_test_protocol.md``, commit ``1ca7af2``) put its
side windows at proper-time offsets +-1/2 and +-1 of the tip separation on the
half-cycle wave (wavelength 2). The +-1/2 windows sit on nodes and read about
zero, the +-1 windows carry little weight, so the corrected reading was about
1.14 times the centre window's own estimate; deleting the side windows or
replacing the links would have passed as well. Here the window offsets are
0, +-h, +-2h and the question is which spacing h makes the side windows and
their links carry the reading, at what variance and cost.

Parts:

1. Geometry, deterministic continuum. Fifth-degree continuum fits of every
   window of every candidate spacing on the benchmark landmark layouts (cases
   0 and 1, run clock covariances), exactly as the multi-window study fits
   them (``intrinsic_weyl_window_study.window_job``), each window mapped to the
   generator frame by its oracle frame. Per spacing: the corrected reading,
   the centre term alone, the per-window contributions, the readings with
   parity-flipped, Haar-random and fixed-angle side links, and the variance
   ratios to one window.
2. Empirical, finite N, post hoc on the frozen run (64 replicates at 2^25).
   Each replicate's operator B is recovered from its saved stencil weights and
   re-used with the stencil of each candidate spacing. The centre-term-only
   reading is then exact for that replicate's centre window. The corrected
   reading's spread uses the recorded windows' noise with the signals of the
   new positions (each linked side estimate shifted by the change of its
   continuum reading times the truth). At h = 1/2 the swap is the identity and
   must reproduce the recorded readings.
3. Selection dry run: selection events only (4096 per window, no measurement
   event, no fit), design seeds disjoint from every test seed, for the
   landmark, shared-landmark and clock counts that set the fit cost. It
   replays the frozen runner's selection steps and is checked bit for bit
   against the frozen engineering pilot's selection records at h = 1/2.
4. Power and cost. Replicates for 90% power of the accuracy gate on its 95%
   two-sided Student-t interval (house rule), the power of the control,
   falsification and link-sign gates, the power of the secondary side-share
   gate G4 against its band (computed last, so that every earlier draw is
   unchanged), and GPU and wall hours from the frozen run's measured
   generation and fit times.
"""

import argparse
import hashlib
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUTS = ROOT / "outputs"  # override with --outputs
WAVELENGTH = 2.0  # diamond heights: the crest-centred half cycle of the benchmark
SPACINGS = {
    "1/4": 0.25,
    "5/16": 0.3125,
    "1/3": 1 / 3,
    "3/8": 0.375,
    "7/16": 0.4375,
    "1/2": 0.5,
}
CHOSEN = "1/3"
FROZEN = "1/2"
OFFSETS = (-2, -1, 0, 1, 2)
CASES = (0, 1)
ORIENTATIONS, ORIENTATION_SEED = 400, 2026112107  # as the frozen design aid
HAAR_DRAWS, HAAR_SEED = 4096, 2026112951
ANGLES_DEG = (5, 10, 20, 30, 45, 60, 90)
ANGLE_DRAWS = 256
# Analysis seeds are 202611295x. Selection dry run: design seeds, disjoint from
# the frozen test (2026112600+),
# its pilot (2026112700, 2026112701) and the proposed test seeds (2026112800+).
DRY_SEEDS = (2026112900, 2026112901, 2026112902, 2026112903)
VALIDATION_SEED = 2026112700  # frozen engineering pilot, replicate 0
SELECTION_N = 4096
TIPS = np.array([[0.0, 0.5, 0.0, 0.0], [1.0, -0.5, 0.0, 0.0]])
GENERATION_ORDER = (2, 1, 3, 0, 4)
INNER = {0: 1, 1: 2, 3: 2, 4: 3}
FROZEN_LEVEL = 33554432
FROZEN_REPLICATES = 64
# Power.
TOLERANCES = (0.04, 0.05, 0.06)
REPLICATES = (
    32,
    40,
    48,
    56,
    64,
    72,
    80,
    88,
    96,
    104,
    112,
    120,
    128,
    160,
    192,
    256,
)
LEVELS = (2**24, 2**25, 2**26)
KAPPAS = (1.0, 1.2)
OFFSETS_MEAN = (0.02, 0.01, 0.0, -0.01, -0.02)
# The design's sizing scenarios (kappa, finite-N offset of the mean reading).
SIZING_SCENARIOS = ((1.2, 0.0), (1.0, 0.01), (1.0, -0.01))
TARGET_POWER = 0.9
DRAWS = 20000
POWER_SEED = 2026112952
# G4 (secondary): the side windows' share a - a0 inside this band; the continuum
# layouts give -0.2018 and -0.1904, widened by 0.02 on each side. Its power is
# also reported at a stress standard deviation and at finite-N mean offsets.
SIDE_SHARE_BAND = (-0.22, -0.17)
SIDE_SHARE_OFFSETS = (-0.01, 0.0, 0.01)
SIDE_SHARE_STRESS_SD = 0.03
# Cost, measured in the frozen run at 2^25 (window_test_results.md): generation
# per window is nearly independent of the landmark count (centre median 94 s
# with 130 landmarks, side 96-98 s with 168); a fit's wall time grows as
# clocks^1.63 (fitted over its 320 windows); the clock moments take 85 s of
# the centre fit (149 s, 3,937 clocks) and 128 s of the side fit (311 s, 6,209).
GENERATION_SECONDS_2_25 = 96.0
CENTRE_FIT = {"clocks": 3936.5, "seconds": 148.7, "moment_seconds": 85.4}
SIDE_FIT = {"clocks": 6209.0, "seconds": 310.7, "moment_seconds": 127.7}
CLOCK_EXPONENT = 1.63
FIT_WORKERS = (3, 4, 6)


def shift_key(shift):
    return f"{shift:+.12f}"


def positions(h):
    return [k * h for k in OFFSETS]


def git_state():
    import subprocess

    def run(*args):
        return subprocess.run(
            ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
        ).stdout.strip()

    return run("rev-parse", "--short=7", "HEAD"), run(
        "status", "--porcelain", "--untracked-files=all"
    ) == ""


def lf_sha(path):
    return hashlib.sha256(Path(path).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


# ---------------------------------------------------------------------------
# Stencil weights for any spacing, and the operator from saved weights
# ---------------------------------------------------------------------------


def stencil_weights(operator, h):
    """Five-point second-order weights W_s = d_s0 I - B d2_s/h^2 + B^2 d4_s/h^4.

    The formula of ``intrinsic_weyl_window_study.correction_specs``, which only
    offers h in (1/4, 1/2, 1); checked equal to it bit for bit at those h.
    """
    from intrinsic_weyl_window_study import STENCILS

    offsets, d2, d4 = STENCILS["five_point"]
    first = [
        (o == 0) * np.eye(10) - operator * w / h**2
        for o, w in zip(offsets, d2, strict=True)
    ]
    return [w + operator @ operator * v / h**4 for w, v in zip(first, d4, strict=True)]


def operator_from_half_weights(weights):
    """B from the five-point weights at h = 1/2: W(-1/2) + 4 W(-1) = -4 B."""
    w = {s: np.asarray(m, dtype=float) for s, m in weights.items()}
    low = -(w[-0.5] + 4 * w[-1.0]) / 4
    high = -(w[0.5] + 4 * w[1.0]) / 4
    centre = w[0.0] - np.eye(10) - 10 * low - 96 * low @ low
    scale = np.max(abs(low))
    return low, {
        "sides_relative_difference": float(np.max(abs(low - high)) / scale),
        "centre_weight_relative_residual": float(np.max(abs(centre)) / scale),
    }


def haar_o3(rng):
    q, r = np.linalg.qr(rng.standard_normal((3, 3)))
    return q * np.sign(np.diag(r))


def weyl_bank(rng, draws):
    from intrinsic_weyl_window_study import weyl_rotation

    return np.array([weyl_rotation(haar_o3(rng)) for _ in range(draws)])


def selfcheck():
    """Stencil formula, operator recovery, parity map, Haar bank."""
    from intrinsic_weyl_window_study import correction_specs, weyl_rotation

    rng = np.random.default_rng(2026112953)
    operator = 0.01 * rng.normal(size=(10, 10))
    alias = rng.normal(size=(10, 100))
    report = {}
    for h in (0.25, 0.5):
        offsets = [round(o * h, 10) for o in OFFSETS]
        specs = correction_specs(
            sorted(set(offsets)),
            alias,
            np.array([1.0, 0, 0, 1]),
            np.array([1.0, 0, 0, 1]),
            operator,
        )
        _, _, _, shifts, weights = next(
            s
            for s in specs
            if s[:3] == ("five_point_second", h, "minimum_norm_time_time")
        )
        ours = stencil_weights(operator, h)
        assert shifts == offsets
        difference = max(
            float(np.max(abs(a - b))) for a, b in zip(ours, weights, strict=True)
        )
        report[f"stencil_equals_correction_specs_h_{h}"] = difference
        assert difference == 0.0
    weights = dict(
        zip((-1.0, -0.5, 0.0, 0.5, 1.0), stencil_weights(operator, 0.5), strict=True)
    )
    recovered, checks = operator_from_half_weights(weights)
    report["operator_recovery_error"] = float(np.max(abs(recovered - operator)))
    assert report["operator_recovery_error"] < 1e-13
    assert max(checks.values()) < 1e-10
    parity = weyl_rotation(-np.eye(3))
    expected = np.diag(np.r_[np.ones(5), -np.ones(5)])
    report["parity_map_is_electric_plus_magnetic_minus"] = float(
        np.max(abs(parity - expected))
    )
    assert report["parity_map_is_electric_plus_magnetic_minus"] < 1e-12
    bank = weyl_bank(np.random.default_rng(2026112954), 2048)
    report["haar_bank_mean_max_abs"] = float(np.max(abs(bank.mean(axis=0))))
    assert report["haar_bank_mean_max_abs"] < 0.1
    return report


# ---------------------------------------------------------------------------
# Part 1: continuum geometry of each spacing
# ---------------------------------------------------------------------------


def fit_jobs(archive, clock_cache, amplitude):
    unique = {}
    for h in SPACINGS.values():
        for s in positions(h):
            unique.setdefault(shift_key(s), s)
    jobs = []
    for case in CASES:
        for key in sorted(unique, key=float):
            jobs.append(
                {
                    "case": case,
                    "amplitude": str(amplitude),
                    "clocks": str(
                        Path(clock_cache) / f"run_clocks_case{case}_n16777216.npz"
                    ),
                    "truth": str(Path(archive) / f"truth_{case}.npz"),
                    "kind": "sinusoid",
                    "wavelength": WAVELENGTH,
                    "shift": unique[key],
                    "orders": (5,),
                }
            )
    return jobs


def run_fit(job):
    from intrinsic_weyl_window_study import window_job

    return window_job(job)


def continuum_fits(archive, clock_cache, amplitude, workers, threads, cache=None):
    jobs = fit_jobs(archive, clock_cache, amplitude)
    results, todo = {}, []
    for job in jobs:
        key = (job["case"], shift_key(job["shift"]))
        path = (
            None if cache is None else Path(cache) / f"fit_case{key[0]}_{key[1]}.json"
        )
        if path is not None and path.exists():
            results[key] = json.loads(path.read_text(encoding="utf-8"))
        else:
            todo.append((key, job, path))
    for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ[name] = str(threads)
    begin = time.time()
    if todo:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for (key, _, path), window in zip(
                todo, pool.map(run_fit, [job for _, job, _ in todo]), strict=True
            ):
                results[key] = window
                if path is not None:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text(json.dumps(window), encoding="utf-8")
                print("fit", key, window["order_5"]["seconds"], "s", flush=True)
    return results, {
        "fits": len(jobs),
        "computed": len(todo),
        "seconds": time.time() - begin,
    }


def window_study_agreement(results, window_study):
    """Shared windows against the multi-window study's recorded fits."""
    recorded = json.loads(
        (Path(window_study) / "summary_windows.json").read_text(encoding="utf-8")
    )
    worst, compared = 0.0, 0
    for window in recorded:
        if window["kind"] != "sinusoid" or window["wavelength"] != WAVELENGTH:
            continue
        key = (window["case"], shift_key(window["shift"]))
        if key not in results:
            continue
        ours = np.asarray(results[key]["order_5"]["coefficients"])
        theirs = np.asarray(window["order_5"]["coefficients"])
        worst = max(worst, float(np.max(abs(ours - theirs)) / np.max(abs(theirs))))
        compared += 1
    return {"windows_compared": compared, "max_relative_coefficient_difference": worst}


def case_geometry(case, fits, clock_cache, amplitude, frozen_power):
    """All continuum quantities of one landmark layout, every spacing."""
    from intrinsic_weyl import second_derivative_basis
    from intrinsic_weyl_joint import TIDAL_GRAM
    from intrinsic_weyl_window_study import (
        alias_map,
        axis_rotation,
        covector,
        minimum_norm_operator,
        plane_wave_operator,
        principal_null_direction,
        second_rotation,
        tidal_norm,
        unit_truth,
        weyl_rotation,
        whitening,
    )
    from scipy.spatial.transform import Rotation

    with np.load(Path(clock_cache) / f"run_clocks_case{case}_n16777216.npz") as saved:
        mask, covariance = saved["mask"], saved["cov"]
    truth = amplitude * unit_truth()
    centre = fits[(case, shift_key(0.0))]
    frame0 = np.asarray(centre["frame"])
    inferred = np.asarray(centre["order_5"]["coordinates"])
    alias_g, consistency, _ = alias_map(inferred, mask, whitening(covariance))
    alias_n = weyl_rotation(frame0.T) @ alias_g @ second_rotation(frame0)
    minimum_norm, rank = minimum_norm_operator(alias_n @ second_derivative_basis())

    def in_normal(window):
        rotation = weyl_rotation(np.asarray(window["frame"]).T)
        fit = window["order_5"]
        return (
            rotation @ np.asarray(fit["coefficients"]),
            rotation @ np.asarray(fit["covariance"]) @ rotation.T,
        )

    windows = {key: in_normal(w) for (c, key), w in fits.items() if c == case}
    n_hat, _ = principal_null_direction(windows[shift_key(0.0)][0])
    plane = alias_n @ plane_wave_operator(covector(n_hat))
    gram = TIDAL_GRAM @ truth
    norm2 = float(truth @ gram)

    def along(vector, target=None):
        if target is None:
            return float(gram @ vector) / norm2
        g = TIDAL_GRAM @ target
        return float(g @ vector) / float(target @ g)

    def variance(weights, covariances, target):
        g = TIDAL_GRAM @ target
        total = sum(
            g @ w @ s @ w.T @ g for w, s in zip(weights, covariances, strict=True)
        )
        return float(total) / float(target @ g) ** 2

    orientations = [
        weyl_rotation(r)
        for r in Rotation.random(
            ORIENTATIONS, random_state=ORIENTATION_SEED
        ).as_matrix()
    ]
    rng = np.random.default_rng([HAAR_SEED, case])
    bank = weyl_bank(rng, HAAR_DRAWS)
    parity = weyl_rotation(-np.eye(3))
    single_var = [
        variance([np.eye(10)], [windows[shift_key(0.0)][1]], o @ truth)
        for o in orientations
    ]
    # Fixed-angle link errors: one bank per angle, shared by all spacings.
    angle_rng = np.random.default_rng([HAAR_SEED, case, 1])
    angle_banks = {
        angle: [
            [
                weyl_rotation(
                    axis_rotation(angle_rng.normal(size=3), np.radians(angle))
                )
                for _ in range(4)
            ]
            for _ in range(ANGLE_DRAWS)
        ]
        for angle in ANGLES_DEG
    }
    haar_sets = HAAR_DRAWS // 5
    result = {
        "case": case,
        "alias_column_consistency": consistency,
        "minimum_norm_slice_rank": rank,
        "operator_along_truth": along(minimum_norm @ truth),
        "single_reading": along(windows[shift_key(0.0)][0]),
        "single_reading_sd_at_2_24_nominal": math.sqrt(
            variance([np.eye(10)], [windows[shift_key(0.0)][1]], truth)
        ),
        "spacings": {},
    }
    for label, h in SPACINGS.items():
        keys = [shift_key(s) for s in positions(h)]
        estimates = [windows[k][0] for k in keys]
        covariances = [windows[k][1] for k in keys]
        weights = stencil_weights(minimum_norm, h)
        plane_weights = stencil_weights(plane, h)
        centre_index = OFFSETS.index(0)
        sides = [i for i in range(5) if i != centre_index]
        inner = [i for i in range(5) if abs(OFFSETS[i]) == 1]
        outer = [i for i in range(5) if abs(OFFSETS[i]) == 2]

        def reading(
            subset,
            transform=None,
            use=weights,
            estimates=estimates,
            centre=centre_index,
        ):
            total = np.zeros(10)
            for i in subset:
                c = estimates[i]
                if transform is not None and i != centre:
                    c = transform(i) @ c
                total = total + use[i] @ c
            return along(total)

        primary = reading(range(5))
        estimate = sum(w @ c for w, c in zip(weights, estimates, strict=True))
        # Haar-random side links. The spread over draws of the projected
        # signals does not shrink with N; the noise part does.
        haar = np.array(
            [
                reading(range(5), transform=lambda i, d=d: bank[5 * d + i])
                for d in range(haar_sets)
            ]
        )
        haar_noise_var = float(
            np.mean(
                [
                    variance(
                        weights,
                        [
                            covariances[i]
                            if i == centre_index
                            else bank[5 * d + i] @ covariances[i] @ bank[5 * d + i].T
                            for i in range(5)
                        ],
                        truth,
                    )
                    for d in range(256)
                ]
            )
        )
        angles = {}
        for angle in ANGLES_DEG:
            values = [
                reading(
                    range(5),
                    transform=lambda i, t=turns, s=sides: t[s.index(i)],
                )
                for turns in angle_banks[angle]
            ]
            angles[str(angle)] = {
                "mean": float(np.mean(values)),
                "min": float(np.min(values)),
                "max": float(np.max(values)),
            }
        ratios = {
            name: [
                variance(w, s, o @ truth) / v
                for o, v in zip(orientations, single_var, strict=True)
            ]
            for name, (w, s) in {
                "primary": (weights, covariances),
                "centre_term_only": (
                    [weights[centre_index]],
                    [covariances[centre_index]],
                ),
                "side_share": (
                    [weights[i] for i in sides],
                    [covariances[i] for i in sides],
                ),
                "five_point_plane_wave": (plane_weights, covariances),
            }.items()
        }
        single_nominal = variance([np.eye(10)], [covariances[centre_index]], truth)
        result["spacings"][label] = {
            "h": h,
            "shifts": positions(h),
            "window_reading_relative_to_centre_truth": {
                str(OFFSETS[i]): along(estimates[i]) for i in range(5)
            },
            "weight_along_truth": {
                str(OFFSETS[i]): along(weights[i] @ truth) for i in range(5)
            },
            "contribution": {
                str(OFFSETS[i]): along(weights[i] @ estimates[i]) for i in range(5)
            },
            "readings": {
                "single": along(estimates[centre_index]),
                "primary_minimum_norm": primary,
                "five_point_plane_wave": reading(range(5), use=plane_weights),
                "centre_term_only": reading([centre_index]),
                "side_share": reading(sides),
                "outer_pair_deleted": reading([centre_index, *inner]),
                "inner_pair_deleted": reading([centre_index, *outer]),
                "parity_flipped_side_links": reading(
                    range(5), transform=lambda i: parity
                ),
                "haar_side_links_mean": float(haar.mean()),
            },
            "primary_full_tidal_relative_error": tidal_norm(estimate - truth)
            / tidal_norm(truth),
            "link_angle_sensitivity_primary": angles,
            "variance_ratio_to_single": {
                name: {
                    "median": float(np.median(v)),
                    "5_95": np.quantile(v, [0.05, 0.95]).tolist(),
                }
                for name, v in ratios.items()
            },
            "haar_side_links_spread": {
                "noise_variance_ratio_to_single": haar_noise_var / single_nominal,
                "signal_sd_independent_of_n": float(haar.std(ddof=1)),
            },
        }
        if label == FROZEN and frozen_power is not None:
            recorded = frozen_power["continuum"][str(case)]
            result["spacings"][label]["frozen_design_aid_agreement"] = {
                "single": abs(
                    recorded["single"]["continuum_reading"]
                    - along(estimates[centre_index])
                ),
                "primary": abs(
                    recorded["five_point_minimum_norm"]["continuum_reading"] - primary
                ),
                "variance_ratio_median": abs(
                    recorded["five_point_minimum_norm"][
                        "variance_ratio_to_single_median"
                    ]
                    - float(np.median(ratios["primary"]))
                ),
            }
    return result


# ---------------------------------------------------------------------------
# Part 2: the frozen run re-weighted (finite N, post hoc)
# ---------------------------------------------------------------------------


def frozen_empirical(run, summary_path, rho_old, rho_new):
    """Per spacing: exact centre-term-only reading and the swapped spread."""
    from intrinsic_weyl_joint import TIDAL_GRAM
    from intrinsic_weyl_window_study import weyl_rotation

    summary = json.loads(Path(summary_path).read_text(encoding="utf-8"))
    recorded = {row["replicate"]: row["readings"] for row in summary["results"]}
    old_shifts = (-1.0, -0.5, 0.0, 0.5, 1.0)
    rng = np.random.default_rng([HAAR_SEED, 7])
    bank = weyl_bank(rng, 1024)
    digests = {}
    rows = {label: [] for label in SPACINGS}
    operator_along, checks, recovery = [], [], []
    for r in range(FROZEN_REPLICATES):
        prefix = f"rep{r:02d}"
        paths = {
            "correction": Path(run) / f"{prefix}_n{FROZEN_LEVEL}_combined_pure.json",
            "evaluation": Path(run) / f"{prefix}_n{FROZEN_LEVEL}_combined.json",
            **{
                f"fit{w}": Path(run) / f"{prefix}_win{w}_n{FROZEN_LEVEL}_pure.json"
                for w in range(5)
            },
        }
        data = {k: json.loads(p.read_text(encoding="utf-8")) for k, p in paths.items()}
        digests.update(
            {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths.values()}
        )
        correction = data["correction"]["correction"]
        truth = np.asarray(
            data["evaluation"]["evaluation"]["truth_coefficients_centre_frame"],
            dtype=float,
        )
        gram = TIDAL_GRAM @ truth
        norm2 = float(truth @ gram)

        def along(vector, gram=gram, norm2=norm2):
            return float(gram @ vector) / norm2

        primary = correction["corrections"]["primary"]
        weights = {
            float(s): np.asarray(m, dtype=float)
            for s, m in zip(primary["shifts"], primary["weights"], strict=True)
        }
        links = {
            s: weyl_rotation(np.asarray(correction["links"]["rotations"][f"{s:+.1f}"]))
            for s in old_shifts
        }
        fits = [
            np.asarray(data[f"fit{w}"]["fit"]["coefficients"], dtype=float)
            for w in range(5)
        ]
        linked = [links[s] @ fits[i] for i, s in enumerate(old_shifts)]
        own = fits
        replayed = sum(along(weights[s] @ linked[i]) for i, s in enumerate(old_shifts))
        checks.append(
            max(
                abs(replayed - recorded[r]["primary"]["value"]),
                abs(along(fits[2]) - recorded[r]["single"]["value"]),
            )
        )
        operator, consistency = operator_from_half_weights(weights)
        recovery.append(max(consistency.values()))
        operator_along.append(along(operator @ truth))
        for label, h in SPACINGS.items():
            new = stencil_weights(operator, h)
            swapped = [
                linked[i] + (rho_new[label][i] - rho_old[i]) * truth for i in range(5)
            ]
            if label == FROZEN:
                # The swap is the identity here: the recorded reading returns.
                checks.append(
                    abs(
                        sum(along(new[i] @ swapped[i]) for i in range(5))
                        - recorded[r]["primary"]["value"]
                    )
                )
            haar = [
                along(
                    new[2] @ swapped[2]
                    + sum(new[i] @ bank[5 * d + i] @ swapped[i] for i in (0, 1, 3, 4))
                )
                for d in range(len(bank) // 5)
            ]
            rows[label].append(
                {
                    "single": along(fits[2]),
                    "centre_term_only": along(new[2] @ fits[2]),
                    "primary_swapped": sum(
                        along(new[i] @ swapped[i]) for i in range(5)
                    ),
                    "side_share_swapped": sum(
                        along(new[i] @ swapped[i]) for i in (0, 1, 3, 4)
                    ),
                    "identity_links_unswapped": sum(
                        along(new[i] @ own[i]) for i in range(5)
                    ),
                    "haar_side_links_swapped_mean": float(np.mean(haar)),
                    "haar_side_links_swapped_sd": float(np.std(haar, ddof=1)),
                    "contribution_swapped": {
                        str(OFFSETS[i]): along(new[i] @ swapped[i]) for i in range(5)
                    },
                }
            )
    listing = "\n".join(f"{name} {digests[name]}" for name in sorted(digests))

    def stats(values):
        values = np.asarray(values, dtype=float)
        return {"mean": float(values.mean()), "sd": float(values.std(ddof=1))}

    out = {
        "post_hoc": True,
        "run": "outputs/intrinsic_weyl_window_test_v1 (frozen test, commit 1ca7af2)",
        "files": len(digests),
        "sha256_of_file_digests": hashlib.sha256(listing.encode()).hexdigest(),
        "largest_replay_difference": float(max(checks)),
        "operator_recovery_largest_relative_residual": float(max(recovery)),
        "operator_along_truth": stats(operator_along),
        "note": (
            "centre_term_only is exact for these centre windows; primary_swapped and "
            "the Haar variant use the recorded windows' noise with the continuum "
            "signals of the new positions, so their spread (not their mean) is the "
            "finite-N input; identity_links_unswapped keeps the old positions' "
            "signals and only its mean is meaningful (random relative frames)."
        ),
        "spacings": {},
    }
    for label in SPACINGS:
        block = rows[label]
        out["spacings"][label] = {
            name: stats([row[name] for row in block])
            for name in (
                "single",
                "centre_term_only",
                "primary_swapped",
                "side_share_swapped",
                "identity_links_unswapped",
                "haar_side_links_swapped_mean",
            )
        }
        out["spacings"][label]["haar_side_links_swapped_within_replicate_sd_mean"] = (
            float(np.mean([row["haar_side_links_swapped_sd"] for row in block]))
        )
        out["spacings"][label]["haar_side_links_total_sd"] = float(
            np.sqrt(
                np.var([row["haar_side_links_swapped_mean"] for row in block], ddof=1)
                + np.mean([row["haar_side_links_swapped_sd"] ** 2 for row in block])
            )
        )
        out["spacings"][label]["contribution_swapped"] = {
            str(k): stats([row["contribution_swapped"][str(k)] for row in block])
            for k in OFFSETS
        }
    return out


# ---------------------------------------------------------------------------
# Part 3: selection dry run (no measurement event, no fit)
# ---------------------------------------------------------------------------


def dry_window(wave, seed, w, shift, neighbour):
    """The frozen runner's selection steps 1-4 for one window."""
    from intrinsic_weyl_window_test import window_selection

    rng = np.random.default_rng([seed, w, 0])
    points = np.vstack((TIPS[:1], wave.sample(SELECTION_N, rng), TIPS[1:]))
    permutation = np.r_[
        0, rng.permutation(np.arange(1, len(points) - 1)), len(points) - 1
    ]
    points = points[permutation]
    candidates, origins = np.empty((0, 4)), np.empty(0, dtype=np.int64)
    record = {"window": w, "shift": shift}
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
        record.update(
            candidates_offered=int(len(rows)), candidates_in_interval=int(len(inside))
        )
    extended = np.vstack((points[:-1], candidates, points[-1:]))
    order = wave.order(extended)
    ids, mask, kept = window_selection(order, SELECTION_N)
    landmarks = extended[ids]
    own_rows = int(np.sum(ids <= SELECTION_N))
    record.update(
        landmarks_including_tips=int(len(landmarks)),
        own_landmarks_including_first_tip=own_rows,
        shared_landmarks=int(len(kept)),
        clocks=int(np.triu(mask, 1).sum()),
    )
    return record, (landmarks, own_rows, shift), ids, origins


def dry_replicate(task):
    """Selection of the five windows of one replicate at spacing h."""
    from intrinsic_weyl_bombelli_benchmark import benchmark_amplitude
    from intrinsic_weyl_window_study import build_window

    seed, h = task
    amplitude = benchmark_amplitude()
    shifts = positions(h)
    waves = {
        w: build_window("sinusoid", amplitude, WAVELENGTH, shifts[w]) for w in range(5)
    }
    info, records, selected = {}, {}, {}
    begin = time.time()
    for w in GENERATION_ORDER:
        neighbour = info.get(INNER[w]) if w != 2 else None
        record, info[w], ids, origins = dry_window(
            waves[w], seed, w, shifts[w], neighbour
        )
        records[w] = record
        selected[w] = (ids.tolist(), origins.tolist())
    return {
        "seed": seed,
        "h": h,
        "windows": [records[w] for w in range(5)],
        "seconds": time.time() - begin,
    }, selected


def dry_run(workers, engineering):
    tasks = [(VALIDATION_SEED, 0.5)] + [(seed, SPACINGS[CHOSEN]) for seed in DRY_SEEDS]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        outcomes = list(pool.map(dry_replicate, tasks))
    validation, selected = outcomes[0]
    agreement = {}
    for w in range(5):
        with np.load(Path(engineering) / f"rep00_win{w}_selection.npz") as saved:
            ids_equal = bool(
                np.array_equal(saved["selected_ids"], np.asarray(selected[w][0]))
            )
            source_equal = bool(
                np.array_equal(saved["candidate_source"], np.asarray(selected[w][1]))
            )
        agreement[str(w)] = {
            "selected_ids_equal": ids_equal,
            "candidate_source_equal": source_equal,
        }
    if not all(all(v.values()) for v in agreement.values()):
        raise ArithmeticError(
            f"Dry run does not replay the pilot's selection: {agreement}"
        )
    design = [outcome for outcome, _ in outcomes[1:]]

    def spread(key, side):
        values = [
            window[key]
            for replicate in design
            for window in replicate["windows"]
            if (window["window"] != 2) == side
        ]
        return {
            "min": int(min(values)),
            "median": float(np.median(values)),
            "max": int(max(values)),
        }

    return {
        "validation_h_1_2_seed": VALIDATION_SEED,
        "validation_against_engineering_pilot": agreement,
        "validation_windows": validation["windows"],
        "design_seeds": list(DRY_SEEDS),
        "spacing": CHOSEN,
        "replicates": design,
        "summary": {
            side_label: {
                key: spread(key, side)
                for key in ("landmarks_including_tips", "shared_landmarks", "clocks")
            }
            for side_label, side in (("centre", False), ("side", True))
        },
    }


# ---------------------------------------------------------------------------
# Part 4: power and cost
# ---------------------------------------------------------------------------


def interval_power(
    mean, sd, replicates, rule, rng, tolerance, draws=DRAWS, band=SIDE_SHARE_BAND
):
    """P(rule) on the 95% two-sided t-interval of the mean of R readings.

    ``sd`` is a scalar or an array of per-replicate standard deviations to
    resample from (selection-by-orientation spread).
    """
    from scipy.stats import t as student

    quantile = student.ppf(0.975, replicates - 1)
    sigma = (
        np.full((draws, replicates), float(sd))
        if np.ndim(sd) == 0
        else rng.choice(np.asarray(sd), size=(draws, replicates))
    )
    values = mean + sigma * rng.normal(size=(draws, replicates))
    centre = values.mean(axis=1)
    half = quantile * values.std(axis=1, ddof=1) / np.sqrt(replicates)
    low, high = centre - half, centre + half
    if rule == "inside":
        hit = (low > 1 - tolerance) & (high < 1 + tolerance)
    elif rule == "outside":
        hit = (high < 1 - tolerance) | (low > 1 + tolerance)
    elif rule == "below":
        hit = high < 1 - tolerance
    elif rule == "positive":
        hit = low > 0
    elif rule == "negative":
        hit = high < 0
    elif rule == "inside_band":
        hit = (low > band[0]) & (high < band[1])
    else:
        raise ValueError(rule)
    return float(hit.mean())


def calibration_spread(calibration, amplitude):
    """One-window reading sd over flat calibration selections x orientations (2^24)."""
    sys.path.insert(0, str(HERE))
    from intrinsic_weyl_window_test_power import single_window_spread
    from scipy.spatial.transform import Rotation

    orientations = list(
        Rotation.random(ORIENTATIONS, random_state=ORIENTATION_SEED).as_matrix()
    )
    return np.sqrt(
        single_window_spread(Path(calibration), amplitude, orientations[:100])
    )


def side_fit_seconds(clocks, level):
    scale = (clocks / SIDE_FIT["clocks"]) ** CLOCK_EXPONENT
    moment = SIDE_FIT["moment_seconds"] / SIDE_FIT["seconds"]
    return SIDE_FIT["seconds"] * scale * (moment * level / 2**25 + (1 - moment))


def centre_fit_seconds(level):
    moment = CENTRE_FIT["moment_seconds"] / CENTRE_FIT["seconds"]
    return CENTRE_FIT["seconds"] * (moment * level / 2**25 + (1 - moment))


def replicate_cost(level, side_clocks):
    generation = 5 * GENERATION_SECONDS_2_25 * level / 2**25
    fits = centre_fit_seconds(level) + 4 * side_fit_seconds(side_clocks, level)
    return {
        "gpu_busy_seconds": generation,
        "fit_cpu_seconds": fits,
        "wall_seconds_by_fit_workers": {
            str(n): max(generation, fits / n) for n in FIT_WORKERS
        },
    }


def power_and_cost(geometry, empirical, dry, calibration_sd, realized):
    """Replicates for 90% power of G1, the other gates' power, and the cost.

    Standard deviations at 2^25 for each reading are the largest of the
    available estimates: the frozen run's realized single-window spread times
    the continuum variance ratio, the calibration's rms (2^24, scaled by
    1/sqrt 2) times the same ratio, and the frozen run re-weighted (Part 2).
    Per-replicate standard deviations are resampled from the calibration's
    selection-by-orientation spread rescaled to that value; kappa inflates the
    variance; the Haar variant's projected-signal spread does not shrink with N.
    """
    rng = np.random.default_rng(POWER_SEED)
    chosen = [geometry[c]["spacings"][CHOSEN] for c in CASES]

    def mean_of(name):
        return float(np.mean([g["readings"][name] for g in chosen]))

    def ratio_of(name):
        return max(g["variance_ratio_to_single"][name]["median"] for g in chosen)

    calibration_rms_2_25 = float(np.sqrt(np.mean(calibration_sd**2))) / math.sqrt(2)
    shape = calibration_sd / float(np.sqrt(np.mean(calibration_sd**2)))
    empirical_chosen = empirical["spacings"][CHOSEN]
    single_candidates = {
        "realized_frozen_run": realized["single_sd"],
        "calibration_rms": calibration_rms_2_25,
    }
    single_sd = max(single_candidates.values())
    ratios = {
        "primary": ratio_of("primary"),
        "centre_term_only": ratio_of("centre_term_only"),
        "haar_noise_part": max(
            g["haar_side_links_spread"]["noise_variance_ratio_to_single"]
            for g in chosen
        ),
    }
    haar_signal_sd = max(
        g["haar_side_links_spread"]["signal_sd_independent_of_n"] for g in chosen
    )
    candidates = {
        "primary": {
            "single_times_continuum_ratio": single_sd * math.sqrt(ratios["primary"]),
            "frozen_run_signal_swap": empirical_chosen["primary_swapped"]["sd"],
        },
        "centre_term_only": {
            "single_times_continuum_ratio": single_sd
            * math.sqrt(ratios["centre_term_only"]),
            "frozen_run_exact": empirical_chosen["centre_term_only"]["sd"],
        },
        "haar_side_links": {
            "continuum": math.sqrt(
                single_sd**2 * ratios["haar_noise_part"] + haar_signal_sd**2
            ),
            "frozen_run_signal_swap": empirical_chosen["haar_side_links_total_sd"],
        },
        "single": single_candidates,
    }
    design_sd = {name: max(values.values()) for name, values in candidates.items()}
    means = {
        "primary": mean_of("primary_minimum_norm"),
        "single": mean_of("single"),
        "centre_term_only": mean_of("centre_term_only"),
        "parity_flipped_side_links": mean_of("parity_flipped_side_links"),
        "haar_side_links": mean_of("haar_side_links_mean"),
    }
    means_frozen_run = {
        "centre_term_only_exact": empirical_chosen["centre_term_only"]["mean"],
        "single_exact": empirical_chosen["single"]["mean"],
    }

    def sds_at(name, level, kappa):
        """Per-replicate sd draws of one reading at a level."""
        scale = math.sqrt(2**25 / level)
        if name == "haar_side_links":
            noise_sq = design_sd[name] ** 2 - haar_signal_sd**2
            total = math.sqrt(kappa * max(noise_sq, 0.0) * scale**2 + haar_signal_sd**2)
            return shape * total
        return shape * design_sd[name] * scale * math.sqrt(kappa)

    result = {
        "spacing": CHOSEN,
        "inputs": {
            "means_continuum_cases_0_1": means,
            "means_frozen_run_reweighted": means_frozen_run,
            "variance_ratio_to_single": ratios,
            "haar_signal_sd_independent_of_n": haar_signal_sd,
            "sd_2_25_candidates": candidates,
            "sd_2_25_design": design_sd,
            "kappa": list(KAPPAS),
            "mean_offsets": list(OFFSETS_MEAN),
            "draws_per_power": DRAWS,
            "seed": POWER_SEED,
        },
        "g1_accuracy": [],
    }
    for level in LEVELS:
        for kappa in KAPPAS:
            sds = sds_at("primary", level, kappa)
            for offset in OFFSETS_MEAN:
                for tolerance in TOLERANCES:
                    powers = {
                        r: interval_power(
                            means["primary"] + offset, sds, r, "inside", rng, tolerance
                        )
                        for r in REPLICATES
                    }
                    needed = next(
                        (r for r in REPLICATES if powers[r] >= TARGET_POWER), None
                    )
                    half = next((r for r in REPLICATES if powers[r] >= 0.5), None)
                    result["g1_accuracy"].append(
                        {
                            "events_per_window": level,
                            "kappa": kappa,
                            "mean_offset": offset,
                            "mean_reading": means["primary"] + offset,
                            "tolerance": tolerance,
                            "power_by_replicates": {
                                str(k): v for k, v in powers.items()
                            },
                            "replicates_for_90pct": needed,
                            "replicates_for_50pct": half,
                        }
                    )
                    print(
                        "G1",
                        level,
                        kappa,
                        round(means["primary"] + offset, 4),
                        tolerance,
                        needed,
                        flush=True,
                    )
    # Proposed size at 2^25 and t = 0.05: the smallest grid R with 90% power in
    # every sizing scenario.
    needed = {}
    for kappa, offset in SIZING_SCENARIOS:
        row = next(
            row
            for row in result["g1_accuracy"]
            if row["events_per_window"] == 2**25
            and row["kappa"] == kappa
            and row["tolerance"] == 0.05
            and row["mean_offset"] == offset
        )
        needed[f"kappa {kappa}, mean offset {offset:+.2f}"] = row[
            "replicates_for_90pct"
        ]
    proposed = None if None in needed.values() else max(needed.values())
    result["proposed_design"] = {
        "events_per_window": 2**25,
        "tolerance": 0.05,
        "replicates": proposed,
        "replicates_for_90pct_by_scenario": needed,
        "sized_for": (
            "90% power of G1 in each scenario: kappa 1.2 at the continuum reading, "
            "and kappa 1 with the mean 0.01 above or below it"
        ),
    }
    sizes = sorted({48, 64, 80, 96, 112, 128} | ({proposed} if proposed else set()))
    others = {}
    for kappa in KAPPAS:
        block = {}
        for r in sizes:
            entry = {
                "g2_single_below": interval_power(
                    means["single"],
                    sds_at("single", 2**25, kappa),
                    r,
                    "below",
                    rng,
                    0.05,
                ),
                "f1_centre_term_only_outside": interval_power(
                    min(
                        means["centre_term_only"],
                        means_frozen_run["centre_term_only_exact"],
                    ),
                    sds_at("centre_term_only", 2**25, kappa),
                    r,
                    "outside",
                    rng,
                    0.05,
                ),
                "f2_haar_side_links_outside_one_draw": interval_power(
                    means["haar_side_links"],
                    sds_at("haar_side_links", 2**25, kappa),
                    r,
                    "outside",
                    rng,
                    0.05,
                ),
                "f2_parity_flipped_outside": interval_power(
                    means["parity_flipped_side_links"],
                    sds_at("primary", 2**25, kappa),
                    r,
                    "outside",
                    rng,
                    0.05,
                ),
            }
            for offset_key in ("1", "-1", "2", "-2"):
                rho = float(
                    np.mean(
                        [
                            g["window_reading_relative_to_centre_truth"][offset_key]
                            for g in chosen
                        ]
                    )
                )
                entry[f"g3_window_{offset_key}_sign"] = interval_power(
                    rho,
                    sds_at("single", 2**25, kappa),
                    r,
                    "positive" if rho > 0 else "negative",
                    rng,
                    0.05,
                )
            block[str(r)] = entry
        others[str(kappa)] = block
    result["other_gates_at_2_25"] = {
        "note": (
            "F2 for one Haar draw; G3 uses the single-window spread for each side "
            "window's linked reading."
        ),
        "by_kappa": others,
    }
    # G4 (secondary), the side share a - a0 against its band, drawn after every
    # other power above so that those draws are unchanged.
    share_means = {
        str(c): g["readings"]["side_share"] for c, g in zip(CASES, chosen, strict=True)
    }
    share_candidates = {
        "single_times_continuum_ratio": single_sd * math.sqrt(ratio_of("side_share")),
        "frozen_run_signal_swap": empirical_chosen["side_share_swapped"]["sd"],
    }
    share_sd = max(share_candidates.values())
    share_power = {}
    for kappa in KAPPAS:
        for case, mean in share_means.items():
            for offset in SIDE_SHARE_OFFSETS:
                share_power[f"kappa {kappa}, case {case}, offset {offset:+.2f}"] = {
                    str(r): interval_power(
                        mean + offset,
                        shape * share_sd * math.sqrt(kappa),
                        r,
                        "inside_band",
                        rng,
                        0.05,
                    )
                    for r in sizes
                }
    stress = {
        f"case {case}, offset {offset:+.2f}": interval_power(
            mean + offset, SIDE_SHARE_STRESS_SD, proposed, "inside_band", rng, 0.05
        )
        for case, mean in share_means.items()
        for offset in SIDE_SHARE_OFFSETS
    }
    result["g4_side_share_secondary"] = {
        "band": list(SIDE_SHARE_BAND),
        "means_continuum_by_case": share_means,
        "mean_frozen_run_reweighted": empirical_chosen["side_share_swapped"]["mean"],
        "variance_ratio_to_single": ratio_of("side_share"),
        "sd_2_25_candidates": share_candidates,
        "sd_2_25_design": share_sd,
        "power_by_scenario": share_power,
        "lowest_power_at_proposed": min(v[str(proposed)] for v in share_power.values()),
        "stress_sd": SIDE_SHARE_STRESS_SD,
        "stress_power_at_proposed": stress,
    }
    print(
        "G4",
        round(share_sd, 4),
        result["g4_side_share_secondary"]["lowest_power_at_proposed"],
        min(stress.values()),
        flush=True,
    )
    side_clocks = dry["summary"]["side"]["clocks"]
    per_replicate = {
        f"{level}_{label}_side_clocks": replicate_cost(level, side_clocks[label])
        for level in LEVELS
        for label in ("median", "max")
    }
    frozen_check = replicate_cost(2**25, SIDE_FIT["clocks"])
    totals = {}
    for level in LEVELS:
        per = replicate_cost(level, side_clocks["max"])
        for r in sizes:
            totals[f"{level}_{r}"] = {
                "events_per_window": level,
                "replicates": r,
                "measurement_events": 5 * level * r,
                "selection_events": 5 * SELECTION_N * r,
                "gpu_busy_hours": per["gpu_busy_seconds"] * r / 3600,
                "fit_cpu_hours": per["fit_cpu_seconds"] * r / 3600,
                "wall_hours_by_fit_workers": {
                    n: s * r / 3600
                    for n, s in per["wall_seconds_by_fit_workers"].items()
                },
            }
    result["cost"] = {
        "model": (
            "generation 96 s per window at 2^25 (frozen run median, independent of "
            "the landmark count there), proportional to N; side-window fit "
            "311 s x (clocks/6209)^1.63 with its clock-moment share (128 of 311 s) "
            "proportional to N; centre fit as measured (149 s, 85 s of it moments). "
            "Wall time per replicate: the larger of generation and fits / workers; "
            "the GPU is reserved for the whole wall time. The frozen run took 8.62 h "
            "for 64 replicates with 3 fit workers; this model gives "
            f"{64 * frozen_check['wall_seconds_by_fit_workers']['3'] / 3600:.2f} h."
        ),
        "side_clocks_from_dry_run": side_clocks,
        "per_replicate": per_replicate,
        "totals_with_max_side_clocks": totals,
    }
    return result


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--outputs",
        type=Path,
        default=OUTPUTS,
        help="folder holding the recorded runs named below (default: outputs/)",
    )
    for name, folder in (
        ("--archive", "intrinsic_weyl_bombelli_benchmark_v1"),
        ("--clock-cache", "intrinsic_weyl_window_study_v1/cache"),
        ("--window-study", "intrinsic_weyl_window_study_v1"),
        ("--calibration", "intrinsic_weyl_calibration_v1"),
        ("--frozen-run", "intrinsic_weyl_window_test_v1"),
        ("--engineering", "intrinsic_weyl_window_test_engineering_v2"),
    ):
        parser.add_argument(name, type=Path, help=f"default: <outputs>/{folder}")
    parser.add_argument(
        "--frozen-summary",
        type=Path,
        default=ROOT / "docs/paper/paper_a/intrinsic_weyl/window_test_summary.json",
    )
    parser.add_argument(
        "--frozen-power",
        type=Path,
        default=ROOT / "docs/paper/paper_a/intrinsic_weyl/window_test_power.json",
    )
    parser.add_argument("--fit-cache", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--threads", type=int, default=3)
    parser.add_argument("--skip-dry-run", action="store_true", help="debug only")
    args = parser.parse_args()
    for name, folder in (
        ("archive", "intrinsic_weyl_bombelli_benchmark_v1"),
        ("clock_cache", "intrinsic_weyl_window_study_v1/cache"),
        ("window_study", "intrinsic_weyl_window_study_v1"),
        ("calibration", "intrinsic_weyl_calibration_v1"),
        ("frozen_run", "intrinsic_weyl_window_test_v1"),
        ("engineering", "intrinsic_weyl_window_test_engineering_v2"),
    ):
        if getattr(args, name) is None:
            setattr(args, name, args.outputs / folder)
    sys.path.insert(0, str(HERE))
    from intrinsic_weyl_bombelli_benchmark import benchmark_amplitude

    commit, clean = git_state()
    began = time.time()
    amplitude = float(benchmark_amplitude())
    checks = selfcheck()
    print(json.dumps(checks), flush=True)
    fits, fit_record = continuum_fits(
        args.archive,
        args.clock_cache,
        amplitude,
        args.workers,
        args.threads,
        args.fit_cache,
    )
    agreement = window_study_agreement(fits, args.window_study)
    frozen_power = json.loads(args.frozen_power.read_text(encoding="utf-8"))
    geometry = {
        case: case_geometry(case, fits, args.clock_cache, amplitude, frozen_power)
        for case in CASES
    }
    print("geometry done", flush=True)

    def rho(label):
        return [
            float(
                np.mean(
                    [
                        geometry[c]["spacings"][label][
                            "window_reading_relative_to_centre_truth"
                        ][str(k)]
                        for c in CASES
                    ]
                )
            )
            for k in OFFSETS
        ]

    empirical = frozen_empirical(
        args.frozen_run,
        args.frozen_summary,
        rho(FROZEN),
        {label: rho(label) for label in SPACINGS},
    )
    print("empirical done", flush=True)
    dry = None if args.skip_dry_run else dry_run(args.workers, args.engineering)
    print("dry run done", flush=True)
    frozen_summary = json.loads(args.frozen_summary.read_text(encoding="utf-8"))
    level = frozen_summary["levels"][str(FROZEN_LEVEL)]["descriptive"]["readings"]
    realized = {
        "single_sd": level["single"]["full"]["sd"],
        "primary_sd": level["primary"]["full"]["sd"],
    }
    calibration_sd = calibration_spread(args.calibration, amplitude)
    power = (
        power_and_cost(geometry, empirical, dry, calibration_sd, realized)
        if dry is not None
        else None
    )
    result = {
        "design_aid": (
            "geometry, power and cost of a finite-N test in which the side windows and "
            "order-only links of the five-window correction carry the reading; no "
            "measurement event drawn; continuum, post hoc and selection-only inputs"
        ),
        "code_version": commit,
        "worktree_clean": clean,
        "source_sha256_lf": lf_sha(__file__),
        "benchmark_amplitude": amplitude,
        "wavelength_diamond_heights": WAVELENGTH,
        "spacings": SPACINGS,
        "chosen_spacing": CHOSEN,
        "selfcheck": checks,
        "fits": fit_record,
        "window_study_agreement": agreement,
        "geometry": {str(c): geometry[c] for c in CASES},
        "frozen_run_reweighted": empirical,
        "realized_frozen_run_sd_2_25": realized,
        "calibration_single_sd_2_24": {
            "selections_x_orientations": int(len(calibration_sd)),
            "rms": float(np.sqrt(np.mean(calibration_sd**2))),
        },
        "selection_dry_run": dry,
        "power_and_cost": power,
        "elapsed_seconds": time.time() - began,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
