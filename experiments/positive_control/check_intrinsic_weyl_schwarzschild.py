"""Independent post-completion audit of the Schwarzschild detection test.

Nothing here calls the runner's judgment or evaluation code. The statistic,
verdicts, radii and the truth comparison are rebuilt from the saved pure fits
through the 3x3 tidal tensors. The truth tensor comes from the covariant Hessian
of the lapse, not from the generator's closed form, and its frame from an O(3)
alignment of all landmarks instead of the two-landmark rule. A sample of the
archived order bits is recomputed with the older Clairaut solver in areal
coordinates, which shares no code with the generator, and the first measurement
events are tested against the proper-volume measure with that solver and the
areal volume element. Every `posthoc` entry was added after the outcomes were
known and decides nothing.
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
    MIRROR,
    close,
    distance,
    framed,
    norm,
    statistic,
    survival10,
    tensors,
)
from check_intrinsic_weyl_weak_accuracy import quantile10, read, sha
from replay_intrinsic_weyl_threads import THREAD_VARIABLES, differences, write_json

ROOT = Path(__file__).resolve().parents[2]
LEVELS = (1048576, 4194304, 16777216, 67108864)
KINDS = ("curved", "curved", "flat", "flat")
SEED, FALSE_POSITIVE = 2026111200, 1e-3
MASS, APEX, PROPER_TIME = 1.0, 6.0, 3.0
ELECTRIC_GRAM = np.einsum("kij,lij->kl", ELECTRIC[:5], ELECTRIC[:5])


def areal(rho, mass):
    return rho * (1 + mass / (2 * rho)) ** 2


def hessian_tidal(mass, step=1e-2):
    """T^2 E at the apex from E_ij = (D_i D_j A) / A for a static metric.

    With spatial metric B^4 delta and phi = ln B, the radial component is
    (A'' - 2 phi' A') / (A B^4) and a transverse one (A'/rho + 2 phi' A') / (A B^4).
    Derivatives are fourth-order central differences in the isotropic radius.
    """
    m = mass / 2
    rho = (APEX - mass + math.sqrt(APEX * APEX - 2 * mass * APEX)) / 2

    def lapse(x):
        return (1 - m / x) / (1 + m / x)

    def phi(x):
        return math.log(1 + m / x)

    def d1(f):
        return (
            f(rho - 2 * step)
            - 8 * f(rho - step)
            + 8 * f(rho + step)
            - f(rho + 2 * step)
        ) / (12 * step)

    first, slope = d1(lapse), d1(phi)
    second = (
        -lapse(rho - 2 * step)
        + 16 * lapse(rho - step)
        - 30 * lapse(rho)
        + 16 * lapse(rho + step)
        - lapse(rho + 2 * step)
    ) / (12 * step**2)
    scale = PROPER_TIME**2 / (lapse(rho) * (1 + m / rho) ** 4)
    radial = (second - 2 * slope * first) * scale
    transverse = (first / rho + 2 * slope * first) * scale
    return np.diag([radial, transverse, transverse])


def electric_coefficients(tensor):
    return np.r_[
        np.linalg.solve(ELECTRIC_GRAM, np.einsum("kij,ij->k", ELECTRIC[:5], tensor)),
        np.zeros(5),
    ]


def aligned_truth(fit, normal, tidal):
    """Truth in the fitted frame by O(3) alignment of all interior landmarks."""
    points = np.array(fit["coordinates_inferred_from_counts"])
    points[:, 1] *= np.sign(points[1, 1])
    points[:, 2] *= np.sign(points[2, 2])
    left, _, right = np.linalg.svd(normal[1:-1, 1:].T @ points[1:-1, 1:])
    rotation = left @ right
    residual = np.c_[normal[:, :1], normal[:, 1:] @ rotation] - points
    return (
        electric_coefficients(rotation.T @ tidal @ rotation),
        float(np.sqrt((residual**2).sum(axis=1).mean())),
    )


def verdict(rows):
    """Frozen gate, restated: rows are (kind, statistic or None, valid)."""
    limit = quantile10(1 - FALSE_POSITIVE)
    if [kind for kind, _, _ in rows] != list(KINDS):
        raise ValueError("All four assigned cases are required, in order")
    if any(q is None or not valid for _, q, valid in rows):
        return "INCOMPLETE"
    if any(q > limit for kind, q, _ in rows if kind == "flat"):
        return "CONTROL_FAILED"
    curved = [q > limit for kind, q, _ in rows if kind == "curved"]
    return "DETECTED" if all(curved) else "MIXED" if any(curved) else "NOT_DETECTED"


def pattern_snr(s, truth):
    """Nominal signal-to-noise along the truth pattern, less favorable mirror."""

    def one(target):
        weight = GRAM2 @ target / norm(target)
        return norm(target) / math.sqrt(float(weight @ s @ weight))

    return min(one(truth), one(MIRROR @ truth))


def audit_row(folder, case, n, recorded, normal, tidal):
    pure = read(folder / f"case_{case}_n{n}_pure.json")
    full = read(folder / f"case_{case}_n{n}.json")
    archive = folder / f"case_{case}_n{n}.npz"
    fit = pure["fit"]
    detection, evaluation = recorded["detection"], recorded["evaluation"]
    c, s, valid = framed(fit)
    q = statistic(fit)
    largest = float(np.linalg.eigvals(s @ GRAM2).real.max())
    radius = (math.sqrt(quantile10(0.95) * largest) + 0.02 * norm(c)) / 0.98
    truth = np.asarray(evaluation["truth_coefficients_xy_positive"])
    electric, magnetic = tensors(c)
    independent, residual = aligned_truth(fit, normal, tidal)
    size = float(np.linalg.norm(tidal))
    row = {
        "case": case,
        "kind": KINDS[case],
        "measurement_n": n,
        "valid_fit": valid,
        "q": q,
        "nominal_p_value": survival10(q),
        "variance_times_N": n * largest,
        "radius_95": radius,
        "error_mod_reflection": distance(c, truth),
        "error_with_independent_truth": distance(c, independent),
        "truth_vectors_differ_by": distance(truth, independent),
        "landmark_coordinate_rms_residual": residual,
        "estimate_tidal_norm": norm(c),
        "electric_eigenvalues": np.linalg.eigvalsh(electric).tolist(),
        "magnetic_norm": float(np.linalg.norm(magnetic)),
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
        close(n * largest, detection["variance_times_N"]),
        close(radius, evaluation["confidence_radius_95"]),
        close(
            row["error_mod_reflection"], evaluation["full_tidal_error_mod_reflection"]
        ),
        bool(row["error_mod_reflection"] <= radius)
        == evaluation["nominal_ball_contains_truth"],
        np.allclose(row["electric_eigenvalues"], detection["electric_eigenvalues"]),
        close(row["magnetic_norm"], detection["magnetic_norm"]),
        close(norm(truth), size) and close(evaluation["truth_norm"], size, 1e-6),
        # A wrong frame would move the truth by the order of its norm.
        row["truth_vectors_differ_by"] <= 0.1 * size,
    ]
    if size:
        snr = pattern_snr(s, truth)
        row["truth_pattern"] = {
            "expected_snr": snr,
            "events_for_expected_snr_one": n / snr**2,
        }
        row["relative_error"] = row["error_mod_reflection"] / size
        row["relative_radius"] = radius / size
        agree += [
            close(snr, evaluation["truth_pattern"]["expected_snr"]),
            close(row["relative_error"], evaluation["relative_full_tidal_error"]),
        ]
    row["recorded_values_reproduced"] = all(bool(v) for v in agree)
    # Post hoc, descriptive. least_squares.cost is half the whitened chi-square.
    freedom = fit["measured_intervals"] - fit["nuisance_parameters"] - 10
    row["posthoc"] = {
        "fit_chi2": 2 * fit["curved_cost"],
        "fit_degrees_of_freedom": freedom,
        "fit_chi2_z": (2 * fit["curved_cost"] - freedom) / math.sqrt(2 * freedom),
        "flat_minus_curved_chi2": 2 * (min(fit["flat_costs"]) - fit["curved_cost"]),
        "noncentrality_of_the_truth": min(
            float(t @ np.linalg.solve(s, t)) for t in (truth, MIRROR @ truth)
        ),
        "q_about_the_truth": min(
            float((c - t) @ np.linalg.solve(s, c - t)) for t in (truth, MIRROR @ truth)
        ),
    }
    return row


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


def measure_check(geometry, mass, events):
    """Moments of the events against the areal volume element and Clairaut heights."""
    offset, measure, height = geometry.interval_nodes()
    position = geometry.base + offset
    rho = np.linalg.norm(position, axis=2)
    own = np.array(
        [
            geometry.half_time - independent_time(geometry.base, x, mass)[0]
            for x in position.reshape(-1, 3)
        ]
    ).reshape(height.shape)
    # sqrt(-g) = r^2 sin(theta) in areal coordinates, pulled back to isotropic ones.
    density = (areal(rho, mass) / rho) ** 2 * (1 - (mass / (2 * rho)) ** 2)
    weight = measure * density * own
    shift = events[:, 1:] - geometry.base
    pairs = (
        (offset[..., 0], shift[:, 0]),
        ((offset**2).sum(axis=2), (shift**2).sum(axis=1)),
        (own**2 / 3, events[:, 0] ** 2),
        (rho, np.linalg.norm(events[:, 1:], axis=1)),
    )
    return {
        "events": len(events),
        "largest_height_difference_from_generator": float(abs(own - height).max()),
        "proper_volume": float(2 * weight.sum()),
        "proper_volume_over_flat": float(
            2 * weight.sum() / (math.pi * PROPER_TIME**4 / 24)
        ),
        "moment_z_scores": [
            float(
                (sample.mean() - (weight * value).sum() / weight.sum())
                / (sample.std() / math.sqrt(len(sample)))
            )
            for value, sample in pairs
        ],
    }


def audit_case(folder, case, normal, pairs=1000):
    """Nested prefixes, a CPU replay, independent order bits and the measure."""
    from schwarzschild_diamond import SchwarzschildDiamond

    kind = KINDS[case]
    mass = MASS if kind == "curved" else 0.0
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
    geometry = SchwarzschildDiamond(mass, APEX, PROPER_TIME)
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


def refit(folder, case, n):
    """Re-estimate one row in a fresh single-threaded process; fields must agree."""
    target = folder.with_name(folder.name + "_audit")
    target.mkdir(exist_ok=True)
    output = target / f"case_{case}_n{n}_{int(time.time())}.json"
    env = {k: v for k, v in os.environ.items() if k not in THREAD_VARIABLES}
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    subprocess.run(
        [
            sys.executable,
            str(Path(__file__).with_name("intrinsic_weyl_schwarzschild.py")),
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
    return {
        "case": case,
        "measurement_n": n,
        "q": statistic(new),
        "fields_identical": new == old,
        "differences_beyond_1e-8": differences(new, old),
    }


def audit(folder, receipt, refits):
    from schwarzschild_diamond import SchwarzschildDiamond
    from scipy.optimize import brentq
    from scipy.stats import chi2, ncx2

    limit = quantile10(1 - FALSE_POSITIVE)
    assert abs(limit - chi2.ppf(1 - FALSE_POSITIVE, 10)) < 1e-10
    summary = read(folder / "summary.json")
    stamp = summary["code_version"]

    def git(*arguments):
        return subprocess.run(
            ["git", *arguments], cwd=ROOT, capture_output=True, text=True
        )

    tidal = {"curved": hessian_tidal(MASS), "flat": hessian_tidal(0.0)}
    closed = PROPER_TIME**2 * MASS / APEX**3 * np.diag([-2.0, 1.0, 1.0])
    truth = {
        "hessian_tidal_eigenvalues": np.diag(tidal["curved"]).tolist(),
        "difference_from_closed_form": float(abs(tidal["curved"] - closed).max()),
        "trace": float(np.trace(tidal["curved"])),
        "flat_tensor_norm": float(np.linalg.norm(tidal["flat"])),
    }
    normal = {}
    for case, kind in enumerate(KINDS):
        with np.load(folder / f"truth_{case}.npz", allow_pickle=False) as file:
            normal[case] = SchwarzschildDiamond(
                MASS if kind == "curved" else 0.0, APEX, PROPER_TIME
            ).normal_coordinates(file["points"])
    recorded = {(r["case"], r["measurement_n"]): r for r in summary["results"]}
    rows = [
        audit_row(folder, case, n, recorded[case, n], normal[case], tidal[KINDS[case]])
        for case in range(4)
        for n in LEVELS
    ]
    verdicts = {
        str(n): verdict(
            [
                (r["kind"], r["q"], r["valid_fit"])
                for r in rows
                if r["measurement_n"] == n
            ]
        )
        for n in LEVELS
    }
    cases = [audit_case(folder, case, normal[case]) for case in range(4)]
    power = {}
    for case in (0, 1):
        top = next(
            r for r in rows if r["case"] == case and r["measurement_n"] == LEVELS[-1]
        )
        rate = top["posthoc"]["noncentrality_of_the_truth"] / LEVELS[-1]
        power[case] = {
            "noncentrality_per_event": rate,
            **{
                f"events_for_power_{int(100 * p)}": brentq(
                    lambda n, p=p, rate=rate: ncx2.sf(limit, 10, rate * n) - p,
                    1e4,
                    1e10,
                )
                for p in (0.5, 0.9)
            },
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
        "threshold_chi2_10": limit,
        "truth": truth,
        "rows": rows,
        "level_verdicts": verdicts,
        "level_verdicts_match_run": verdicts == summary["level_verdicts"],
        "cases": cases,
        # Unbiased-estimate assumption: the model bias is not measured in this run.
        "posthoc_power_if_unbiased": power,
        "refits": [refit(folder, case, n) for case, n in refits],
    }
    result["passed"] = bool(
        result["run_stamp_is_ancestor_of_head"]
        and result["level_verdicts_match_run"]
        and truth["difference_from_closed_form"] < 1e-7
        and abs(truth["trace"]) < 1e-7
        and truth["flat_tensor_norm"] == 0
        and all(
            r["valid_fit"]
            and r["recorded_values_reproduced"]
            and r["blind_boundary_kept"]
            and r["archive_hash_matches"]
            for r in rows
        )
        and all(
            c["truth_file_labelled"]
            and c["tips_at_half_unit_time_in_normal_coordinates"] < 1e-9
            and c["nested_prefixes_bit_exact"]
            and c["tips_see_every_event"]
            and c["cpu_replay_of_first_4096_events_bit_exact"]
            and c.get("closed_form_bits_match", True)
            and c.get("independent_bits_match", 0)
            == c.get("independent_pairs_checked", 0)
            - c.get("independent_pairs_inside_the_error_bound", 0)
            and c["measure"]["largest_height_difference_from_generator"] < 1e-7
            and max(map(abs, c["measure"]["moment_z_scores"])) < 4.5
            for c in cases
        )
        and all(r["fields_identical"] for r in result["refits"])
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
                    "posthoc_power_if_unbiased",
                    "refits",
                    "passed",
                )
            }
        )
    )
    sys.exit(0 if outcome["passed"] else 1)
