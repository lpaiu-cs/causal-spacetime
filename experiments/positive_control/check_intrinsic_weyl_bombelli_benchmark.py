"""Independent post-completion audit of the Ibrahim-Bombelli benchmark run.

Nothing here calls the runner's judgment or evaluation code. The statistic,
threshold, verdicts, reflection, tidal distances, radii and pattern
signal-to-noise ratios are rebuilt from the saved pure fits through the 3x3
tidal tensors, and a sample of the archived order bits is recomputed with a
general-purpose ODE integrator that shares no code with the generator. The
truth vectors are reused from the run after their invariants are checked. Every
`posthoc` entry was added after the outcomes were known and decides nothing.
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
from check_intrinsic_weyl_weak_accuracy import quantile10, read, sha
from intrinsic_weyl import WEYL
from replay_intrinsic_weyl_threads import THREAD_VARIABLES, differences, write_json

ROOT = Path(__file__).resolve().parents[2]
LEVELS = (1048576, 4194304, 16777216, 67108864)
KINDS = ("wave", "wave", "flat", "flat")
SEED, H0, FALSE_POSITIVE = 2026111100, 0.01, 1e-3
EPSILON = np.zeros((3, 3, 3))
for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
    EPSILON[i, j, k], EPSILON[i, k, j] = 1, -1
# E_ij = C_0i0j and B_ij = eps_ikl C_kl0j / 2, taken from the full tensor basis.
ELECTRIC = WEYL[:, 0, 1:, 0, 1:]
MAGNETIC = 0.5 * np.einsum("ikl,nklj->nij", EPSILON, WEYL[:, 1:, 1:, 0, 1:])
GRAM2 = sum(np.einsum("nij,mij->nm", part, part) for part in (ELECTRIC, MAGNETIC))


def tensors(c):
    """Electric and magnetic 3x3 parts of a ten-coefficient vector."""
    return tuple(np.einsum("n,nij->ij", c, part) for part in (ELECTRIC, MAGNETIC))


def reflection(signs):
    """10x10 action of diag(1, signs) on the four-index tensor itself."""
    d = np.r_[1.0, signs]
    turned = WEYL * d[:, None, None, None] * d[:, None, None] * d[:, None] * d
    flat = WEYL.reshape(10, -1).T
    matrix = np.linalg.lstsq(flat, turned.reshape(10, -1).T, rcond=None)[0]
    if not np.allclose(flat @ matrix, turned.reshape(10, -1).T, atol=1e-13):
        raise ArithmeticError("The tensor basis is not closed under this reflection")
    return matrix


MIRROR = reflection([1, 1, -1])


def norm(c):
    return math.sqrt(float(c @ GRAM2 @ c))


def distance(a, b):
    """Full-tidal distance modulo the unobservable z reflection."""
    return min(norm(a - c) for c in (b, MIRROR @ b))


def survival10(x):
    return math.exp(-x / 2) * sum((x / 2) ** j / math.factorial(j) for j in range(5))


def statistic(fit):
    c = np.asarray(fit["coefficients"])
    lower = np.linalg.cholesky(np.asarray(fit["nominal_component_covariance"]))
    return float(np.sum(np.linalg.solve(lower, c) ** 2))


def framed(fit):
    """Coefficients and covariance with the fitted x and y landmark axes positive."""
    points = np.asarray(fit["coordinates_inferred_from_counts"])
    turn = reflection([np.sign(points[1, 1]), np.sign(points[2, 2]), 1])
    c = turn @ np.asarray(fit["coefficients"])
    s = turn @ np.asarray(fit["nominal_component_covariance"]) @ turn.T
    valid = (
        fit["optimizer_success"]
        and fit["boundary_excess"] <= 0
        and fit["coordinate_rank"] == fit["coordinate_parameters"]
        and fit["projected_curvature_rank"] == 10
        and fit["nuisance_rank"] == fit["nuisance_parameters"]
        and min(abs(points[1, 1]), abs(points[2, 2])) > 1e-12
        and np.isfinite(c).all()
        and np.linalg.eigvalsh(s)[0] > 0
    )
    return c, s, bool(valid)


def verdict(rows):
    """Frozen gate, restated: rows are (kind, statistic or None, valid)."""
    limit = quantile10(1 - FALSE_POSITIVE)
    if [kind for kind, _, _ in rows] != list(KINDS):
        raise ValueError("All four assigned cases are required, in order")
    if any(q is None or not valid for _, q, valid in rows):
        return "INCOMPLETE"
    wave = [q > limit for kind, q, _ in rows if kind == "wave"]
    flat = [q > limit for kind, q, _ in rows if kind == "flat"]
    if any(flat):
        return "CONTROL_FAILED"
    return "DETECTED" if all(wave) else "MIXED" if any(wave) else "NOT_DETECTED"


def patterns(c, s, truth, n):
    """Signal-to-noise along the truth's tidal patterns, less favorable mirror."""
    masks = {
        "electric": np.r_[np.ones(5), np.zeros(5)],
        "magnetic": np.r_[np.zeros(5), np.ones(5)],
        "full": np.ones(10),
    }
    candidates = (truth, MIRROR @ truth)
    matched = min(candidates, key=lambda t: norm(c - t))

    def read_out(target):
        weight = GRAM2 @ target / norm(target)
        return weight, math.sqrt(float(weight @ s @ weight))

    rows = {}
    for name, mask in masks.items():
        snr = min(norm(mask * t) / read_out(mask * t)[1] for t in candidates)
        weight, sigma = read_out(mask * matched)
        rows[name] = {
            "expected_snr": snr,
            "events_for_expected_snr_one": n / snr**2,
            "posthoc_matched_mirror_realized_z": float(weight @ c) / sigma,
        }
    return rows


def close(left, right, tolerance=1e-8):
    return abs(left - right) <= tolerance * max(1.0, abs(left), abs(right))


def audit_row(folder, case, n, recorded):
    pure = read(folder / f"case_{case}_n{n}_pure.json")
    full = read(folder / f"case_{case}_n{n}.json")
    archive = folder / f"case_{case}_n{n}.npz"
    fit, continuum = pure["fit"], full["full"]["continuum_fit"]
    detection, evaluation = recorded["detection"], recorded["evaluation"]
    c, s, valid = framed(fit)
    mean = framed(continuum)[0]
    q = statistic(fit)
    largest = float(np.linalg.eigvals(s @ GRAM2).real.max())
    radius = (math.sqrt(quantile10(0.95) * largest) + 0.02 * norm(c)) / 0.98
    truth = np.asarray(evaluation["truth_coefficients_xy_positive"])
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
        "continuum_bias_mod_reflection": distance(mean, truth),
        "continuum_tidal_norm": norm(mean),
        "estimate_tidal_norm": norm(c),
        "electric_eigenvalues": np.linalg.eigvalsh(tensors(c)[0]).tolist(),
        "magnetic_eigenvalues": np.linalg.eigvalsh(tensors(c)[1]).tolist(),
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
        close(
            row["continuum_bias_mod_reflection"],
            evaluation["continuum_bias_mod_reflection"],
        ),
        close(row["continuum_tidal_norm"], evaluation["continuum_tidal_norm"]),
        bool(row["error_mod_reflection"] <= radius)
        == evaluation["nominal_ball_contains_truth"],
    ]
    if KINDS[case] == "wave":
        row["patterns"] = patterns(c, s, truth, n)
        agree += [
            close(row["patterns"][name][key], evaluation["patterns"][name][key])
            for name in row["patterns"]
            for key in ("expected_snr", "events_for_expected_snr_one")
        ]
    row["recorded_values_reproduced"] = all(agree)
    # Post hoc, descriptive. least_squares.cost is half the whitened chi-square.
    noncentrality = statistic(continuum)
    nearest = min((mean, MIRROR @ mean), key=lambda m: norm(c - m))
    row["posthoc"] = {
        "fit_chi2": 2 * fit["curved_cost"],
        "fit_degrees_of_freedom": fit["measured_intervals"]
        - fit["nuisance_parameters"]
        - 10,
        "flat_minus_curved_chi2": 2 * (min(fit["flat_costs"]) - fit["curved_cost"]),
        "added_curvature_parameters": fit["nuisance_parameters"]
        - fit["coordinate_parameters"]
        + 10,
        "continuum_noncentrality": noncentrality,
        "expected_q": 10 + noncentrality,
        "expected_q_sd": math.sqrt(2 * (10 + 2 * noncentrality)),
        "q_about_continuum_refit": float(
            (c - nearest) @ np.linalg.solve(s, c - nearest)
        ),
        "continuum_refit_chi2": 2 * continuum["curved_cost"],
    }
    return row


def independent_margin(first, last, amplitude):
    """Causal margin from a direct ODE solution of both transverse Jacobi equations.

    Along u the extremal of (x'^2 + a x^2)/2 has action [x x']/2, and a pair is
    related iff the drop in v exceeds the summed action of x'' = a x, y'' = -a y.
    """
    from scipy.integrate import solve_ivp

    def rhs(u, y):
        a = amplitude * math.sin(math.pi * u)
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


def audit_case(folder, case, pairs=1000):
    """Nested prefixes, a CPU replay, and integrator-independent order bits."""
    from profile_plane_wave import constant, half_cycle

    kind = KINDS[case]
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
        points, amplitude = file["points"], str(file["amplitude"].item())
        labelled = str(file["kind"].item()) == kind
    value = float(amplitude)
    expected = math.pi**2 * H0 / 2 if kind == "wave" else 0.0
    wave = (half_cycle if kind == "wave" else constant)(amplitude)
    events = wave.sample(4096, np.random.default_rng([SEED + case, 1]))
    stored = [
        np.unpackbits(top[key][:, :512], axis=1).astype(bool)
        for key in ("after", "before")
    ]
    replayed = wave.relations(points[:, None], events[None])
    left, right = points[:, None], events[None]
    du = right[..., 0] - left[..., 0]
    span = np.where(du == 0, 1.0, abs(du))
    flat = -np.sign(du) * (right[..., 1] - left[..., 1]) - (
        (right[..., 2:] - left[..., 2:]) ** 2
    ).sum(axis=-1) / (2 * span)
    report = {
        "case": case,
        "kind": kind,
        "landmarks_including_tips": len(points),
        "truth_file_labelled": bool(labelled and close(value, expected, 1e-14)),
        "nested_prefixes_bit_exact": bool(nested),
        "tips_see_every_event": bool(
            (top["after"][0] == 255).all() and (top["before"][-1] == 255).all()
        ),
        "cpu_replay_of_first_4096_events_bit_exact": all(
            np.array_equal(a, b) for a, b in zip(replayed, stored, strict=True)
        ),
        "first_chunk_pairs_decided_differently_from_flat": int(
            ((stored[0] | stored[1]) != ((flat > 0) & (du != 0))).sum()
        ),
    }
    if kind == "flat":
        clear = (abs(flat) > 1e-12) & (du != 0)
        related = (flat > 0) & clear
        report["closed_form_pairs_checked"] = int(clear.sum())
        report["closed_form_bits_match"] = bool(
            np.array_equal(stored[0][clear], (related & (du > 0))[clear])
            and np.array_equal(stored[1][clear], (related & (du < 0))[clear])
        )
        return report
    # The pairs nearest the flat cone are the ones curvature can decide.
    order = np.argsort(abs(flat), axis=None)
    chosen = np.unique(
        np.r_[
            order[:pairs],
            np.random.default_rng([SEED + case, 2]).choice(flat.size, pairs, False),
        ]
    )
    matched, smallest, curved = 0, math.inf, 0
    for landmark, event in zip(*np.unravel_index(chosen, flat.shape), strict=True):
        ahead = du[landmark, event] > 0
        first, last = (
            (points[landmark], events[event])
            if ahead
            else (events[event], points[landmark])
        )
        margin = independent_margin(first, last, value)
        smallest = min(smallest, abs(margin))
        bits = (stored[0][landmark, event], stored[1][landmark, event])
        matched += bits == ((margin > 0) and ahead, (margin > 0) and not ahead)
        curved += (margin > 0) != (flat[landmark, event] > 0)
    report.update(
        independent_pairs_checked=len(chosen),
        independent_bits_match=int(matched),
        independent_pairs_decided_differently_from_flat=int(curved),
        smallest_independent_margin=smallest,
    )
    return report


def truth_invariants(truth):
    """A plane-wave tidal pair has eigenvalues (-A0, 0, A0) twice and EB + BE = 0."""
    e, b = tensors(np.asarray(truth))
    target = math.pi**2 * H0 / 2 * np.array([-1.0, 0.0, 1.0])
    return {
        "electric_eigenvalue_error": float(abs(np.linalg.eigvalsh(e) - target).max()),
        "magnetic_eigenvalue_error": float(abs(np.linalg.eigvalsh(b) - target).max()),
        "anticommutator_norm": float(np.linalg.norm(e @ b + b @ e)),
    }


def album_density(total, thickness, exponent=0.11, width=2.0, separation=0.5):
    """rho H^4 at which the box of arXiv:2611.09056, App. B, holds `total` elements.

    The box is (a_w d)^2 (L + a_w d) H with L = a_L lambda = 2 a_L H, and the
    tube thickness is d/H = A (pi rho H^4 / 24)^-beta. Defaults are the paper's
    most liberal a_w, a_L and its four-dimensional exponent.
    """
    from scipy.optimize import brentq

    def excess(log_density):
        density = 10**log_density
        d = thickness * (math.pi * density / 24) ** -exponent
        return density * (width * d) ** 2 * (2 * separation + width * d) - total

    return 10 ** brentq(excess, 6, 16)


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
            str(Path(__file__).with_name("intrinsic_weyl_bombelli_benchmark.py")),
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
    from scipy.optimize import brentq
    from scipy.stats import chi2, ncx2

    limit = quantile10(1 - FALSE_POSITIVE)
    assert abs(limit - chi2.ppf(1 - FALSE_POSITIVE, 10)) < 1e-10
    # Coefficients 0-4 are purely electric and 5-9 purely magnetic in this basis.
    assert not ELECTRIC[5:].any() and not MAGNETIC[:5].any()
    for signs in ([1, 1, -1], [-1, 1, 1], [1, -1, 1], [-1, -1, -1]):
        turn = reflection(signs)
        assert np.allclose(turn @ turn, np.eye(10))
        assert np.allclose(turn.T @ GRAM2 @ turn, GRAM2)
    summary = read(folder / "summary.json")
    stamp = summary["code_version"]

    def git(*arguments):
        return subprocess.run(
            ["git", *arguments], cwd=ROOT, capture_output=True, text=True
        )

    recorded = {(r["case"], r["measurement_n"]): r for r in summary["results"]}
    rows = [
        audit_row(folder, case, n, recorded[case, n])
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
    cases = [audit_case(folder, case) for case in range(4)]
    invariants = {
        case: truth_invariants(
            recorded[case, LEVELS[-1]]["evaluation"]["truth_coefficients_xy_positive"]
        )
        for case in (0, 1)
    }
    power = {}
    for case in (0, 1):
        top = next(
            r for r in rows if r["case"] == case and r["measurement_n"] == LEVELS[-1]
        )
        rate = top["posthoc"]["continuum_noncentrality"] / LEVELS[-1]
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
        "rows": rows,
        "level_verdicts": verdicts,
        "level_verdicts_match_run": verdicts == summary["level_verdicts"],
        "cases": cases,
        "truth_invariants": invariants,
        "posthoc_power_from_continuum_refit": power,
        # A diamond of height H holds pi rho H^4 / 24 events.
        "posthoc_density_rho_H4": {
            "this_run_by_level": {str(n): 24 * n / math.pi for n in LEVELS},
            "ibrahim_bombelli_eq20": {
                f"m={m},A={a}": album_density(4.5e11 / m**3.5, a)
                for m in (11, 5, 4)
                for a in (0.3, 0.32)
            },
        },
        "refits": [refit(folder, case, n) for case, n in refits],
    }
    result["passed"] = bool(
        result["run_stamp_is_ancestor_of_head"]
        and result["level_verdicts_match_run"]
        and all(
            r["valid_fit"]
            and r["recorded_values_reproduced"]
            and r["blind_boundary_kept"]
            and r["archive_hash_matches"]
            for r in rows
        )
        and all(
            c["truth_file_labelled"]
            and c["nested_prefixes_bit_exact"]
            and c["tips_see_every_event"]
            and c["cpu_replay_of_first_4096_events_bit_exact"]
            and c.get("closed_form_bits_match", True)
            and c.get("independent_bits_match", 0)
            == c.get("independent_pairs_checked", 0)
            for c in cases
        )
        and all(
            max(v["electric_eigenvalue_error"], v["magnetic_eigenvalue_error"]) < 1e-12
            and v["anticommutator_norm"] < 1e-12
            for v in invariants.values()
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
        print(json.dumps({k: v for k, v in line.items() if k != "patterns"}))
    for line in outcome["cases"]:
        print(json.dumps(line))
    print(
        json.dumps(
            {
                key: outcome[key]
                for key in (
                    "level_verdicts",
                    "level_verdicts_match_run",
                    "truth_invariants",
                    "posthoc_power_from_continuum_refit",
                    "posthoc_density_rho_H4",
                    "refits",
                    "passed",
                )
            }
        )
    )
    sys.exit(0 if outcome["passed"] else 1)
