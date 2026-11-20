"""Independent post-completion check of the covariance calibration stage.

Nothing here calls the runner's verdict or evaluation code. The statistic is
rebuilt from the saved pure fits by a Cholesky solve, the chi-square quantiles
of the interval by bisection on the regularized incomplete gamma function
(cross-checked against Wilson-Hilferty), and the truth of the secondary
constant-wave population from the closed-form inverse exponential map and a
Gram-Schmidt frame on the first two landmarks, instead of the generator's map
and the run's O(3) alignment. Archives are hashed, nested prefixes compared bit
for bit, seeds and amplitudes compared with the frozen protocol, and the bound
sources with the commit that ran. Every `posthoc` entry was added after the
outcomes were known and decides nothing.
"""

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path
from statistics import NormalDist

import numpy as np
from check_intrinsic_weyl_bombelli_benchmark import (
    ELECTRIC,
    GRAM2,
    MAGNETIC,
    MIRROR,
    close,
    framed,
    statistic,
    survival10,
)
from check_intrinsic_weyl_weak_accuracy import read, sha
from replay_intrinsic_weyl_threads import THREAD_VARIABLES, differences, write_json

ROOT = Path(__file__).resolve().parents[2]
# Restated from docs/paper/paper_a/intrinsic_weyl/calibration_protocol.md.
LEVELS = (1048576, 16777216)
POPULATIONS = (("flat", 64, 2026112100), ("constant", 16, 2026113100))
ORDERS = (5, 6)
LOW, HIGH = 0.8, 1.25
COMPLETE_FRACTION = 60 / 64
H0 = 0.01


def chi2_quantile(probability, dof):
    """Bisection on the regularized lower incomplete gamma; no scipy.stats."""
    from scipy.special import gammainc

    lo, hi = 0.0, dof + 50 * math.sqrt(2 * dof) + 100
    for _ in range(200):
        mid = (lo + hi) / 2
        if gammainc(dof / 2, mid / 2) < probability:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def wilson_hilferty(probability, dof):
    z = NormalDist().inv_cdf(probability)
    return dof * (1 - 2 / (9 * dof) + z * math.sqrt(2 / (9 * dof))) ** 3


def judge(values, planned):
    """The frozen verdict for one jet and level, restated from the protocol."""
    values = np.asarray(values, dtype=float)
    out = {"completed": int(len(values)), "planned": planned}
    if len(values) < COMPLETE_FRACTION * planned:
        out["verdict"] = "INCOMPLETE"
        return out
    total, dof = float(values.sum()), 10 * len(values)
    upper, lower = chi2_quantile(0.975, dof), chi2_quantile(0.025, dof)
    interval = [total / upper, total / lower]
    out.update(
        kappa=total / dof,
        interval_95=interval,
        quantiles_vs_wilson_hilferty=max(
            abs(upper / wilson_hilferty(0.975, dof) - 1),
            abs(lower / wilson_hilferty(0.025, dof) - 1),
        ),
        verdict="CALIBRATED"
        if LOW <= interval[0] and interval[1] <= HIGH
        else "NOT_CALIBRATED",
        direction="UNDERSTATED"
        if interval[0] > 1
        else "OVERSTATED"
        if interval[1] < 1
        else "CONSISTENT_WITH_ONE",
        exceedances_of_q_star=int((values > chi2_quantile(0.999, 10)).sum()),
    )
    return out


def constant_truth(points, amplitude):
    """Ten coefficients of the constant wave in the gauge frame of the truth.

    In the midpoint frame (t, x, y, z) of the Brinkmann chart, R_uxux = -A and
    R_uyuy = A give T^2 E = diag(-A, A, 0) and T^2 B_xy = T^2 B_yx = -A, with
    E_ij = C_0i0j and B_ij = eps_ikl C_kl0j / 2.
    """
    from intrinsic_weyl_joint_audit import normal_coordinates

    normal = normal_coordinates(points, amplitude)
    first = normal[1, 1:] / np.linalg.norm(normal[1, 1:])
    second = normal[2, 1:] - (normal[2, 1:] @ first) * first
    second /= np.linalg.norm(second)
    axes = np.vstack((first, second, np.cross(first, second)))
    electric = np.diag([-amplitude, amplitude, 0.0])
    magnetic = np.zeros((3, 3))
    magnetic[0, 1] = magnetic[1, 0] = -amplitude
    e, b = axes @ electric @ axes.T, axes @ magnetic @ axes.T
    rhs = np.einsum("nij,ij->n", ELECTRIC, e) + np.einsum("nij,ij->n", MAGNETIC, b)
    c = np.linalg.solve(GRAM2, rhs)
    residual = max(
        abs(np.einsum("n,nij->ij", c, ELECTRIC) - e).max(),
        abs(np.einsum("n,nij->ij", c, MAGNETIC) - b).max(),
    )
    return c, float(residual)


def check_row(folder, replicate, kind, n, recorded, truth):
    pure = read(folder / f"case_{replicate}_n{n}_pure.json")
    archive = folder / f"case_{replicate}_n{n}.npz"
    row = {
        "replicate": replicate,
        "kind": kind,
        "measurement_n": n,
        "blind_boundary_kept": pure["oracle_modules_before_return"] == [],
        "archive_hash_matches": sha(archive)
        == recorded["archive_sha256"]
        == pure["archive_sha256"],
        "measurement_n_matches": pure["measurement_n"] == n,
    }
    agree = []
    for order in ORDERS:
        name = f"order_{order}"
        fit = pure["fits"][name]
        c, s, valid = framed(fit)
        q = statistic(fit)
        candidates = [np.zeros(10)] if truth is None else [truth, MIRROR @ truth]
        values = [float((c - t) @ np.linalg.solve(s, c - t)) for t in candidates]
        row[name] = {"valid_fit": valid, "q": q, "truth_statistic": min(values)}
        evaluation = recorded["evaluation"][name]
        agree += [
            valid,
            close(q, evaluation["q"]),
            close(min(values), evaluation["truth_statistic"], 1e-7),
        ]
    row["recorded_values_reproduced"] = all(bool(v) for v in agree)
    return row


def nested(folder, replicate):
    """The lower level is a bit-exact prefix of the upper one."""
    low, high = LEVELS
    with (
        np.load(folder / f"case_{replicate}_n{low}.npz", allow_pickle=False) as a,
        np.load(folder / f"case_{replicate}_n{high}.npz", allow_pickle=False) as b,
    ):
        same = set(a.files) == set(b.files)
        same &= int(a["measurement_n"]) == low and int(b["measurement_n"]) == high
        for key in ("after", "before"):
            same &= np.array_equal(a[key], b[key][:, : low // 8])
        for key in ("selection_order", "selection_design", "selected_ids"):
            same &= np.array_equal(a[key], b[key])
    return bool(same)


def refit(folder, replicate, n):
    """Re-estimate one row in a fresh single-threaded process; fields must agree."""
    target = folder.with_name(folder.name + "_audit")
    target.mkdir(exist_ok=True)
    output = target / f"case_{replicate}_n{n}_{int(time.time())}.json"
    env = {k: v for k, v in os.environ.items() if k not in THREAD_VARIABLES}
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    subprocess.run(
        [
            sys.executable,
            str(Path(__file__).with_name("intrinsic_weyl_calibration.py")),
            "--worker",
            str(folder / f"case_{replicate}_n{n}.npz"),
            "--truth",
            str(folder / f"truth_{replicate}.npz"),
            "--output",
            str(output),
        ],
        check=True,
        env=env,
    )
    new = read(output.with_name(output.stem + "_pure.json"))["fits"]
    old = read(folder / f"case_{replicate}_n{n}_pure.json")["fits"]
    for fits in (new, old):
        for fit in fits.values():
            fit.pop("fit_seconds", None)
    found = differences(new, old)
    return {
        "replicate": replicate,
        "measurement_n": n,
        "fields_identical": new == old,
        "differences_beyond_1e-8": found,
    }


def posthoc(values, kappa):
    """Descriptive: KS distance of Q / kappa to chi2_10, mean, median."""
    x = np.sort(np.asarray(values) / kappa)
    cdf = np.array([1 - survival10(v) for v in x])
    k = np.arange(1, len(x) + 1)
    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "ks_distance_scaled": float(
            max((k / len(x) - cdf).max(), (cdf - (k - 1) / len(x)).max())
        ),
    }


def audit(folder, receipt, refits):
    summary = read(folder / "summary.json")
    manifest = read(folder / "manifest.json")
    stamp = summary["code_version"]

    def git(*arguments):
        return subprocess.run(
            ["git", *arguments], cwd=ROOT, capture_output=True, text=False
        )

    sources = {}
    for name, digest in manifest["source_sha256_lf"].items():
        blob = git("show", f"{stamp}:{name}").stdout.replace(b"\r\n", b"\n")
        sources[name] = hashlib.sha256(blob).hexdigest() == digest
    amplitude = math.pi**2 * H0 / 2
    expected = []
    for kind, count, seed in POPULATIONS:
        expected += [(kind, seed + r) for r in range(count)]
    assigned = [(a["kind"], a["seed"]) for a in manifest["assignments"]]
    amplitudes = [float(a["amplitude"]) for a in manifest["assignments"]]
    plan = {
        "assignments_match_protocol": assigned == expected,
        "amplitudes_match_protocol": all(
            (a == 0) if k == "flat" else abs(a - amplitude) < 1e-15
            for (k, _), a in zip(assigned, amplitudes, strict=True)
        ),
        "levels_match_protocol": list(manifest["measurement_levels"]) == list(LEVELS),
        "orders_match_protocol": list(manifest["orders"]) == list(ORDERS),
        "tolerance_matches_protocol": manifest["tolerance_on_kappa"] == [LOW, HIGH],
    }
    recorded = {(r["replicate"], r["n"]): r for r in summary["rows"]}
    rows, truths, failures = [], {}, []
    for replicate, (kind, _) in enumerate(assigned):
        truth = None
        if kind == "constant":
            with np.load(folder / f"truth_{replicate}.npz", allow_pickle=False) as f:
                points, value = f["points"], float(str(f["amplitude"].item()))
            truth, residual = constant_truth(points, value)
            truths[replicate] = residual
        for n in LEVELS:
            row = recorded.get((replicate, n))
            if row is None or "failure" in row or "evaluation" not in row:
                failures.append({"replicate": replicate, "n": n, "row": row})
                continue
            rows.append(check_row(folder, replicate, kind, n, row, truth))
    prefixes = {r: nested(folder, r) for r in range(len(assigned))}
    verdicts, descriptive = {}, {}
    for kind, count, _ in POPULATIONS:
        for order in ORDERS:
            for n in LEVELS:
                key = f"{kind}/order_{order}/n{n}"
                values = [
                    r[f"order_{order}"]["truth_statistic"]
                    for r in rows
                    if r["kind"] == kind and r["measurement_n"] == n
                ]
                verdicts[key] = judge(values, count)
                if "kappa" in verdicts[key]:
                    descriptive[key] = posthoc(values, verdicts[key]["kappa"])
    run = summary["calibration"]
    match = {
        key: verdict["verdict"] == run[key]["verdict"]
        and verdict.get("direction") == run[key].get("direction")
        and (
            "kappa" not in verdict
            or (
                close(verdict["kappa"], run[key]["kappa"], 1e-9)
                and all(
                    close(a, b, 1e-9)
                    for a, b in zip(
                        verdict["interval_95"], run[key]["interval_95"], strict=True
                    )
                )
                and verdict["exceedances_of_q_star"]
                == run[key]["exceedances_of_q_star"]
            )
        )
        for key, verdict in verdicts.items()
    }
    result = {
        "audited_folder": str(folder),
        "audit_code_version": git("rev-parse", "HEAD").stdout.decode().strip(),
        "audit_worktree_clean": not git("status", "--porcelain").stdout.strip(),
        "run_code_version": stamp,
        "run_stamp_is_ancestor_of_head": git(
            "merge-base", "--is-ancestor", stamp, "HEAD"
        ).returncode
        == 0,
        "bound_sources_match_stamp": sources,
        "summary_sha256": sha(folder / "summary.json"),
        "plan": plan,
        "failed_rows": failures,
        "constant_truth_basis_residual": max(truths.values(), default=0.0),
        "rows": rows,
        "nested_prefixes_bit_exact": prefixes,
        "verdicts": verdicts,
        "verdicts_match_run": match,
        "posthoc": {
            "descriptive": descriptive,
            "kappa_ratio_sixth_to_fifth": {
                f"{kind}/n{n}": verdicts[f"{kind}/order_6/n{n}"]["kappa"]
                / verdicts[f"{kind}/order_5/n{n}"]["kappa"]
                for kind, _, _ in POPULATIONS
                for n in LEVELS
                if "kappa" in verdicts[f"{kind}/order_5/n{n}"]
                and "kappa" in verdicts[f"{kind}/order_6/n{n}"]
            },
        },
        "refits": [refit(folder, r, n) for r, n in refits],
    }
    result["passed"] = bool(
        result["run_stamp_is_ancestor_of_head"]
        and all(sources.values())
        and all(plan.values())
        and result["constant_truth_basis_residual"] < 1e-12
        and all(
            r["blind_boundary_kept"]
            and r["archive_hash_matches"]
            and r["measurement_n_matches"]
            and r["recorded_values_reproduced"]
            for r in rows
        )
        and all(prefixes.values())
        and all(match.values())
        and all(
            v["quantiles_vs_wilson_hilferty"] < 1e-3
            for v in verdicts.values()
            if "kappa" in v
        )
        and all(not r["differences_beyond_1e-8"] for r in result["refits"])
    )
    write_json(receipt, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", type=Path)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--refit", action="append", default=[], metavar="REPLICATE:N")
    args = parser.parse_args()
    outcome = audit(
        args.folder.resolve(),
        args.receipt or args.folder.with_name(args.folder.name + ".audit.json"),
        [tuple(int(part) for part in item.split(":")) for item in args.refit],
    )
    print(json.dumps({k: outcome[k] for k in ("verdicts", "verdicts_match_run")}))
    print(json.dumps(outcome["posthoc"]))
    print(json.dumps({"refits": outcome["refits"], "passed": outcome["passed"]}))
    sys.exit(0 if outcome["passed"] else 1)
