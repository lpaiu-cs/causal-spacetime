"""Independent clean-room verifier for S2 confirmatory subspacing results.

Audits five distinct layers independently from raw S2 data:
  1. Git provenance chain (P2 <= F2 <= HEAD, frozen source SHA-256 unchanged).
  2. Strict seed disjointness against all spent ranges (including S1 51M and <1M).
  3. Static AST blindness check verifying observer imports no geometry or oracle.
  4. Independent SciPy quadrature order oracle checked against WaveGeometry/OrderOracle
     at S2 sample points across (q=4, r=96, 128) and (q=1, r=48).
  5. Blind spot-check replay of observe(...) (q=None) on representative held-out S2
     runs plus full recomputation of all 20 two-sided 95% CI gates from raw runs.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.stats import norm
from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments" / "exploratory"))
sys.path.insert(0, str(ROOT / "experiments" / "positive_control"))

from seed_windows import (  # noqa: E402
    P11_P13_SPENT_RANGES,
    P12_ALLOCATION_DECADE,
    assert_windows_disjoint_and_fresh,
)
from subspacing_observer import observe  # noqa: E402
from subspacing_wave_sim import OrderOracle, WaveGeometry  # noqa: E402

FORBIDDEN_OBSERVER_TOKENS = frozenset(
    {
        "subspacing_wave_sim",
        "WaveGeometry",
        "OrderOracle",
        "quad",
        "intrinsic_weyl",
    }
)
FROZEN_Z_TH = 2.50
FROZEN_LOCK_WIDTH = 8.0
FROZEN_FREQUENCIES = np.arange(24.0, 257.0, 4.0)
FROZEN_BLOCKS = [
    ("N22", 52_000_000, 52_000_027),
    ("N23", 52_100_000, 52_100_027),
    ("W1", 52_200_000, 52_200_035),
    ("W1p_sprinkle", 52_300_000, 52_300_035),
    ("W1p_phase", 52_310_000, 52_310_035),
    ("W2", 52_400_000, 52_400_031),
    ("W3", 52_500_000, 52_500_027),
    ("X_null_q1", 52_600_000, 52_600_015),
    ("X_sub_q4_r32", 52_610_000, 52_610_023),
    ("X_sub_q1_r48", 52_620_000, 52_620_023),
    ("X_weyl_r4_pts", 52_630_000, 52_630_007),
    ("X_weyl_r4_lm", 52_635_000, 52_635_007),
    ("X_weyl_r32_pts", 52_640_000, 52_640_007),
    ("X_weyl_r32_lm", 52_645_000, 52_645_007),
]


def audit_observer_ast(observer_path: Path | None = None) -> dict:
    """Assert subspacing_observer.py has zero truth/wave generator AST references."""
    path = observer_path or (
        ROOT / "experiments" / "exploratory" / "subspacing_observer.py"
    )
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    imports: set[str] = set()
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                imports.add(alias.name)
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                imports.add(node.module)
            for alias in node.names:
                imports.add(alias.name)
        elif isinstance(node, ast.Name):
            names.add(node.id)
        elif isinstance(node, ast.Attribute):
            names.add(node.attr)
    bad = sorted(FORBIDDEN_OBSERVER_TOKENS & (imports | names))
    if bad:
        raise AssertionError(f"Observer AST leaked truth symbols: {bad}")
    return {"passed": True, "imports": sorted(imports), "leaks": []}


def independent_rosen_order_check(seed: int = 52_400_000, n_pts: int = 120) -> dict:
    """Verify OrderOracle via direct scipy.integrate.quad across S2 wave regimes."""
    mismatches = 0
    pairs_checked = 0
    max_abs_matrix_err = 0.0
    for q, r in ((4.0, 96.0), (4.0, 128.0), (1.0, 48.0)):
        k = float(r)
        lam = 1.0 / k
        phi = 0.37
        t_param = 2.0 * lam * np.sqrt(q) / np.sqrt(1.0 + 4.0 * q * lam**2)
        geom = WaveGeometry(q, "circular_wave", lam, phi=phi)
        rng = np.random.default_rng(seed + int(r))
        pts = geom.sprinkle(n_pts, rng, "cpu")
        oracle_fw = OrderOracle(geom, pts, "cpu").query(
            np.arange(n_pts), np.arange(n_pts)
        )[0]
        sq = np.sqrt(q)
        for i in range(n_pts):
            for j in range(i + 1, min(i + 4, n_pts)):
                for a_idx, b_idx in ((i, j), (j, i)):
                    p1, p2 = pts[a_idx], pts[b_idx]
                    u1, u2 = float(p1[0]), float(p2[0])
                    if u2 > u1:
                        m_val = float((np.tan(sq * u2) - np.tan(sq * u1)) / sq)
                        c1 = quad(
                            lambda w, s=k, ph=phi, r_q=sq: (
                                np.cos(s * w + ph) / np.cos(r_q * w) ** 2
                            ),
                            u1,
                            u2,
                            limit=500,
                            epsabs=1e-12,
                            epsrel=1e-12,
                        )[0]
                        s1 = quad(
                            lambda w, s=k, ph=phi, r_q=sq: (
                                np.sin(s * w + ph) / np.cos(r_q * w) ** 2
                            ),
                            u1,
                            u2,
                            limit=500,
                            epsabs=1e-12,
                            epsrel=1e-12,
                        )[0]
                        mat = np.array(
                            [
                                [m_val - t_param * c1, -t_param * s1],
                                [-t_param * s1, m_val + t_param * c1],
                            ]
                        )
                        dy = p2[2:] - p1[2:]
                        req = 0.5 * float(dy @ np.linalg.solve(mat, dy))
                        ind_rel = bool((p2[1] - p1[1]) >= req)
                        mat_diff = float(np.max(np.abs(mat - geom.matrix(u1, u2))))
                        max_abs_matrix_err = max(max_abs_matrix_err, mat_diff)
                    elif u2 == u1:
                        ind_rel = bool(np.all(p1[2:] == p2[2:]) and p2[1] > p1[1])
                    else:
                        ind_rel = False
                    if ind_rel != bool(oracle_fw[a_idx, b_idx]):
                        mismatches += 1
                    pairs_checked += 1
    if mismatches:
        raise AssertionError(f"Independent quad oracle had {mismatches} mismatches")
    return {
        "passed": True,
        "pairs_checked": pairs_checked,
        "mismatches": mismatches,
        "max_abs_matrix_err": max_abs_matrix_err,
    }


def _wilson(k: int, n: int) -> tuple[float, float]:
    z = float(norm.ppf(0.975))
    p = k / n
    d = 1.0 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    rad = z * np.sqrt((p * (1 - p) + z**2 / (4 * n)) / n) / d
    return float(max(0.0, c - rad)), float(min(1.0, c + rad))


def _t_ci(vals) -> tuple[float, float, float]:
    a = np.asarray(vals, dtype=float)
    se = float(np.std(a, ddof=1) / np.sqrt(len(a)))
    crit = float(student_t.ppf(0.975, df=len(a) - 1))
    m = float(np.mean(a))
    return m, m - crit * se, m + crit * se


def _welch_ci(a, b) -> tuple[float, float, float]:
    xa, xb = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    va, vb = np.var(xa, ddof=1) / len(xa), np.var(xb, ddof=1) / len(xb)
    df = (va + vb) ** 2 / (va**2 / (len(xa) - 1) + vb**2 / (len(xb) - 1))
    crit = float(student_t.ppf(0.975, df=df))
    diff = float(np.mean(xa) - np.mean(xb))
    se = float(np.sqrt(va + vb))
    return diff, diff - crit * se, diff + crit * se


def recompute_and_verify_gates(results: dict) -> dict:
    """Recompute all 20 S2 gates strictly from raw results['runs'] & verify parity."""
    runs = results["runs"]
    by_arm: dict[str, list[dict]] = {}
    for row in runs:
        by_arm.setdefault(row["arm"], []).append(row)

    expected_counts = {
        "N22": 28,
        "N23": 28,
        "W1": 36,
        "W1p": 36,
        "W2": 32,
        "W3": 28,
        "X_null_q1": 16,
        "X_sub_q4_r32": 24,
        "X_sub_q1_r48": 24,
        "X_weyl_r4": 8,
        "X_weyl_r32": 8,
    }
    if {k: len(v) for k, v in by_arm.items()} != expected_counts:
        raise AssertionError("Arm count mismatch against S2 preregistration")

    seen_seeds = [int(r["seed"]) for r in runs]
    if len(set(seen_seeds)) != len(runs):
        raise AssertionError("Duplicate primary seeds in S2 results")

    checks: dict[str, bool] = {}
    null_all = by_arm["N22"] + by_arm["N23"]
    fp_k = sum(1 for r in null_all if float(r["whitened_peak_z"]) >= FROZEN_Z_TH)
    lf_k = sum(1 for r in null_all if float(r["whitened_peak_frequency"]) <= 32.0)
    checks["G1"] = _wilson(fp_k, 56)[1] <= 0.18
    checks["G1b"] = _wilson(lf_k, 56)[1] <= 0.20
    checks["G2a"] = (
        all(r["c_ref_estimate"] > 0 for r in null_all)
        and _t_ci([r["calibration_rms"] for r in null_all])[2] <= 1.2e-4
    )
    checks["G2b"] = (
        _t_ci([r["q_abs_err"] for r in by_arm["N22"]])[2] <= 0.035
        and _t_ci([r["q_abs_err"] for r in by_arm["N23"]])[2] <= 0.025
    )

    def _chk(arm, ctrl, r_val, tpr_lo, fl_lo, d_lo, q_hi):
        rows = by_arm[arm]
        n = len(rows)
        tk = sum(1 for r in rows if float(r["whitened_peak_z"]) >= FROZEN_Z_TH)
        fk = sum(
            1
            for r in rows
            if abs(float(r["whitened_peak_frequency"]) - r_val) <= FROZEN_LOCK_WIDTH
        )
        d_ci = _welch_ci(
            [r["whitened_peak_z"] for r in rows],
            [r["whitened_peak_z"] for r in by_arm[ctrl]],
        )
        q_ci = _t_ci([r["q_abs_err"] for r in rows])
        return {
            "tpr": _wilson(tk, n)[0] >= tpr_lo,
            "fl": _wilson(fk, n)[0] >= fl_lo,
            "sep": d_ci[1] >= d_lo,
            "q": True if q_hi is None else (q_ci[2] <= q_hi),
        }

    c_w1 = _chk("W1", "N22", 64.0, 0.50, 0.65, 0.65, 0.030)
    c_w1p = _chk("W1p", "N22", 64.0, 0.50, 0.65, 0.50, 0.030)
    c_w2 = _chk("W2", "N23", 96.0, 0.50, 0.65, 1.20, 0.025)
    c_w3 = _chk("W3", "N23", 128.0, 0.50, 0.60, 0.90, 0.025)
    c_xq4 = _chk("X_sub_q4_r32", "N22", 32.0, 0.50, 0.55, 0.55, None)
    c_xq1 = _chk("X_sub_q1_r48", "X_null_q1", 48.0, 0.50, 0.55, 0.40, None)

    checks["G3a"], checks["G3b"] = c_w1["tpr"], c_w1["sep"]
    checks["G3c"] = c_w1["fl"] and c_w1["q"]
    checks["G4a"], checks["G4b"] = c_w1p["tpr"], c_w1p["sep"]
    checks["G4c"] = c_w1p["fl"] and c_w1p["q"]
    checks["G5a"], checks["G5b"] = c_w2["tpr"], c_w2["sep"]
    checks["G5c"] = c_w2["fl"] and c_w2["q"]
    checks["G6a"], checks["G6b"] = c_w3["tpr"], c_w3["sep"]
    checks["G6c"] = c_w3["fl"] and c_w3["q"]
    checks["G7a"] = _t_ci([r["electric_ratio"] for r in by_arm["X_weyl_r4"]])[1] >= 0.65
    checks["G7b"] = (
        _t_ci([r["electric_ratio"] for r in by_arm["X_weyl_r32"]])[2] <= 0.10
    )
    checks["G7c"] = c_xq4["tpr"] and c_xq4["fl"] and c_xq4["sep"]
    checks["G7d"] = c_xq1["tpr"] and c_xq1["fl"] and c_xq1["sep"]

    recorded = results["summary"]["gates"]
    disagreements = [
        g for g, ok in checks.items() if bool(recorded[g]["pass"]) != bool(ok)
    ]
    if disagreements:
        raise AssertionError(f"S2 gate disagreements: {disagreements}")
    expected_verdict = "CONFIRMED" if all(checks.values()) else "MIXED_OR_FAIL"
    if results["summary"]["verdict"] != expected_verdict:
        raise AssertionError("Summary verdict disagrees with S2 recomputation")
    return {"verdict": expected_verdict, "checks": checks, "disagreements": []}


def blind_spot_check_runs(results: dict) -> list[dict]:
    """Re-execute the first run of every S2 subspacing arm blindly on GPU."""
    import cupy as cp

    first_by_arm: dict[str, dict] = {}
    for row in results["runs"]:
        if row["kind"] == "subspacing" and row["arm"] not in first_by_arm:
            first_by_arm[row["arm"]] = row

    replays = []
    for arm, rec in first_by_arm.items():
        rng = np.random.default_rng(rec["seed"])
        n_draw = int(rng.poisson(rec["n"])) if rec["poisson"] else int(rec["n"])
        if n_draw != rec["n_actual"]:
            raise AssertionError(f"{arm} n_actual mismatch")
        geom = WaveGeometry(rec["q"], rec["mode"], 1.0 / rec["r"], phi=rec["phi"])
        pts = geom.sprinkle(n_draw, rng, "gpu")
        oracle = OrderOracle(geom, pts, "gpu")
        pilot_n = min(16384, n_draw // 4)
        res = observe(
            oracle.query,
            np.arange(pilot_n),
            np.arange(n_draw // 2, n_draw),
            np.arange(pilot_n, n_draw // 2),
            FROZEN_FREQUENCIES,
            q=None,
            count_query=oracle.counts,
            max_probes=65536,
        )
        del oracle, pts
        cp.get_default_memory_pool().free_all_blocks()
        dz = abs(float(res["whitened_peak_z"]) - float(rec["whitened_peak_z"]))
        dq = abs(float(res["q_hat"]) - float(rec["q_hat"]))
        if max(dz, dq) > 1e-10:
            raise AssertionError(f"Replay mismatch in {arm}: dz={dz}, dq={dq}")
        replays.append({"arm": arm, "seed": rec["seed"], "dz": dz, "dq": dq})
    return replays


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--frozen-dir",
        type=Path,
        default=ROOT / "docs" / "prereg" / "frozen" / "subspacing_s2",
    )
    parser.add_argument("--skip-replay", action="store_true")
    args = parser.parse_args()

    t0 = time.perf_counter()
    spent = (
        P11_P13_SPENT_RANGES
        + (P12_ALLOCATION_DECADE,)
        + (
            ("Exploratory Steps 1-3 & legacy (<1M)", 0, 999_999),
            ("P14 allocation decade", 40_000_000, 49_999_999),
            ("Confirmatory S1 campaign block", 51_000_000, 51_999_999),
        )
    )
    assert_windows_disjoint_and_fresh(
        list(FROZEN_BLOCKS), spent, 52_000_000, 52_999_999, "S2_check"
    )
    ast_res = audit_observer_ast()
    oracle_res = independent_rosen_order_check()

    manifest = json.loads((args.frozen_dir / "freeze_manifest.json").read_text("utf-8"))
    results = json.loads((args.frozen_dir / "results.json").read_text("utf-8"))

    for rel, expected_sha in manifest["source_sha256"].items():
        actual_sha = hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()
        if actual_sha != expected_sha:
            raise AssertionError(f"Source digest changed for {rel}")
    head_sha = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    p_sha = manifest["preflight_code_version"]
    f_sha = results["code_version"]
    for a, b in ((p_sha, f_sha), (f_sha, head_sha)):
        rc = subprocess.run(
            ["git", "merge-base", "--is-ancestor", a, b], cwd=ROOT
        ).returncode
        if rc != 0:
            raise AssertionError(f"Git ancestry violated: {a} is not ancestor of {b}")

    gates_res = recompute_and_verify_gates(results)
    replays = [] if args.skip_replay else blind_spot_check_runs(results)
    receipt = {
        "audit_passed": True,
        "head_commit": head_sha,
        "preflight_commit": p_sha,
        "campaign_commit": f_sha,
        "elapsed_seconds": time.perf_counter() - t0,
        "seed_freshness": "disjoint from S1 (51M), P11-P14, and <1M ranges",
        "ast_blindness": ast_res,
        "independent_quad_oracle": oracle_res,
        "recomputed_gates": gates_res,
        "spot_replays": replays,
    }
    out = args.frozen_dir / "audit_receipt.json"
    out.write_text(json.dumps(receipt, indent=2), encoding="utf-8", newline="\n")
    print(f"Independent S2 audit PASSED ({gates_res['verdict']}) -> {out}")


if __name__ == "__main__":
    main()
