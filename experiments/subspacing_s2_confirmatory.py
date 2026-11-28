"""Preregistered confirmatory S2 runner for w^-2 + w^-4 whitened detection.

Preserves strict P2 -> F2 -> R2 commit provenance on a clean worktree, enforces
seed disjointness against all spent ranges (including exploratory and S1 51M),
runs all 268 held-out replicates on GPU, and evaluates 20 two-sided 95% CI gates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import norm
from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments" / "exploratory"))
sys.path.insert(0, str(ROOT / "experiments" / "positive_control"))

import intrinsic_weyl_joint as iwj  # noqa: E402
from seed_windows import (  # noqa: E402
    P11_P13_SPENT_RANGES,
    P12_ALLOCATION_DECADE,
    assert_windows_disjoint_and_fresh,
)
from subspacing_observer import observe  # noqa: E402
from subspacing_wave_sim import OrderOracle, WaveGeometry  # noqa: E402

Z_95 = float(norm.ppf(0.975))
FREQUENCIES = np.arange(24.0, 257.0, 4.0)
Z_TH = 2.50
LOCK_HALF_WIDTH = 8.0
MAX_PROBES = 65536
PILOT_CAP = 16384

S2_WINDOW_FLOOR = 52_000_000
S2_WINDOW_CEIL = 52_999_999

EXTRA_SPENT_RANGES = (
    ("Exploratory Steps 1-3 & legacy (<1M)", 0, 999_999),
    ("P14 allocation decade", 40_000_000, 49_999_999),
    ("Confirmatory S1 campaign block", 51_000_000, 51_999_999),
)
ALL_SPENT_RANGES = P11_P13_SPENT_RANGES + (P12_ALLOCATION_DECADE,) + EXTRA_SPENT_RANGES

SEED_BLOCKS = [
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

FROZEN_SOURCE_PATHS = (
    "docs/prereg/subspacing_s2_confirmatory.md",
    "experiments/exploratory/subspacing_observer.py",
    "experiments/exploratory/subspacing_wave_sim.py",
    "experiments/positive_control/intrinsic_weyl_joint.py",
    "experiments/positive_control/seed_windows.py",
    "experiments/subspacing_s2_confirmatory.py",
    "experiments/check_subspacing_s2_confirmatory.py",
)

DEFAULT_FROZEN_DIR = ROOT / "docs" / "prereg" / "frozen" / "subspacing_s2"
DEFAULT_CKPT = ROOT / "outputs" / "subspacing_s2_ckpt.json"


def check_seed_freshness() -> list[tuple[str, int, int]]:
    """Verify every S2 block is disjoint from spent history and within 52M."""
    return assert_windows_disjoint_and_fresh(
        list(SEED_BLOCKS),
        ALL_SPENT_RANGES,
        S2_WINDOW_FLOOR,
        S2_WINDOW_CEIL,
        "S2",
    )


def git_head_and_status(repo: Path = ROOT) -> tuple[str, str]:
    """Return (full_sha, porcelain_status) from git."""
    sha = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=repo, text=True
    ).strip()
    status = subprocess.check_output(
        ["git", "status", "--porcelain"], cwd=repo, text=True
    ).strip()
    return sha, status


def source_digests(repo: Path = ROOT) -> dict[str, str]:
    """Compute SHA-256 hex digests for every frozen protocol and code file."""
    return {
        rel: hashlib.sha256((repo / rel).read_bytes()).hexdigest()
        for rel in FROZEN_SOURCE_PATHS
    }


def wilson_ci(k: int, n: int, z: float = Z_95) -> tuple[float, float]:
    """Exact two-sided Wilson score confidence interval for k / n."""
    if n <= 0:
        raise ValueError("n must be positive")
    p = k / n
    den = 1.0 + z**2 / n
    ctr = (p + z**2 / (2.0 * n)) / den
    rad = z * np.sqrt((p * (1.0 - p) + z**2 / (4.0 * n)) / n) / den
    return float(max(0.0, ctr - rad)), float(min(1.0, ctr + rad))


def student_mean_ci(values) -> dict[str, float]:
    """Exact two-sided 95% Student-t CI for a single-sample mean."""
    arr = np.asarray(values, dtype=float)
    n = len(arr)
    mean = float(np.mean(arr))
    sd = float(np.std(arr, ddof=1))
    se = sd / float(np.sqrt(n))
    crit = float(student_t.ppf(0.975, df=n - 1))
    return {
        "n": int(n),
        "mean": mean,
        "sd": sd,
        "se": se,
        "df": int(n - 1),
        "ci95_lo": float(mean - crit * se),
        "ci95_hi": float(mean + crit * se),
    }


def welch_diff_ci(treat, control) -> dict[str, float]:
    """Exact two-sided 95% Welch-Satterthwaite t-CI for mean(treat) - mean(control)."""
    a = np.asarray(treat, dtype=float)
    b = np.asarray(control, dtype=float)
    va = float(np.var(a, ddof=1)) / len(a)
    vb = float(np.var(b, ddof=1)) / len(b)
    se = float(np.sqrt(va + vb))
    if se == 0.0:
        df = float(len(a) + len(b) - 2)
    else:
        df = float((va + vb) ** 2 / (va**2 / (len(a) - 1) + vb**2 / (len(b) - 1)))
    crit = float(student_t.ppf(0.975, df=df))
    diff = float(np.mean(a) - np.mean(b))
    return {
        "n_treat": int(len(a)),
        "n_control": int(len(b)),
        "treat_mean": float(np.mean(a)),
        "control_mean": float(np.mean(b)),
        "diff": diff,
        "se": se,
        "df": df,
        "ci95_lo": float(diff - crit * se),
        "ci95_hi": float(diff + crit * se),
    }


def build_jobs() -> list[dict]:
    """Deterministic ordered list of all 268 S2 confirmatory replicates."""
    check_seed_freshness()
    jobs: list[dict] = []

    for i in range(28):
        jobs.append(
            dict(
                arm="N22",
                kind="subspacing",
                n=2**22,
                q=4.0,
                mode="background",
                r=64,
                phi=0.37,
                poisson=False,
                seed=52_000_000 + i,
            )
        )
    for i in range(28):
        jobs.append(
            dict(
                arm="N23",
                kind="subspacing",
                n=2**23,
                q=4.0,
                mode="background",
                r=96,
                phi=0.37,
                poisson=False,
                seed=52_100_000 + i,
            )
        )
    for i in range(36):
        jobs.append(
            dict(
                arm="W1",
                kind="subspacing",
                n=2**22,
                q=4.0,
                mode="circular_wave",
                r=64,
                phi=0.37,
                poisson=False,
                seed=52_200_000 + i,
            )
        )
    for i in range(36):
        phi_rng = np.random.default_rng(52_310_000 + i)
        phi = float(phi_rng.uniform(0.0, 2.0 * np.pi))
        jobs.append(
            dict(
                arm="W1p",
                kind="subspacing",
                n=2**22,
                q=4.0,
                mode="circular_wave",
                r=64,
                phi=phi,
                phi_seed=52_310_000 + i,
                poisson=True,
                seed=52_300_000 + i,
            )
        )
    for i in range(32):
        jobs.append(
            dict(
                arm="W2",
                kind="subspacing",
                n=2**23,
                q=4.0,
                mode="circular_wave",
                r=96,
                phi=0.37,
                poisson=False,
                seed=52_400_000 + i,
            )
        )
    for i in range(28):
        jobs.append(
            dict(
                arm="W3",
                kind="subspacing",
                n=2**23,
                q=4.0,
                mode="circular_wave",
                r=128,
                phi=0.37,
                poisson=False,
                seed=52_500_000 + i,
            )
        )
    for i in range(16):
        jobs.append(
            dict(
                arm="X_null_q1",
                kind="subspacing",
                n=2**22,
                q=1.0,
                mode="background",
                r=48,
                phi=0.37,
                poisson=False,
                seed=52_600_000 + i,
            )
        )
    for i in range(24):
        jobs.append(
            dict(
                arm="X_sub_q4_r32",
                kind="subspacing",
                n=2**22,
                q=4.0,
                mode="circular_wave",
                r=32,
                phi=0.37,
                poisson=False,
                seed=52_610_000 + i,
            )
        )
    for i in range(24):
        jobs.append(
            dict(
                arm="X_sub_q1_r48",
                kind="subspacing",
                n=2**23,
                q=1.0,
                mode="circular_wave",
                r=48,
                phi=0.37,
                poisson=False,
                seed=52_620_000 + i,
            )
        )
    for i in range(8):
        jobs.append(
            dict(
                arm="X_weyl_r4",
                kind="weyl_jet",
                n=2**21,
                q=0.25,
                mode="circular_wave",
                r=4.0,
                phi=0.37,
                poisson=False,
                seed=52_630_000 + i,
                lm_seed=52_635_000 + i,
            )
        )
    for i in range(8):
        jobs.append(
            dict(
                arm="X_weyl_r32",
                kind="weyl_jet",
                n=2**21,
                q=0.25,
                mode="circular_wave",
                r=32.0,
                phi=0.37,
                poisson=False,
                seed=52_640_000 + i,
                lm_seed=52_645_000 + i,
            )
        )
    return jobs


def _run_subspacing_one(job: dict) -> dict:
    import cupy as cp

    t0 = time.perf_counter()
    rng = np.random.default_rng(job["seed"])
    n_draw = int(rng.poisson(job["n"])) if job["poisson"] else int(job["n"])
    geom = WaveGeometry(job["q"], job["mode"], 1.0 / job["r"], phi=job["phi"])
    pts = geom.sprinkle(n_draw, rng, "gpu")
    oracle = OrderOracle(geom, pts, "gpu")
    gen_s = time.perf_counter() - t0

    pilot_n = min(PILOT_CAP, n_draw // 4)
    pilot = np.arange(pilot_n)
    probes = np.arange(pilot_n, n_draw // 2)
    references = np.arange(n_draw // 2, n_draw)
    cache: dict[tuple[bytes, ...], np.ndarray] = {}

    def counts(refs, at, early=(), late=(), _oracle=oracle):
        key = tuple(
            hashlib.blake2b(
                np.asarray(g, dtype=np.int64).tobytes(), digest_size=16
            ).digest()
            for g in (refs, at, early, late)
        )
        if key not in cache:
            cache[key] = _oracle.counts(refs, at, early, late)
        return cache[key]

    try:
        res = observe(
            oracle.query,
            pilot,
            references,
            probes,
            FREQUENCIES,
            q=None,
            count_query=counts,
            max_probes=MAX_PROBES,
        )
        pk_freq = float(res["whitened_peak_frequency"])
        pk_z = float(res["whitened_peak_z"])
        return {
            **job,
            "n_actual": n_draw,
            "alpha": float(np.log(n_draw) / np.log(job["r"])),
            "generation_seconds": gen_s,
            "total_seconds": time.perf_counter() - t0,
            "whitened_peak_z": pk_z,
            "whitened_snr": float(res["whitened_snr"]),
            "whitened_peak_frequency": pk_freq,
            "detected": bool(pk_z >= Z_TH),
            "freq_locked": bool(abs(pk_freq - job["r"]) <= LOCK_HALF_WIDTH),
            "low_freq_peak": bool(pk_freq <= 32.0),
            "q_hat": float(res["q_hat"]),
            "q_abs_err": float(abs(res["q_hat"] - job["q"])),
            "calibration_rms": float(res["calibration_rms"]),
            "c_ref_estimate": float(res["c_ref_estimate"]),
            "c1_estimate": float(res["c1_estimate"]),
            "c2_steep_estimate": float(res["c2_steep_estimate"]),
            "shot_sd": float(res["shot_sd"]),
            "accepted": int(res["accepted"]),
        }
    finally:
        del oracle, pts
        cp.get_default_memory_pool().free_all_blocks()


def _run_weyl_one(job: dict) -> dict:
    import cupy as cp

    t0 = time.perf_counter()
    geom = WaveGeometry(job["q"], job["mode"], 1.0 / job["r"], phi=job["phi"])
    sel_int = geom.sprinkle(2048, np.random.default_rng(job["lm_seed"]), "cpu")
    sel_pts = np.vstack(([-0.5, 0.0, 0.0, 0.0], sel_int, [0.5, 1.0, 0.0, 0.0]))
    sel_ord = OrderOracle(geom, sel_pts, "cpu").query(
        np.arange(len(sel_pts)), np.arange(len(sel_pts))
    )[0]
    ids, mask = iwj.selection(sel_ord, design="two")
    lm = sel_pts[ids]
    pts = geom.sprinkle(job["n"], np.random.default_rng(job["seed"]), "gpu")
    all_pts = np.vstack((lm, pts))
    oracle = OrderOracle(geom, all_pts, "gpu")
    try:
        d_lm = np.arange(len(lm))
        d_ev = np.arange(len(lm), len(all_pts))
        counts = np.zeros((len(lm), len(lm)), dtype=np.float64)
        tile = 32768
        for s in range(0, job["n"], tile):
            sub = d_ev[s : s + tile]
            fw, bw = oracle.query(d_lm, sub)
            counts += fw.astype(np.float64) @ bw.T.astype(np.float64)
        lm_ord = sel_ord[np.ix_(ids, ids)]
        causal = np.where(lm_ord, counts, 0.0)
        fit = iwj.fit_joint(causal + causal.T, clock_mask=mask, starts=1, model_order=5)
        e2 = float(fit["electric_norm_squared"])
        sqrt_e2 = float(np.sqrt(max(e2, 0.0)))
        amp = float(np.sqrt(job["q"] * (job["r"] ** 2 + 4.0 * job["q"])))
        ratio = float(sqrt_e2 / (np.sqrt(2.0) * amp))
        return {
            **job,
            "total_seconds": time.perf_counter() - t0,
            "e2_hat": e2,
            "sqrt_e2_hat": sqrt_e2,
            "amp_true": amp,
            "electric_ratio": ratio,
            "c_norm": float(np.linalg.norm(fit["coefficients"])),
            "optimizer_success": bool(fit["optimizer_success"]),
        }
    finally:
        del oracle, pts, all_pts
        cp.get_default_memory_pool().free_all_blocks()


def evaluate_gates(runs: list[dict]) -> dict:
    """Evaluate all 20 frozen S2 two-sided 95% CI gates from raw run dicts."""
    by_arm: dict[str, list[dict]] = {}
    for r in runs:
        by_arm.setdefault(r["arm"], []).append(r)

    n22, n23 = by_arm["N22"], by_arm["N23"]
    w1, w1p = by_arm["W1"], by_arm["W1p"]
    w2, w3 = by_arm["W2"], by_arm["W3"]
    xq1_null = by_arm["X_null_q1"]
    xq4_r32, xq1_r48 = by_arm["X_sub_q4_r32"], by_arm["X_sub_q1_r48"]
    xw_r4, xw_r32 = by_arm["X_weyl_r4"], by_arm["X_weyl_r32"]

    pooled_null = n22 + n23
    fp_k = sum(1 for r in pooled_null if r["detected"])
    fp_ci = wilson_ci(fp_k, len(pooled_null))
    lf_k = sum(1 for r in pooled_null if float(r["whitened_peak_frequency"]) <= 32.0)
    lf_ci = wilson_ci(lf_k, len(pooled_null))
    null_rms_ci = student_mean_ci([r["calibration_rms"] for r in pooled_null])
    null_c0 = sum(1 for r in pooled_null if r["c_ref_estimate"] <= 0.0)
    n22_qerr = student_mean_ci([r["q_abs_err"] for r in n22])
    n23_qerr = student_mean_ci([r["q_abs_err"] for r in n23])

    gates: dict[str, dict] = {}
    gates["G1"] = {
        "criterion": "Pooled Null FPR upper95 <= 0.18",
        "k": fp_k,
        "n": len(pooled_null),
        "rate": fp_k / len(pooled_null),
        "ci95": list(fp_ci),
        "pass": bool(len(pooled_null) == 56 and fp_ci[1] <= 0.18),
    }
    gates["G1b"] = {
        "criterion": "Pooled Null low-freq peak P(pf <= 32) upper95 <= 0.20",
        "k": lf_k,
        "n": len(pooled_null),
        "rate": lf_k / len(pooled_null),
        "ci95": list(lf_ci),
        "pass": bool(len(pooled_null) == 56 and lf_ci[1] <= 0.20),
    }
    gates["G2a"] = {
        "criterion": "Pooled Null valid fits == 56, cal_rms upper95 <= 1.2e-4",
        "c_ref_zero_count": null_c0,
        **null_rms_ci,
        "pass": bool(
            len(pooled_null) == 56 and null_c0 == 0 and null_rms_ci["ci95_hi"] <= 1.2e-4
        ),
    }
    gates["G2b"] = {
        "criterion": "N22 |q_err| upper95 <= 0.035 AND N23 |q_err| upper95 <= 0.025",
        "N22_q_err": n22_qerr,
        "N23_q_err": n23_qerr,
        "pass": bool(n22_qerr["ci95_hi"] <= 0.035 and n23_qerr["ci95_hi"] <= 0.025),
    }

    def _wave(arm_rows, ctrl_rows, tpr_lo, fl_lo, diff_lo, q_hi):
        n = len(arm_rows)
        tk = sum(1 for r in arm_rows if r["detected"])
        tci = wilson_ci(tk, n)
        fk = sum(1 for r in arm_rows if r["freq_locked"])
        fci = wilson_ci(fk, n)
        diff = welch_diff_ci(
            [r["whitened_peak_z"] for r in arm_rows],
            [r["whitened_peak_z"] for r in ctrl_rows],
        )
        qerr = student_mean_ci([r["q_abs_err"] for r in arm_rows])
        return {
            "tpr": {
                "k": tk,
                "n": n,
                "rate": tk / n,
                "ci95": list(tci),
                "pass": bool(tci[0] >= tpr_lo),
            },
            "freq_lock": {
                "k": fk,
                "n": n,
                "rate": fk / n,
                "ci95": list(fci),
                "pass": bool(fci[0] >= fl_lo),
            },
            "z_sep": {**diff, "pass": bool(diff["ci95_lo"] >= diff_lo)},
            "q_err": {
                **qerr,
                "pass": bool(q_hi is None or qerr["ci95_hi"] <= q_hi),
            },
        }

    w1_s = _wave(w1, n22, 0.50, 0.65, 0.65, 0.030)
    w1p_s = _wave(w1p, n22, 0.50, 0.65, 0.50, 0.030)
    w2_s = _wave(w2, n23, 0.50, 0.65, 1.20, 0.025)
    w3_s = _wave(w3, n23, 0.50, 0.60, 0.90, 0.025)
    xq4_s = _wave(xq4_r32, n22, 0.50, 0.55, 0.55, None)
    xq1_s = _wave(xq1_r48, xq1_null, 0.50, 0.55, 0.40, None)

    gates["G3a"] = w1_s["tpr"]
    gates["G3b"] = w1_s["z_sep"]
    gates["G3c"] = {
        "freq_lock": w1_s["freq_lock"],
        "q_err": w1_s["q_err"],
        "pass": bool(w1_s["freq_lock"]["pass"] and w1_s["q_err"]["pass"]),
    }
    gates["G4a"] = w1p_s["tpr"]
    gates["G4b"] = w1p_s["z_sep"]
    gates["G4c"] = {
        "freq_lock": w1p_s["freq_lock"],
        "q_err": w1p_s["q_err"],
        "pass": bool(w1p_s["freq_lock"]["pass"] and w1p_s["q_err"]["pass"]),
    }
    gates["G5a"] = w2_s["tpr"]
    gates["G5b"] = w2_s["z_sep"]
    gates["G5c"] = {
        "freq_lock": w2_s["freq_lock"],
        "q_err": w2_s["q_err"],
        "pass": bool(w2_s["freq_lock"]["pass"] and w2_s["q_err"]["pass"]),
    }
    gates["G6a"] = w3_s["tpr"]
    gates["G6b"] = w3_s["z_sep"]
    gates["G6c"] = {
        "freq_lock": w3_s["freq_lock"],
        "q_err": w3_s["q_err"],
        "pass": bool(w3_s["freq_lock"]["pass"] and w3_s["q_err"]["pass"]),
    }

    wr4_ci = student_mean_ci([r["electric_ratio"] for r in xw_r4])
    wr32_ci = student_mean_ci([r["electric_ratio"] for r in xw_r32])
    gates["G7a"] = {**wr4_ci, "pass": bool(wr4_ci["ci95_lo"] >= 0.65)}
    gates["G7b"] = {**wr32_ci, "pass": bool(wr32_ci["ci95_hi"] <= 0.10)}
    gates["G7c"] = {
        "tpr": xq4_s["tpr"],
        "freq_lock": xq4_s["freq_lock"],
        "z_sep": xq4_s["z_sep"],
        "pass": bool(
            xq4_s["tpr"]["pass"]
            and xq4_s["freq_lock"]["pass"]
            and xq4_s["z_sep"]["pass"]
        ),
    }
    gates["G7d"] = {
        "tpr": xq1_s["tpr"],
        "freq_lock": xq1_s["freq_lock"],
        "z_sep": xq1_s["z_sep"],
        "pass": bool(
            xq1_s["tpr"]["pass"]
            and xq1_s["freq_lock"]["pass"]
            and xq1_s["z_sep"]["pass"]
        ),
    }

    arm_summaries = {
        arm: (
            {
                "n": len(rows),
                "z_max": student_mean_ci([r["whitened_peak_z"] for r in rows]),
                "detected_k": sum(1 for r in rows if r["detected"]),
                "freq_locked_k": sum(1 for r in rows if r["freq_locked"]),
                "low_freq_peak_k": sum(
                    1 for r in rows if float(r["whitened_peak_frequency"]) <= 32.0
                ),
                "q_abs_err": student_mean_ci([r["q_abs_err"] for r in rows]),
                "cal_rms": student_mean_ci([r["calibration_rms"] for r in rows]),
            }
            if rows[0]["kind"] == "subspacing"
            else {
                "n": len(rows),
                "electric_ratio": student_mean_ci([r["electric_ratio"] for r in rows]),
                "c_norm": student_mean_ci([r["c_norm"] for r in rows]),
            }
        )
        for arm, rows in by_arm.items()
    }
    all_pass = all(g["pass"] for g in gates.values())
    return {
        "verdict": "CONFIRMED" if all_pass else "MIXED_OR_FAIL",
        "gates": gates,
        "arm_summaries": arm_summaries,
    }


def run_preflight(frozen_dir: Path) -> Path:
    """Verify clean checkout P2, seed disjointness, GPU, and write freeze_manifest."""
    import cupy as cp

    check_seed_freshness()
    sha, status = git_head_and_status()
    if status:
        raise SystemExit(f"Refusing preflight on dirty worktree:\n{status}")
    if cp.cuda.runtime.getDeviceCount() < 1:
        raise SystemExit("CUDA GPU device required for S2 preflight")
    dev = cp.cuda.runtime.getDeviceProperties(0)["name"].decode("utf-8", "replace")
    digests = source_digests()
    jobs = build_jobs()
    manifest = {
        "prereg_doc": "docs/prereg/subspacing_s2_confirmatory.md",
        "preflight_code_version": sha,
        "gpu_device": dev,
        "z_th": Z_TH,
        "lock_half_width": LOCK_HALF_WIDTH,
        "max_probes": MAX_PROBES,
        "frequencies": FREQUENCIES.tolist(),
        "total_runs": len(jobs),
        "seed_blocks": [
            {"label": lb, "lo": lo, "hi": hi} for lb, lo, hi in SEED_BLOCKS
        ],
        "source_sha256": digests,
    }
    manifest_bytes = json.dumps(manifest, indent=2, sort_keys=True).encode("utf-8")
    manifest["manifest_sha256"] = hashlib.sha256(manifest_bytes).hexdigest()

    frozen_dir.mkdir(parents=True, exist_ok=True)
    out = frozen_dir / "freeze_manifest.json"
    out.write_text(json.dumps(manifest, indent=2), encoding="utf-8", newline="\n")
    print(f"Wrote S2 freeze manifest for HEAD={sha[:7]} to {out}")
    return out


def run_campaign(frozen_dir: Path, ckpt_path: Path) -> Path:
    """Execute full held-out S2 campaign on clean frozen commit F2."""
    import cupy as cp

    check_seed_freshness()
    sha, status = git_head_and_status()
    if status:
        raise SystemExit(f"Refusing campaign on dirty worktree:\n{status}")
    manifest_path = frozen_dir / "freeze_manifest.json"
    if not manifest_path.exists():
        raise SystemExit(f"Missing freeze manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest["source_sha256"] != source_digests():
        raise SystemExit("Frozen source SHA-256 mismatch; code changed after preflight")
    res_path = frozen_dir / "results.json"
    if res_path.exists():
        raise SystemExit(f"Refusing to overwrite existing confirmatory {res_path}")

    p_sha = manifest["preflight_code_version"]
    anc = subprocess.run(["git", "merge-base", "--is-ancestor", p_sha, sha], cwd=ROOT)
    if anc.returncode != 0:
        raise SystemExit(f"Preflight commit {p_sha} is not an ancestor of {sha}")

    jobs = build_jobs()
    runs_by_key: dict[tuple[str, int], dict] = {}
    if ckpt_path.exists():
        saved = json.loads(ckpt_path.read_text(encoding="utf-8"))
        for row in saved.get("runs", []):
            runs_by_key[(row["arm"], int(row["seed"]))] = row
        print(f"Resuming S2 from checkpoint with {len(runs_by_key)}/{len(jobs)} done")

    ckpt_path.parent.mkdir(parents=True, exist_ok=True)
    t_all = time.perf_counter()
    for idx, job in enumerate(jobs):
        key = (job["arm"], int(job["seed"]))
        if key in runs_by_key:
            continue
        row = (
            _run_subspacing_one(job)
            if job["kind"] == "subspacing"
            else _run_weyl_one(job)
        )
        runs_by_key[key] = row
        ordered = [
            runs_by_key[(j["arm"], int(j["seed"]))]
            for j in jobs
            if (j["arm"], int(j["seed"])) in runs_by_key
        ]
        tmp = ckpt_path.with_suffix(".tmp")
        tmp.write_text(json.dumps({"runs": ordered}, indent=2), encoding="utf-8")
        tmp.replace(ckpt_path)

        if row["kind"] == "subspacing":
            print(
                f"[{idx + 1:3d}/{len(jobs)} {row['arm']:12s} s={row['seed']}] "
                f"Z={row['whitened_peak_z']:.2f}@{row['whitened_peak_frequency']:.0f} "
                f"qh={row['q_hat']:.3f} rms={row['calibration_rms']:.1e} "
                f"dt={row['total_seconds']:.1f}s",
                flush=True,
            )
        else:
            print(
                f"[{idx + 1:3d}/{len(jobs)} {row['arm']:12s} s={row['seed']}] "
                f"ratio={row['electric_ratio']:.4f} dt={row['total_seconds']:.1f}s",
                flush=True,
            )

    ordered_runs = [runs_by_key[(j["arm"], int(j["seed"]))] for j in jobs]
    summary = evaluate_gates(ordered_runs)
    dev = cp.cuda.runtime.getDeviceProperties(0)["name"].decode("utf-8", "replace")
    out = {
        "code_version": sha,
        "preflight_code_version": p_sha,
        "freeze_manifest_sha256": manifest["manifest_sha256"],
        "gpu_device": dev,
        "elapsed_seconds": time.perf_counter() - t_all,
        "summary": summary,
        "runs": ordered_runs,
    }
    res_path.write_text(json.dumps(out, indent=2), encoding="utf-8", newline="\n")
    print(f"Wrote S2 confirmatory results ({summary['verdict']}) to {res_path}")
    return res_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode", choices=("preflight", "campaign", "evaluate"), required=True
    )
    parser.add_argument("--frozen-dir", type=Path, default=DEFAULT_FROZEN_DIR)
    parser.add_argument("--ckpt", type=Path, default=DEFAULT_CKPT)
    args = parser.parse_args()

    if args.mode == "preflight":
        run_preflight(args.frozen_dir)
    elif args.mode == "campaign":
        run_campaign(args.frozen_dir, args.ckpt)
    else:
        data = json.loads((args.frozen_dir / "results.json").read_text("utf-8"))
        summary = evaluate_gates(data["runs"])
        print(
            json.dumps(
                {"verdict": summary["verdict"], "gates": summary["gates"]}, indent=2
            )
        )


if __name__ == "__main__":
    main()
