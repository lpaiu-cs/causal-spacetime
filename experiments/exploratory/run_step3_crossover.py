"""Exploratory local-Weyl-jet vs subspacing crossover benchmark.

Evaluates both exact continuum volumes (CPU) and finite-N causal order relations
(GPU) across wave frequency r = 1/lambda for:
  1. Fixed coarse-grained energy parameter q (circular_wave vs background control),
  2. Fixed crest Weyl norm C0 (exact Ricci-flat family).

Purely exploratory seeds (760000+); disjoint from preregistered confirmatory seeds.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exploratory"))
sys.path.insert(0, str(ROOT / "experiments" / "positive_control"))

import intrinsic_weyl_joint as iwj  # noqa: E402
from subspacing_observer import observe  # noqa: E402
from subspacing_wave_sim import OrderOracle, WaveGeometry  # noqa: E402


def exact_diamond_matrix(geom, landmarks, n_quad=64):
    """Exact 1D Gauss-Legendre diamond volumes for landmarks (u, v, y0, y1)."""
    m = len(landmarks)
    nodes, weights = np.polynomial.legendre.leggauss(n_quad)
    vols = np.zeros((m, m), dtype=float)
    for i in range(m):
        u1, v1, y1 = landmarks[i, 0], landmarks[i, 1], landmarks[i, 2:]
        for j in range(i + 1, m):
            u2, v2, y2 = landmarks[j, 0], landmarks[j, 1], landmarks[j, 2:]
            du = u2 - u1
            if du <= 1e-12:
                continue
            m12 = geom.matrix(u1, u2)
            dy = y2 - y1
            gap = (v2 - v1) - 0.5 * float(dy @ np.linalg.solve(m12, dy))
            if gap <= 1e-12:
                continue
            u = 0.5 * (u1 + u2) + 0.5 * du * nodes
            left = np.linalg.inv(geom.matrix(np.full_like(u, u1), u))
            right = np.linalg.inv(geom.matrix(u, np.full_like(u, u2)))
            det_total = np.linalg.det(left + right)
            integrand = (np.cos(geom.k * u) ** 2) / np.sqrt(
                np.maximum(det_total, 1e-300)
            )
            vols[i, j] = vols[j, i] = float(
                np.pi * gap**2 * 0.5 * du * np.sum(weights * integrand)
            )
    return vols


def _continuum_fit(geom, sel_seed=2026):
    t0 = time.perf_counter()
    interior = geom.sprinkle(2048, np.random.default_rng(sel_seed), "cpu")
    pts = np.vstack(([-0.5, 0.0, 0.0, 0.0], interior, [0.5, 1.0, 0.0, 0.0]))
    oracle = OrderOracle(geom, pts, "cpu")
    all_idx = np.arange(len(pts))
    ord_mat = oracle.query(all_idx, all_idx)[0]
    ids, mask = iwj.selection(ord_mat, design="two")
    lm = pts[ids]
    vols = exact_diamond_matrix(geom, lm, n_quad=64)
    fit = iwj.fit_joint(
        1e9 * vols / geom.volume, clock_mask=mask, starts=1, model_order=5
    )
    e2 = float(fit["electric_norm_squared"])
    c_norm = float(np.linalg.norm(fit["coefficients"]))
    g_norm = float(np.linalg.norm(fit["gradient_coefficients"]))
    rms_clk = float(np.sqrt(2.0 * fit["curved_cost"] / fit["measured_intervals"]))
    return {
        "landmarks": int(fit["landmarks"]),
        "clocks": int(fit["measured_intervals"]),
        "c_norm": c_norm,
        "g_norm": g_norm,
        "e2_hat": e2,
        "sqrt_e2_hat": float(np.sqrt(max(e2, 0.0))),
        "rms_clock": rms_clk,
        "optimizer_success": bool(fit["optimizer_success"]),
        "seconds": time.perf_counter() - t0,
    }


def run_continuum():
    results = {"fixed_q": [], "fixed_crest_c0": []}
    for q in (0.01, 0.25):
        for r in (0.10, 0.25, 0.50, 1.0, 2.0, 4.0, 8.0, 16.0, 24.0, 32.0, 48.0, 64.0):
            g = WaveGeometry(q, "circular_wave", 1.0 / r)
            fit = _continuum_fit(g)
            amp = float(np.sqrt(q * (r**2 + 4 * q)))
            fit.update(
                q=q,
                mode="circular_wave",
                r=r,
                amp_true=amp,
                electric_ratio=float(fit["sqrt_e2_hat"] / (np.sqrt(2.0) * amp)),
            )
            results["fixed_q"].append(fit)
            print(
                f"[cont fixed_q q={q:g} r={r:5.2f}] ratio={fit['electric_ratio']:.4f} "
                f"||c||={fit['c_norm']:.4f} dt={fit['seconds']:.2f}s",
                flush=True,
            )
        bg = WaveGeometry(q, "background", 1.0)
        fit_bg = _continuum_fit(bg)
        fit_bg.update(q=q, mode="background", r=None, amp_true=0.0, electric_ratio=None)
        results["fixed_q"].append(fit_bg)
        print(
            f"[cont fixed_q q={q:g} BG] ||c||={fit_bg['c_norm']:.4f} "
            f"sqrt(E2)={fit_bg['sqrt_e2_hat']:.4f}",
            flush=True,
        )

    c0 = 0.5
    for r in (0.50, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0):
        q = float((-(r**2) + np.sqrt(r**4 + 16 * c0**2)) / 8.0)
        g = WaveGeometry(q, "circular_wave", 1.0 / r)
        fit = _continuum_fit(g)
        fit.update(
            c0=c0,
            q=q,
            mode="circular_wave",
            r=r,
            amp_true=c0,
            electric_ratio=float(fit["sqrt_e2_hat"] / (np.sqrt(2.0) * c0)),
        )
        results["fixed_crest_c0"].append(fit)
        print(
            f"[cont fixed_C0={c0} r={r:5.2f} q={q:.5f}] "
            f"ratio={fit['electric_ratio']:.4f} dt={fit['seconds']:.2f}s",
            flush=True,
        )
    return results


def _finite_weyl(q, mode, r, n=2**21, seed=760000):
    import cupy as cp

    t0 = time.perf_counter()
    geom = WaveGeometry(q, mode, 1.0 / r)
    rng = np.random.default_rng(seed)
    sel_int = geom.sprinkle(2048, np.random.default_rng(seed + 10000), "cpu")
    sel_pts = np.vstack(([-0.5, 0.0, 0.0, 0.0], sel_int, [0.5, 1.0, 0.0, 0.0]))
    sel_ord = OrderOracle(geom, sel_pts, "cpu").query(
        np.arange(len(sel_pts)), np.arange(len(sel_pts))
    )[0]
    ids, mask = iwj.selection(sel_ord, design="two")
    lm = sel_pts[ids]
    pts = geom.sprinkle(n, rng, "gpu")
    # Strict Tip[0] -> Interior -> Tip[-1] order using GPU OrderOracle
    all_pts = np.vstack((lm, pts))
    oracle = OrderOracle(geom, all_pts, "gpu")
    d_lm = np.arange(len(lm))
    d_ev = np.arange(len(lm), len(all_pts))
    # Accumulate causal interval counts I(lm[i], lm[j]) tile by tile.
    counts = np.zeros((len(lm), len(lm)), dtype=np.float64)
    tile = 32768
    for s in range(0, n, tile):
        sub = d_ev[s : s + tile]
        fw, bw = oracle.query(d_lm, sub)
        counts += fw.astype(np.float64) @ bw.T.astype(np.float64)
    lm_ord = sel_ord[np.ix_(ids, ids)]
    causal = np.where(lm_ord, counts, 0.0)
    sym_counts = causal + causal.T
    fit = iwj.fit_joint(sym_counts, clock_mask=mask, starts=1, model_order=5)
    e2 = float(fit["electric_norm_squared"])
    del oracle
    cp.get_default_memory_pool().free_all_blocks()
    return {
        "n": n,
        "q": q,
        "mode": mode,
        "r": r,
        "seed": seed,
        "sqrt_e2_hat": float(np.sqrt(max(e2, 0.0))),
        "c_norm": float(np.linalg.norm(fit["coefficients"])),
        "seconds": time.perf_counter() - t0,
    }


def _finite_subspacing(q, mode, r, n=2**22, seed=761000, max_probes=65536):
    import cupy as cp

    t0 = time.perf_counter()
    geom = WaveGeometry(q, mode, 1.0 / r)
    pts = geom.sprinkle(n, np.random.default_rng(seed), "gpu")
    oracle = OrderOracle(geom, pts, "gpu")
    pilot_n = min(16384, n // 4)
    pilot = np.arange(pilot_n)
    probes = np.arange(pilot_n, n // 2)
    refs = np.arange(n // 2, n)
    cache = {}

    def counts(rf, at, early=(), late=(), _oracle=oracle):
        key = tuple(
            hashlib.blake2b(
                np.asarray(g, dtype=np.int64).tobytes(), digest_size=16
            ).digest()
            for g in (rf, at, early, late)
        )
        if key not in cache:
            cache[key] = _oracle.counts(rf, at, early, late)
        return cache[key]

    freqs_wide = np.arange(8.0, 257.0, 4.0)
    res = observe(
        oracle.query,
        pilot,
        refs,
        probes,
        freqs_wide,
        q=None,
        count_query=counts,
        max_probes=max_probes,
    )
    mask_std = freqs_wide >= 24.0
    w_std = res["whitened_amplitude"][mask_std]
    f_std = freqs_wide[mask_std]
    pk_std = int(np.argmax(w_std))
    del oracle, pts
    cp.get_default_memory_pool().free_all_blocks()
    return {
        "n": n,
        "q": q,
        "mode": mode,
        "r": r,
        "seed": seed,
        "q_hat": float(res["q_hat"]),
        "calibration_rms": float(res["calibration_rms"]),
        "wide_peak_freq": float(res["whitened_peak_frequency"]),
        "wide_peak_z": float(res["whitened_peak_z"]),
        "std_peak_freq": float(f_std[pk_std]),
        "std_peak_z": float(w_std[pk_std]),
        "c_ref_estimate": float(res["c_ref_estimate"]),
        "seconds": time.perf_counter() - t0,
    }


def run_finite():
    rows = {"weyl_finite": [], "subspacing_finite": []}
    # Clean finite-N comparison at q=0.25 and q=1.0 across overlap band
    for idx, r in enumerate((1.0, 2.0, 4.0, 8.0, 16.0, 32.0)):
        out = _finite_weyl(0.25, "circular_wave", r, n=2**21, seed=760000 + idx)
        amp = float(np.sqrt(0.25 * (r**2 + 1.0)))
        out["electric_ratio"] = float(out["sqrt_e2_hat"] / (np.sqrt(2.0) * amp))
        rows["weyl_finite"].append(out)
        print(
            f"[finite weyl q=0.25 r={r:4.1f}] ratio={out['electric_ratio']:.3f} "
            f"dt={out['seconds']:.1f}s",
            flush=True,
        )
    out_bg = _finite_weyl(0.25, "background", 1.0, n=2**21, seed=760099)
    rows["weyl_finite"].append(out_bg)

    for q, seed_base in ((0.25, 761000), (1.0, 761100)):
        for idx, r in enumerate((4.0, 8.0, 16.0, 24.0, 32.0, 48.0, 64.0, 96.0, 128.0)):
            out = _finite_subspacing(
                q, "circular_wave", r, n=2**22, seed=seed_base + idx
            )
            rows["subspacing_finite"].append(out)
            print(
                f"[finite subsp q={q:g} r={r:5.1f}] "
                f"wide={out['wide_peak_z']:.2f}@{out['wide_peak_freq']:.0f} "
                f"std={out['std_peak_z']:.2f}@{out['std_peak_freq']:.0f} "
                f"qh={out['q_hat']:.3f} dt={out['seconds']:.1f}s",
                flush=True,
            )
        bg_sub = _finite_subspacing(q, "background", 1.0, n=2**22, seed=seed_base + 99)
        rows["subspacing_finite"].append(bg_sub)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode", choices=("continuum", "finite", "all"), default="continuum"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "outputs" / "exploratory_step3_crossover.json",
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = {}
    if args.output.exists():
        payload = json.loads(args.output.read_text(encoding="utf-8"))
    if args.mode in ("continuum", "all"):
        payload["continuum"] = run_continuum()
    if args.mode in ("finite", "all"):
        payload["finite"] = run_finite()
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
