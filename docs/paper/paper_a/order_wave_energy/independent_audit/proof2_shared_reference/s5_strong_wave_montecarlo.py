# ruff: noqa: E501, UP031
"""Monte Carlo validation of the EXACT wave moment formulas used in Sections 3-4.

The population inputs of the proof,

    p = alpha_lambda(a,u) eta^2 / gamma_lambda,  f = alpha_lambda(u,b) zeta^2 / gamma_lambda,
    c = Kbar C_3 / (3 C_2^2)   (pair probability of a sub-diamond),

are checked against brute-force sampling of the exact wave at LARGE amplitude
(sinh(2 beta) = 2 lambda sqrt(q) of order one), where the wave and background
formulas differ by many standard errors.  The sampler uses only the pairwise
causal criterion  dv >= dy^T Mt^{-1} dy / 2  and the density a(u)^2; it does not
use the volume/chain hierarchy.  Also checks that the sampled relation is a
strict partial order (the nestedness of Section 1 needs transitivity under BOTH
hypotheses).

q = 6 is outside the theorem's range 0 < q <= 1; it is used only because the
pair formula is far from the flat value 1/10 there.

Deterministic (fixed seeds).  Prints PASS/FAIL lines; exits nonzero on any FAIL.
Run:  python s5_strong_wave_montecarlo.py      (system Python: numpy only; ~1.5 min)
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_geometry import Wave, p0  # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ((" | " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)


def volume_check(tag, q, lam, n, seed):
    wave = Wave(q, lam, phase=0.37)
    back = Wave(q, lam, background=True)
    event = np.array([0.10, 0.55, 0.15, -0.10])
    theta = wave.moments(event, factor=3.0)
    theta_bg = back.moments(event)
    rng = np.random.default_rng(seed)
    pts = wave.sample(n, rng)
    p_hat = float(wave.precedes(pts, event).mean())
    f_hat = float(wave.follows(pts, event).mean())
    zp = (p_hat - theta[0]) / math.sqrt(theta[0] * (1 - theta[0]) / n)
    zf = (f_hat - theta[1]) / math.sqrt(theta[1] * (1 - theta[1]) / n)
    zp_bg = (p_hat - theta_bg[0]) / math.sqrt(theta_bg[0] * (1 - theta_bg[0]) / n)
    zf_bg = (f_hat - theta_bg[1]) / math.sqrt(theta_bg[1] * (1 - theta_bg[1]) / n)
    check(
        f"{tag} q={q}, lambda={lam} (eps={wave.eps:.3f}): MC volumes match the wave formulas",
        abs(zp) < 4.5 and abs(zf) < 4.5 and max(abs(zp_bg), abs(zf_bg)) > 10,
        f"p {theta[0]:.6f} vs MC {p_hat:.6f} (z={zp:.2f}); f {theta[1]:.6f} vs MC {f_hat:.6f} "
        f"(z={zf:.2f}); background formula at the same event is off by z={zp_bg:.1f}, {zf_bg:.1f}",
    )
    return wave, pts


def order_check(tag, wave, pts):
    sub = pts[:500]
    n = len(sub)
    rel = np.zeros((n, n), dtype=bool)
    for i in range(n):
        rel[i] = wave.follows(sub, sub[i])  # rel[i, j]: sub[i] < sub[j]
    antisym = not np.any(rel & rel.T)
    two_step = (rel.astype(np.int32) @ rel.astype(np.int32)) > 0
    violations = int(np.count_nonzero(two_step & ~rel))
    check(
        f"{tag} sampled exact relation is a strict partial order (500 events)",
        antisym and violations == 0,
        f"relations {int(rel.sum())}, transitivity violations {violations}",
    )


def pair_check():
    q, lam = 6.0, 0.3
    wave = Wave(q, lam, phase=0.37)
    formula = wave.pair_probability(-0.5, 0.5, factor=3.0)
    background = p0(q)
    values = []
    for seed in range(6):
        rng = np.random.default_rng(700 + seed)
        pts = wave.sample(20000, rng)
        values.append(wave.comparable_pairs(pts) / math.comb(len(pts), 2))
    values = np.array(values)
    mean = float(values.mean())
    se = float(values.std(ddof=1) / math.sqrt(len(values)))
    check(
        "W3 q=6, lambda=0.3: MC pair probability matches Kbar C3/(3 C2^2), not the background",
        abs(mean - formula) < 0.0025 and abs(mean - background) > 0.0035,
        f"wave formula {formula:.5f}, background p0 {background:.5f}, "
        f"MC {mean:.5f} (empirical se {se:.5f}, 6 x 20000 events)",
    )


if __name__ == "__main__":
    w1, pts1 = volume_check("W1", 1.0, 0.5, 4000000, 11)
    order_check("W1b", w1, pts1)
    w2, pts2 = volume_check("W2", 6.0, 0.3, 3000000, 12)
    order_check("W2b", w2, pts2)
    small = Wave(1.0, 1.0 / 64, phase=0.37)
    order_check(
        "W4 q=1, lambda=1/64:", small, small.sample(500, np.random.default_rng(13))
    )
    pair_check()
    print(f"SUMMARY s5_strong_wave_montecarlo: {len(FAILS)} FAIL")
    sys.exit(1 if FAILS else 0)
