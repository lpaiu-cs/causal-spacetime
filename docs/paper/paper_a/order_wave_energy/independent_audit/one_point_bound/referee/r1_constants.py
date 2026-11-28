# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of one_point_strip_bound.md.
"""Referee r1: recompute the constants table of one_point_strip_bound.md
independently (scipy.quad on the defining integrals, m from quad of f, no closed
forms used for W or Vhat), and the closed forms of Step 7 for comparison.
Also: the lambda -> 0 limit L(q) of TV/lambda^2 from a leading-order expansion
(derived by the referee, see report), to compare with the note's numerics.
"""
import math
import numpy as np
from scipy.integrate import quad

QS = (0.25, 1.0, 4.0, 9.0, 9.8)


def run(q):
    s = math.sqrt(q)
    f = lambda u: 1.0 / math.cos(s * u) ** 2
    a2 = lambda u: math.cos(s * u) ** 2
    fmax = f(0.5)
    CE = 2 * s * (4 * fmax - 2)
    d = (math.pi / 2 - s / 2) / s
    pts = [0.5 - d * k for k in (0.5, 0.25, 0.1, 0.05)] + [-0.5 + d * k for k in (0.5, 0.25, 0.1, 0.05)]
    pts = sorted(p for p in pts if -0.5 < p < 0.5)

    def mm(u):
        return quad(f, -0.5, u, epsabs=0, epsrel=1e-13, limit=400)[0]

    def mp(u):
        return quad(f, u, 0.5, epsabs=0, epsrel=1e-13, limit=400)[0]

    M = mm(0.5)
    vhat = quad(lambda u: a2(u) * mp(u) * mm(u) / M, -0.5, 0.5, points=pts, epsabs=0, epsrel=1e-11, limit=400)[0]
    W = quad(lambda u: a2(u) * (mp(u) ** 2 + mm(u) ** 2) / M ** 2, -0.5, 0.5, points=pts, epsabs=0, epsrel=1e-11, limit=400)[0]
    # closed forms of Step 7
    Th = math.tan(s / 2)
    Ic = 0.5 + math.sin(s) / (2 * s)
    Is = 0.5 - math.sin(s) / (2 * s)
    vhat_cf = (Th * Th * Ic - Is) / (2 * s * Th)
    W_cf = (Is + Th * Th * Ic) / (2 * Th * Th)
    t1 = 2 * s / math.sqrt(1 + 4 * q)
    r1 = (1 + t1) / (1 - t1)
    C10 = CE * W / vhat
    C1 = r1 * C10
    # leading-order limit L(q) = (2 s / (pi Vhat)) int a^2 <|f(u) - fmax e^{i chi}|>_chi (m_+^2+m_-^2)/M^2 du
    chi = (np.arange(4096) + 0.5) * 2 * math.pi / 4096

    def avgabs(u):
        return float(np.mean(np.abs(f(u) - fmax * np.exp(1j * chi))))

    Lint = quad(lambda u: a2(u) * avgabs(u) * (mp(u) ** 2 + mm(u) ** 2) / M ** 2, -0.5, 0.5, points=pts + [0.0], epsabs=0, epsrel=1e-9, limit=400)[0]
    L = 2 * s * Lint / (math.pi * vhat)
    # tighter constant: (B-A)^+ <= e/(m(m+e)) <= e/m^2 gives 1/(1-t) instead of (1+t)/(1-t)
    print(
        f"q={q}: f_max={fmax:.6g} C_E={CE:.6g} M={M:.12g} (cf {2*Th/s:.12g}) "
        f"Vhat={vhat:.12g} (cf {vhat_cf:.12g}) W={W:.12g} (cf {W_cf:.12g}) "
        f"W/Vhat={W/vhat:.6g} C_1^0={C10:.6g} t_1={t1:.6f} (1+t1)/(1-t1)={r1:.6g} C_1={C1:.6g} "
        f"| L(q)=lim TV/lam^2={L:.6g}  C_1^0/L={C10/L:.4g}  dist_to_pole={d:.4g}"
    )


for q in QS:
    run(q)
