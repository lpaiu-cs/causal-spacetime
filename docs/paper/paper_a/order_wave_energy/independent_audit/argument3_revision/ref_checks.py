# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these
# three lines. Absolute paths point to the worktree it ran in.
"""Referee checks for gram_fourier_noise_r2.md. Deterministic, no sampling."""
import json, math, pathlib, sys
from fractions import Fraction as Fr
import numpy as np

ROOT = pathlib.Path(r"E:\lab\causal-spacetime\.claude\worktrees\upbeat-newton-3e92fe\docs\paper\paper_a\order_wave_energy")
sys.path.insert(0, str(ROOT / "independent_audit" / "proof3b_wave_signal"))
import wavegeom as wg  # noqa

CERT = json.loads((ROOT / "gram_fourier_evidence" / "anchor_guard_certificate.json").read_text())

print("== 1. a_2x a_2y supremum, exact rational-ish")
G11, G22, G12 = 0.00638, 0.00642, 2e-5
sup = G12 * math.sqrt(G22 / G11 - G12**2 / G11**2)
print("sup", sup)
# is it monotone in |G12| on [0,2e-5]? derivative of x sqrt(r - x^2/G11^2)
xs = np.linspace(0, 2e-5, 1001)
vals = xs * np.sqrt(G22 / G11 - xs**2 / G11**2)
print("monotone in |G12|:", bool(np.all(np.diff(vals) > 0)))

print("== 2. certificate reconstruction")
ell = CERT["ell_magnitude_upper"]; a3 = CERT["third_anchor_coordinate_lower"]
first = CERT["D_imaginary_lower"] + ell * 2e-5 / 0.348
print("first term", first)
for lab, e in (("cert ell 1.00943", ell), ("text ell 1.02", 1.02)):
    print(lab, "-> Im D_A lower", first - e * sup / 0.348)
print("grid coefficient", first - ell * sup / 0.348 - CERT["D_absolute_upper"] / 80)

print("== 3. alpha0 closed form and monotonicity, q in (0,1]")
for q in (1e-6, 0.01, 0.25, 1.0):
    bg = wg.Background(q)
    L = np.linspace(1e-3, 1.0, 2000)
    a = bg.alpha0(L)
    x = math.sqrt(q) * L
    closed = (1 - x / np.tan(x)) / (2 * q) if q > 1e-5 else L**2 / 6
    print(f"q={q}: max rel err closed form {np.max(np.abs(a/closed-1)):.1e}; increasing {bool(np.all(np.diff(a)>0))}")
    blow = CERT["anchor_pair_slack_lower"] ** 2 * float(bg.alpha0(0.096)) / bg.gamma0
    bpr = CERT["active_probe_slack_lower"] ** 2 * float(bg.alpha0(0.248)) / bg.gamma0
    print(f"   b_low {blow:.4e}  probe bound {bpr:.4e} ratio {bpr/blow:.1f}")

print("== 4. anchor-pair slack lower bound from the box alone (dm>=du)")
dv = Fr(96, 1000); zz = 2 * Fr(642, 100000) + Fr(4, 100000)
print("pair(1,2) slack >=", float(dv - zz / (2 * dv)), "cert", CERT["anchor_pair_slack_lower"])

print("== 5. states: eta, zeta of probes in supp w and of anchors on the outer box")
for q in (1e-6, 0.25, 1.0):
    bg = wg.Background(q)
    M = bg.M
    def ez(u, v, r2):
        m = bg.m(u); t = m / M
        eta = v - r2 / (2 * M * t)
        zeta = 1 - v - r2 / (2 * M * (1 - t))
        return eta, zeta
    U, V, R = np.meshgrid(np.linspace(0.1, 0.3, 41), np.linspace(0.65, 0.75, 21), np.linspace(0.0025, 0.04, 21))
    e, z = ez(U, V, R)
    print(f"q={q}: probe min eta {e.min():.4f} min zeta {z.min():.4f} max eta+zeta {(e+z).max():.4f}")
    mins = []
    for (u0, v0, r0) in ((-0.35, 0.15, 0.0064), (-0.25, 0.25, 0.0064), (-0.15, 0.35, 0.0128)):
        uu, vv, rr = np.meshgrid(u0 + np.array([-2e-3, 2e-3]), v0 + np.array([-2e-3, 2e-3]), r0 + np.array([-2e-5, 2e-5]))
        e, z = ez(uu, vv, rr)
        mins.append((round(float(e.min()), 4), round(float(z.min()), 4), round(float((e + z).max()), 4)))
    print("   anchors (min eta, min zeta, max eta+zeta):", mins)

print("== 6. A_reg is not compactly supported: H vs b_3X at a background config")
q = 1.0; bg = wg.Background(q)
# frozen anchors, d = 0.08
ua = np.array([-0.35, -0.25, -0.15]); va = np.array([0.15, 0.25, 0.35])
ya = np.array([[0.08, 0.0], [0.0, 0.08], [0.08, 0.08]])
ma = bg.m(ua)
GA = ya[:2] @ ya[:2].T
g3 = ya[:2] @ ya[2]
ellv = np.append(-np.linalg.solve(GA, g3), 1.0)
ux, vx, yx = 0.2, 0.7, np.array([0.05, -0.03]); mx = bg.m(ux)
r2x = yx @ yx
def gIX(i, b):
    eta0 = math.sqrt(b * bg.gamma0 / float(bg.alpha0(ux - ua[i])))
    return 0.5 * (ya[i] @ ya[i] + r2x) - (mx - ma[i]) * (vx - va[i] - eta0)
def slack(i):
    return vx - va[i] - np.sum((yx - ya[i])**2) / (2 * (mx - ma[i]))
b0 = [slack(i)**2 * float(bg.alpha0(ux - ua[i])) / bg.gamma0 for i in range(3)]
H0 = sum(ellv[i] * gIX(i, b0[i]) for i in range(3))
print("background b_iX", [f"{b:.4f}" for b in b0], "H at population", f"{H0:.2e}")
for b3 in (b0[2], 1.0, 10.0, 100.0):
    H = sum(ellv[i] * gIX(i, b0[i] if i < 2 else b3) for i in range(3))
    print(f"   b_3X={b3:8.3f}  H={H:.4f}")

print("== 7. |beta(s,t)| / (lam min(du, lam)) bounded")
for r in (64, 256, 1024):
    W = wg.Wave(1.0, r, 0.37)
    s = np.linspace(-0.5, 0.45, 60)
    worst = 0
    for du in (1e-4, 1e-3, 1e-2, 0.05, 0.2, 0.5):
        t = np.minimum(s + du, 0.5)
        ratio = np.abs(W.beta(s, t)) / (W.lam * np.minimum(t - s, W.lam))
        worst = max(worst, float(ratio.max()))
    print(f"r={r}: sup ratio {worst:.3f}")

print("== 8. elementary check: full exact variance, n=1, m=2, by quadrature")
x, w = np.polynomial.legendre.leggauss(200); x = 0.5 * (x + 1); w = 0.5 * w
for omega in (2 * math.pi, 7.3, 20.0):
    # cross term E_R |g(R)|^2, g(R) = int_0^1 (1[R<=t] - t) e^{-i omega t} dt, exact
    def g(R):
        return (np.exp(-1j*omega*R) - np.exp(-1j*omega)) / (1j*omega) - ((np.exp(-1j*omega)*(1+1j*omega) - 1) / (omega**2))
    # check closed form of int_0^1 t e^{-i w t}
    tint = np.sum(w * x * np.exp(-1j*omega*x))
    assert abs(tint - ((np.exp(-1j*omega)*(1+1j*omega) - 1)/omega**2)) < 1e-12
    cross = np.sum(w * np.abs(g(x))**2)
    diag = np.sum(w * x * (1 - x))
    n, m = 1, 2
    exact = (1 - 1/m) * cross / n + diag / (n * m)
    s_ = math.sin(omega/2)/(omega/2)
    formula = (1 - 1/m) * (1 - s_**2) / (n * omega**2) + 1 / (6 * n * m)
    print(f"omega={omega:.3f}: exact {exact:.10f} formula {formula:.10f}")

print("== 9. scan-range domination and exponent algebra")
for K in (1, 2, 3, 5):
    for n in (1e6, 1e12, 1e24):
        wmax = n ** (K / (2*K + 2))
        for om in (1.0, wmax ** 0.5, wmax):
            terms = [om**(k-2) * n**(-k/2) for k in range(1, K+1)]
            rem = om**K / n**((K+1)/2)
            ok = all(t <= max(1/(om*n**0.5), 1/n) * (1+1e-12) for t in terms) and rem <= 1/(om*n**0.5)*(1+1e-12)
            if not ok:
                print("FAIL", K, n, om)
    # K=1: is 1/n absorbed? sum has only k=1; remainder om/n >= 1/n
print("domination checked; K=1 sum has no k=2 term, 1/n <= omega/n remainder")
