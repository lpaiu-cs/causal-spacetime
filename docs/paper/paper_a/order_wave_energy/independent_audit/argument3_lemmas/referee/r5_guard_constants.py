# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of gram_fourier_lemmas.md.
"""Referee check r5: the anchor-guard constants used by Proposition L, re-derived independently.

Exact rationals with outward-rounded square roots (Fractions, 60-digit isqrt), on the outer
box: u_i, v_i within 0.002 of (-0.35,-0.25,-0.15), (0.15,0.25,0.35); Gram entries within 2e-5
of G* = [[.0064,0,.0064],[0,.0064,.0064],[.0064,.0064,.0128]].
Monotonicity arguments are written per quantity. Also: a dense floating scan of the box
corners/edges as a sanity check that the bounds are not violated, and the r2 number 0.0106978
recomputed with |ell_2| <= 1.009434 (not the rounded 1.0094) and m' <= 1.30.
"""
import itertools
import math
from fractions import Fraction as Fr

S = 10**60


def sq_lo(x):
    return Fr(math.isqrt(x.numerator * S * S // x.denominator), S)


def sq_hi(x):
    return sq_lo(x) + Fr(1, S)


t = Fr(2, 10**5)
du = Fr(2, 10**3)
g = Fr(64, 10**4)
G11, G22, G33, G13, G23 = (g - t, g + t), (g - t, g + t), (2 * g - t, 2 * g + t), (g - t, g + t), (g - t, g + t)
G12 = t  # |G12| <= t

# frame a1 = (sqrt G11, 0); a2 = (G12/sqrt G11, sqrt(G22 - G12^2/G11)); a3 = (G13/sqrt G11, (G23 - G12 G13/G11)/a2y)
a3x_lo = G13[0] / sq_hi(G11[1])  # increasing in G13, decreasing in G11
a3y_num_lo = G23[0] - G12 * G13[1] / G11[0]
a2y_hi = sq_hi(G22[1])  # sqrt(G22 - G12^2/G11) <= sqrt(G22)
a3y_lo = a3y_num_lo / a2y_hi
# ell_1 = -(G13 G22 - G23 G12)/det, ell_2 = -(G23 G11 - G13 G12)/det, det = G11 G22 - G12^2
det_lo = G11[0] * G22[0] - G12**2
# crude: numerator and denominator bounded separately
ell1_crude = (G13[1] * G22[1] + G23[1] * G12) / det_lo
ell2_crude = (G23[1] * G11[1] + G13[1] * G12) / det_lo
# sharper and still valid: G23 G11/(G11 G22 - G12^2) = G23/(G22 - G12^2/G11), worst at G11 low
ell2_hi = G23[1] / (G22[0] - G12**2 / G11[0]) + G13[1] * G12 / det_lo
ell1_hi = G13[1] / (G11[0] - G12**2 / G22[0]) + G23[1] * G12 / det_lo
a2xy_hi = G12 * sq_hi(G22[1] / G11[0])  # |G12| sqrt(G22/G11 - G12^2/G11^2) <= |G12| sqrt(G22/G11)
# sec^2(1/2) upper bound: cos(1/2) >= 1 - 1/8 + 1/384 - 1/46080
cos_lo = 1 - Fr(1, 8) + Fr(1, 384) - Fr(1, 46080)
sec2_hi = 1 / cos_lo**2
u_lo, u_hi = Fr(1, 10), Fr(3, 10)
d3_hi = sec2_hi * (u_hi - (Fr(-15, 100) - du))
d2_lo = u_lo - (Fr(-25, 100) + du)
d1_lo = u_lo - (Fr(-35, 100) + du)
d3_lo = u_lo - (Fr(-15, 100) + du)
imDA_lo = a3x_lo * a3y_lo / d3_hi - ell2_hi * a2xy_hi / d2_lo
absDA_hi = (ell1_hi * G11[1] / d1_lo + ell2_hi * G22[1] / d2_lo + G33[1] / d3_lo) / 2
grid_lo = imDA_lo - absDA_hi * Fr(1, 80)

out = []
out.append(f"a3x >= {float(a3x_lo):.10f}   a3y >= {float(a3y_lo):.10f}")
out.append(f"crude |ell| <= {float(ell2_crude):.7f}; ")
out.append(f"|ell1| <= {float(ell1_hi):.10f}   |ell2| <= {float(ell2_hi):.10f}   |a2x a2y| <= {float(a2xy_hi):.6e}")
out.append(f"sec^2(1/2) <= {float(sec2_hi):.10f}")
out.append(f"Im D_A >= {float(imDA_lo):.10f}   |D_A| <= {float(absDA_hi):.10f}   grid >= {float(grid_lo):.10f}")

# r2's chain with 1.30 and its quoted |ell_2| 1.0094, and with the true 1.009434
for ell2 in (Fr(10094, 10**4), ell2_hi, ell2_crude):
    v = a3x_lo * a3y_lo / (Fr(130, 100) * Fr(452, 1000)) - ell2 * a2xy_hi / d2_lo
    out.append(f"r2 chain (m' <= 1.30), |ell2| = {float(ell2):.7f}: Im D_A >= {float(v):.10f}  (>= 0.0106978: {v >= Fr(106978, 10**7)})")
v = Fr(106978, 10**7) - Fr(424, 10**4) * Fr(1, 80)
out.append(f"grid with 0.0106978 and 0.0424: {float(v):.8f}")

# slacks
pairs = {(1, 2): (G11[1] + G22[1] + 2 * G12, Fr(1, 10) - 2 * du, Fr(1, 10) - 2 * du),
         (2, 3): (G22[1] + G33[1] - 2 * G23[0], Fr(1, 10) - 2 * du, Fr(1, 10) - 2 * du),
         (1, 3): (G11[1] + G33[1] - 2 * G13[0], Fr(2, 10) - 2 * du, Fr(2, 10) - 2 * du)}
ps = {kk: dv - dz2 / (2 * dmm) for kk, (dz2, dv, dmm) in pairs.items()}
out.append("anchor-pair slack lower bounds: " + ", ".join(f"{kk}: {float(v):.8f}" for kk, v in ps.items()))
probe = []
for vi, ui, Gii in ((Fr(15, 100), Fr(-35, 100), G11[1]), (Fr(25, 100), Fr(-25, 100), G22[1]),
                    (Fr(35, 100), Fr(-15, 100), G33[1])):
    dz = (sq_hi(Fr(4, 100)) + sq_hi(Gii)) ** 2
    probe.append(Fr(65, 100) - (vi + du) - dz / (2 * (u_lo - (ui + du))))
out.append(f"active-probe slack >= {float(min(probe)):.8f}")


# floating scan: vertices of the Gram box x u-extremes; true D_A from an actual frame
def frame(G):
    a1 = (math.sqrt(G[0][0]), 0.0)
    a2x = G[0][1] / a1[0]
    a2y = math.sqrt(G[1][1] - a2x**2)
    a3x = G[0][2] / a1[0]
    a3y = (G[1][2] - a2x * a3x) / a2y
    return [a1, (a2x, a2y), (a3x, a3y)]


def mfun(u, q):
    k = math.sqrt(q)
    return (math.tan(k * u) + math.tan(k / 2)) / k


worst_im, worst_abs = 1.0, 0.0
for s in itertools.product((-1, 1), repeat=6):
    G = [[0.0064 + s[0] * 2e-5, s[1] * 2e-5, 0.0064 + s[2] * 2e-5],
         [0, 0.0064 + s[3] * 2e-5, 0.0064 + s[4] * 2e-5],
         [0, 0, 0.0128 + s[5] * 2e-5]]
    G[1][0], G[2][0], G[2][1] = G[0][1], G[0][2], G[1][2]
    det = G[0][0] * G[1][1] - G[0][1] ** 2
    l1 = -(G[0][2] * G[1][1] - G[1][2] * G[0][1]) / det
    l2 = -(G[1][2] * G[0][0] - G[0][2] * G[0][1]) / det
    ell = (l1, l2, 1.0)
    a = frame(G)
    for q in (1e-4, 0.25, 1.0):
        for uu in (-0.352, -0.348), (-0.252, -0.248), (-0.152, -0.148):
            pass
        for us in itertools.product((-0.352, -0.348), (-0.252, -0.248), (-0.152, -0.148)):
            for u in [0.1 + 0.2 * j / 40 for j in range(41)]:
                D = 0
                for i in range(3):
                    zc = complex(a[i][0], -a[i][1])
                    D += -0.5 * ell[i] * zc**2 / (mfun(u, q) - mfun(us[i], q))
                worst_im = min(worst_im, D.imag)
                worst_abs = max(worst_abs, abs(D))
out.append(f"float scan over box vertices, q in {{1e-4,0.25,1}}: min Im D_A = {worst_im:.6f}, max |D_A| = {worst_abs:.6f}")
print("\n".join(out))
