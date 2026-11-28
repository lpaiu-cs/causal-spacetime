"""Deterministic checks for gram_fourier_noise_r2.md (argument 3 written out).

No sampling. Run: python checks.py
"""

from __future__ import annotations

import json
import math
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "proof3b_wave_signal"))
import wavegeom as wg  # noqa: E402

CERT = json.loads(
    (HERE.parents[1] / "gram_fourier_evidence" / "anchor_guard_certificate.json")
    .read_text()
)
failed = 0


def check(name, ok, detail=""):
    global failed
    print(("PASS " if ok else "FAIL ") + name + (" | " + detail if detail else ""))
    failed += 0 if ok else 1


# 1. Corrected Gram identity: a_2x a_2y = G12 sqrt(G22/G11 - G12^2/G11^2).
rng = np.random.default_rng(1)
worst = 0.0
for _ in range(1000):
    a = rng.normal(size=(2, 2))
    g = a @ a.T + 0.1 * np.eye(2)
    a1 = math.sqrt(g[0, 0])
    a2x = g[0, 1] / a1
    a2y = math.sqrt(g[1, 1] - a2x**2)
    formula = g[0, 1] * math.sqrt(g[1, 1] / g[0, 0] - (g[0, 1] / g[0, 0]) ** 2)
    worst = max(worst, abs(a2x * a2y - formula))
check("a_2x a_2y formula on 1000 random Gram matrices", worst < 1e-12, f"{worst:.1e}")

# 2. Its supremum on the stated box, G11, G22 in [0.00638, 0.00642], |G12| <= 2e-5.
sup = 2e-5 * math.sqrt(0.00642 / 0.00638 - (2e-5 / 0.00638) ** 2)
check("supremum of |a_2x a_2y| on the box", abs(sup - 2.0063e-5) < 1e-9, f"{sup:.6e}")

# 3. The certificate's Im D_A bound with the corrected supremum.
first = CERT["D_imaginary_lower"] + CERT["ell_magnitude_upper"] * 2e-5 / 0.348
corrected = first - CERT["ell_magnitude_upper"] * sup / 0.348
check(
    "certificate first term reproduced",
    abs(first - 0.0107560) < 5e-7,
    f"first term {first:.7f}",
)
check(
    "corrected Im D_A lower bound",
    abs(corrected - 0.0106978) < 5e-8 and corrected > 0.01,
    f"{corrected:.7f}, overstatement {CERT['D_imaginary_lower'] - corrected:.2e}",
)
grid = corrected - CERT["D_absolute_upper"] / 80
check("grid coefficient with the corrected bound", grid > 0.009, f"{grid:.7f}")

# 4. Positive b thresholds from the certified slacks: b = slack^2 alpha0(du) / gamma.
rows = []
for q in (1.0, 0.25, 0.01):
    bg = wg.Background(q)
    pair_slack, probe_slack = (
        CERT["anchor_pair_slack_lower"],
        CERT["active_probe_slack_lower"],
    )
    pair = pair_slack**2 * float(bg.alpha0(0.096)) / bg.gamma0
    probe = probe_slack**2 * float(bg.alpha0(0.248)) / bg.gamma0
    rows.append((q, pair, probe))
    check(
        f"b lower bounds positive at q={q}",
        pair > 0 and probe > 0,
        f"anchor pair {pair:.3e}, anchor to active probe {probe:.3e}",
    )
b_star = min(min(r[1], r[2]) for r in rows)
print(f"b_* (smallest over q in {{1, 0.25, 0.01}}) = {b_star:.3e}")

# 5. Elementary check for general omega, by deterministic quadrature.
#    S = (1/m) sum_i Hhat(X_i) e^{-i omega X_i}, Hhat = F_n - id on [0, 1].
#    Claim: E|S|^2 = (1 - 1/m)(1 - sinc^2(omega/2)) / (n omega^2) + 1/(6 n m).
x, w = np.polynomial.legendre.leggauss(400)
x, w = 0.5 * (x + 1), 0.5 * w


def kernel(r, omega):
    """g(R) = int_0^1 (1[R <= t] - t) e^{-i omega t} dt, the indicator part
    integrated exactly on [R, 1] by Gauss panels (no discontinuity inside)."""
    nodes = r + (1 - r) * x
    head = (1 - r) * np.sum(w * np.exp(-1j * omega * nodes))
    return head - np.sum(w * x * np.exp(-1j * omega * x))


for omega in (2 * math.pi, 7.3, 20.0):
    g = np.array([kernel(r, omega) for r in x])
    var_g = np.sum(w * np.abs(g - np.sum(w * g)) ** 2)
    sinc = math.sin(omega / 2) / (omega / 2)
    closed = (1 - sinc**2) / omega**2
    check(
        f"cross term at omega={omega:.3f}",
        abs(var_g - closed) < 2e-4 * closed + 1e-9,
        f"quadrature {var_g:.6e}, closed form {closed:.6e}",
    )

print(f"SUMMARY failed={failed}")
raise SystemExit(1 if failed else 0)
