"""Stage B constraint (1), part one: the dS_2 diamond volume, exactly.

Section 5 forbids relying on the small-diamond expansion, because its
coefficient convention is "exactly the kind of unverified constant that
killed Stage C v1". So the relation is derived in closed form and
checked here against brute-force quadrature before anything is frozen.

DERIVATION (flat slicing, ds^2 = (ell^2/eta^2)(-d eta^2 + dx^2),
R = 2/ell^2, null coordinates U = eta + x, V = eta - x):

The causal diamond of a timelike pair p < q is the coordinate rectangle
U in [U_p, U_q], V in [V_p, V_q]. Its proper volume is

    Vol = int int Omega^2 d eta dx,  Omega^2 = ell^2/eta^2,
        = 2 ell^2 int int dU dV / (U+V)^2         [eta = (U+V)/2,
                                                  d eta dx = dU dV / 2]
        = 2 ell^2 ln[ (U_q+V_p)(U_p+V_q) / ((U_p+V_p)(U_q+V_q)) ].

With a = eta_p + eta_q and b = x_q - x_p, the numerator is a^2 - b^2 =
(eta_p+eta_q)^2 - (Delta x)^2 = 2 eta_p eta_q (Z + 1) where Z is the de
Sitter invariant of Section 3, and the denominator is 4 eta_p eta_q. So

    Vol = 2 ell^2 ln((Z+1)/2) = 4 ell^2 ln cosh(tau / (2 ell)),

using Z = cosh(tau/ell) and (cosh 2t + 1)/2 = cosh^2 t. Position
drops out entirely, as maximal symmetry requires, and the flat limit is

    Vol -> tau^2/2 - tau^4/(48 ell^2) = (tau^2/2)(1 - R tau^2/48),

so Section 5's "V(tau) = (tau^2/2)(1 - c R tau^2 + ...)" has c = 1/48
in this convention. The campaign uses the exact form; c is printed only
to retire the expansion's ambiguity explicitly.

THE INVERSION. Let g = Vol / tau^2. Then with s = tau/(2 ell),

    g = ln cosh(s) / s^2 =: G(s),

monotone decreasing from G(0+) = 1/2 (flat) to 0, so G is invertible on
(0, 1/2). Given an order-plus-count measurement of g the instrument
recovers s, hence ell = tau/(2s), hence

    R_hat = 2 / ell^2 = 8 s^2 / tau^2,   s = G^-1(g).

Run from the repository root. No experimental seeds are touched: this
file computes deterministic quadratures and contains no sprinkling.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "experiments/positive_control")
from p11_metric import _preflight_clean  # noqa: E402

#: Writes to the gitignored outputs/ tree, not straight into frozen/.
#: The preflight counts untracked files as dirty, so a script writing
#: into the repository would make the NEXT script's stamp dirty --
#: which is exactly what happened on the first attempt. Freezing is a
#: separate copy step, as it is for every runner in this programme.
OUT = Path("outputs/p12_stage_b_volume_check.json")
TOLERANCE = 1e-9          # Section 3's pinned tolerance, reused here

#: Review C12: the first version of this file used a midpoint rule at
#: n = 4000, measured a worst relative error of 4.1e-9 against a pinned
#: 1e-9, RECORDED it, and continued -- a guard that does not guard,
#: while the addendum called the relation verified at the tolerance. The
#: integrand 1/(U+V)^2 is smooth over the diamond, so Gauss-Legendre at
#: modest order reaches machine precision, and the tolerance is now
#: ASSERTED rather than reported.
GL_NODES = 128


def tau_curved(eta_p, x_p, eta_q, x_q, ell=1.0):
    z = ((eta_p ** 2 + eta_q ** 2 - (x_p - x_q) ** 2)
         / (2.0 * eta_p * eta_q))
    return ell * math.acosh(z)


def volume_closed_form(tau, ell=1.0):
    """4 ell^2 ln cosh(tau / 2 ell) -- the relation Stage B freezes."""

    return 4.0 * ell ** 2 * math.log(math.cosh(tau / (2.0 * ell)))


def volume_quadrature(eta_p, x_p, eta_q, x_q, ell=1.0, n=GL_NODES):
    """Integrate Omega^2 over the diamond in (U, V) with no use of the
    closed form and no use of tau. Tensor Gauss-Legendre; the integrand
    is smooth over the rectangle and bounded away from the lightcone
    tip because eta stays strictly negative, so this converges to
    machine precision at moderate order."""

    u_p, v_p = eta_p + x_p, eta_p - x_p
    u_q, v_q = eta_q + x_q, eta_q - x_q
    nodes, weights = np.polynomial.legendre.leggauss(n)
    u = 0.5 * (u_q - u_p) * nodes + 0.5 * (u_q + u_p)
    v = 0.5 * (v_q - v_p) * nodes + 0.5 * (v_q + v_p)
    wu = 0.5 * (u_q - u_p) * weights
    wv = 0.5 * (v_q - v_p) * weights
    eta = 0.5 * (u[:, None] + v[None, :])
    # Omega^2 d eta dx = Omega^2 * dU dV / 2
    return float(np.einsum("i,j,ij->", wu, wv,
                           ell ** 2 / eta ** 2) * 0.5)


def g_of_s(s):
    """G(s) = ln cosh(s) / s^2, with the s -> 0 limit filled in."""

    if s < 1e-8:
        return 0.5 - s ** 2 / 12.0
    return math.log(math.cosh(s)) / s ** 2


def invert_g(g):
    """s = G^-1(g) by bisection. G is strictly decreasing on (0, inf),
    so the inversion is unique wherever g < 1/2. Returns None when the
    measurement lands at or above the flat value, which is the failure
    mode Section 5 item 4 demands be visible: noise pushes g above 1/2
    and no positive curvature is representable."""

    if g >= 0.5:
        return None
    lo, hi = 1e-8, 1.0
    while g_of_s(hi) > g:
        hi *= 2.0
        if hi > 1e6:
            return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g_of_s(mid) > g:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main() -> None:
    results: dict = {
        "code_version": _preflight_clean(),
        "relation": "Vol = 4 ell^2 ln cosh(tau / (2 ell))",
        "expansion_coefficient_c": 1.0 / 48.0,
        "tolerance": TOLERANCE,
        "quadrature": f"tensor Gauss-Legendre, {GL_NODES} nodes per axis",
    }

    # (a) the closed form against quadrature, at positions chosen so the
    # SAME tau occurs at different places in the patch -- if the formula
    # were position-dependent this is where it would show.
    checks = []
    for eta_p, x_p, deta, dx in ((-1.9, 0.0, 0.6, 0.1),
                                 (-1.5, 0.7, 0.4, 0.05),
                                 (-1.2, -0.4, 0.15, 0.02),
                                 (-6.0, 2.0, 3.0, 0.5),
                                 (-3.0, -1.0, 1.2, 0.3)):
        eta_q, x_q = eta_p + deta, x_p + dx
        tau = tau_curved(eta_p, x_p, eta_q, x_q)
        closed = volume_closed_form(tau)
        quad = volume_quadrature(eta_p, x_p, eta_q, x_q)
        checks.append({
            "eta_p": eta_p, "x_p": x_p, "eta_q": eta_q, "x_q": x_q,
            "tau": tau, "closed_form": closed, "quadrature": quad,
            "rel_error": abs(closed - quad) / closed,
        })
        print(f"  eta {eta_p:+.2f}->{eta_q:+.2f} x {x_p:+.2f}->{x_q:+.2f}"
              f" | tau {tau:.4f} | closed {closed:.8f}"
              f" quad {quad:.8f} | rel {checks[-1]['rel_error']:.2e}")
    results["closed_form_vs_quadrature"] = checks
    results["worst_quadrature_rel_error"] = max(
        c["rel_error"] for c in checks)

    # (b) position independence at FIXED tau, stated as its own check
    target = 0.45
    same_tau = []
    for eta_p in (-1.9, -1.6, -1.3):
        # solve for deta at dx = 0 giving tau = target:
        # tau = ell * arccosh(1 + (deta^2)/(2 eta_p (eta_p+deta))) ...
        # easier: bisect on deta
        lo_d, hi_d = 1e-6, abs(eta_p) - 1e-6
        for _ in range(200):
            mid = 0.5 * (lo_d + hi_d)
            t = tau_curved(eta_p, 0.0, eta_p + mid, 0.0)
            if t < target:
                lo_d = mid
            else:
                hi_d = mid
        deta = 0.5 * (lo_d + hi_d)
        same_tau.append({
            "eta_p": eta_p, "deta": deta,
            "tau": tau_curved(eta_p, 0.0, eta_p + deta, 0.0),
            "quadrature": volume_quadrature(eta_p, 0.0, eta_p + deta, 0.0),
        })
    vols = [c["quadrature"] for c in same_tau]
    results["position_independence"] = {
        "target_tau": target, "samples": same_tau,
        "spread": (max(vols) - min(vols)) / float(np.mean(vols)),
        "closed_form": volume_closed_form(target),
    }
    vols_txt = " / ".join(f"{v:.8f}" for v in vols)
    print(f"  position independence at tau={target}: volumes"
          f" {vols_txt} | spread"
          f" {results['position_independence']['spread']:.2e}")

    # (c) the inversion is exact on exact inputs, and its amplification
    # is what Section 5 item 4 wants budgeted
    inv = []
    for tau_over_ell in (0.30, 0.60, 1.00, 1.50):
        tau = tau_over_ell            # ell = 1
        g_true = volume_closed_form(tau) / tau ** 2
        s = invert_g(g_true)
        r_hat = 8.0 * s ** 2 / tau ** 2
        # amplification: d(R_hat)/R_hat per unit relative error in g
        eps = 1e-4
        s_p = invert_g(g_true * (1.0 + eps))
        amp = abs((8.0 * s_p ** 2 / tau ** 2) / r_hat - 1.0) / eps
        inv.append({
            "tau_over_ell": tau_over_ell, "g_true": g_true,
            "one_minus_2g": 1.0 - 2.0 * g_true,
            "s_recovered": s, "R_hat": r_hat, "R_true": 2.0,
            "rel_error": abs(r_hat - 2.0) / 2.0,
            "amplification_dRhat_over_dg": amp,
        })
        print(f"  tau/ell={tau_over_ell:.2f}: g {g_true:.6f}"
              f" (1-2g {1.0 - 2.0 * g_true:.6f}) -> R_hat {r_hat:.9f}"
              f" | rel err {inv[-1]['rel_error']:.2e}"
              f" | amplification {amp:.1f}x")
    results["inversion"] = inv

    # Review C12: the tolerance is a GATE, not a field. Nothing is
    # written unless every check meets it.
    failures = [
        f"quadrature worst rel error "
        f"{results['worst_quadrature_rel_error']:.3e} > {TOLERANCE:.1e}"
    ] if results["worst_quadrature_rel_error"] > TOLERANCE else []
    if results["position_independence"]["spread"] > TOLERANCE:
        failures.append(
            f"position-independence spread "
            f"{results['position_independence']['spread']:.3e}"
            f" > {TOLERANCE:.1e}")
    worst_inv = max(c["rel_error"] for c in results["inversion"])
    if worst_inv > TOLERANCE:
        failures.append(f"inversion worst rel error {worst_inv:.3e}"
                        f" > {TOLERANCE:.1e}")
    results["gate_passed"] = not failures
    if failures:
        raise SystemExit("volume check FAILED the pinned tolerance: "
                         + "; ".join(failures))

    # Review C20: outputs/ is gitignored, so on a fresh checkout it does
    # not exist and write_text does not create parents -- the script
    # would finish its work and then raise FileNotFoundError.
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2))
    print(f"\ngate passed at {TOLERANCE:.0e}; wrote {OUT}")


if __name__ == "__main__":
    main()
