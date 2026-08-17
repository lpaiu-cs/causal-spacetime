"""S6 mass-ladder rung geometry -- ONE exact path, general in M.

THE LADDER AXIS (PI ruling). The certified shell and the anchors stay
FIXED in absolute coordinates -- shell r in [10, 20], anchors at
r = 12 and r = 18, radially aligned -- and ONLY the mass changes.
Nothing is rescaled with M, so a rung is never an isometric copy: the
dimensionless compactness MU = 2M/r_c at the anchor midpoint
r_c = 15 genuinely differs between rungs, and with it the optical
curvature, the box shape and the diamond volume.

THE TIME WINDOW (PI ruling). dt is NOT free per rung: the frozen
convention is

    dt(M) = 8.5 * T_min(M) / T_min(1)      (KAPPA = 8.5 / T_min(1)),

the central rung's dimensionless slack carried to every rung. The
ratio form is used, not KAPPA * T_min(M): x / x == 1 exactly in
binary64, so dt(1.0) reproduces the executed central configuration's
8.5 EXACTLY, which a rounded KAPPA product would miss by an ulp.

ONE EXACT PATH. Every derived constant (box, B_box, SCALE) follows
`o4_sizing`'s certified path verbatim -- containment certificate,
outward binary64 widening, certified shell-volume enclosure, centre
-- only parametrized by M. `o4_sizing` itself is NOT modified: it is
part of the O4b recovery's frozen decision surface and stays the
M = 1 instance; this module must reproduce its SCALE bit-exactly at
M = 1, and the contract tests pin that.

WHAT IS AND IS NOT INHERITED (the certification addendum, S6
section of docs/theory/schwarzschild_volume_oracle_certification.md).
The patch-level lemmas are mass-GENERIC only under conditions that
must be re-certified per rung; `lemma_table` evaluates every one of
them as a certified interval comparison and refuses no rung silently:

    exterior            2M < R_MIN            (horizon below shell)
    K < 0 on shell      1.5M < R_MIN          (L3, negative optical
                                               curvature)
    w monotone          3M < R_MIN            (L1/L2c, photon sphere
                                               below shell)
    Q > 0 on arcs       perihelion floor > 3M (L2b via u_t <= 1/r_p)
    perihelion in shell reported vs R_MIN cos(1) (L2a patch bound)
    L4 margins          all four certainly positive (entry gate)
    L5 winding check    w_glob (2 pi - psi_max) > dt + w(r_hi) psi_max
                        (the review-supplied cross-check of the
                        alternative reduction, re-verified per rung)

M = 0 is OUT OF SCOPE here: its angle-cost constant is a different
formula (2 R_MIN / pi) and no S6 rung is massless.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

import certified_flight_time as cft  # noqa: E402
from certified_interval import Iv, float_down, float_up  # noqa: E402

R_IN, R_OUT = 12.0, 18.0
R_SHELL = (cft.R_MIN, cft.R_MAX)          # (10, 20), absolute
R_C = 0.5 * (R_IN + R_OUT)                # 15.0, the MU anchor point

#: The frozen ladder (PI ruling): central rung executed (O3'/O5),
#: two deeper rungs to freeze. MU = 2M/r_c is THE pre-frozen
#: dimensionless curvature indicator.
LADDER = (1.0, 1.4, 1.8)


def t_min(m: float) -> float:
    """Radial null flight time between the anchors at mass m:
    rho(R_OUT) - rho(R_IN) with rho(r) = r + 2M ln(r/2M - 1)."""

    if not 0.0 < m < R_IN / 2.0:
        raise ValueError(f"t_min: mass {m} outside (0, {R_IN / 2})")
    return (R_OUT - R_IN) + 2.0 * m * math.log(
        (R_OUT - 2.0 * m) / (R_IN - 2.0 * m))


#: KAPPA, reported for the record; dt() uses the ratio form.
KAPPA = 8.5 / t_min(1.0)


def dt(m: float) -> float:
    """The frozen time window: 8.5 * (T_min(M) / T_min(1)). Exactly
    8.5 at M = 1 (x / x == 1 in binary64)."""

    return 8.5 * (t_min(m) / t_min(1.0))


def mu(m: float) -> float:
    """The pre-frozen dimensionless curvature indicator: compactness
    2M/r_c at the anchor midpoint r_c = 15 (absolute)."""

    return 2.0 * m / R_C


def rung_geometry(m: float) -> dict:
    """The o4_sizing exact path, parametrized by mass: containment
    certificate -> outward-widened binary64 box -> certified
    shell-volume enclosure -> centres. Raises (CertificationError)
    if L4 does not certify."""

    window = dt(m)
    cert = cft.containment_certificate(R_IN, R_OUT, window, m=m)
    r_lo = float_down(cert["r_box"][0].lo)
    r_hi = float_up(cert["r_box"][1].hi)
    psi_max = float_up(cert["psi_max"].hi)
    cap = Iv(1) - Iv(psi_max).cos()
    cube = (Iv(r_hi) * Iv(r_hi) * Iv(r_hi)
            - Iv(r_lo) * Iv(r_lo) * Iv(r_lo))
    b_box_iv = Iv(2) * Iv.pi() * cap * cube / Iv(3)
    scale_iv = Iv(window) * b_box_iv

    def centre(iv: Iv) -> float:
        return 0.5 * (iv.lo_float() + iv.hi_float())

    return {
        "m": m, "dt": window, "t_min": t_min(m), "mu": mu(m),
        "r_lo": r_lo, "r_hi": r_hi, "psi_max": psi_max,
        "b_box": centre(b_box_iv), "scale": centre(scale_iv),
        "l_max_ub": float_up(cert["margins"]["nonempty"].hi),
        "margins": {k: v.lo_float()
                    for k, v in cert["margins"].items()},
        "perihelion_floor": cert["perihelion_floor"].lo_float(),
    }


def lemma_table(m: float) -> dict:
    """Every mass-genericity condition of the certification, evaluated
    as a certified comparison for THIS rung. `all_pass` gates any S6
    freeze; each row carries its margin so a failure is a coordinate,
    not a bare refusal."""

    g = rung_geometry(m)          # raises if L4 itself fails
    m_i = Iv(m)
    r_min_i = Iv(R_SHELL[0])
    rows = {}

    def row(name, lhs: Iv, rhs: Iv, statement: str):
        rows[name] = {
            "statement": statement,
            "certified": bool(lhs.certainly_gt(rhs)),
            "margin": (lhs - rhs).lo_float(),
        }

    row("exterior", r_min_i, Iv(2) * m_i,
        "horizon 2M below the shell floor R_MIN")
    row("k_negative", r_min_i, Iv(1.5) * m_i,
        "optical curvature K < 0 on the shell (r > 1.5M, L3)")
    row("w_monotone", r_min_i, Iv(3) * m_i,
        "photon sphere below the shell (r > 3M, L1/L2c)")
    row("q_positive", Iv(g["perihelion_floor"]), Iv(3) * m_i,
        "u_t <= 1/r_p < 1/(3M) on every oracle arc (L2b)")
    # L2a patch bound, reported with its own certified margin
    row("perihelion_patch", Iv(g["perihelion_floor"]),
        r_min_i * Iv(1.0).cos(),
        "box orbits stay above the patch perihelion bound "
        "R_MIN cos(1) (L2a)")
    # L5 winding cross-check of the alternative reduction: any
    # winding path costs at least w_glob (2 pi - psi_max); the direct
    # class is bounded by dt + w(r_hi) psi_max
    w_glob = Iv(3) * Iv(3).sqrt() * m_i
    w_hi = (Iv(g["r_hi"])
            / (Iv(1) - Iv(2) * m_i / Iv(g["r_hi"])).sqrt())
    winding_lb = w_glob * (Iv(2) * Iv.pi() - Iv(g["psi_max"]))
    direct_ub = Iv(g["dt"]) + w_hi * Iv(g["psi_max"])
    row("l5_winding", winding_lb, direct_ub,
        "winding paths certifiably longer than the direct class "
        "(review-supplied cross-check, per rung)")
    for name, margin in g["margins"].items():
        rows[f"l4_{name}"] = {
            "statement": f"L4 containment margin `{name}`",
            "certified": margin > 0.0,
            "margin": margin,
        }
    return {
        "m": m, "mu": g["mu"],
        "rows": rows,
        "all_pass": all(r["certified"] for r in rows.values()),
        "geometry": g,
    }


#: The frozen per-rung constants, PINNED as literals and re-derived
#: at import through the exact path above -- a drift in either
#: direction fails the import, so no consumer can read a stale or a
#: hand-rounded value. (The PI's independent check values: dt(1.4) =
#: 9.07056519..., SCALE(1.4) ~ 18,141.08; dt(1.8) = 9.72524817...,
#: SCALE(1.8) ~ 13,679.09.)
RUNG_CONSTANTS = {
    1.0: {"dt": 8.5, "scale": 28546.5704734038,
          "mu": 0.13333333333333333},
    1.4: {"dt": 9.070565190742672, "scale": 18141.08004658374,
          "mu": 0.18666666666666665},
    1.8: {"dt": 9.725248174609407, "scale": 13679.093767152488,
          "mu": 0.24000000000000002},
}

for _m, _want in RUNG_CONSTANTS.items():
    _g = rung_geometry(_m)
    if (_g["dt"] != _want["dt"] or _g["scale"] != _want["scale"]
            or _g["mu"] != _want["mu"]):
        raise AssertionError(
            f"s6_rungs: frozen constants for M={_m} drifted: "
            f"derived (dt={_g['dt']!r}, scale={_g['scale']!r}, "
            f"mu={_g['mu']!r}) != pinned {_want}")
del _m, _want, _g


def main() -> None:
    print(f"KAPPA = {KAPPA!r}   (dt(1.0) == 8.5: {dt(1.0) == 8.5})")
    for m in LADDER:
        tab = lemma_table(m)
        g = tab["geometry"]
        print(f"\nM = {m}: mu = {g['mu']:.6f}  dt = {g['dt']!r}")
        print(f"  box r in [{g['r_lo']!r}, {g['r_hi']!r}]  "
              f"psi_max {g['psi_max']!r}")
        print(f"  SCALE = {g['scale']!r}  all_pass = "
              f"{tab['all_pass']}")
        for name, r in tab["rows"].items():
            print(f"    {name:18s} {'PASS' if r['certified'] else 'FAIL'}"
                  f"  margin {r['margin']:+.6f}")


if __name__ == "__main__":
    main()
