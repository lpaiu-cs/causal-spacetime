"""Isolated diagnostic for campaign v2's bottom-rung offset.

Section 11.2 left one question live. The whole of `Delta_13 = +0.0243`
is a `tau/ell = 0.30` offset -- curved rungs 2, 3 and 4 agree to 0.003
dex -- and the FLAT twin, whose true contrast is zero by construction,
carries the same offset at `+0.0179 +/- 0.0136`. Of the three
candidates 11.2 lists, the `m`-distribution one is already closed by
the frozen artifact (`sd_m` flat to 2%). This script separates the
remaining two.

Post hoc, quarantined, labelled: seeds live in the design-check space
`9000000+` that Section 10 reserves, never in an experimental window,
and nothing here re-reads a gate. Run from the repository root.

THREE PARTS
-----------
D1, control template at precision. The flat arm at all four rungs,
1000 samples at the two rungs that carry the step and 500 at the
plateau. The campaign measured the twin at `n = 179`, giving a step
standard error of 0.0136; at 1000 it is about 0.0058, which resolves
a `+0.018` offset at roughly 3 sigma. This decides whether the flat
arm's offset is real at all.

D2, curved template at matching precision. Same sizes on the curved
arm, so the difference-of-differences that 11.2 could only quote at
1.1 sigma gets a standard error near 0.008 instead of 0.017.

D3, the geometry separation, which is the actual discriminator. The
suspected quantity is packing pressure, `K` box areas over the
patch's `(u,v)` area, which the frozen constants leave UNMATCHED
across the ladder at 0.168 / 0.222 / 0.263 / 0.263. Rung 0.30 is the
loosest. So move that ratio at FIXED `tau/ell` and see whether `y`
follows it:

    V-LOOSE   rung 0.60, X 4.0 -> 5.28, ratio 0.222 -> 0.168
              (given rung 0.30's packing pressure)
    V-TIGHT   rung 0.30, X 1.8 -> 1.36, ratio 0.168 -> 0.222
              (given rung 0.60's packing pressure)

Both variants keep `eta_lo`, the band, `K`, the eligibility pair, the
rejection cap and the estimator exactly as frozen; only `X` moves, and
`rho` is recalibrated so the realized mean `m` still matches the
campaign's grand mean 75.755. Without that recalibration the variant
would trade a curvature question for a discreteness one, which is the
confound Section 5 gate 1 exists to prevent.

PREDICTIONS, STATED BEFORE THE RUN
----------------------------------
Under candidate 2 (packing pressure drives the offset):
  V-LOOSE reads like rung 0.30, near -0.761, NOT like rung 0.60's
  -0.737; and V-TIGHT reads like rung 0.60, near -0.737, NOT like
  rung 0.30's -0.761. The two variants CROSS.
Under candidate 3 (a real curvature effect, not of (tau/ell)^2 form):
  each variant stays at its own rung's value -- V-LOOSE near -0.737,
  V-TIGHT near -0.761 -- because `X` is not a curvature parameter.
The two candidates therefore predict opposite orderings of the same
two numbers, on both arms, which is why this is worth 4 minutes.

A third possibility is that neither moves cleanly and the offset is
some other property of the bottom rung's constants. The script reports
completion rate and realized `m` per configuration so that outcome is
visible rather than silent.
"""

from __future__ import annotations

import json
import sys
from contextlib import contextmanager
from pathlib import Path

import numpy as np

sys.path.insert(0, "experiments/positive_control")
import p13_tau_ell as P  # noqa: E402

OUT = Path("docs/prereg/frozen/p13v2/p13v2_rung1_diag_results.json")

#: Design-check space, Section 10: 9000000+. Blocks are 400000 apart
#: and no block consumes more than 1000 * STRIDE = 200000 seeds, so
#: they are pairwise disjoint. Every experimental window ever used
#: sits below 8800000 (v1 6000000-6475999, v2 8000000-8795999) and
#: the v1 design space is 7000000-7999999.
D1_BASE = {0.30: 9_000_000, 0.60: 9_400_000,
           1.00: 9_800_000, 1.50: 10_200_000}
D2_BASE = {0.30: 10_600_000, 0.60: 11_000_000,
           1.00: 11_400_000, 1.50: 11_800_000}
D3_BASE = {("V-LOOSE", False): 12_200_000, ("V-LOOSE", True): 12_600_000,
           ("V-TIGHT", False): 13_000_000, ("V-TIGHT", True): 13_400_000}
CAL_BASE = {("V-LOOSE", False): 13_800_000, ("V-LOOSE", True): 14_000_000,
            ("V-TIGHT", False): 14_200_000, ("V-TIGHT", True): 14_400_000}

N_TEMPLATE = {0.30: 1000, 0.60: 1000, 1.00: 500, 1.50: 500}
N_VARIANT = 500
N_CALIBRATE = 150
M_TARGET_DIAG = 75.755081300813      # campaign v2's grand mean
M_CAL_TOLERANCE = 0.01
CAL_ITERATIONS = 4
#: A 150-sample probe at STRIDE = 200 consumes 30000 seeds, so
#: successive calibration iterations must be spaced by more than that
#: or they would re-measure the same ensembles; 4 * 40000 also fits
#: inside the 200000 gap between CAL_BASE blocks.
CAL_SPACING = 40_000

#: The two variants: (rung, new X). The X values solve
#: packing_ratio(rung, X) = the other rung's frozen ratio.
VARIANTS = {"V-LOOSE": (0.60, 5.28), "V-TIGHT": (0.30, 1.36)}

LOWEST_DESIGN_SEED = 9_000_000


def packing_ratio(tau_c: float, xhalf: float) -> float:
    """K box areas over the patch's (u,v) area. The box side is the
    rung's measured mean sqrt(dU dV) -- the same quantity Section 5
    uses to centre the twin's band -- and the patch's (u,v) area is
    twice its (eta, x) area because the map has Jacobian 2."""

    eta_lo = P.PATCH[tau_c][0]
    box_area = P.TWIN_BAND_CENTRE[tau_c] ** 2
    patch_uv = 4.0 * (abs(eta_lo) - 1.0) * xhalf
    return P.K_PAIRS * box_area / patch_uv


@contextmanager
def patched_rung(tau_c: float, xhalf: float, rho: float, rho_twin: float):
    """Move X and both intensities for one rung, then restore. The
    band, eta_lo, K, both eligibility conditions, the rejection cap
    and the estimator are untouched, so this is the campaign's code
    path with one geometric parameter moved."""

    saved_patch, saved_twin = P.PATCH[tau_c], P.TWIN_RHO[tau_c]
    eta_lo = saved_patch[0]
    P.PATCH[tau_c] = (eta_lo, xhalf, rho)
    P.TWIN_RHO[tau_c] = rho_twin
    try:
        yield
    finally:
        P.PATCH[tau_c] = saved_patch
        P.TWIN_RHO[tau_c] = saved_twin


def block(tau_c: float, base: int, count: int, flat: bool) -> dict:
    ys, ms, complete = [], [], 0
    for k in range(count):
        record, ok = P.run_sample(tau_c, base + P.STRIDE * k, flat=flat)
        if ok:
            complete += 1
            ys.append(record["y"])
            ms.append(record["mean_m"])
    y = np.asarray(ys)
    return {
        "n_attempted": count, "n_complete": complete,
        "completion": complete / count,
        "mean_y": float(y.mean()), "sd_y": float(y.std(ddof=1)),
        "se_y": float(y.std(ddof=1) / np.sqrt(y.size)),
        "mean_m": float(np.mean(ms)),
    }


def calibrate_rho(tau_c: float, xhalf: float, rho0: float,
                  base: int, flat: bool) -> dict:
    """Bring the realized mean m back to the campaign's grand mean
    after X moves. m is a count in a box whose size the band fixes, so
    it is very nearly proportional to rho and one or two Newton steps
    suffice; the loop records every step rather than just the answer."""

    rho, steps = rho0, []
    for iteration in range(CAL_ITERATIONS):
        curved_rho = rho if not flat else P.PATCH[tau_c][2]
        twin_rho = rho if flat else P.TWIN_RHO[tau_c]
        with patched_rung(tau_c, xhalf, curved_rho, twin_rho):
            probe = block(tau_c, base + CAL_SPACING * iteration,
                          N_CALIBRATE, flat)
        steps.append({"iteration": iteration, "rho": float(rho),
                      "mean_m": probe["mean_m"],
                      "completion": probe["completion"]})
        if abs(probe["mean_m"] / M_TARGET_DIAG - 1.0) <= M_CAL_TOLERANCE:
            return {"rho": float(rho), "converged": True, "steps": steps}
        rho = rho * M_TARGET_DIAG / probe["mean_m"]
    return {"rho": float(rho), "converged": False, "steps": steps}


def contrast(low: dict, high: dict) -> dict:
    delta = high["mean_y"] - low["mean_y"]
    se = float(np.sqrt(low["se_y"] ** 2 + high["se_y"] ** 2))
    return {"delta": float(delta), "se": se,
            "ci95": [float(delta - 1.96 * se), float(delta + 1.96 * se)],
            "sigma": float(delta / se) if se > 0 else float("nan")}


def main() -> None:
    for base in (list(D1_BASE.values()) + list(D2_BASE.values())
                 + list(D3_BASE.values()) + list(CAL_BASE.values())):
        assert base >= LOWEST_DESIGN_SEED, (
            f"seed block {base} is below the design-check space")

    results: dict = {
        "code_version": P._preflight_clean(),
        "purpose": "campaign v2 bottom-rung offset: control template, "
                   "curved template, and the packing-pressure variants",
        "quarantine": "design-check seeds 9000000+, no experimental "
                      "window touched, no gate re-read",
        "m_target": M_TARGET_DIAG,
        "frozen_packing_ratio": {
            str(t): packing_ratio(t, P.PATCH[t][1]) for t in P.RUNGS},
    }

    print("=== D1: flat template (the control at precision) ===",
          flush=True)
    results["d1_flat_template"] = {}
    for tau_c in P.RUNGS:
        out = block(tau_c, D1_BASE[tau_c], N_TEMPLATE[tau_c], flat=True)
        results["d1_flat_template"][str(tau_c)] = out
        print(f"  tau/ell={tau_c}: n {out['n_complete']:4d} | mean_y "
              f"{out['mean_y']:+.4f} +/- {out['se_y']:.4f} | m "
              f"{out['mean_m']:.2f} | completion {out['completion']:.3f}",
              flush=True)

    print("=== D2: curved template at matching precision ===", flush=True)
    results["d2_curved_template"] = {}
    for tau_c in P.RUNGS:
        out = block(tau_c, D2_BASE[tau_c], N_TEMPLATE[tau_c], flat=False)
        results["d2_curved_template"][str(tau_c)] = out
        print(f"  tau/ell={tau_c}: n {out['n_complete']:4d} | mean_y "
              f"{out['mean_y']:+.4f} +/- {out['se_y']:.4f} | m "
              f"{out['mean_m']:.2f} | completion {out['completion']:.3f}",
              flush=True)

    flat, curved = results["d1_flat_template"], results["d2_curved_template"]
    results["template_contrasts"] = {
        "flat_step_030_060": contrast(flat["0.3"], flat["0.6"]),
        "curved_step_030_060": contrast(curved["0.3"], curved["0.6"]),
        "flat_endpoint": contrast(flat["0.3"], flat["1.5"]),
        "curved_endpoint": contrast(curved["0.3"], curved["1.5"]),
    }
    tc = results["template_contrasts"]
    for name in ("step_030_060", "endpoint"):
        c_curved, c_flat = tc[f"curved_{name}"], tc[f"flat_{name}"]
        did = c_curved["delta"] - c_flat["delta"]
        se = float(np.sqrt(c_curved["se"] ** 2 + c_flat["se"] ** 2))
        tc[f"did_{name}"] = {
            "delta": float(did), "se": se,
            "ci95": [float(did - 1.96 * se), float(did + 1.96 * se)],
            "sigma": float(did / se)}
        print(f"  {name}: curved {c_curved['delta']:+.4f}"
              f" +/- {c_curved['se']:.4f} | flat {c_flat['delta']:+.4f}"
              f" +/- {c_flat['se']:.4f} | DiD {did:+.4f} +/- {se:.4f}",
              flush=True)

    print("=== D3: packing-pressure variants at fixed tau/ell ===",
          flush=True)
    results["d3_variants"] = {}
    for name, (tau_c, xhalf) in VARIANTS.items():
        eta_lo, x_frozen, rho_frozen = P.PATCH[tau_c]
        entry = {
            "rung": tau_c, "x_frozen": x_frozen, "x_variant": xhalf,
            "eta_lo": eta_lo,
            "packing_ratio_frozen": packing_ratio(tau_c, x_frozen),
            "packing_ratio_variant": packing_ratio(tau_c, xhalf),
            "arms": {},
        }
        for flat in (False, True):
            arm = "flat" if flat else "curved"
            rho0 = P.TWIN_RHO[tau_c] if flat else rho_frozen
            cal = calibrate_rho(tau_c, xhalf, rho0,
                                CAL_BASE[(name, flat)], flat)
            curved_rho = P.PATCH[tau_c][2] if flat else cal["rho"]
            twin_rho = cal["rho"] if flat else P.TWIN_RHO[tau_c]
            with patched_rung(tau_c, xhalf, curved_rho, twin_rho):
                out = block(tau_c, D3_BASE[(name, flat)],
                            N_VARIANT, flat=flat)
            out["calibration"] = cal
            entry["arms"][arm] = out
            print(f"  {name} {arm}: rho {cal['rho']:.1f}"
                  f" (converged {cal['converged']}) | n"
                  f" {out['n_complete']:4d} | mean_y {out['mean_y']:+.4f}"
                  f" +/- {out['se_y']:.4f} | m {out['mean_m']:.2f}"
                  f" | completion {out['completion']:.3f}", flush=True)
        results["d3_variants"][name] = entry

    #: The discriminator, evaluated against the predictions in the
    #: docstring: does the variant move toward the rung whose packing
    #: pressure it was given, or stay at its own rung?
    results["d3_readings"] = {}
    for name, (tau_c, _xhalf) in VARIANTS.items():
        other = 0.30 if tau_c == 0.60 else 0.60
        for flat in (False, True):
            arm = "flat" if flat else "curved"
            template = (results["d1_flat_template"] if flat
                        else results["d2_curved_template"])
            variant = results["d3_variants"][name]["arms"][arm]
            own, target = template[str(tau_c)], template[str(other)]
            gap = target["mean_y"] - own["mean_y"]
            moved = variant["mean_y"] - own["mean_y"]
            results["d3_readings"][f"{name}-{arm}"] = {
                "own_rung": tau_c, "own_rung_mean_y": own["mean_y"],
                "pressure_donor_rung": other,
                "donor_mean_y": target["mean_y"],
                "variant_mean_y": variant["mean_y"],
                "rung_gap": float(gap), "variant_shift": float(moved),
                "fraction_of_gap_traversed":
                    float(moved / gap) if gap != 0 else float("nan"),
                "shift_se": float(np.sqrt(variant["se_y"] ** 2
                                          + own["se_y"] ** 2)),
            }
            r = results["d3_readings"][f"{name}-{arm}"]
            print(f"  {name}-{arm}: shift {r['variant_shift']:+.4f}"
                  f" +/- {r['shift_se']:.4f} of a {gap:+.4f} rung gap"
                  f" -> {r['fraction_of_gap_traversed']:+.2f} of the way",
                  flush=True)

    OUT.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
