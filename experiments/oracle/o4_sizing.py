"""O4 sizing: every frozen constant on ONE path, and an ANALYTIC
worst-case certification of the G1 equivalence gate's power.

The review ruling fixed what this module must do:

  * derive `B_box`, `B_out`, `SCALE`, the Z ceiling and the sample
    sizes from the sampler's OWN boundaries in a single path -- no
    hand-copied geometry, so a boundary edit cannot silently
    desynchronise the sizing from the run;
  * pin `L_max` as an UPWARD bound of the certified margin, never a
    nearest-rounded regeneration (a downward L_max would understate
    the variance ceiling and inflate the power claim);
  * prove the worst case over the CONTINUUM `V_true in [V_lo, V_hi]`.
    A finite grid is not a certification, so the proof here is an
    interval-arithmetic branch and bound over the whole segment.

Why the power argument closes analytically at all (no distributional
model, per the ruling's choice of certification route (a)):

  L1. Z_i^2 <= c Z_i for Z_i in [0, c], hence
      V_n <= c * Zbar * n/(n-1) DETERMINISTICALLY. The empirical
      Bernstein half-width is therefore a monotone function of the
      observed mean and needs no separate concentration argument.
      The bound is nearly lossless here: (c - m)/c ~ 0.989.
  L2. With H(.) monotone, `pass` is the event that Chat lies in a
      fixed window [c_lo, c_hi], obtained by solving Chat +/- H(Chat)
      against the equivalence band. Using L1's UPPER bound for H
      shrinks the window, so P(Chat in window) is a LOWER bound on
      the true pass probability.
  L3. Var(Z) <= m (c - m) for Z in [0, c] with mean m (two-point
      extremal distribution) -- the only distributional input, and it
      is an inequality, not a model.
  L4. Bernstein's inequality bounds P(Chat outside the window). The
      normal approximation is reported alongside as a DIAGNOSTIC and
      is never the basis of the frozen size.

The Z ceiling `c = L_max/dt` holds under the defect-free alternative,
which is exactly the alternative the power is computed at. Coverage
never uses it: the empirical Bernstein interval only needs
`Z in [0, 1]`, i.e. `L <= dt`, which the runner asserts fail-closed.

Run:  python experiments/oracle/o4_sizing.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import certified_flight_time as cft  # noqa: E402
from certified_interval import Iv, float_down, float_up  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
_ARTIFACT = _REPO / "docs" / "prereg" / "p14_o4_sizing.json"

# --- the frozen decision constants (the gate itself) ------------------

TAU = 0.030                 # equivalence margin, relative to V_ref
DELTA_G1 = 0.02             # per side -> two-sided 96% EB interval
ALPHA_G2 = 0.005            # per end -> G2 spends 0.01 in total
ALPHA_G3 = 0.05             # one-sided cluster CP
LEAK_BUDGET = 0.0025        # G2 budget, as a fraction of V_ref
G3_CLUSTERS = 100_000
POWER_TARGET = 0.90

#: Frozen sample sizes. G1 is rounded UP from the smallest size whose
#: worst-case Bernstein power reaches 0.90 (26,012,722) so that the
#: continuum proof closes with a machine-checkable margin instead of
#: the 8e-7 headroom the exact threshold leaves.
N_G1 = 26_200_000
MS_PER_POINT = 1.536e-3     # 2 flight_time calls at 768 us each

# --- geometry, derived once from the certified containment box -------

_CERT = cft.containment_certificate(*cft.FROZEN_ANCHORS)
R_IN, R_OUT, DT = cft.FROZEN_ANCHORS

#: The sampling box: the certified box widened OUTWARD to binary64.
#: Sampling a slightly larger box is safe -- it still contains the
#: diamond, and the extra shell contributes L = 0, which costs a
#: little variance and no bias.
R_LO = float_down(_CERT["r_box"][0].lo)
R_HI = float_up(_CERT["r_box"][1].hi)
PSI_MAX = float_up(_CERT["psi_max"].hi)

#: L_max = dt - T_cert(P, Q), the certified ceiling on L. Taken as an
#: UPWARD bound of the certified margin: a larger ceiling only makes
#: the variance bound (L3) more conservative.
L_MAX_UB = float_up(_CERT["margins"]["nonempty"].hi)

#: The S1 frozen patch, for the G2 leakage stratum.
PATCH_R = (cft.R_MIN, cft.R_MAX)
PATCH_CAP = 1.0             # s1_schwarzschild_cost.CAP_HALF_ANGLE


def _shell_volume(r_lo: float, r_hi: float, psi: float) -> Iv:
    """2 pi (1 - cos psi) (r_hi^3 - r_lo^3) / 3, the r^2 sin(theta)
    measure of a polar-cap shell, as a certified enclosure."""

    cap = Iv(1) - Iv(psi).cos()
    cube = (Iv(r_hi) * Iv(r_hi) * Iv(r_hi)
            - Iv(r_lo) * Iv(r_lo) * Iv(r_lo))
    return Iv(2) * Iv.pi() * cap * cube / Iv(3)


B_BOX_IV = _shell_volume(R_LO, R_HI, PSI_MAX)
B_PATCH_IV = _shell_volume(PATCH_R[0], PATCH_R[1], PATCH_CAP)
B_OUT_IV = B_PATCH_IV - B_BOX_IV
SCALE_IV = Iv(DT) * B_BOX_IV          # V = SCALE * E[Z]

def _centre(iv: Iv) -> float:
    """The binary64 centre of a certified enclosure. These are known
    constants of the sampling region, not estimates: the enclosures
    are ~1e-25 wide, twenty orders below the statistical width, and
    the contract test pins each centre inside its own enclosure."""

    return 0.5 * (iv.lo_float() + iv.hi_float())


B_BOX = _centre(B_BOX_IV)
B_OUT = _centre(B_OUT_IV)
SCALE = _centre(SCALE_IV)
C_CEIL = L_MAX_UB / DT                # Z ceiling under no defect


def oracle_interval() -> tuple[float, float]:
    """[V_lo, V_hi] from the published O3 campaign artifact."""

    art = json.loads(
        (_REPO / "docs" / "prereg" / "p14_o3_volume.json")
        .read_text(encoding="utf-8"))["result"]
    return art["v_lo"], art["v_hi"]


V_LO, V_HI = oracle_interval()
V_REF = 0.5 * (V_LO + V_HI)           # frozen band scale
BAND = TAU * V_REF


# --- the accept window (lemmas L1 + L2) ------------------------------


def eb_half_width_bound(chat: float, n: int) -> float:
    """H-bar(Chat): the empirical Bernstein half-width in V units,
    with V_n replaced by its deterministic bound c*Zbar*n/(n-1)."""

    ln = math.log(2.0 / DELTA_G1)
    kappa = n / (n - 1.0)
    return (math.sqrt(2.0 * B_BOX * L_MAX_UB * chat * kappa * ln / n)
            + SCALE * 7.0 * ln / (3.0 * (n - 1.0)))


def accept_window(n: int) -> tuple[float, float]:
    """Counts of Chat for which the identified discrepancy interval
    [C_lo - V_hi, C_hi - V_lo] lies inside +/- BAND."""

    a, b = V_HI - BAND, V_LO + BAND
    lo, hi = 0.0, 4.0 * V_HI
    for _ in range(200):                      # Chat + H increasing
        mid = 0.5 * (lo + hi)
        if mid + eb_half_width_bound(mid, n) <= b:
            lo = mid
        else:
            hi = mid
    c_hi = lo
    lo, hi = 0.0, 4.0 * V_HI
    for _ in range(200):                      # Chat - H increasing
        mid = 0.5 * (lo + hi)
        if mid - eb_half_width_bound(mid, n) >= a:
            hi = mid
        else:
            lo = mid
    return hi, c_hi


# --- rigorous exp(-x) upper bound without an interval exp ------------

_TAYLOR_TERMS = 64


def exp_neg_upper(x: Iv) -> Iv:
    """An upper bound on exp(-x) for x >= 0, using only +, *, /.

    Every term of exp(x)'s Taylor series is positive for x >= 0, so a
    truncated sum is a valid LOWER bound on exp(x); its reciprocal is
    therefore a valid UPPER bound on exp(-x). This keeps the proof
    inside the certified interval module's existing operation set --
    `certified_interval.py` is digest-pinned by the O3 freeze and is
    not touched to add a transcendental."""

    lo = x.lo if x.lo > 0 else Iv(0).lo
    xl = Iv(lo)
    term, total = Iv(1), Iv(1)
    for k in range(1, _TAYLOR_TERMS + 1):
        term = term * xl / Iv(k)
        total = total + term
    return Iv(1) / Iv(total.lo)


def failure_upper(v: Iv, n: int, c_lo: float, c_hi: float) -> Iv:
    """Bernstein upper bound on P(Chat outside the accept window) for
    every V_true in the interval `v`, at the defect-free alternative.

    Both tails use the same sigma^2 enclosure, so the competition
    between them (the upper tail grows with v while the lower tail
    shrinks) is resolved by the enclosure itself rather than assumed
    away."""

    scale = SCALE_IV
    m = v / scale
    sigma2 = m * (Iv(C_CEIL) - m)                       # lemma L3
    two_c_over_3 = Iv(2) * Iv(C_CEIL) / Iv(3)
    out = Iv(0)
    for edge, sign in ((Iv(c_hi), 1), (Iv(c_lo), -1)):
        t = (edge - v) / scale if sign > 0 else (v - edge) / scale
        if not t.lo > 0:
            return Iv(1)                    # window not resolved here
        phi = Iv(n) * t * t / (Iv(2) * sigma2 + two_c_over_3 * t)
        out = out + exp_neg_upper(phi)
    return out


def certify_power(n: int, max_boxes: int = 2_000_000) -> dict:
    """Branch and bound proving sup_{V_true in [V_lo, V_hi]} failure
    <= 1 - POWER_TARGET over the CONTINUUM, not a grid."""

    c_lo, c_hi = accept_window(n)
    target = 1.0 - POWER_TARGET
    work = [(V_LO, V_HI)]
    boxes = worst = 0
    worst_ub = 0.0
    while work:
        lo, hi = work.pop()
        boxes += 1
        if boxes > max_boxes:
            raise RuntimeError("branch and bound did not converge")
        ub = failure_upper(Iv(lo, hi), n, c_lo, c_hi).hi_float()
        if ub <= target:
            worst_ub = max(worst_ub, ub)
            continue
        if hi - lo <= 1e-12:
            raise RuntimeError(
                f"cannot certify near V_true={lo!r}: bound {ub}")
        mid = 0.5 * (lo + hi)
        work.append((lo, mid))
        work.append((mid, hi))
        worst = worst + 1
    ends = {f"{name}": failure_upper(Iv(x, x), n, c_lo,
                                     c_hi).hi_float()
            for name, x in (("V_lo", V_LO), ("V_ref", V_REF),
                            ("V_hi", V_HI))}
    return {
        "n": n,
        "accept_window": [c_lo, c_hi],
        "eb_half_width_at_V_hi": eb_half_width_bound(V_HI, n),
        "failure_upper_bound": worst_ub,
        "power_lower_bound": 1.0 - worst_ub,
        "failure_at_endpoints": ends,
        "worst_endpoint": max(ends, key=ends.get),
        "boxes": boxes,
        "splits": worst,
    }


def normal_diagnostic(n: int) -> float:
    """Normal-approximation power -- DIAGNOSTIC ONLY, never the basis
    of a frozen size."""

    def phi(z: float) -> float:
        return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))

    c_lo, c_hi = accept_window(n)
    worst = 1.0
    for v in (V_LO, V_REF, V_HI):
        m = v / SCALE
        sd = math.sqrt(m * (C_CEIL - m) / n) * SCALE
        worst = min(worst, phi((c_hi - v) / sd) - phi((c_lo - v) / sd))
    return worst


# --- G2 and G3 sizing ------------------------------------------------


def g2_points(budget_frac: float = LEAK_BUDGET,
              alpha: float = ALPHA_G2) -> int:
    """Points needed so that a ZERO-count Clopper-Pearson upper bound
    certifies leaked volume <= budget.

    A zero count gives p <= 1 - alpha^(1/n) <= ln(1/alpha)/n, and
    V_leak <= dt * B_out * p because L_out <= dt."""

    return math.ceil(DT * B_OUT * math.log(1.0 / alpha)
                     / (budget_frac * V_REF))


def g3_upper(clusters: int = G3_CLUSTERS,
             alpha: float = ALPHA_G3) -> float:
    """One-sided CP upper bound on the per-CLUSTER mismatch rate from
    zero mismatching clusters. The cluster, not the call, is the
    independent unit: the two stress calls at one spatial point share
    that point's flight times."""

    return 1.0 - alpha ** (1.0 / clusters)


def summary() -> dict:
    n_g2 = g2_points()
    cert = certify_power(N_G1)
    calls = 2 * N_G1 + 2 * n_g2 + 4 * G3_CLUSTERS
    hours = (N_G1 + n_g2) * MS_PER_POINT / 3600 \
        + 4 * G3_CLUSTERS * 0.768e-3 / 3600
    return {
        "stage": ("O4 sizing: analytic worst-case certification of "
                  "the G1 equivalence gate, with G2/G3 companions"),
        "cost_convention": (
            "one flight_time call = 768 us (p14_s1_cost.json, tol "
            "1e-8); one spatial point = 2 calls = 1.536 ms; one G3 "
            "cluster = 4 causal_relation calls"),
        "geometry": {
            "source": "certified_flight_time.containment_certificate",
            "anchors": {"r_in": R_IN, "r_out": R_OUT, "dt": DT},
            "sampling_box": {"r_lo": R_LO, "r_hi": R_HI,
                             "psi_max": PSI_MAX},
            "patch": {"r_lo": PATCH_R[0], "r_hi": PATCH_R[1],
                      "cap_half_angle": PATCH_CAP},
            "B_box": B_BOX, "B_out": B_OUT, "scale_dt_B_box": SCALE,
            "L_max_upper_bound": L_MAX_UB,
            "Z_ceiling_c": C_CEIL,
        },
        "oracle": {"v_lo": V_LO, "v_hi": V_HI, "v_ref": V_REF,
                   "certified_half_width": (V_HI - V_LO)
                   / (V_HI + V_LO),
                   "tau_floor_full_width": (V_HI - V_LO) / V_REF},
        "gate": {"tau": TAU, "band_abs": BAND,
                 "delta_g1_per_side": DELTA_G1,
                 "alpha_g2_per_end": ALPHA_G2,
                 "alpha_g3": ALPHA_G3,
                 "leak_budget_frac": LEAK_BUDGET,
                 "familywise": DELTA_G1 * 2 + ALPHA_G2 * 2,
                 "total_sentence": TAU + LEAK_BUDGET},
        "g1": cert | {"power_normal_diagnostic": normal_diagnostic(
            N_G1)},
        "g2": {"points": n_g2,
               "leak_bound_abs": LEAK_BUDGET * V_REF},
        "g3": {"clusters": G3_CLUSTERS,
               "cp_upper_zero_mismatch": g3_upper()},
        "budget": {"calls": calls, "hours": hours,
                   "cap_calls": 80_000_000, "cap_wall_s": 86_400.0,
                   "cap_calls_ratio": 80_000_000 / calls,
                   "cap_wall_ratio": 86_400.0 / (hours * 3600)},
    }


def main() -> None:
    s = summary()
    print(f"box r in [{R_LO:.10f}, {R_HI:.10f}], "
          f"psi <= {PSI_MAX:.12f}")
    print(f"B_box = {B_BOX:.6f}   B_out = {B_OUT:.6f}   "
          f"scale = {SCALE:.4f}")
    print(f"L_max upper bound = {L_MAX_UB:.13f}  ->  c = {C_CEIL:.12f}")
    print(f"V_ref = {V_REF:.6f}  band = {BAND:.6f}  "
          f"tau floor (full width) = "
          f"{s['oracle']['tau_floor_full_width']:.6%}")
    c = s["g1"]
    print(f"\nG1 n = {c['n']:,}")
    print(f"  accept window = [{c['accept_window'][0]:.7f}, "
          f"{c['accept_window'][1]:.7f}]")
    print(f"  EB half-width at V_hi = "
          f"{c['eb_half_width_at_V_hi']:.7f}")
    print(f"  CONTINUUM proof: {c['boxes']} boxes, "
          f"{c['splits']} splits")
    print(f"  failure upper bound = {c['failure_upper_bound']:.8f}"
          f"  ->  power >= {c['power_lower_bound']:.8f}")
    print(f"  endpoints: {c['failure_at_endpoints']}")
    print(f"  worst endpoint = {c['worst_endpoint']}")
    print(f"  normal diagnostic (NOT a basis) = "
          f"{c['power_normal_diagnostic']:.6f}")
    print(f"\nG2 points = {s['g2']['points']:,}  "
          f"(leak <= {s['g2']['leak_bound_abs']:.6f})")
    print(f"G3 clusters = {s['g3']['clusters']:,}  "
          f"CP upper = {s['g3']['cp_upper_zero_mismatch']:.6e}")
    b = s["budget"]
    print(f"\ntotal {b['calls']:,} calls, {b['hours']:.2f} h  "
          f"| caps {b['cap_calls']:,} ({b['cap_calls_ratio']:.2f}x) "
          f"/ 24 h ({b['cap_wall_ratio']:.2f}x)")
    with open(_ARTIFACT, "w", encoding="utf-8", newline="\n") as f:
        json.dump(s, f, ensure_ascii=False, indent=1)
        f.write("\n")
    print(f"artifact: {_ARTIFACT}")


if __name__ == "__main__":
    main()
