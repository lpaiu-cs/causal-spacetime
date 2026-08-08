"""P14 §8 P2: the volume prediction, under the exact error model.

EXPLORATORY. No gate, no threshold, no verdict, no confirmation seed
window appears in this file, and none may be added to it (design
§8.2). The `equivalent` / `discriminates` labels below are properties
of confidence intervals, never program verdicts.

The question: do the curved arm's anchor-diamond cardinalities
reproduce the §4.4 quadrature `r = V_A/V_0`? Scope is AXIS anchors
only -- the one diamond with an independent analytic prediction
(off-axis curved volumes are known only through MC, which would make
the check MC-vs-MC circular).

Protocol, frozen in the P2 design review (2026-08-08, in-session,
three rounds):

- **Primary** -- equivalence, per the accuracy house rule (AGENTS.md:
  interval, 95%, two-sided). Normalized error `theta = E[Z] / lam_0`
  with `Z_i = N_A,i - r N_0,i`; the Student-t 95% two-sided CI on
  `theta` must lie inside `[-tau, +tau]` with the chosen effect-scale
  margin `tau = delta = r - 1`. The licensed sentence is "the error
  from the quadrature prediction is smaller than the flat-to-
  prediction gap" -- no coefficient-grade recovery is claimed.
- **Secondary** -- discrimination: the same CI excludes the flat-truth
  value `theta = -delta`.
- **Marginal GOF per arm** -- the rate ratio (observed pooled count /
  predicted pooled mean) with its EXACT Garwood Poisson 95% CI must
  lie inside `[1 - tau_m, 1 + tau_m]`, `tau_m = delta`. This is the
  check that catches common-mode error, which `Z` is exactly blind to
  (P1 review R7.1); it is sized separately and BINDS the sprinkling
  counts at slice-a1.0.
- **Var(Z) diagnostic** -- the empirical between-sprinkling variance
  (bootstrap percentile CI) against the analytic `r rho V_dis` (CP
  bracket from a 2e6-sample MC). Reported, never judged: it is a
  cross-check of the error model, not an estimand.
- **Coverage** -- pointwise 95% per operating point; any joint
  sentence over the three points must be Bonferroni-adjusted and
  reported separately.

Error model (§8 P2): the anchors are EXTERNAL -- the fattest eligible
axis diamond of the frozen geometry, chosen before any realization,
never present among the sprinkled points -- and `sprinkle` draws
`N ~ Poisson(rho V_box)`, so `N_A` and `N_0` are exactly Poisson over
deterministic regions and the covariance identity holds as written.
Post-selecting anchors from a realization would void all of it.

Seeds: consumption unit is ONE rng stream per operating point (the
stream's first consumer is the campaign's first sprinkling); the MC
cross-check and the bootstrap derive from `MC_SEED` sub-streams,
disjoint from every campaign stream. Burned seeds (P1 campaign,
design checks) may not reappear. `assert_seed_layout` enforces the
whole layout and a test pins it.

Sizing: n = 21_000 / 400 / 100 from the design review -- primary TOST
at the V_dis 95% UCB (2e6 MC), marginal equivalence at `tau_m = delta`
(the binding constraint at slice-a1.0 and high-a2.0), and a floor of
100 for the variance estimate itself. Verified coverage/power of the
frozen t-CI by compound-Poisson simulation: 0.9495/0.9492/0.9473
coverage, 1.0000 equivalence power. The P1 inputs behind delta and
sd_Z live in `docs/prereg/p14_probe_p1_sizing.json` and are
cross-checked at run time.

Run:  python experiments/positive_control/p14_probe_p2.py
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
from p14_plane_wave import (
    PlaneWaveGeometry,
    Slab,
    arms,
    class_c_eligible,
    conjugate_du,
    sprinkle,
)
from p14_probe_p1 import (
    _between,  # the ONE membership predicate: counts and MC volumes
    axis_volume_ratio,
    diamond_volumes_mc,
    fattest_axis_diamond,
)

# ---------------------------------------------------------------------
# Frozen protocol constants (design review, 2026-08-08)
# ---------------------------------------------------------------------

#: (label, w, du, dv, dx, dy, n_sprinklings). n from the review's
#: sizing: max(primary TOST at the V_dis UCB, marginal tau_m = delta,
#: variance floor 100) -> 20_939 / 349 / 13 exact minima, frozen with
#: margin as 21_000 / 400 / 100.
OPERATING_POINTS = (
    ("slice-a1.0", 1.0, 1.0, 0.2, 6.0, 6.0, 21_000),
    ("high-a2.0", 1.0, 2.0, 0.5, 2.0, 2.0, 400),
    ("edge-a2.4", 1.0, 2.4, 0.2, 1.2, 0.8, 100),
)

#: Expected elements per box -- E[N_box], not a fixed count: sprinkle
#: draws N ~ Poisson(rho V_box), which is what makes the marginal laws
#: exactly Poisson. External-anchor counting is O(N), so N is free
#: where P1's O(N^2) pair census pinned it at 300.
N_BOX = 30_000

#: MC samples for the Var(Z) analytic cross-check -- 10x P1's 2e5 so
#: V_dis resolves at every point (120/79/737 hits at the design seed).
MC_SAMPLES = 2_000_000

#: One campaign stream per operating point; the MC / bootstrap streams
#: derive from MC_SEED sub-streams. 20260808 (P1 campaign) and 777
#: (design checks) are burned and may not reappear.
CAMPAIGN_SEEDS = {
    "slice-a1.0": 20260811,
    "high-a2.0": 20260812,
    "edge-a2.4": 20260813,
}
MC_SEED = 20260814
BURNED_SEEDS = (20260808, 777)

#: Bootstrap replicates for the Var(Z) percentile CI (diagnostic).
B_BOOT = 4_000

_T_LEVEL = 0.975  # two-sided 95%: the house-rule interval

#: Bonferroni per-point level for the ONE joint sentence over the
#: three operating points (protocol: any joint claim is adjusted and
#: reported separately; the primary record stays pointwise).
_T_LEVEL_BONF3 = 1.0 - 0.05 / 6.0


def assert_seed_layout() -> None:
    """The frozen seed layout: campaign seeds and MC_SEED pairwise
    distinct, and none of them burned. Runs at campaign start; a test
    pins the literals so an edit here cannot pass silently."""

    seeds = [*CAMPAIGN_SEEDS.values(), MC_SEED]
    assert len(set(seeds)) == len(seeds), "seed collision in the layout"
    clash = set(seeds) & set(BURNED_SEEDS)
    assert not clash, f"burned seeds reused: {sorted(clash)}"


# ---------------------------------------------------------------------
# Exact distributions, dependency-free (the repo already hand-builds
# Clopper-Pearson the same way; scipy is not a dependency)
# ---------------------------------------------------------------------


def _beta_cf(a: float, b: float, x: float) -> float:
    """Continued fraction for the regularized incomplete beta
    (modified Lentz), the standard construction."""

    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-15:
            return h
    raise RuntimeError("incomplete beta did not converge")


def _beta_inc(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta `I_x(a, b)`."""

    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    ln_front = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                + a * math.log(x) + b * math.log1p(-x))
    front = math.exp(ln_front)
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _beta_cf(a, b, x) / a
    return 1.0 - front * _beta_cf(b, a, 1.0 - x) / b


def student_t_cdf(t: float, df: int) -> float:
    """Student-t CDF via the incomplete beta."""

    x = df / (df + t * t)
    p = 0.5 * _beta_inc(0.5 * df, 0.5, x)
    return 1.0 - p if t >= 0.0 else p


def student_t_crit(df: int, level: float = _T_LEVEL) -> float:
    """`t` with `CDF(t) = level` by bisection -- the frozen CI's
    critical value (95% two-sided at `level = 0.975`)."""

    lo, hi = 0.0, 200.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if student_t_cdf(mid, df) < level:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _gamma_p(a: float, x: float) -> float:
    """Regularized lower incomplete gamma `P(a, x)`: series for
    `x < a + 1`, continued fraction otherwise (both standard)."""

    if x <= 0.0:
        return 0.0
    if x < a + 1.0:
        term = 1.0 / a
        total = term
        ap = a
        for _ in range(10_000):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * 1e-16:
                break
        else:
            raise RuntimeError("gamma series did not converge")
        return total * math.exp(-x + a * math.log(x) - math.lgamma(a))
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 10_000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-15:
            break
    else:
        raise RuntimeError("gamma continued fraction did not converge")
    return 1.0 - h * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gamma_quantile(a: float, p: float) -> float:
    """`x` with `P(a, x) = p` by bisection (a >= 1 here)."""

    lo = 0.0
    hi = a + 20.0 * math.sqrt(a) + 20.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _gamma_p(a, mid) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def poisson_mean_ci(total: int, level: float = 0.95) -> tuple[float, float]:
    """Exact (Garwood) two-sided CI for a Poisson mean given an
    observed `total`: `[qgamma(alpha/2, T), qgamma(1 - alpha/2,
    T + 1)]`, zero lower at zero counts -- the frozen marginal CI."""

    half = (1.0 - level) / 2.0
    lower = 0.0 if total == 0 else _gamma_quantile(float(total), half)
    upper = _gamma_quantile(float(total + 1), 1.0 - half)
    return lower, upper


# ---------------------------------------------------------------------
# Counting
# ---------------------------------------------------------------------


def count_diamond(geometry: PlaneWaveGeometry, p: np.ndarray,
                  q: np.ndarray, pts: np.ndarray) -> int:
    """Anchor-diamond cardinality over sprinkled points.

    The membership predicate is `_between` -- the SAME function the MC
    volume estimator classifies with, so counts and volumes share one
    definition (P2 design review). The `u` window is strict, so the
    anchors themselves, or any point at exactly their `u`, are never
    counted -- the anchors are external and excluded by construction.
    """

    c = 0
    p0, q0 = float(p[0]), float(q[0])
    for i in range(len(pts)):
        u = pts[i, 0]
        if u <= p0 or u >= q0:
            continue
        if _between(geometry, p, q, pts[i, 0], pts[i, 1],
                    pts[i, 2], pts[i, 3]):
            c += 1
    return c


# ---------------------------------------------------------------------
# The campaign
# ---------------------------------------------------------------------


def bootstrap_var_ci(z: np.ndarray, rng: np.random.Generator,
                     b: int = B_BOOT) -> tuple[float, float]:
    """Percentile 95% bootstrap CI for `Var(Z)` -- the diagnostic's
    empirical side. A chi-square interval would assume normal `Z`,
    which a weighted Poisson difference is not (design review)."""

    n = len(z)
    out = np.empty(b)
    step = max(1, 50_000_000 // max(n, 1))  # cap the index block size
    done = 0
    while done < b:
        k = min(step, b - done)
        idx = rng.integers(0, n, size=(k, n))
        out[done:done + k] = z[idx].var(axis=1, ddof=1)
        done += k
    lo, hi = np.percentile(out, [2.5, 97.5])
    return float(lo), float(hi)


def run_point(label: str, w: float, du: float, dv: float, dx: float,
              dy: float, n_sprinklings: int, point_index: int) -> dict:
    """One operating point's full P2 record: raw sufficient statistics
    plus every derived quantity, so a test can recompute the derived
    fields from the raw ones exactly (the PR #41 lesson applied from
    the start -- nothing hand-transcribed anywhere)."""

    slab = Slab(du=du, dv=dv, dx=dx, dy=dy)
    curved, flat = arms(slab, w)
    assert slab.du < conjugate_du(w)
    p, q = fattest_axis_diamond(curved)
    assert class_c_eligible(curved, p, q), f"{label}: anchor ineligible"

    a = w * slab.du
    r_pred = axis_volume_ratio(a)
    delta = r_pred - 1.0
    dv_win = float(p[1]) - float(q[1])
    v0_exact = math.pi * slab.du ** 2 * dv_win ** 2 / 6.0
    rho = N_BOX / slab.coordinate_volume
    lam0 = rho * v0_exact
    lam_a = r_pred * lam0
    tau = delta          # the chosen effect-scale margin (primary)
    tau_m = delta        # and the marginal's, same yardstick

    # MC cross-check volumes on a MC_SEED sub-stream, disjoint from
    # every campaign stream (design review: MC seed separation)
    mc_rng = np.random.default_rng([MC_SEED, point_index])
    vols = diamond_volumes_mc(curved, flat, p, q, MC_SAMPLES, mc_rng)
    floor = delta * v0_exact
    var_lo = r_pred * rho * max(vols.disagree_lcb, floor)
    var_hi = r_pred * rho * max(vols.disagree_ucb, floor)

    rng = np.random.default_rng(CAMPAIGN_SEEDS[label])
    na = np.empty(n_sprinklings, dtype=np.int64)
    n0 = np.empty(n_sprinklings, dtype=np.int64)
    t0 = time.perf_counter()
    for i in range(n_sprinklings):
        pts = sprinkle(curved, rho, rng)
        na[i] = count_diamond(curved, p, q, pts)
        n0[i] = count_diamond(flat, p, q, pts)
        if (i + 1) % 1000 == 0:
            el = time.perf_counter() - t0
            print(f"    {label}: {i + 1}/{n_sprinklings} "
                  f"({el:.0f}s)", flush=True)
    seconds = time.perf_counter() - t0

    z = na.astype(float) - r_pred * n0.astype(float)
    n = n_sprinklings
    sum_na, sum_n0 = int(na.sum()), int(n0.sum())
    sum_z2 = float((z * z).sum())
    zbar = (sum_na - r_pred * sum_n0) / n
    var_emp = (sum_z2 - n * zbar * zbar) / (n - 1)
    se = math.sqrt(var_emp / n)
    t_crit = student_t_crit(n - 1)
    ci_lo, ci_hi = zbar - t_crit * se, zbar + t_crit * se
    theta_lo, theta_hi = ci_lo / lam0, ci_hi / lam0
    equivalent = bool(-tau < theta_lo and theta_hi < tau)
    discriminates = bool(theta_lo > -delta or theta_hi < -delta)
    # the Bonferroni-widened CI feeds ONLY the joint sentence
    t_bonf = student_t_crit(n - 1, _T_LEVEL_BONF3)
    b_lo = (zbar - t_bonf * se) / lam0
    b_hi = (zbar + t_bonf * se) / lam0
    theta_ci_bonf = [b_lo, b_hi]
    equivalent_bonf = bool(-tau < b_lo and b_hi < tau)
    discriminates_bonf = bool(b_lo > -delta or b_hi < -delta)

    marginal = {}
    for arm_name, total, lam in (("curved", sum_na, lam_a),
                                 ("flat", sum_n0, lam0)):
        m_lo, m_hi = poisson_mean_ci(total)
        r_lo, r_hi = m_lo / (n * lam), m_hi / (n * lam)
        marginal[arm_name] = {
            "total": total,
            "ratio_ci": [r_lo, r_hi],
            "equivalent": bool(1.0 - tau_m < r_lo and r_hi < 1.0 + tau_m),
        }

    boot_rng = np.random.default_rng([MC_SEED, 100 + point_index])
    boot_lo, boot_hi = bootstrap_var_ci(z, boot_rng)

    return {
        "label": label, "w": w, "du": du, "dv": dv, "dx": dx, "dy": dy,
        "n": n, "seed": CAMPAIGN_SEEDS[label],
        "a": a, "r_pred": r_pred, "delta": delta, "rho": rho,
        "lam0": lam0, "lam_a": lam_a, "tau": tau, "tau_m": tau_m,
        "raw": {"sum_na": sum_na, "sum_n0": sum_n0, "sum_z2": sum_z2},
        "zbar": zbar, "var_emp": var_emp, "se": se, "t_crit": t_crit,
        "theta_ci": [theta_lo, theta_hi],
        "equivalent": equivalent, "discriminates": discriminates,
        "theta_ci_bonf": theta_ci_bonf,
        "equivalent_bonf": equivalent_bonf,
        "discriminates_bonf": discriminates_bonf,
        "marginal": marginal,
        "var_diag": {"emp": var_emp, "boot_ci": [boot_lo, boot_hi],
                     "analytic_ci": [var_lo, var_hi],
                     "v_dis_hits": int(vols.disagree_hits),
                     "mc_samples": MC_SAMPLES},
        "seconds": seconds,
    }


# ---------------------------------------------------------------------
# Artifact and rendered tables (the doc embeds these verbatim)
# ---------------------------------------------------------------------

_ARTIFACT = (Path(__file__).resolve().parents[2]
             / "docs" / "prereg" / "p14_probe_p2_results.json")
_P1_SIZING = (Path(__file__).resolve().parents[2]
              / "docs" / "prereg" / "p14_probe_p1_sizing.json")


def check_p1_link() -> None:
    """The sizing inputs behind this protocol are P1's committed
    artifact; the quadrature both scripts share must agree exactly."""

    art = json.loads(_P1_SIZING.read_text(encoding="utf-8"))
    by_label = {r["label"]: r for r in art["points"]}
    for label, w, du, *_ in OPERATING_POINTS:
        rec = by_label[label]
        assert abs(rec["delta"] - (axis_volume_ratio(w * du) - 1.0)) \
            <= 1e-12 * rec["delta"], label


def _sig3(x: float) -> str:
    return f"{x:.3g}"


def primary_table(art: dict) -> str:
    """Primary + secondary, rendered from the artifact."""

    lines = [
        "| point | n | θ 95% CI (t) | τ = δ | equivalent | "
        "discriminates |",
        "|---|---|---|---|---|---|",
    ]
    for r in art["points"]:
        lo, hi = r["theta_ci"]
        lines.append(
            f"| {r['label']} | {r['n']} | [{lo:+.5f}, {hi:+.5f}] "
            f"| {r['tau']:.3e} | {'yes' if r['equivalent'] else 'no'} "
            f"| {'yes' if r['discriminates'] else 'no'} |")
    joint_eq = all(r["equivalent_bonf"] for r in art["points"])
    joint_di = all(r["discriminates_bonf"] for r in art["points"])
    lines.append(
        "\nJoint sentence over the three points (Bonferroni, per-point "
        "t at 99.17%): equivalent "
        f"{'yes' if joint_eq else 'no'}, discriminates "
        f"{'yes' if joint_di else 'no'}. The primary record stays "
        "pointwise.")
    return "\n".join(lines)


def marginal_table(art: dict) -> str:
    """Marginal rate-ratio equivalence per arm, rendered."""

    lines = [
        "| point | arm | rate ratio 95% CI (Garwood) | 1 ± τ_m | "
        "equivalent |",
        "|---|---|---|---|---|",
    ]
    for r in art["points"]:
        for arm_name in ("curved", "flat"):
            m = r["marginal"][arm_name]
            lo, hi = m["ratio_ci"]
            lines.append(
                f"| {r['label']} | {arm_name} | [{lo:.5f}, {hi:.5f}] "
                f"| 1 ± {r['tau_m']:.3e} "
                f"| {'yes' if m['equivalent'] else 'no'} |")
    return "\n".join(lines)


def diagnostic_table(art: dict) -> str:
    """Var(Z) empirical vs analytic, rendered. Reported, not judged."""

    lines = [
        "| point | Var(Z) emp [boot 95%] | analytic r·ρ·V_dis "
        "[CP 95%] | ratio |",
        "|---|---|---|---|",
    ]
    for r in art["points"]:
        d = r["var_diag"]
        b_lo, b_hi = d["boot_ci"]
        a_lo, a_hi = d["analytic_ci"]
        mid = 0.5 * (a_lo + a_hi)
        lines.append(
            f"| {r['label']} | {d['emp']:.3f} [{b_lo:.3f}, {b_hi:.3f}] "
            f"| [{a_lo:.3f}, {a_hi:.3f}] | {d['emp'] / mid:.2f} |")
    return "\n".join(lines)


def write_results_json(art: dict) -> Path:
    with open(_ARTIFACT, "w", encoding="utf-8", newline="\n") as f:
        json.dump(art, f, indent=2, sort_keys=True)
        f.write("\n")
    return _ARTIFACT


def main() -> None:
    assert_seed_layout()
    check_p1_link()
    print(f"P14 probe P2 -- exploratory, E[N_box] = {N_BOX}, "
          f"mc = {MC_SAMPLES}, seeds {sorted(CAMPAIGN_SEEDS.values())}"
          f" + MC {MC_SEED}\n")

    art = {
        "script": "experiments/positive_control/p14_probe_p2.py",
        "protocol": {
            "n_box": N_BOX, "mc_samples": MC_SAMPLES,
            "mc_seed": MC_SEED, "campaign_seeds": CAMPAIGN_SEEDS,
            "burned_seeds": list(BURNED_SEEDS),
            "ci": "student-t, 95% two-sided",
            "tau_rule": "delta", "tau_m_rule": "delta",
            "bootstrap_b": B_BOOT,
        },
        "points": [],
    }
    for k, (label, w, du, dv, dx, dy, n) in enumerate(OPERATING_POINTS):
        print(f"  {label}: n = {n}")
        art["points"].append(run_point(label, w, du, dv, dx, dy, n, k))
        print(f"  {label} done ({art['points'][-1]['seconds']:.0f}s)")

    print("\n== primary: theta equivalence (tau = delta) ==")
    print(primary_table(art))
    print("\n== marginal rate-ratio equivalence (tau_m = delta) ==")
    print(marginal_table(art))
    print("\n== Var(Z) diagnostic (reported, not judged) ==")
    print(diagnostic_table(art))
    print(f"\nresults artifact written: {write_results_json(art)}")


if __name__ == "__main__":
    main()
