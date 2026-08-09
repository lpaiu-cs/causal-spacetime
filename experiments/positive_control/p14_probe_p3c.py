"""P14 §8 P3-C: confirmation/termination stage -- PREFLIGHT + runner.

This module is the OPERATIVE certifier and campaign runner for P3-C;
it supersedes the normal-model preliminary record in `p14_probe_p3e`
(kept there as the historical design-stage candidate). The campaign
itself is NOT run by the preflight: per the review, execution waits
until this machinery -- the three-branch rules, both-branch empirical
certification with the campaign's own pipeline, and the raw-integer
artifact -- is committed and green.

**Three-branch verdict (frozen).** Each of the three metrics gets a
95% CI by the FROZEN construction below; the composite verdict is

- ``confirmed``   -- every CI lies ENTIRELY outside its equivalence
  band in the frozen direction (curved above flat at aniso-a1.0):
  `s_lo > +eps_s`, `auc_lo > 0.5 + eps_auc`, `ba_lo > 0.5 + eps_ba`;
- ``equivalent``  -- every CI lies entirely INSIDE its band;
- ``inconclusive`` -- everything else. Equivalence failure alone is
  never separation.

The positive sentence, frozen: "aniso-a1.0의 고정된 유한 상자·밀도에서
global relation fraction이 평탄 ensemble과 분리된다" -- never a
general-Weyl or box-independent claim. The termination sentence is
the three-rule closure frozen in the P3 design review.

**Frozen CI constructions** (identical in certification and
campaign; this identity is the review's requirement):

- `s`:   `s_hat ± z * sqrt(2/n + s_hat^2 / (4(n-1)))` -- the full
  variance; the null-form `sqrt(2/n)` is ANTI-conservative at large
  observed `s`;
- `AUC`: `auc ± z * sqrt(S10/m + S01/n)` -- the DELONG variance from
  the placement values (ties via midplacements). The null-form
  `(2n+1)/(12 n^2)` used by the first preflight is the variance
  under the identical-continuous-distributions null ONLY, not an
  upper bound under alternatives (PR #47 review: flat constant with
  curved two-valued at 0.55/0.45 gives AUC = 0.55 with asymptotic
  variance `0.2475/n`, ~48% above `1/(6n)` -- an anti-conservative
  CI that can manufacture `confirmed`); DeLong tracks the
  alternative and the certification was re-run with it;
- `BA`:  `ba ± z / (2 sqrt(n))` -- the binomial variance maximum
  `p(1-p) <= 1/4` (test half is `n/2` per arm), genuinely
  conservative. The threshold and direction are learned on the first
  half of each arm's stream and applied to the second half, in
  certification replicates and in the campaign alike.

**Certification (both branches, empirical distributions).** Powers
are measured on bootstrap draws from P3-E's stored aniso-a1.0
per-sprinkling 3a samples, the two arms resampled INDEPENDENTLY
(index sharing would certify a paired design, and P3-C is unpaired):
the null draws both arms from the flat sample, the effect branch
draws each arm from its own sample. Certification rule: the exact CP
95% LOWER bound of the branch's pass rate must reach 0.90 --
equivalence branch at the null, confirmed branch at the effect.
`n = 4800/arm` (the equivalence branch binds; the effect branch is
saturated). Margins are `p14_probe_p3e.P3C_MARGINS`, frozen -- they
never move to meet a power result.

**Campaign (held until preflight approval).** Two UNPAIRED arms at
aniso-a1.0, `E[N] = 300`, fresh seed streams (curved 20260831, flat
20260832; one rng stream per arm), `n = 4800` sprinklings per arm,
3a per sprinkling, then the three CIs and the branch verdict, all
into a raw-integer artifact.

Run:  python experiments/positive_control/p14_probe_p3c.py preflight
      python experiments/positive_control/p14_probe_p3c.py campaign
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from p14_plane_wave import Slab, arms, sprinkle
from p14_probe_p1 import clopper_pearson, relation_census
from p14_probe_p3e import P3C_MARGINS

# ---------------------------------------------------------------------
# Frozen constants
# ---------------------------------------------------------------------

N_ARM = 4_800
POINT = ("aniso-a1.0", 1.0, 1.0, 1.0, 2.0, 6.0)
E_N = 300

#: One rng stream per arm; unpaired by construction. 20260831/32 fed
#: the FIRST campaign, whose verdict was downgraded to exploratory in
#: the PR #48 review (the AUC interval degenerated at complete
#: separation, so the frozen valid-CI requirement was unmet at run
#: time) -- they are burned, and the fresh confirmation block runs on
#: new streams under the corrected boundary rule.
CAMPAIGN_SEEDS = {"curved": 20260841, "flat": 20260842}
BURNED_SEEDS = (20260808, 777, 778, 779, 780, 781,
                20260811, 20260812, 20260813, 20260814,
                20260821, 20260822, 20260823, 20260824,
                20260831, 20260832)

PREFLIGHT_NULL_REPS = 20_000
PREFLIGHT_EFFECT_REPS = 4_000
#: Substreams of the burned design lineage (781) -- documented, never
#: campaign streams.
PREFLIGHT_NULL_SEED = (781, 50)
PREFLIGHT_EFFECT_SEED = (781, 51)

_Z = 1.959964

_P3E_ARTIFACT = (Path(__file__).resolve().parents[2]
                 / "docs" / "prereg" / "p14_probe_p3e_results.json")
_PREFLIGHT_ARTIFACT = (Path(__file__).resolve().parents[2]
                       / "docs" / "prereg"
                       / "p14_probe_p3c_preflight.json")
_CAMPAIGN_ARTIFACT = (Path(__file__).resolve().parents[2]
                      / "docs" / "prereg" / "p14_probe_p3c_results.json")


def assert_seed_layout() -> None:
    seeds = list(CAMPAIGN_SEEDS.values())
    assert len(set(seeds)) == len(seeds)
    clash = set(seeds) & set(BURNED_SEEDS)
    assert not clash, f"burned seeds reused: {sorted(clash)}"


# ---------------------------------------------------------------------
# Frozen CI constructions and the three-branch verdict
# ---------------------------------------------------------------------


def ci_s(s_hat: float, n: int) -> tuple[float, float]:
    se = math.sqrt(2.0 / n + s_hat * s_hat / (4.0 * (n - 1)))
    return s_hat - _Z * se, s_hat + _Z * se


def auc_delong(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    """(AUC, ci_lo, ci_hi) with the DeLong variance.

    Placements: `V10_i = P(b < a_i) + 0.5 P(b = a_i)` over the flat
    sample, `V01_j = P(a > b_j) + 0.5 P(a = b_j)` over the curved
    one; `AUC = mean(V10) = mean(V01)` and
    `Var = var(V10)/m + var(V01)/n`. Valid under alternatives and
    ties -- the null-form `1/(6n)` is not an upper bound off the
    null (PR #47 review) and appears nowhere in the verdict path.
    At the boundaries the frozen rule switches to the exact
    placement bound (PR #48 review; see inline).
    """

    m, n = len(a), len(b)
    b_sorted = np.sort(b)
    left = np.searchsorted(b_sorted, a, side="left")
    right = np.searchsorted(b_sorted, a, side="right")
    v10 = (left + 0.5 * (right - left)) / n
    a_sorted = np.sort(a)
    hi_idx = np.searchsorted(a_sorted, b, side="right")
    lo_idx = np.searchsorted(a_sorted, b, side="left")
    v01 = ((m - hi_idx) + 0.5 * (hi_idx - lo_idx)) / m
    auc = float(v10.mean())
    # At the boundaries the Wald/DeLong interval DEGENERATES to zero
    # width -- complete separation makes every placement 1 (or 0), the
    # influence variance vanishes, and [1, 1] is not a valid
    # finite-sample CI for the POPULATION AUC (PR #48 review: the
    # empirical AUC being 1 shows the sample separates, not that the
    # population overlap is exactly zero). The frozen boundary rule is
    # the exact placement bound: all m values of V10 equal 1 means the
    # one-sided 97.5% CP lower bound on P(V10 = 1) is 0.025^(1/m), and
    # AUC = E[V10] >= P(V10 = 1), so [0.025^(1/m), 1] is a valid 95%
    # interval; symmetric at 0.
    if auc >= 1.0:
        return 1.0, 0.025 ** (1.0 / m), 1.0
    if auc <= 0.0:
        return 0.0, 0.0, 1.0 - 0.025 ** (1.0 / m)
    var = (v10.var(ddof=1) / m if m > 1 else 0.0) \
        + (v01.var(ddof=1) / n if n > 1 else 0.0)
    se = math.sqrt(max(var, 0.0))
    return auc, auc - _Z * se, auc + _Z * se


def ci_ba(ba: float, n: int) -> tuple[float, float]:
    se = 1.0 / (2.0 * math.sqrt(n))
    return ba - _Z * se, ba + _Z * se


def classify(s_ci, auc_ci, ba_ci, margins=P3C_MARGINS) -> str:
    """The frozen three-branch verdict; the positive direction is
    curved-above-flat."""

    es, ea, eb = margins["eps_s"], margins["eps_auc"], margins["eps_ba"]
    confirmed = (s_ci[0] > es and auc_ci[0] > 0.5 + ea
                 and ba_ci[0] > 0.5 + eb)
    equivalent = (-es < s_ci[0] and s_ci[1] < es
                  and 0.5 - ea < auc_ci[0] and auc_ci[1] < 0.5 + ea
                  and 0.5 - eb < ba_ci[0] and ba_ci[1] < 0.5 + eb)
    if confirmed:
        return "confirmed"
    if equivalent:
        return "equivalent"
    return "inconclusive"


# ---------------------------------------------------------------------
# The metric pipeline (one dataset -> the three CIs), shared verbatim
# by certification replicates and the campaign
# ---------------------------------------------------------------------


def metric_cis(a: np.ndarray, b: np.ndarray) -> dict:
    """(curved sample, flat sample) -> the three frozen CIs plus the
    point statistics. Threshold and direction are learned on the
    first half of each arm's STREAM ORDER and applied to the second
    half; ties classify as flat."""

    n = len(a)
    assert len(b) == n
    sd = math.sqrt(0.5 * (a.var(ddof=1) + b.var(ddof=1)))
    s_hat = float((a.mean() - b.mean()) / sd) if sd > 0 else 0.0
    auc, auc_lo, auc_hi = auc_delong(a, b)
    h = n // 2
    ta, tb = a[:h].mean(), b[:h].mean()
    thr = 0.5 * (ta + tb)
    if ta > tb:
        acc_a = float((a[h:] > thr).mean())
        acc_b = float((b[h:] <= thr).mean())
    elif ta < tb:
        acc_a = float((a[h:] < thr).mean())
        acc_b = float((b[h:] >= thr).mean())
    else:
        acc_a, acc_b = 0.0, 1.0
    ba = 0.5 * (acc_a + acc_b)
    return {"s": s_hat, "auc": auc, "ba": ba,
            "ci_s": ci_s(s_hat, n), "ci_auc": (auc_lo, auc_hi),
            "ci_ba": ci_ba(ba, n)}


# ---------------------------------------------------------------------
# Preflight certification (both branches, empirical distributions)
# ---------------------------------------------------------------------


def _branch_counts(f_a: np.ndarray, f_b: np.ndarray, n_arm: int,
                   reps: int, seed_path) -> dict:
    """Replicates of the FULL pipeline on independent bootstrap draws
    (arm A from `f_a`, arm B from `f_b`, separate index draws) and
    the three-branch verdict for each."""

    rng = np.random.default_rng(list(seed_path))
    counts = {"confirmed": 0, "equivalent": 0, "inconclusive": 0}
    for _ in range(reps):
        a = f_a[rng.integers(0, len(f_a), n_arm)]
        b = f_b[rng.integers(0, len(f_b), n_arm)]
        m = metric_cis(a, b)
        counts[classify(m["ci_s"], m["ci_auc"], m["ci_ba"])] += 1
    return counts


def preflight(p3e_art: dict) -> dict:
    rec = next(r for r in p3e_art["ladder_g"]
               if r["label"] == POINT[0])
    fa = np.asarray(rec["raw"]["f_curved"])
    f0 = np.asarray(rec["raw"]["f_flat"])

    t0 = time.perf_counter()
    null = _branch_counts(f0, f0, N_ARM, PREFLIGHT_NULL_REPS,
                          PREFLIGHT_NULL_SEED)
    eff = _branch_counts(fa, f0, N_ARM, PREFLIGHT_EFFECT_REPS,
                         PREFLIGHT_EFFECT_SEED)
    null_lo, null_hi = clopper_pearson(null["equivalent"],
                                       PREFLIGHT_NULL_REPS)
    eff_lo, eff_hi = clopper_pearson(eff["confirmed"],
                                     PREFLIGHT_EFFECT_REPS)
    out = {
        "n_arm": N_ARM,
        "margins": dict(P3C_MARGINS),
        "ci_constructions": "s: z*sqrt(2/n + s^2/(4(n-1))); "
                            "auc: DeLong placement variance "
                            "z*sqrt(S10/m + S01/n), midplacement "
                            "ties; ba: z/(2 sqrt(n)); z = 1.959964",
        "resampling": "arms drawn INDEPENDENTLY (unpaired); null = "
                      "both arms from the flat sample",
        "source_artifact": "p14_probe_p3e_results.json ladder_g "
                           "aniso-a1.0 raw samples",
        "null": {"reps": PREFLIGHT_NULL_REPS,
                 "seed_path": list(PREFLIGHT_NULL_SEED),
                 "counts": null,
                 "equivalent_rate": null["equivalent"]
                 / PREFLIGHT_NULL_REPS,
                 "equivalent_ci95_exact": [null_lo, null_hi]},
        "effect": {"reps": PREFLIGHT_EFFECT_REPS,
                   "seed_path": list(PREFLIGHT_EFFECT_SEED),
                   "counts": eff,
                   "confirmed_rate": eff["confirmed"]
                   / PREFLIGHT_EFFECT_REPS,
                   "confirmed_ci95_exact": [eff_lo, eff_hi]},
        "certification_rule": "exact CP 95% LOWER bound >= 0.90 for "
                              "BOTH branches",
        "certified": bool(null_lo >= 0.90 and eff_lo >= 0.90),
        "campaign_seeds": dict(CAMPAIGN_SEEDS),
        "seconds": time.perf_counter() - t0,
    }
    assert out["certified"], "preflight certification failed"
    return out


# ---------------------------------------------------------------------
# P4 from the stored samples (pair-level joint bootstrap)
# ---------------------------------------------------------------------

P4_BOOT_B = 4_000
P4_BOOT_SEED = (781, 41)


def p4_block(p3e_art: dict) -> dict:
    """Paired vs independent-arm variance of the aniso 3a difference,
    with a PAIR-LEVEL joint bootstrap (sprinkling pairs resampled
    together -- the pairing is the object being measured)."""

    rec = next(r for r in p3e_art["ladder_g"]
               if r["label"] == POINT[0])
    fa = np.asarray(rec["raw"]["f_curved"])
    f0 = np.asarray(rec["raw"]["f_flat"])
    n = len(fa)
    paired = float((fa - f0).var(ddof=1))
    indep = float(fa.var(ddof=1) + f0.var(ddof=1))
    rng = np.random.default_rng(list(P4_BOOT_SEED))
    g_b = np.empty(P4_BOOT_B)
    c_b = np.empty(P4_BOOT_B)
    for b in range(P4_BOOT_B):
        i = rng.integers(0, n, n)          # ONE index vector: pairs
        a, z = fa[i], f0[i]
        g_b[b] = (a.var(ddof=1) + z.var(ddof=1)) / (a - z).var(ddof=1)
        c_b[b] = np.corrcoef(a, z)[0, 1]
    return {
        "point": POINT[0], "n_sprinklings": n,
        "paired_variance": paired,
        "independent_arm_variance_sum": indep,
        "variance_gain": indep / paired,
        "arm_correlation": float(np.corrcoef(fa, f0)[0, 1]),
        "bootstrap": {"B": P4_BOOT_B, "seed_path": list(P4_BOOT_SEED),
                      "resampling": "pair-level joint (one index "
                                    "vector for both arms)",
                      "gain_ci95": [float(np.percentile(g_b, 2.5)),
                                    float(np.percentile(g_b, 97.5))],
                      "corr_ci95": [float(np.percentile(c_b, 2.5)),
                                    float(np.percentile(c_b, 97.5))]},
        "sizing_note": "paired design needs 1/gain of the unpaired "
                       "sprinklings for equal mean-difference "
                       "precision",
        "source_artifact": "p14_probe_p3e_results.json ladder_g "
                           "aniso-a1.0 raw samples",
    }


# ---------------------------------------------------------------------
# The campaign (HELD until preflight approval; code committed only)
# ---------------------------------------------------------------------


def arm_samples(geom, rho: float, seed: int,
                count: int) -> tuple[np.ndarray, int, int]:
    """One unpaired arm's per-sprinkling 3a values from a single rng
    stream, plus the arm's TOTAL ambiguous and escalated pair counts.

    `related` leaves an undecided pair False, so summing it without
    reading `ambiguous` silently counts undecided as unrelated --
    the "undecided is never a silent False" contract (PR #48 review
    R2). The campaign stores both totals and asserts ambiguity is
    zero; a nonzero count would require bracketing the relation
    fraction instead. Module-level so a test can rerun a PREFIX of
    the stream against the committed raw record.
    """

    rng = np.random.default_rng(seed)
    out = np.empty(count)
    amb = esc = 0
    for i in range(count):
        pts = sprinkle(geom, rho, rng)
        n = len(pts)
        cen = relation_census(geom, pts)
        out[i] = cen.related.sum() / (n * (n - 1) / 2)
        amb += cen.ambiguous
        esc += cen.escalated
        if (i + 1) % 500 == 0:
            print(f"    {seed}: {i + 1}/{count}", flush=True)
    return out, amb, esc


def run_campaign() -> dict:
    """Two unpaired arms, fresh seed streams, n = 4800 sprinklings
    per arm; 3a per sprinkling; the frozen pipeline and verdict."""

    assert_seed_layout()
    label, w, du, dv, dx, dy = POINT
    slab = Slab(du=du, dv=dv, dx=dx, dy=dy)
    rho = E_N / slab.coordinate_volume
    curved, flat = arms(slab, w)

    print("  curved arm ...", flush=True)
    f_curved, amb_a, esc_a = arm_samples(curved, rho,
                                         CAMPAIGN_SEEDS["curved"], N_ARM)
    print("  flat arm ...", flush=True)
    f_flat, amb_0, esc_0 = arm_samples(flat, rho,
                                       CAMPAIGN_SEEDS["flat"], N_ARM)
    # undecided is never a silent False: with any ambiguity the
    # relation fraction would need bracketing, not a point value
    assert amb_a == 0 and amb_0 == 0, (amb_a, amb_0)
    m = metric_cis(f_curved, f_flat)
    verdict = classify(m["ci_s"], m["ci_auc"], m["ci_ba"])
    art = {
        "script": "experiments/positive_control/p14_probe_p3c.py",
        "point": label, "e_n": E_N, "n_arm": N_ARM,
        "seeds": dict(CAMPAIGN_SEEDS),
        "margins": dict(P3C_MARGINS),
        "ambiguity": {"curved": {"ambiguous": amb_a, "escalated": esc_a},
                      "flat": {"ambiguous": amb_0, "escalated": esc_0}},
        "raw": {"f_curved": [float(v) for v in f_curved],
                "f_flat": [float(v) for v in f_flat]},
        "metrics": m,
        "verdict": verdict,
        "positive_sentence_scope": "aniso-a1.0의 고정된 유한 상자·밀도에서 "
                                   "global relation fraction이 평탄 "
                                   "ensemble과 분리된다",
    }
    with open(_CAMPAIGN_ARTIFACT, "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(art, f, indent=2, sort_keys=True)
        f.write("\n")
    print(f"verdict: {verdict}; artifact: {_CAMPAIGN_ARTIFACT}")
    return art


def campaign_table(art: dict) -> str:
    """The campaign verdict table the results doc embeds verbatim."""

    m = art["metrics"]
    mg = art["margins"]
    rows = [
        ("s", m["s"], m["ci_s"], f"±{mg['eps_s']:.4f}"),
        ("AUC", m["auc"], m["ci_auc"], f"0.5 ± {mg['eps_auc']:.4f}"),
        ("BA", m["ba"], m["ci_ba"], f"0.5 ± {mg['eps_ba']:.4f}"),
    ]
    lines = ["| metric | value | 95% CI | band |",
             "|---|---|---|---|"]
    for name, v, ci, band in rows:
        lines.append(f"| {name} | {v:.6f} | [{ci[0]:.6f}, "
                     f"{ci[1]:.6f}] | {band} |")
    lines.append(f"\n**Verdict: `{art['verdict']}`** — every CI "
                 "entirely outside its band in the frozen direction.")
    return "\n".join(lines)


def main(mode: str) -> None:
    assert_seed_layout()
    p3e_art = json.loads(_P3E_ARTIFACT.read_text(encoding="utf-8"))
    if mode == "preflight":
        art = {
            "script": "experiments/positive_control/p14_probe_p3c.py",
            "preflight": preflight(p3e_art),
            "p4": p4_block(p3e_art),
        }
        with open(_PREFLIGHT_ARTIFACT, "w", encoding="utf-8",
                  newline="\n") as f:
            json.dump(art, f, indent=2, sort_keys=True)
            f.write("\n")
        p = art["preflight"]
        print(f"null equivalent {p['null']['counts']['equivalent']}"
              f"/{p['null']['reps']} CP "
              f"{p['null']['equivalent_ci95_exact']}")
        print(f"effect confirmed {p['effect']['counts']['confirmed']}"
              f"/{p['effect']['reps']} CP "
              f"{p['effect']['confirmed_ci95_exact']}")
        print(f"certified: {p['certified']}")
        print(f"p4 gain {art['p4']['variance_gain']:.3f} "
              f"{art['p4']['bootstrap']['gain_ci95']}")
        print(f"artifact: {_PREFLIGHT_ARTIFACT} ({p['seconds']:.0f}s)")
    elif mode == "campaign":
        run_campaign()
    else:
        raise SystemExit("mode: preflight | campaign")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "preflight")
