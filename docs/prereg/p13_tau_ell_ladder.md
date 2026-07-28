# P13: how far does curvature let the flat-normalized reading go? — power-first preregistration

Status: **DESIGN v1.0 (2026-07-28), for in-session review. Nothing
below has been run.** Dates are local (UTC+9); commit timestamps carry
their +09:00 offset.

Lineage and the question. P12 showed that the unchanged P11 chain
estimator reads CURVED proper time in `dS_2` at `tau / ell ~ 0.3`,
improving with density at the longest-chain rate (slope `-0.304` over
four rungs in the archived diagnostic, `Delta = -0.193` against a
derived `-0.2007`). P12's own Section 7 named the boundary it could
not see: the chain law is asymptotic in `interval size / curvature
radius`, so the open question is **how large `tau / ell` may grow
before a flat-normalized reading degrades, and whether the `-1/3`
exponent survives there.**

**The design crux: vary `tau / ell` at FIXED discreteness.** P12
swept density at fixed `tau / ell`; P13 does the opposite, and the
comparison is only clean if the interval count `m` is held constant
across rungs — otherwise a curvature effect and a discreteness effect
are the same number. Every rung therefore gets its own patch and its
own density, chosen so the realized `m` matches. That matching is
itself gated (Section 5), because a design that silently drifts in `m`
would answer a different question than the one asked.

---

## 1. Power design (first, by house rule)

### 1.1 Ladder, statistic, contrast

Ladder: `tau / ell in {0.30, 0.60, 1.00, 1.50}`, each a band of
`+/- 10%`, at a common target `m = 76` (P12's top-rung count).
`ell = 1` throughout, so `R = 2` is fixed and only the interval's
size relative to the curvature radius moves.

Per-sample statistic, unchanged in form from P11 and P12:
`y = log10( median over K = 6 pairs of
|tau_hat - tau_curved| / tau_curved )`. Rung statistic: the MEAN of
`y`. Primary contrast

    Delta_13 = mean_y(tau/ell = 1.50) - mean_y(tau/ell = 0.30).

### 1.2 The polarity is inverted, and the verdict table says so

In P11 and P12 the interesting outcome was a NEGATIVE contrast. Here
the null is `Delta_13 = 0` — the flat-normalized reading is
curvature-robust at fixed discreteness — and the interesting outcome
is POSITIVE: the reading degrades as the interval approaches the
curvature radius. The frozen table, evaluated by precedence so an
observed CI has exactly one outcome (`lo`, `hi` = the 95% bootstrap
interval, `delta_eq = 0.05` dex):

| # | condition | verdict |
|---|---|---|
| 1 | `lo > 0` | **CURVATURE-DEGRADES** — the boundary is located; a positive finding |
| 2 | `hi < 0` | **CURVATURE-HELPS** — anomaly, representable, would need explaining |
| 3 | `-delta_eq < lo` and `hi < delta_eq` | **CURVATURE-ROBUST** — the reading survives to `tau/ell = 1.5`; available only if the pilot declared it affordable (1.3) |
| 4 | otherwise | **INCONCLUSIVE** if row 3 was affordable, else **UNRESOLVED** |

`delta_eq = 0.05` dex is 12% in relative error: the threshold below
which "the reading did not notice the curvature" is a fair
description. It is a design choice, stated as one, not a derived
quantity.

### 1.3 Sample size, and affording both verdicts

From the pilot's calibrated Bonett bounds (P11 Section 1.2 machinery,
applied verbatim, including publishing the calibrated bound that
power consumes):

    n_sup = ceil( S^2_90 * (1.960 + 1.282)^2 / 0.15^2 )   # 90% power
                                                          # at Delta = 0.15
    n_eq  = ceil( S^2_90 * (1.960 + 1.645)^2 / 0.05^2 )   # 90% P(ROBUST)
                                                          # at Delta = 0

`n_per_rung = clamp(max(n_sup, n_eq), 12, 60)` when `n_eq` fits the
cap, else `clamp(n_sup)` with CURVATURE-ROBUST declared unavailable
ex ante. The detectable effect `0.15` dex is a design choice too: it
is roughly the size a `(tau/ell)^2` correction with an order-unity
coefficient would produce at the top rung (Section 7), stated as an
order-of-magnitude expectation and NOT frozen as a target — the
coefficient is precisely the unverified-constant class that killed
P11's Stage C v1.

Stage P-13 pilots the two ENDPOINT rungs (0.30 and 1.50), 200 samples
each, cross-rung statistics forbidden.

---

## 2. Ensemble and the frozen per-rung constants

Same `dS_2` flat slicing, same causal rule, same genuine
inhomogeneous Poisson sprinkling by thinning as P12 Section 2 — only
the patch and density change per rung. Each patch is
`eta in [eta_lo, -1]`, `|x| <= X`, and `rho` is set so the realized
interval count matches the target.

| `tau/ell` | `eta_lo` | `X` | `rho` | design-check completion | realized mean `m` | typical `N` |
|---|---|---|---|---|---|---|
| 0.30 | -1.8 | 1.8 | 1707.0 | 100% | 76.6 | 2733 |
| 0.60 | -2.6 | 4.0 | 434.8 | 100% | 76.7 | 2137 |
| 1.00 | -4.0 | 8.0 | 163.0 | 100% | 77.0 | 1959 |
| 1.50 | -7.0 | 17.0 | 76.1 | 100% | 75.8 | 2218 |

Two things in that table are the P12 lessons applied before freezing
rather than after. **The `X` values are large and grow steeply**
because a pair at `tau/ell = 1.5` spans a coordinate `eta` interval
comparable to the whole patch, so six disjoint boxes can only be
stacked along `x`; the first candidate constants (`X` of 1.5 to 4.5)
completed 100% / 10% / 0% / 0%, the same packing collapse P12 v1.0
suffered, and were rejected in design check. **The `rho` values are
calibrated, not analytic**: at fixed `rho tau^2` the realized `m`
drifted downward by ~10% across the ladder — that drift IS the
curvature correction to the diamond volume — so each `rho` was
adjusted until the realized `m` matched, which is why the fourth
column is flat to 1.6%.

Design-check seed space: `7000000-7999999` (the runs above used
`7000000` and `7100000`); experimental windows are `6000000-6999999`
(Section 6). No experimental data has been produced.

## 3. Truth and estimator

Truth is P12's exact `dS_2` geodesic proper time via the de Sitter
invariant, already pinned against an independent geodesic integration
at `1e-9` (P12 tests). The estimator is P11's, CALLED through
`estimate_tau_from_longest_chain_1p1` with the rung's `rho` — never
re-derived (the P12 review finding), so a recalibration of the shared
definition reaches all three experiments at once.

## 4. Protocol

P12's, verbatim: both frozen eligibility conditions (band membership
AND causal box inside the patch, both decided on coordinates before
any measurement), greedy disjoint-support draws with the 200-rejection
cap, first-`n`-complete fill with reserve slots and published skip
identities, the `>= 1998 of 2000` verification pin per rung with
measured wall times, clean-worktree preflight, stamp-equality within a
chain, and one-significant-figure quoting of every measured timing.
Cross-stage gate: P13 runs only after P12 Stage A's frozen IMPROVES
with a reachable stamp, since P13's instrument is that same estimator.

## 5. Controls that can void the contrast (frozen)

The design's whole claim rests on discreteness being held fixed, so
two controls gate the verdict rather than decorating it:

1. **`m`-matching.** Realized mean `m` per rung must lie within
   `+/- 5%` of the campaign's grand mean. Outside that, the record
   reports **CONFOUNDED** instead of a verdict: the contrast would be
   mixing curvature with discreteness. (Design check: 75.8 to 77.0,
   a spread of 1.6%.)
2. **Flat twin.** The same protocol run in FLAT Minkowski at matched
   `m` and matched coordinate box sizes, where the true `Delta` is
   zero by construction. If the flat twin's contrast exceeds
   `delta_eq` in magnitude, the protocol itself has a size-dependent
   artifact and the record reports **CONFOUNDED**. This is the
   control P12's Stage B addendum asked for, promoted to a gate.

Labelled, never gating: the estimator's error against the flat
template, and the geometric gap. P12's lesson stands — the first of
these is density- and geometry-dependent, so no level is frozen for
it this time; it is reported, not scored.

## 6. Seed windows (frozen; above every range used)

| block | base | fills | window span |
|---|---|---|---|
| verification-13 (non-experimental) | 6000000 | 2000 per rung | 6000000-6007999 |
| Stage P-13 pilot, tau/ell = 0.30 | 6020000 | 200 of 220 | 6020000-6063999 |
| Stage P-13 pilot, tau/ell = 1.50 | 6064000 | 200 of 220 | 6064000-6107999 |
| Stage A-13, tau/ell = 0.30 | 6120000 | <= 60 of 80 | 6120000-6135999 |
| Stage A-13, tau/ell = 0.60 | 6136000 | <= 60 of 80 | 6136000-6151999 |
| Stage A-13, tau/ell = 1.00 | 6152000 | <= 60 of 80 | 6152000-6167999 |
| Stage A-13, tau/ell = 1.50 | 6168000 | <= 60 of 80 | 6168000-6183999 |
| flat-twin control | 6200000 | <= 60 of 80 per rung | 6200000-6263999 |

All spans are disjoint from P11's (to 705999), P12's experimental
(1000000-1999999) and design (2000000-5999999) spaces, and from
P13's own design space (7000000+).

## 7. Theory appendix: what a curvature effect should look like

The chain law reads the longest chain as the geodesic proper time in
the regime where the interval is small compared with the curvature
radius; the correction is controlled by `(tau / ell)^2`, the same
combination that governs the small-diamond volume expansion. Across
this ladder `(tau/ell)^2` grows 25-fold, from 0.09 to 2.25, while
discreteness is pinned — so if the correction carries an order-unity
coefficient it should be plainly visible at the top rung, and if the
coefficient is small the ladder should read CURVATURE-ROBUST and the
instrument's usable range extends further than the theory's comfort
zone. Both are publishable readings of the same experiment, which is
the property a good gate has.

What this design cannot see: any effect that is not monotone in
`tau / ell`, and anything beyond `tau / ell = 1.5`, where the patch
would have to grow faster than the `O(N^2)` pool scan can afford. The
`m = 76` operating point is also a choice; a floor that appears only
at larger `m` would be invisible here, and Section 9 will say so
rather than generalize.

## 8. Implementation plan (after this design is approved)

Extend `experiments/positive_control/p12_curved.py` into a
per-rung-parameterized module (patch and `rho` become rung
parameters rather than module constants) and add `p13_tau_ell.py`
with the stage runners; the estimator, truth, eligibility, packing,
and gate machinery are imported, not re-implemented. Tests before any
run: the per-rung packing pin (>= 99% completion on the design-check
seeds, all four rungs); the `m`-matching pin; the flat-twin's
zero-contrast pin; verdict-table precedence including the new
polarity and the ROBUST-unavailable branch; and window privacy and
freshness against every documented range. Stage P-13 runs only from a
clean commit containing all of it, after the verification pin passes.
