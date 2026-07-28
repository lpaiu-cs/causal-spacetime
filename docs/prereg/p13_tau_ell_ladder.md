# P13: how far does curvature let the flat-normalized reading go? — power-first preregistration

Status: **DESIGN v1.1 (2026-07-28), for in-session review. Nothing
below has been run.** Dates are local (UTC+9); commit timestamps
carry their +09:00 offset.

v1.1 corrections, pre-freeze and pre-data, from design review (the
central one reproduced independently before it was applied): **the
sample cap rises to 200 for P13, because at cap 60 this design could
not afford its own most likely outcome.** Two independent
design-check runs put the endpoint contrast at zero
(`+0.006 +/- 0.025` and `+0.0074 +/- 0.0296`) with per-rung
`s(y) ~ 0.11`, i.e. `S^2 ~ 0.025-0.027`, which makes `n_eq ~ 143`
before calibration and ~170 after — far above 60. Under v1.0 the
pilot would therefore have declared CURVATURE-ROBUST unavailable ex
ante and a 12-sample campaign would have produced intervals too wide
for rows 1 and 2, so the near-certain reading would have been
UNRESOLVED by construction. Cap 60 is a P11 inheritance from an era
when a sample cost minutes; a P13 sample costs about 0.2 s, so 200
per rung across four rungs plus the flat twin is a few minutes. Seed
windows are re-cut for 200-of-220 blocks (Section 6), and the
design-check contrast is disclosed here rather than discovered later
(the P12 precedent for quarantined design numbers).

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

Rows 1 and 2 can fire on a significant contrast SMALLER than
`delta_eq` — significance and size are different questions, and the
table answers only the first. Section 9's narration must therefore
state the magnitude against `delta_eq` alongside the verdict word,
never the word alone.

### 1.3 Sample size, and affording both verdicts

From the pilot's calibrated Bonett bounds (P11 Section 1.2 machinery,
applied verbatim, including publishing the calibrated bound that
power consumes):

    n_sup = ceil( S^2_90 * (1.960 + 1.282)^2 / 0.15^2 )   # 90% power
                                                          # at Delta = 0.15
    n_eq  = ceil( S^2_90 * (1.960 + 1.645)^2 / 0.05^2 )   # 90% P(ROBUST)
                                                          # at Delta = 0

`n_per_rung = clamp(max(n_sup, n_eq), 12, 200)` when `n_eq` fits the
cap, else `clamp(n_sup)` with CURVATURE-ROBUST declared unavailable
ex ante. **The cap is 200 here, not P11's 60** (v1.1): with the
design-check variance, `n_eq ~ 170` after calibration, so a cap of 60
would have made the experiment's most likely verdict unpurchasable —
a design that cannot buy its own expected conclusion is not a design.
The frozen power computation is unchanged; only the affordability
ceiling moves, and it moves because a P13 sample costs about 0.2 s
where P11's cost minutes. The detectable effect `0.15` dex is a design choice too: it
is roughly the size a `(tau/ell)^2` correction with an order-unity
coefficient would produce at the top rung (Section 7), stated as an
order-of-magnitude expectation and NOT frozen as a target — the
coefficient is precisely the unverified-constant class that killed
P11's Stage C v1.

Stage P-13 pilots the two ENDPOINT rungs (0.30 and 1.50), 200 samples
each, cross-rung statistics forbidden.

**What this experiment is really testing, stated plainly.** The
design checks say the contrast is zero to within `+/- 0.03`, so the
likely reading is CURVATURE-ROBUST: at fixed discreteness the
flat-normalized chain reading does not notice curvature out to
`R tau^2 = 4.5`. That is the strong form of this programme's thesis —
curvature entering the order only through the counting measure — and
P13 exists to establish or refute it with frozen power, not to
discover it. Those design numbers are quarantined (seeds `7000000+`
and `9000000+`, disjoint from every experimental window) and are
disclosed here precisely so that a later ROBUST verdict cannot be
mistaken for a surprise.

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
`estimate_tau_from_longest_chain_1p1` — never re-derived (the P12
review finding), so a recalibration of the shared definition reaches
all three experiments at once. The `rho` it is given is the P12
convention, pinned here to remove the ambiguity (v1.1): **the
realized event count divided by the rung's frozen patch proper
volume**, `ell^2 (1 - 1/|eta_lo|) * 2X`, not the nominal `rho` of the
Section 2 table (which sets the sprinkler's intensity).

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
2. **Flat twin**, operationally defined (v1.1, so the control is
   reproducible rather than gestural). For each rung: the SAME
   coordinate rectangle (`eta_lo`, `X`); a UNIFORM intensity
   `rho_twin` calibrated by design check so the realized mean `m` is
   76; a band on `tau_flat = sqrt(dU dV)` of `+/- 10%`, centred on
   the curved rung's design-check mean box area, so the
   box-to-patch ratio — the packing pressure — matches; and the same
   eligibility conditions, `K`, rejection cap, fill rule and
   `m`-gate. The true `Delta` is then zero by construction, so a
   twin contrast exceeding `delta_eq` in magnitude means the
   protocol has a size-dependent artifact and the record reports
   **CONFOUNDED**. What is NOT matched, and is published rather than
   glossed: the VARIANCE of box areas, which is larger on the curved
   side.

Labelled, never gating: a **`(tau/ell)^2` trend fit across all four
rungs** with a bootstrap CI on its slope (v1.1) — the endpoint
contrast discards the two middle rungs, and in a small-coefficient
regime the trend fit is what rescues an INCONCLUSIVE endpoint
reading; the per-rung dispersion of `m` (so the mean-matching
assumption of Section 5 is auditable rather than asserted); and the
estimator's error against the flat template together with the
geometric gap. P12's lesson stands — the first of
these is density- and geometry-dependent, so no level is frozen for
it this time; it is reported, not scored.

## 6. Seed windows (frozen; above every range used)

| block | base | fills | window span |
|---|---|---|---|
| verification-13 (non-experimental) | 6000000 | 2000 per rung | 6000000-6007999 |
| Stage P-13 pilot, tau/ell = 0.30 | 6020000 | 200 of 220 | 6020000-6063999 |
| Stage P-13 pilot, tau/ell = 1.50 | 6064000 | 200 of 220 | 6064000-6107999 |
| Stage A-13, tau/ell = 0.30 | 6120000 | <= 200 of 220 | 6120000-6163999 |
| Stage A-13, tau/ell = 0.60 | 6164000 | <= 200 of 220 | 6164000-6207999 |
| Stage A-13, tau/ell = 1.00 | 6208000 | <= 200 of 220 | 6208000-6251999 |
| Stage A-13, tau/ell = 1.50 | 6252000 | <= 200 of 220 | 6252000-6295999 |
| flat twin, tau/ell = 0.30 | 6300000 | <= 200 of 220 | 6300000-6343999 |
| flat twin, tau/ell = 0.60 | 6344000 | <= 200 of 220 | 6344000-6387999 |
| flat twin, tau/ell = 1.00 | 6388000 | <= 200 of 220 | 6388000-6431999 |
| flat twin, tau/ell = 1.50 | 6432000 | <= 200 of 220 | 6432000-6475999 |

(Re-cut in v1.1 for the raised cap: 220 slots of stride 200 per
block, 44000 seeds each, all inside `6000000-6999999`.)

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
`tau / ell`, and anything beyond `tau / ell = 1.5`. The stopping
point is a scope decision, not a cost one (v1.1 — the `O(N^2)` scan
would only roughly double at `tau/ell = 2`): `R tau^2 = 4.5` already
violates the small-diamond premise by more than a factor of two, so
the ladder answers "outside the theory's comfort zone" and closes.
If the reading is ROBUST there, the natural successor sweeps
`tau/ell in {1.5, 2.0, 2.5, 3.0}`, where the conjugate point at
`pi ell ~ 3.14` is a hard ceiling — a different experiment with a
different failure mode, not an extension of this one. The
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

---

## 9. Stage records (campaign v1)

### 9.1 Verification-13 and Stage P-13 (2026-07-28)

- **Verification-13**: 2000 / 2000 / 2000 / **1998** of 2000 complete
  across the four rungs — the pin met exactly at the top rung, so the
  v1.1 packing constants hold at pin resolution where v1.0's would
  have collapsed. Measured wall times ~0.1 s per sample.
- **Stage P-13**: 200 samples at each endpoint rung, zero skips. The
  calibration fired hard at the top rung (coverage 0.833 at nominal,
  `z = 2.86` — heavier tails than any earlier pilot), and the
  realized variance came in ABOVE the design check: `n_sup = 21`,
  **`n_eq = 231 > cap 200`**, so CURVATURE-ROBUST was declared
  unavailable ex ante and the campaign ran at `n = 21` with the
  reduced vocabulary. The cap was not raised at this point: the
  pilot's variance is data, and changing a frozen rule after seeing
  it is precisely what preregistration forbids.

### 9.2 Stage A-13 (2026-07-28): **CONFOUNDED**

21 samples per rung on both arms, zero skips.

    Delta_13 = +0.014,  95% CI [-0.050, +0.077]
    flat twin = +0.085, 95% CI [+0.032, +0.136]   <- gate 2 fails

| gate | reading | result |
|---|---|---|
| `m`-matching | 73.9 / 74.4 / 75.8 / 79.5, grand 75.9 | PASS (top rung +4.7%, inside the 5% band) |
| flat twin | `|Delta_twin| = 0.085 > delta_eq = 0.05` | **FAIL** |

The frozen rule returns **CONFOUNDED**: a control whose true contrast
is zero by construction read larger than the equivalence margin, and
the curved contrast (+0.014) is SMALLER than the artifact the control
reports, so nothing can be said about curvature from this campaign.
The verdict stands as run; gates are not re-read.

Labelled, for the record: the curved rungs' median relative errors
were 0.200 / 0.182 / 0.184 / 0.196 — flat to within 10% while
`(tau/ell)^2` moved 25-fold — and the `(tau/ell)^2` trend fit read
`+0.021`, CI `[-0.008, +0.051]`, consistent with zero. The raw signal
still points at curvature-insensitivity; the campaign simply has no
standing to say so.

### 9.3 Isolated diagnostic: the control was innocent (post hoc, quarantined)

Run in design-check space (`7500000+`, `7600000+`), never touching an
experimental window.

**Pair level (720 pairs per rung).** Relative error
0.201 / 0.192 / 0.186 / 0.192; `m = 74.5 +/- 12.2` identically at
every rung; mean chain length 16.0-16.1; aspect-ratio distributions
overlapping. **Conditioned on `m in [70, 82]`**: 0.192 / 0.184 /
0.184 / 0.185. This is what theory demands and is stronger than a
measurement: the LIS distribution of a uniform Poisson box is
invariant under axis-wise monotone rescaling, so given the count it
cannot know the box's size or aspect. The twin's implementation is
therefore not broken, and the packing-selection worry closes too,
since the `m` distributions coincide across rungs.

**Sample level (150 samples per rung, the campaign's own code path).**
Mean `y` = -0.759 / -0.739 / -0.751 / -0.744, giving
**`Delta_twin = +0.015 +/- 0.030`** — consistent with zero.

**Where +0.085 came from.** Resampling `n = 21` subsets of those 150
gives a contrast SD of 0.037 and

    P(|Delta_twin| > delta_eq | true value 0) = 0.218,
    P(Delta_twin > 0.085 | true value 0) = 0.031.

So a perfect protocol trips this control roughly one campaign in
five, and the observed value is a 3%-tail draw. The campaign's CI
excluded zero, which is why the CONFOUNDED verdict was correct as
frozen; the grounds for now calling it a fluctuation are the
150-per-rung diagnostic and the invariance theorem, not a re-reading
of the campaign.

**The defect, named.** v1.1 specified the twin's protocol but never
its SAMPLE SIZE, so the control inherited the curved arm's `n_sup`,
which was sized to detect `Delta = 0.15` — not to demonstrate that a
control sits inside `delta_eq = 0.05`. An equivalence-grade
requirement was given a superiority-grade sample. That is the same
class as P11 v1.5's "the design promised a verdict it had not
priced", recurring one level out: this time in the control rather
than the gate.

---

## 10. Design v2 (2026-07-28): the control gets its own power and semantics

Campaign v1's record above stands. v2 changes three frozen things and
nothing else, each because of a fact now on the record rather than a
preference.

**(1) The twin gate becomes an equivalence test, evaluated by
precedence.** A point-estimate threshold is a point test at any `n`,
which is why a perfect control failed it. With `lo_t`, `hi_t` the
twin's 95% bootstrap interval:

| # | condition | control result |
|---|---|---|
| 1 | `-delta_eq < lo_t` and `hi_t < delta_eq` | **CONTROL-CLEAN** — equivalence demonstrated; the campaign's verdict stands |
| 2 | `lo_t > delta_eq` or `hi_t < -delta_eq` | **CONFOUNDED** — a real artifact larger than the margin |
| 3 | otherwise | **UNDERPOWERED-CONTROL** — a record, not a verdict: the campaign reports its contrast with the control's imprecision stated |

Row 3 is the honest home for the reading v1 produced, and it cannot
fire spuriously the way a point threshold does.

**(2) The twin is piloted and sized like an arm.** Stage P-13B pilots
BOTH twin endpoint rungs (200 samples each), and

    n_twin = ceil( S^2_90(twin) * (1.960 + 1.645)^2 / delta_eq^2 ),

the same equivalence formula the curved arm uses, clamped to the same
floor and cap. The twin's measured `s(y) ~ 0.11-0.14` puts this near
`n ~ 170-230`, which is the point: demonstrating a control is clean
costs as much as demonstrating equivalence anywhere else.

**(3) The cap rises to 300, ex ante.** This is legitimate now and was
not before: `n_eq = 231` is disclosed pre-freeze information for v2,
where in v1 it was pilot data arriving after the rule was frozen. At
300 the curved arm can again buy CURVATURE-ROBUST, and the total cost
is about `2 * 231 * 4 * 0.1 s`, a few minutes.

**(4) The `m`-gate's false-firing rate is stated rather than
assumed** (same class check, applied to the other control): at
`n = 21` the sample-level spread of per-rung mean `m` made the
top rung's `+4.7%` a near miss of the `5%` band; at the v2 sample
sizes the standard error of a rung's mean `m` is under `1%`, so the
gate fires on drift rather than on noise. The realized figures are
published per rung either way.

**Fresh windows (v2, frozen).** Campaign v1's seeds are spent, so v2
runs entirely in `8000000-8799999`; design-check space remains
`7000000-7999999` and `9000000+`.

| block | base | slots | span |
|---|---|---|---|
| verification-13B | 8000000 | 2000 per rung | 8000000-8007999 |
| pilot, curved 0.30 | 8020000 | 320 | 8020000-8083999 |
| pilot, curved 1.50 | 8084000 | 320 | 8084000-8147999 |
| pilot, twin 0.30 | 8148000 | 320 | 8148000-8211999 |
| pilot, twin 1.50 | 8212000 | 320 | 8212000-8275999 |
| Stage A, curved 0.30 | 8280000 | 320 | 8280000-8343999 |
| Stage A, curved 0.60 | 8344000 | 320 | 8344000-8407999 |
| Stage A, curved 1.00 | 8408000 | 320 | 8408000-8471999 |
| Stage A, curved 1.50 | 8472000 | 320 | 8472000-8535999 |
| twin 0.30 | 8540000 | 320 | 8540000-8603999 |
| twin 0.60 | 8604000 | 320 | 8604000-8667999 |
| twin 1.00 | 8668000 | 320 | 8668000-8731999 |
| twin 1.50 | 8732000 | 320 | 8732000-8795999 |

Everything else — the ladder, the patch constants, the bands, the
estimator, both eligibility conditions, the verdict table of Section
1.2, `delta_eq`, and the detection target — is unchanged from v1.1
and is NOT re-opened by this revision.
