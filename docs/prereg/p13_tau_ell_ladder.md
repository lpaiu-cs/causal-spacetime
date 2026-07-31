# P13: how far does curvature let the flat-normalized reading go? — power-first preregistration

Status: **CLOSED. Three campaigns have run and the question is
answered.** The result is Section 13:
**CURVATURE-ROBUST with the control CLEAN**, at `delta_eq = 0.02` dex,
`Delta_13 = +0.0072` with CI `[-0.0004, +0.0148]`, and a flat twin that
reproduces the curved arm's entire residual drift to `+0.00007`. At
fixed discreteness `m ~ 76` the flat-normalized chain reading does not
notice curvature out to `tau/ell = 1.5`, i.e. `R tau^2 = 4.5`.

Getting there took three campaigns, and each verdict stands at its own
stamp, none re-read. v1 **CONFOUNDED** (Section 9): the control had no
sample size of its own. v2 **CURVATURE-DEGRADES** (Section 11): the
control was sized, but the verdict's trigger was scale-free and fired
at `+0.024` — under half the margin then in force — and 11.3's
quarantined diagnostic at four times the precision does not reproduce
that contrast, so what the word may be read to MEAN is bounded there.
v3 (Sections 12 and 13) repaired the trigger to key on the margin,
tightened `delta_eq` to 0.02, raised the equivalence power target to
99%, and ran on fresh windows at `15000000+`.

Dates are local (UTC+9); commit timestamps carry their +09:00 offset.

Sections 1 through 8 describe the design as frozen for campaign v1.
Where v2 changes a frozen rule it says so in Section 10; everything
those sections state that Section 10 does not touch remains in force
verbatim.

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

**What row 3 does to the curved verdict word, pinned so it cannot be
read two ways** (the ambiguity v1 paid for elsewhere): on row 3 **the
verdict word IS issued** and carries UNDERPOWERED-CONTROL beside it,
exactly as P11's SELECTION-CAVEAT rides alongside a verdict rather
than replacing it, and Section 9's narration must print the twin's
interval next to the verdict word every time. Withholding the word
would kill a campaign twice over for a control that is merely
imprecise; with `n_twin` sized as in (2), row 3 fires only when the
realized variance badly overruns the pilot, which is rare and is
itself worth reporting.

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

---

## 11. Stage records (campaign v2)

Artifacts: `docs/prereg/frozen/p13v2/`. All four are stamped
`d916f34`, so stamp-equality across the chain holds. Campaign v1's
artifacts stay in `docs/prereg/frozen/p13/` and the two directories
are never mixed in a quotation.

### 11.1 Verification-13B and Stage P-13B (2026-07-28)

- **Verification-13B**: 2000 / 2000 complete at every rung against
  the 1998 pin, in the fresh `8000000+` windows. Measured wall times
  ~0.2 / ~0.09 / ~0.08 / ~0.1 s per sample (one significant figure,
  per 9.6).
- **Stage P-13B**: 200 samples in each of FOUR pilot blocks — both
  curved endpoint rungs and, new in v2, both twin endpoint rungs.
  Curved: `s^2_90 = 0.0473` after calibration (bottom coverage 0.8805
  -> `z = 2.756`; top 0.929 -> `z = 1.813`), giving `n_sup = 23` and
  **`n_eq = 246`, inside the raised cap of 300**, so
  CURVATURE-ROBUST was purchasable ex ante and `n_per_rung = 246`.
  Twin: `s^2_90 = 0.0343`, `n_twin = 179`, equivalence affordable.
  Projected Stage A 0.05 h. FEASIBLE.

The cap is what bought the vocabulary. At v1's cap of 200 this
pilot's `n_eq = 246` would again have declared ROBUST unavailable and
the campaign would again have run at `n_sup = 23`.

### 11.2 Stage A-13 v2 (2026-07-28): **CURVATURE-DEGRADES**, control CLEAN, effect at half the margin

246 curved samples and 179 twin samples per rung, zero skips on
either arm, no selection caveat.

    Delta_13  = +0.0243,  95% CI [+0.0044, +0.0456]  ->  CURVATURE-DEGRADES
    flat twin = +0.0063,  95% CI [-0.0194, +0.0328]  ->  CONTROL-CLEAN

| gate | reading | result |
|---|---|---|
| `m`-matching | 75.57 / 75.58 / 75.98 / 75.89, grand 75.755 | PASS, widest rung at +0.30% |
| flat twin (v2 equivalence test) | interval inside `+/- delta_eq` | **CONTROL-CLEAN** |

Both v2 controls did what v2 was for. The `m`-gate now fires on drift
rather than on noise — the widest rung sits at `+0.30%` where v1's
`n = 21` put it at `+4.7%`, which is 10 (4) discharged as promised.
And the twin, piloted and sized as an arm, demonstrated equivalence
instead of failing a point threshold.

**The magnitude, printed beside the verdict word because 1.2 requires
it.** `Delta = +0.0243` is under HALF of `delta_eq = 0.05`, and even
the interval's upper end (`+0.0456`) lies inside the margin. The
campaign therefore says: the effect is statistically detected — zero
is excluded — and it is SMALLER than the threshold below which "the
reading did not notice the curvature" is a fair description. Quoting
the verdict word alone inverts the finding. This is the exact
situation 1.2 anticipated when it pinned that rows 1 and 2 can fire
on a contrast smaller than `delta_eq`.

**Where the contrast comes from, and it is not a trend.** The
labelled checks: median relative error `0.1795 / 0.1916 / 0.1880 /
0.1860`, a 6% spread while `(tau/ell)^2` moves 25-fold; the
`(tau/ell)^2` trend fit reads `+0.0073`, CI `[-0.0013, +0.0158]`,
consistent with zero. Both hold at once because of the rung profile:

| `tau/ell` | mean `y`, curved (n = 246) | mean `y`, twin (n = 179) |
|---|---|---|
| 0.30 | -0.7608 +/- 0.0075 | -0.7552 +/- 0.0097 |
| 0.60 | -0.7367 +/- 0.0071 | -0.7373 +/- 0.0095 |
| 1.00 | -0.7392 +/- 0.0071 | -0.7508 +/- 0.0078 |
| 1.50 | -0.7364 +/- 0.0071 | -0.7489 +/- 0.0088 |

(`+/-` is one standard error. The curved rung means are in the
summary artifact; the twin rung means and every standard error are
computed from `p13_stage_a.csv`, the same frozen rows the gate
scored. A labelled addition, not a re-reading of the gate.)

Curved rungs 2, 3 and 4 agree to 0.003 dex. The `0.30 -> 0.60` step
is `+0.0241 +/- 0.0103` and the primary contrast is `+0.0243`:
**essentially the entire contrast is a bottom-rung offset.** A
`(tau/ell)^2` correction does not have that shape.

**And the control has the same offset.** The twin's `0.30 -> 0.60`
step is `+0.0179 +/- 0.0136` — same sign, comparable size, in an
ensemble whose true contrast is ZERO by construction. The twin's
endpoint contrast stays small (`+0.0063`) only because its rungs 3
and 4 come back down; the feature that produces the curved arm's
signal is present in the control. Campaign v1's own 9.3 twin
diagnostic, at 150 samples per rung in the disjoint design-check
space `7600000+`, recorded `-0.759 / -0.739 / -0.751 / -0.744`: the
same shape, so the offset replicates across two independent seed
spaces.

**The control is clean at its own resolution and silent at the
effect's.** `n_twin = 179` was sized by the equivalence formula
against `delta_eq = 0.05`, exactly as 10 (2) specifies, and its
realized interval half-width is 0.026. The verdict, however, fires on
`lo > 0`, which has no size floor. So the gate rang at `+0.0243`
while the control can only certify the absence of artifacts at the
`0.05` scale, and the twin's interval `[-0.0194, +0.0328]` contains
`+0.0243` comfortably. **CONTROL-CLEAN is a true statement about
0.05-scale artifacts and no statement at all about 0.024-scale
ones.**

That is this programme's recurring defect class, one level further
out again. P11 v1.5: the design promised a verdict it had not priced.
P13 v1: an equivalence-grade requirement was given a
superiority-grade sample. P13 v2: **the verdict fires at a resolution
finer than the control was sized to certify.** v2 gave the control
power against `delta_eq` and left the trigger scale-free, so the
mismatch moved rather than closed. Naming it here is not a
re-reading of the gate; it is a bound on what the gate's word can be
taken to mean.

**Netting out the control** — a labelled calculation, not the frozen
statistic, since P13 freezes the twin as a gate and not as a
subtrahend:

    endpoint difference-of-differences    +0.0180 +/- 0.0167   (1.1 sigma)
    rung-1 step difference-of-differences +0.0062 +/- 0.0171   (0.4 sigma)

Nothing survives either subtraction. The difference construction is
in-programme — P12 Section 5 item 5 makes exactly this subtraction
Stage B's primary, for exactly this reason — but P13 did not freeze
it, so it is reported and not scored.

**What the campaign establishes.** The frozen gate says
CURVATURE-DEGRADES and it stands: the contrast excludes zero at 246
samples per rung with both controls passing as frozen, and the
verdict is not re-read. What the campaign does NOT establish is that
curvature is the cause of it, and the reason is on the record above
rather than in an objection: the whole contrast is a bottom-rung
offset, and a zero-effect ensemble shows the same offset at
comparable size.

**Three candidates for the bottom-rung offset, and what the frozen
artifact already decides.**

1. **The band's absolute width, entering through the `m`
   distribution.** The `+/- 10%` band is proportional, so its
   absolute width grows with `tau`. **Excluded by the artifact**:
   `sd_m = 11.52 / 11.74 / 11.53 / 11.59`, flat to 2% across the
   ladder. And quantitatively — P12's measured `dy / dlog10 m ~
   -0.304` makes a 2% shift in `m` worth about 0.003 dex, an order
   below `+0.024`. No run was required to close this one.
2. **The 0.30 rung's geometry at fixed `tau/ell`.** The candidate the
   twin supports, and the quantity is packing pressure, which is not
   flat across the ladder. With `K = 6` boxes of `(u,v)` area
   `centre^2` in a patch of `(u,v)` area `4 (|eta_lo| - 1) X`, the
   frozen constants give **0.168 / 0.222 / 0.263 / 0.263**. Two
   notes, both corrections to the natural guess: the bottom rung has
   the LOOSEST packing, not the tightest (Section 2's steeply growing
   `X` is about the top rungs needing room, and the bottom rung ended
   up with a surplus); and the relation to `y` is a step, not a
   proportionality, since rung 2's ratio is intermediate while its
   `y` already sits on the plateau. A correlation across rungs
   therefore cannot settle it — only moving the ratio at fixed
   `tau/ell` can, which is 11.3.
3. **A real curvature effect that is not of `(tau/ell)^2` form.**
   Still available. It has to beat candidate 2, and it has to explain
   why a flat ensemble shows the same step.

Candidate 2 has direct control-arm support, so it gets an isolated
diagnostic, on the v1 precedent of 9.3: quarantined seed space, the
campaign's own code path, labelled throughout, and no gate re-read.

### 11.3 Isolated diagnostic: the offset does not replicate, and the mechanism has no effect

Post hoc, quarantined, labelled. Design-check seeds `9000000+` per
Section 10, no experimental window touched, no gate re-read. The
script is `docs/prereg/frozen/p13v2/p13v2_rung1_diag.py`, committed at
`7b1442d` with its predictions stated **before** it ran, and the
results at that same stamp are in
`p13v2_rung1_diag_results.json`. It calls `run_sample`, so
eligibility, packing, the `rho` convention and the estimator are the
campaign's. Every configuration completed 1.000 and every `rho`
calibration converged, so nothing below is a selection or a
discreteness artifact.

**D1 and D2: the ladder at four times the campaign's precision.**
1000 samples at the two rungs that carry the step, 500 at the
plateau, both arms.

| `tau/ell` | curved, campaign (n=246) | curved, diagnostic | flat, campaign (n=179) | flat, diagnostic |
|---|---|---|---|---|
| 0.30 | -0.7608 +/- 0.0075 | **-0.7432 +/- 0.0036** | -0.7552 +/- 0.0097 | -0.7533 +/- 0.0036 |
| 0.60 | -0.7367 +/- 0.0071 | **-0.7514 +/- 0.0039** | -0.7373 +/- 0.0095 | -0.7485 +/- 0.0039 |
| 1.00 | -0.7392 +/- 0.0071 | -0.7368 +/- 0.0049 | -0.7508 +/- 0.0078 | -0.7520 +/- 0.0053 |
| 1.50 | -0.7364 +/- 0.0071 | -0.7376 +/- 0.0051 | -0.7489 +/- 0.0088 | -0.7455 +/- 0.0055 |

The plateau rungs agree between the two samples to within a standard
error. **The two bottom rungs swapped.** In the campaign rung 0.30 sat
0.024 dex BELOW rung 0.60; in the diagnostic it sits 0.008 dex ABOVE
it.

| contrast | campaign | diagnostic | separation |
|---|---|---|---|
| curved, `0.30 -> 1.50` (the primary) | +0.0243 +/- 0.0105 | **+0.0056 +/- 0.0062**, CI `[-0.0066, +0.0179]` | 1.5 sigma |
| curved, `0.30 -> 0.60` (the step) | +0.0241 +/- 0.0103 | **-0.0082 +/- 0.0054** | 2.8 sigma |
| flat, `0.30 -> 1.50` | +0.0063 +/- 0.0133 | +0.0079 +/- 0.0066 | 0.1 sigma |
| flat, `0.30 -> 0.60` | +0.0179 +/- 0.0136 | +0.0048 +/- 0.0053 | 0.9 sigma |
| endpoint difference-of-differences | +0.0180 +/- 0.0167 | **-0.0022 +/- 0.0091** | — |

So both offsets 11.2 reported were sampling noise at their own sample
sizes: the twin's `+0.0179` falls to `+0.0048 +/- 0.0053`, and the
curved arm's `+0.0241` — which was the entire primary contrast —
reverses to `-0.0082 +/- 0.0054`. The endpoint contrasts are formally
compatible between the two samples (1.5 sigma); the steps are in mild
tension (2.8 sigma). Either way the higher-precision sample puts the
primary contrast at about a quarter of the campaign's value with an
interval that excludes nothing.

Across the diagnostic's curved ladder the total spread is **0.0146
dex** while `(tau/ell)^2` moves 25-fold, from 0.09 to 2.25.

**D3: the suspected mechanism, measured directly.** Packing pressure
moved at FIXED `tau/ell`, `rho` recalibrated so realized `m` stayed on
the campaign's grand mean 75.755 (`m` turned out to be nearly
`X`-insensitive on the curved arm — both curved calibrations converged
at their frozen `rho` on the first probe, which is what a count in a
band-fixed box should do).

| variant | move | packing ratio | curved shift | flat shift |
|---|---|---|---|---|
| V-LOOSE | rung 0.60, `X` 4.0 -> 5.28 | 0.2217 -> 0.1679 | **+0.0008 +/- 0.0064** | +0.0061 +/- 0.0064 |
| V-TIGHT | rung 0.30, `X` 1.8 -> 1.36 | 0.1680 -> 0.2224 | **-0.0002 +/- 0.0065** | +0.0070 +/- 0.0063 |

Realized `m`: 76.33 / 74.86 / 75.32 / 74.52. Completion 1.000
everywhere, including the narrowed patch.

On the curved arm, a `+/- 32%` move in packing pressure shifts `y` by
under 0.001 dex, bounding the mechanism at `+/- 0.013` dex at 95%. On
the flat arm both shifts are about `+1 sigma` and — decisively — they
have the SAME sign under OPPOSITE moves of the ratio, which is what
noise looks like and is not what a mechanism looks like. The
docstring's two-sided prediction (the variants CROSS if packing
pressure drives the offset) is refuted on both arms.

**All three candidates are now closed.**

1. The `m` distribution: closed by the frozen artifact before any run
   (`sd_m` flat to 2%; P12's measured slope makes a 2% `m` shift worth
   0.003 dex).
2. Packing pressure: closed by D3 above. Note this closes the
   mechanism independently of the offset — even had the offset been
   real, this is not what caused it.
3. A real curvature effect not of `(tau/ell)^2` form: nothing is left
   to explain. The high-precision curved ladder is flat to 0.0146 dex
   across a 25-fold change in `(tau/ell)^2`.

**What this means for the campaign's word.** The gate stands exactly
as frozen and is not re-read: `Delta_13 = +0.0243`, CI
`[+0.0044, +0.0456]`, CURVATURE-DEGRADES, control CLEAN, `m`-gate
passed. What the diagnostic removes is the grounds for reading that
word as a curvature effect. This is the P12 9.2 pattern repeating with
a different number: a campaign-sized signature that a quarantined
high-`n` diagnostic dissolves, and the honest record says so rather
than letting the verdict word carry a mechanism it cannot support.

**The substantive reading, labelled and not a verdict.** At fixed
discreteness `m ~ 76`, no degradation of the flat-normalized chain
reading is detected out to `tau/ell = 1.5`, i.e. `R tau^2 = 4.5`, with
the endpoint contrast bounded inside `[-0.0066, +0.0179]` dex at 95%
— well inside `delta_eq = 0.05`. That interval satisfies the
numerical condition of row 3 of the 1.2 table, and it is the strong
form of this programme's thesis: curvature entering the order only
through the counting measure. **But a diagnostic cannot award a
verdict.** It carries no verification pin, its `n` is not
pilot-derived, and it ran after the campaign's data were seen.
CURVATURE-ROBUST remains unawarded, and awarding it requires a
campaign.

**The class, stated precisely this time.** Rows 1 and 2 are
significance tests with no size floor, so a ladder with no true effect
fires one of them in about one campaign in twenty by construction —
that is the `alpha` of the test, not a defect in it. The campaign drew
one. What the design did not require is that the trigger clear the
resolution of the campaign's own controls: `n_per_rung = 246` gives
`se ~ 0.010` on the contrast, `n_twin = 179` gives `+/- 0.026` on the
control, and row 1 fired at `+0.024` — below the control's resolution
and at twice the arm's. So the three recurrences read: P11 v1.5
promised a verdict it had not priced; P13 v1 gave an equivalence-grade
requirement a superiority-grade sample; P13 v2 let a scale-free
trigger fire below the scale its controls could certify. Section
1.2's magnitude pin is what kept this from being published as a
mechanism, and it earned its place.

**Named for the next design review, adopted by nothing here.** Two
options, neither frozen: size the campaign on an equivalence
requirement at a `delta_eq` near the effect scale that actually
matters (~0.02 rather than 0.05), which at this variance costs
`n ~ 1500` per rung and a few minutes; or give rows 1 and 2 a size
floor, so a significant contrast must also exceed some stated fraction
of `delta_eq` before it earns the word DEGRADES. The first buys
precision, the second buys honesty in the vocabulary. Choosing between
them is a design decision, and this record does not make it.

---

## 12. Design v3 (2026-07-31): one margin, and the trigger keyed to it

Approved in session. 11.3 answered P13's scientific question at
diagnostic grade and left the design owing one repair and offering one
upgrade; both are taken, and the power target moves with them.
Campaign v2's record stands and is not re-scored. Everything below is
frozen before any v3 datum exists, and v3 runs on fresh windows.

### 12.1 The repair: the trigger and the margin become the same number

Section 10 (1) replaced the twin's point threshold with a three-way
equivalence test, because a point test is a point test at any `n`. It
did not carry that fix to the PRIMARY verdict table, and the module
shows the asymmetry with the two functions side by side: the control
asks `lo_t > delta_eq`, the verdict asks `lo > 0`. That gap IS 11.3's
third recurrence, so the repair is to propagate a fix already made
rather than to invent one.

Frozen v3 table, by precedence as before:

| # | condition | verdict |
|---|---|---|
| 1 | `lo > delta_eq` | **CURVATURE-DEGRADES** |
| 2 | `hi < -delta_eq` | **CURVATURE-HELPS** |
| 3 | `-delta_eq < lo` and `hi < delta_eq` | **CURVATURE-ROBUST** (available only if the pilot declared it affordable) |
| 4 | otherwise | **INCONCLUSIVE** if row 3 was affordable, else **UNRESOLVED** |

The vocabulary now matches the question. A contrast that does not clear
the margin cannot earn the word DEGRADES, because "does not clear the
margin" is the definition of the reading not noticing the curvature.
Row 4 widens to hold every interval that straddles a margin edge,
which is the honest home for an imprecise campaign.

**1.2's magnitude pin STAYS.** It is now belt and braces rather than
the only guard, and 11.2 is the reason to keep it: it is what stopped
a scale-free trigger from being published as a mechanism.

**What this table may not be used for.** Applied to campaign v2's
frozen interval `[+0.0044, +0.0456]` it returns CURVATURE-ROBUST. That
is a demonstration that the repair targets the real defect, and it is
**not** a re-scoring. v2's verdict was issued under v2's rules and
stands as CURVATURE-DEGRADES. A rule change that flips one's own
verdict may only touch data not yet seen; that is why v3 gets fresh
windows and why this paragraph is here rather than discovered later.

### 12.2 The upgrade: `delta_eq` falls from 0.05 to 0.02 dex

A margin should sit near the scale at which the answer would change.
When v1 chose `delta_eq = 0.05` dex -- 12% in relative error -- that
was where the answer would have changed, and it was stated as a design
choice rather than a derived quantity. 11.3 moved that scale: the
endpoint contrast is now bounded inside `[-0.0066, +0.0179]` dex, so a
0.05 margin is three times wider than the best available bound and
would certify something no longer in doubt.

`delta_eq = 0.02` dex is 4.7% in relative error. Since the estimator's
own error is 18-19% at `m ~ 76`, ROBUST at this margin says curvature
moves that error by less than 4.7% OF ITSELF while `(tau/ell)^2` moves
25-fold. That is a materially stronger claim than v2 could have made,
and it is the claim the diagnostic suggests is true.

**Affordability is not the justification.** A sample costing 0.1 s is
what makes the stronger claim purchasable; it is not what makes 0.02
the right margin. The argument for the number is the one above, and it
would hold at any cost per sample.

**Pre-freeze disclosure, per the v1.1 and P12 precedent.** The
quarantined numbers that make ROBUST the likely v3 reading are
11.3's, in design-check space `9000000+`: curved endpoint
`+0.0056 +/- 0.0062`, interval `[-0.0066, +0.0179]`, already inside
`+/- 0.02`; the curved ladder flat to 0.0146 dex; packing pressure
bounded at `+/- 0.013` dex. They are disclosed here so that a v3
ROBUST cannot later be mistaken for a surprise, and so that this
design is on record as having been built to establish a reading it
expects rather than to discover one.

### 12.3 Power: the equivalence target rises from 90% to 99%

`n_eq` is defined so that `P(ROBUST | Delta = 0) = 90%`. That 10% is
not slack -- it is exactly the probability that a null ladder's
interval pokes outside the margin. Campaign v2 drew a `+2.3 sigma`
fluctuation and its upper end landed 0.0044 short of the margin:
inside, but by a hair. **That thinness is the `beta`, and it scales
with `delta_eq`, so tightening the margin does not fix it.** Only power
does.

    n_sup = ceil( S^2_90 * (1.960 + 1.282)^2
                  / (DELTA_DETECT - delta_eq)^2 )    # 90% power
    n_eq  = ceil( S^2_90 * (1.960 + 2.326)^2
                  / delta_eq^2 )                     # 99% P(ROBUST)
    n_twin = ceil( S^2_90(twin) * (1.960 + 2.326)^2 / delta_eq^2 )

`n_sup` now measures against `DELTA_DETECT - delta_eq`, because row 1
requires clearing the margin rather than clearing zero. At this
variance it is about 30 and is not binding; it is corrected for
accuracy, not for effect.

**The cap rises from 300 to 3000, ex ante.** This has the same
standing v2's raise from 200 to 300 had: the variance that sets `n_eq`
is disclosed pre-freeze rather than arriving as pilot data after the
rule is frozen. Projecting from v2's calibrated pilot bounds,
`S^2_90 = 0.0473` gives `n_eq ~ 2172` and
`S^2_90(twin) = 0.0343` gives `n_twin ~ 1575`; the cap absorbs a
realized variance up to about `S^2_90 = 0.065` before ROBUST becomes
unpurchasable again. Projected campaign cost is roughly 24,000
samples, well under an hour.

**One margin, three consumers.** `delta_eq` now sets the verdict
trigger, the control's gate, the curved arm's size and the twin's size
-- and that identity is the point. v2's resolution mismatch existed
precisely because the trigger keyed to 0 while the control keyed to
0.05. A test asserts the single-margin invariant, so the mismatch
cannot reappear by editing one constant and forgetting the other.

### 12.4 Fresh windows (frozen)

Campaign v1 spent `6000000-6475999` and v2 spent `8000000-8795999`;
11.3's diagnostic consumed design-check space up to `14550000`. So v3
runs entirely at `15000000+`, and design-check space for anything after
this is `22000000+`.

| block | base | slots | span |
|---|---|---|---|
| verification-13C | 15000000 | 2000 per rung, consecutive | 15000000-15007999 |
| pilot, curved 0.30 | 15100000 | 320 | 15100000-15163999 |
| pilot, curved 1.50 | 15200000 | 320 | 15200000-15263999 |
| pilot, twin 0.30 | 15300000 | 320 | 15300000-15363999 |
| pilot, twin 1.50 | 15400000 | 320 | 15400000-15463999 |
| Stage A, curved 0.30 | 16000000 | 3020 | 16000000-16603999 |
| Stage A, curved 0.60 | 16700000 | 3020 | 16700000-17303999 |
| Stage A, curved 1.00 | 17400000 | 3020 | 17400000-18003999 |
| Stage A, curved 1.50 | 18100000 | 3020 | 18100000-18703999 |
| twin 0.30 | 18800000 | 3020 | 18800000-19403999 |
| twin 0.60 | 19500000 | 3020 | 19500000-20103999 |
| twin 1.00 | 20200000 | 3020 | 20200000-20803999 |
| twin 1.50 | 20900000 | 3020 | 20900000-21503999 |

Slots are `cap + 20` at stride 200, the v2 pattern scaled to the new
cap, so a Stage A rung can fill 3000 completions with reserve.

### 12.5 What v3 does NOT change

The ladder and its four rungs; every per-rung patch constant; the
`+/- 10%` bands; the twin's construction, band centres and intensities;
the estimator, called and never re-derived; both eligibility
conditions; `K = 6`; the rejection cap; the first-`n`-complete fill
with published skip identities; the `>= 1998 of 2000` verification pin;
clean-worktree preflight and stamp equality; the `m`-gate at `+/- 5%`
of the grand mean; `DELTA_DETECT = 0.15`; the three-way control
semantics of 10 (1) including UNDERPOWERED-CONTROL riding beside a
verdict; and 1.2's requirement that magnitude be printed with the
word.

### 12.6 A consequence, stated once

`delta_eq` and the cap are module constants, and the module's comment
calls the literals the frozen spec. Changing them means v1's and v2's
verdicts are reproducible only at their own stamps -- `3f46313` and
`d916f34` -- which is exactly why every artifact carries
`code_version`. Nothing in Sections 9 or 11 is re-run, re-read or
re-scored by this amendment.

---

## 13. Stage records (campaign v3)

Artifacts: `docs/prereg/frozen/p13v3/`. Every file, and every row of
the CSV, carries `code_version = 33e371f`, so stamp equality holds
across the whole chain. The three campaigns keep three directories --
`p13/`, `p13v2/`, `p13v3/` -- and are never mixed in a quotation.

### 13.1 Verification-13C and Stage P-13C (2026-07-31)

- **Verification-13C**: 2000 / 2000 / **1999** / 2000 of 2000 complete
  against the 1998 pin, in the fresh `15000000+` windows. Measured wall
  times ~0.2 / ~0.1 / ~0.08 / ~0.1 s per sample (one significant
  figure, per 9.6).
- **Stage P-13C**: 200 samples in each of four pilot blocks. Curved
  `S^2_90 = 0.0377` after calibration, which fired on both endpoints
  and hard at the top: coverage 0.9335 -> `z = 1.755` at the bottom,
  **0.8055 -> `z = 2.869`** at the top — the largest under-coverage
  the programme has recorded, past v1's 0.833 and P12's 0.865, and
  the reason P11 v1.9 replaced the nominal bound with a calibrated
  one. `n_sup = 24`, **`n_eq = 1731`, inside the cap of 3000**, so
  CURVATURE-ROBUST was purchasable and `n_per_rung = 1731`. Twin
  `S^2_90 = 0.0359`, `n_twin = 1650`, equivalence affordable.
  Projected Stage A 0.4 h. FEASIBLE.

The realized variance came in BELOW 12.3's projection — `n_eq = 1731`
where v2's bounds projected 2172 — so the 99% target and the tightened
margin were bought with room to spare rather than at the cap.

### 13.2 Stage A-13C (2026-07-31): **CURVATURE-ROBUST**, control CLEAN

1731 curved samples and 1650 twin samples per rung.

    Delta_13  = +0.0072,  95% CI [-0.0004, +0.0148]  ->  CURVATURE-ROBUST
    flat twin = +0.0071,  95% CI [-0.0007, +0.0153]  ->  CONTROL-CLEAN

| gate | reading | result |
|---|---|---|
| `m`-matching | 75.64 / 75.79 / 76.12 / 76.00, grand 75.89 | PASS, widest rung +0.31% |
| flat twin, equivalence at `delta_eq = 0.02` | interval inside `+/- 0.02` | **CONTROL-CLEAN** |
| skips | 0 / 0 / 2 / 1 | selection caveat rides along |

**The selection caveat, with its size.** Three samples of 6924 failed
to pack and were replaced from reserve slots; their seeds are published
in the artifact (`17454000`, `17600600`, `18342600`). The flag fires on
any skip at all, so its magnitude is stated for the same reason 1.2
requires the contrast's: 3 in 6924 is 0.04%.

**This is the verdict the programme has been trying to buy since v1,
and v3 is the first campaign that could afford it.** v1 could not
(cap 200 against `n_eq = 231`); v2 could, and a scale-free trigger
pre-empted it. Here the interval sits inside a margin two and a half
times tighter than either used.

**What makes it strong is the control, not the verdict word.**

| `tau/ell` | mean `y`, curved (n = 1731) | mean `y`, twin (n = 1650) |
|---|---|---|
| 0.30 | -0.7499 +/- 0.0028 | -0.7462 +/- 0.0030 |
| 0.60 | -0.7470 +/- 0.0029 | -0.7459 +/- 0.0029 |
| 1.00 | -0.7443 +/- 0.0028 | -0.7453 +/- 0.0029 |
| 1.50 | -0.7427 +/- 0.0027 | -0.7391 +/- 0.0028 |

    curved endpoint contrast   +0.00720
    twin endpoint contrast      +0.00713
    difference-of-differences  +0.00007 +/- 0.00562   (0.01 sigma)

Both ladders drift monotonically upward by about `+0.007` dex, and the
flat ensemble — where curvature cannot act — reproduces the curved
arm's drift to **seven parts in a hundred thousand of a dex**. The
residual is therefore a property of the protocol, not of the geometry,
and this time the control has the precision to say so: `n_twin = 1650`
against v2's 179. (`+/-` is one standard error; rung means are in the
summary artifact, standard errors computed from the frozen CSV.)

**How close v2's rule came to the wrong word again.** The contrast's
lower bound is `-0.00039`. Under v2's table row 1 fired whenever
`lo > 0`, so a shift of four ten-thousandths of a dex would have
returned CURVATURE-DEGRADES on a contrast the flat control reproduces
to `+0.00007`. The twin's own lower bound, `-0.00068`, sat just as
close to the line. That is 12.1's repair vindicated by a near miss
rather than by an argument: under the old rule this campaign was
roughly a coin flip from publishing the same false word a second time,
and what separated them was less noise than a single sample carries.

**Labelled checks.** The `(tau/ell)^2` trend fit reads `+0.0030`, CI
`[-0.0001, +0.0062]` — at the edge of zero, an order below `delta_eq`,
and consistent with the same protocol drift the twin carries. Median
relative error `0.1839 / 0.1858 / 0.1878 / 0.1865`, a 2% spread while
`(tau/ell)^2` moves 25-fold. `sd_m` = `11.74 / 11.71 / 11.65 / 11.16`,
flat to 5%.

**Consistency with the disclosure.** 12.2 published 11.3's quarantined
bound, `[-0.0066, +0.0179]`, precisely so a ROBUST reading could not be
mistaken for a surprise. The campaign returned `[-0.0004, +0.0148]`.
The design predicted its own outcome and then bought it under frozen
rules, which is the only sense in which this result is not news.

### 13.3 What P13 concludes

Over the tested ladder, at fixed discreteness `m ~ 76`:

**The flat-normalized chain reading does not notice curvature out to
`tau/ell = 1.5`, that is `R tau^2 = 4.5`.** The endpoint contrast is
bounded inside `[-0.0004, +0.0148]` dex — under 3.5% in relative
error — against a margin of 0.02 dex, with a flat control that
reproduces the entire residual to `+0.00007` and an `m`-gate flat to
0.31%.

This is the strong form of the programme's thesis: curvature enters the
order only through the counting measure. P12 established that the
unchanged P11 estimator reads CURVED proper time at `tau/ell ~ 0.3` and
improves with density at the longest-chain rate. P13 establishes that
the same reading survives a 25-fold growth in `(tau/ell)^2`, well past
the point where the small-diamond premise the chain law rests on has
itself failed by more than a factor of two.

Scope, carried and not smoothed:

- the claim binds to this estimator family at this operating point,
  `m ~ 76`. A floor appearing only at larger `m` would be invisible
  here, as Section 7 said in advance;
- it binds to `tau/ell <= 1.5`. The successor sweeps
  `{1.5, 2.0, 2.5, 3.0}`, where the conjugate point at `pi ell ~ 3.14`
  is a hard ceiling — a different experiment with a different failure
  mode, not an extension of this one;
- the `+0.007` residual drift is real at the edge of significance and
  is NOT explained here beyond "the flat arm has it identically". It is
  a protocol property of the ladder; naming its mechanism would take
  its own design, and 11.3 already closed the two candidates that were
  on the table;
- three samples of 6924 were replaced from reserve slots.

The three campaigns' verdicts stand as issued, each at its own stamp:
v1 **CONFOUNDED** at `3f46313`; v2 **CURVATURE-DEGRADES** at
`d916f34`, with 11.3 bounding what that word may be read to mean; v3
**CURVATURE-ROBUST** at `33e371f`. None is re-read, and the reason
three campaigns were needed is written down in each: a control that
was not sized, then a trigger that was not scaled, then neither.
