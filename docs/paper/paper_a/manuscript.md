# An operational reconstruction ladder for spacetime quantities from causal order

**Draft v0.4.** Causal Spacetime Lab. Every quantitative result is tied to a
producing experiment and expected output path. Five figure-source summaries are
committed; the remaining cited legacy tables require regeneration before
submission, while the Section 6 capstone cites only committed,
provenance-locked artifacts (Section 10 and `artifact_manifest.md`). No number
is from memory.

## Abstract

We study, in controlled 1+1D (and, for dimension, higher-D) Minkowski models,
which spacetime quantities can be operationally reconstructed from a causal
(accessibility) order once a minimal, explicitly declared set of additional
ingredients is supplied. We organize the reconstructions as a ladder. In flat
Poisson-sprinkled Alexandrov intervals, order statistics estimate dimension,
and raw longest-chain length is an uncalibrated timelike statistic; converting
it to proper time requires density and dimension-dependent normalization.
Adding a global event-density calibration turns Alexandrov interval cardinality
into a timelike proper-time estimate whose error is consistent with
finite-sampling noise at the tested settings. Adding an observer chain with clock
labels yields radar time and unsigned radar distance; adding an orientation
reference lifts the reflection degeneracy to signed coordinates and lets one
recover the Lorentz map between two inertial protocols; adding overlapping
charts yields an observer atlas with approximately consistent Poincare
transition maps; adding supplied measure information, implemented in exp19 as
local weights, makes volume reconstruction possible under the conformal
ambiguity, and density-rescaled reconstruction is stable in the tested
random-thinning protocol. We also give a Rindler horizon analogue, in which an
accelerated observer's two-way radar reconstruction is confined to the expected
wedge, and a finite-speed lattice counterexample showing that finite signal
speed alone does not produce Lorentzian structure. As a capstone we validate,
in a preregistered two-claim design on a 3+1D vacuum plane wave whose
coordinate volume form and sprinkling law are exactly flat, that the
conformal-class information order carries is detectable at finite density:
both the paired mean shift of the global relation fraction and a single-poset
classifier confirm against frozen margins. The contribution is not any
single reconstruction — most are standard — but the explicit accounting of what
each rung requires and the negative results that bound it. We make no claim
that spacetime is reducible to causal order; the reconstructions are controlled
validations inside known models.

## 1. Introduction

A conservative reading of the causal-set and order-first programs is that, for
future- and past-distinguishing Lorentzian spacetimes, causal or chronological
structure determines the conformal class. Supplying a volume element fixes the
remaining conformal factor, up to diffeomorphism and unit conventions
[@blms1987; @malament1977; @hkm1976; @kronheimer1967; @sorkin2005; @surya2019].
This paper treats that statement operationally and quantitatively. Rather than
ask whether spacetime "is" a causal set, we ask: fix a causal (or
null-accessibility) order in a known model; then, supplying an explicitly
declared ingredient set, which spacetime quantities can a finite procedure
actually reconstruct, with what error, and where does each reconstruction stop?

The value of posing it this way is an explicit *ledger*. Many individual
reconstructions here are standard — proper time from longest chains
[@brightwell1991], volume from interval cardinality [@blms1987], dimension from
ordering fractions [@myrheim1978; @meyer1988], radar coordinates from an
observer protocol. What is easy to lose, and what we make explicit, is exactly
which extra structure each requires and what remains underdetermined without it.
The ladder also isolates three instructive negative results: causal order alone
does not fix conformal scale; a single observer yields only unsigned distance
(a reflection degeneracy); and finite signal speed alone does not imply
Lorentzian structure.

Contributions: (1) a reconstruction ladder that indexes recoverable spacetime
quantities by the minimal supplied ingredient set, with per-rung error behavior in
controlled models; (2) grounded validations of dimension, proper time, radar
decomposition, Lorentz-map and atlas consistency, and measure-dependent volume;
(3) three bounding negative results (conformal scale, reflection degeneracy,
finite-speed lattice); (4) a Rindler reconstruction-inaccessibility
analogue; and (5) a preregistered 3+1D capstone validation that a pure-Weyl
deformation of the light cones — at identical coordinate volume form and
sprinkling law — measurably moves the causal order at finite density
(Section 6). All results are controlled validations, not evidence that geometry
reduces to order.

## 2. The reconstruction ladder

We separate the *primitive* from *supplied ingredients* and index rungs by the
minimal ingredient each requires.

Primitive:

- **Causal / accessibility order.** A strict partial order on events (in the
  continuum models, Minkowski causal precedence; null-inclusive).

Supplied ingredient symbols (these form branches, not one linear hierarchy):

- **M — global event-density calibration.** The conversion between event count
  and physical volume in the declared sampling model.
- **O — observer chain with clock labels.** A timelike chain carrying ordered
  tick labels (an operational clock), defining a radar protocol.
- **R — oriented synchronized reference.** A second beacon chain with calibrated
  separation, fixing a spatial side.
- **A — observer atlas.** Multiple overlapping observer charts.
- **W — local measure profile / weights.** The physical volume element relative
  to the support sampling measure.

Rungs (minimal ingredient -> what is reconstructed):

| Rung | Ingredient | Reconstructed | Bounded by |
| --- | --- | --- | --- |
| R0 | order only | flat-sprinkling dimension estimate; uncalibrated chain statistics | manifold/model assumptions; no metric scale |
| R1 | order + M | timelike proper time (interval cardinality -> volume -> tau) | needs a density |
| R2 | order + O | radar time; unsigned radar distance | sign undetermined |
| R3 | order + O + R | signed coordinates; Lorentz map between protocols | supplied orientation |
| R4 | order + O + R + A | atlas transition maps (Poincare); invariant agreement | oriented, calibrated charts |
| R5 | order + M + W | volume under conformal ambiguity; coarse-graining stability | global and local measure supplied |

Three negative results bound the ladder from below (Section 5): conformal scale
is not fixed by order alone (motivating M and W); a single observer gives only
unsigned distance (motivating R); and finite signal speed alone does not give
Lorentzian structure. Figure 1 summarizes the ladder.

![Figure 1](figures/fig1_ladder.png)

*Figure 1. The reconstruction dependency ledger. Each row names a minimal
ingredient set and the quantity it unlocks; the observer and measure rows are
branches rather than one cumulative linear chain. The "bounded by" line records
what that ingredient set still does not supply.*

## 3. Methods (verified foundation layer)

Natural units, c = 1. The foundation modules were checked for numerical
correctness against known results (sprinkling measure, Minkowski causal
precedence, longest-chain and interval-cardinality normalizations,
Myrheim-Meyer inversion, radar/Lorentz/Rindler formulas); this section states
the conventions the results depend on.

- **Sprinkling.** Events are sampled uniformly with respect to the Minkowski
  volume in a causal diamond (rejection sampling on the ball-volume slice;
  equivalently uniform in null coordinates). The 1+1D and D-dimensional
  diamonds and the 1+1D forward cone are provided.
- **Causal order.** For events (t, x), i precedes j iff t_j - t_i > 0 and
  (t_j - t_i)^2 - |x_j - x_i|^2 >= -atol (null-inclusive J+; strict in time,
  so irreflexive).
- **Longest chain.** Topological-sort longest-path DP on the order; the 1+1D
  proper-time normalization uses the Brightwell-Gregory constant
  (L ~ sqrt(2 rho) tau, i.e. L ~ 2 sqrt(N) in a diamond of N events,
  consistent with Section 4.2 and Appendix A), with an
  endpoint-inclusive counting convention and acknowledged finite-size
  corrections.
- **Interval cardinality.** The open Alexandrov interval count K between two
  events estimates volume V = K/rho; in 1+1D V = tau^2/2, so tau_est =
  sqrt(2K/rho).
- **Myrheim-Meyer dimension.** The ordering fraction (related pairs over all
  pairs) is inverted against the known f(d) curve to estimate spacetime
  dimension in flat Alexandrov intervals.
- **Radar coordinates.** For a supplied observer chain, tau_minus and tau_plus
  are the latest preceding and earliest succeeding tick labels; radar time is
  their mean and radar distance their half-difference. A single chain gives
  unsigned distance; a second oriented chain supplies the sign. Lorentz maps
  between two inertial protocols are fit on their overlap [@perlick2008].
- **Rindler.** An accelerated (Rindler) observer in flat spacetime, with
  analytic two-way radar ticks; the reconstructible region is the Rindler
  wedge, with the horizon appearing as a reconstruction-inaccessibility
  boundary [@rindler1966].
- **Conformal / measure.** Positive conformal rescalings preserve causal order
  while changing volume and clock scale; volume reconstruction then requires
  supplied measure information (implemented in exp19 as local weights). In the
  tested random-thinning protocol, reconstruction is stable after density
  rescaling.

## 4. Results

Each result is a controlled validation in a known model; we report the grounded
number and its producing script.

### 4.1 R0 — order alone

**Dimension.** In the declared flat Poisson-sprinkling model, the
Myrheim-Meyer order statistic estimates spacetime dimension: for true
dimension 2 the estimate moves 1.95 -> 2.03 -> 2.01 -> 1.99 at N = 300, 600,
1200, 2400; for dimension 3, 2.96 -> 2.95 -> 3.00 -> 3.00; for dimension 4,
4.07 -> 3.97 -> 3.93 -> 3.97. Endpoint RMSE is lower at N = 2400 than at
N = 300, with non-monotonic finite-sample fluctuations for dimensions 3 and 4
(exp10; `outputs/data/dimension_reconstruction_summary.csv`).

**Uncalibrated timelike statistic (longest chain).** The longest chain follows the
Brightwell-Gregory scaling L ~ sqrt(2 rho) tau, approached from below at finite
density: at rho = 300 the normalized length L/(sqrt(rho) tau) is 1.325 with
endpoints included and 1.085 with endpoints removed, versus the asymptotic
sqrt(2) ~ 1.414, giving a finite-size low bias (mean chain proper-time error
-0.127) (exp09; `outputs/data/longest_chain_calibration_summary.csv`). The
approach toward the asymptotic value with N is shown by the chain estimator's
error falling with N in Section 4.2. Raw longest-chain length is
order-theoretic. Interpreting it as proper time, including the normalization
quoted here, belongs to R1 because it requires density and a dimension-specific
asymptotic law; absolute scale is not fixed by order alone.

### 4.2 R1 — order + measure: timelike proper time

With event density supplied, Alexandrov interval cardinality reconstructs
timelike proper time between internal pairs, tau_est = sqrt(2K/rho). The
volume-estimator relative RMSE falls from 0.273 at N = 300 to 0.105 at
N = 2400. Its separately aggregated RMSE is below the chain estimator's
(0.382 -> 0.231) at every N, but the volume statistic uses 500 pairs per N
whereas the chain statistic uses 125, so this is not a paired comparison
(exp07; `outputs/data/timelike_pair_reconstruction_summary.csv`). In a
fixed-interval sanity check, the interval-cardinality formula returns
tau = 2.0 exactly because the full sampling diamond uses K = N and rho = N/V;
this is a normalization identity, not an independent accuracy test. The chain
estimate converges toward it
(0.303 at N = 200 -> 0.0125 at N = 2000) (exp03;
`outputs/data/timelike_reconstruction_summary.csv`). The observed errors are
consistent with finite-sampling noise at the tested settings: the ratio of
reconstruction RMSE to the predicted Poisson standard deviation is 0.93-1.07
across N (binomial 1.00-1.14). Smaller bias or other estimator errors are not
excluded (exp08;
`outputs/data/probe_pair_statistical_calibration_summary.csv`).

### 4.3 R2 — order + observer: radar time and unsigned distance

Given an observer chain with clock labels, radar time and unsigned radar
distance are reconstructed from causal accessibility; at N = 300 the error
falls roughly by half per doubling of tick resolution — radar-time RMSE from
0.027 at 16 ticks to 0.0033 at 128 ticks, radar-distance RMSE from 0.072 to 0.0085, with
the accessible fraction 1.0 throughout (exp11;
`outputs/data/discrete_radar_reconstruction_summary.csv`). A single observer
determines only unsigned distance (the reflection degeneracy of Section 5).

### 4.4 R3 — + orientation: signed coordinates and the Lorentz map

A second synchronized beacon chain with known separation supplies orientation,
lifting the degeneracy to signed 1+1D coordinates. The affine map recovered
between two oriented inertial protocols approaches the Lorentz boost, with the
fitted-beta RMSE falling with tick density. Averaged over N as in Figure 2, it
falls from 0.00495 at 32 ticks to 0.000975 at 128 ticks for beta = 0.3, and
from 0.01345 to 0.00431 for beta = 0.6 (exp13;
`outputs/data/oriented_radar_lorentz_summary.csv`).

Figure 2 collects the convergence behavior of the first four rungs: dimension
recovery (a), timelike proper time (b), observer radar (c), and Lorentz-map
recovery (d).

![Figure 2](figures/fig2_convergence.png)

*Figure 2. Reconstruction accuracy in controlled models. (a) Myrheim-Meyer
dimension estimates lie near D = 2, 3, 4; endpoint RMSE is lower at N = 2400
than N = 300 despite non-monotonic finite-sample fluctuations (exp10).
(b) Timelike proper-time relative RMSE falls with N, the interval-volume
estimator below the separately aggregated longest-chain estimator (exp07).
(c) Observer radar time and distance RMSE fall about by half per doubling of
clock ticks (exp11).
(d) Recovered Lorentz beta RMSE falls with clock ticks for two boosts (exp13).
Panels (c, d) are means over N; axes are logarithmic where noted.*

### 4.5 R4 — + oriented atlas: transition-map consistency

With overlapping observer charts that retain the synchronized, calibrated
orientation reference from R3, affine Lorentz/Poincare transition maps fit
on chart overlaps show approximately consistent composition and invariant
agreement, improving with tick density: the mean transition-map beta error
falls 0.0047 -> 0.0012 and the invariant-interval RMSE 0.050 -> 0.010 (ticks
32 -> 128), while loop closure (A -> B -> C versus direct A -> C) has a
beta-composition error of 0.0072 -> 0.0017 (exp14;
`outputs/data/observer_atlas_transition_summary.csv`,
`outputs/data/observer_atlas_loop_summary.csv`). On exact analytic input the same
transition/loop machinery recovers the maps to machine precision (beta error
~5.6e-17, RMSE ~4e-17) (exp15;
`outputs/data/exact_poincare_map_sanity.csv`).

### 4.6 R5 — + conformal measure: volume and coarse-graining

Positive conformal rescalings preserve causal order while changing volume and
clock scale: under constant (1.0/1.5/2.0) and sinusoidal rescalings the causal
matrix is unchanged and the reconstructed dimension is identical (2.020), while
the proper-time ratio tracks 1.0/1.5/2.0 and the volume ratio 1.0/2.25/4.0
(exp18; `outputs/data/conformal_order_ambiguity_summary.csv`). With supplied
measure weights, weighted volume reconstruction has negligible observed bias
and decreasing relative RMSE (0.234 -> 0.118 across the tested N), whereas,
for exp19's coordinate-volume support sampling, the unweighted estimate
retains a large low bias with no material RMSE improvement
(0.569 -> 0.557, bias ~ -0.29); the analytic
volume/proper-time formulas are verified to ~1e-7 (exp19, exp20;
`outputs/data/weighted_conformal_volume_summary.csv`,
`outputs/data/conformal_volume_exact_sanity.csv`). Under random thinning, density-rescaled
reconstruction is stable (volume RMSE 0.018 -> 0.009, dimension steady ~2.0)
while the uncorrected estimate blows up to RMSE 0.270 (bias -0.183) at 25%
retention (exp23; `outputs/data/thinning_coarse_graining_summary.csv`).
Figure 3 shows the measure dependence directly.

![Figure 3](figures/fig3_measure.png)

*Figure 3. R5: volume reconstruction requires supplied measure information. For
the coordinate-volume support sampling used in exp19, the weighted estimate
has negligible observed bias and decreasing RMSE over the tested N, while
omitting the supplied local weights gives the reported persistent bias
(constant-1.5 conformal profile). Equivalent measure
information can conceptually be encoded in the sampling law rather than as
post-hoc weights; that alternative is not part of Figure 3's evidence package.
Order alone fixes no absolute scale.*

### 4.7 Horizon analogue — Rindler reconstruction-inaccessibility

For an accelerated (Rindler) observer in flat spacetime, two-way radar
reconstruction is confined to the Rindler wedge. The ideal-wedge classification
is exact — precision = recall = 1.0, zero false positives and false negatives
in every configuration — with the wedge covering about a quarter of events
(accessible fraction ~0.25) and the in-wedge radar-time RMSE falling with tick
resolution (e.g. 0.0098 -> 0.0023) (exp16;
`outputs/data/rindler_horizon_reconstruction_summary.csv`). Comparing observers
directly, all events are accessible to the inertial observer while only ~0.25
(ideal wedge) / ~0.21 (finite clock coverage) are accessible to the Rindler
observer, and no event is Rindler-accessible but inertial-inaccessible — a
strict subset (exp17;
`outputs/data/inertial_vs_rindler_accessibility.csv`). This is a controlled
flat-spacetime horizon analogue, not a black-hole simulation.

## 5. Negative results that bound the ladder

- **Conformal scale is not fixed by order alone.** Positive conformal
  rescalings leave the causal order invariant while changing physical volume
  and clock scale (Section 4.6). Absolute scale therefore requires a supplied
  measure (rung M); R0 therefore supplies uncalibrated order statistics, not a
  proper-time scale.
- **A single observer gives only unsigned distance.** One chain's radar
  distance is |x|: two targets at x = +0.1 and x = -0.1 return the identical
  single-observer distance 0.1, while the two-chain oriented protocol recovers
  the signed positions +0.1 and -0.1. Signed coordinates therefore require an
  orientation reference (rung R) (exp12;
  `outputs/data/single_observer_reflection_degeneracy.csv`).
- **Finite signal speed alone does not give Lorentzian structure.** After
  calibrating density at the final time, a regular finite-speed lattice shares
  the continuum model's leading quadratic count growth, while finite-t counts
  differ (t = 5: 21 vs 14.75; t = 30: 496 vs 496). Yet its edges lie only along the
  two lightcone diagonals (465 each), so it has a discrete symmetry, not the
  continuous Lorentz symmetry of a sprinkled causal set (exp05;
  `outputs/data/finite_speed_lattice_growth.csv`). Finite speed is necessary
  but not sufficient; statistical Lorentz compatibility is the additional
  structure.

An exploratory spacelike-distance proxy (common-past / common-future /
enclosing-interval counts) is reported as boundary-dependent and *not* a
validated spacelike estimator: the counts vary strongly with the sampling
region rather than tracking spacelike distance cleanly (exp06;
`outputs/data/spacelike_distance_proxy_summary.csv`).

## 6. Capstone validation: finite-density detectability of the conformal class

The ladder reconstructs quantities from supplied ingredients inside flat
models; Section 5 bounds it from below. This section caps it from the other
side, with the converse of the first negative result. Section 4.6 showed that
*within* a conformal class, order carries no scale: positive rescalings leave
the causal matrix unchanged. By the theorems that motivate the program, order
determines the conformal class itself — so *across* conformal classes order
must differ. Whether that difference survives finite density, for a curvature
channel that touches nothing but the order, is the question this capstone
answers affirmatively, as a preregistered experiment (P14 in the program's
internal numbering; within this program, the first finite-density crossing of
the conformal-flatness ceiling, via a 4D type N pure-Weyl construction).

### 6.1 The 1+1D conformal ceiling

Every 1+1D metric is conformally flat, so in the 1+1D stages of this program
(and its curvature-focused companions P12-P13) curvature could reach the order
only through the volume-coefficient channel: those experiments tested
curvature entering through the measure/counting channel, and the channel in
which order directly carries a conformal-class difference could not be tested
there at all. This was an intended dimensional limit of the models, not a
confounded experiment: Weyl curvature — the conformal-class content of
curvature — is identically zero below four spacetime dimensions, so a direct
test needs a 4D construction built so that the measure channel is exactly
frozen while the conformal class moves.

### 6.2 Pure-Weyl control construction

The construction is a vacuum plane wave in Brinkmann coordinates with profile
`A(u)(x^2 - y^2)` on a constant-`A` slab: Ricci-flat, Petrov type N, with
`det g = -1` identically. Because `det g = -1`, the coordinate volume form and
the Poisson sprinkling law are *identical* to flat spacetime, and the control
arm is the *same sprinkled point set* read with `A = 0`: any difference
between the arms is a functional of the relation change alone — it cannot be
a sampling artifact, an intensity mismatch, or a calibration residual. What
is *not* identical is the causal-diamond volume: the light-cone boundary
tilts, and at equal proper time the curved-arm diamond volume sits 0.4%
(`wT = 1`) to 7.3% (`wT = 2`) above flat — a deterministic design check pins
both numbers. "Same measure" here means the volume form and sprinkling law
only, never interval volumes; conflating the two was the first design draft's
central error, caught in review, and the distinction is load-bearing.

### 6.3 Finite-density detection design

An exploratory probe chain (four probes and a confirmation stage, run on
seed streams disjoint from everything that follows) selected the operating
point — slab `(1.0, 1.0, 2.0, 6.0)`, `w = 1.0`, expected 300 events per
sprinkling — and the primary statistic, the global relation fraction (the one
candidate computable from a single poset alone). The preregistered stage then
froze two claims with deliberately different boundaries, never merged into
one sentence:

- **C1 (paired ensemble mean).** The estimand is
  `theta_Delta = E[f_A - f_0]`, the ensemble mean of the *net* change of the
  relation fraction over paired readings of the same points — not any single
  point set's displacement, and not the gained-plus-lost total motion.
  Confirmation requires the 95% Student-t interval of the mean paired
  difference (n = 3000 sprinklings) to sit entirely above a frozen margin
  `epsilon_Delta = 3.579e-4`, an operational margin anchored to the probe
  chain's confirmed block.
- **C2 (single-poset classifier, independent replication).** A frozen
  classifier reading only the single-poset relation fraction must separate
  curved from flat ensembles (n = 4800 per unpaired arm) under three frozen
  interval rules (standardized separation, AUC with an exact boundary bound
  at complete separation, balanced accuracy), each with its own equivalence
  band. This claim replicates the probe chain's confirmation stage on fresh
  seeds; its failure would have been reported as a conflicting replication,
  never a retroactive cancellation.

The stage verdict POSITIVE requires both claims confirmed, certified as a
joint event before execution (4000/4000 joint-effect replicates; the null
branches certify each claim's equivalence power separately at an exact
Clopper-Pearson 95% lower bound of at least 0.90; a joint equivalence verdict
is deliberately not defined, because its joint null rate sits below that
floor at the frozen sizes).

### 6.4 Preregistered result

The campaign ran once, from a clean checkout of the freeze commit, on frozen
seeds, and recorded the executing commit in the artifact. Both claims
confirmed; the stage is POSITIVE.

| Claim | Verdict | Result (95% CI) |
| --- | --- | --- |
| C1 paired ensemble mean, n = 3000 | confirmed | mean 0.0502929 [0.0501046, 0.0504812]; lower end 140x the margin 3.579e-4 |
| C2 classifier replication, n = 4800/arm | confirmed | separation s = 11.199 [11.035, 11.362]; AUC = 1.0 [0.999232, 1.0]; balanced accuracy = 1.0 [0.986, 1.014] |

The two arms' relation-fraction samples are completely separated (the minimum
curved-arm value exceeds the maximum flat-arm value in 4800 draws per arm);
the AUC interval at complete separation is the frozen exact bound, not the
degenerate Wald interval. Every relation census in both claims recorded zero
ambiguous and zero escalated pairs. Frozen sentences are carried verbatim:
C1 "the paired ensemble mean shift exceeds epsilon_Delta"; C2 "the probe
chain's separation is independently reproduced"
(`docs/prereg/p14_prereg_results.json`).

### 6.5 What this establishes

At a fixed box, density, and profile, an order-only statistic separates flat
from curved ensembles: the conformal-class information that order carries by
theorem is *detectable at finite density*, in the one construction where the
sprinkling measure is exactly frozen. Together with Section 4.6 this closes
the conformal story in both directions — within a class, order is blind to
scale; across classes, at least at this operating point, order visibly moves.

### 6.6 What this does not establish

It is not a general-Weyl discriminator (the plane wave is Petrov type N and
its exact volume form is unrepresentative); it is not Weyl-tensor recovery
(detection and recovery are different claims; no reconstructed quantity is
produced, which is why this is a capstone boundary experiment and not a new
ladder rung); it does not establish box- or density-independence (the licensed
sentence is bound to the frozen operating point); and it does not open the
Schwarzschild generalization. A separate cost measurement (S1) prices only the
causal-predicate component of that path — about 0.77 ms per pair on the tested
solver, patch, and tolerance, roughly 360x the plane-wave predicate — and even
an affordable predicate leaves the Schwarzschild diamond-volume oracle
unresolved; an unaffordable one would hold only that solver-domain-budget
path (`docs/prereg/p14_s1_cost.json`).

The full execution provenance — the preregistration's freeze ordering, its
mechanical gates, and the commit-ancestry contract — is in Appendix B.

## 7. Discussion

The ladder is a branched accounting device, not a single cumulative chain.
Across its dependency rows it recovers a substantial fraction of operational
Lorentzian geometry — dimension, timelike duration, radar decomposition,
Lorentz and Poincare consistency, volume — from a causal order plus short,
explicit ingredient sets. It is equally a list of what each reconstruction
*costs*: absolute scale costs a
measure; a signed spatial coordinate costs an orientation; a metric-not-merely-
conformal statement costs a volume element; Lorentzian structure is not bought
by finite speed alone. Stating both directions together is the point: it turns
"causal structure fixes the conformal class, and a volume element fixes the
remaining conformal factor" from a slogan into a per-quantity operational
ledger with measured error behavior.

This framing also clarifies what the program does *not* show and sets up the
open question. Everything here is reconstruction inside known models: the
geometry is put in (by sprinkling from Minkowski or by a supplied protocol) and
recovered. It does not address whether geometry could *emerge* from an order
that was not built from a geometry — the representability question — which
requires a validated discriminator and is taken up separately. The Section 6
capstone sharpens the boundary from the other side: the one piece of geometry
that order itself carries — the conformal class — is not merely present by
theorem but detectable at finite density, exactly where the measure channel is
frozen and only the light cones move.

## 8. Claim boundary

We claim, as controlled validations in known 1+1D (and higher-D for dimension)
models: dimension is estimated from order statistics in flat sprinklings;
timelike proper time is
recoverable from interval cardinality once a density is supplied, with
finite-sampling-consistent error; radar time and unsigned distance are
recoverable from an observer protocol, signed coordinates and the Lorentz map
with an orientation reference, and atlas transition maps with overlapping
oriented, calibrated charts; volume is recoverable with supplied global
density and local measure information and is stable
under density-rescaled coarse-graining; and a Rindler wedge is the
reconstructible region for an accelerated observer. As the preregistered
capstone (Section 6), at the frozen 3+1D operating point: the paired ensemble
mean of the relation-fraction change under a pure-Weyl deformation exceeds its
frozen margin (C1), and a frozen single-poset classifier separates curved from
flat ensembles, independently replicating the probe chain's confirmation (C2).

We do not claim: that spacetime is reducible to, or emerges from, causal order;
that causal order alone yields absolute scale, the conformal factor, signed
coordinates, or a unique atlas; that finite signal speed implies relativity; or
that any of these finite validations establish a physical theory. Reconstructing
a geometry that was put into the model is not deriving geometry from order.
For the capstone we additionally do not claim: a general-Weyl discriminator,
Weyl-tensor recovery, box- or density-independence of the separation, or any
verdict on the Schwarzschild path (S1 supplies a price for one solver path,
not a verdict; the diamond-volume oracle remains open). Priority is claimed
only within this program; no literature-wide first is asserted before a
priority search.

## 9. Limitations and future work

Results are controlled and mostly 1+1D (dimension is checked to 4D, and the
Section 6 capstone is 3+1D). The
constants are convention-dependent (chain endpoint convention, null-inclusive
causal relation, Myrheim-Meyer normalization); we state each where it is used.
The spacelike proxy is exploratory. The capstone is bounded by its
construction: Petrov type N with an exactly flat volume form — the property
that makes the experiment clean also makes it unrepresentative of generic
`Weyl != 0` spacetimes — and by its single frozen operating point;
generalizing beyond type N sits behind the unresolved Schwarzschild
diamond-volume oracle, with only the causal-predicate component priced so
far. The natural next question — whether
observer-relative distance *order* can be validated as recovering latent
geometry, as opposed to being reconstructed from a supplied one — is the
subject of a companion study that builds a preregistered discriminator on this
foundation, measures its dose-response to geometry dilution and its dimension
selection in 2+1D, and then carries it to orders produced by growth dynamics
and by an action-weighted 2D-order ensemble. At N = 600, random-start
post-burn-in configurations at beta = 2 and beta = 8 pass the frozen
instrument, while the beta = 32 bipartite-start control blocks structurally.
These are not certified equilibrium draws, and the companion study did not
establish an equilibrium transition or finite-size scaling.

## 10. Reproducibility

Foundation-layer baseline commit `325df55`. Every number in Sections 4-5 has a
cited `experiments/exp*.py` producer and expected summary path. Five summaries
used by the figures are committed under `figures/data/`; the other named
`outputs/data/` tables are gitignored generated outputs and are not present in
a clean checkout. `artifact_manifest.md` records this boundary and the hashes
of the committed summaries. The submission gate is the 19 legacy tables plus
the finalized Section 6 evidence bundle, provenance-locked. Conventions
(sprinkling measure, causal relation, chain and interval normalizations,
Myrheim-Meyer inversion) are fixed in the foundation modules and stated in
Section 3.

Section 6 carries no regeneration debt: every number it cites lives in a
committed artifact, and the repository's test suite recomputes the capstone's
metrics, verdicts, and sentences from the stored raw samples, reproduces a
prefix of every frozen seed stream, and asserts the commit-ancestry contract
of Appendix B. The stage runner (`experiments/positive_control/p14_prereg.py`)
exposes three modes — `preflight`, `manifest`, `campaign` — and the campaign
mode refuses to run unless the freeze manifest exists, its recorded preflight
digest matches the committed artifact, and the certified source digests match
the checked-out files.

## Appendix A: conventions and normalizations

| Quantity | Convention used |
| --- | --- |
| Causal relation | null-inclusive J+: i precedes j iff t_j - t_i > 0 and (dt)^2 - (dx)^2 >= -atol, atol = 1e-12; strict in time (irreflexive) |
| Longest chain | endpoint-inclusive count; 1+1D normalization L ~ sqrt(2 rho) tau (asymptotic normalized value sqrt(2) ~ 1.414), finite-size corrected |
| Interval volume (1+1D) | open Alexandrov count K; V = K/rho; V = tau^2/2, so tau_est = sqrt(2K/rho) |
| Myrheim-Meyer | ordering fraction r = (related ordered pairs)/(N(N-1)) inverted against f(d) = Gamma(d+1)Gamma(d/2)/(4 Gamma(3d/2)); equals the standard unordered form up to a factor 2 on both sides |
| Radar coordinates | tau_minus = latest preceding tick, tau_plus = earliest succeeding tick; radar time = (tau_plus + tau_minus)/2, radar distance = (tau_plus - tau_minus)/2 |
| Density (finite-N) | fixed-N sprinkle with empirical density; Poisson and fixed-N binomial sampling models both reported |

These are stated so the reported constants (e.g. the sqrt(2) chain
normalization, the interval-volume factor) are unambiguous; each is fixed in
the foundation modules.

## Appendix B: capstone execution provenance

The Section 6 result is reproducible not merely in the sense that its code is
committed, but in the sense that the *order of operations* — freeze before
execution, execution before recording — is provable from the repository's
commit graph and enforced by mechanical gates.

**Commit chain.** Design and probe chain: PR #38-#48. Preregistration:
preparation commit (PR #49), S1 cost measurement (PR #50), then
`P' = 51875a2` (a one-source seed-ledger fix) -> preflight executed on a clean
checkout of `P'`, its artifact recording `code_version = P'` and the SHA-256
digests of every execution-relevant source file -> final-freeze commit
`F = b858b08` carrying the certification artifact and the freeze manifest
(preflight digest, `preflight_code_version = P'`, the S1 price sentence, the
frozen sizes and seeds; PR #51) -> the campaign executed on a checkout of
exactly `F` (not a descendant merge commit, so the recorded executing HEAD is
`F` itself) -> results commit `R = 5126ddf`, a direct child of `F` (PR #52).
A test asserts `P' ≺ F ≺ R` with `git merge-base --is-ancestor`.

**Gates.** Stochastic modes refuse a dirty worktree at entry. The campaign
mode refuses to run unless (i) the freeze manifest exists, (ii) its recorded
preflight digest matches the committed preflight artifact, and (iii) the
source digests certified at preflight entry match the checked-out files —
captured at entry and re-checked at exit, so a mid-run edit aborts the run
with nothing recorded. Interrupted runs restart the same seed; an ambiguity
violation blocks the stage rather than permitting a seed swap.

**Artifacts.** The complete Section 6 evidence bundle, with digests, is
enumerated in `artifact_manifest.md`; the per-sentence mapping from manuscript
text to evidence file is recorded there as the citation-to-artifact inventory.

## References

Verified bibliography in `citations/references.bib` (shared verified core with
Paper B), including the radar-method and uniformly accelerated-frame sources
used here.
