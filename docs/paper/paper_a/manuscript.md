# An operational reconstruction ladder for spacetime quantities from causal order

**Draft v0.1.** Causal Spacetime Lab. Quantitative results are grounded in the
experiment output CSVs (regenerable; each cited with its producing script);
cells awaiting a number are marked `[NUM: exp##]`.

## Abstract

We study, in controlled 1+1D (and, for dimension, higher-D) Minkowski models,
which spacetime quantities can be operationally reconstructed from a causal
(accessibility) order once a minimal, explicitly declared set of additional
ingredients is supplied. We organize the reconstructions as a ladder. From
causal order alone, order statistics recover the spacetime dimension and the
longest antichain-free chain fixes the *shape* of timelike separation, but not
its absolute scale. Adding a counting measure (event density) turns Alexandrov
interval cardinality into a timelike proper-time estimate whose error is
consistent with finite-sampling noise. Adding an observer chain with clock
labels yields radar time and unsigned radar distance; adding an orientation
reference lifts the reflection degeneracy to signed coordinates and lets one
recover the Lorentz map between two inertial protocols; adding overlapping
charts yields an observer atlas with approximately consistent Poincare
transition maps; adding conformal measure weights makes volume reconstruction
possible under the conformal ambiguity, and coarse-graining is stable only when
density is rescaled. We also give a Rindler horizon analogue, in which an
accelerated observer's two-way radar reconstruction is confined to the expected
wedge, and a finite-speed lattice counterexample showing that finite signal
speed alone does not produce Lorentzian structure. The contribution is not any
single reconstruction — most are standard — but the explicit accounting of what
each rung requires and the negative results that bound it. We make no claim
that spacetime is reducible to causal order; the reconstructions are controlled
validations inside known models.

## 1. Introduction

A conservative reading of the causal-set and order-first programs is that
causal order carries much, but not all, of Lorentzian geometry: with a volume
element it fixes geometry up to the well-understood conformal freedom
[@blms1987; @malament1977; @hkm1976; @kronheimer1967; @sorkin2005; @surya2019].
This paper treats that statement operationally and quantitatively. Rather than
ask whether spacetime "is" a causal set, we ask: fix a causal (or
null-accessibility) order in a known model; then, supplying one additional
ingredient at a time, which spacetime quantities can a finite procedure
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
quantities by the minimal supplied ingredient, with per-rung error behavior in
controlled models; (2) grounded validations of dimension, proper time, radar
decomposition, Lorentz-map and atlas consistency, and measure-dependent volume;
(3) three bounding negative results (conformal scale, reflection degeneracy,
finite-speed lattice); and (4) a Rindler reconstruction-inaccessibility
analogue. All results are controlled validations, not evidence that geometry
reduces to order.

## 2. The reconstruction ladder

We separate the *primitive* from *supplied ingredients* and index rungs by the
minimal ingredient each requires.

Primitive:

- **Causal / accessibility order.** A strict partial order on events (in the
  continuum models, Minkowski causal precedence; null-inclusive).

Supplied ingredients, in increasing strength:

- **M — counting measure / event density.** A volume element, supplied as a
  density or a sampling process.
- **O — observer chain with clock labels.** A timelike chain carrying ordered
  tick labels (an operational clock), defining a radar protocol.
- **R — orientation reference.** A second synchronized chain with known
  separation, fixing a spatial side.
- **A — observer atlas.** Multiple overlapping observer charts.
- **W — conformal measure weights.** A supplied local volume weighting.

Rungs (minimal ingredient -> what is reconstructed):

| Rung | Ingredient | Reconstructed | Bounded by |
| --- | --- | --- | --- |
| R0 | order only | dimension; timelike *shape* (longest chain); conformal/causal structure | no absolute scale |
| R1 | order + M | timelike proper time (interval cardinality -> volume -> tau) | needs a density |
| R2 | order + O | radar time; unsigned radar distance | sign undetermined |
| R3 | order + O + R | signed coordinates; Lorentz map between protocols | supplied orientation |
| R4 | order + O + A | atlas transition maps (Poincare); invariant agreement | supplied charts |
| R5 | order + M + W | volume under conformal ambiguity; coarse-graining stability | conformal factor supplied |

Three negative results bound the ladder from below (Section 5): conformal scale
is not fixed by order alone (motivating M and W); a single observer gives only
unsigned distance (motivating R); and finite signal speed alone does not give
Lorentzian structure.

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
  (L ~ 2 sqrt(rho) tau, i.e. L ~ 2 sqrt(N/2) in a diamond), with an
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
  between two inertial protocols are fit on their overlap.
- **Rindler.** An accelerated (Rindler) observer in flat spacetime, with
  analytic two-way radar ticks; the reconstructible region is the Rindler
  wedge, with the horizon appearing as a reconstruction-inaccessibility
  boundary.
- **Conformal / measure.** Positive conformal rescalings preserve causal order
  while changing volume and clock scale; volume reconstruction then requires
  supplied measure weights, and coarse-graining (random thinning) is stable
  only with density rescaling.

## 4. Results

Each result is a controlled validation in a known model; we report the grounded
number and its producing script.

### 4.1 R0 — order alone

**Dimension.** The Myrheim-Meyer estimator recovers spacetime dimension in flat
Alexandrov intervals and converges toward the true value as N grows: for true
dimension 2 the estimate moves 1.95 -> 2.03 -> 2.01 -> 1.99 at N = 300, 600,
1200, 2400; for dimension 3, 2.96 -> 2.95 -> 3.00 -> 3.00; for dimension 4,
4.07 -> 3.97 -> 3.93 -> 3.97, with RMSE decreasing with N
(exp10; `outputs/data/dimension_reconstruction_summary.csv`).

**Timelike shape (longest chain).** The longest chain reproduces the
Brightwell-Gregory scaling L ~ 2 sqrt(rho) tau, with L/sqrt(N) approaching the
constant from below as N grows: `[NUM: exp09 L/sqrt(N) trend]`
(exp09; `outputs/data/longest_chain_calibration_summary.csv`). This fixes the
*shape* of timelike separation from order plus a density normalization, with
O(1) endpoint and finite-size corrections; absolute scale is not fixed by order
alone.

### 4.2 R1 — order + measure: timelike proper time

With event density supplied, Alexandrov interval cardinality reconstructs
timelike proper time between internal pairs, tau_est = sqrt(2K/rho). The
reconstruction error decreases with interval size / N: `[NUM: exp07 RMSE or
relative error vs N]` (exp07;
`outputs/data/timelike_pair_reconstruction_summary.csv`). The full-diamond
sanity check agrees: `[NUM: exp03 error]`
(exp03; `outputs/data/timelike_reconstruction_summary.csv`). The observed
errors are consistent with finite-sampling (Poisson / fixed-N binomial)
expectations rather than bias or estimator error: `[NUM: exp08 consistency]`
(exp08; `outputs/data/probe_pair_statistical_calibration_summary.csv`).

### 4.3 R2 — order + observer: radar time and unsigned distance

Given an observer chain with clock labels, radar time and unsigned radar
distance are reconstructed from causal accessibility; the error decreases as
the observer tick density increases: `[NUM: exp11 error vs ticks]`
(exp11; `outputs/data/discrete_radar_reconstruction_summary.csv`). A single
observer determines only unsigned distance (the reflection degeneracy of
Section 5).

### 4.4 R3 — + orientation: signed coordinates and the Lorentz map

A second synchronized beacon chain with known separation supplies orientation,
lifting the degeneracy to signed 1+1D coordinates. The affine map recovered
between two oriented inertial protocols approaches the Lorentz boost, with the
recovered rapidity/beta residual decreasing with tick density: `[NUM: exp13
beta residual vs ticks]`
(exp13; `outputs/data/oriented_radar_lorentz_summary.csv`).

### 4.5 R4 — + atlas: transition-map consistency

With overlapping observer charts, affine Lorentz/Poincare transition maps fit
on chart overlaps show approximately consistent composition and invariant
agreement, improving with tick density: `[NUM: exp14 transition/loop
consistency]` (exp14; `outputs/data/observer_atlas_transition_summary.csv`,
`observer_atlas_loop_summary.csv`); an exact Poincare-map sanity check confirms
the machinery `[NUM: exp15]`.

### 4.6 R5 — + conformal measure: volume and coarse-graining

Positive conformal rescalings preserve causal order while changing volume and
clock scale (exp18; `outputs/data/conformal_order_ambiguity_summary.csv`):
`[NUM: exp18 order preserved / scale change]`. With supplied measure weights,
weighted volume is reconstructed (exp19-20;
`outputs/data/weighted_conformal_volume_summary.csv`,
`conformal_volume_exact_sanity.csv`): `[NUM: exp19/20 volume error]`. Under
random thinning, corrected (density-rescaled) reconstruction is stable while
uncorrected density fails (exp23;
`outputs/data/thinning_coarse_graining_summary.csv`): `[NUM: exp23 corrected vs
uncorrected]`.

### 4.7 Horizon analogue — Rindler reconstruction-inaccessibility

For an accelerated (Rindler) observer in flat spacetime, two-way radar
reconstruction is confined to the Rindler wedge; events beyond the horizon are
reconstruction-inaccessible, and finite observer-chain coverage produces the
expected false-positive/negative behavior near the boundary: `[NUM: exp16/exp17
accessibility]` (exp16, exp17;
`outputs/data/rindler_horizon_reconstruction_summary.csv`,
`inertial_vs_rindler_accessibility.csv`). This is a controlled flat-spacetime
horizon analogue, not a black-hole simulation.

## 5. Negative results that bound the ladder

- **Conformal scale is not fixed by order alone.** Positive conformal
  rescalings leave the causal order invariant while changing physical volume
  and clock scale (Section 4.6). Absolute scale therefore requires a supplied
  measure (rung M); this is why R0 recovers only timelike *shape*.
- **A single observer gives only unsigned distance.** Radar distance from one
  chain is |x|; the side is undetermined. Signed coordinates require an
  orientation reference (rung R): `[NUM: exp12 reflection degeneracy]`
  (exp12; `outputs/data/single_observer_reflection_degeneracy.csv`).
- **Finite signal speed alone does not give Lorentzian structure.** A
  finite-speed lattice with a bounded propagation speed produces
  causal-cone-like growth without Lorentz invariance: `[NUM: exp05 growth]`
  (exp05; `outputs/data/finite_speed_lattice_growth.csv`). Finite speed is
  necessary but not sufficient; statistical Lorentz compatibility (as in
  sprinkled causal sets) is the additional structure.

An exploratory spacelike-distance proxy (common-past / common-future /
enclosing-interval counts) is reported as boundary-dependent and *not* a
validated spacelike estimator (exp06;
`outputs/data/spacelike_distance_proxy_summary.csv`): `[NUM: exp06 note]`.

## 6. Discussion

The ladder is an accounting device. Read top-down it recovers a substantial
fraction of operational Lorentzian geometry — dimension, timelike duration,
radar decomposition, Lorentz and Poincare consistency, volume — from a causal
order plus a short, explicit list of supplied ingredients. Read bottom-up it is
equally a list of what each reconstruction *costs*: absolute scale costs a
measure; a signed spatial coordinate costs an orientation; a metric-not-merely-
conformal statement costs a volume element; Lorentzian structure is not bought
by finite speed alone. Stating both directions together is the point: it turns
"causal order plus a volume element fixes geometry up to conformal factor" from
a slogan into a per-quantity operational ledger with measured error behavior.

This framing also clarifies what the program does *not* show and sets up the
open question. Everything here is reconstruction inside known models: the
geometry is put in (by sprinkling from Minkowski or by a supplied protocol) and
recovered. It does not address whether geometry could *emerge* from an order
that was not built from a geometry — the representability question — which
requires a validated discriminator and is taken up separately.

## 7. Claim boundary

We claim, as controlled validations in known 1+1D (and higher-D for dimension)
models: dimension is recoverable from order statistics; timelike proper time is
recoverable from interval cardinality once a density is supplied, with
finite-sampling-consistent error; radar time and unsigned distance are
recoverable from an observer protocol, signed coordinates and the Lorentz map
with an orientation reference, and atlas transition maps with overlapping
charts; volume is recoverable with supplied conformal weights and is stable
under density-rescaled coarse-graining; and a Rindler wedge is the
reconstructible region for an accelerated observer.

We do not claim: that spacetime is reducible to, or emerges from, causal order;
that causal order alone yields absolute scale, the conformal factor, signed
coordinates, or a unique atlas; that finite signal speed implies relativity; or
that any of these finite validations establish a physical theory. Reconstructing
a geometry that was put into the model is not deriving geometry from order.

## 8. Limitations and future work

Results are controlled and mostly 1+1D (dimension is checked to 4D). The
constants are convention-dependent (chain endpoint convention, null-inclusive
causal relation, Myrheim-Meyer normalization); we state each where it is used.
The spacelike proxy is exploratory. The natural next question — whether
observer-relative distance *order* can be validated as recovering latent
geometry, as opposed to being reconstructed from a supplied one — is the
subject of a companion study that builds a preregistered discriminator on this
foundation.

## 9. Reproducibility

Foundation-layer baseline commit `325df55`. Every number in Section 4/6 is
produced by the cited `experiments/exp*.py` script (run with `PYTHONPATH=src`)
and read from the named summary CSV under `outputs/data/`. Conventions
(sprinkling measure, causal relation, chain and interval normalizations,
Myrheim-Meyer inversion) are fixed in the foundation modules and stated in
Section 3.

## Appendix A: conventions and normalizations

`[to be finalized: table of the exact constants — Brightwell-Gregory chain
constant and endpoint convention; eta_D interval-volume constants; Myrheim-Meyer
f(d) form and inversion; radar tau+/tau- definition; null-inclusive atol.]`

## References

Verified bibliography in `citations/references.bib` (shared verified core with
Paper B). Additional radar-coordinate and Rindler references to be added in the
bibliography pass.
