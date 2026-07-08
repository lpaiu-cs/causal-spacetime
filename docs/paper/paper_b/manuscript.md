# A validated discriminator for latent geometry in discrete causal order, and its response to geometry dilution

**Draft v0.1.** Causal Spacetime Lab. All quantitative results trace to frozen
preregistration artifacts (`docs/prereg/frozen/`); see Section 9.

## Abstract

Causal-set and order-first programs ask when metric geometry can be treated as
an effective representation of a discrete causal order rather than as primitive
structure. A recurring difficulty is methodological: a procedure that fits
low-dimensional coordinates to order-derived comparisons can appear to
"recover geometry" even when none is present, because embeddability is not the
same as recovering true spatial structure. We address this directly. We build
a preregistered pipeline that measures order-level bracket-width response
profiles for a set of observer reference chains on a causal set, converts them
to a reference-centered (parallax) dissimilarity, and tests whether the
implied distance order admits a stable low-dimensional ordinal embedding, under
fixed acceptance gates. We first validate the pipeline as a *discriminator*
(experiment PC-V1): on 1+1D Minkowski-sprinkled causal sets it passes every
gate on 20/20 held-out confirmatory seeds and recovers the true spatial order
of targets, while on density-matched geometry-free order and on
consistency-destroyed controls it blocks on every seed (10/10 for each of two
control families). Establishing the discriminator required catching, before
freezing, two confounds that would otherwise have produced false positives: a
shared per-observer scalar (removed by parallax centering) and a degenerate
near-complete control order (fixed by matching post-closure density). We then
make geometry the manipulated variable (experiment P1): rewiring a fraction
epsilon of the geometric covering relation toward random time-respecting edges
at fixed relation density. Geometry recovery degrades monotonically in epsilon
(Spearman rho >= 0.85 on 20/20 confirmatory seeds), with an identifiable
graded transition near epsilon ~ 0.31 and a robust "false-pass" window in which
profiles still embed in one dimension while no longer encoding true space. We
conclude that, on this controlled family at finite scale, effective metric
representability responds to the *amount* of geometric consistency in a causal
order, not to its density; and that embeddability alone is an unsafe indicator
of geometry, which has direct methodological implications for order-first
reconstruction claims. We make no claim of spacetime emergence.

## 1. Introduction

The hypothesis motivating this work is that operational time and distance can
be reconstructed from primitive causal or null information-accessibility
relations, with metric scale requiring additional structure such as event
density or observer protocols, rather than from velocity as displacement over
time. Within causal-set theory this is a familiar stance: causal order plus a
counting measure fixes a great deal of Lorentzian geometry in suitable
continuum limits [@blms1987; @sorkin2005; @surya2019]. A weaker, finite, and more
operational question is whether *observer-relative distance order*, derived
from causal accessibility and an observer protocol, admits a stable
low-dimensional metric representation, and under what conditions.

This question is easy to ask and easy to answer badly. Any procedure that fits
coordinates to minimize violations of order constraints will return *some*
low-dimensional configuration. Low fit error ("embeddability") does not
establish that the fitted coordinates track true spatial structure: a
dominant, geometry-free common mode in the input can be perfectly embeddable
while carrying no spatial information. A program that reads embeddability as
evidence of geometry can therefore accumulate results that are internally
consistent but uninterpretable, because it has never established that its
instrument distinguishes geometry from its absence.

We take the discriminator seriously as the primary object. Before asking
whether any causal order "has geometry," we ask whether our measurement-and-fit
pipeline can be shown to (i) pass on data with known latent geometry and (ii)
block on matched data without it. Only a pipeline that does both is a valid
instrument, and only then are its verdicts on ambiguous inputs meaningful. This
is the content of experiment PC-V1 (Section 4). We then use the validated
instrument to measure a dose-response: how recovery changes as geometric
consistency is titrated away at fixed density (experiment P1, Section 5).

Our contributions are: (1) a preregistered, provenance-tracked pipeline from
causal order to a representability verdict, with an instrument-integrity
protocol that fixes the fit budget and forbids a saturating stability metric;
(2) a positive/negative control (PC-V1) demonstrating sensitivity and
specificity, including two confounds that had to be removed before freezing;
(3) a geometry-dilution dose-response (P1) showing monotone degradation, an
identifiable transition, and a false-pass window; and (4) the methodological
consequence that embeddability alone over-reports geometry, so a
truth-recovery check is load-bearing.

## 2. Related work and positioning

Causal-set theory reconstructs timelike distance from longest chains and
volume from interval cardinality [@blms1987; @brightwell1991], and dimension
from ordering fractions [@myrheim1978; @meyer1988]. Order-theoretic results
fix conformal structure from causal order in the continuum
[@malament1977; @hkm1976; @kronheimer1967]. Ordinal-embedding methods from
machine learning fit low-dimensional coordinates to quadruplet ("is d(i,j) <
d(k,l)") comparisons [@shepard1962; @kruskal1964; @kleindessner2014;
@agarwal2007]. Our pipeline sits at the intersection: it derives quadruplet
comparisons operationally from causal accessibility and observer protocols,
then applies ordinal embedding, but treats the fitted coordinates as a
representation diagnostic rather than as physical coordinates. The novel
methodological emphasis is the explicit discriminator validation and the
geometry-dilution dose-response, neither of which, to our knowledge, is
standard in order-first reconstruction studies.

Inline citations use pandoc-style keys (`@blms1987`, ...) matching the
verified `references.bib`; see `citations/citation_verification_report.md` for
the source confirming each entry.

## 3. Methods

All quantities are in natural units with c = 1. The pipeline reuses a
foundation layer (sprinkling, causal order, discrete radar, ordinal embedding)
whose numerical correctness against known causal-set results was verified
independently; this paper's contributions are the profile measurement, the
discriminator protocol, and the dilution experiment built on that layer.

### 3.1 Scenes: causal sets with observer chains

A scene is an N-event sprinkle in a 1+1D Minkowski causal diamond
(N = 900, diamond duration T = 2.0) plus K = 6 supplied stationary observer
"reference chains" at fixed spatial positions x0 in {-0.25, ..., 0.25}, each
carrying 96 tick events. The causal order is the standard 1+1D Minkowski
precedence relation (null-inclusive). Targets are bulk events in a thin band
(|t| <= 0.10, |x| <= 0.25) that are two-sided bracketed by all six chains, so
that missing-data effects are separated from representability; a scene must
yield >= 30 such targets (subsampled to at most 40). Every scene records
content hashes of its event array and causal matrix for provenance.

### 3.2 Bracket-width echo profiles and parallax dissimilarity

For target j and reference chain r, using tick indices as clocks, let p be the
last tick index that causally precedes j and s the first tick index that j
precedes on chain r. The bracket width W[j, r] = s - p is an order-level
"echo" observable; in the continuum it scales as 2|x_j - x0_r|/dtau. Clock
values are never used, only tick order.

Bracket widths carry a shared per-target component (a global/temporal common
mode present in any time-respecting order) that is not spatial. We therefore
center each target's widths across the chains that reach it,
P[j, r] = W[j, r] - mean_r W[j, r], keeping only the cross-observer *parallax*.
This operationalizes a simple principle: a single scalar shared across
observers is not a distance structure; only cross-observer disagreement is.
The profile dissimilarity is the RMS parallax difference over common reachable
chains, D[i, j] = RMS_r (P[i, r] - P[j, r]).

### 3.3 Representability gates

From D we sample quadruplet constraints "D(i,j) + delta < D(k,l)" with the
margin delta derived from the measured dissimilarity resolution (25th
percentile of positive gaps), split at the *pair* level (no target pair
contributes to both train and held-out, avoiding leakage), and fit a
low-dimensional ordinal embedding by hinge-loss gradient descent at dimensions
{1, 2, 3}. A scene is scored at dimension d = 1 (the spatial dimension of the
model) by four quantities: held-out constraint violation; the separation
between held-out violation on the real profiles and on a destructive,
marginal-preserving column-shuffle null ("null gap"); restart stability
measured as sign-discordance of pair distances across independent restarts;
and, where the true coordinates are known, the truth-order error (sign
discordance between fitted and true pair distances).

### 3.4 Instrument-integrity protocol

Two failure modes are precluded by construction. First, the fit budget is
fixed (>= 1500 steps, 5 restarts) and the runner asserts the executed budget
equals the requested one; an under-converged structured fit is otherwise
indistinguishable from a null fit. Second, a positional argsort-mismatch
stability metric (which saturates near one for dense continuous distances under
any restart jitter) is forbidden; sign discordance is used instead. Both
choices follow from a prior program in which a silently reduced budget and a
saturating metric made every verdict uninterpretable.

### 3.5 Preregistration

Each experiment is preregistered with fixed hypotheses, an exploratory
calibration stage that proposes thresholds by a mechanical rule, and a
confirmatory stage on disjoint fresh seeds under frozen thresholds. Instrument
repairs are permitted only before the freeze and are logged. Thresholds are
never retuned after the freeze; under-covered or scene-invalid cells are
recorded, never silently dropped. Frozen thresholds and decision registries
live in version control (`docs/prereg/frozen/`), not in the (gitignored)
outputs tree, so no decision input is ephemeral.

## 4. Experiment PC-V1: the pipeline is a discriminator

**Design.** Calibration on seeds 0-9 (exploratory) proposes gate thresholds;
confirmatory sensitivity on seeds 100-119 and specificity on seeds 200-209
apply the frozen thresholds. Sensitivity uses geometric scenes (H-SENS: they
should pass and recover x). Specificity uses two geometry-free control
families (H-SPEC: they should block): a density-matched random time-respecting
order, and a column-shuffle of the real profiles that preserves marginals but
destroys cross-observer consistency.

**Two confounds fixed before freezing.** The first calibration passed, but its
specificity controls exposed two instrument defects. (i) A geometry-free order,
once made density-comparable, still passed the embeddability gate: its profiles
were dominated by a shared per-target scalar (mean cross-chain column
correlation ~0.7 versus ~0.0 for geometric scenes) that embeds in one
dimension without encoding space (fitted 1D coordinate correlated with true x
at ~0.1-0.5 versus ~0.98 for geometric scenes). Centering to parallax profiles
(Section 3.2) removes it. (ii) The random control matched pre-closure edge
density but percolated under transitive closure to a near-complete order
(density ~0.997), an unfair degenerate foil; we instead match the geometric
post-closure density by bisection. Both are logged as pre-freeze repairs; the
frozen thresholds derive only from the post-repair calibration.

**Frozen thresholds** (gate dimension d = 1): held-out violation <= 0.05, null
gap >= 0.15, restart stability <= 0.15, truth-order error <= 0.15. These were
set mechanically from the calibration distributions; both hard floors passed
(median structured held-out 0.025 <= 0.35; effect size between structured and
the destructive null, Cohen's d = 7.2).

**Results.** On calibration (d = 1), the two arms are cleanly separated with no
overlap:

| condition | role | held-out (min-max) | mean |
| --- | --- | --- | --- |
| structured (geometric) | should pass | 0.014-0.032 | 0.023 |
| relabel-symmetry (control) | should mirror structured | 0.012-0.047 | 0.023 |
| column-shuffled | should block | 0.210-0.376 | 0.277 |
| density-matched random | should block | 0.202-0.348 | 0.261 |

*Table 1. PC-V1 calibration: held-out violation by condition at the gate
dimension d = 1 (seeds 0-9). The two arms do not overlap.*

Structured scenes recover the true spatial order (truth-order error mean
0.097); the symmetry control mirrors structured (mean 0.023), confirming
label-equivariance. On the confirmatory runs under frozen thresholds:

- **H-SENS supported.** 20/20 seeds pass all four gates (rule >= 16/20); every
  individual gate passed 20/20, including the tightest (stability).
- **H-SPEC supported.** Both control families block on every seed (rule
  >= 8/10): density-matched random order 10/10 and column-shuffled 10/10.

The pipeline is therefore a validated discriminator: it passes on measured
geometric order and recovers true spatial order, and blocks on matched
geometry-free order. Figure 1 shows the confirmatory separation: the geometric
held-out violations cluster near 0.02, far below the 0.05 gate, while both
geometry-free families sit at 0.15-0.38, with no overlap.

![Figure 1](figures/fig1_discriminator.png)

*Figure 1. PC-V1 confirmatory held-out violation by condition at the gate
dimension d = 1. Structured (geometric) scenes pass the frozen 0.05 gate on
every seed; the two geometry-free control families block on every seed. Bars
are per-condition means.*

## 5. Experiment P1: geometry-dilution dose-response

**Design.** We reuse the frozen PC-V1 instrument unchanged and add one element:
an epsilon-diluted order generator. For dilution fraction epsilon, we keep each
non-chain geometric covering edge (transitive-reduction edge) with probability
1 - epsilon, add random time-respecting edges at a bisection-tuned probability
that holds the post-closure relation density at the geometric value, always
keep chain-internal edges, and re-close. epsilon = 0 is the geometric order;
epsilon = 1 is the density-matched geometry-free order of PC-V1. Rewiring
covering edges (rather than closed edges) makes dilution effective, because
removing a covering edge genuinely reduces reachability. Calibration is on
seeds 0-9; confirmation on fresh seeds 300-319 over the grid
epsilon in {0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.0}. The frozen test constant
is the monotonicity floor rho_min = 0.85; a seed enters the test only with >= 6
density-held cells spanning both endpoints (under-covered seeds are recorded,
not counted).

**Results.** Density is held across the grid (mean achieved density 0.566-0.584
at every epsilon), so the response is to geometry, not density. The
dose-response (confirmatory means over density-held seeds):

| epsilon | truth-order error | held-out violation | achieved density |
| --- | --- | --- | --- |
| 0.00 | 0.098 | 0.023 | 0.572 |
| 0.15 | 0.120 | 0.025 | 0.584 |
| 0.30 | 0.148 | 0.026 | 0.581 |
| 0.45 | 0.203 | 0.036 | 0.571 |
| 0.60 | 0.302 | 0.097 | 0.566 |
| 0.75 | 0.446 | 0.203 | 0.575 |
| 0.90 | 0.504 | 0.251 | 0.579 |
| 1.00 | 0.496 | 0.259 | 0.576 |

*Table 2. P1 dose-response: confirmatory means over density-held seeds
(300-319). Relation density is held near 0.57 throughout, so the response is to
geometric consistency, not density.*

All four preregistered criteria are supported on the fresh confirmatory seeds:

- **H-MONO (monotone degradation): supported.** Every one of 20/20 covered
  seeds has Spearman(epsilon, truth-order error) >= the frozen rho_min 0.85
  (per-seed rho 0.929-1.0, median 0.976).
- **H-THRESH (identifiable transition): supported.** The geometry-recovery
  crossing (truth-order error crossing 0.15) is estimable for 20/20 seeds with
  median epsilon* = 0.31. The curve rises smoothly across epsilon ~ 0.3-0.75:
  a graded transition, not a cliff.
- **H-LAG (false-pass window): supported.** The held-out crossing (embeddability
  crossing 0.05) sits at a larger epsilon than the truth crossing; the median
  gap is 0.19 (embeddability crosses at ~0.50, truth at ~0.31). Between these,
  profiles still embed in one dimension while no longer recovering true space.
- **Endpoint reproduction: supported.** 19/20 seeds reproduce both endpoints
  (epsilon = 0 passes the frozen PC-V1 gate; epsilon = 1 blocks), above the
  18/20 rule.

Figure 2 shows the dose-response. Truth-order error rises monotonically from
0.10 to 0.50; held-out violation stays below its 0.05 gate until epsilon ~ 0.5,
crossing later than truth. The shaded band between the two crossings is the
false-pass window.

![Figure 2](figures/fig2_dose_response.png)

*Figure 2. P1 geometry-dilution dose-response (confirmatory seeds 300-319,
density-held cells; bands are +/-1 SD). Relation density is held at ~0.57 across
all epsilon, so the response is to geometric consistency, not density.
Truth-order error (space recovery) crosses its 0.15 gate near epsilon ~ 0.31;
held-out violation (embeddability) crosses its 0.05 gate near epsilon ~ 0.50.
The gap is the false-pass window, where profiles still embed in one dimension
but no longer recover true space.*

## 6. Discussion

On this controlled family at finite scale, effective metric representability of
observer-relative distance order responds monotonically to the amount of
geometric consistency in the underlying causal order, at fixed relation
density. The validated discriminator does not merely separate the two extremes;
it tracks intermediate geometric consistency, with a graded transition near
epsilon ~ 0.31.

The most consequential result is methodological: the false-pass window (H-LAG).
Across a non-trivial range of dilution, the profiles remain embeddable in one
dimension (held-out violation below the 0.05 gate) while the fitted coordinates
have already stopped recovering true space (truth-order error above 0.15). A
procedure that reads embeddability alone as evidence of geometry would, in this
window, confidently report geometry that is not there. This is not a
hypothetical: it explains why an order-first pipeline can produce internally
consistent "representability" verdicts that are nonetheless uninterpretable,
and it shows why a truth-recovery (or equivalent independent) check is
load-bearing rather than optional. In the absence of ground truth on real
data, the discriminator's specificity arm (blocking on matched geometry-free
order) plays the same role: it certifies that a pass is not achievable by the
common-mode route alone.

The two pre-freeze confounds in PC-V1 reinforce the point. The shared-scalar
confound is exactly a common mode masquerading as structure; the degenerate
control is a fairness failure that would have flattered the gate. Both were
invisible to a single pass/looks-fine reading and were caught only because the
protocol required a matched negative control and a density-held foil before any
threshold was frozen.

## 7. Claim boundary

We claim: (a) the described pipeline is, on 1+1D Minkowski causal sets at the
stated scale, a discriminator that passes on geometric order and blocks on
matched geometry-free order; (b) its recovery degrades monotonically as
geometry is diluted at fixed density, with an identifiable graded transition;
and (c) embeddability alone over-reports geometry, so an independent
recovery/specificity check is required.

We do not claim: spacetime emergence; that causal order alone yields metric
scale, signed coordinates, or a manifold; a continuum limit or dynamics; or
that these finite results establish any physical theory. The experiments are
controlled dose-response studies inside a fixed generator family; "pass" means
the representability gate is satisfiable by measured geometric order, and
"block" is the expected outcome on geometry-free order, not a falsification of
any theory.

## 8. Limitations and future work

The study is 1+1D, single diamond geometry, a fixed observer-chain layout, and
a single dilution family (covering-edge rewiring at held density). The absolute
gate values are scale-dependent; the *contrast* (pass vs block, and the
monotone curve) is the transferable content, not the specific thresholds. The
epsilon* estimate is a property of this generator and scale. Natural extensions
are: 2+1D scenes (where one spatial dimension no longer suffices and dimension
selection becomes a live test); boosted or mixed observer layouts; finer
epsilon grids near the transition; alternative dilution families (e.g.,
link-swap at fixed count) as robustness checks; and, most importantly, orders
produced by a candidate dynamics rather than by hand-tuned dilution, which is
where an emergence question could eventually be posed.

## 9. Reproducibility

All results derive from frozen artifacts under `docs/prereg/frozen/`:
PC-V1 thresholds and Stage A/B/C registries; P1 test constants and P1-A/P1-B
registries; and the per-seed CSVs. Preregistrations:
`docs/prereg/pc_v1_positive_control.md` (FROZEN) and
`docs/prereg/p1_epsilon_sweep.md` (FROZEN). Provenance commits: PC-V1
preregistration `4a2d0fa`, instrument repair `9162e8e`, freeze `b77f588`,
confirmatory `891498f`; P1 skeleton `f1bb7d5`, coverage refinement `6b21bb7`,
confirmatory logic `8bc65ae`, freeze `a218d9a`, confirmatory `4c05cf2`.
Calibration runs are deterministic: re-running the calibration script at the
cited commit reproduces the frozen constants. Every output row carries scene
content hashes, seed, stage, and the requested-vs-executed fit budget.

## Appendix A: Fixed parameters

All parameters are frozen in the preregistrations; the pipeline reuses one
scene/measurement/fit configuration across both experiments (only the order
generator differs).

| Group | Parameter | Value |
| --- | --- | --- |
| Scene | events N | 900 |
| | diamond duration T | 2.0 |
| | reference chains K | 6 at x0 = linspace(-0.25, 0.25, 6) |
| | ticks per chain | 96 |
| | target band | \|t\| <= 0.10, \|x\| <= 0.25 |
| | targets per scene | >= 30, subsampled to <= 40 |
| | required bracketing chains | 6 (two-sided) |
| Dissimilarity | centering | parallax (per-target mean over chains removed) |
| | min common columns | 4 |
| | margin | 25th percentile of positive \|D(i,j)-D(k,l)\| gaps |
| Split | pair train fraction | 0.8 (pair-level, leak-free) |
| | constraints | 4000 train / 1000 held-out |
| Fit | embedding dims | {1, 2, 3}; gate dim 1 |
| | steps / restarts | 1500 / 5 (asserted, no silent clamp) |
| | learning rate / batch | 0.05 / full |
| | stability | sign discordance over 4 restart fits |
| PC-V1 gates | held-out / null-gap / stability / truth | <= 0.05 / >= 0.15 / <= 0.15 / <= 0.15 |
| PC-V1 seeds | A / B / C | 0-9 / 100-119 / 200-209 |
| P1 generator | dilution | covering-edge rewiring, re-closed |
| | density hold | post-closure density matched, tol 0.02 |
| | epsilon grid | {0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.0} |
| P1 constants | rho_min / mono-pass / endpoint | 0.85 / 0.8 of covered / >= 18/20 |
| | crossing levels (truth / held-out) | 0.15 / 0.05 |
| | test coverage | >= 6 density-held cells spanning endpoints |
| P1 seeds | A / B | 0-9 / 300-319 |

*Table 3. Fixed parameters, frozen in the preregistrations.*

## References

Verified bibliography: `citations/references.bib` (each entry confirmed
against an authoritative source; see `citations/citation_verification_report.md`).
Anchor sources: Bombelli, Lee, Meyer & Sorkin (1987); Brightwell & Gregory
(1991); Myrheim (1978); Meyer (1988); Malament (1977); Hawking, King & McCarthy
(1976); Kronheimer & Penrose (1967); Sorkin (2005); Surya, Living Reviews in
Relativity (2019); Shepard (1962); Kruskal (1964); Agarwal et al. (2007);
Kleindessner & von Luxburg (2014).
