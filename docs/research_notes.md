# Research Notes

## Milestone 18: State-Change Causal Order Refactor

Milestone 18 is a theory refactor. It reorganizes the earlier reconstruction
layers around locally finite state-changing events and causal trigger order.
Observer chains, radar-time ranks, and observer slices become observer-derived
protocol structures. Observer-slice-relative distance order is then the spatial
order layer, while metric geometry, seconds, meters, ratios, velocity,
curvature values, and quantitative dynamics remain calibrated effective
representations.

This milestone does not add a state-change simulation or quantum model.
Quantum compatibility is marked as future work requiring an additional
amplitude or Hilbert-space layer.

## Milestone 19: Minimal State-Change Toy Model

Milestone 19 implements the first finite state-change causal trigger network.
Events are local state transitions, and immediate trigger edges form a finite
DAG. Causal-order queries use the transitive closure.

The validation checks are deliberately order-theoretic: acyclicity,
irreflexivity, transitivity, finite intervals, local chain lengths, and a local
reference-chain diagnostic. No metric reconstruction, observer extraction,
finite-speed spatial geometry, quantum amplitude layer, or curved-spacetime
model is introduced.

## Milestone 20: Reference-Chain Selection

Milestone 20 studies reference-chain utility inside finite state-change
trigger networks. A local-system chain can serve as a candidate reference
protocol, and order-only greedy, longest-chain, and random baselines are also
tested.

The diagnostics are coverage, two-sided bracketing, interval-profile
regularity, local-system purity, top-score gaps, and candidate overlap. These
are finite reference-chain utility diagnostics. They do not assign calibrated
time, identify a unique observer, or reconstruct metric distance.

## Milestone 21: Reference-Chain Bracket Diagnostics

Milestone 21 computes order-level brackets from selected reference chains. For
each target event, the code records predecessor and successor reference
positions, two-sided accessibility, radar-time rank, bracket-width rank, and
rank slices.

These are rank-level diagnostics for a chosen reference protocol. The
bracket-width rank is not metric distance, and the rank slices are not global
time slices.

## Milestone 1: Timelike Reconstruction In 1+1D

This milestone works in 1+1D Minkowski spacetime with natural units `c = 1`.
Events are represented as `(t, x)`.

The causal relation is:

```text
p precedes q iff dt > 0 and dt^2 >= dx^2
```

where `dt = q.t - p.t` and `dx = q.x - p.x`.

The experiment samples events inside the causal diamond bounded by:

```text
p = (-T/2, 0)
q = ( T/2, 0)
```

The diamond condition is:

```text
abs(x) <= T/2 - abs(t)
```

For this 1+1D interval, the continuum volume is:

```text
V = tau^2 / 2
```

so an event-density estimate `rho` and an interval count `N` give:

```text
tau_hat = sqrt((N / rho) / (1/2))
```

The density is not derived from causal order alone. It is extra metric-scale
structure supplied to the reconstruction procedure.

## Interpretation

The current experiment is a controlled check that causal interval cardinality
and longest chains behave consistently in a known spacetime. It is not a claim
that causal order alone fixes all metric quantities, nor that finite signal
speed by itself implies Lorentzian spacetime.

## Milestone 2: Internal Timelike Pair Validation

The original full-diamond endpoint experiment is useful as a sanity check, but
it can be close to tautological if the same full-diamond event count is used to
set the density and then reconstruct the full-diamond proper time.

The milestone 2 validation experiment instead:

```text
1. generates one larger 1+1D Minkowski causal diamond,
2. treats global event density as supplied scale information,
3. samples many internal timelike pairs using the causal order,
4. counts each pair's Alexandrov interval cardinality,
5. compares reconstructed proper time to the hidden Minkowski value.
```

For an internal pair `(i, j)`, causal order contributes the set of sampled
events `k` satisfying:

```text
i precedes k and k precedes j
```

Event density contributes the scale needed to convert that count into a volume
estimate. In 1+1D, the operational reconstruction being tested is:

```text
tau_hat = sqrt((interval_count / rho) / (1/2))
```

Longest-chain estimates are also reported with the asymptotic 1+1D relation:

```text
L ~ sqrt(2 * rho) * tau
```

The implementation marks this normalization as finite-size sensitive. These
experiments are validation in known 1+1D Minkowski spacetime, not a proof of a
new spacetime theory.

## Milestone 3: Probe-Pair Statistical Calibration

Milestone 3 separates support-event sampling from endpoint selection. The
support events are sprinkled in a known 1+1D causal diamond, while probe
timelike endpoint pairs are sampled independently and are not inserted into the
support set.

For each probe pair, the code counts support events in the probe pair's
Alexandrov interval and uses:

```text
tau_hat = sqrt(2K / rho)
```

where `rho` is supplied from the known global support density. The hidden probe
coordinates are used to validate the result and to compute expected
finite-sampling error scales.

For Poisson support-event fluctuations:

```text
E[K] = rho * tau^2 / 2
std(tau_hat) ~= 1 / sqrt(2 rho)
```

For the fixed-count sprinkling used by the code, the interval count is better
approximated as binomial with:

```text
p = V_interval / V_global
K ~ Binomial(N_support, p)
```

The probe-pair experiment compares observed RMSE against both Poisson and
fixed-N binomial delta-method predictions. This helps distinguish ordinary
finite-sampling noise from possible bias, boundary effects, estimator mistakes,
or convention errors. It does not imply that causal order alone determines
metric scale.

## Milestone 4: Dimension As An Order-Statistical Observable

Milestone 4 adds D-dimensional flat Alexandrov intervals and the
Myrheim-Meyer dimension estimator. The validation question is whether dimension
can be estimated from causal-order relation statistics in controlled known
intervals.

For a finite causal matrix `C`, the ordered relation fraction is:

```text
r_obs = (# ordered relations) / (N * (N - 1))
```

The Myrheim-Meyer curve used here is:

```text
r(D) = Gamma(D + 1) * Gamma(D / 2) / (4 * Gamma(3D / 2))
```

The experiment inverts this curve by bisection to estimate `D`. This tests
dimension as an order-statistical observable. It does not prove spacetime
emergence, and it still assumes the events were generated inside known flat
Minkowski Alexandrov intervals.

The theory-facing hypothesis is a mathematical reconstruction program:

```text
primitive causal/information-accessibility structure
  + counting measure / event density
  + observer protocol
  -> operational time, distance, dimension, and spacetime-like geometry
```

The project separates established mathematics, implemented numerical models,
reconstruction hypotheses, speculation, and open proof obligations in
`docs/theory/`.

## Milestone 5: Observer-Chain Radar Reconstruction

Milestone 5 moves radar decomposition away from direct coordinate formulas. A
stationary observer chain with clock labels is supplied, and the reconstruction
uses only:

```text
causal matrix + observer chain indices + observer clock labels
```

For a target event, the protocol identifies:

```text
tau_minus = latest observer clock tick preceding the target
tau_plus  = earliest observer clock tick following the target
```

Then:

```text
radar_time = (tau_plus + tau_minus) / 2
radar_distance = (tau_plus - tau_minus) / 2
```

Hidden coordinates are used only for validation against the known stationary
observer values `t` and `abs(x)`. This supports the observer-dependent spatial
decomposition layer of the reconstruction program, while keeping the observer
chain and clock labels explicit as additional protocol structure.

## Milestone 6: Oriented Radar And Lorentz-Map Recovery

Milestone 6 makes explicit that one observer chain reconstructs an unsigned
radar distance. For a stationary chain at `x = 0`, the events `(t, x)` and
`(t, -x)` have the same radar time and radar distance, so a single chain cannot
recover a signed spatial coordinate.

The oriented protocol supplies a second synchronized comoving beacon chain at
rest-frame separation `a`. Radar distances to the primary chain and beacon
chain are combined as:

```text
x_hat = (R0^2 - Ra^2 + a^2) / (2a)
```

The beacon separation and synchronization are additional protocol structure.
They are not derived from causal order alone.

The Lorentz-map experiment builds two oriented protocols for the same causal
set: a stationary lab protocol and a moving inertial protocol. Each protocol
uses causal-order-based radar ticks and supplied clock labels. Hidden
coordinates are used only to validate the reconstructed coordinates and the
fitted boost parameter. The validation question is whether the compatibility
map between the two reconstructed coordinate decompositions approaches the
known Lorentz map as clock tick density increases.

## Milestone 7: Observer Atlas Consistency

Milestone 7 moves from isolated oriented protocols to a small atlas of
overlapping reconstructed charts. Each chart is built from the same causal
matrix plus supplied protocol chains, clock labels, beacon separation, and
observer origin data.

For chart pairs, the experiment fits an affine Lorentz/Poincare transition map:

```text
chart_B ~= L(beta_AB) chart_A + translation_AB
```

It compares the fitted relative beta to the value expected from the controlled
hidden protocol construction. It also checks whether the two charts assign
approximately the same Minkowski interval squared `dt^2 - dx^2` to sampled
pairs of shared events.

For the chart loop `A -> B -> C`, the experiment composes fitted maps and
compares that composed result with the direct fitted map `A -> C`. This is an
atlas-consistency diagnostic in known 1+1D Minkowski intervals. It does not
show that causal order alone gives a manifold, a unique atlas, signed space, or
affine translations.

## Milestone 8: Rindler Horizon And Reconstruction Inaccessibility

Milestone 8 adds a uniformly accelerated Rindler observer in flat 1+1D
Minkowski spacetime. The observer chain and proper-time labels are supplied
protocol structure. The reconstruction uses only the causal matrix, observer
indices, and clock labels to determine whether each support event has both
radar ticks.

The validation compares three notions:

```text
1. analytic infinite-wedge accessibility,
2. finite-chain coverage accessibility,
3. reconstructed two-way radar accessibility.
```

Events inside the ideal Rindler wedge may still be inaccessible for a finite
observer chain if their analytic `tau_minus` or `tau_plus` falls outside the
simulated clock range. This is reported separately so finite chain duration is
not mistaken for the true Rindler horizon boundary.

This is a controlled flat-spacetime horizon analogue. It is not a black hole
simulation and does not prove spacetime emergence.

## Milestone 9: Conformal Ambiguity And Measure Dependence

Milestone 9 demonstrates a limitation of causal order. For any positive
conformal factor `Omega(t, x)`, the metric

```text
ds_Omega^2 = Omega(t, x)^2 * (dt^2 - dx^2)
```

has the same null cones and the same causal order as the flat metric. The
physical volume element and clock rates change:

```text
dV_Omega = Omega^2 dt dx
d tau_Omega = Omega dt
```

The constant-scale experiment shows that relation fraction and Myrheim-Meyer
dimension remain unchanged while proper time and volume change under the same
causal matrix. The weighted-volume experiment then supplies `Omega^2` as
additional measure information and checks whether weighted interval counts
recover conformal physical volume better than unweighted coordinate counts.

This supports the reconstruction program by clarifying its assumptions:

```text
causal/accessibility order -> conformal/light-cone structure
causal/accessibility order + measure/density -> volume-scaled reconstruction
```

It does not derive the conformal factor from causal order.

## Milestone 10: Measure Encoding And Coarse-Graining Stability

Milestone 10 extends the measure story in two controlled ways. First, support
events can be sprinkled uniformly with respect to conformal physical volume
rather than coordinate volume. For a positive conformal factor in 1+1D,

```text
dV_Omega = Omega^2 dt dx
```

so the coordinate-space event density is proportional to `Omega^2`. If the
physical density `rho_physical = N / V_global_Omega` is supplied, unweighted
interval counts estimate conformal physical volume. If a coordinate density is
used instead, non-flat physical measures produce a predictable bias.

Second, the milestone tests random thinning as a coarse-graining diagnostic. If
events are retained with probability `p`, the density must be rescaled:

```text
rho_thinned = p * rho_original
```

Corrected or realized density estimates should remain approximately stable
under thinning, with larger finite-sampling noise at lower keep probability.
Using the original unthinned density after thinning should underestimate
volume.

The local measure-profile experiment also demonstrates a limitation: nonuniform
relative `Omega^2` shape can be recovered statistically from event counts, but a
constant conformal scale has the same normalized coordinate distribution as the
flat profile. Absolute scale still requires additional density or volume
information.

## Milestone 11: Order-First Reformulation

Milestone 11 reinterprets the previous reconstruction experiments as
representation-layer checks. The order-first thesis is:

```text
causal order -> primitive temporal order
observer protocol + causal order -> observer-relative distance order
order structures + calibration/dynamics -> metric representation
```

In this framing, seconds, meters, ratios, metric tensors, and curvature values
are not primitive. They are effective objects that become useful when order
structures admit stable calibration, observer-atlas consistency, and
representability conditions.

The new ordinal experiments check:

- radar-return distance order from causal order and observer tick order,
- preservation of order under positive monotone transformations,
- failure of arbitrary monotone transformations to preserve ratios,
- ratio stabilization under equal-step calibration,
- distance-order preservation in oriented radar charts,
- finite diagnostics for order cycles, triangle inequality failure, and
  Euclidean candidate embeddings.

These checks do not show that distance order alone gives metric geometry. They
separate order preservation from metric-value accuracy.

## Milestone 12: Ordinal Embedding And Effective Metric Representation

Milestone 12 asks when distance-order constraints admit a useful
low-dimensional metric representation. The controlled ordinal embedding
experiments use constraints of the form:

```text
d(i,j) < d(k,l)
```

and fit low-dimensional coordinates by minimizing a finite quadruplet hinge
loss. The fitted coordinates are evaluated by violation rate, order-preservation
loss, and Procrustes-aligned RMSE against hidden validation coordinates.

This milestone treats metric geometry as a low-complexity compression of rich
order data. More constraints should generally improve a consistent embedding,
candidate dimensions below the effective dimension should show higher ordinal
stress, and noisy flipped comparisons should create irreducible violations.

The observer-derived experiment uses oriented radar reconstruction to generate
distance-order constraints from a supplied observer protocol, then asks whether
those constraints support a 1D effective spatial representation. This remains a
finite diagnostic; it does not prove that distance order alone gives full
metric geometry.

## Milestone 13: Representation Stability And Null Models

Milestone 13 adds held-out order validation, shuffled and random null-model
baselines, controlled noisy constraints, and subset-stability diagnostics. A
metric representation is treated as useful only if it generalizes to held-out
order constraints and is reasonably stable across independent subsets of the
same structured order data.

The null models are baselines, not physical alternatives. They test whether the
ordinal embedding optimizer can partially fit constraints that do not encode a
consistent geometric order. Structured geometric and observer-derived order
data should outperform shuffled or random constraints on held-out validation if
they support a low-complexity effective metric representation.

These checks strengthen the finite diagnostic standard. They do not prove
spacetime emergence, derive meters, or make metric geometry fundamental.

## Milestone 14: Simultaneity-Sliced Distance Order

Milestone 14 refines observer-relative spatial order by adding an explicit
slice-selection protocol. Radar tick brackets define an order-level radar-time
rank, and radar-time bins define same-slice target sets. Spatial distance-order
comparisons are then made only within those observer-derived slices.

This avoids mixing events from different observer times when constructing
spatial distance-order constraints. Hidden coordinates are used only for
validation, not for choosing slices or generating observer-derived constraints.

## Milestone 15: Cross-Slice Identification And Transport

Milestone 15 separates slice-local spatial order from cross-slice predicates.
Without transport, anchors, persistence, calibration, or dynamics, statements
about same position, same direction, velocity, constant velocity, or spatial
evolution are undefined. With a chosen protocol, those statements become
transport-relative.

Anchor and persistence experiments test how supplied structure constrains
independent per-slice translation, reflection, orientation, and scale freedom.
Noisy-transport diagnostics show that derived cross-slice quantities degrade
when the supplied protocol is uncertain.

## Milestone 16: Relational Spatial Evolution

Milestone 16 asks which cross-slice statements can be made without
same-position transport. With supplied persistent object labels and
slice-local distance order, pair-distance order histories define ordinal shape
signatures across slices. Changes in those signatures are relational spatial
evolution: persistence-dependent but transport-independent.

This does not make persistence primitive, and it does not define velocity or
metric dynamics. It identifies weaker transport-gauge content that survives
independent per-slice affine/reflection transformations.

## Milestone 17: Persistence Ambiguity And Identity Matching

Milestone 17 studies object persistence itself. Without supplied persistence,
cross-slice object identity and pair-distance order histories are undefined.
If a matching criterion is used, the resulting relational history is defined
only relative to that persistence hypothesis.

The experiments enumerate finite identity matchings, show that symmetric
configurations can have multiple equally good matchings, test when
relational-continuity matching recovers hidden validation labels, and check how
partial labels restrict ambiguity. Crossing histories are included as failure
cases. These diagnostics do not derive object identity from causal order.

## Radar Coordinates

For the stationary observer worldline `O(tau) = (tau, 0)`, the radar
decomposition of an event `(t, x)` is:

```text
tau_minus = t - abs(x)
tau_plus  = t + abs(x)
radar_time = (tau_plus + tau_minus) / 2
radar_distance = (tau_plus - tau_minus) / 2
```

For an inertial observer moving with velocity `v`, the event is first Lorentz
boosted into that observer's rest frame:

```text
t' = gamma * (t - v x)
x' = gamma * (x - v t)
```

The same stationary radar decomposition is then applied to `(t', x')`. This is
an operational coordinate protocol relative to a chosen observer.

## Lorentz Length Contraction

For a rod at rest in `S'` with endpoints `x'=0` and `x'=L0`, the transform to a
lab frame `S` where the rod moves at velocity `beta` is:

```text
t = gamma * (t' + beta x')
x = gamma * (x' + beta t')
```

A lab-frame length measurement compares endpoint events with the same lab time.
At lab time `t=0`, the selected rest-frame endpoint events are:

```text
left:  t'=0,        x'=0
right: t'=-beta L0, x'=L0
```

Transforming these two events gives a lab-frame separation:

```text
L = L0 / gamma
```

This is the usual special-relativistic length contraction as an event-selection
result, not a claim about a new spacetime model.

## Finite-Speed Lattice Counterexample

A regular 1+1D lattice graph can impose a finite propagation rule:

```text
(t, x) -> (t + 1, x - 1)
(t, x) -> (t + 1, x + 1)
```

This graph has a finite causal cone, but it also has preferred lattice
directions, a preferred time slicing, and a parity structure. Comparing the
lattice cone with a Poisson sprinkle in a continuum causal cone is a simple
counterexample to the claim that finite signal speed alone implies
Lorentz-invariant spacetime behavior.

The comparison is intentionally modest. It shows that extra structure matters;
it does not claim that a Poisson sprinkle is a complete physical model.

## Exploratory Spacelike-Distance Proxies

Spacelike distance is harder to reconstruct from causal order than timelike
proper time in this toy setting. A timelike interval has a direct Alexandrov
interval cardinality, while a spacelike pair has no causal interval between the
two events themselves.

The current exploratory implementation records three simple proxy counts for
spacelike-separated event pairs:

```text
common future overlap count
common past overlap count
minimal enclosing Alexandrov interval count
```

These quantities are compared against the true 1+1D Minkowski spacelike
separation:

```text
sqrt(dx^2 - dt^2)
```

The proxies are strongly affected by the finite sampling region, event density,
and boundary placement. They are useful diagnostics for failure modes and
possible correlations, but they are not validated estimators of spacelike
distance.
