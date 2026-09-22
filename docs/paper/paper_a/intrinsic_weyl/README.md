# Local Weyl reconstruction from interval cardinalities

Development record, 2026-09-22. The active research objective is to determine
whether Weyl tidal information can actually be recovered from causal order,
without supplied coordinates, a wave profile, or a paired flat control.
**That objective is not yet achieved.**

The proposed lever is to retain the geometry of multiple intervals instead of
averaging it into a single chain statistic. Their first-order incompatibility
with flat geometry contains signed Weyl information. This record establishes
a local algebraic route and tests its continuum bias. The original checkpoint
left finite-order recovery, order-only landmark selection, and calibrated
uncertainty open; the update below records which parts have since progressed.
This is not a new theory of spacetime, a universal global formula, or a
preregistered accuracy result. Earlier frozen verdicts are unchanged.

Update: the subsequent [joint finite-order pilot](joint_pilot_results.md)
implements order-derived landmark selection and a joint geometry/Weyl fit.
It produced quantitative tidal spectra from six saved finite orders, with a
separately replayed estimator that receives no generator data. Independent
accuracy confirmation and varying-profile validation are still outstanding.
The [three-band follow-up](three_band_pilot_results.md) adds six further
orders and exposes a selection failure: the central band was mostly removed.
Reliable recovery of the magnetic components is not established.
The [connected selection pilot](connected_pilot_results.md) repaired central
retention in six further independent orders. Both spectra remain noisy, and
an expanded evaluation checks all ten components and relative directions.
On the four curved cases, relative full-tidal errors are 31--112% from finite
orders versus 0.24--0.62% from continuum volumes on the same landmarks.
These are development point errors, not confidence or accuracy verdicts.
The subsequent [independent coverage protocol](coverage_protocol.md) fixes
256 new orders, an internal frame modulo reflection, and a simultaneous
ten-component confidence ball with an explicit bias assumption. The
[completed run](coverage_results.md), generated from `e057d10`, returned
254 valid-and-covered cases out of 256; its exact 95% two-sided interval
is [0.97207, 0.99905], meeting the frozen coverage gate. Median relative
full-tidal errors remain 89.54% (weak) and 32.31% (strong), and the balls
are too wide for 25% accuracy. Eight cases exceeded the provisional later
variance envelope, so its old sample-price assumption is not validated.
The [rotating-profile finite-order pilot](rotating_pilot_results.md) then
completed eight independent orders using the unchanged estimator. All four
nonconstant curved cases exceeded the old 2% continuum-bias allowance
(observed 2.70--4.74%), despite broad nominal balls containing the truth.
The next lever is the magnetic information in the landmark geometry,
together with an explicit varying-profile bias budget; increasing event
count alone does not address both limitations.
The [fixed-event paired design comparison](design_pilot_results.md) tested
wider centers on eight fresh shared orders. It did not improve magnetic
uncertainty and increased actual tensor error in seven pairs; that selector
is not adopted. A separate [gradient diagnosis](gradient_model_diagnosis.md)
reduced continuum bias through oracle subtraction of the missing fifth-degree
term. An order-only local information audit retains rank ten when all 24
vacuum gradient components are freed, at a nominal variance cost of
13.5--48.4% in four checked configurations.
The subsequent [blind curvature-and-gradient pilot](jet_pilot_results.md)
implements that enlarged inverse and tests eight fresh orders. On all four
varying-profile cases its continuum bias is 1.19--1.83%, below the unchanged
2% allowance; paired quartic biases are 2.14--4.31%. Actual finite-order
errors and confidence widths remain large. The enlarged-model coverage,
regularity scope and powered accuracy demonstration are still unfinished.
The [independent fifth-degree coverage protocol](jet_coverage_protocol.md)
now freezes 256 fresh cases from an iid constant/rotating-wave mixture,
including selection, fit and bias failures in its confidence gate. Its
protocol is not a result; no coverage or accuracy success is claimed here.
The [local error derivation and fixed-scale check](local_consistency_scope.md)
explain the possible rate benefit: under explicit regularity, conditioning
and local-solution assumptions the physical-curvature bias becomes quadratic
in window size. A four-scale continuum audit observes successive bias
reductions approaching four when the window is halved. These conditions
are not yet proved for the random selector and implemented optimizer.
The [curvature-squared correction](square_correction.md) adds known distance
and volume-clock terms without adding fitted components. On three selected
coverage failures it reduces continuum bias to 1.25--1.32%, but actual finite
errors remain 28--53% and do not improve consistently. The
[fresh paired pilot](square_pilot_protocol.md) freezes eight independent
orders before testing this correction. Neither selected-case diagnosis nor
its protocol establishes coverage or accuracy.

The [finite-resolution obstruction](finite_resolution_obstruction.md) gives
an explicit oscillating vacuum family: midpoint tidal strength stays fixed
while finite-order laws approach the flat law as frequency grows. Thus a
uniform finite-sample guarantee for unrestricted smooth profiles is impossible.
Positive reconstruction on a stated regularity class, including varying
profiles and tested uncertainty, remains to be established.

## 1. What the input can determine

Assume a smooth four-dimensional Ricci-flat spacetime and uniform sampling with
respect to its physical volume. All intervals used below must lie in a regular
convex normal neighborhood; their shape stays bounded as its scale shrinks.
The estimator sees interval cardinalities or the order from which to count
them. It does not see the embedding, the metric profile, or a flat sample.
The original Schur construction below takes landmark IDs from the caller;
the subsequent joint pipeline selects its own landmarks from a separate order.

Absolute dimensional curvature is not identifiable without an extra scale:
under `g -> lambda^2 g`, the causal order and normalized uniform sampling law
are unchanged, while orthonormal curvature components scale as `lambda^-2`.
Consequently the target is `tau_ref^2 C_abcd`, where the reference time is
itself determined by an interval. Unknown common density cancels in ratios.
An arbitrary unknown spatially varying sampling density does not cancel and
is outside this model.

The observer is the initial unit tangent from landmark 0 to landmark 1.
Three further reference vectors fix spatial axes by Lorentz Gram-Schmidt.
Thus the tensor is specified in an internal, physically defined frame. Its
electric eigenvalues and electric norm can also be reported. A tensor without
a stated observer/frame is not the target.

## 2. Eliminate the unknown flat embedding

Let `epsilon` denote the size of the neighborhood in a bounded local frame,
`eta = diag(-1,1,1,1)`, and `S_ij = -tau_ij^2` for timelike pairs. Normal
coordinates are used for this derivation only. With the curvature convention
fixed by `C_uiuj = -H_ij` for the Brinkmann test geometries,

```text
S(x,y) = (x-y)^T eta (x-y) - C(x,y,x,y)/3 + O(epsilon^5).
G_xy   = (S_0x + S_0y - S_xy)/2
       = x^T eta y + C(x,y,x,y)/6 + O(epsilon^5).
```

This is the standard world-function expansion, not a new curvature identity;
see Brewin [1], section 9. For four independent reference vectors `a_i`, put
`G_A = (G_ai,aj)`, `g_x = (G_ai,x)`, and `c_x = G_A^-1 g_x`.
Factor the observed Lorentzian Gram matrix and reconstruct
`x_hat = sum_i c_xi a_hat_i`. Lorentzian Gram embedding is also established
methodology [2]. Its error relative to a fixed normal frame is `O(epsilon^3)`.

The observable Schur residual is

```text
D_xy = G_xy - g_x^T G_A^-1 g_y.
```

Expanding the inverse Gram matrix gives the linear response

```text
6 D_xy = C(x,y,x,y)
       - sum_i c_yi C(a_i,x,a_i,x)
       - sum_i c_xi C(a_i,y,a_i,y)
       + sum_ij c_xi c_yj C(a_i,a_j,a_i,a_j)
       + O(epsilon^5).
```

The quartics use reconstructed vectors in the implementation. Their
`O(epsilon^3)` coordinate errors change this equation only at order six.
The ten unknowns are the five independent components of the electric and
five of the magnetic Weyl tensor. No plane-wave algebraic type or polarization
is imposed on the estimator.

In vacuum, a small diamond has

```text
V_ij = (pi/24) tau_ij^4 [1 + O(epsilon^4)],
sqrt(V_ij/(pi/24)) = tau_ij^2 + O(epsilon^6).
```

The Ricci term is zero and a single scalar volume has no linear Weyl term.
Replacing `S` by `-sqrt(n_ij/n_04)` therefore preserves the quartic response,
in squared reference-time units. This does not assert that curved diamonds
have the same volume as flat ones. The order-six clock error is retained in
the bias. Nonvacuum Ricci terms enter earlier and require a different fit.

## 3. Identifiability, including a failed reduction

`check_identifiability.py` computes the response in exact rational arithmetic.
For the ten-event, all-timelike configuration saved in
`identifiability_certificate.json`, the 15-by-10 matrix has rank ten. One
10-by-10 minor is exactly

```text
1713617480891173 / 315179916321896792064 != 0.
```

Because this minor is a nonzero polynomial after clearing denominators, it
also proves local generic identifiability within the open set of admissible
configurations. It is not an accuracy or conditioning guarantee.

For the SAME four reference vectors, keeping only probe-to-reference clocks
leaves a rank-eight polynomial map even with arbitrarily many probe locations.
The two exact null vectors are in the certificate. Cross-probe intervals
remove this obstruction. This rank-eight result is not a theorem about all
possible reference configurations or methods.

## 4. Continuum checks and development changes

`intrinsic_weyl_development.py` uses exact plane-wave interval volumes as a
labelled oracle. It reuses the existing transfer-matrix integrator. A separate
fit using exact proper times helps distinguish volume-clock bias from the
curvature truncation. Only a matrix of volumes enters the estimator.

The present 33-event configuration uses four references and 28 probes with
wider spatial spread. At the flat limit its 224 usable residual rows have
rank ten and a design condition number of approximately 21.28.

| Neighborhood scale | Constant wave: relative tensor error | Rotating, varying wave: relative tensor error |
|---|---:|---:|
| 0.04 | 0.029564 | 0.627712 |
| 0.02 | 0.007776 | 0.313227 |
| 0.01 | 0.001969 | 0.155238 |
| 0.005 | 0.000494 | 0.077135 |

The constant wave shows approximately quadratic bias decay; the varying wave
shows approximately linear decay, as the local expansion predicts. Flat
absolute tensor errors are at most `1.1e-12`. These are deterministic checks
on the displayed configurations, not confidence intervals or validation on
finite causal sets. Source hashes and full outputs are in
`continuum_diagnostic.json`.

Development history matters: an initial ten-event thin configuration had
condition number about 364 and constant-wave relative errors approximately
22.88, 6.221, 1.579, 0.3962 on the same scale ladder. The wider arrangement
was chosen AFTER seeing that failure. No preregistration is claimed for this
choice or the numerical table above. The exact rank certificate deliberately
keeps the thin configuration: identifiable does not mean statistically useful.

Independent mechanical checks cover the Weyl symmetries and trace, exact
rank, interval extraction from an order, invariance under relabelling, and
invariance under common rescaling of all counts. The combinatorial order
test checks count extraction only; it does not recover curvature.

## 5. Noise is the next decision, not an afterthought

For a fixed nondegenerate landmark shape, a Poisson intensity `rho` gives
counts of order `rho epsilon^4`. Squared-time noise is of order
`rho^-1/2`, while the unnormalized curvature response is quartic. Conditional
on valid independently chosen landmarks and controlled shape,

```text
curvature bias  = O(epsilon)       [O(epsilon^2) in the constant test],
curvature noise = O_p(rho^-1/2 epsilon^-4).
```

Thus the candidate consistency regime is `epsilon -> 0` and
`rho epsilon^8 -> infinity`. This scaling argument is not a completed
finite-sample theorem: landmark selection and overlapping-interval covariance
must be handled. Absolute error of the dimensionless tensor is insufficient
as a recovery claim, because the dimensionless target itself shrinks.

### Prospective development cost diagnostic

The following is fixed before this diagnostic's random data are generated.
It estimates an optimistic local noise price and has NO scientific accuracy
gate. It does not run a curved recovery campaign.

- Keep the 33 landmarks above and ordinary least squares unchanged.
- Use the existing uniform 4D Minkowski-diamond sampler, outer endpoints
  `(0,0,0,0)` and `(11,0,0,0)`, containing all landmark intervals.
- Four independent seeds `2026092200..2026092203`, each with `262144`
  events; process chunks of `16384`. No extension after observing results.
- Estimate the exact flat-limit influence covariance from shared event
  membership in ALL intervals. Do not replace overlapping interval counts by
  independent Poisson draws. Preserve the four replicate covariance results.
- Derive the count-to-fit Jacobian analytically at the flat limit and compare
  deterministic directional derivatives with finite differences before
  drawing random points. Common count scaling must lie in its nullspace.
- Report predicted full-tensor RMS at `N=1e6` for scales 0.04 and 0.02.
  Also report an optimistic normal-model price for a SINGLE nonzero tensor
  projection: a 95% two-sided interval contained in a relative +/-25% band
  with 90% probability under zero bias. Its factor is
  `(z_0.975 + z_0.95)^2 / 0.25^2`, not the break-even `z_0.975^2` factor.
  This omits truncation bias and is not a full-tensor certification price.
- Reconstruct the flat tensor from each replicate's actual counts, retaining
  signature/rank failures as failures. This checks the noisy input path; it
  does not establish curved recovery or confidence coverage.
- Bind the run to a clean, reachable Git commit and source hashes. Save all
  outputs, including unfavorable prices. No flat-control sample enters any
  curved estimator; the flat computation here is a separate design study.

If this price is prohibitive, the next lever is the statistical design
(landmark geometry, joint fitting, and covariance weighting), not a larger
uncontrolled sampling run. A confirmatory study would separately freeze
landmark selection, estimand, bias treatment, 95% two-sided construction,
90% power, independent validation, and a realistic budget before data.

## 6. Relationship to prior work

Checked against primary sources on 2026-09-22; this is a bounded search, not
a proof of novelty or an exhaustive priority audit.

1. [Brewin, *Riemann Normal Coordinate expansions using Cadabra*, v3](https://arxiv.org/html/0903.2087),
   section 9: the squared-geodesic-distance expansion used here.
2. [Clough and Evans, *Embedding graphs in Lorentzian spacetime*](https://arxiv.org/html/1602.03103),
   section 2.2: indefinite Gram embedding is established.
3. [Puetzfeld, Obukhov and Laemmerzahl, *Gravitational clock compass in General Relativity*](https://arxiv.org/html/1805.10673):
   curvature from appropriately arranged clocks is established. The abstract
   and framework were checked; no exhaustive equation-by-equation comparison
   has yet been completed.
4. [Shuman, *Partial Derivatives on Causal Sets*](https://arxiv.org/html/2405.10574v1),
   and [*The Geometry of Tangent Spaces on Causal Sets*](https://arxiv.org/html/2405.11805v1),
   sections I.2 and II: already combine order-derived coordinates, proper
   times, metric estimation and curvature. The reported numerical curvature
   examples are flat 1+1D, where Weyl is identically zero. We must not claim
   the first inverse-geometric or coordinate-reconstruction approach.
5. [Ibrahim Ahmed and Bombelli, APS 2026 abstract](https://meetings-archive.aps.org/smt/2026/apr-r93/1/):
   reports Weyl indications through discrete geodesic deviation in plane
   gravitational waves. It blocks a broad first-Weyl-detection claim; an
   abstract alone does not permit a full estimator or precision comparison.
6. For chain-volume reconstruction and the Alfyorov--Shnyukov observable,
   retain the estimand-level comparisons in
   [the existing literature audit](../window_response/plane_wave_chain_literature_and_limits.md).

The prospective contribution is therefore specific: an internally framed,
signed 4D Weyl reconstruction from multiple timelike interval counts, with
finite-order uncertainty and cost demonstrated. The last two requirements
are still missing; this is not yet a submission-ready claim.

The subsequent [noise-cost decision record](noise_design.md) finds an
unaffordable cost for the current ordinary least-squares geometry and a
315--355-fold reduction in a separately validated, linearized covariance-
weighting diagnostic. The remaining price is still about two billion events
for one projection at scale 0.04 under optimistic assumptions. This is a
design finding, not an actual curved recovery result.

## Reproduction

From the repository root, without adding dependencies:

```powershell
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_development.py
rtk proxy C:/Users/lpaiu/AppData/Local/Programs/Python/Python311/python.exe docs/paper/paper_a/intrinsic_weyl/check_identifiability.py
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_noise.py
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_weighting.py
```

The second command uses the already installed SymPy environment. All new
code is isolated from frozen experiments and requires no package changes.
The noise diagnostics require a clean worktree and record their generating
commit. To reproduce the saved outputs' source versions, use respectively
`5dc22a9` and `a63e490` in a separate checkout. Historical source hashes are
checked against those reachable commits, not against later edited documents.
