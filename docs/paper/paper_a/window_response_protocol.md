# Window response: prospective calibration, pilot and confirmation protocol

2026-09-21. This refines the accepted supplementary plan before new pilot data.
G1 instrument checks passed at `f3d4427056f327465fc45758aed281d262cb2850`.
The explicit box coefficients justify retaining the narrow G0 contribution;
they do not establish priority, general reconstruction, or a physical discovery.

## Development and calibration

- Regions: R+, R0, R-, D exactly as the theory note; A=lambda, Du=1.
- Grid A: 0, 1/64, 1/16, 1/4, 1, 2. Fixed pilot N: 200,800.
- Grid quadrature: 16 independent scrambled Sobol replicates, 2^17 triple
  backbones each. All v coordinates are integrated analytically. Use all Sobol
  points, 52-bit scrambling, no thinning. Numerical seeds are independent of
  poset seeds. This is an approximate numerical error analysis, not an interval
  arithmetic certification or an exact finite-replicate coverage theorem.
- Namespace: G1 `[20260921,6101,...]`; grid quadrature `[20260921,6201,region,A]`;
  final numerical forecasts `[20260921,6202,cell,level]`; pilot
  `[20260921,6301,region,N,A,rep]`; confirmation `[20260921,6401,region,N,arm,rep]`.
  The full list suffixes and parameters are recorded in artifacts.
- Fixed pilot: 80 independently generated posets per region/N/A cell, 48 cells,
  3,840 posets total. No cap extension after inspection. Both f and 3-chain
  fraction are retained, but only f enters the primary decision.
- Samplers, raw statistics and per-cell process CPU time are saved. The R+,A=1
  curved-minus-flat mean at each pilot N is the established positive control;
  its two-sided 95% Welch lower limit must be positive (an existence/sign gate).
- Normal forecasts use the exact U-statistic variance formula and numerical
  p,zeta. Flat means use their analytic values. For the flat diamond,
  p=1/10 and zeta=1/100 exactly: integrate the squared normalized future-plus-past
  interval volume over the flat diamond. This is an analytic control, not novelty.
- The pilot assesses numerical errors, variance and tail plausibility. It cannot
  certify the true FPR as 1% or the normal approximation bias as negligible.

## Freeze selection and numerical uncertainty

Use the grid normal forecast at N=300 to find the first positive-A interval
crossing beta=0.5. Its geometric-mean A is used at both N=300 and 1200. If none
crosses, choose the closest positive boundary and report no transition located.
No new A is selected by looking at confirmation data. Across regions, different
A values are not an isolated observation-window contrast.

The final forecast calculation uses 32 independent Sobol scrambles at power 19.
One prespecified refinement at power 20 is allowed only if the propagated beta
error half-width exceeds 0.01; no further numerical extension. Bonferroni t
intervals use per-parameter error 0.05/64, covering the two fixed precision levels
and the distinct final p,zeta parameters conservatively under the t approximation.
The same selected level is used for the entire forecast artifact.
The baseline threshold is the numerical point estimate, frozen as an actual
number. Curved p/zeta uncertainty is propagated through the two normal tails
by bounding each tail over the parameter box. Add the diamond inverse-CDF
discretization comparison to that error budget. Reject numerical feasibility
if the resulting beta half-width is above 0.01. No claim of rigorous coverage
is made for this scrambled-quadrature t approximation.

The diamond quadrature inverse CDF is compared at grids 4097 and 8193 on the
same development unit points. Actual posets use exact-density rejection sampling
and therefore do not inherit that inverse-CDF approximation.

## Power, cost and stopping

The 16 primary intervals remain two-sided Clopper–Pearson with per-interval
alpha=0.05/16: eight flat upper limits <=0.05 and eight beta-difference intervals
entirely inside [-0.10,+0.10]. Include numerical forecast error by widening each
difference interval. A mismatch is a scientific result, never a reason to refit
on the same holdout.

Power is necessarily conditional on a design population. Before pilot data,
declare sensitivity scenarios true FPR in {0.01,0.02,0.03} and curved probability
within +/-{0,0.02,0.04} of the forecast (clipped to [0,1]). These are planning
assumptions, not confidence bounds established by 80 pilot samples. Use the
least favourable endpoint in each scenario; the probability of an integer
binomial acceptance interval is unimodal in its probability parameter.
Size m in {1200,1600,2000,2400,3200,4800,6400,9600} for a joint-power lower bound
of 0.90 at the worst declared scenario, using exact binomial probabilities and
a union bound. This is stricter than the original illustrative +/-0.02 example.
If a pilot interval excludes the declared scenario range, or a core sampler/
mean/variance discrepancy is unexplained, G2 remains unresolved and confirmation
does not start. A pilot interval overlapping a scenario is not its validation.
Report the full sensitivity table, including the limitations of these assumptions.

The resource limits remain pilot/development 8 CPU-hours and confirmation
48 CPU-hours. Price the full sampling, relation and 3-chain pipeline, with 2x
the measured per-cell CPU cost. Choose the smallest m meeting the stated power
and CPU conditions; otherwise stop as infeasible. Confirmation runs once from a
clean committed state, with forecasts, gates, m, A values and seed layout saved
before any of its data exist. Execution authorization is already present.

No claim requires all eight curved conditions to be detectable. Report a cell
as detectable only if its CP lower limit is at least 0.8 and the matching flat
upper limit is at most 0.05. Prediction of low power can pass the accuracy gate.

The numerical precision rule and power scenarios above operationalize choices
left open by the accepted plan; they do not alter its primary tolerances,
confidence family, target power, populations, pilot cap or resource limits.
