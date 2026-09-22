# Prospective finite-order rotating-profile pilot

2026-09-22. Commit this protocol, generator and runner before generating
any new selection or measurement orders. This is a bounded development
experiment, not an accuracy or coverage certification. Its purpose is to
measure whether the unchanged inverse procedure retains the actual
midpoint Weyl information when polarization changes across the diamond,
and to separate finite sampling error from continuum model bias.

## Fixed cases and input boundary

Use the bounded [rotating generator](rotating_generator.md). The ordered
list `(A,w)` is fixed as follows:

```text
(0,0), (0,.5), (.25,0), (.75,0),
(.25,.5), (.75,.5), (.25,-.5), (.75,-.5).
```

These are eight separate orders. The two flat cases check the production
path, including the rotating coordinate convention. The two constant
curved cases check transfer of the generator. The four nonconstant curved
cases test both senses of polarization rotation; H rotates by one radian
in its spin-two components across the full diamond. No flat case, profile,
rotation rate, true coordinates or tensor enters any curved fit.

- Data seed for case i is `2026092600+i` (i=0..7).
- Selection and measurement streams are `default_rng([seed,0])` and
  `default_rng([seed,1])`, respectively.
- Each case has 4096 independent selection events, 262144 measurement
  events, and the supplied outer tips. Total event budget is **2,129,920**
  plus tips. Generate in chunks of 4096; no enlargement or seed replacement.
- Randomize interior labels before constructing the selection order. The
  generator preserves those labels; tip IDs alone are distinguished.
- Keep the connected selector, selected clock mask, quartic joint fitter,
  correlated-membership covariance, gauge, optimizer and ten-component
  confidence calculation exactly as in `e057d10`.
- The archive has only `selection_order`, `selection_design`, `selected_ids`,
  `after`, `before`, and `measurement_n`. It contains neither interval volumes
  nor profile information. The estimator reselects landmark IDs from the order.
- Store generator coordinates and parameters in a separate truth archive.
  The evaluator opens it only after `fit_order_archive` returns.

## Evaluation and failure handling

The target is all ten components of `T^2 C` at the midpoint, in the internal
frame fixed by the first two landmarks and modulo residual z reflection.
Here T=sqrt(2), and the true tidal norm is 2A. Use the true inverse
exponential map for evaluation only. Report the full signed-component
distance, relative error where A>0, nominal confidence radius, covariance
eigenvalue, selected ranks and optimizer status. Flat relative error is
undefined, so report its absolute error.

Compute continuum interval volumes for the SAME selected landmarks, clock
mask and covariance weights. A 32/64-node relative quadrature discrepancy
>=1e-10 is a failure. Refit those mean clocks and report the resulting
full-tensor model bias. Record whether it satisfies the preceding relative
2% bias allowance (absolute tolerance 1e-10 for numerical flat zero).
Also record whether the nominal ball contains the truth. A failure of the
bias condition is evidence against transferring that allowance; do not
inflate it, reduce the tested frequency, add parameters, or repair the
fit within this pilot after observing results.

All selection failures, unresolved floating causal margins, empty selected
measurement clocks, invalid covariance/ranks, and optimizer/boundary
failures remain in the eight-case record. Continue the remaining fixed
cases after a recorded case failure. Preserve full outputs, order/truth
hashes, runtime versions, source hashes and a reachable pre-data commit.
No result from this small development sample is an independent estimate
that coverage is at least 95%, and no point error is an accuracy gate.

The next decision follows the results: if varying-profile continuum bias
exceeds the allowance, establish a regularity/window restriction or a new
local model before an accuracy campaign. If it is controlled, address the
observed order-derived variance and finite-count noise with a separately
powered design. The old provisional variance cap was exceeded by 8 of the
254 completed coverage fits, so it is not a validated cost guarantee.

```powershell
rtk proxy .venv/Scripts/python.exe experiments/positive_control/rotating_plane_wave.py
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_rotating_pilot.py
```
