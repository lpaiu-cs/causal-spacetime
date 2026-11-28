# Missing Weyl-gradient term: bias diagnosis and order-derived information

2026-09-22. The wider-center comparison failed to improve the dominant
uncertainty. Two separate calculations now test a specific next mechanism.
The first uses oracle information to diagnose bias. The second uses only
order-derived fitted quantities to assess whether freeing the missing
derivatives destroys local identification. Neither is an operational
gradient-reconstruction result or an accuracy verdict.

## The omitted term is established geometry

For normal-coordinate vectors x,y, the squared geodesic interval contains

```text
S(x,y) = (x-y)^T eta (x-y) - C(x,y,x,y)/3
         - (nabla_e C)(x,y,x,y) (x+y)^e/12 + O(length^6).
```

The fifth-degree term follows the normal-coordinate expansion in
[Brewin, *Riemann Normal Coordinate expansions using Cadabra*, section 9 and Eq. 11.21c](https://arxiv.org/html/0903.2087).
The current joint fitter keeps only the quartic term. In vacuum the scalar
volume-clock correction begins at degree six, so it does not remove this
leading derivative term. This is an existing expansion used to diagnose
our implementation, not a new physical identity or priority claim.

For the tested rotating wave at the midpoint,
`H'=2Aw [[0,1],[1,0]]` and `nabla_u C_uiuj=-H'_ij`. The Brinkmann connection
vanishes on the central geodesic at that point. Transform the tensor and
derivative to the same tip-time-normalized frame used by the evaluator.
The resulting fifth-degree clock correction is symmetric in its endpoints.

An independent deterministic check compares exact plane-wave proper-time
intervals with the quartic and fifth-degree expressions. Shrinking the
dimensionless parameters as `A=.75*s^2`, `w=.5*s` gives:

| s | Maximum quartic error | Maximum error including gradient |
|---:|---:|---:|
| 1 | 6.37946e-5 | 1.35017e-5 |
| .5 | 7.67754e-6 | 8.82393e-7 |
| .25 | 9.37888e-7 | 5.63763e-8 |
| .125 | 1.15778e-7 | 3.55857e-9 |

The corrected errors decrease approximately sixteenfold per halving, as
expected after removing the leading derivative contribution in the
tip-normalized interval. This is a fixed-configuration expansion check.

## Oracle subtraction on the completed varying-profile pilot

`intrinsic_weyl_gradient_diagnostic.py`, committed at
`aa20b01fa6888d3b2cc214b2a8d356bba2ff4c0a`, subtracts the TRUE fifth-degree
correction from the same-landmark continuum volume clocks of the earlier
rotating pilot. It retains the original clock mask and covariance weights,
then refits the quartic model. True normal coordinates and A,w are used
here, so this is explicitly an oracle diagnostic and cannot satisfy the
order-only objective or overturn any previous failure.

| Previous case | A | w | Original continuum relative bias | After oracle subtraction |
|---|---:|---:|---:|---:|
| 4 | .25 | .5 | 2.7006% | 1.4306% |
| 5 | .75 | .5 | 4.7356% | 1.2598% |
| 6 | .25 | -.5 | 3.7357% | 1.4592% |
| 7 | .75 | -.5 | 2.8441% | 1.3914% |

The reduction supports this omitted term as a material cause of the bias.
It does not establish that estimating it from finite orders will reproduce
the oracle gain. Residual second-derivative/curvature-squared effects remain.
Full diagnostic output is `oracle_gradient_diagnostic.json`; the referenced
original summary hash is preserved in it.

## Free gradients do not erase local Weyl information in four checked cases

Let `q[e,k]` be the four derivatives of the ten Weyl basis coefficients.
The vacuum differential Bianchi identity

```text
nabla_a C_bcde + nabla_b C_cade + nabla_c C_abde = 0
```

gives a 24-row integer constraint matrix on 40 entries. Its rank is **16**,
leaving **24** gradient degrees of freedom. This count was checked both by
the numerical null space and independently by exact SymPy rank/null-space
calculation. This is the general vacuum derivative space, not a plane-wave
polarization or frequency ansatz.

`intrinsic_weyl_gradient_information.py`, committed at
`8bd7e54`, reads only the existing order archives for rotating-pilot cases
2,3,4,5 and reruns the pure quartic fitter. It does not import a generator
or open truth files. At those order-derived coordinates it appends all 24
gradient-response columns, whitens by the actual correlated clock covariance,
and projects the ten curvature columns off both position and gradient
nuisance columns. The projected-curvature threshold is fixed at
`1e-9 * norm(whitened_curvature_design,2)`.

| Existing case | Nuisance rank before -> after | Projected Weyl rank after | Largest tidal variance inflation |
|---|---|---:|---:|
| 2 | 253 -> 277 | 10 | 1.4841 |
| 3 | 253 -> 277 | 10 | 1.3090 |
| 4 | 253 -> 277 | 10 | 1.1352 |
| 5 | 253 -> 277 | 10 | 1.2799 |

The baseline covariance was independently reconstructed by the projected
information matrix and matched the original fitter to relative differences
below `5.5e-14`. Output, source and order hashes are in
`gradient_information.json`.

This is a local Jacobian calculation at gradient zero and the existing noisy
quartic fit. It does not solve the enlarged nonlinear inverse problem, check
its optimizer, certify coverage, bound omitted higher derivatives, or supply
a validated sample-size price. The observed 13.5--48.4% variance increase
is an encouraging feasibility estimate for these configurations only.

The next operational model must infer positions, C and the allowed gradient
coefficients together from order clocks. Its bias, finite-sample covariance,
coverage and powered accuracy must then be tested on fresh orders. Providing
the true derivative or imposing the known rotating-wave profile would not
meet that requirement.

```powershell
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_gradient_diagnostic.py --check
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_gradient_information.py --output outputs/new_information_diagnostic.json
```
