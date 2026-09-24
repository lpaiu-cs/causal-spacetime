# A constructive improvement of the sufficient resolution region

2026-09-24. The sharp full-order boundary is **still unresolved**. This step
supplies a more informative observable, not just another global chain count:
reconstruct a background null coordinate from conditional order statistics,
and reject the background if it decreases along an observed causal pair.

For fixed known `0<q<=1` in the exact circular family, the
[construction and proof](null_coordinate_challenge.md) give the sufficient
region

\[
N\lambda^4/\log N\longrightarrow\infty.
\]

The former pair-count sufficient region was
`N lambda^8 log^2(1/lambda) -> infinity`. The new argument uses the arrangement
of conditional counts across events; it therefore avoids the earlier
fixed-global-chain obstruction. The full-order impossibility statement
`N lambda -> 0` is unchanged. The exponents 1 and 4 do not match, and 4 is
**not** claimed optimal. Both mathematical arguments await independent audit.

## 1. Frozen evidence and the cases that failed

Design and implementation were frozen at
`ab16d6a3f96a5a462e976e6efb4b3e9059fc336e`. A clean-tree run completed 12
constructed pairs and one refinement, using 13 geometry ODE solves in
24.65 seconds on CPU. Raw SHA256:
`0924f97d4bb3a44dd0d7c1f405605d006f807e7f48ce307058d10e9a75a11357`.
The [manifest](null_coordinate_evidence/manifest.json),
[raw rows](null_coordinate_evidence/rows.jsonl), and
[summary](null_coordinate_evidence/summary.json) are retained together.

For every constructed pair the exact wave causal margin was positive.
But a generator-coordinate violation by itself would be insufficient.
The observer-side reconstruction received only the exact population past
fraction, future fraction, and past-subdiamond pair probability, plus q.
Positive entries below mean that its reconstructed coordinate decreased
despite the pair being causally related. Negative entries are failures of
this construction at that rung, retained without changing the design.

| r=1/lambda | q=0.25: reconstructed decrease / lambda^2 | q=1: reconstructed decrease / lambda^2 |
|---:|---:|---:|
| 64 | -3.49210 | -3.36643 |
| 128 | -1.32457 | -1.21805 |
| 256 | -0.32420 | -0.23373 |
| 512 | -0.03450 | 0.06063 |
| 1024 | 0.06138 | 0.14222 |
| 2048 | 0.09071 | 0.17118 |

Thus 5 of 12 constructed population witnesses succeed. This is **not** a
classifier accuracy or a success probability for random sprinklings.
The lower-frequency failures show that the reconstruction bias cannot be
discarded in a finite-frequency interpretation. At the largest rung, the
q=1 reconstructed decrease is `4.08123e-8`, while its actual wave causal
margin is `4.37452e-8`. The refinement changed the saved moments and decrease
by zero at stored double precision. This is consistency, not a certified
ODE/quadrature enclosure or an independent replication.

The exact-rational inverse-calibration control bounded
`d log(p0)/dz <= -0.03555319047 < -1/50` across all of `[0,1]`, using
intervals rather than pointwise sampling. Its background-coordinate control
enclosed the exact background v in a box of width `1.14516e-10` when supplied
moment intervals of radius `1e-13`.

An audit clarified the kernel control: the pre-run scalar correction and
matrix formula shared the complex-primitive implementation. A **post-run**
control was therefore added against the original independent five-state
Rosen ODE; it passed the `2e-12` absolute comparison. It does not retrospectively
make the original shared-primitive check independent.

The largest q=1, r=2048 witness was then recomputed at `7904c10` using
the five-state Rosen ODE for all scalar kernels and direct Gauss integration
for its action matrices. This replaces both the primitive and stable-change
routes used by the primary run. The [independent-formulation check](null_coordinate_evidence/independent_formulation_check.json)
found decrease `4.08122519113e-8`, differing from the primary value by
`3.42504e-14` (relative `8.39217e-7`). It passed the pre-execution numerical
check of positive sign and less than 1% relative change. This is a distinct
mathematical formulation implemented by the same author, not an independent
researcher audit or a rigorous numerical enclosure.

## 2. Why the violating set has positive measure

This expands the measure argument in the design; the numerical pair search
does not prove it. Represent a trace-free symmetric matrix by a complex
number. Write `D(u)=m(u)[M-m(u)]`, `a=-1/2`, `b=1/2`, and
`f(u)=cos(sqrt(q)u)^-2`. On two separated compact interior intervals I1,I2,
integration by parts in the exact primitive gives

\[
\frac{F_\lambda(u)}{\lambda^2}
=\frac{2i\sqrt q}{D(u)}\left[
 f(u)e^{iu/\lambda}-(1-t(u))f(a)e^{ia/\lambda}
 -t(u)f(b)e^{ib/\lambda}\right]+O(\lambda).
\tag{1}
\]

A global phase only rotates this expression. The local oscillatory terms in
the squared difference have integral

\[
C_0=4q\left[|I_2|\int_{I_1}(f/D)^2du
                 +|I_1|\int_{I_2}(f/D)^2du\right]>0.
\]

The squared endpoint contribution is nonnegative. Mixed terms between a
local phase and smooth endpoint coefficients are `O(lambda)` by integration
by parts; the two-local-phase cross term is `O(lambda^2)`. Hence the squared
norm integral is at least `C0/2` for sufficiently small lambda, uniformly
in endpoint phases. The same expression has a uniform upper bound.
An elementary second-moment bound then gives a fixed positive measure of
u pairs where `||F_lambda(u1)-F_lambda(u2)||>=c lambda^2`.

For those pairs, a fixed fraction of transverse directions has a positive
quadratic form. Keep the first transverse point in an interior annulus and
allow the second within a ball of radius `c' lambda` about it. This costs a
factor `lambda^2` in measure and only `O(lambda^2)` in the causal action.
Put the second reconstructed null time in a negative interval of width
`c'' lambda^2` relative to the first. Taking the constants small enough
leaves the wave causal margin positive. This costs another `lambda^2`.
The slice transformation has uniformly positive bounded Jacobian on this
compact set. Therefore an independent physical-volume pair has probability
at least `c''' lambda^4` of such a violation.

This is why order N disjoint probe pairs suffice in the proof. Their failure
probability is at most `exp(-c N lambda^4)`. The additional `log N` in the
sufficient condition comes from simultaneous uncertainty bounds when the
null coordinate is reconstructed from reference events.

## 3. What was constructed, and what was not experimentally shown

The source contains an executable `order_test(order,q,rng)`. Its inputs are
a boolean transitive order matrix, a declared exact q calibration, and an
independent random generator. It randomly partitions the events, estimates
three conditional probabilities, uses rational interval calibration, and
checks coordinate monotonicity along disjoint probe pairs. It has no
wavelength, phase, or coordinate argument. Full transitive comparabilities
must be supplied, not merely an edge list lacking its closure.

Under the specified iid physical-volume protocol, the analytical
false-positive bound is `N^-2`. Under the alternative the proof uses uniform
compact-set chain asymptotics to show a coordinate bias
`O(lambda^4 log(1/lambda))`, smaller than the violation of order lambda^2.
The reference-sampling uncertainty is `O(sqrt(log N/N))`. Their comparison,
together with the positive-measure argument above, supplies the sufficient
condition. None of these error probabilities is estimated from the twelve
constructed pairs.

The direct implementation has cubic worst-case reference-chain counting
and exact rational interval overhead. It is an executable construction,
not a demonstrated practical large-N algorithm. No powered finite-poset
campaign has been run, and the potentially large constants have not been
sized for finite 90% power. Unknown q, q>1, other wave classes, and caustics
are not covered by this new test's proof. The exact supplied q is not an
estimated parameter with an ignored measurement error.

## 4. A separation that now changes the research decision

Set `r=1/lambda` and, for example, `N=r^6`. The new test's sufficient ratio
is `r^2/log(r^6)`, which diverges. Its two error probabilities tend to zero
under the argument above. But any fixed list of global chain fractions has
normalized mean shift of order
`sqrt(N) lambda^4 log(1/lambda)=log(r)/r`, which vanishes.
The [positive-definite projection result](full_order_chain_results.md)
therefore gives no asymptotic power beyond size to regular Gaussian/Wald
tests of that fixed vector, while the new conditional-order test separates.

This is a genuine difference in what the two observation procedures retain
within the stated model, rather than a better fitted constant. The comparison
is restricted to the regular fixed-chain tests covered by that result;
weak convergence does not rule out every irregular test of their exact
discrete distributions. Nor does it establish a minimax boundary for the
full order. The remaining goal is to close the exponent gap, audit the proof,
and identify an observable attaining the actual boundary. Novelty has not
been established by a complete literature review.

Verification of the portable evidence:

```powershell
E:/lab/causal-spacetime/.venv/Scripts/python.exe experiments/positive_control/order_wave_null_coordinate.py --self-check --summarize docs/paper/paper_a/order_wave_energy/null_coordinate_evidence
```

Subsequent update, 2026-09-24: the [quantized-coordinate lower-bound proof](quantized_order_bound.md)
improves the necessity side quoted above from exponent 1 to every fixed power
below 4/3. It leaves this instrument and its sufficient condition unchanged.
Neither the endpoint at 4/3 nor a matching upper/lower boundary is established.
