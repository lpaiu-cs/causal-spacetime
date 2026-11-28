# A stronger order-only test: failure of a reconstructed null time

2026-09-24. Analytic challenge and deterministic-check design, frozen before
new numerical evidence. This attacks the full-order resolution problem. The
proposed sufficient region improves its exponent from 8 to 4 (with a logarithmic
loss). It does **not** match the `N lambda -> 0` impossibility region.
Proof and implementation are by the same author and need independent audit.

## 1. Observation contract and target theorem

Keep the exact circular vacuum family and fixed-N physical-volume protocol
of [the probe](full_order_resolution_probe.md). For this construction assume
fixed known `0<q<=1`, clock/axial window length 1 and longitudinal slack h=1.
The restriction supplies an elementary certified inverse calibration; no
claim is made here for larger q. Phase and wavelength are not test inputs.

The target sufficient-side theorem is

\[
\boxed{N\lambda^4/\log N\longrightarrow\infty
\quad\Longrightarrow\quad
\operatorname{TV}(P_{N,\lambda},P_{N,0})\longrightarrow1.}
\tag{1}
\]

The witness is not a fixed collection of global motif counts. It estimates
three conditional quantities at a growing collection of observed events,
then checks a necessary monotonicity property along observed causal pairs.
The test below uses no matching background coordinates or true wave phases.

## 2. An exactly reconstructible background null coordinate

In the effective background put `s=u+1/2`, `k=sqrt(q)` and

\[
m(u)=\frac{\tan(ku)+\tan(k/2)}k,\quad M=2\tan(k/2)/k,
\quad t(u)=m(u)/M.
\]

The conformally related metric in coordinates `(m,v,x)` is Minkowski.
In particular **v is nondecreasing on every future causal relation**.

For a segment of u-length ell define

\[
K(\ell)=\sin(k\ell)/k,\quad
\alpha(\ell)=C_2(\ell)/K(\ell)
=\frac{1-k\ell\cot(k\ell)}{2q},\qquad\gamma=\alpha(1).
\]

For an interior event X let p and f be its past/future diamond volumes as
fractions of the full diamond. Write eta and zeta for its two null slacks.
The exact chain/volume hierarchy gives

\[
p=\alpha(s)\eta^2/\gamma,\quad
f=\alpha(1-s)\zeta^2/\gamma,\quad
c=P_2(I(p_{\rm tip},X))=p_0(qs^2),
\tag{2}
\]

where `p0(z)` is the known isotropic pair probability. Recover s from c, then
eta and zeta from p,f. Eliminating the transverse radius yields

\[
\boxed{V(p,f,c;q)=1-t+t\sqrt{p\gamma/\alpha(s)}
 -(1-t)\sqrt{f\gamma/\alpha(1-s)}=v.}
\tag{3}
\]

Thus `X precedes Y` requires `V(X)<=V(Y)`. This is an order-only population
constraint: all three inputs in (2) are probabilities of relations among
uniform events. An absolute geometric coordinate is used to *prove* the
identity, not supplied to its observer-side evaluation.

### Invertibility is checked rather than presumed

Put `C_r(z)=0F1(r+1/2;-z/4)/(2r-1)!` for a unit segment. Their derivatives
obey `C_r'=-r C_(r+1)` and they are positive for `0<=z<=1`. Therefore

\[
(\log p_0)'=-C_2/C_1-3C_4/C_3+4C_3/C_2.
\]

The source checks this expression with exact rational alternating-series
enclosures on all sixteen intervals `[j/16,(j+1)/16]`, not just at their
endpoints. Monotonicity of every C_r makes the interval upper bound
`-C_2(b)/C_1(a)-3C_4(b)/C_3(a)+4C_3(a)/C_2(b)` valid. The certificate requires
this upper bound to be below `-1/50` and `p0(1)>0.09`.
Consequently p0 has a strictly negative derivative bounded away from zero
throughout the range used here. Also alpha is strictly increasing: its
derivative has the sign of `x-sin(x)cos(x)>0` for `0<x<=1`.

## 3. A finite, randomized observable with a declared error bound

Randomly permute the unlabeled order independently of its content. Use half
the vertices as references, the rest as probes, and pair probes disjointly.
For each probe count its reference ancestors and descendants. Among its
reference ancestors count all comparable pairs, using transitive relations.
These give empirical versions of p,f,c. Fewer than two ancestors means skip.

Set the whole-test false-positive bound to `alpha_N=N^-2`. Choose integer B
with `2^B>=6N^3`. A sample mean based on n references gets radius
`sqrt(B/(2n))`; a relation fraction among n ancestors gets radius
`sqrt(B/(2 floor(n/2)))`. Use outward enclosures for these square roots.
The usual bounded-variable and order-two U-statistic Hoeffding inequalities
give two-sided failure at most `2 exp(-B)<=2^(1-B)` per quantity. Conditional
on a probe and its ancestor count, its reference ancestors are independent
uniform events in that subdiamond. Thus a union bound over the at most N
probes and three quantities gives simultaneous coverage at least `1-N^-2`.
Shared references do not invalidate this union bound.

Push each triple of probability intervals through (2)-(3) by interval
arithmetic, retaining every background-compatible s. Skip noninformative
boxes touching s=0 or s=1. Reject only if an observed causal pair X<Y has
`upper V(Y)<lower V(X)`. For the reverse relation reverse the comparison.
Under the background all exact v values are covered on the simultaneous
event, so the false-positive probability is at most `N^-2`.
This is a **detection** rule; no accuracy gate is replaced by a point estimate.

The executable construction uses exact rational alternating-series bounds,
dyadic outward square roots and interval arithmetic for `0<q<=1`. Its
calibration routines take counts/probabilities and q, not embedding data.
No Gaussian variance estimate is required. The large finite-N cost constants
are not evaluated by the deterministic check below; no finite 90%-power
sampling campaign is claimed.

## 4. Why the wave violates this coordinate at order lambda squared

Use the transverse rescaling from the probe, so its matrix primitive is

\[
\widetilde M(a,u)=m(u)I+E_\lambda(u),\qquad
E_\lambda(u)=-(2\lambda\sqrt q/c)J(a,u),\quad
c=\sqrt{1+4q\lambda^2}.
\]

For interior u, `E_lambda=O(lambda^2)` is trace free. Let

\[
F_\lambda(u)=\frac{E_\lambda(u)-t(u)E_\lambda(b)}{m(u)[M-m(u)]}.
\tag{4}
\]

Map each wave slice to the background slice while preserving its two
slacks. Explicitly, with
`Q_lambda=Mtilde(a,u)^-1+Mtilde(u,b)^-1` and
`W0=m(M-m)/M`, define `z=(W0 Q_lambda)^(1/2)y`. Set
`v_tilde=eta+|z|^2/(2m)`. This is a coordinate device for the proof, not an
input to the observable. On every compact interior patch, matrix expansion
gives

\[
v_{\rm wave}=\widetilde v-\tfrac12 z^TF_\lambda(u)z+O(\lambda^4),
\qquad y=z+O(\lambda^2).
\tag{5}
\]

Take two u values in separated compact interior intervals, with the same
z and the same v_tilde. Their exact wave causal margin is

\[
v_2-v_1-\tfrac12(y_2-y_1)^T\widetilde M(u_1,u_2)^{-1}(y_2-y_1)
=\tfrac12 z^T[F_\lambda(u_1)-F_\lambda(u_2)]z+O(\lambda^4).
\tag{6}
\]

The matrix on the right is trace free. Whenever its norm exceeds
`c1 lambda^2`, a set of transverse angles of fixed positive measure makes
the margin positive by at least `c2 lambda^2`. Decrease v_tilde of the second
point by a smaller positive multiple of lambda^2. The points remain wave
causally related but their background reconstructed null time decreases.

Such u pairs occupy a uniformly positive fraction, not an isolated tuned
sequence. Integration by parts gives
`E_lambda/lambda^2` as a nonzero bounded-amplitude rotating matrix at phase
`u/lambda`, plus smooth coefficients multiplying the endpoint phase and a
uniform `O(lambda)` remainder. On the two separated u intervals, integration
of the squared norm of the difference in (4) has a positive lower limit:
the local oscillatory variances add and the mixed oscillatory integrals
vanish. Uniform boundedness then gives a fixed positive-measure set where
the norm exceeds a fixed threshold. This argument is uniform in global phase.

Allow the second transverse vector to differ from the first by at most
`c3 lambda`, and its v_tilde to occupy an interval of width `c4 lambda^2`
inside the negative-time gap. The wave action costs `O(lambda^2)` there and
the constants can be reduced to retain (6). All other coordinates occupy
fixed interior sets with physical-volume density bounded above and below.
Thus an independent pair has probability at least `c5 lambda^4` to be a
causal pair whose v_tilde decreases by at least `c6 lambda^2`.

## 5. Removing the coordinate oracle and obtaining the sufficient region

In the wave, the same exact volume hierarchy gives p,f from its actual
slacks and perturbed scalar convolutions. On compact subdiamonds the
[fixed-q chain expansion](full_order_chain_asymptotics.md) is uniform.
Its constant kernel factor cancels from volume ratios. Therefore the wave's
three population inputs to (3) differ from their background expressions at
the same `(u,eta,zeta)` by `O(lambda^4 log(1/lambda))`.
The inverse calibration is smooth on the chosen compact set, so

\[
V_{\rm wave\ population}(X)=\widetilde v(X)
                   +O(\lambda^4\log(1/\lambda)).
\tag{7}
\]

This error is smaller than the lambda-squared violation. Ancestor and
descendant volume fractions are bounded away from zero on that set.
Standard binomial concentration makes their observed reference counts
proportional to N with probability tending to one. The simultaneous
probability boxes, and therefore their V images there, have widths
`O(sqrt(log N/N))`. Under (1) these widths are `o(lambda^2)`.

There are order N disjoint independent probe pairs. The probability that
none lies in the violating set is at most `exp(-c N lambda^4)`.
Combining this with interval coverage and (7) gives false negatives tending
to zero under (1), while false positives are bounded by `N^-2`.
This proves the proposed sufficient side *if the preceding geometric and
uniform-remainder arguments withstand scrutiny*. It is not a sharp lower
bound. In particular no argument yet excludes an even better full-order test.

## 6. Frozen deterministic attack on the delicate step

Source: `experiments/positive_control/order_wave_null_coordinate.py`.
Freeze this design and source before generating new rows. Use CPU only.

- Exact q values `{0.25,1}` and inverse wavelengths
  `{64,128,256,512,1024,2048}`, phase 0.37.
- At each model, generator-side u candidates occupy one full phase cycle
  around 0 and 1/4, with 16 equally spaced phases in each interval. Fix
  `|z|=0.2` and initial `v_tilde=0.5`. Use the top eigenvector of the matrix
  in (6), and retain the candidate with largest *exact* causal margin.
  This is a predeclared construction of a counterexample, not a blind
  classification score or an estimate of its probability.
- Lower the second v_tilde by half that margin. Evaluate exact past/future
  volume ratios and past pair probabilities, then pass only those three
  moments and q to the observer-side inverse. Record its recovered V gap,
  the exact wave causal margin, and each divided by lambda squared.
- Use `max(96,ceil(r*ell/2)+64)` Gauss nodes for each subdiamond. At q=1,
  r=2048 retain the same chosen candidate, double quadrature order, and
  tighten the primitive ODE. Record every discrepancy; do not retune the
  construction or remove failed lower-frequency rows.
- Controls certify monotonicity with rational bounds, check background
  reconstruction/interval containment and observer rejection logic, and
  compare a wave matrix determinant with the independent existing kernel.

The decisive numerical question is whether the **moment-reconstructed** V
decrease survives, not merely whether the generator's canonical coordinate
decreases. Failure at every increasing rung would directly challenge (7)
or the claimed signal size. Success supplies a deterministic consistency
check only. No sample-complexity exponent, minimax law, novelty claim, or
finite-sample power is inferred from these twelve constructed pairs.
