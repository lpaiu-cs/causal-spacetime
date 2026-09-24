# Using shared-reference cancellation to improve the sufficient region

2026-09-24. Analytic construction and deterministic-check design. Freeze
before running new geometry cases. The target remains the sharp full-order
law and an attaining observable; this step does not complete that target.
All derivations await independent researcher audit.

For the same fixed known `0<q<=1` circular vacuum family and fixed-N iid
physical-volume protocol, the proposed improvement is

\[
\boxed{N\lambda^3/\log N\longrightarrow\infty
\quad\Longrightarrow\quad
\operatorname{TV}(P_{N,\lambda},P_{N,0})\longrightarrow1.}
\tag{1}
\]

The previous sufficient condition used exponent 4. The difference here is
to bound the **difference** between two reconstructed coordinates with their
shared sampling noise retained, and inspect all probe pairs. No wavelength,
phase, embedding, or generator-selected pair is supplied to the observer.

## 1. Use unconditional moments with nested kernels

Split an independently permuted order into M=floor(N/2) references and the
remaining probes. For each probe X estimate

\[
p_X=\Pr(Z<X),\quad f_X=\Pr(X<Z),\quad
a_X=\Pr(Z_1<X,Z_2<X, Z_1\text{ and }Z_2\text{ comparable}).
\]

The empirical a is the number of comparable reference-ancestor pairs divided
by `choose(M,2)`. Its denominator is **not** the number of ancestor pairs.
At population level `a_X=p_X^2 c_X`, where c is the past-subdiamond pair
probability used by the previous null-coordinate construction. Thus set

\[
\mathcal F(p,f,a)=V(p,f,a/p^2;q).
\tag{2}
\]

Under the background, F equals the exact null coordinate v.
For an observed causal probe pair X<Y, the three increments
`p_Y-p_X`, `f_X-f_Y`, and `a_Y-a_X` have nonnegative kernels bounded by 1:
past sets are nested, and future sets are nested in the reverse direction.
Their kernel variances are therefore bounded by their respective means.
The same references must be used at X and Y; treating the counts as two
independent coordinate estimates would lose the improvement.

## 2. Simultaneous difference confidence bounds

Choose integer B with `2^B>=32N^4`. For a mean with effective sample size k,
the initial Hoeffding radius is `r_H=sqrt(B/(2k))`. Here k=M for p,f, and
k=floor(M/2) for a. For each nonnegative pair increment with empirical mean
d_hat, use the upper bound `d_U=min(1,d_hat+r_H)` and radius

\[
r_D=\min\left(r_H,
\sqrt{2d_U B/k}+2B/(3k)\right).
\tag{3}
\]

The second term is Bernstein's bound with variance at most the true increment
mean, on the initial Hoeffding coverage event. For an order-two U-statistic,
average over all permutations of disjoint reference pairs. Jensen's inequality
for the exponential transfers the iid-pair Bernstein bound to the complete
U-statistic with k=floor(M/2). This does not assume that overlapping reference
pairs are independent. See the iid U-statistic version in
[Theorem 1.6 of the 2022 concentration paper](https://www.sciencedirect.com/science/article/abs/pii/S0167715222000992).

The implementation also applies (3) to each individual moment, whose
nonnegative kernel has variance at most its mean. There are at most 3N
individual Hoeffding and 3N individual Bernstein events, 3N^2 increment
Hoeffding events and 3N^2 increment Bernstein events. Each has two-sided
failure at most `2 exp(-B)<=2^(1-B)`. A union bound gives simultaneous
coverage with probability at least `1-N^-2`. Dependence among these events,
reuse of references, and subsequent selection of a pair do not invalidate it.

For a pair put `b=theta_X`, `d=theta_Y-theta_X`, where theta=(p,f,a), so the
f component of d is negative. Let B_b be the individual-moment box for b,
and B_d the signed increment box from (3). Use the rectangular domain
`R=B_b+hull(0,B_d)`, which contains all the straight paths used next.
Skip the pair if the calibration is not defined throughout R. Obtain bounds
L on every component of gradient F, and H on every entry of Hessian F,
valid on R. Let e_b bound `|b_hat-b|_1`, e_d bound `|d_hat-d|_1`, and
D bound `|d|_1`, from those boxes. The deterministic error bound is

\[
|[\mathcal F(\widehat b+\widehat d)-\mathcal F(\widehat b)]
 -[\mathcal F(b+d)-\mathcal F(b)]|
\le L e_d+H D e_b.
\tag{4}
\]

First replace d_hat by d at fixed b_hat, using L. Then differentiate
`F(b+d)-F(b)` with respect to b; the Hessian bound gives H|d|_1. The domain
R contains b_hat+d as well as the actual and estimated endpoints, so no
uncovered mixed endpoint enters this argument.

The implementation retains the componentwise version of (4):
`sum_i L_i e_d_i + sum_(i,j) H_ij e_b_i D_j`, with bounds on each derivative
entry. This avoids multiplying every error by the largest derivative, while
the coarser displayed inequality suffices for the rate proof.

Evaluate the central difference with certified intervals and add the radius
in (4). Reject if its upper endpoint is negative for any observed X<Y.
Under H0, `F(theta_Y)-F(theta_X)=v_Y-v_X>=0`; the whole-test false-positive
probability is at most N^-2. This is a detection guarantee, not an accuracy
gate decided by a point estimate.

### Executable derivative certification

The new source reuses the rational intervals, alternating-series calibration,
and inverse monotonicity certificate. For c=a/p^2, write
`z=p0^-1(c)` and `s=sqrt(z/q)`. Then

\[
\mathcal F=1-t(s)+A(s)\sqrt p-A(1-s)\sqrt f,
\qquad A(s)=t(s)\sqrt{\gamma/\alpha(s)}.
\]

The source differentiates this expression twice, including the implicit
inverse derivatives `z_c=1/p0'(z)` and `z_cc=-p0''(z)/p0'(z)^3`.
All values and derivatives are enclosed on R, not sampled at its centre.
Any box touching p=0, f=0, or an inverse-calibration endpoint is skipped.
On compact interior population sets these derivatives are uniformly bounded
for fixed q, and the empirical boxes eventually remain in such a set.
The constants may be large; no finite-N 90%-power price is claimed.

## 3. Short null-neighbour pairs have a larger signal-to-noise ratio

Use the canonical slice map and F_lambda matrix from
[the previous construction](null_coordinate_challenge.md), Section 4.
Take u in a fixed compact interval and `u'=u+tau lambda`, with tau in a
fixed subinterval of `(pi/2,3pi/2)`. Integration by parts in the primitive
now gives, in complex trace-free notation,

\[
\frac{F_\lambda(u')-F_\lambda(u)}{\lambda^2}
=\frac{2i\sqrt q f(u)}{m(u)[M-m(u)]}
e^{i(u/\lambda+\phi)}(e^{i\tau}-1)+O(\lambda).
\tag{5}
\]

Smooth coefficients multiplying the endpoint phase change by O(lambda).
Thus the leading matrix has norm bounded below uniformly on the chosen
u,tau sets. A fixed fraction of transverse directions makes the relevant
quadratic form positive by `c lambda^2`.

For equal canonical transverse z and equal reconstructed v_tilde, the actual
transverse difference is O(lambda^2). The connecting matrix inverse is now
O(1/lambda), so its action is O(lambda^3), still smaller than the leading
lambda-squared positive margin. This replaces the O(lambda^4) action term
for widely separated u endpoints in the earlier proof.

Allow the second canonical transverse position in a ball of radius
`c' lambda^(3/2)` about the first, and decrease its v_tilde by a small
multiple of lambda^2. The transverse action remains O(lambda^2), with a
coefficient that can be made smaller than the positive margin. An open
interval of negative v_tilde gaps of width `c'' lambda^2` is available.
All other variables stay in compact interior sets of positive density.
The resulting independent-pair mass is at least

\[
c\underbrace{\lambda}_{\text{u separation}}
 \underbrace{\lambda^3}_{\text{two transverse coordinates}}
 \underbrace{\lambda^2}_{\text{null-time gap}}
=c\lambda^6.
\tag{6}
\]

The conditional mass given either endpoint is at most C lambda^6, by the
same coordinate support sizes and bounded density. For the number W of
such pairs among order N probes,

\[
\mathbb EW\ge cN^2\lambda^6,\qquad
\operatorname{Var}W\le C(N^2\lambda^6+N^3\lambda^{12}).
\tag{7}
\]

Consequently `N lambda^3 -> infinity` makes `Pr(W=0)->0`. Using only N
disjoint pairs would miss this gain. The actual observer checks every
causal probe pair and does not receive the generator's witness selection.

## 4. Remove the remaining population coordinates

The prior uniform compact-set expansion still gives
`F(theta_X)=v_tilde(X)+O(lambda^4 log(1/lambda))`. It applies to each of
these endpoints: their past and future subdiamonds have fixed nonzero
longitudinal sizes, despite the short separation between the endpoints.
The population-coordinate violation is therefore still `c lambda^2`.

For these pairs, `|theta_Y-theta_X|=O(lambda)` uniformly. One way to check
this is to differentiate the exact prefix/suffix volume and chain formulas:
their kernels and first endpoint derivatives are uniformly bounded; the
endpoint-to-tip segments stay bounded below, and the endpoint coordinates
differ by O(lambda). The rapidly growing second curvature derivatives are
not used. Nestedness then bounds each increment kernel's variance by
its O(lambda) mean.

Individual moment errors are O(sqrt(log N/N)). The initial increment
Hoeffding error is o(lambda) under (1), so d_U=O(lambda). Formula (3)
and then (4) give coordinate-difference width

\[
O_q\left(\sqrt{\frac{\lambda\log N}{N}}
          +\frac{\log N}{N}
          +\lambda\sqrt{\frac{\log N}{N}}\right)
=o(\lambda^2)
\tag{8}
\]

under (1). Together with (6)-(7), simultaneous coverage and the negligible
population bias, this proves the proposed sufficient condition. The test
needs no frequency-dependent threshold or phase estimate.

## 5. Frozen deterministic challenge and limits

Implementation: `experiments/positive_control/order_wave_local_null.py`.
Reuse the exact WaveFixture, its scalar hierarchy, and the old independent
five-state Rosen ODE cross-check; do not fork a second geometry generator.

- Freeze q in `{0.25,1}`, r in `{64,128,256,512,1024,2048}`, phase 0.37.
  At each rung set `u1=1/5`, `u2=1/5+pi/r`, canonical v=0.5 and transverse
  radius 0.2. Use the positive eigenvector of `F(u1)-F(u2)`. This is a
  fixed construction, with no phase/u search or reselection.
- Construct the canonical pair, record its exact causal margin, and lower
  the second v_tilde by half that margin when it is positive. Retain failures.
  Pass only the population (p,f,c) values to the inverse calibration, and
  record the recovered decrease divided by lambda^2. Record the three
  unconditional moment increments divided by lambda as a variance control.
- At q=1,r=2048 retain the exact same pair, double quadrature order and
  tighten the primitive. Then check that pair again using the existing
  independent five-state Rosen kernel and direct action quadrature. Freeze
  the independent-formulation criterion as positive recovered decrease and
  less than 1% relative change, as in the previous witness check.
- Rational controls test first and second derivative enclosures against
  independently evaluated finite differences of the calibrated value,
  using whole perturbation boxes. Check the signed difference interval on
  nearby background population moments, and observer input/count invariants.
- At q=1,r=2048 also feed rounded population moments, with an explicitly
  synthetic reference count `10^30`, into the certified difference routine.
  Record its interval without treating this as an observed order or a power
  result. This checks the algebraic decision path at nearly exact moments;
  it supplies no evidence that such a sample was collected or is affordable.

These are constructed population witnesses and numerical consistency checks.
They are not random-sprinkling classification accuracy, a finite-N powered
campaign, or an independent researcher proof audit. The implementation has
cubic reference-chain counting and potentially many rational pair bounds;
no efficient large-N runtime is established.

The new exponent 3, if this argument survives audit, still does not match
the [necessary powers below 4/3](quantized_order_bound.md). No claim is made
at either critical endpoint, for unknown q, or for arbitrary wave families.
