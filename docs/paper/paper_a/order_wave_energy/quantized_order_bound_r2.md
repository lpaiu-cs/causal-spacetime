# A stronger impossibility region after hiding exact coordinates, revision 2: constants linked and the smoothing lemma written out

2026-11-26. **Revision of [quantized_order_bound.md](quantized_order_bound.md).**
- **The original is unchanged.** Its Sections 1 to 5 (the proof) were
  frozen at `5dda03d` and are byte-identical to that freeze. The later
  commits `ac2cd7d` and `9205b3a` added its Section 7, edited one preamble
  sentence, and replaced the bracket paragraph of its Section 6. SHA-256 of
  the original as it stands:
  `39217a16c9770ff9b56b326984262b1424f4aaa4e1a2bd934238d658bed0564b`.
- **What this revision does.** It carries the independent audit's findings
  for argument 1 into the text ([AUDIT.md](independent_audit/AUDIT.md)
  Section 3, argument 1, items 1 to 6; referee report
  [proof1_lower_bound/report.md](independent_audit/proof1_lower_bound/report.md),
  findings F1 to F9):
  - the constants are linked, so that the non-edge claim holds;
  - the hypothesis of (3) is stated;
  - the common orientation of an edge is argued;
  - the smoothing lemma is written as a lemma with its constant and proof;
  - the Hellinger convention is fixed;
  - the bound is shown to hold for every `N >= 1` and `0 < lambda <= 1`;
  - `K = 1` is included;
  - the strict-separation sentence names its unaudited dependency.

  Section 0 maps every change. The remaining text is the original's.
- **Authorship and review.** AI-written. The original derivation is by one
  AI author. It was audited by an AI referee, whose re-derivations this
  revision follows and cites. A second AI referee, in a fresh context,
  reviewed this revision:
  [independent_audit/argument1_revision/review.md](independent_audit/argument1_revision/review.md).
  Section 9 lists the corrections made after it.
- **What it is not.** It is not a human-reviewed or published theorem. It
  is asymptotic only. It covers one family: the exact circular vacuum plane
  wave with fixed `0 < q < pi^2` and fixed-N i.i.d. physical-volume
  sampling.
- **Status.** The original's preamble said it "has not received independent
  researcher review". Since then the AI-only audit graded it CORRECT, with
  minor findings only (AUDIT.md Section 2). Section 8 lists what is now
  written in full and what remains imported or cited.

For every fixed `0<q<pi^2` in the exact circular family, and every fixed
integer `K >= 1`, the result is

\[
\boxed{N^{3K+1}\lambda^{4K}\longrightarrow0
\quad\Longrightarrow\quad
\operatorname{TV}(P_{N,\lambda},P_{N,0})\longrightarrow0.}
\tag{1}
\]

Quantitatively (Section 5): there is `C_K(q)` such that for **all**
`N >= 1` and `0 < lambda <= 1`, uniformly in the wave's phase,

\[
\operatorname{TV}(P_{N,\lambda},P_{N,0})
\le C_K\left[N\lambda^2+(N\lambda^4)^{1/(2K+1)}
 +(N^{3K+1}\lambda^{4K})^{1/(2K+1)}\right].
\tag{8}
\]

`K = 1` gives `N lambda -> 0`. `K = 2` gives `N lambda^(8/7) -> 0`.
Letting K be a sufficiently large **fixed** integer proves
indistinguishability at every fixed power `N=O(lambda^(-alpha))` with
`alpha<4/3`. This is **not** the statement `N lambda^(4/3) -> 0`: constants
below depend on K, and no growing-K estimate or logarithmic critical
boundary is claimed. Nothing is proved at `alpha = 4/3` itself.

The [coordinate-oracle result](coordinate_oracle_barrier.md) has a sharp
boundary at `N lambda~1`. So for any fixed `1<alpha<4/3`, at
`N~lambda^(-alpha)`, exact coordinates plus order distinguish consistently,
while the full finite order alone is asymptotically indistinguishable.
**This separation of observation protocols is conditional.** Its oracle
side rests on `kappa(q) > 0` and the other statements of
`coordinate_oracle_barrier.md`, which the audit did not audit (Section 8).
It is not evidence that the order-only sufficient exponent is sharp. The
later [shared-reference construction](local_null_challenge.md) improves
that exponent from 4 to 3, still leaving a gap.

## 0. What changed against the frozen original

| where | change | audit item / finding |
|---|---|---|
| preamble | replaced. The original's status sentence is summarised under "Status". The statement is extended to `K >= 1` and given in the quantitative form (8) for all `N, lambda`. The strict-separation sentence names its dependency on `coordinate_oracle_barrier.md` | items 5, 6; F5, F8, F9 |
| Section 1 | `rho_min`, `rho_max`, `rho'_max` and `f_max` defined; the phase-uniformity stated; the one-point constant renamed `C_1` and its source cited | notation |
| Section 2 | rewritten: Lemma 2.1 (operator bound with `C_E`), Lemma 2.2 (the sandwich (4) with `c0 = 2 sqrt 2 C_E`, proved for every ordered pair), Corollary 2.3 ((3) with its hypothesis, the crossing set), Lemma 2.4 (`T` is a common time function); `b` is the signed `Delta m` | items 1, 2; F1, F3 |
| Section 3 | `L` fixed with `L > 1 + sqrt 3`, and `h >= h_min = c0 delta / (L - 1 - sqrt 3)`; the non-edge argument written, with the explicit `(1 + sqrt 3) h` shift replacing "the null-margin functions are uniformly Lipschitz"; the orientation argument written via Lemma 2.4 | items 1, 2; F2, F3 |
| Section 3 | the cubes stated to be half-open on one aligned grid, so that T-intervals coincide or are disjoint (the orientation argument needs this) | item 2; F3 |
| Section 3 | the discard rule changed from "time orientation not fixed throughout its two cubes" to "every edge between cubes with the same T-interval", which contains it | item 2; F3 |
| Section 3 | the boundary-layer reason rewritten (cubes within `2h` of the boundary); the conditional density bounds made explicit; the uniform bound and the tree count written with their constants | item 1 |
| Section 4 | rewritten: Lemma 4.1 (Lipschitz gradient, with the constant `Lambda`), Lemma 4.2 (`sqrt p` Lipschitz across zeros), Lemma 4.3 (component Hellinger bound). The lemma's Lipschitz constant is renamed `Lambda` so that it does not clash with the edge constant `L` | item 3; F4 |
| Section 5 | the Hellinger convention stated; "components with an edge <= uncertain pairs" stated; (7) with its range of `h`; (8) shown for every `N >= 1`, `0 < lambda <= 1`; the chosen `h` is not the minimiser; `K = 1` | items 4, 5; F5, F6, F9 |
| Section 5 | the original's "Under (1), h tends to zero ... This proves (1)" replaced by the paragraphs "(1) follows from (8)" and "(1) is exactly what (7) delivers"; a new "Exponents" paragraph with the formula for K and a table of examples | item 5; F5, F9 |
| Sections 6, 7 | the original's text, unchanged apart from an italic note at their head | |
| Section 8 | new: status, scope, constants, and what remains | F7, F8 |
| Section 9 | new: the review, and the corrections made after it | review |
| (after review) preamble, Sections 0, 2, 4, 5, 6, 8 | the corrections of Section 9 | review |

Deterministic checks of the new constants:
[independent_audit/argument1_revision/checks.py](independent_audit/argument1_revision/checks.py).

## 1. Put the background in a fixed Minkowski chart

Use the y chart of the coordinate-oracle derivation, set
`m(u)=integral_(-1/2)^u a^-2`, and define
`T=(m+v)/sqrt(2)`, `Z=(m-v)/sqrt(2)`, `z=(Z,y1,y2)`.
The background relation is exactly

\[
X\prec_0Y\quad\Longleftrightarrow\quad
T_Y-T_X\ge |z_Y-z_X|.
\tag{2}
\]

The reference sampling domain D is a fixed bounded Minkowski diamond.
Its density in this chart is proportional to `a(u(m))^4`: on a neighbourhood
of the closed diamond it is positive, bounded, and continuously differentiable
with bounded derivative. Write `rho_min`, `rho_max` and `rho'_max` for the
infimum and supremum of the normalised density and the supremum of
`|d rho / d T|` there. With `a = cos(sqrt(q) u)` and `|u| <= 1/2`, write
`f = a^-2` and `f_max = sec^2(sqrt(q)/2)`, finite exactly because
`q < pi^2`. Constants can depend on the fixed q and positive patch margin.
They are not uniform as q approaches pi^2. Every constant below is uniform
in the wave's phase.

Write `delta=lambda^2`. The normalized wave and reference one-point laws
differ by at most `C_1 delta` in TV, by the existing domain-strip bound
([coordinate_oracle_barrier.md](coordinate_oracle_barrier.md) Section 1;
re-derived with an explicit constant in the referee report, Section 3.3).
Replacing the wave sampling locations by N iid reference locations therefore
costs at most `C_1 N delta` after applying any order map. For the remaining
argument **both models use exactly the same reference product measure**.

## 2. A uniform causal-boundary sandwich, including short segments

The exact wave relation is imported from the coordinate-oracle derivation
(re-derived in the referee report, Section 3.1). Put
`b = m(u1,u2) = m(u2) - m(u1)`, the signed `Delta m`, so that `b > 0`
exactly when `u1 < u2`. Put `Delta v = v2 - v1`, `Delta y = y2 - y1`, and,
for `u1 < u2`, `E = E_lambda(u1,u2)`, a symmetric trace-free 2x2 matrix
with eigenvalues `+/- e`. Then

- `X <_lam Y` iff `2 Delta v >= Delta y^T (bI+E)^-1 Delta y`, when `u1 < u2`;
- `X <_lam Y` iff `Delta y = 0` and `Delta v >= 0`, when `u1 = u2`;
- `X <_lam Y` never holds when `u1 > u2`, because `u` is non-decreasing
  along causal curves.

**Lemma 2.1 (operator bound).** With
`C_E = 2 sqrt(q) (4 f_max - 2)`, for every `lambda` and every `u1 < u2`
with `l = u2 - u1`,

\[
e\le C_E\,\lambda\min(l,\lambda)\le C_E\,\delta,
\qquad e<b .
\]

*Proof.* `e = t |integral_(u1)^(u2) f e^(iu/lambda) du|` with
`t = sinh(2 beta)/cosh(2 beta) < 1` and `sinh(2 beta) = 2 lambda sqrt(q)`, so
`t <= 2 lambda sqrt(q)`. One integration by parts bounds the integral by
`lambda (f(u1) + f(u2) + TV f) <= lambda (4 f_max - 2)`. Trivially it is at
most `integral f = b <= f_max l`, and `f_max <= 4 f_max - 2`. So
`e <= C_E lambda min(l, lambda)`, and `e <= t b < b`. []

**Lemma 2.2 (sandwich).** Let `c0 = 2 sqrt 2 C_E` and
`g = Delta T - |Delta z|`. For every ordered pair `(X, Y)`,

\[
\{g\ge c_0\delta\}
\subseteq\{X\prec_\lambda Y\}
\subseteq\{g\ge-c_0\delta\}.
\tag{4}
\]

*Proof.* `Delta T = (b + Delta v)/sqrt 2` and
`Delta Z = (b - Delta v)/sqrt 2`, so

\[
g\,(\Delta T+|\Delta z|)=\Delta T^2-|\Delta z|^2=2b\,\Delta v-|\Delta y|^2 .
\]

Put `Q = Delta y^T (bI+E)^-1 Delta y`. By Lemma 2.1, `bI + E` is positive
definite and `|Delta y|^2/(b+e) <= Q <= |Delta y|^2/(b-e)`.

- **Case `u1 < u2`, first inclusion.** Suppose `g >= c0 delta`. Then
  `Delta T > |Delta Z|`, which with `b > 0` gives `Delta v > 0`. So
  `2 b Delta v - |Delta y|^2 >= c0 delta Delta T >= c0 delta Delta v / sqrt 2`.
  Since `c0 delta >= 2 sqrt 2 e` (Lemma 2.1), this is at least
  `2 e Delta v`. Hence `|Delta y|^2 <= 2 Delta v (b - e)`, so
  `Q <= 2 Delta v` and `X <_lam Y`.
- **Case `u1 < u2`, second inclusion.** Suppose `X <_lam Y`. Then
  `2 Delta v >= Q >= 0`, so `Delta v >= 0` and
  `|Delta y|^2 <= 2 Delta v (b + e)`. So
  `g (Delta T + |Delta z|) >= -2 e Delta v`. If `g < 0`, then
  `g >= -2 e Delta v / Delta T = -2 sqrt 2 e Delta v/(b + Delta v)`, which
  is at least `-2 sqrt 2 e >= -c0 delta`.
- **Case `u1 = u2`.** Then `b = 0` and
  `|Delta z| >= |Delta Z| = |Delta T|`, so `g <= 0`, and the first set is
  empty. A wave relation forces `Delta y = 0` and `Delta v >= 0`, so
  `g = 0`, and the second inclusion holds.
- **Case `u1 > u2`.** Then `b < 0`. `g > 0` would need
  `b + Delta v > |b - Delta v| >= Delta v - b`, which is false, so the
  first set is empty. There is no forward wave relation, so the second
  inclusion is vacuous.

The reverse orientation follows by exchanging the endpoints. []

The smallest valid `c0` is the supremum of `|g| / delta` over wave-null
pairs in D and `lambda` in `(0, 1]`.
- The audit referee's sampled lower bounds of it lie below `c0`: 2.96, 7.16
  and 36.6, against 6.39, 18.1 and 132 at q = 0.25, 1 and 4.
- This revision's referee found `|g| <= 0.407 c0 delta` on the wave null
  cone of its grid.

`c0` is a valid constant, not the sharp one.

**Corollary 2.3 (the original's (3), with its hypothesis).** Label the pair
so that `Delta T >= 0`. On the **crossing set**, where the two relations
of the pair differ in either orientation,

\[
\left|\Delta T-|\Delta z|\right|\le c_0\delta .
\tag{3}
\]

*Proof.* By Lemma 2.4 below, the reverse relation fails under both kernels
when `Delta T >= 0`. So the relations can differ only in the forward
orientation. There, by (4) and (2), `g >= c0 delta` forces both relations
and `g < -c0 delta` forbids both. []

The bound does not divide by a short `u` separation: Lemma 2.1 is uniform
in `l`, and Lemma 2.2 uses only `e <= C_E delta` and `e < b`. This uniform
sandwich, not a fixed-gap expansion, is the geometry needed below.

**Lemma 2.4 (a common time function).** Under both kernels,
`X < Y` with `X != Y` implies `T_X < T_Y`.

*Proof.* For the background this is (2). For the wave, a relation gives
`Delta u >= 0` and `Delta v >= 0` (from the case analysis above), so
`Delta m >= 0` and `Delta T = (Delta m + Delta v)/sqrt 2 >= 0`. Equality
forces `Delta u = Delta v = 0`, hence `Delta y = 0` and `X = Y`. []

## 3. Reveal cells, retaining the uncertainty inside them

Partition the `(T,z)` chart into fixed axis-aligned half-open cubes of side
h on a grid. Then the T-intervals of two cubes either coincide or are
disjoint. Fix an edge constant `L > 1 + sqrt 3`, for example `L = 4`. Take

\[
h_{\min}:=\frac{c_0\,\delta}{L-1-\sqrt3}\ \le\ h\ \le\ 1 .
\tag{3a}
\]

This is the precise form of the original's `c0 delta << h`. The audit's
control with `L = 1` violates the non-edge claim (154 violations).

Consider the stronger experiment revealing the cube index of every event
together with the order. Its cube-index distribution is identical in both
reference-sampling models. Conditioned on these indices, the locations
remain independent in their respective cubes.

Discard as a bad mark configuration any cube touching the boundary of D.
The probability is at most `C_bd N h`: such cubes lie within distance `2h`
of the boundary, a layer of volume `O(h)`. Every remaining cube is wholly
inside D. In such a cube, conditioned also on the three spatial
coordinates of an event, its T coordinate has a density on an interval of
length h. That density lies between `rho_min/(rho_max h)` and
`rho_max/(rho_min h)`, and its derivative is at most
`rho'_max/(rho_min h)`.

**Uncertain pairs.** A pair of cubes is an edge if the forward or backward
null margin of their centres, `+/- Delta T_c - |Delta z_c|`, has magnitude
at most `L h`.

- **Every non-edge has one deterministic relation, the same under both
  kernels, at every location in its two cubes.** Moving two events inside
  their cubes changes `Delta T` by at most h and `|Delta z|` by at most
  `sqrt 3 h`. So each true margin differs from its centre value by at most
  `(1 + sqrt 3) h`. At a non-edge both true margins therefore keep the
  sign of their centre values and have magnitude more than
  `(L - 1 - sqrt 3) h >= c0 delta`, by (3a). By (4) and (2), a margin
  above `c0 delta` forces the relation in that orientation under both
  kernels, and one below `-c0 delta` forbids it under both.
- **Orientation.** If an edge joins cubes with disjoint T-intervals, then
  by Lemma 2.4 only the earlier-to-later relation is possible under either
  kernel. The pair is one binary variable with a fixed orientation.

  If an edge joins cubes with the same T-interval, then `Delta T_c = 0` and
  `|Delta z_c| <= L h`. The two events are then within h in T and within
  `(L + sqrt 3) h` in z. Such a pair has probability at most `C h^4`. All
  such pairs are discarded as bad marks, at a cost of at most
  `C_amb N^2 h^4`. Two events in the same cube are included in this case.
- **The uniform bound.** Fix a first event X anywhere in D and any grid.
  If Y's cube forms an edge with X's cube, then
  `| |T_Y - T_cX| - |z_Y - z_cX| | <= (L + (1+sqrt 3)/2) h`. For each
  `z_Y`, this confines `T_Y` to two intervals of total length
  `4 (L + (1+sqrt 3)/2) h`. Hence

  \[
  \Pr(\text{Y's cube forms an edge with X's cube}\mid X)\le C_{\rm unc}\,h,
  \]

  with `C_unc = 4 (L + (1+sqrt 3)/2) vol_3(pi_z D) rho_max`, uniformly in
  X and the grid. Consequently `E U <= C(N,2) C_unc h`, where U is the
  number of uncertain pairs of events.
- **Components.** Discard mark configurations whose uncertain-pair graph
  on the events has a component with more than K vertices. Such a
  component contains a set of K+1 events spanned by a tree of uncertain
  pairs.
  - Fix the set and a labelled tree on it. Order the vertices so that each
    vertex after the first has its parent earlier, and peel the last one,
    a leaf.
  - Its location is independent of the others. The only tree edge it
    enters joins it to its parent, which costs at most `C_unc h`
    conditionally on everything else.
  - By induction the tree costs at most `(C_unc h)^K`. No independence
    between edges is assumed, and cycles are allowed.
  - Cayley's `(K+1)^(K-1)` labelled trees and at most `N^(K+1)/(K+1)!`
    sets then give

\[
\Pr(\text{bad marks})
\le C_{\rm bd}Nh+C_{\rm amb}N^2h^4+\frac{(K+1)^{K-1}}{(K+1)!}N^{K+1}(C_{\rm unc}h)^K .
\tag{5}
\]

## 4. The finite-component smoothing lemma

This section writes the original's sketch as three lemmas. The proof of
Lemma 4.1 is the referee's slicing argument (report, Section 3.5), written
out here. The original gave the structure of this argument without a
constant.

**Lemma 4.1 (Lipschitz gradient).** Let `k >= 2`, and let
`rho(s) = prod_i rho_i(s_i)` on `[0,1]^k`, with `0 <= rho_i <= R_i` and each
`rho_i` Lipschitz with constant `R'_i`. Let rows
`a_e = +/- (e_j - e_i)`, `e = 1..m`, be indexed by distinct unordered pairs
`{i, j}`, and let

\[
p(\tau)=\int_{[0,1]^k}\rho(s)\prod_{e=1}^m\mathbf 1\{a_e\cdot s\ge\tau_e\}\,ds,
\qquad \tau\in\mathbb R^m .
\tag{6}
\]

Then `p` is `C^1` on `R^m`, and `grad p` is Lipschitz in the Euclidean norm
with constant

\[
\Lambda=(m+1)\prod_i R_i+\max_j R'_j\prod_{l\ne j}R_l .
\]

*Proof.* Fix a row `e = (i -> j)`, say `a_e . s = s_j - s_i`. Integrating
`s_j` first,

\[
p(\tau+\eta\,\mathbf e_e)-p(\tau)=-\int_0^\eta G_e(\tau+\sigma\mathbf e_e)\,d\sigma,
\]

where `G_e(tau)` is the integral of `F = rho x (the other indicators) x
(the cube)` over the hyperplane `s_j = s_i + tau_e`, parametrised by the
`k - 1` coordinates other than `s_j`, with unit Jacobian. So
`partial_e p = -G_e`. It remains to bound how `G_e` changes.

- **(i) Another threshold `tau_f`, `f != e`, moved by `eta`.** The
  indicator of row f changes only on a slab of width `eta` in the
  hyperplane. Row f is a different unordered pair. On the hyperplane,
  `a_f . s` therefore depends on at least one free coordinate other than
  `s_i` with coefficient `+/- 1`. This holds whether f shares no vertex
  with e, shares `j` (then `s_j = s_i + tau_e` and the free coordinate is
  f's other vertex, which is not `i`), or shares `i`. So the slab has
  `(k-1)`-measure at most `eta`, and `|Delta G_e| <= eta prod R`.
- **(ii) Its own threshold, moved by `eta`.** For fixed values of the other
  `k - 2` coordinates, the admissible `s_j` form an interval `[lo, hi]`.
  It is cut out by the cube and by the other rows through `j`. None of
  them contains `i`, since the pairs are distinct, so `lo` and `hi` do not
  depend on `s_i`.
  - The indicator of `s_i + tau_e` lying in `[lo, hi]` changes only on two
    `s_i`-intervals of length `eta`. This contributes at most
    `2 eta prod R`.
  - The factor `rho_j(s_i + tau_e)` changes by at most `R'_j eta`. This
    contributes at most `eta R'_j prod_(l != j) R_l`.
- **(iii) Assembly.** So `|G_e(tau) - G_e(tau')| <= sum_f c_ef |tau_f - tau'_f|`,
  with `c_ef = prod R` for `f != e` and
  `c_ee = 2 prod R + R'_j prod_(l != j) R_l`. The matrix `(c_ef)` is
  symmetric and non-negative. Its Euclidean operator norm is at most its
  largest row sum, which is at most `Lambda`. So `G` is continuous, `p` is
  `C^1`, and `grad p` is `Lambda`-Lipschitz.

Linear dependence among three or more rows never enters. []

*Scaling to cubes of side h.* Put `t = h s` and let
`p_h(theta) = p(theta/h)`. Then `grad p_h` is `(Lambda/h^2)`-Lipschitz, with
`R_i = sup (h rho_i^(h))` and `R'_i = sup (h^2 |rho_i^(h)'|)`. Here
`rho_i^(h)` are the conditional T densities of Section 3, so
`R_i <= R := rho_max/rho_min` and, for `h <= 1`,
`R'_i <= R' := rho'_max/rho_min`. Since `R >= 1`, the constant is monotone
in `k <= K` and `m <= M := K(K-1)/2`. One constant serves every component:

\[
\Lambda_K=(M+1)R^K+R'R^{K-1}.
\]

**Lemma 4.2 (square root).** Let `p >= 0` on `R^m` have a
`Lambda`-Lipschitz gradient. Then `|grad p|^2 <= 2 Lambda p`, and
`sqrt p` is `sqrt(Lambda/2)`-Lipschitz on all of `R^m`, including across
the zero set.

*Proof.*
- **The gradient bound.** By the descent inequality at
  `tau - grad p(tau)/Lambda`, `0 <= p(tau) - |grad p(tau)|^2/(2 Lambda)`.
- **Where `p > 0` on a segment,** `|grad sqrt p| = |grad p|/(2 sqrt p)`,
  which is at most `sqrt(Lambda/2)`.
- **Across zeros.** Suppose the segment from `x` to `y` meets the zero set.
  Let `x'` be its first zero and `y'` its last. Then
  `sqrt p(x) = |sqrt p(x) - sqrt p(x')| <= sqrt(Lambda/2) |x - y|`, and
  likewise at `y`. Since both values are non-negative,
  `|sqrt p(x) - sqrt p(y)| <= max(sqrt p(x), sqrt p(y))`. []

No common-support or positive-cell-probability assumption is hidden here.
Patterns with zero reference probability are covered.

**Lemma 4.3 (one component).** On good marks, consider a component with
`k <= K` events and `m <= M` uncertain pairs, each with its fixed
orientation from Section 3. Let `S in {0,1}^m` record which uncertain pairs
are related. Every other relation inside the component is a non-edge, so it
is deterministic and common to both kernels. Then

\[
H^2_{\rm comp}:=\sum_S\Big(\sqrt{p_\lambda(S)}-\sqrt{p_0(S)}\Big)^2
\le 2^{M+1}M\,\Lambda_K\,c_0^2\,\frac{\delta^2}{h^2}=:C_H\,\frac{\delta^2}{h^2}.
\]

*Proof.*
- **Condition on the spatial coordinates.** Condition on the component's
  spatial coordinates z and translate each T-interval to `[0, h]`.
  - The T coordinates are then independent, with the densities of
    Section 3.
  - A background bit 1 for the pair `i -> j` is
    `1{t_j - t_i >= theta_e}`, with `theta_e` fixed by z and the cube
    offsets.
  - A bit 0 is the reversed row with the negated threshold. The
    difference between strict and non-strict inequality has probability
    zero.
  - So `p_0(S | z) = p_h(tau)` for a signed row set of the form (6).
- **Squeeze the wave probability.** By (4), every wave indicator lies
  pointwise between the background indicators with its threshold moved by
  `+c0 delta` and by `-c0 delta`. In signed-row form this gives
  `p_h(tau + c0 delta 1) <= p_lam(S | z) <= p_h(tau - c0 delta 1)`, since
  `p_h` is non-increasing in each threshold.
- **Apply Lemma 4.2.** `sqrt p_lam(S | z)` and `sqrt p_0(S | z)` both lie in
  the interval between `sqrt p_h(tau +/- c0 delta 1)`. By Lemma 4.2 that
  interval has length at most `2 sqrt(Lambda_K/2) c0 delta sqrt(m) / h`.
- **Sum over patterns.** There are at most `2^m` patterns. Joint convexity
  of `H^2` removes the conditioning on z. Forgetting z cannot increase the
  bound. []

The spatial coordinates are an auxiliary device in this lemma, not extra
inputs to the final order-only experiment. The potentially large constants,
`2^(M+1)` and `Lambda_K`, are why K is kept fixed in (1).

## 5. Assemble the bound and choose h

**Convention.** `H^2(P, Q) = sum (sqrt p - sqrt q)^2`, without a factor
1/2. Then `TV = (1/2) sum |p - q| <= H`, by Cauchy-Schwarz with
`sum (sqrt p + sqrt q)^2 <= 4`. For product measures,
`H^2 <= sum_i H_i^2`, because `1 - H^2/2` is multiplicative and
`1 - prod (1 - x_i) <= sum x_i`.

On good marks:
- components use disjoint events, so they are independent under both
  kernels;
- every relation outside the uncertain pairs is deterministic and common to
  both kernels;
- every component with an edge contains at least one uncertain pair, and
  no pair lies in two components. So the number of such components is at
  most U, and components without an edge contribute nothing.

Hence, on good marks,
`TV(order | marks) <= H(order | marks) <= sqrt(U C_H) delta / h`.

The mark law is common to both models. Revealing the marks can only
increase TV. So, averaging over marks and using
`E sqrt U <= sqrt(E U) <= N sqrt(C_unc h / 2)`, for every `h` in (3a),

\[
\operatorname{TV}(P_{N,\lambda},P_{N,0})
\le C_1N\delta+C_{\rm bd}Nh+C_{\rm amb}N^2h^4
+\frac{(K+1)^{K-1}}{(K+1)!}N^{K+1}(C_{\rm unc}h)^K
+\sqrt{C_HC_{\rm unc}/2}\;\frac{N\delta}{\sqrt h}.
\tag{7}
\]

The first term restores the actual wave sampling law. The proof then
forgets all cube indices; the observer for P is still supplied only the
order.

**Choice of h.** Take `h = (delta/N^K)^(2/(2K+1))`. Then exactly

\[
N^{K+1}h^K=\frac{N\delta}{\sqrt h}=(N^{3K+1}\delta^{2K})^{1/(2K+1)},
\qquad Nh=(N\delta^2)^{1/(2K+1)},
\]

and `N^2 h^4 <= (N h)^3` (since `h <= 1 <= N`), which is at most `N h`
when `N h <= 1`. When `N h > 1`, the second term of (8) exceeds
`C_K >= 1`, so (8) is trivial. Substituting `delta = lambda^2` turns (7)
into (8) whenever this h lies in (3a).

**(8) holds for every `N >= 1` and `0 < lambda <= 1`,** after enlarging
`C_K` (finding F5). TV is at most 1, so it is enough that the right side
of (8) is at least 1 whenever h leaves (3a). Since `N >= 1` and
`delta <= 1`, `h <= 1` always holds, so only `h < h_min` needs care.

Suppose `h < h_min = c1 delta`, with `c1 = c0/(L - 1 - sqrt 3)`. Then
`N^(2K) > c1^-(2K+1) delta^(1-2K)`. Hence

\[
N^{3K+1}\delta^{2K}>c_1^{-(2K+1)(3K+1)/(2K)}\,\delta^{-(2K+1)(K-1)/(2K)}
\ge c_1^{-(2K+1)(3K+1)/(2K)},
\]

since `delta <= 1` and `K >= 1`. So the third term of (8) is at least
`C_K c1^-((3K+1)/(2K))`. Taking `C_K >= max(1, c1^((3K+1)/(2K)))`,
together with the constants of (7), covers this case.

**(1) follows from (8).** Let `N >= 1` and `0 < lambda <= 1` vary so that
`N^(3K+1) lambda^(4K) -> 0`.
- The third term of (8) tends to zero.
- Since `(N lambda^2)^(2K) = N^(2K) lambda^(4K) <= N^(3K+1) lambda^(4K)`,
  also `N lambda^2 -> 0`.
- `N lambda^4 <= N lambda^2`, so the second term tends to zero as well.

So the right side of (8) tends to zero, and with it TV. This proves (1).

**(1) is exactly what (7) delivers.** For every h,
`(N^(K+1) h^K) (N delta / sqrt h)^(2K) = N^(3K+1) delta^(2K)`. So one of
these two terms of (7) is at least a fixed positive constant times
`(N^(3K+1) delta^(2K))^(1/(2K+1))`, and
the right side of (7) can tend to zero, for any choice of h, only under
(1).

**This h is not the minimiser of (7).** It balances the tree term against
the Hellinger term. When `N h` dominates, a smaller h gives a better rate:
up to a factor 1840 on the audit referee's grid. This costs rate only, not
the condition (1).

**K = 1.** Then `M = 0`, and components on good marks have one vertex, so
`U = 0` and the last term of (7) is absent. (8) with `K = 1` reads
`TV <= C_K [N lambda^2 + (N lambda^4)^(1/3) + (N lambda)^(4/3)]`, which
returns the earlier condition `N lambda -> 0`.

**Exponents.** At `N = lambda^-alpha`, the third term of (8) decays with
exponent `(4K - (3K+1) alpha)/(2K+1)`, which is positive iff
`alpha < 4K/(3K+1)`. The other two terms vanish for every `alpha < 2`. So
every fixed `alpha < 4/3` is covered by `K = floor(alpha/(4 - 3 alpha)) + 1`.
The coverage is not uniform as `alpha -> 4/3`. Examples of sufficient
vanishing conditions:

| K | condition |
|---:|---|
| 2 | `N lambda^(8/7) -> 0` |
| 3 | `N lambda^(6/5) -> 0` |
| 4 | `N lambda^(16/13) -> 0` |

These are necessary-resource bounds, not detection tests.

## 6. Core lemma check and what remains

*Sections 6 and 7 are the original's text, unchanged, and partly stale.*
- *"This proof" in them means the frozen original at `5dda03d`. The check
  they describe exercised the original's lemma.*
- *Section 6's bracket sentence is the state on 2026-09-24. The current
  statement is AUDIT.md Section 2's
  `4/3 <= alpha_low <= alpha_high <= 2`.*
- *"Sections 1-5 here" in Section 7 means the original's Sections 1 to 5.*
- *How Lemma 4.1 here was tested:*
  - *by the audit referee, with exact rationals for uniform densities up to
    k = 5 (report, Section 3.5);*
  - *by this revision's referee, with exact piecewise polynomials for
    non-uniform densities at k = 3
    ([ref_checks.py](independent_audit/argument1_revision/ref_checks.py),
    R5);*
  - *by the author's
    [checks.py](independent_audit/argument1_revision/checks.py), only for
    step (ii) at k = 2 and for Lemma 4.2's ratio on the uniform triangle.*

Before running its small deterministic check, freeze this proof and
`experiments/positive_control/order_wave_cluster_check.py` in Git.
The check enumerates vertices of three-dimensional polytopes and evaluates
their volumes with ConvexHull, without using a smoothing or square-root bound.
For independent uniform times in [0,1], the cyclic relaxed inequalities
`T2-T1>=-eps`, `T3-T2>=-eps`, `T1-T3>=-eps` have exact probability
`(9/2)eps^2-4eps^3` for `0<=eps<=1/2`. This deliberately exercises a
reference-impossible pattern becoming possible, where assuming a positive
reference probability would invalidate a Hellinger argument.

Freeze eps in `{0,0.001,0.01,0.1,0.4}` and numerical absolute tolerance
`2e-12` for that identity and for the sum of all eight pattern probabilities.
Also record the full-pattern Hellinger distance from eps=0. This is a
deterministic algebra/geometry control, not a finite-event power campaign,
an independent researcher audit, or a proof for arbitrarily large K.

The full-order bracket after the later [shared-reference improvement](local_null_challenge.md)
is: indistinguishable at every fixed power below 4/3, and distinguishable when
`N lambda^3/log N -> infinity` for known q<=1. A matching law,
an observer attaining it, and novelty verification remain open. The
coordinate-oracle formula is not substituted for that missing result.

## 7. Frozen checks completed; provenance and limits

The coordinate-oracle design and implementation were committed at
`b1bf8801a8416baf93ec01bcf85d42db13c5acb2` before the 18 primary cases and
one refinement. The CPU run took 47.65 seconds. Its
[manifest](oracle_barrier_evidence/manifest.json),
[raw rows](oracle_barrier_evidence/rows.jsonl), and
[summary](oracle_barrier_evidence/summary.json) are preserved. Raw SHA256:
`faa777d3b665bb7172c591e3ab7f2dd6e6eaf18a37235a871564a66dd9816fe0`.

| q | Predicted kappa(q) | At r=1024: twice exact reference loss / (lambda^2 kappa) |
|---:|---:|---:|
| 0.25 | 0.3221254984 | 1.0000181590 |
| 1 | 0.6316283583 | 1.0000475553 |
| 4 | 1.1375203435 | 1.0000828696 |

Doubling all quadrature orders and tightening the ODE at q=1,r=1024 changed
the scaled reference-loss integral by `2.0638319e-7` relatively, within the
pre-execution numerical tolerance `1e-5`. Refining the limiting coefficient
changed it by `1.1623540e-9` relatively. The lower r=32 case at q=0.25 lies
below its limiting coefficient; all rows are retained, and no monotonicity
or exact finite-frequency scaling was assumed. These are deterministic
consistency checks, not confidence intervals or certified integration bounds.

The quantized-order proof and polytope control were then frozen at
`5dda03d9786a96d4e8b6891cff8186e0e08d704c`. The
[separate cluster check](oracle_barrier_evidence/cluster_check.json) passed
all five frozen polynomial and probability-sum checks. At eps=0.001 the
new cyclic-pattern mass was `4.496e-6`, agreeing with the exact polynomial;
the full-pattern squared Hellinger distance divided by eps^2 was 13.44677.
This pattern is used to challenge a mathematical bound; no cyclic order
is asserted to arise from a wave spacetime. The calculation does not prove
the arbitrary-K lemma, whose proof is Section 4.

All four portable evidence files were copied byte-for-byte from the original
output. Both stamps precede their respective executions and remain reachable
in the actual branch history. The evidence does not replace a separate
researcher's audit of either mathematical derivation.

The old coordinate-oracle design records the pre-existing necessity exponent
1 in its final section. Sections 1-5 here are the subsequent improvement,
not a retroactive change to that design or its frozen numerical predictions.

## 8. Status, scope and what remains

- **What this revision completes.** The audit's items 1 to 6 for argument
  1 are now in the text:
  - the constant link (3a);
  - the hypothesis of (3), as Corollary 2.3, and the common orientation, as
    Lemma 2.4;
  - the smoothing lemma with its constant and proof (Lemmas 4.1 to 4.3);
  - the Hellinger convention and the component count;
  - `K = 1`, and the remark that the chosen h is not the minimiser;
  - the strict separation's dependency.

  Also written:
  - the sandwich (4) for every ordered pair, with an explicit `c0`;
  - the tree count with its constants;
  - the validity of (8) for every `N >= 1` and `0 < lambda <= 1`;
  - the step from (8) to (1);
  - why (1) is exactly what (7) delivers (Section 5).
- **What is imported or cited, not derived here:**
  - the model itself: `P_{N,lambda}`, `P_{N,0}`, the effective background,
    the diamond and its tips, and the unlabelled order
    ([full_order_resolution_probe.md](full_order_resolution_probe.md));
  - the background relation (2) and the `a^4` density in the `(T, z)`
    chart, asserted in Section 1 and derived in the audit referee's report,
    Section 3.1;
  - the exact causal relation of the wave in the Rosen chart, and the form
    of `E_lambda` (from `coordinate_oracle_barrier.md`, re-derived in the
    same report, Section 3.1);
  - the one-point bound `TV(mu_lam, mu_0) <= C_1 lambda^2`, uniform in
    phase and in `lambda <= 1`. This is the domain-strip bound, re-derived
    in that report, Section 3.3. Every `alpha > 1` depends on it. This
    revision's referee names it the step it trusts least, since it is
    neither written here nor given a constant.
- **The separation of protocols** in the preamble rests on `kappa(q) > 0`
  and the oracle-side statements of `coordinate_oracle_barrier.md`, which
  the audit did not audit. The order-only bound (1) does not depend on
  them.
- **Constants (finding F7).** (8) is stated for every N, but `C_K` is not
  evaluated, and the bound has no content at reachable sizes.
  - `C_K` contains the tree coefficient `(K+1)^(K-1)/(K+1)!`. Through
    `sqrt(C_H)` it also contains `2^((M+1)/2)`, with `M = K(K-1)/2`: about
    7.0e13 at K = 14.
  - `c0 = 2 sqrt 2 C_E` is 6.39, 18.1, 132 and 13532 at q = 0.25, 1, 4
    and 9.
  - `alpha = 1.3` needs `K >= 14`. At K = 14 the main-term exponent is
    0.00345. It grows with K toward `(4 - 3 alpha)/2 = 0.05`: 0.0171 at
    K = 20 and 0.0366 at K = 50. `C_K` grows with K too.
- **What it does not change.** The bracket `4/3 <= alpha_low <= alpha_high <= 2`
  is the audit's statement (AUDIT.md Section 2), not something new here.
  Nothing is proved at `alpha = 4/3`. Whether 4/3 is the true boundary is
  open.

## 9. Review and the corrections made after it

[independent_audit/argument1_revision/review.md](independent_audit/argument1_revision/review.md).
An AI referee in a fresh context was given this revision as committed at
`a032f0e`, together with:
- the frozen original;
- AUDIT.md;
- the audit referee's report;
- `coordinate_oracle_barrier.md`;
- the author's checks.

It ran only short deterministic checks.

**Verdict: ERROR (E1, E2)**, both in the revision's account of itself.
Neither touches an inequality: the referee found no false step in
Sections 2 to 5, and the proof of (8) does not change. With E1 and E2
fixed, its grade is CORRECT WITH MINOR FIXES. In more detail:
- AUDIT.md's items 1 to 6 for argument 1 are all addressed.
- The SHA-256 recorded for the original is right.
- The original's Sections 1 to 5 are byte-identical to the freeze.

Corrected after the review, in the text above:

- **E1.** Section 8 said that "(1) implies (8) tends to zero" was written.
  It was not: the revision had dropped the original's "Under (1) ... This
  proves (1)" without replacing it. Section 5 now has "(1) follows from
  (8)". It also has "(1) is exactly what (7) delivers", which replaces a
  citation.
- **E2.** Section 0 did not list every change. The missing rows are added:
  - the density constants and `C_1`;
  - the half-open aligned grid, which the orientation argument needs;
  - the changed discard rule;
  - the boundary-layer reason;
  - the removed "uniformly Lipschitz" sentence;
  - the replaced closing sentences of Section 5;
  - the Exponents paragraph.
- **M1, M2, M10.** Section 8's constants are restated:
  - the Hellinger factor is `2^((M+1)/2)`, not `2^(M+1)`;
  - the tree coefficient is `(K+1)^(K-1)/(K+1)!`;
  - the `alpha = 1.3` margin is 0.00345 at K = 14, growing with K;
  - "no finite-N statement" becomes "`C_K` is not evaluated; no content at
    reachable sizes".
- **M3.** The `K = 1` display uses `C_K`, not `C_1`.
- **M4.** The preamble says what the later commits to the original changed.
- **M5.** The preamble points to the review record.
- **M6.** Section 8 lists two further imports: the model definitions, and
  the background relation (2) with its `a^4` density.
- **M7.** `b` is the signed `Delta m` throughout Lemma 2.2.
- **M8.** The sharp-constant remark says over which set the supremum runs,
  and that the quoted values are sampled lower bounds.
- **M9.** The note at the head of Section 6 says what each check of
  Lemma 4.1 covers. "The same constant structure" now reads "without a
  constant". The author's Lemma 2.4 check, which could not fail on its
  grid, now ranges over negative `Delta v`.
- **M11.** The same note says which sentences of Sections 6 and 7 are
  stale.
- **M12.** The vacuous case `h > 1` is removed: `h <= 1` always holds.

**The step the referee trusts least** is the imported one-point bound
`TV(mu_lam, mu_0) <= C_1 lambda^2`, uniform in phase and in `lambda <= 1`.
Every `alpha > 1` depends on it.
- This revision neither writes it nor gives its constant.
- The referee confirmed it numerically only at three values of q and one
  phase, with TV/lambda^2 from 1.6 to 8.1.

Among the written steps, it trusts least Lemma 4.1 for general k. It
re-derived the lemma in full but tested it numerically only at k = 3.
