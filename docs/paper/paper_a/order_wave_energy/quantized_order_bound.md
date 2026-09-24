# A stronger impossibility region after hiding exact coordinates

2026-09-24. Analytic result after the coordinate-oracle calculation. The
sharp order-only boundary remains unresolved. This proof improves the
necessary side; it does not supply an attaining test. It has not received
independent researcher review.

For every fixed `0<q<pi^2` in the exact circular family, and every fixed
integer K>=2, the result is

\[
\boxed{N^{3K+1}\lambda^{4K}\longrightarrow0
\quad\Longrightarrow\quad
\operatorname{TV}(P_{N,\lambda},P_{N,0})\longrightarrow0.}
\tag{1}
\]

In particular K=2 gives `N lambda^(8/7) -> 0`. Letting K be a sufficiently
large **fixed** integer proves indistinguishability at every fixed power
`N=O(lambda^(-alpha))` with `alpha<4/3`. This is **not** the statement
`N lambda^(4/3) -> 0`: constants below depend on K, and no growing-K estimate
or logarithmic critical boundary is claimed.

The [coordinate-oracle result](coordinate_oracle_barrier.md) has a sharp
boundary at `N lambda~1`. Thus, for any fixed `1<alpha<4/3`, at
`N~lambda^(-alpha)` exact coordinates plus order distinguish consistently,
while the full finite order alone is asymptotically indistinguishable.
This is a strict separation of observation protocols, conditional on the
proofs. It is not evidence that the order-only sufficient exponent 4 is sharp.

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
with bounded derivative. Constants can depend on the fixed q and positive
patch margin. They are not uniform as q approaches pi^2.

Write `delta=lambda^2`. The normalized wave and reference one-point laws
differ by at most `C delta` in TV, by the existing domain-strip bound.
Replacing the wave sampling locations by N iid reference locations therefore
costs at most `CN delta` after applying any order map. For the remaining
argument **both models use exactly the same reference product measure**.

## 2. A uniform causal-boundary displacement, including short segments

For ordered u endpoints let `b=m(u1,u2)` and `E=E_lambda(u1,u2)`.
The exact family satisfies

\[
\|E\|\le C\lambda\min(u_2-u_1,\lambda)\le C\delta,
\quad \|E\|/b\le C\lambda.
\]

The wave action is `Delta y^T(bI+E)^-1 Delta y/2`, compared with
`|Delta y|^2/(2b)`. If their causal relations differ, the longitudinal
difference is nonnegative and comparable to `|Delta y|^2/b`. Consequently

\[
\left|\Delta T-|\Delta z|\right|
=\frac{|2b\Delta v-|\Delta y|^2|}{\Delta T+|\Delta z|}
\le C\|E\|\le C\delta.
\tag{3}
\]

The estimate uses `|Delta y|^2/(b+Delta v)<=Cb` on the crossing set;
it does not divide by an uncontrolled short u-separation. Equivalently,
with a fixed constant c0,

\[
\{\Delta T-|\Delta z|\ge c_0\delta\}
\subseteq\{X\prec_\lambda Y\}
\subseteq\{\Delta T-|\Delta z|\ge-c_0\delta\}.
\tag{4}
\]

For the first inclusion, start in the background cone, apply the same action
bound, and note `b>0, Delta v>0`. For the second, start in the wave cone;
ellipticity supplies the same estimate. Reverse orientations are treated
by interchanging endpoints. This uniform sandwich, rather than a fixed-gap
expansion, is the geometry needed below.

## 3. Reveal cells, retaining the uncertainty inside them

Partition the `(T,z)` chart into fixed axis-aligned cubes of side h, with
`c0 delta << h << 1`. Consider the stronger experiment revealing the cube
index of every event together with the order. Its cube-index distribution
is identical in both reference-sampling models. Conditioned on these
indices, the locations remain independent in their respective cubes.

Discard as a bad mark configuration any cube touching the boundary of D.
The probability is at most `CNh`: D has a boundary layer of volume O(h).
Every remaining cube is wholly inside D. In such a cube, conditioned also
on the three spatial coordinates of an event, its T coordinate has a density
on an interval of length h between `c/h` and `C/h`, with derivative at most
`C/h`. These bounds follow directly from the fixed positive smooth density.

Make a graph of **uncertain pairs of cubes**. A pair is an edge if a forward
or backward null margin of its centres has magnitude at most Lh, for a fixed
large L. The null-margin functions are uniformly Lipschitz. Thus every
nonedge has a single deterministic relation, shared by both kernels on
every location in those cubes, by (4).

Uniformly in the first location and in the grid, another event lies in an
uncertain cube with probability at most Ch. This follows by integrating T
over a strip of width O(h), then spatial coordinates over a bounded region.
Hence the expected number of uncertain pairs is at most `CN^2h`.
If an uncertain pair's time orientation is not fixed throughout its two
cubes, the two centres are within O(h) in all four coordinates. Discarding
all such pairs costs at most `CN^2h^4`.

Discard mark configurations whose uncertain-pair graph has a component
with more than K vertices. Such a component contains a tree on K+1 vertices.
For each fixed labelled tree, integrate leaf vertices in turn; every required
edge costs a factor Ch by the preceding uniform bound. Counting finitely
many trees for fixed K gives

\[
\Pr(\text{bad marks})
\le C_K[N h+N^2h^4+N^{K+1}h^K].
\tag{5}
\]

No edge-independence assumption is made for a component. Cycles and all
other graphs on at most K vertices are allowed.

## 4. The finite-component smoothing lemma

For k<=K independent T coordinates in intervals of length h as above,
consider m distinct unordered pairs, with one fixed orientation per pair.
Condition on the spatial coordinates. Any background binary relation pattern
has probability of the form

\[
p(\theta)=\int_{[0,h]^k}\rho(t)
\prod_{e=1}^m\mathbf1\{a_e\cdot t\ge\theta_e\}\,dt,
\tag{6}
\]

after translating the intervals. Each row a_e is `+/- (e_j-e_i)`; false
bits reverse the row and threshold. All rows concern distinct unordered
pairs, so any two are linearly independent. The probability density is the
product of the k bounded, positive, continuously differentiable densities.

**Claim:** globally in theta, `p` is nonnegative and has Lipschitz gradient
with constant `C_K/h^2`. To see the bound, first smooth the indicators.
Mixed second derivatives for two different rows insert two delta kernels.
There are two integration variables with Jacobian determinant +/-1, even
when the edges share a vertex. The remaining integration volume, multiplied
by the bounded density, gives `C_K/h^2`.

For a repeated derivative, integrate the derivative of its delta kernel
by parts in one of that row's endpoint variables. The two cube faces give
a face plus one delta constraint; these are independent because the row
contains two endpoints. They also contribute at most `C_K/h^2`. Derivatives
on the other indicators give the already bounded mixed terms; the derivative
on the density gives at most `C_K/h`, hence `C_K/h^2` for h<=1. All bounds
are uniform in the smoothing width. Passing to the limit establishes the
claim. Equivalently, the first derivatives are continuous boundary integrals:
every other row and every single cube face intersects that boundary in
lower dimension. Linear dependencies among three or more rows do not spoil
the two-row estimates. Only boundedness, not a formula for each polytope's
volume, is needed.

A nonnegative function with L-Lipschitz gradient satisfies
`|gradient p|^2<=2Lp`: apply the descent inequality at
`theta-gradient(p)/L` and use nonnegativity there. Therefore `sqrt(p)` is
globally Lipschitz with constant `sqrt(L/2)`, including at p=0 by continuity.
This deals with relation patterns that have zero reference probability;
no common-support or positive-cell-probability assumption is hidden here.

For any requested pattern, (4) sandwiches its wave probability between
the background probabilities with all its signed thresholds moved by
`+/- c0 delta`. The square-root bound consequently gives

\[
|\sqrt{p_\lambda(S)}-\sqrt{p_0(S)}|\le C_K\delta/h.
\]

Summing over the at most `2^(K(K-1)/2)` patterns shows that the squared
Hellinger distance of the entire component is at most `C_K delta^2/h^2`.
Average over the common conditional spatial distribution, then forget
those spatial coordinates; the bound cannot increase. They are an auxiliary
device in this lemma, not extra inputs to the final order-only experiment.
The potentially large constants are why K is kept fixed in (1).

## 5. Assemble the bound and choose h

On good marks, components use disjoint event locations and hence are
independent under both kernels. All remaining relations are deterministic.
Squared Hellinger distances of product distributions are at most their sum,
and TV is bounded by Hellinger distance with this convention. Averaging over
the common mark law and using Cauchy-Schwarz and the uncertain-pair count,

\[
\operatorname{TV}(P_{N,\lambda},P_{N,0})
\le C_K\left[N\delta+Nh+N^2h^4+N^{K+1}h^K
                    +\frac{N\delta}{\sqrt h}\right].
\tag{7}
\]

The first term restores the actual wave sampling law. The proof then forgets
all cube indices; the observer for P is still supplied only the order.

Choose `h=(delta/N^K)^(2/(2K+1))`. Under (1), h tends to zero, h/delta tends
to infinity, and Nh tends to zero. The `N^2h^4` term is then dominated by Nh.
Substituting `delta=lambda^2` gives

\[
\boxed{\operatorname{TV}(P_{N,\lambda},P_{N,0})
\le C_K\left[N\lambda^2+(N\lambda^4)^{1/(2K+1)}
 +(N^{3K+1}\lambda^{4K})^{1/(2K+1)}\right].}
\tag{8}
\]

This proves (1). Examples of sufficient vanishing conditions are
`N lambda^(8/7)->0` for K=2, `N lambda^(6/5)->0` for K=3, and
`N lambda^(16/13)->0` for K=4. Increasing K approaches the exponent 4/3
from below. These are necessary-resource bounds, not detection tests.

## 6. Core lemma check and what remains

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

The present full-order bracket is now: indistinguishable at every fixed
power below 4/3, and distinguishable when `N lambda^4/log N -> infinity`
under the earlier known-q<=1 null-coordinate construction. A matching law,
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
