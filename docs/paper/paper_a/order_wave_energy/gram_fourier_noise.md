# Shared-reference Fourier noise: proposed sufficient theorem

2026-09-24. This extends the [frozen population challenge](gram_fourier_challenge.md).
The proof below is a same-author derivation, not an independent audit. It gives
an explicit mathematical test family. The corresponding complete, precision-
controlled finite-order implementation and a finite-data power assessment are
not yet delivered. Accordingly the established programme record remains
**insufficient below 4/3, sufficient at 3**. Do not silently promote its upper
endpoint using the population calculation or this proof draft.

The proposed result is: for each fixed integer K>=1 there is a randomized,
order-only, wavelength-free test, with known fixed 0<q<=1, separating the exact
circular vacuum family from its effective background whenever

\[
N\asymp\lambda^{-\alpha},\qquad \alpha>2+2/K.                 \tag{1}
\]

Consequently the family would cover each fixed power above 2. It does not
establish a critical-window law at 2, show that 2 is necessary, or give a single
test uniform in a growing K. The impossibility exponent is still below 4/3.

## 1. Raw moments and the pointwise zero

Condition on a good, compact anchor triple, selected using disjoint pilot
events. References R_j and probes X_i are independent of that selection.
Use the bounded order kernels

\[
p_X=E1_{R<X},\quad f_X=E1_{X<R},\quad
a_X=E[1_{R<X}1_{R'<X}1_{R\sim R'}],\quad
b_{AX}=E1_{A<R<X}.
\]

Here ~ means comparable, not a link or a null relation. Estimate a by the
unbiased order-two U-statistic over **all reference pairs**, before forming
`c_hat=a_hat/p_hat^2`. In particular, a conditional ancestor-pair fraction
has a different finite-n denominator and must not be silently substituted
in this expansion. Anchor p,f,a and anchor-to-anchor b use the same references.

Collect these raw probabilities into a finite vector mu(X), and let epsilon(X)
be its centered empirical error. All coordinates of epsilon are averages or
order-two U-statistics. Define A(mu)=w(u,v,R^2)H, with H and the state calibration
in the challenge. Under the background

\[
A(\mu_0(X))=0\quad\hbox{for every }X.                       \tag{2}
\]

Work with a smooth compact extension of A and u in raw-moment space. This is
possible because the active probes are strictly inside the diamond, all
anchor-to-probe gaps and slacks are positive, the inverse chain calibration
has nonzero derivative for fixed q>0, and the anchor Gram block is positive
definite. Explicitly: choose inner and outer compact neighbourhoods of these
regular population vectors, use a smooth cutoff equal to one on the inner
set and zero outside the outer set, and extend the calibration smoothly
outside the outer set. A zero extension then preserves (2) for every
background probe; only the regular part can carry a nonzero weight.

The construction must use computable, fixed cutoffs in a complete instrument;
"discard errors after looking at the data" is not an acceptable replacement.
Every derivative up to a fixed order is bounded, with constants depending on
q, K and these compact sets. There is no assertion of uniformity as q->0.

The score is the probe average of

\[
F_\omega(\mu_0(X)+\epsilon(X))
=A(\mu_0(X)+\epsilon(X))e^{-i\omega u(\mu_0(X)+\epsilon(X))}.
\]

## 2. A bounded-variation lemma specific to this background

Fix (v,y) and vary u. In the background this curve is a future null generator.
Each past indicator increases at most once; each future indicator decreases
at most once. An anchor-interval indicator is a past indicator multiplied by
the constant 1_{A<R}. The raw ancestor-pair kernel increases at most once.
All therefore have uniformly bounded variation in u.

For the order-two kernel h, write its Hoeffding decomposition

\[
U_nh-Eh=\frac2n\sum_jh_1(R_j)+
\binom n2^{-1}\sum_{j<k}h_2(R_j,R_k).
\]

Taking an expectation of a bounded-variation function does not increase its
variation bound. Thus h_1 and the canonical h_2 retain uniform variation,
and h_2 is centered in each reference argument. Anchor kernels are constant
in the probe u. Products of any fixed number of these kernels also have
bounded variation. Multiplying by a smooth coefficient and the background
density preserves the bound. Integration by parts as a Stieltjes integral
gives, for every fixed reference-index configuration,

\[
\left|\int e^{-i\omega u}B(u)\,du\right|
\le\frac{|B(\hbox{left})|+|B(\hbox{right})|+\operatorname{Var}B}{\omega}.
                                                               \tag{3}
\]

All integrations in v,y are over a fixed compact set. The probe weight is
strictly inside the diamond, so moving domain boundaries introduce no
uncontrolled surface term. This is why independent noise attached separately
to every probe would be the wrong error model: it loses (3).

## 3. Taylor terms, repeated indices, and the remainder

In the order-k Taylor term of F around mu_0, the coefficient has degree at
most k-1 in omega. The degree-k term would multiply A(mu_0), which is zero.
For any monomial, expand each empirical error into its order-one projection
and its canonical order-two part. Suppose this produces L reference slots
in the product, k<=L<=2k, with normalization comparable to n^(-L).

To bound its squared expectation, use two copies, with L+L' slots. If a
reference index appears only once across those slots, conditioning on all
other indices kills the term, even before the u integrals. Consequently a
surviving equality pattern has at most (L+L')/2 distinct indices. There are
only finitely many patterns at fixed k. Their normalization and counts give
at most C_k n^(-(L+L')/2), hence at most C_k n^(-k).

For each surviving configuration, (3) applies separately in the two probe
integrals. Thus the L2 norm of the integrated degree-k term is bounded by

\[
C_k\omega^{k-2}n^{-k/2},\quad\omega\ge1.                    \tag{4}
\]

This argument includes repeated indices within a product, bias terms and
cycles involving canonical order-two kernels. It does not assume independent
edges or discard the repeated-index terms. For k=1 it is 1/(omega sqrt(n));
for k=2 it is 1/n.

The remainder requires care: directly bounding the (K+1)-st derivative of
F would appear to cost omega^(K+1). Instead, on the segment mu_0+t epsilon,

\[
A(\mu_0+t\epsilon)
=t\int_0^1 DA(\mu_0+st\epsilon)[\epsilon]ds=tJ(t).
\]

Expand J(t) exp(-i omega u(mu_0+t epsilon)) through degree K-1, then multiply
by t. Its order-K remainder is bounded by
`C_K omega^K |epsilon|^(K+1)`. Bounded-kernel U-statistic moment bounds give
`E|epsilon|^(2K+2)<=C_K n^(-(K+1))`; this can also be obtained by the same
index-centering count without invoking an asymptotic normal approximation.

Finally the conditional variance from n_P~n independent probes is at most
C/n^2: by (2), |A(mu_0+epsilon)|<=C|epsilon|. Combining terms yields

\[
\|S_n(\omega)\|_2\le C_K\left[
\sum_{k=1}^K\omega^{k-2}n^{-k/2}
+\frac1n+\frac{\omega^K}{n^{(K+1)/2}}\right].                \tag{5}
\]

For 1<=omega<=n^(K/(2K+2)), the sum is bounded by a constant times its first
term, and the remainder is no larger than that term. The three-term candidate
bound in the challenge is therefore valid on the frequency range the test
uses. No claim of that simplified bound for arbitrary omega is needed.

As an elementary check of the mechanism, let X,R be uniform on [0,1],
H_hat(X)=F_n(X)-X, and omega=2 pi j, j a nonzero integer. With n references
and m independent probes the exact variance is
`(1-1/m)/(n omega^2)+1/(6nm)`. This follows by integrating the empirical-CDF
error once by parts. It explains both leading terms of (5), but is not a
validation of the nonlinear causal-set instrument.

## 4. Canonical density and kernel coupling, including the tips

Let s be distance in u to the nearer tip. For the positive trace-free
primitive B, write M_left=mI-B_left and similarly on the right. The previous
oscillatory integral bound implies a relative matrix error at most

\[
e(s)=C\lambda\min(1,\lambda/s).
\]

The canonical slice map y=T_lambda(u)z has
`T=(W_0(M_left^-1+M_right^-1))^(-1/2)`. Its first-order perturbation is
trace-free. Matrix Taylor expansion, uniformly where e<=C lambda, gives

\[
|\det T-1|\le C e(s)^2.                                    \tag{6}
\]

The background slice volume is proportional to W_0~s. For s<lambda the
integrated error is O(lambda^2 int_0^lambda s ds)=O(lambda^4); for s>lambda
it is O(lambda^4 int_lambda^1 ds/s). Normalization changes by the same order.
Thus the canonical one-point laws satisfy

\[
\|\nu_\lambda-\mu_0\|_{TV}\le C\lambda^4\log(1/\lambda).
                                                               \tag{7}
\]

This is stronger than the common-physical-chart O(lambda^2) estimate because
the canonical map preserves both tip slacks. It does not follow from an
interior-only expansion.

For later order-kernel estimates, the transverse shift obeys
`|delta y|<=C lambda min(sqrt(s),lambda/sqrt(s))`. Its L1 average is
O(lambda^2). The longitudinal shift is uniformly O(lambda^2). For the latter,
use either
`delta v=(z^T T M_left^-1 T z-|z|^2/m)/2`
or, from T(M_left^-1+M_right^-1)T=W_0^-1 I, the corresponding right expression.
Choose the expression with the long denominator; |z|^2<=Cs removes the tip
singularity. These are uniform bounds, not merely formal first variations.

Two independently sampled reference points have background-cone disagreements
under this coordinate change on O(lambda^2) mass. Indeed the Minkowski null
margin is Lipschitz in both points, and its conditional strip probability is
uniformly at most C times the strip width. Split the crossing event into the
strip attributable to either endpoint, condition on that endpoint, and
integrate its displacement. The wave-versus-background cone difference in
the common physical chart adds O(lambda^2), by the existing uniform sandwich.

The same estimate holds with one endpoint fixed in the compact anchor/probe
region. Split the other endpoint into an interior region, where all shifts
are O(lambda^2), and the tips. At a tip its u-separation from the fixed point
is bounded below; the action is then uniformly Lipschitz in transverse
coordinates and longitudinal slack. Its strip width is bounded by
`C lambda^2/sqrt(s)`, which is integrable. Hence, uniformly for active X,

\[
\Pr\{K_\lambda(X,R,R')\ne K_0(X,R,R')\}\le C\lambda^2       \tag{8}
\]

for each raw indicator kernel. For ancestor pairs, apply a union bound to
the two ancestor indicators and the reference-reference comparison. Couple
the reference locations using (7); its additional mismatch probability is
smaller than (8).

## 5. Alternative noise without a whole-sample density coupling

Let epsilon_lambda and epsilon_0 be centered errors from these coupled
reference lists. For every fixed p>=2, (8), boundedness and the independent-
pair representation of an order-two U-statistic imply, when n lambda^2>=1,

\[
\|\epsilon_\lambda(X)-\epsilon_0(X)\|_p\le C_p\lambda/\sqrt n.
                                                               \tag{9}
\]

For completeness, a bounded centered difference kernel with variance
O(lambda^2), averaged over independent pairs, has p-th norm bounded by
`C_p[lambda/sqrt(n)+lambda^(2/p)/n^(1-1/p)]`. The second term is no larger
when n lambda^2>=1. Averaging over pairings and Jensen transfers this bound
to the U-statistic. The order-one case is the same calculation.

On the active compact set, mu_lambda-mu_0=O(lambda^2), A(mu_lambda)=O(lambda^2),
and the one-point state errors are O(lambda^4 log(1/lambda)). At
omega=1/lambda+O(1), derivatives satisfy

\[
|DF|\le C(1+\omega|A|),\qquad
|D^2F|\le C(1+\omega+\omega^2|A|).
\]

Compare the two centered plug-in errors by first replacing epsilon_lambda
by epsilon_0, then moving the population vector from mu_lambda to mu_0.
The first change is bounded in L2 by (9); the second by
`C omega |mu_lambda-mu_0|/sqrt(n)`. Higher-moment versions of (9) control the
terms containing omega epsilon. Since omega/sqrt(n)->0, both changes are
O(lambda/sqrt(n)). This avoids treating the empirical coordinate error as
independent noise at each probe.

Changing the probe integration law from mu_0 to nu_lambda costs at most
`C lambda^4 log(1/lambda)/sqrt(n)` for the background empirical error, using
(7) and its uniform pointwise L2 bound. Conditional probe averaging contributes
O(1/n+lambda^2/sqrt(n)). There is **no** assumption that all N coupled sample
coordinates agree, and hence no unnecessary restriction N lambda^4 log ->0.
Together with (5), at the target frequency the total score error is

\[
O_{L^2}\left(\lambda/\sqrt n+1/n+
\omega^K n^{-(K+1)/2}\right).                              \tag{10}
\]

## 6. The signal survives intrinsic anchor selection

For general anchor vectors a_i and ell_i with sum ell_i a_i=0, define

\[
D_A(u)=-\frac12\sum_{i=1}^3
\frac{\ell_i(a_{i1}-ia_{i2})^2}{m(u)-m_i}.
\]

This is only a proof coordinate system. Rotate it away using the anchor Gram
matrix; the test never receives it. For the frozen triple, it reduces to D
in the challenge, whose imaginary part is strictly positive. Its complex
integral is therefore nonzero. The same holds uniformly in a sufficiently
small open neighbourhood, including its rotations and reflections. Reflection
reverses the complex helicity representation; it changes which conjugate
coefficient a positive frequency selects, while preserving its magnitude.

Such a neighbourhood is specified intrinsically by u_i,v_i,R_i^2 and g_ij.
For a concrete compact neighbourhood, use outer tolerances 0.002 for the
three u_i and v_i, and 0.00002 for all Gram entries, around the frozen triple.
In the orientation a_1=(sqrt(G11),0), a_2 has positive second coordinate.
The elementary Gram formulas bound a_3's two coordinates below 0.079,
|a_2x a_2y|=|G12|<=0.00002 and |ell_2|<1.02. For q<=1,
`m(u)-m_3<1.30*0.452` on 0.1<=u<=0.3 and `m(u)-m_2>0.348`.
Thus Im D_A(u)>0.01 throughout that interval. The anchor block stays positive
definite, anchor pair slacks stay positive, and active probes are timelike
after all three anchors. These bounds can be checked with rational interval
arithmetic in the full instrument; they are not numerical results of this run.

Take disjoint candidate triples from O(log n) candidate events. Pilot moment
estimates use disjoint O(n) references. Accept the first triple in an inner
neighbourhood with half these tolerances. A still smaller neighbourhood has
positive physical probability. Bounded U-statistic concentration and a union
bound over O(log n) candidates imply: some triple is accepted with probability
tending to one, and every accepted triple is in the outer good set with
probability tending to one. Wave population moments converge uniformly to
the canonical background moments there, so the same conclusion holds under
the alternative. This argument does not supply a practical sample size: the
good-set probability and inverse-calibration constants may be very small.

The population expansion in the challenge, including angular integration,
now yields uniformly over these anchors and wave phase

\[
|E_\lambda[wH e^{-i\omega u}]|\ge c_q\lambda^2,
\quad |\omega-1/\lambda|\le1/8,                             \tag{11}
\]

for sufficiently small lambda. The O(1) frequency offset cannot cancel the
integral: |D_A|<0.05 on the same outer set, and centering the Fourier phase at
u=0.2 bounds its change by (0.1/8)*0.05<0.001, below the imaginary lower bound
0.01. Smooth terms involving
anchor phases vanish on integration by parts; one-point calibration errors
are o(lambda^2). Taking the magnitude respects the previous even-likelihood
symmetry and does not resurrect a forbidden signed, orientation-fixed score.

## 7. A scan with no wavelength input or unknown threshold constant

Here a simpler asymptotic construction replaces the optional Monte Carlo
calibration discussed in the challenge. It is a new derivation, not a
re-scoring of that population experiment.

Randomly relabel events and split them into L=ceil(16 log N) equal blocks,
discarding the remainder. In each block reserve 3 ceil(log m) candidates,
then divide the other m events between pilot references, estimation references
and probes. Put n equal to the smallest of these three list sizes. If the
block is too small or has no accepted anchor triple, its vote is zero.

For fixed K use the same deterministic grid in every block:

\[
\Omega_n=\{j/4:4\le j\le\lfloor4n^{K/(2K+2)}\rfloor\}.
\]

At each frequency let a block vote if its score magnitude exceeds

\[
t_n(\omega)=\log n\left[
1/(\omega\sqrt n)+1/n+\omega^K/n^{(K+1)/2}\right].            \tag{12}
\]

Reject the background if a strict majority votes at any common grid point.
All quantities in this definition are order counts, known q, deterministic
functions or independent relabeling randomness. There is no true wavelength,
phase, embedding, labelled generator anchor, or coordinate oracle input.

Under the null, (5) and Markov's inequality bound each fixed-frequency block
vote by C_K/(log n)^2 plus the vanishing pilot failure probability. Uniformly
on the grid it is eventually below 1/4. Blocks are independent, so Hoeffding
gives a majority probability at most exp(-L/8)<=N^-2. The grid has fewer than
4 sqrt(N) points. A union bound therefore makes the **whole scan** size tend
to zero. The logarithm in (12) eventually dominates the unknown fixed C_K;
this is an asymptotic statement, not a finite-N significance guarantee.

Under (1), 1/lambda is eventually strictly below the scan cap. A grid point
is within 1/8 of it. Equations (10)--(12) give signal/noise tending to infinity,
including the logarithmic threshold and block splitting. Its block vote tends
to one, and the majority rejects with probability tending to one. This proves
the proposed implication (1), conditional on the lemmas and regularized
calibration being realized as stated.

Finite precision must also be part of the instrument: compute the regularized
score with a certified absolute error, for example n^-4, and vote only when its
lower magnitude bound exceeds (12). This error is negligible compared with
the threshold throughout the scan range. Fixed double precision and a fixed
root tolerance do not implement this asymptotic requirement.

Before using the concrete geometric guards, freeze this document and the
source's `--guard-output` check. It will verify the Section 6 bounds using
Fraction interval arithmetic and the existing certified cosine series at
q=1, whose derivative bound controls every 0<q<=1. It must retain a failure
without enlarging the neighbourhood or weakening its target inequalities.
This algebraic control is not a noise-rate or finite-data-power test.

## 8. Remaining decisive work and prior art

The next required deliverable is the complete order-only count implementation
with fixed regularization, anchor guards, certified increasing precision and
the block scan. Its nontrivial positive controls must exercise the same data
path; a zero score on a small sample that never found anchors is not evidence
of detection power. Any numerical accuracy/power campaign must be sized and
frozen before sampling. The sharp lower bound remains a separate obstacle:
this proof does not exclude more informative geometry-fitting instruments
between powers 4/3 and 2.

Order-derived embeddings, distances and Gram algebra are established tools,
not our novelty claim: see [Henson (2006)](https://arxiv.org/abs/gr-qc/0601069)
and [Johnston (2025)](https://arxiv.org/abs/2502.09701). The possible contribution
would be the oscillatory statistical bound in this exact high-frequency
experiment and eventually a matching resolution law. A bounded search on
2026-09-24 did not retrieve this complete result, which is not a novelty proof.
For general U-statistic background see
[Gine, Latala and Zinn (2000)](https://arxiv.org/abs/math/0003228); no theorem
there has been substituted for the model-specific variation and coupling
arguments above.
