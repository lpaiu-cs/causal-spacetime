# Fixed-energy separation: a sufficient side, not the sharp resolution law

2026-09-24. Same-author derivation and a frozen deterministic check design.
The [full-order probe](full_order_resolution_probe.md) defines the experiment.
This document removes the exchange of weak-amplitude and high-frequency limits
from its pair-count sufficient-side heuristic. It does **not** close the gap to
the full-order impossibility bound or complete the research goal.

## 1. Fixed model and observation contract

Use the exact circular vacuum family in the probe, with fixed
`0 < q < pi^2`, axial tips at `u=-1/2,1/2`, fixed positive longitudinal slack,
and `lambda -> 0`. Constants below may depend on q and its positive-patch
margin. The comparison is its effective null-matter metric with the **same q**.
There are N independent uniform physical-volume events, conditional on N.
The observer receives the order, N, and the specified q/clock calibration.
It receives no coordinates, wavelength, polarization axis, or matched sample.
An unknown-q composite test is a separate problem.

A phase change in this circular family is a constant transverse rotation.
It fixes the axial tips and preserves the sampling measure. Thus the *entire
order law*, not just chain means, is independent of phase. No phase average
or phase inference is needed for this particular family.

Write `a(u)=cos(sqrt(q)*u)`, `f=a^-2`, `ell=t-s`, and

\[
K_0(s,t)=\frac{\sin(\sqrt q\,\ell)}{\sqrt q},\quad
m(s,t)=\int_s^t f(u)du=\frac{K_0(s,t)}{a(s)a(t)},\quad
I_\lambda(s,t)=\int_s^t f(u)e^{iu/\lambda}du.
\]

For `c=sqrt(1+4q lambda^2)`, divide the exact Jacobi kernel by c:

\[
\bar K_\lambda
=K_0\sqrt{1-\frac{4q\lambda^2}{1+4q\lambda^2}
                  \frac{|I_\lambda|^2}{m^2}}.
\tag{1}
\]

This constant division preserves every normalized chain probability. It
does not preserve the physical volume or remove microscopic geometry from
the full order. Let `delta K=bar K_lambda-K_0`.

## 2. The short-segment logarithm

For a continuously differentiable weight h and any fixed nonzero segment,

\[
\int_s^b h(t)\,\delta K(s,t)dt
=-4q\,h(s)\lambda^4\log(1/\lambda)+O(\lambda^4).
\tag{2}
\]

There is an analogous right-end formula. For a smooth two-variable weight,
integrating over the triangle gives the corresponding one-sided diagonal
measure: the leading term is
`-4q lambda^4 log(1/lambda) integral h(u,u)du`.
These are weighted integral statements, not pointwise delta functions.

Here is a remainder argument. Uniformly on the fixed positive patch,
`|I| <= C min(ell,lambda)`, `m` is comparable to ell, and

\[
I_\lambda=\frac{\lambda}{i}
 [f(t)e^{it/\lambda}-f(s)e^{is/\lambda}]+O(\lambda^2).
\tag{3}
\]

For `ell >= lambda`, expand the square root in (1). The leading term is

\[
-\frac{2q\lambda^4}{K_0(\ell)}
 \left[\frac{a(t)^2}{a(s)^2}+\frac{a(s)^2}{a(t)^2}
             -2\cos(\ell/\lambda)\right].
\tag{4}
\]

The bracket differs from `2[1-cos(ell/lambda)]` by `O(ell^2)`, and
`1/K_0-1/ell=O(ell)`. The integrated error from these replacements is
`O(lambda^4)`. The integration-by-parts remainder contributes at most
`O(lambda^5 log(1/lambda))`. The square-root remainder has integral
`O(lambda^6)`: bound it by
`C ell lambda^4 min(1,lambda^4/ell^4)` and split at ell=lambda.
The denominator `1+4q lambda^2` only adds a smaller-order correction.

The omitted strip `0 < ell < lambda` contributes `O(lambda^4)`, since
`|delta K| <= C lambda^2 ell` there. Finally,

\[
\int_\lambda^L\frac{1-\cos(\ell/\lambda)}{\ell}d\ell
=\log(1/\lambda)+O_L(1).
\]

Replacing h(t) by h(s) costs `O(lambda^4)` because h is Lipschitz.
This proves (2). For the triangular version the bound can be integrated
over s; the endpoint remainder `|log(b-s)|` is integrable. Fixed-segment
`delta K(s,t)=O(lambda^4)` has no logarithmic term.

## 3. Every fixed chain, then a strictly nonzero pair signal

Let `C_r^0=K_0^{*r}(-1/2,1/2)` and let `bar C_r` use `bar K_lambda`.
Expanding the r factors and applying (2) to each collapsed segment yields,
for every **fixed** integer r>=2,

\[
\boxed{\bar C_r=C_r^0-4qr C_{r-1}^0
                    \lambda^4\log(1/\lambda)+O_{q,r}(\lambda^4).}
\tag{5}
\]

There are r choices of the collapsed segment, each leaving `C_(r-1)^0`.
Terms with two or more perturbations are smaller: the integrated absolute
perturbation is `O(lambda^4 log(1/lambda))`, and the remaining bounded
convolution kernels give products of these bounds. For fixed r at least one
segment has macroscopic length, so no all-collapsed constraint is inserted.
There is no assertion of uniformity as r grows with 1/lambda or N.

The earlier exact hierarchy gives
`P_k=k! 2^k/(2k)! K^(k-1) C_(k+1)/C_2^k`. Consequently

\[
P_{k,\lambda}=P_{k,0}\left\{1+4q\left[
 \frac{2kK_0}{C_2^0}-\frac{(k+1)C_k^0}{C_{k+1}^0}\right]
 \lambda^4\log(1/\lambda)\right\}+O_{q,k}(\lambda^4).
\tag{6}
\]

For pairs this simplifies substantially. Put
`p0=P_(2,0)=K_0 C_3^0/[3(C_2^0)^2]` and

\[
\boxed{p_\lambda-p_0=-A(q)\lambda^4\log(1/\lambda)+O_q(\lambda^4),
\qquad A(q)=4q\frac{K_0}{C_2^0}(1-4p_0)>0.}
\tag{7}
\]

The sign is analytic, not a fitted coefficient. With `k=sqrt(q) in (0,pi)`,
`D=sin k-k cos k`, and `E=(3-k^2)sin k-3k cos k`,

\[
K_0=\sin k/k,\quad C_2^0=D/(2k^3),\quad
C_3^0=E/(8k^5),\quad p_0=\sin k E/(6D^2).
\]

Both kernels and convolutions are positive. Also
`6D^2-4 sin(k) E=2[3k^2-(3+k^2)sin^2 k]>0`.
Indeed Euler's sine product implies
`(sin k/k)^2 <= exp(-k^2/3) < 1/(1+k^2/3)`.
Thus `p0<1/4`, giving `A(q)>0` throughout the stated patch range.
At weak energy `A(q)/q -> 72/5`. This agrees with the separately derived
weak-wave transfer's logarithmic term; (7) itself holds at **fixed q**.

## 4. An intrinsic, frequency-independent sufficient test

Let `U_N` be the fraction of unordered comparable pairs. It is a bounded
order-two U-statistic, not a binomial average of independent edges. For
`m_N=floor(N/2)` and t>0, the standard Hoeffding bound gives

\[
\Pr(U_N-\mathbb E U_N\ge t),\quad
\Pr(U_N-\mathbb E U_N\le-t)\ \le\ e^{-2m_Nt^2}.
\tag{8}
\]

One proof averages the means of disjoint pairs over all permutations.
Convexity bounds the U-statistic moment generating function by that of
one such independent-pair mean; the bounded-variable Hoeffding lemma then
gives (8). The tool is established, not a new concentration theorem.
See [Hoeffding (1963)](https://doi.org/10.1080/01621459.1963.10500830), section 5.

Fix the one-sided false-positive bound **alpha=0.025** before data. Reject
the effective-matter model when

\[
U_N<p_0-\sqrt{\frac{\log(1/\alpha)}{2m_N}}.
\tag{9}
\]

This needs neither lambda nor phase. Its false-positive probability is at most
alpha for every N>=2. Equations (7)-(8) imply power tending to one when

\[
\boxed{N\lambda^8\log^2(1/\lambda)\longrightarrow\infty.}
\tag{10}
\]

This is a detection test, not an accuracy gate. No estimated confidence
interval or Gaussian approximation enters it. At an exact population gap
`Delta=p0-p_lambda>0`, a conservative **90% power** sufficient size is

\[
m_N\ge\frac{[\sqrt{\log 40}+\sqrt{\log 10}]^2}{2\Delta^2}.
\tag{11}
\]

For simple-versus-simple total variation, a threshold at the exact midpoint
has both errors tending to zero under (10); hence the full-order TV tends to
one in that region. This midpoint argument can depend on lambda. It must not
be confused with a single wavelength-free test having both errors vanish
uniformly over every arbitrarily slowly separating sequence. Test (9) has
fixed size alpha and power tending to one instead.

Combining this sufficient region with the probe yields a **bracket**, not a
boundary: `N lambda -> 0` is indistinguishable, whereas (10) is distinguishable.
No assertion is made about the large gap between them. Adding finitely many
fixed chain means retains the same order of the population signal; proving
an optimal noise rate or a full-order lower bound requires further work.

## 5. Frozen deterministic falsification check

Source: `experiments/positive_control/order_wave_resolution.py`.
This design and source are committed together before the diagnostic runs.
Outputs bind the commit and file SHA256 values; a dirty tree is rejected.
No GPU, random training sample, publication submission, or finite-N campaign
is involved.

The solver integrates the complex oscillatory primitive directly with the
analytic scale factor. It computes changes in C2, C3, and P2 algebraically,
avoiding subtraction of two almost equal absolute probabilities. Controls
compare the kernel and P2 with the earlier independent five-state Rosen ODE,
and compare the q->0 coefficient with the existing transfer formula.

Primary grid: q in `{0.25,1,4}` and r=1/lambda in
`{32,64,128,256,512,1024}`. The Gauss rule has
`max(128,ceil(r/2)+64)` nodes per dimension. The q=1, r=128 and r=1024 cases
also use twice that order and tighter ODE tolerances. Record the coefficient
`-Delta P2/[lambda^4 log r]`, its predicted limit A(q), and the bounded-scale
remainder `Delta P2/lambda^4+A(q)log r`. The last quantity need not converge:
an O(lambda^4) endpoint oscillation is permitted by the theorem.

Do not demand rapid convergence in `1/log r`, fit a power and relabel it a
proof, or treat a refinement difference as a certified numerical error bound.
If the sign or remainder behavior contradicts the derivation, retain all rows
and revisit the proof. These numerical results do not validate finite-sample
power, a sharp information boundary, or novelty. The explicit gaps in the
full-order probe remain the active research task.
