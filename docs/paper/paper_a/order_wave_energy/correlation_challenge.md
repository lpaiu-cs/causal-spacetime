# A risky target: causal-order spectroscopy, and its first principled pivot

2026-09-24. Derivation and diagnostic design, before the new numerical run.
The user explicitly prefers a discovery-sized conclusion even if it probably
fails, and permits moving to another direction after failure. The previous
526-solve forward result is not treated as sufficient scientific payoff.

## The strong claim to attack

**Can scale-dependent causal relation counts determine the hidden gravitational
wave's entire energy spectrum, without fitting a supplied waveform?**

The unrestricted finite-window form of this claim already has a mathematical
counterexample, below. We therefore move to a precise related discovery target:

> For weak stationary plane waves with both polarizations, does the *single*
> scale-dependent pair relation probability determine the complete derivative-
> strain structure function inside the observed lag range, without supplying
> the spectrum, mode count, polarization, or energy to the inverse?

This is a nonlocal correlation target, not another calibration of mean Ricci
or a finite list of fitted Fourier peaks. Its scope is still classical plane
waves and the leading weak-amplitude coefficient. Novelty is undecided.

## 1. Observation contract and physical quantity

Let a weak transverse strain be
\[
B(u)=\epsilon\begin{pmatrix}b_+(u)&b_\times(u)\\
b_\times(u)&-b_+(u)\end{pmatrix},\qquad S=e^{2B}.
\]
The Rosen scale a solves the exact vacuum equation in the earlier design.
Fix the freely falling central observer and its proper-time unit, with axial
diamonds of duration L. Assume a stationary phase ensemble with enough mean
square derivatives for the operations below. Define the amplitude-coefficient
correlations (the factors epsilon^2 are suppressed in these quantities)
\[
C(L)=\langle b_+(u)b_+(u+L)+b_\times(u)b_\times(u+L)\rangle,
\quad R(L)=-C''(L),\quad q=R(0),\quad S_R(L)=q-R(L).
\]
Here R is the two-point correlation of **derivatives of strain**, whose zero-
lag value enters the usual averaged wave energy. It is not the fourth-order
correlation of local energy densities. S_R is half the mean squared difference
of the derivative-strain vectors, and is nonnegative for a stationary process.
The trace over polarizations does not recover helicity or separate polarization
components. An absolute frequency unit needs the stated observer clock.

The inverse receives only L and the population relation fraction P2(L), averaged
over the stationary ensemble. In the mathematical identity it receives the
coefficient of epsilon^2 in that curve; a finite nonzero amplitude has higher-
order corrections. How to extract it from a single unknown-amplitude, unknown-
phase finite poset is deliberately **not** assumed solved.

## 2. Exact identity at quadratic order

Let D(L) be the phase mean of the quadratic perturbation of the Jacobi kernel,
and let Y(L) be the quadratic coefficient of log(10 P2(L)). The earlier
scalar-strain derivation extends by taking a trace: at quadratic order the
two polarization contributions add, including for noncommuting B(u).
\[
D(L)=2LC(0)-\frac4L\int_0^L(L-s)C(s)ds-\frac{qL^3}{6},
\tag{1}
\]
\[
Y(L)=\frac{D(L)}L+\frac{60}{L^5}\int_0^L(L-s)^3D(s)ds
-\frac{24}{L^3}\int_0^L(L-s)D(s)ds.
\tag{2}
\]
Define e(L) as the quadratic coefficient of
\(105[1-10P_2(L)]/(4L^2)\), so
\(Y=-4L^2e/105\). Let theta = L d/dL and d=D/L. Evaluating the Volterra
operator on powers, or differentiating its integrals, gives
\[
Y=A(\theta)d,\qquad
A(z)=\frac{z(z-1)(z^2+15z+62)}{(z+2)(z+3)(z+4)(z+5)}.
\tag{3}
\]
This rational identity follows from
\(A(z)=1+360/[(z+2)(z+3)(z+4)(z+5)]-24/[(z+2)(z+3)]\).
The integral realization fixes its boundary conditions at zero. Vacuum and
smooth stationary covariance imply d=O(L^4). The homogeneous solutions
1, L, and \(L^{(-15\pm i\sqrt{23})/2}\) do not satisfy that condition,
so they cannot be inserted as invisible physical modes.

Four derivatives of (1) remove its affine covariance ambiguity:
\[
(LD)''''=4[R(L)-q]=-4S_R(L).
\tag{4}
\]
Combining (3)--(4) yields the candidate inversion identity
\[
\boxed{(\theta^2+19\theta+96)S_R
=\frac{(\theta+3)(\theta+4)^2(\theta+5)(\theta+6)(\theta+7)}{105}\,e.}
\tag{5}
\]
The solution bounded at zero is unique: the two homogeneous roots on the
left are \((-19\pm i\sqrt{23})/2\), both with negative real part.
This proves the claimed injectivity of the map to S_R within the stated
quadratic-order regularity class, if (1)--(2) hold. It does not prove stability
in the raw observation norm.

An explicit inverse requiring four rather than six raw derivatives is useful.
Set a=19/2, b=sqrt(23)/2. Polynomial division gives
\[
S_R(L)=\frac1{105}\left[
(\theta^4+10\theta^3+59\theta^2+74\theta+384)e(L)
-\int_0^1 x^{a-1}\left(864\cos(b\log(1/x))
+\frac{18576}{b}\sin(b\log(1/x))\right)e(Lx)dx\right].
\tag{6}
\]
All values required by the integral lie inside the observed window. A positive
minimum observed size leaves additional extrapolation uncertainty; replacing
it by an analytic polynomial is an assumption of the numerical prototype.

Four derivatives are a genuine noise penalty. The leading symbol grows as
|k|^4 in log-scale frequency; injectivity does not make this inverse a useful
finite-event instrument. The numerical challenge includes this failure mode.

## 3. Why arbitrary full-spectrum recovery over a finite window fails

Let Lmax be finite. Choose a nonzero real even smooth bump b(t), supported in
|t| in (Lmax+1,Lmax+2), and let h(omega) be its inverse cosine transform with
\(b(t)=\int_0^\infty h(\omega)\cos(\omega t)d\omega\).
It is real, signed and rapidly decreasing. Define
\[
\delta\rho(\omega)=\omega^2h(\omega),\quad
\rho_\pm(\omega)=|\delta\rho(\omega)|\pm\delta\rho(\omega).
\tag{7}
\]
These are distinct nonnegative derivative-energy spectra. Their total
energies are equal because \(\int\delta\rho=-b''(0)=0\). Their strain
spectra rho/omega^2 are integrable, including at zero, and all required
derivative moments exist. Their strain covariance difference is 2b(L), hence
both covariance and derivative-covariance differences vanish on [0,Lmax].
Equations (1)--(2) therefore give **exactly identical quadratic-order pair-
response curves throughout the available window**, despite distinct spectra.
Stationary Gaussian processes with these spectra provide realizations.

This defeats nonparametric full-spectrum uniqueness from this observable on a
finite window. It does not claim identical exact nonlinear metrics or full
finite-poset laws. A known finite spectral band changes the uniqueness problem:
analytic continuation can restore formal uniqueness while remaining unstable.
The counterexample is a standard Fourier-window mechanism applied here, not
itself a claim of new Fourier theory. It is a reason to drop that strong target,
not to relabel failure as a successful spectrum instrument.

## 4. Frozen second challenge

Freeze this document and source before calculating the test curves. Record the
reachable commit, source hashes, inputs and all outputs. The following are
deterministic diagnostics, not a statistical accuracy pass/fail campaign.

- Reuse the prior exact Rosen ODE and chain integrals through a shared
  coefficient entry point. Preserve the earlier result bytes and provenance.
- Verify (3) and the polynomial division with exact rational/integer arithmetic.
- Numerical inverse: 65 Chebyshev size nodes in [0,1]; fit
  `0.1-P2 = L^4 sum c_j T_j(2L^2-1)`. Recover e and its theta derivatives from
  that fit and evaluate (6) with 128 Gaussian nodes. The inverse receives only
  the scale nodes and P2 values, not the spectrum, q or target structure function.
- Primary polynomial degree 18; degrees 12 and 24 are retained as prespecified
  discretization challenges. Evaluate at 29 lags in [0.15,0.85]. No result-based
  choice of degree or lag range is permitted.
- Analytic input controls use the earlier F integral for three fixed spectra,
  comparing the inverse with `q sum weights*(1-cos(omega*L))`. They check
  the derivation/implementation; exact ODE cases prevent relying only on the
  same forward approximation used to derive the inverse.
- Exact vacuum cases: physical q=0.0025 and 0.00125, 24 uniform phase shifts,
  base angular frequency 2 and the following harmonics. All share an initial
  averaged scale at u=-1/2. One ODE per phase covers all smaller diamonds.

| Case | Harmonics | Energy weights | Ellipticities | Polarization rotations |
|---|---|---|---|---|
| separated | 1,4,9 | .25,.35,.40 | 0,0,0 | 0,0,0 |
| rotating | 3,7,11 | .20,.35,.45 | .7,-.4,1 | .2,.8,1.3 |
| broad | 1,2,3,4,5,6,7,8 | normalized 1/n | alternating +.6,-.6 | .17*n |

Mode phase offsets are .31*n. For each mode, strain is a rotated vector
`sqrt(2*q*weight/(1+eta^2))/omega * (sin(z),eta*cos(z))`. Its mean derivative
energy is q*weight. The shared phase enters as n*phi, making 24 samples exact
for the quadratic phase average but not for all higher orders. S=e^(2B) and
Q=tr[(S^-1 S')^2]/8 are used **exactly**, not replaced by the quadratic sum.

- Use 64 nodes per chain quadrature for main curves and ODE tolerances as in
  the previous tighter solver. Repeat rotating q=.00125 with 48 phases and
  128 chain nodes. Keep that result alongside the main curve.
- Report normalized RMSE and maximum error of S_R/q for every degree and
  amplitude; these descriptive errors are not confidence intervals. Prediction:
  approximation residuals should shrink with amplitude until numerical
  differentiation dominates. Failure to do so stops the inverse claim.
- Determine the linear observation-to-S_R matrix of the frozen inverse. Report
  its noise amplification. For **independent diamonds at each scale**, use
  the exact flat P2 variance as a weak-field planning approximation. A 95%
  two-sided interval with 90% power inside an absolute tolerance .25*q uses
  `(z_.975+z_.95)` in the price. Do not quote break-even or treat the resulting
  estimate as a proved lower bound for other statistics or correlated windows.

## 5. Interpretation and literature boundary

If the identity and no-oracle inverse survive, the precise candidate result is
that one multiscale order probability determines a nonlocal derivative-strain
structure function in a weak stationary plane wave. This is neither a new
gravitational energy definition nor a first detection of Weyl curvature.
If finite-data conditioning is disastrous, this prototype is stopped as a
practical reconstruction route; more fits do not repair an information problem.

Background: [Roy--Sinha--Surya](https://arxiv.org/abs/1212.0631) already relate
chains to local Ricci, [Harte](https://arxiv.org/abs/1502.03658) treats nonlinear
plane-wave optics, and [Huneau--Luk](https://arxiv.org/abs/2404.07659) review
high-frequency backreaction. [Ahmed--Bombelli's 2026 APS abstract](https://meetings-archive.aps.org/smt/2026/apr-r93/1/)
reports causal-set Weyl detection through geodesic deviation; the retrieved
abstract does not provide an energy-correlation inversion formula. Searches
for causal-set wave spectrum/correlation and causal-diamond spectroscopy have
not retrieved the same identity. This limited non-retrieval is not evidence of
priority, and no independent expert has audited this derivation.
