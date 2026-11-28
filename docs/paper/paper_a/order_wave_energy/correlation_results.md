# Discovery search: abandon unrestricted spectroscopy; test nonlocal curvature information

2026-09-24. The user authorized a high-risk direction with a substantial possible
conclusion, and switching direction when it fails. This record preserves the
failures rather than treating every narrower result as success. No submission
or discovery claim follows from these calculations.

## Decision

1. **Stop the unrestricted full-spectrum claim.** On a finite window, distinct
   nonnegative spectra can have the same entire quadratic-order pair-response
   curve. The counterexample is in the [frozen derivation](correlation_challenge.md).
2. **Retain a mathematical candidate:** multiscale causal pair probabilities
   determine a nonlocal derivative-strain structure function at quadratic order
   for weak stationary plane waves. This survived a no-spectrum-input continuum
   calculation with two polarizations. It is not a finite-event recovery result.
3. **Stop the direct pointwise numerical inverse as a practical route.** Its
   fourth-derivative noise penalty leads to prohibitive sample prices under the
   tested protocol. More precision runs on that estimator are not the next step.

The next discovery-sized question is whether losing the wave signal is a
consequence of compressing the order into pair fractions, or a limitation of
the full finite causal order. A successful answer would give a quantitative
event-count versus wavelength distinguishability law, and identify an intrinsic
statistic attaining it. That question is not answered by the present forward
formula or its inverse. No full-order rate experiment has been run here.
The first analytic attack on it is saved in
[full_order_resolution_probe.md](full_order_resolution_probe.md): for the exact
circular family, a coupling argument gives a candidate bound proportional to
N^2*lambda^2 and hence indistinguishability when N*lambda tends to zero. It is
not a matching rate or an independently audited theorem.

## 1. The conditional mathematical result

Write the two strain components as b_plus and b_cross. For a fixed central
observer with a calibrated clock, set

\[
C(L)=\langle\mathbf b(u)\cdot\mathbf b(u+L)\rangle,\quad
R(L)=-C''(L),\quad q=R(0),\quad S_R(L)=q-R(L).
\]

S_R is half the mean squared increment of the derivative-strain vector. It
measures how wave changes at two times differ. Calling it an energy-density
autocorrelation would be wrong: the latter is a fourth-order field statistic.

Let e be the quadratic amplitude coefficient of
\(105[1-10P_2(L)]/(4L^2)\), and theta=L*d/dL. The derived identity is

\[
\boxed{
(\theta^2+19\theta+96)S_R
=\frac{(\theta+3)(\theta+4)^2(\theta+5)(\theta+6)(\theta+7)}{105}\,e.
}
\]

Its solution bounded at zero is unique. The explicit integral inverse requires
up to four derivatives of e, rather than a fitted list of spectral lines. The
derivation uses the prior exact chain hierarchy and a trace expansion of the
Rosen metric; exact rational algebra and separate matrix differentiation check
the new factors and the two-polarization vacuum generator. This same-agent
derivation is not an independent mathematical audit.

There is a useful physical reading. In the leading parallel transverse frame,
the vacuum plane-wave tidal profile is
\(H=\epsilon\begin{pmatrix}b_+''&b_\times''\\b_\times''&-b_+''\end{pmatrix}
+O(\epsilon^2)\). Thus its trace correlation coefficient is

\[
W(L)=\tfrac12\langle\operatorname{tr}(H(u)H(u+L))\rangle/\epsilon^2
=C''''(L)=S_R''(L).
\]

Equivalently S_R is the double integral
\(\int_0^L(L-s)W(s)ds\). The candidate therefore encodes a nonlocal **tidal
curvature correlation**, integrated twice, rather than merely a local curvature
norm. A raw W reconstruction would require two more derivatives and has not
been executed. Comparing tidal tensors at separated points requires the stated
transported frame; no frame-free tensor correlation is claimed.

The assumptions are consequential: weak-amplitude coefficient, stationary
ensemble or a justified phase average, a positive plane-wave patch, and a
chosen observer clock. There is no inference of absolute units from a bare
poset and no recovery of polarization-resolved spectra or helicity.

## 2. What the new calculation actually did

The inverse function receives only the scale nodes and 0.1-P2 at those nodes.
The generator's energy, Fourier frequencies, weights, polarization rotations,
and ellipticities are not arguments to it. Those values are read by the separate
evaluation function after inversion. These are preselected mathematical test
cases, not an independently administered blind campaign or sampled posets.

The design and source were committed at
`fe685a8747ce37979e664406acae3eac151b1b3c` before execution. The CPU run took
305.97 seconds and completed:

- Three analytic-input controls of the derived inverse.
- Six exact-vacuum cases: three spectra/polarizations at two amplitudes, each
  with 24 phases and 64 nonzero size nodes. These use 144 ODE solutions and
  9,216 continuum pair-probability evaluations.
- One prespecified refinement with 48 phases and 128 chain quadrature nodes:
  48 ODE solutions and 3,072 pair-probability evaluations.
- All three prespecified polynomial degrees, including the primary degree 18,
  and their linear noise amplification matrices.

All 13 raw rows, including all 192 phase curves, were retained. No GPU was used.

The errors below are divided by q, **not** by S_R at each lag. They are
deterministic discrepancies against a quadratic-order target, not statistical
confidence intervals. Primary degree 18, 29 lags in [0.15,0.85]:

| Case | q | RMSE/q | Maximum error/q |
|---|---:|---:|---:|
| separated spectral lines | 0.0025 | 0.03428% | 0.08462% |
| separated spectral lines | 0.00125 | 0.01714% | 0.04232% |
| rotating / elliptical polarizations | 0.0025 | 0.04623% | 0.09314% |
| rotating / elliptical polarizations | 0.00125 | 0.02311% | 0.04656% |
| broad eight-mode profile | 0.0025 | 0.09832% | 0.14926% |
| broad eight-mode profile | 0.00125 | 0.04913% | 0.07456% |

Halving q multiplies RMSE by 0.49967--0.49994, consistent with the leading
higher-order correction to the normalized response. This supports the tested
weak-wave identity; it does not prove an all-profile remainder bound.

The largest analytic-control error/q with degree 18 was 6.741e-8. Degree 12
fails to resolve the highest-frequency rotating case as well as degrees 18/24;
its exact-input maximum error/q is about 1.76% at q=.00125. It was retained,
not selected away. Degree 24 does not remove the physical amplitude correction.

Doubling the phase count and chain quadrature for the rotating q=.00125 case
changed P2 by at most 1.122e-12 but changed the recovered S_R/q by up to
3.260e-5. The refined maximum error/q was 0.04362%. This is itself a visible
example of numerical differentiation amplifying tiny input differences.

## 3. Why the practical inverse is stopped

The primary inverse's linear map from measured P2 values to S_R has row norms
between 2.291e4 and 1.571e7 over the evaluated lags. For independent diamonds
at each size, the flat-limit relation-fraction standard deviation is
approximately 0.2/sqrt(N), where N is the number of events **per size**.

The frozen planning calculation uses a pointwise 95% two-sided interval with
90% probability of lying wholly within an absolute error band 0.25*q at
q=.00125, neglecting approximation bias. It gives approximately

\[
2.79\times10^{15}\ \text{to}\ 1.31\times10^{21}
\quad\text{events per measured size, depending on the target lag.}
\]

This is a normal/weak-field planning approximation for this particular linear
inverse and independence protocol. It is **not** a proved minimax lower bound,
not a simultaneous confidence band over all lags, and not a bound on estimators
using the full causal order or correlated windows. Even the favorable end is
far beyond a sensible campaign. No such campaign was launched.

The practical claim is therefore stopped. Quietly replacing pointwise recovery
by a strongly regularized few-mode fit would change the question and introduce
spectral assumptions; it is not counted as success of the original inverse.

## 4. A further interpretation check: even q needs a finite-window assumption

This extension was derived during interpretation after the frozen run; it is
not one of its preregistered numerical tests.

The spectrum counterexample in the design preserves q. A related construction
shows why S_R on a finite interval need not determine q. Choose an even smooth
compactly supported covariance perturbation b(t) equal to c*t^2 throughout
|t|<=Lmax, with c nonzero. Let h be its inverse cosine transform and again set
delta_rho=omega^2*h and rho_plus/minus=|delta_rho| +/- delta_rho. These are
nonnegative spectra of mean-square smooth processes, with finite strain
variance and moments.
Their total derivative energies differ by -2*b''(0)=-4c, but

\[
\Delta S_R(L)=-2b''(0)+2b''(L)=0\qquad(0\le L\le L_{\max}).
\]

Substitution into the D formula also gives exactly zero pair-response
difference at quadratic order. Thus a finite-window response identifies the
increment structure, not an unrestricted total derivative energy. A known
short-wave band or a measured decorrelation plateau can supply the additional
assumption. The construction does not refute that constrained high-frequency
regime and does not identify derivative energy at arbitrary wavelength with
a generally covariant local gravitational stress tensor.

## 5. Novelty comparison and next falsifiable direction

The following primary sources delimit the claim:

- [Roy--Sinha--Surya](https://arxiv.org/abs/1212.0631): chain statistics already
  estimate local Ricci quantities. Our target is a nonlocal wave correlation.
- [Harte](https://arxiv.org/abs/1502.03658) and
  [Huneau--Luk](https://arxiv.org/abs/2404.07659): nonlinear plane-wave optics
  and effective high-frequency backreaction are established surrounding theory.
- [Ahmed--Bombelli, APS 2026](https://meetings-archive.aps.org/smt/2026/apr-r93/1/):
  causal-set Weyl detection via geodesic deviation is already an active project.
- [Alfyorov--Shnyukov, author-uploaded 2026 preprint](https://www.researchgate.net/publication/403462568_Weyl_curvature_from_the_Hasse_diagram_a_parameter-free_bridge_formula_for_causal_sets):
  a conditional local electric-Weyl-squared proposal uses paired flat/curved
  realizations and coordinate-based strata. The accessed definition is not a
  single-poset, multiscale correlation inverse. Its claimed accuracy and
  derivation have not been independently audited here. SSRN access failed;
  the author-uploaded full text was used for this limited comparison.

No matching formula was retrieved in the bounded spectrum/correlation searches.
That is not evidence of priority. The current result remains a conditional
identity with numerical support and an impractical direct inverse, not a
declared new physical discovery.

The next target should attack the compression bottleneck instead of adding
more spectra to this table: for exact rapidly rotating vacuum waves versus
their effective null-matter geometry, determine how the **whole order law**
separates when event count grows jointly with frequency. Fixed N is already
blocked; raw pair fractions provide only an inefficient upper construction.
Success requires a matching rate or a demonstrably better intrinsic statistic
on withheld frequencies/phases, with no matched flat coordinates. The first
step is a coupling/information bound and an affordable separation protocol,
not a large simulation or a new publication package.

## Evidence and reproduction

Original output:
`E:/lab/causal-spacetime/outputs/order_wave_energy/correlation_challenge_20260924`.

Portable copies: [manifest](correlation_evidence/manifest.json),
[raw rows](correlation_evidence/rows.jsonl),
[summary](correlation_evidence/summary.json).

Raw SHA-256:
`ac0d0e1e99b5d09f7539e37294390cb1d3370aa55fb2d8f5a74589807a885c0d`.
The manifest binds the design, new runner and reused generator to the reachable
pre-execution commit. The summary verifier checks those Git bytes, the raw hash,
completion status and row count; Windows manifest paths are normalized only
for Git lookup. The recorded data are not rewritten.

```powershell
rtk proxy E:/lab/causal-spacetime/.venv/Scripts/python.exe experiments/positive_control/order_wave_correlation.py --self-check
rtk proxy E:/lab/causal-spacetime/.venv/Scripts/python.exe experiments/positive_control/order_wave_correlation.py --output NEW_OUTPUT_DIRECTORY
rtk proxy E:/lab/causal-spacetime/.venv/Scripts/python.exe experiments/positive_control/summarize_order_wave_correlation.py NEW_OUTPUT_DIRECTORY
```

The original energy-run analytical controls also passed after extracting the
shared Rosen coefficient entry point. The changed Python files pass Ruff.
These checks validate this implementation; they do not constitute independent
physics review or a statistical reconstruction campaign.
