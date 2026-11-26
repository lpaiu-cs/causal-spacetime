# Is exponent 2 sharp? Formal statement, attack routes, stop rules

2026-11-24. Written by the audit orchestrator after the independent audit
recorded in [AUDIT.md](AUDIT.md). Nothing in this file is a result. Section 3
is heuristic accounting by one author and has not been reviewed by anyone.
No design in Section 5 is frozen and none has been run.

## 1. The question

Exact circular vacuum family, fixed known `0<q<=1`, fixed-N iid physical-volume
sampling, order-only observation, `N = lambda^-alpha`. Define

- `alpha_low  = sup{alpha : TV(P_{N,lambda}, P_{N,0}) -> 0}`,
- `alpha_high = inf{alpha : TV(P_{N,lambda}, P_{N,0}) -> 1}`.

Status after the audit (see AUDIT.md for the conditions attached to each
number): all three arguments were found correct, so
`4/3 <= alpha_low <= alpha_high <= 2`. The upper 2 rests on the audited
draft, whose regularised statistic still has to be written into its text;
the upper bound with a fully written proof is 3.

*Note added 2026-11-26.* The regularised statistic is now written into a
revision of the draft, [../gram_fourier_noise_r2.md](../gram_fourier_noise_r2.md).
An AI referee in a fresh context graded it CORRECT WITH MINOR FIXES, and the
fixes are applied ([AUDIT.md](AUDIT.md) Section 10). The upper 2 is not yet
written out in every step. Its index count, its tip estimates and the
expansion behind (11) are still at outline level, and the referee reports
supply them. The numerical work on U1 and U2 is stopped, with its readings
left as recorded ([STATUS_AND_LEVERS.md](STATUS_AND_LEVERS.md) Section 7).
Later the same day the audit's corrections to the lower bound were written
into [../quantized_order_bound_r2.md](../quantized_order_bound_r2.md) and
reviewed ([AUDIT.md](AUDIT.md) Section 11). Its one-point bound has since
been proved in [../one_point_strip_bound.md](../one_point_strip_bound.md)
(AUDIT.md Section 12). So has the wave's exact causal relation, in
[../wave_causal_relation.md](../wave_causal_relation.md) (AUDIT.md
Section 13). The lower end now rests only on the model's definitions.
Later still, the upper end's outline steps were written out in
[../gram_fourier_lemmas.md](../gram_fourier_lemmas.md) (AUDIT.md
Section 14), so both ends are now written as proofs with AI-level review.
Question Q is unchanged.

**Question Q.** Is `alpha_low = alpha_high = 2`? Equivalently: can the
necessary condition be raised from 4/3 to 2, or is there an order-only test
below 2?

In units of the discreteness length `l = N^(-1/4)` of the unit diamond,
`alpha = 2` is `lambda ~ l^2` and `alpha = 4/3` is `lambda ~ l^3`. Every
number in the bracket lies below `alpha = 4` (`lambda ~ l`): in this
plane-symmetric family all N events constrain one null coordinate, so the
order resolves wavelengths far below the discreteness length. That statement
is about this family with known q. It is not a claim about generic waves.

## 2. Why each endpoint is a property of a method, not yet of the problem

**Why 2 appears.** The Gram-Fourier statistic is first order in the wave
amplitude and uses one global Fourier phase `exp(-i omega u_hat)`. Two
things in it give the same number.

- Its leading noise is the shared-reference term `1/(omega sqrt(n))`. At
  `omega = 1/lambda` that is `lambda/sqrt(n)`, and the signal `lambda^2`
  exceeds it exactly when `lambda sqrt(n) -> infinity`, that is `alpha > 2`.
  The audit found that the wave-side error terms sit at the same order with
  no slack, and that the Taylor remainder is not what binds: a scan capped at
  `sqrt(n)` would still be controlled. Raising K only removes the `2/K`.
- It needs the order-reconstructed null coordinate to carry the wave phase
  across the whole probe window, that is `|u_hat - u| << lambda`. A
  coordinate built from reference counts has error `n^(-1/2)`.

So 2 is the threshold for *global phase lock by count statistics*, and it is
intrinsic to this statistic.

**Why 4/3 appears.** The impossibility proof gives the observer the cube
index of every event at side `h ~ N^(-1-1/K)`. In that stronger experiment
every event's position is known to about `1/N`, which is far below `lambda`
whenever `alpha > 1`. The stronger observer therefore knows every event's
wave phase for free. So 4/3 is the threshold for *an observer who is handed
the phase*.

The interval `(4/3, 2)` is therefore one question: **how much of the wave
phase, and of the fine null-margin information, does the order itself
carry?**

A constraint already in the record narrows the form of any answer: by the
quarter-turn symmetry in `cone_constraint_challenge.md` the order likelihood
is even in the wave amplitude, so the linear score vanishes and the leading
order-only information is quadratic. A test below 2 must be an even, joint
statistic. The modulus of the Gram-Fourier coefficient is one such statistic.

## 3. Heuristic accounting (unreviewed; orientation only)

Write `delta = lambda^2` for the displacement of a causal boundary. Two
structural unknowns decide the answer.

- **U1, margin resolution.** For a near-null pair whose relation could be
  flipped without violating transitivity (a link or a critical pair), how
  sharply does the rest of the order predict that relation? Call the width of
  the null-margin window over which it is unpredictable `w`.
  Cell-limited: `w ~ 1/N`, the spacing of the N other light cones along a
  line. Count-limited: `w ~ N^(-1/2)`, the accuracy of any count-based
  coordinate. The argument for count-limited is that the *sequence* of cone
  crossings near an event is not witnessed by other events at scale `1/N`
  (the lens between two cones holds about `1/N` events there). The argument
  for cell-limited is that the impossibility proof finds exactly `N^3 delta^2`
  of information once cells are revealed.
- **U2, phase transport.** How well does the order fix the difference of the
  null coordinate `u` between two events at separation `d`? Count increments
  with shared references give `sqrt(d/N)` (the nested-kernel bound in the
  shared-reference proof). If nothing does better, phase is coherent only
  over `d_c ~ N lambda^2`, which is below 1 exactly when `alpha < 2`.

Accounting: items with relative tilt `delta / w`, second-order information
`(number of coherent item pairs) x (tilt)^4`.

| U1 | U2 | accounting | exponent |
|---|---|---|---|
| cell-limited | phase known | `N (N delta)^2` | 4/3 (reproduces the proved bound) |
| cell-limited | `d_c = N lambda^2` | `N^2 d_c (N delta)^4` | 10/7 |
| cell-limited | none, per event only | `N (N delta)^4` | 8/5 |
| count-limited | phase known | `N^(3/2) (sqrt(N) delta)^2` | 8/5 |
| count-limited | none, per event only | `sqrt(N) (N^(3/4) delta)^2` | 2 |
| count-limited | windowed Fourier filter | needs `1/(N lambda^2) << d << N lambda^2` | 2 |

Reading of the table, with its limits stated:

1. In the count-limited scenario three unrelated accountings return the
   same number, 2: global phase lock, a windowed filter that must be longer
   than `1/(N lambda^2)` to suppress the reference staircase noise and
   shorter than `N lambda^2` to stay coherent, and the phase-free per-event
   quadrupole of an event's roughly `sqrt(N)` near-null partners. That
   coincidence is the best available reason to think 2 might be sharp.
2. Every cell-limited row is below 2. If the `1/N`-scale information is
   accessible to the order, 2 is not sharp, and the attaining statistic is
   not a smooth functional of counts. It may also not be efficiently
   computable; a statistical-computational gap is possible and the question
   posed here is the statistical one.
3. These are scalings from an independent-items picture. No row is a
   theorem, the rows are not exhaustive, and two of them (10/7, 8/5) rest on
   a tilt per item that no explicit statistic in this folder extracts.

So Q reduces to U1, with U2 deciding the value if U1 is cell-limited.

## 4. Attack routes and stop rules

Each route has a result that would end it. A stopped route is recorded as a
result. It is not converted into a smaller claim by default.

### Route A. Decide U1 on the background alone (cheapest, first)

No wave is needed. The object is the conditional law of one relation bit
given the rest of a background order.

- A1 (analytic). For the background, compute the expected number, per event,
  of pairs whose single-relation flip is realisable by a background
  embedding, and the conditional flip probability of each. Cell-limited
  predicts O(1) pairs per event with probability of order one. Count-limited
  predicts about `sqrt(N)` pairs per event with probability about
  `N^(-1/2)` each.
- A2 (numerical, small N, CPU). Section 5, design D2.

Stop rule: if A1 or A2 shows the count-limited scaling, record "order
resolution of null margins is count-limited", drop Routes C and D below, and
put all effort on Route B. If cell-limited, drop Route B as stated (2 cannot
be necessary) and go to Routes C and D.

### Route B. Raise the necessary condition toward 2

Target: `TV -> 0` for every fixed `alpha < alpha_1` with `alpha_1 > 4/3`.
The cell argument cannot be reused with larger cells: for `h > 1/N` the
uncertain-pair graph percolates and the finite-component smoothing lemma is
unavailable. Two candidate replacements:

- B1. A second-moment (chi-square) bound on the likelihood ratio of the
  order laws after averaging over the wave phase and over the conformal
  symmetries of the background diamond, which preserve the background order
  exactly and move `u` by `N^(-1/2)`. The quarter-turn symmetry removes the
  linear term; the bound must control the quadratic term.
- B2. A sequential (event-by-event) Lindeberg replacement in which the
  conditional law of each new event's relation vector is compared under the
  two kernels given only the previous order, not the previous coordinates.

The prior-art report ([prior_art/report.md](prior_art/report.md), Section
3.2) identifies the nearest worked technique for B2: the vertex-exposure
chain rule with a posterior (cavity) analysis and a truncated chi-square,
which is the method that has produced sharp exponents for deterministic
hard-threshold relations in random geometric graphs. Every worked case there
has a product-measure null. Here both hypotheses are geometric, so the method
would have to compare two posteriors. No example of that was found.

Stop rules: (i) if the quadratic term in B1 has a contribution that does not
vanish for some `alpha < 2`, that term is a candidate statistic, Route B
stops and the term goes to Route C; (ii) if two sessions produce no bound
beyond 4/3, record the method barrier and stop.

### Route C. An explicit even statistic below 2

Only if U1 is cell-limited or Route B exposes a surviving quadratic term.
Candidate: the per-event light-cone quadrupole. The wave gives each event a
spin-2 quantity `B(u)` of size `lambda^2` that deforms its past and future
cones identically for all partners. The statistic correlates, for pairs of
flippable relations that share an event, their signs with the relative
transverse direction of the two partners. It needs no phase.

First step is a population challenge in the style of
`gram_fourier_challenge.md`: an analytic leading coefficient, then a frozen
quadrature. Stop rules: (i) the leading coefficient vanishes by symmetry, as
the linear cone residual did; (ii) the noise is not conditionally independent
across events and the exponent returns to 2; (iii) the statistic cannot be
written without a tightness oracle.

### Route D. Phase transport (U2)

Only if U1 is cell-limited. Measure or bound the posterior spread of `u`
differences. Stop rule: if it is `sqrt(d/N)` the heuristic value is 10/7 and
the work is the statistic of Route C with windows; if it is smaller, the
target moves toward 4/3.

### What would count as the discovery-sized outcome

Either a proof that `alpha_low = alpha_high = 2` (the order's resolution of a
plane wave is exactly the count limit `N lambda^2 ~ 1`), or an order-only
statistic proved to separate at some `alpha < 2` (the order stores null-margin
information below the count limit). A new upper or lower exponent strictly
inside the bracket, without the matching side, is progress and not that
outcome.

## 5. Design drafts. Not frozen. Not run.

These are drafts for discussion. Thresholds, grids and sample sizes below
are placeholders until a design is argued, frozen and stamped.

### D1. Finite-n check of the Gram-Fourier null noise law (optional)

Purpose: test the mechanism behind exponent 2, eq. (5) of
`gram_fourier_noise.md`, on the real four-dimensional background
calibration, at sizes a CPU can reach. It is not a power test of the
complete instrument. That test is out of reach: the frozen population
challenge has weight mass `4.4e-5` and normalised coefficient about `0.019`,
and its asymptotic regime starts near `1/lambda = 512`, so signal over noise
of order one needs `n` of order `1e9` times an unknown squared constant, with
an order-two U-statistic per probe.

- Estimand: `s(omega, n) = sqrt(E|S_n(omega)|^2)` under the background, fixed
  good anchor triple chosen by the generator (declared, not order-selected),
  probes drawn from the weight support (declared importance sampling), fixed
  smooth regularisation.
- Claims: slope of `log s` in `log omega` at fixed n is `-1` on
  `2 <= omega <= n^(1/3)`; slope in `log n` at fixed omega is `-1/2`.
- Gate type: both are accuracy claims, so the house rule applies. Pass only
  if the 95% two-sided interval of the fitted slope lies inside
  `target +/- t`. Draft `t = 0.15`.
- Sizing for 90% power: with R independent replicates per grid cell and
  approximately complex-Gaussian `S`, `Var(log s_hat) ~ 1/(4R)`. For G
  log-spaced grid points over a natural-log range D,
  `SE(slope) = (1/(2 sqrt(R))) / sqrt(D^2 G (G+1) / (12 (G-1)))`.
  The requirement is `t - |bias allowance| >= (1.960 + 1.282) SE`.
  Draft numbers: `G = 5`, `D = ln 16`, bias allowance `0.05`, giving
  `SE <= 0.0308` and `R >= 55`; use `R = 64`. Break-even at the same
  allowance is `R ~ 20`. That is 50% power and is not the price.
- Cost: `n` in `{1024, 2048, 4096}`; the ancestor-pair moment for all probes
  is one `n^3` matrix product, about 20 s at `n = 4096`; about one CPU hour
  per `n`.
- Retain every row. A failed slope is a result against the draft.

### D2. Background-only flip-predictability probe (Route A2)

- Sample background orders at `N` in `{64, 128, 256, 512}`.
- For each event list its links and critical pairs (the transitivity-safe
  flips). Record their number per event.
- With generator coordinates (declared: this is a diagnostic, not an
  observer) record for each such pair its null margin and whether flipping it
  is realisable by moving that one event inside the arrangement of the other
  cones.
- Gates are existence and ordering claims, not accuracy claims: is the count
  of realisable flips per event flat in N or growing, and does the margin
  scale of the realisable ones follow `1/N` or `N^(-1/2)`. If a slope is
  quoted as a recovered exponent it becomes an accuracy gate and must be
  sized as in D1 before freezing.
- Limit: single-event realisability is a necessary condition for a joint
  embedding and conditions on the other coordinates. It bounds U1 from the
  cell-limited side only. A joint-posterior sampler is the stronger and much
  more expensive version and is not drafted here.
