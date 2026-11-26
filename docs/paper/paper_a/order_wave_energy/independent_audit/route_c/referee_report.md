# Referee report on the Route C leading-coefficient note

2026-11-26. **AI-level review; no human review.** An AI referee in a fresh
context, which had not seen this work, was given the first version of
[leading_coefficient.md](leading_coefficient.md) (in the history), its
checks, the plane-wave geometry, next_question.md and the cone-constraint
record. It was told to find what is wrong rather than to confirm. It was
also told not to run any Monte Carlo of the note's Section 6, so that a
later preregistered test of that claim stays clean. Its deterministic
checks A to E are [referee_checks.py](referee_checks.py), with output in
[referee_checks_out.txt](referee_checks_out.txt). The report below is the
referee's, lightly formatted. The note was corrected after it (its
Section 9).

It reproduced `algebra_checks.py`: every number the note quotes. The
note's Section 4 table shows the first of three events. No order was
sampled and no Monte Carlo of Section 6 was run.

## Verdicts

1. **Eq. (1): correct.**
   - Checked: `Mt^-1 = (dI+B)/(d^2-|beta|^2)`, the sign, the `1/2`,
     `omega = Delta zeta / d` (later minus earlier), and `+phi_P` for a
     future partner, `-phi_P` for a past one.
   - The remainder `-(1/2)|omega|^2|beta|^2/d` is `O(lambda^4)` at fixed d
     but only `O(lambda^3 |omega|^2)` uniformly. Check E gives
     `err/(lambda^3 |omega|^2) = 3.4` at `d ~ lambda`. This is harmless.
   - C cancels exactly from every event-to-event relation.
   - **The gauge statement is incomplete.**
     - (a) `s(m) = C0 + C1 m + C2 m^2` (six real parameters) is pure gauge
       at first order. The map is
       `zeta -> zeta + (a0 + a1 m) conj(zeta)`,
       `v -> v + Re(a1 conj(zeta)^2)/2`, with Jacobian `1 - |a|^2`. Check A
       gives residual / eps^2 constant at 9.245.
     - (b) C is not gauge relative to the diamond tips (see 5c).
2. **Section 3: correct.**
   - `b(u)`, the `(1+4q lambda^2)^-1/2` factor and the phase all hold.
   - Check 2's constant 0.72 is the `lambda f'/f` term (`2 tan 0.35`).
   - "As if at `v_P - phi_P(omega)`" holds for both cones.
3. **Section 4: correct.**
   - The six functions span the quadratics in `(omega1, omega2)`. Position
     moves give `1`, `omega`, `conj(omega)` and `|omega|^2`.
   - Check 3 is set up correctly: the signs, the partners placed on the
     exact wave cone via `inv_action`, and the right-hand side `-mb`.
   - Partner shears enter only through their weighted `u`-average, so the
     fit has no gauge problem.
   - Wording: "six parameters of (2) plus the shear" should read four plus
     two.
   - Its coefficient of 1 rests on oracle partner `v` coordinates (see 5a).
4. **Section 5: (3) is correct**, and no discrete symmetry removes
   `|b|^2`.
   - Because of 1(a), "`<s>_loc` = local mean" must become a local quadratic
     detrend in `m`. A symmetric mean leaves `C2 Var_W(m)`.
   - In a joint fit, shear trends and global position moves share an exact
     six-dimensional null direction.
   - (3) stands for the high-passed power.
5. **Section 6: the conclusion is right, but parts of the argument are
   wrong.**
   - **(a) A stronger proof exists.**
     - The exact wave margin is `Delta v - F_lam(u_P, y_P, u_X, y_X)`. The
       map `v_X -> v_X - [F_lam - F_0](P; u_X, y_X)` has unit Jacobian,
       and the density does not depend on `v`.
     - So P's whole relation vector, together with all others' `(u, y)`,
       has the background law to all orders. This holds apart from the
       boundary, and with fixed N under Palm conditioning.
     - Links and critical pairs add Z-X relations along one generator. On
       an isolated sliver the `b_P`-dependent residual integrates to zero
       at first order. Check B gives `1e-11`, against `1.1` if the shift is
       not absorbed.
   - **(b) Step (c), "disjoint slivers", is false.**
     - Links and critical pairs have margin `eps ~ N^-1/2`, since sliver
       volume is proportional to `eps^2`; not `1/N`.
     - The velocity width is therefore about `N^-1/4`. With about
       `2.6 sqrt(N)` partners (Route A), neighbouring slivers overlap
       generically.
     - On overlaps the `b_P` residual survives, suppressed by `r^2 ~ eps/L`
       relative to the unabsorbed tilt. Check B gives ratio / r^2 =
       0.56 to 0.57. That is relative `O(lambda^2)` per pair, with no gain
       in N.
     - The survivors list should say `N^-1/2`, not "cell scale".
   - **(c) The boundary is not a boundary layer.**
     - `vol(D(P)\D(X)) = pi eps^2 int_0^S s(L+s)/L ds` is dominated by its
       far end, so every critical pair of every bulk event reaches the
       boundary.
     - In the frame with P held fixed, the past exit condition tightens by
       `phi_P(w)` and the future one loosens by `phi_P(w)`, plus terms in
       the tips' shears.
     - This is a per-event mode-2 response with the phase of `b(u_P)`,
       plus a tip-coherent one. Check C gives `d log V / d delta` of
       +0.16, -0.11, +0.17 and -0.33, against 200 (`2/eps`) if unabsorbed.
       Links are insensitive.
     - So "None of these is the per-event quadrupole" is false. What is
       true: the N-enhanced coefficient vanishes.
     - Opinion, by the same accounting: this channel needs `alpha > 4` per
       event, or `alpha > 8/3` tip-coherent, so it does not threaten 2.
   - **(d) Second order holds**, apart from (b) and (c). A
     `cos 2 (theta1 - theta2)` term needs a product of two first-order
     `b_P` dependences, and the second-order margin terms are isotropic in
     `|omega|^2`.
   - **(e) "Why the joint fit escapes" gives the wrong reason.**
     - A `v`-shift depending on `(u, y)` is one measure-preserving map that
       absorbs every direction of P's star. Counting angular modes rules
       out only shifts that depend on `u` alone.
     - The real obstruction is the relations between partners in
       different directions (Weyl curvature).
     - The argument does not prove too much: the Gram-Fourier statistic
       uses interval counts that span many stars.
   - **(f)** `T_w` under fixed N, Palm conditioning on P, and stationarity
     in `v`: fine in the bulk. The vertex holds about `N eps^4` events, so
     it is negligible.
6. **Section 7 overclaims.**
   - **7.2** asks that positions be known "better than `lambda^2`", which
     "only the cell-resolution posterior can" supply.
     - Since `1/N < lambda^2` iff `alpha > 2`, nothing meets this anywhere
       in the interval of interest. It also contradicts 7.3's 8/5, where
       `N lambda^2 -> 0`.
     - The correct requirement is aggregate:
       `N (lambda^2 / w)^4 -> infinity`, with conditional independence.
   - **7.3 attaches the wrong U2 assumption to the 8/5 row.**
     - P's wall partners sit at macroscopic distance along its cones:
       shell measure grows with L, so the smallest margins come from the
       largest L.
     - Pinning them against each other at cell scale is transport over
       `d ~ 1`, which is U2.
     - Under "U2 none" the joint fit is count-like: `N^3 lambda^8`,
       `alpha > 8/3` (the referee's scaling).
     - Under cell-like transport the phase is known, so the problem's row
       is 4/3, and 8/5 is only this statistic's value.
     - Opinion: six unknowns against about 9.5 effective walls (a 4-D
       measurement) may leave the 6-D cell much wider in the shear
       directions.
   - **Unrecorded consequence:** Section 6 removes the premise of both
     "per event only" rows. That includes the count-limited row, one of
     the three legs of next_question.md's "coincidence at 2".
   - **Stop rule (i)** is defensible for the N-enhanced coefficient only.
     - Route C asked for an analytic coefficient and then a frozen
       quadrature. For the local form the zero was argued, not computed,
       so it should be recorded conditionally.
     - The joint fit is a new candidate with its own stop rules, not
       "Route C survives".
   - **7.4 / (iii):** the "posterior of position and shear" has no
     likelihood yet. Free per-event shears in the first-order rule break
     transitivity. Check D gives `mu(Y,Z) = mu(Z,W) = +1e-9` but
     `mu(Y,W) = -7.5e-7` at `|s| = 1e-3`.

## Ranked weak points

1. 7.2 and 7.3: an unmeetable requirement; the 8/5 row with the wrong U2
   assumption; the consequence for the coincidence at 2 left unrecorded.
2. Section 6: step (c) (disjointness), and the reason given for why the
   joint fit escapes.
3. The boundary channel is miscast; an unenhanced per-event quadrupole
   survives.
4. The gauge is incomplete: six real parameters, not two.
5. The posterior is undefined: there is no transitive model behind it.
6. Stop-rule wording.

## A preregistered finite-N test of Section 6 (design only, nothing run)

- **Estimand.** The first-order response with the phase taken from
  generator coordinates. This also bounds the quadratic statistic.
  - For central events P and each link or critical pair X, regress the
    indicator on `x = Re(b(u_P) conj(omega_PX)^2)`. Generator coordinates
    are declared as a diagnostic.
  - Stratify links against critical pairs, and isolated against
    overlapping partners.
  - Report `R = kappa_hat / kappa_H`, where `kappa_H` is the unabsorbed
    tilt from a deterministic quadrature.
- **Frozen predictions.** Section 6 predicts `R_links` inside ±0.25 at
  every N, and `R_crit` matching the boundary-channel quadrature, `O(eps)`.
  The heuristic predicts `R ~ 1`.
- **Grid.** N from 256 to 2048 with `lambda^2 sqrt(N)` held fixed, several
  phases, the exact wave relation and exact wave-diamond membership.
- **Positive control.** Shift only P's relations by `phi_P(omega)` against
  fixed partners and recompute; it must return R in `1 ± 0.25`.
- **Null control.** Background orders with a fictitious phase.
- **House rule.** Decide on the 95% two-sided interval of R. For 90% power,
  `SE(R) <= t/3.24` (0.077 at `t = 0.25`), with variances from
  background-only pilots.
- **Secondary.** The order-only `cos 2 Delta-theta` statistic, predicted
  null.

## Not checked

- The generator code.
- The `a(u)^2` density of the exact wave (taken from `wavegeom.py`).
- The exact wave-diamond geometry (check C is flat and first order).
- Orders beyond first in `b` on sliver volumes, and fixed-N corrections.
- Any bound on the survivors, and the size of the 6-D cell.
- Routes A and D were not audited.

The referee's exponents (4, 8/3, 8/3) are scalings, and opinion.
