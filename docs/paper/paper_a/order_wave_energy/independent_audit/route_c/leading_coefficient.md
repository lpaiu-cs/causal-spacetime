# Route C, first step: the leading coefficient of the light-cone quadrupole

2026-11-26. Route C of [next_question.md](../next_question.md) Section 4,
started after the user approved Routes C and D ("C·D로 진행", relayed by
the session that started this one). AI-written. Reviewed by an AI referee
in a fresh context ([referee_report.md](referee_report.md)). This text was
corrected after that review; Section 9 says what changed. Generator
coordinates appear only in the deterministic checks of
[algebra_checks.py](algebra_checks.py). Nothing here is a measurement, and
no order was sampled.

**Summary.**

- The wave gives every event a spin-2 deformation of its light cones, of
  amplitude `2 sqrt(q) lambda^2 f(u)`. No symmetry of the problem removes
  its square, so the per-event quadrupole power has a nonzero leading
  coefficient (Sections 2 to 5).
- The statistic as written in next_question.md correlates the signs of an
  event's flippable relations with the relative transverse direction of
  its partners. **The response of that statistic that grows with N is
  absent at leading order** (Section 6). This is argued, not computed.
  - An event's own relations never see the wave: a unit-Jacobian shift of
    `v` absorbs them exactly.
  - The flippable relations add only relations among the event's partners
    near single null lines. There the deformation is again absorbed,
    except for overlap and boundary channels that do not grow with N.
- **Stop rule (i) therefore fires for the local form**, recorded
  conditionally (Section 7).
- The same argument removes the premise of the two "per event only" rows of
  next_question.md's heuristic table. One of them is one of the three
  accountings behind the "coincidence at 2".
- A joint per-event fit of position and shear is a new candidate, not a
  surviving form of Route C. It has no likelihood yet. Which resolution it
  would have is a question about U2, the transport of position across the
  diamond.

## 1. Setting

Conventions of [wavegeom.py](../proof3b_wave_signal/wavegeom.py) and
`gram_fourier_challenge.md`:

- `a(u) = cos(k u)`, `k = sqrt(q)`, `f = a^-2`, `m(u) = int f`, and
  `d = Delta m(u1, u2)`.
- The transverse vector is the complex number `zeta = y1 + i y2`.
- A trace-free symmetric matrix `[[b1, b2], [b2, -b1]]` is stored as
  `b1 + i b2`.

For an ordered pair, earlier event 1 and later event 2, the exact wave
margin is

```
mu_lam = Delta v - (1/2) Delta y^T Mt^-1 Delta y,   Mt = d I - B(u1, u2),
B(u1, u2) <-> beta(u1, u2) = fac e^{i phi} int_{u1}^{u2} f e^{i u / lam} du,
fac = 2 lam k / sqrt(1 + 4 q lam^2),
```

and the background margin is `mu_0 = Delta v - |Delta zeta|^2 / (2 d)`. The
background sampling density is proportional to `a(u)^2 du dv d^2y` on the
diamond. It depends on `u` alone, apart from the diamond boundary.

## 2. The first-order margin shift splits into per-event terms

Put `beta_E = beta(-1/2, u_E)`, so that `beta(u1, u2) = beta_2 - beta_1`.
Because `B` is trace-free, `Mt^-1 = (d I + B) / (d^2 - |beta|^2)`. Write
`omega = Delta zeta / d`: for a null pair this is the transverse velocity of
the null line through the two events. Define, for every event E,

```
phi_E(omega) = (1/2) Re( beta_E conj(omega)^2 ).
```

Then

```
mu_lam - mu_0 = phi_1(omega) - phi_2(omega) - (1/2) |omega|^2 |beta|^2 / d + ...   (1)
```

For a fixed event P and a partner X this reads
`delta(P, X) = sigma (phi_P(omega) - phi_X(omega))`, where `sigma = +1` for
a future partner and `-1` for a past partner.

The remainder is `O(lambda^4)` at fixed `d`, but only `O(lambda^3 |omega|^2)`
uniformly. At `d ~ lambda` the referee found about `3.4 lambda^3 |omega|^2`.

Check 1, on 4000 fixed pairs:

- The error of the first-order terms divided by `lambda^4` levels off.
  Its maximum is 187, 326, 426, 383 and 400 at `r = 1/lambda` = 64 to 1024,
  and its 99th percentile is 97 to 142.
- The shift itself is about `10 lambda^2`.

## 3. The per-event shear, and its gauge

Integration by parts gives

```
beta_E = b(u_E) + C + O(lambda^3),
b(u) = -2 i k lambda^2 f(u) e^{i (u / lambda + phi)} / sqrt(1 + 4 q lambda^2),
```

where `C` is a constant from the lower endpoint `u = -1/2`. Check 2:
differences of `beta_E` match differences of `b` to relative error
`0.72 lambda`, at every `r` from 64 to 1024. The `0.72` is the term
`lambda f'/f` at the ends of the range (`2 tan 0.35`).

**The gauge is larger than `C`.** The referee pointed out, and this record
re-derived, that every shear field quadratic in `m`,

```
s(m) = C0 + C1 m + C2 m^2,
```

(six real parameters) changes the first-order relations exactly as a
coordinate change of unit Jacobian to first order does. The coordinate
change is `zeta -> zeta + (a0 + a1 m) conj(zeta)` with
`v -> v + Re(a1 conj(zeta)^2) / 2`, and `C1 = 2 a0`, `C2 = a1`. The
referee's check A agrees to second order. These changes move the diamond
boundary, so relative to the tips they are not gauges (Section 6).

Every event therefore carries a spin-2 quantity `b(u_E)`:

- **Amplitude:** `2 k lambda^2 f(u_E)`, up to `(1 + 4 q lambda^2)^-1/2`.
- **Phase:** it turns once per wavelength in `u`.
- **Identifiability:** only up to a quadratic trend in `m`.

Along each null direction `omega`, (1) moves both cones of P by the same
amount: P behaves as if it sat at `v_P - phi_P(omega)`.

## 4. The event's own coordinates cannot absorb it

Moving P by `(dm, dv, dzeta)` changes the margin with every partner by

```
sigma ( -dv + Re(conj(omega) dzeta) - |omega|^2 dm / 2 ).      (2)
```

These are the angular modes 0, ±1 and the isotropic `|omega|^2`, all
independent of the distance `d`. The shear in (1) is the pair of modes
`±2`. Together the six span all quadratics in `(omega1, omega2)`: four
position parameters and two shear parameters per event. The shear cannot be
traded for a move of the event itself.

Against any weight that is smooth along P's generators, the partners'
shears and any partner coordinate errors that oscillate in `u_X` average to
`O(lambda^3)` (plus a slow trend, which is gauge). Check 3 fits P's exact
wave cone with the four parameters of (2) plus the two shear parameters.
It uses quadrature along both cones and **exact (oracle) partner
coordinates**. The table shows the first of three events; the other two
agree to the digits shown in the real part:

| r | fitted shear / `b(u_P)` | worst relative error, three events | largest position adjustment / `lambda^2` |
|---|---|---|---|
| 64 | 1.0230 - 0.0012i | 2.3e-2 | 2.7e-4 |
| 128 | 1.0026 - 0.0007i | 3.1e-3 | 6.7e-5 |
| 256 | 1.0000 - 0.0004i | 9.4e-4 | 1.7e-5 |
| 512 | 1.0000 - 0.0002i | 4.7e-4 | 4.2e-6 |
| 1024 | 1.0000 - 0.0001i | 2.4e-4 | 1.0e-6 |

## 5. The leading coefficient, and why no symmetry removes it

High-pass the fitted shears by subtracting a local quadratic trend in `m`
over a window much longer than `lambda`. A plain mean would leave a residual
`C2 Var(m)`. The per-event power of what is left has the population value

```
|b(u_P)|^2 = 4 q lambda^4 f(u_P)^2 / (1 + 4 q lambda^2).          (3)
```

The symmetries of the problem act on `b` as follows. None of them changes
`|b|^2`:

- the quarter-turn `(y1, y2) -> (-y2, y1)`, which removed the linear cone
  residual (`cone_constraint_challenge.md`), sends `b -> -b`;
- the reflection `y2 -> -y2` sends `b -> conj(b)`;
- a rotation by `a` sends `b -> e^{2ia} b`.

The referee found no exact cancellation at order `lambda^4`.

**This is the coefficient with the other events held at their true
positions.** In a joint fit two further things happen:

- shear trends and global position moves share an exact six-dimensional
  null direction;
- the other events' positions absorb part of the shear.

The frozen cone-constraint campaign (`cone_constraint_challenge.md`,
2026-09-24) already measured the second effect at population level. Its
least-squares residual is exactly the part of the first-order forcing
`sigma (phi_P - phi_X)` on near-null pairs that refitting **all** events'
coordinates cannot remove. On that campaign's coordinate-selected edges,
the residual was:

| N | residual / forcing RMS | residual RMS / `lambda^2` | mean degree |
|---|---|---|---|
| 128 | 0.45 to 0.50 | 0.88 to 1.03 | 11 to 14 |
| 256 | 0.57 to 0.63 | 0.96 to 1.08 | 16 to 17 |
| 512 | 0.61 to 0.67 | 0.95 to 1.04 | 19 to 22 |

The fully joint coefficient is therefore not 1. It is a substantial
fraction, and it rose with N over that range.

## 6. The local form has no response that grows with N

next_question.md describes the statistic as correlating, "for pairs of
flippable relations that share an event, their signs with the relative
transverse direction of the two partners". The flippable relations are the
transitivity-safe flips: links and critical pairs. The argument below is
first order in the shear and is not a bound.

**(a) An event's own relations never see the wave.** Fix P. The map
`v_X -> v_X - [F_lam - F_0](P; u_X, y_X)` moves the other events. Here
`F` is the `v`-offset of P's cone at `(u_X, y_X)`, so that
`mu = Delta v - F`. The map has unit Jacobian, and the density does not
depend on `v`. So P's whole relation vector, together with everyone
else's `(u, y)`, has the background law **to all orders**, apart from the
diamond boundary.

**(b) Flippable relations add partner-to-partner relations, and only near
single null lines.** Take a partner X near P's future cone along the null
line with transverse velocity `w`, at margin `-eps` and distance `L`. Take
Z on P's past cone at distance `s` with velocity `u`, at `t` below that
cone. The background chart gives exactly (check 4):

```
mu_0(Z, X) = t - eps + L s |w - u|^2 / (2 (L + s)).                (4)
```

- For a critical pair (P, X), adding `P < X`, the sets `D(P) \ D(X)` and
  `U(X) \ U(P)` are slivers along the same null line. Their angular width
  is `O(sqrt(eps / min(s, L)))`. A link's interval is a sliver of the same
  kind.
- The margins of links and critical pairs are typically `eps ~ N^-1/2`,
  not `1/N`, because a sliver's volume goes like `eps^2`. So the widths are
  about `N^-1/4`. The smallest margins, the walls, come from the largest
  distances.

What remains of P's own shear in these slivers:

- **Isolated sliver.** Within it `omega ≈ w`, and the shear acts as the
  shift `v -> v - G_w(u)`, with `G_w(u) = (1/2) Re(beta(-1/2, u) conj(w)^2)`.
  The part of the sliver's response that depends on `b(u_P)` integrates to
  zero at first order. (Referee's check B gives `1e-11`, against `1.1` for
  a shift that is not absorbed.)
- **Overlapping slivers.** With about `2.6 sqrt(N)` flippable partners per
  event (Route A) and widths `N^-1/4`, overlaps are the rule, not the
  exception. On an overlap the `b(u_P)` part survives, smaller than the
  unabsorbed tilt by about `0.56 r^2`, where `r^2 ~ eps / L` (check B).
  That is a relative `O(lambda^2)` per pair, with no gain in N.
- **Boundary channel.** The volume of `D(P) \ D(X)` is dominated by its
  far end, so every critical pair of every bulk event reaches the diamond
  boundary.
  - With P held fixed, the past exit condition tightens by `phi_P(w)` and
    the future one loosens by `phi_P(w)`, plus terms from the tips' shears.
  - This is a per-event mode-2 response carrying the phase of `b(u_P)`,
    plus a part common to all events set by the tips.
  - Its size: `d log V / d delta` of 0.1 to 0.3 in the referee's check C,
    against `2 / eps` (200 there) for an unabsorbed shift.
  - Links are not affected.
- **Second order.** A `cos 2 (theta_1 - theta_2)` correlation needs two
  first-order dependences on `b(u_P)` multiplied together. The
  second-order margin terms depend on `omega` only through `|omega|^2`. So
  the quadrupole correlation inherits the suppressions above.

So the response that the heuristic counted is absent at leading order:
each flippable item tilted by `delta / w`, growing with N. Responses with
no gain in N survive. The referee's scaling opinion is that the boundary
channel would need `alpha > 4` per event, and `alpha > 8/3` for its
tip-coherent part, so it does not threaten 2. These are scalings, not
bounds.

**(c) What the wave does change.** A single shift of `v` that depends on
`(u, y)` absorbs P's whole star, as in (a). It cannot also keep every
relation between partners that lie in different directions from P. That
residue of the partners' mutual relations is the Weyl curvature. The
argument does not prove too much: the Gram-Fourier statistic, which
separates above exponent 2, uses interval counts that span many events'
cones.

## 7. What this means for Route C

1. **Stop rule (i) fires for the local form, conditionally.** Record:

   > For the sign correlation of an event's flippable relations, the
   > coefficient that grows with N is zero at leading order. An event's
   > own relations absorb the wave exactly. The partner relations that
   > make a relation flippable absorb it, except in overlap and boundary
   > channels that have no gain in N.

   This was argued, not computed. The referee sketched a preregistered
   finite-N test ([referee_report.md](referee_report.md)). By
   next_question.md's rule the stopped route is a result. It is not
   converted into a smaller claim.
2. **Consequence for next_question.md Section 3.** Both "per event only"
   rows of the heuristic table assume that an event's near-null partners
   each carry the tilt `delta / w`: the cell-limited 8/5 and the
   count-limited 2. At leading order they do not, so both rows lose their
   premise. The count-limited per-event row is one of the three
   accountings behind the "coincidence at 2". That leg is gone. The other
   two legs, global phase lock and the windowed filter, are not touched by
   this.
3. **The joint per-event fit is a new candidate, with its own open
   questions.**
   - **It has no likelihood yet.** Free per-event shears in the
     first-order rule break transitivity: the referee's check D gives
     `mu(Y, Z) = mu(Z, W) = +1e-9` but `mu(Y, W) = -7.5e-7`. So "the
     posterior of an event's position and shear" is not defined. The
     defined objects are the wave likelihood itself, and how much of the
     wave forcing the background posterior absorbs (the cone-constraint
     residual of Section 5).
   - **Its requirement is an aggregate one**, `N (lambda^2 / w)^4 -> oo`
     with conditional independence across events, where `w` is the
     per-event shear resolution. The first version of this note asked for
     partner positions finer than `lambda^2`. That is impossible inside
     the bracket, since `1/N < lambda^2` only when `alpha > 2`.
   - **Its resolution is a U2 question.** An event's walls sit at
     macroscopic distance along its cones, so pinning them against each
     other is transport over distance about 1. The two cases, by scaling
     only:
     - if transport is count-like, the joint fit behaves like count
       resolution (the referee's scaling: `N^3 lambda^8`, `alpha > 8/3`);
     - if transport is cell-like, the phase is effectively known, the
       relevant row of the table is 4/3, and 8/5 would be only this
       statistic's own value.
   - **Opinion (the referee's):** six unknowns per event against about 9.5
     effective walls (Route A) may leave the cell much wider in the shear
     directions than in position.
4. **So the open question for Routes C and D is the same: U2**, the
   posterior spread of relative position across the diamond. Route D's
   stop rule frames it: `sqrt(d / N)` gives the heuristic 10/7; anything
   smaller moves the target toward 4/3.

## 8. Limits

- AI-level work: the derivation and its review are both by AI, with no
  human review.
- Section 6 is a first-order argument. Its residual channels are
  estimated by scalings and by the referee's deterministic sliver checks,
  not bounded, and nothing about their exponents is proved.
- Sections 2 to 5 use the plane-wave action of `wavegeom.py`, the
  referee-side geometry of the earlier audit. The generator code was not
  re-audited.
- Check 3 uses oracle partner coordinates and a fixed radial weight.
- The cone-constraint fractions of Section 5 come from a frozen campaign
  with its own coordinate-selected edges. They are not a statistic of the
  order.
- Nothing here bears on the bracket `4/3 <= alpha_low <= alpha_high <= 2`.

## 9. Review and what was corrected

[referee_report.md](referee_report.md): an AI referee in a fresh context,
given this note, the checks, the geometry and next_question.md. It was
told not to run any Monte Carlo of Section 6. The first version of this
note (in the history) was corrected as follows:

- **Section 3:** the gauge is the six-parameter quadratic family, not the
  constant `C`. The high-pass in Section 5 changed from a local mean to a
  quadratic detrend.
- **Section 6, rewritten.**
  - The first version argued that slivers in different directions are
    disjoint. They are not, at the true margin scale `N^-1/2`.
  - It said that no surviving channel is a per-event quadrupole. The
    boundary channel is one, without gain in N.
  - It gave the wrong reason why a joint statistic escapes.
  - It wrote "cell scale" where the margin scale is `N^-1/2`.
  - The conclusion about the response that grows with N stands, on the
    stronger argument (a) and the sliver analysis (b).
- **Section 7, rewritten.**
  - The first version said "the joint fit survives (i)" and asked for
    partner positions finer than `lambda^2`, which is impossible inside
    the bracket.
  - It attached the 8/5 row to the joint fit without the U2 assumption
    behind it.
  - It treated a per-event posterior of shear as defined.
  - It did not record the consequence for the "coincidence at 2".
- **Section 5:** added that the coefficient 1 holds only with the other
  events held fixed, together with the cone-constraint fractions.
- Wording of check 3, and the note that the table shows one of three
  events.

## 10. Files

- [algebra_checks.py](algebra_checks.py): checks 1 to 4 (CPU, about two
  minutes). Output as quoted above.
- [referee_report.md](referee_report.md): the review, with the referee's
  checks A to E summarised.
