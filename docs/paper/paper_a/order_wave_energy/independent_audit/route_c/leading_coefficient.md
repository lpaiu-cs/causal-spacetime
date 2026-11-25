# Route C, first step: the leading coefficient of the light-cone quadrupole

2026-11-26. Route C of [next_question.md](../next_question.md) Section 4,
started after the user approved Routes C and D ("C·D로 진행", relayed by
the session that started this one). AI-written. **Not yet reviewed**; an
independent review in a fresh context is recorded in Section 8 when it
exists. Generator coordinates appear only in the deterministic checks of
[algebra_checks.py](algebra_checks.py); nothing here is a measurement and no
order was sampled.

**In one paragraph.** The wave gives every event a spin-2 deformation of
its light cones. Its amplitude is `2 sqrt(q) lambda^2 f(u)`, and no symmetry
of the problem removes its square. So the per-event quadrupole power has a
nonzero leading coefficient (Sections 2 to 5). But the statistic as
written in next_question.md, correlating the signs of an event's flippable
relations with the partners' relative transverse direction, has **no
leading response in the bulk**. Along any one null direction the
deformation acts as a translation in `v`, and a measure-preserving map of
the sampling law absorbs it (Section 6). Stop rule (i) therefore fires for
that form. It does not fire for a joint fit of each event's position and
shear given the positions of the others. That fit needs the cell-resolution
posterior that Route D studies, so Routes C and D now meet at one object
(Section 7).

## 1. Setting

Conventions of [wavegeom.py](../proof3b_wave_signal/wavegeom.py) and
`gram_fourier_challenge.md`:

- `a(u) = cos(k u)`, `k = sqrt(q)`, `f = a^-2`, `m(u) = int f`, and
  `d = Delta m(u1, u2)`.
- The transverse vector is the complex number `zeta = y1 + i y2`.
- A trace-free symmetric matrix `[[b1, b2], [b2, -b1]]` is stored as
  `b1 + i b2`.

For an ordered pair, earlier event 1 and later event 2, the exact wave margin is

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
mu_lam - mu_0 = phi_1(omega) - phi_2(omega) + O(|beta|^2 |omega|^2 / d).   (1)
```

For a fixed event P and a partner X this reads
`delta(P, X) = sigma (phi_P(omega) - phi_X(omega))`, where `sigma = +1` for
a future partner and `-1` for a past partner. Check 1:

- On 4000 fixed pairs the error of (1) divided by `lambda^4` levels off.
  Its maximum is 187, 326, 426, 383 and 400 at `r = 1/lambda` = 64 to 1024,
  and its 99th percentile is 97 to 142.
- The shift itself is about `10 lambda^2`.

## 3. The per-event shear

Integration by parts gives

```
beta_E = b(u_E) + C + O(lambda^3),
b(u) = -2 i k lambda^2 f(u) e^{i (u / lambda + phi)} / sqrt(1 + 4 q lambda^2),
```

where `C` is a constant from the lower endpoint `u = -1/2`. `C` cancels from
every pair in (1), so it is a gauge: a common shear added to every event
changes no relation. Check 2: differences of `beta_E` match differences of
`b` to relative error `0.72 lambda`, at every `r` from 64 to 1024.

Every event therefore carries a spin-2 quantity `b(u_E)`:

- **Amplitude:** `2 k lambda^2 f(u_E)`, up to `(1 + 4 q lambda^2)^-1/2`.
- **Phase:** it turns once per wavelength in `u`.

Along each null direction `omega`, (1) moves both cones of P by the same
amount: P behaves as if it sat at `v_P - phi_P(omega)` (minus the gauge).
This is the identical deformation of past and future cones that
next_question.md describes.

## 4. The event's own coordinates cannot absorb it

Moving P by `(dm, dv, dzeta)` changes the margin with every partner by

```
sigma ( -dv + Re(conj(omega) dzeta) - |omega|^2 dm / 2 ).      (2)
```

These are the angular modes 0, ±1 and the isotropic `|omega|^2`, all
independent of the distance `d`. The shear in (1) is the pair of modes
`±2`. So each event has a six-dimensional perturbation space: four
position directions and two shear directions. The shear cannot be traded
for a move of the event itself.

What the partners contribute to P's cone:

- **Partners' shears.** The partner's `phi_X` depends on `u_X`, which varies
  along every generator of P's cone. Against any weight that is smooth along
  the generators, its oscillating part averages to `O(lambda^3)`. Its
  constant part is the gauge `C`, which the fit returns as an offset.
- **Partners' coordinate errors.** If the partners' coordinates carry errors
  that oscillate in `u_X`, those errors average out in the same way.

Check 3 fits P's exact wave cone with the six parameters of (2) plus the
shear, by quadrature along both cones with exact partner coordinates. At
three events:

| r | fitted shear / `b(u_P)` | worst relative error | largest position adjustment / `lambda^2` |
|---|---|---|---|
| 64 | 1.0230 - 0.0012i | 2.3e-2 | 2.7e-4 |
| 128 | 1.0026 - 0.0007i | 3.1e-3 | 6.7e-5 |
| 256 | 1.0000 - 0.0004i | 9.4e-4 | 1.7e-5 |
| 512 | 1.0000 - 0.0002i | 4.7e-4 | 4.2e-6 |
| 1024 | 1.0000 - 0.0001i | 2.4e-4 | 1.0e-6 |

From `r = 256` on, the error is `0.24 lambda`.

## 5. The leading coefficient, and why no symmetry removes it

The gauge-invariant per-event power is `|s_P - <s>_loc|^2`. Here `<s>_loc`
is the mean of the fitted shears over a `u` window much longer than
`lambda`, which removes `C`. Its population value is

```
|b(u_P)|^2 = 4 q lambda^4 f(u_P)^2 / (1 + 4 q lambda^2).          (3)
```

The symmetries of the problem act on `b` as follows. None of them changes
`|b|^2`:

- the quarter-turn `(y1, y2) -> (-y2, y1)`, which removed the linear cone
  residual (`cone_constraint_challenge.md`), sends `b -> -b`;
- the reflection `y2 -> -y2` sends `b -> conj(b)`;
- a rotation by `a` sends `b -> e^{2ia} b`;
- a shift of `u` by one wavelength sends `b -> b` up to the slowly varying
  `f`.

So for the per-event quadrupole power of a **joint fit**, stop rule (i)
does not fire. The normalised coefficient of a weighted sum is 1:

```
sum_P w(P) |s_P - <s>_loc|^2  /  ( 4 q lambda^4 sum_P w(P) f(u_P)^2 )  ->  1.
```

## 6. The local form of the statistic has no leading response

next_question.md describes the statistic as correlating, "for pairs of
flippable relations that share an event, their signs with the relative
transverse direction of the two partners". The flippable relations are the
transitivity-safe flips: links and critical pairs. This section argues that
at first order in the shear, the joint law of an event's flippable
relations is unchanged in the bulk. So no statistic of that local form has
a leading response.

**Sector locality.** Take a partner X near P's future cone along the
null line with transverse velocity `w`, at margin `-eps` and distance `L`.
Take Z on P's past cone at distance `s` with transverse velocity `u`, at
`t` below that cone. Then the identity of check 4 holds, exactly in the
background chart:

```
mu_0(Z, X) = t - eps + L s |w - u|^2 / (2 (L + s)).                (4)
```

Consequences:

- For a critical pair (P, X), adding `P < X`, the set `D(P) \ D(X)` is
  therefore a sliver along P's past generator with the same velocity `w`.
  Its angular width is `O(sqrt(eps / min(s, L)))`. The same holds for
  `U(X) \ U(P)` beyond X.
- For a link `P < X`, the interval is a sliver of the same kind around the
  segment.
- So whether P has a link or a critical pair with a partner in direction
  `w` depends only on events near the null line through P with velocity
  `w`. The exception is the vertex region `s, L <~ eps`, where all
  directions meet.

**Along one line the shear is a measure-preserving map.**

- For two events Y earlier and Z later near that line, `omega = w + O(width)`.
  So by (1), `mu_lam(Y, Z) = mu_0(Y, Z) + G_w(u_Y) - G_w(u_Z)`, with
  `G_w(u) = (1/2) Re(beta(-1/2, u) conj(w)^2)`, up to corrections
  proportional to the width.
- That is the background margin of the configuration mapped by
  `T_w : (u, v, y) -> (u, v - G_w(u), y)`. This map has unit Jacobian. It
  moves `v` by a function of `u` alone, and the background density depends
  on `u` alone.
- So away from the diamond boundary, `T_w` maps the background sampling
  law to itself.
- P is moved along with the rest. By stationarity in `v`, the law of the
  configuration relative to P is unchanged.
- Hence, given `u_P`, the law of P's flippable relations in direction `w`
  is the background law.

**Different directions use disjoint slivers**, so they are independent.
The joint law of P's flippable relations across all directions, with
their partners' directions estimated from counts, is the background joint
law. In particular the common shear `b(u_P)` induces no correlation between
the signs in two directions. The `|b(u_P)|^2 cos 2(theta_1 - theta_2)`
term that the heuristic expected is absent at leading order.

**What survives, by scaling only (not bounded):**

- the variation of the shear across a sliver, relative size
  `O(sqrt(eps / L))`, which is a negative power of N at the cell scale;
- the vertex region;
- the diamond boundary, which `T_w` moves by `O(lambda^2)`;
- partner directions estimated from order counts, whose `O(lambda^2)`
  wave-induced shifts enter only at `O(lambda^4)` with weight `O(1)` per pair
  and no gain in N;
- the `O(lambda^4 / d)` remainder of (1).

None of these is the per-event quadrupole, and none was shown to beat the
exponent 2.

**Why the joint fit escapes.** The map `T_w` depends on `w`. A single map
for all directions would be needed to absorb the shear everywhere, and
none exists: the shear is a mode ±2 in `w`, while a single `v` shift is
mode 0. So the shear shows only where the order compares directions. It
shows as an inconsistency between P's position read through its relations
in different directions, relative to partners whose own positions are
pinned against each other more finely than `lambda^2`.

## 7. What this means for Route C

1. **Stop rule (i) fires for the local form.** Record: "the sign
   correlation of an event's flippable relations has no leading response to
   the wave in the bulk; along each null direction the per-event cone
   deformation is a translation that a measure-preserving map of the
   sampling law absorbs." By next_question.md's rule this is a result. It
   is not converted into a smaller claim.
2. **The joint per-event fit survives (i).** Its population coefficient is
   1 with exact partner coordinates (check 3). Two things are open:
   - With partner coordinates estimated from the order, a smooth
     `O(lambda^2)` distortion of the reconstructed frame can leak into the
     fitted shear. Examples are the anchor-phase terms of the Gram
     reconstruction, or the fitted coordinate adjustment of a posterior
     embedding. Such a distortion changes P's cone by terms that grow with
     `d`, which project onto the shear column. It is smooth in `u_P`,
     whereas `b(u_P)` turns once per wavelength, so a local mean
     subtraction should separate them. This has to be computed, not assumed.
   - The fit needs the other events' positions to better than the shear
     itself, about `lambda^2` in margin units. Count-based coordinates
     (`N^-1/2`) cannot supply that whenever `alpha < 4`. Only the
     cell-resolution joint posterior can, which is the object of
     Routes A and D.
3. **Stop rule (ii) is untested.** The exponent of the joint statistic is
   set by how finely the order fixes a per-event shear given the rest:
   - cell-like (`~1/N` in margin units): the heuristic row "cell-limited,
     per event only" of next_question.md (8/5);
   - count-like: back to 2 or above.

   Nothing here estimates which.
4. **Stop rule (iii).** A posterior functional is order-only. Computing it
   needs an embedding found from the order. Starting a chain from
   generator coordinates, as Route D did, is legitimate for the statistical
   question only if the chain mixes; Route D's non-ergodicity caveat
   carries over unchanged.

So Routes C and D now meet at one object: **the posterior of an event's
position and shear given the order**.

## 8. Limits

- AI-level derivation, not yet reviewed.
- Section 6 is a first-order argument. The sector independence and the
  listed corrections are scalings, not bounds, and there is no theorem
  that the corrections stay below the exponent 2.
- Sections 2 to 5 use the plane-wave action of `wavegeom.py`. That module
  is the referee-side geometry of the earlier audit, and the generator code
  was not re-audited here.
- Check 3 uses generator coordinates for the partners and a fixed radial
  weight. It verifies (1) to (3) and nothing about finite N.
- Nothing here bears on the bracket `4/3 <= alpha_low <= alpha_high <= 2`.

## 9. Files

- [algebra_checks.py](algebra_checks.py): checks 1 to 4 (CPU, about two
  minutes). Output as quoted above, at the commit that adds this note.
