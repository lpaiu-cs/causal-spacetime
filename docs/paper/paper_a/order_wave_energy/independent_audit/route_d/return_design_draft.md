# Route D, second measurement: does the chain forget a far start? DRAFT

2026-11-25. **Draft. Not frozen. Not implemented. Not run.** It becomes a
design only when the user approves it; freezing is then a separate commit.
Generator coordinates throughout: a diagnostic of the background order law,
not an observer and not a wave test. AI-written; not reviewed. It replaces
the withdrawn joint-slack draft ([RECORD.md](RECORD.md) Section 2).

## 1. The question

The first measurement ([results.md](results.md)) found that the residual
of central events, after the seven-parameter family is fitted out, stays
near `12.4 / N` along the chain. The doubt left is whether that is the
size of the whole set of typical configurations with the given order, or
only of a part the single-event chain stays in.

Test: push the configuration, with its order kept exactly, along a
collective direction that is not a family motion; start the chain there;
watch whether the push decays.

- If the conditional law is concentrated at the cell scale along that
  direction, a pushed start is atypical and the chain is pulled back: the
  push is **forgotten**.
- If the direction is free over a range larger than the push, or the
  pushed start lies in a part the chain cannot leave, nothing pulls it
  back: the push is **remembered**.

This is why it does not suffer from the fault of the withdrawn draft. The
far configuration is only a starting point; what is measured is where the
chain goes, which is a property of typical configurations, not of an
extreme one.

## 2. Quantities

`x0`: sampled configuration, order C, central set and group fit as in
[design.md](design.md). `s1`: one-sweep displacement of central events
along the reference chain. For any configuration y with order C, `r(y)` is
the residual field on the central events after the group element fitted
from `x0` to y is removed.

**Push fields.** With X the chart position relative to the midpoint of the
tips: *shear*, `(0, X_y1, X_Z, 0)`; *quadratic*, `(X_y1 X_y2, 0, 0, 0)`.
Neither is a conformal Killing field. Each is evaluated at the central
events of `x0`, made orthogonal there to the seven generator fields, and
normalised to unit root mean square: `F`.

**Memory.** `P(y) = mean over central events of F_k . r_k(y)`, and for a
chain started at a pushed configuration `x*`, `M(t) = P(y_t) / P(x*)`.
`M = 1` at the start; a chain that has forgotten the push has `M` near 0.

**Distance.** `R(t)`: root mean square of `r(y_t)`, with the correction
for seven fitted parameters, divided by `s1`.

## 3. Making the pushed start

A path of linear programmes in which all events move, started at `x0`.
Rows are the links and critical pairs of C, including the tip relations.
In each round, at the current configuration:

- maximise `sum over central events of F_k . xi_k`;
- every row keeps at least 0.3 of its current margin (linearised);
- the cumulative displacement from `x0` satisfies the seven family pins
  (`sum over central events of G_a(x0,k) . (x - x0)_k = 0`);
- `|xi| <= rho` per coordinate, `rho` starting at `2 / N`.

The step is scaled by the largest of `1, 1/2, 1/4, 1/8` for which the
order matrix equals C exactly and every row margin exceeds `1e-9`. If none
works `rho` is halved (at most four times in a row, then the path stops);
after an accepted step `rho` grows by 1.5 up to `8 / N`. The path stops
when the distance of the current configuration reaches 3 times the
reference plateau `12.4 / N`, after 60 rounds, or after three rounds that
each add less than 1%. Solver: HiGHS interior point. The sign of the push
is drawn at random per order and field.

Every accepted configuration has order exactly C, so the result is a
legitimate start whatever the quality of the linear model. If the path
stops short, the start is less far; it is still used and its distance is
reported.

## 4. Chains

Per order: one reference chain from `x0`, and one chain from each pushed
start. The Route A chain, unchanged in law. Length `N^2` sweeps. States
stored after `N 2^j` sweeps (`j = 0, 1, ...` up to `N^2`), after `N^2/4`,
`N^2/2`, `N^2`, and one sweep after the start of the reference chain for
`s1`. Exact order comparison at every stored state.

For the reference chain the same `P` is computed for both fields: it has
mean zero and gives the noise against which `M` is read.

## 5. Controls with answers fixed in advance

Both use an exact group element, so no linear programme is involved, and
both are scored on the fitted group parameter, not on the residual.

1. **A free direction must be remembered.** Rotate `x0` about the wave
   axis by the angle that moves central events by `30 s1` root mean
   square. The rotation preserves the order and the density exactly, so
   nothing pulls it back; the chain can only diffuse. Memory is the fitted
   rotation parameter over the applied one. Requirement at lag `N^2`: the
   95% interval lies above 0.5.
2. **A stiff direction must be forgotten.** Push `x0` along the stiffest
   eigenvector of the family second-derivative matrix by 8 times its
   predicted equilibrium spread. The order is preserved exactly; the
   density is not, and pulls it back. (In the first measurement the stiff
   parameters reached their spread by about `N^2 / 4`.) Memory is the
   fitted parameter along that eigenvector over the applied one.
   Requirement at lag `N^2`: the 95% interval lies below 0.3.

These show that the instrument returns "remembered" for a collective
motion that is free and "forgotten" for one that is held, at the sizes and
chain length of the measurement. If either fails, the main chains are not
read.

## 6. Campaign

CPU only; twelve single-thread workers at below-normal priority, recorded;
clean-tree preflight, commit stamp, no overwriting. Seeds
`SeedSequence([20261129, 4, N, r])`. Sizes N = 128 and 256.

| stage | chains | orders | estimated wall time |
|---|---|---|---|
| 0 | the two controls | 12 at N = 128, 6 at N = 256 | 20 min |
| 1 | reference, shear start, quadratic start | 12 at each size | 45 min |

Building a pushed start is estimated at under two minutes at N = 256.

## 7. Descriptions fixed in advance. No verdict on U2

Per size and field, means over orders with 95% t intervals: the starting
distance over the reference distance at lag `N^2`; `M(t)` at every stored
lag; `R(t)` of the pushed chain over `R(t)` of the reference chain.

- **Not far enough.** The interval of the starting distance over the
  reference distance does not lie above 2. No further description; the
  numbers are reported. (That a collective push cannot be made far is
  itself what a stiff cell would give, and is not evidence by this rule.)
- **Forgotten.** At lag `N^2` the interval of `M` lies below 0.2 and the
  interval of the distance ratio has an upper end of at most 1.10.
  *Establishes:* along this field, at this distance, the chain is pulled
  back to the neighbourhood of `x0`; the plateau of the first measurement
  is not an artefact of the chain being unable to leave a displaced
  configuration of this kind. *Does not establish:* the same for other
  directions, for larger displacements, or for configurations no path of
  this kind reaches.
- **Remembered.** At lag `N^2` the interval of `M` lies above 0.2.
  *Establishes:* the chain did not forget a push of this size in `N^2`
  sweeps. *Does not establish* which of three things caused it: a free
  direction, a part of the set the chain cannot leave, or a slow return.
  *Consequence:* the first measurement's plateau cannot be read as the
  size of the typical set along this direction.
- **Neither** otherwise.

The first lag at which the mean of `M` is below 0.2 is reported as the
forgetting time.

## 8. Limits

- Two push fields, one distance (about 3 times the plateau), two sizes
  with the larger at 256, central events, one background.
- A push is collective and smooth by construction. Freedom that is
  collective but not smooth is not probed.
- The pushed start is a corner of the cell (many margins are small after
  the path). A corner is atypical in ways other than the push; the
  reference chain and the control say how fast the chain leaves it.
- "Forgotten" within `N^2` sweeps says the direction is held. It does not
  say how tightly beyond the plateau value.
- Nothing here is about a wave, a statistic or the exponent.

## 9. What has been seen that bears on this

The Route A and Route D records and their checks, and the review of the
withdrawn draft: straight steps along family generators keep the order
only to a few `1/N` of central displacement, and a default linear-programme
solver is too slow at N of 512 and above. No pushed start has been built
and no chain has been started away from a sampled configuration.
