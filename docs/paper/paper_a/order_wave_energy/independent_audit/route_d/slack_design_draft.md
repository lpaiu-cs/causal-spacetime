# Route D, second measurement: joint slack without a chain. DRAFT

2026-11-25. **Withdrawn draft. Not to be frozen. Never implemented or
run.** An independent review
([slack_referee_report.md](slack_referee_report.md)) found that its
positive description would be returned for a cell that is typically rigid,
that its positive control cannot pass, and that its certification
under-reports slack that exists. The reasons are summarised in
[RECORD.md](RECORD.md). The text below is kept as it was reviewed.
Generator coordinates throughout: a diagnostic of the background order law,
not an observer and not a wave test. AI-written.

## 1. Why

The first measurement ([results.md](results.md)) found that along an
order-preserving chain the residual of central events, after the
seven-parameter family is fitted out, stays near `12 / N`. The largest
doubt left by its independent check is non-ergodicity: the chain moves one
event at a time, and freedom that needs many events to move together might
not be reached in `2 N^2` sweeps.

This measurement does not use a chain. At a sampled configuration it asks a
linear programme, in which **all events move at once**, how far a central
event can be displaced while every relation is kept and the family is held
fixed. It then checks the answer against the exact order.

What it can and cannot give:

- A displacement that passes the exact order comparison is a **certificate**:
  a configuration with the same order exists in which that event sits that
  far away after the family is removed. Certified slack is a rigorous lower
  bound on the joint slack.
- The linear programme's own optimum is a model value. It is not a bound in
  either direction (Section 3).
- Slack is a maximum, not a typical spread. Large certified slack would not
  show that positions are typically that uncertain; it would show that the
  chain's plateau cannot be read as rigidity. Small model slack would
  support rigidity and not prove it.

## 2. Quantities

Sampled configuration `x0` with order C; rows = links and critical pairs of
C including the relations with the two tips (the irredundant set of Route
A's first design: if all rows hold, the order is C). `m_c` is the oriented
null margin of row c and `g_c` its gradient.

**Family pin.** `G_a(x)` is the generator field of the a-th parameter of the
order-preserving group of [design.md](design.md) Section 2, by central
differences of the group map. A displacement `xi` of all events is pinned
when `sum over central events of G_a(x0,k) . xi_k = 0` for all seven a.
This is the first-order form of "the best-fitting group element between the
two configurations is the identity", with the same central set as before.

**Directions.** For a tested central event i: the wave coordinate
`(T + Z) / sqrt(2)` and one transverse coordinate `y1`, each in both signs.

**Three programmes,** each maximising the displacement of event i along the
direction, inside a box `|xi| <= 4 / sqrt(N)` per coordinate:

- *single*: only event i moves. Rows: its relations with every other event.
- *joint*: every event moves; all rows; seven pins.
- *joint, one pin released* (positive control, Section 4): as joint, with
  the dilation pin removed.

Lengths are multiplied by N inside the programmes so that solver
tolerances are small against the margins.

**Certification, the same for all three.** Rounds, at most six. In a round
the programme is linearised at the current configuration `x_r` and
maximises the displacement of event i from `x0`; the pins are applied to
the cumulative displacement from `x0`. The step is scaled by the largest
`theta` in `1, 1/2, ..., 2^-10` for which the order matrix of the stepped
configuration equals C exactly; `x_(r+1)` is that configuration. Rounds
stop when a round adds less than 5% or no `theta` works. For the final
configuration the group element is fitted on the central events as in the
first measurement (non-linear least squares), and the **certified
displacement** of event i is its residual along the direction.

**Reported per tested event and coordinate,** summing the two signs:

- `w_single`: certified width with only event i moving;
- `w_model`: the joint programme's optimum in the first round;
- `w_cert`: certified joint width;
- `J_model = w_model / w_single`, `J_cert = w_cert / w_single`;
- whether the box was reached, the rounds used, the largest and
  root-mean-square displacement of the other events in the certificate.

## 3. What the linear model is, and is not

A link row requires `m_c > 0`. The margin is concave in the displacement,
so its tangent over-estimates it: the linearised row admits displacements
the true cone does not. A critical row requires `m_c < 0`; there the
tangent is the stricter condition. So the model region is neither inside
nor outside the true one, and `w_model` is not a bound. It is accurate when
the displacements it uses are small against the separations of the pairs
involved, which holds for cell-scale moves and fails for box-scale ones.
That is why every reading below leans on `w_cert` where it can.

## 4. Controls

1. **The programme finds order-preserving joint motion when it exists.**
   With the dilation pin released, the family direction it frees is an
   exact first-order symmetry of every row, so the model optimum must reach
   the box, and the certified displacement must be large. Requirement, at
   every size: the box is reached in every released programme, and the
   median over tested events of `w_cert(released) / w_cert(joint)` has a
   95% interval lower end above 2. This is the positive control for
   count-scale joint motion that the chain instruments lacked in 3+1: the
   motion exists, is collective, and must be found and certified.
2. **Certificates are exact.** Every certified configuration is compared
   with C by the full order matrix, inside the diamond; a mismatch is an
   error.
3. **Single against the chain.** `w_single` should scale like `1 / N`. Its
   mean times N is reported per size against the chain's one-sweep
   displacement (`2.9 / N`) and the frozen-cell radius found by the Route A
   check (`9.3 / N`). Reported, not gated: the three are different
   functionals of the same cell.

There is no control with a known *small* joint slack. Section 6 keeps the
wording that follows from that.

## 5. Campaign

CPU only; twelve single-thread workers at below-normal priority, recorded;
clean-tree preflight, commit stamp, no overwriting. Seeds
`SeedSequence([20261128, 4, N, r])`.

| N | orders | tested central events per order | estimated wall time |
|---|---|---|---|
| 128, 256, 512 | 12 each | 16 | 15 min together |
| 1024 | 12 | 16 | 30 min |

Released-pin control: 4 of the tested events per order. Estimates assume 1
to 3 seconds per joint programme at N = 1024; the first minutes of the run
will show whether that holds, and a size that would take more than three
times its estimate is stopped and reported, not shortened.

## 6. Descriptions fixed in advance. No verdict on U2

Per size: means over orders of `J_model` and `J_cert` for each coordinate,
with 95% t intervals; the fraction of joint programmes that reached the
box; the widths in chart units times N and times `sqrt(N)`.

Across sizes, for each coordinate, the least-squares fit of the per-order
values to `J^2 = A + B N`, with intervals, for `J_cert` and for `J_model`.

- **Count-scale joint slack certified.** The interval of B for `J_cert`
  lies above zero. *Establishes:* configurations with the same order exist
  in which a central event is displaced, beyond the family, by an amount
  that grows relative to the single-event width. *Does not establish:* that
  such displacements are typical. *Consequence for the first measurement:*
  its plateau cannot be read as rigidity; the chain either does not reach
  these configurations or they carry little weight.
- **No count-scale joint slack found.** The intervals of B for both
  `J_cert` and `J_model` contain zero, and fewer than 5% of joint
  programmes reached the box. Reported with the smallest B the data could
  have told from zero. *Does not establish* rigidity: the model is not a
  bound, the search is local to the sampled configuration, and three sizes
  a factor 8 apart with one in the pre-asymptotic range make a modest
  lever.
- **Model slack not certified.** B for `J_model` lies above zero and B for
  `J_cert` does not. Reported as that: the linear model has room the exact
  order did not confirm.
- Anything else is reported as numbers without a description.

## 7. Limits

- Maximum, not typical (Section 1). Linear model, not a bound (Section 3).
- Local: the programme starts at the sampled configuration and its
  certificates stay connected to it through order-preserving steps.
  Disconnected sets of configurations with the same order are invisible.
- The pin is first order and on the central set. A motion that looks like
  a family element on the central events is removed by it whether or not
  it is one elsewhere; the rows involving the other events are what stops
  such a motion, and they are all in the programme.
- Central events, 16 per order, two coordinates, one background, one
  diamond, sizes to 1024.
- Nothing here is about a wave, a statistic or the exponent.

## 8. What has been seen that bears on this

The Route A and Route D records and their checks. In the pilot of Route
A's first design (N = 96) a walk built from programmes of this kind
certified many flips that single-event moves did not; at that size the cell
and count scales coincide. No joint programme with pinned family has been
run.
