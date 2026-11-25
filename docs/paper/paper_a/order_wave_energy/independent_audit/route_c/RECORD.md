# Route C: progress record, and the proposed next step for Routes C and D

2026-11-26. AI-run work; no human review. Routes C and D were approved by
the user ("C·D로 진행"), relayed by the session that started this one. No
design has been drafted in full, frozen or run under this approval.

## 1. First step of Route C: done, and stopped by its own rule

[leading_coefficient.md](leading_coefficient.md),
[referee_report.md](referee_report.md).

- **The per-event quadrupole exists.**
  - The wave gives every event a spin-2 deformation of its light cones, of
    amplitude `2 sqrt(q) lambda^2 f(u)`, with a phase that turns once per
    wavelength in `u`.
  - It is identifiable only up to a shear field quadratic in `m`, which is
    a six-parameter gauge.
  - No symmetry removes its square. With the other events held at their
    true positions, a fit of one event's exact cone returns it with
    coefficient 1.
  - With all events refitted, the part that cannot be absorbed is the
    cone-constraint residual already on record: 0.45 to 0.67 of the forcing
    RMS, rising with N from 128 to 512.
- **The local form of the statistic has no response that grows with N.**
  - That form correlates the signs of an event's flippable relations with
    the relative direction of the partners.
  - An event's own relations absorb the wave exactly, through a
    unit-Jacobian shift in `v`.
  - The partner relations that make a relation flippable absorb it too,
    except in overlap and boundary channels that have no gain in N.
  - **Stop rule (i) fires for that form, conditionally.** The zero was
    argued, not computed. The referee sketched a preregistered finite-N
    test that would compute it.
- **Consequence for next_question.md.**
  - Both "per event only" rows of the heuristic table lose their premise.
  - One of them, the count-limited row, is one of the three accountings
    behind the "coincidence at 2". The other two are untouched.
- **The joint per-event fit is a new candidate, not a surviving form.**
  - It has no likelihood yet: free per-event shears break transitivity.
  - Its resolution is set by U2. An event's walls sit at macroscopic
    distance along its cones, so pinning them against each other is
    transport across the diamond.
  - If transport is count-like, the scaling is about `alpha > 8/3`. If it
    is cell-like, the phase is known and the relevant row is 4/3. These are
    scalings, not results.

By next_question.md's rule a stopped route is a result. It is not converted
into a smaller claim.

## 2. Where Routes C and D now meet

Both routes now turn on U2: how finely the order fixes the **relative
position of events far apart**, beyond the global seven-parameter group.

Route D's stop rule already frames it:

- a spread of `sqrt(d / N)` gives the heuristic 10/7;
- anything smaller moves the target toward 4/3.

Route D's first measurement and its extension measured the residual of
**central** events. That is how finely positions are fixed locally. It is
not transport.

## 3. Proposed next step for Route D (a proposal, not a design)

**Why not the repaired far-start test.** By
[../route_d/RECORD.md](../route_d/RECORD.md) Sections 3 and 4, at the sizes
it can afford it would tell a width of one plateau from three or more. It
needs rebuilt controls and 24 to 48 orders. Its "forgotten" outcome means
narrow relative to the push, not cell scale. It addresses the chain's reach
rather than U2.

**Proposed instead: a U2 transport measurement with Route D's frozen
instrument.**

- **Instrument.** Unchanged: the same order-preserving chain, the same
  seven-parameter group fit, the same state pairs and the same unit.
- **New readout.** The posterior spread of the **relative** position of
  event pairs, after the group fit, as a function of their separation `d`,
  from neighbours to across the diamond.
- **Separate the coordinate that carries the wave phase.** The wave
  coordinate `(T + Z)/sqrt(2)` was not shown flat in lag by the first
  measurement.
- **Descriptions fixed in advance**, in the style of the first measurement:
  - "flat in `d`" (cell-like transport);
  - "growing with `d`", with the `sqrt(d / N)` coefficient bounded as the
    count-scale term was;
  - size dependence over N = 128 to 512, or 1024 if affordable;
  - fresh seeds;
  - a power section under the house rule if any accuracy claim is made.
- **Known limits, carried over.**
  - The reading is one-sided: a chain can under-report collective motion.
    The 1+1 control showed about three quarters at `N^2/2`.
  - The fit absorbs a common shift of interior events.
  - It says nothing about the wave itself.
- **Cost.** Same order as the extension: about two hours per order at
  N = 512 with twelve workers. A size ladder with 12 to 24 orders per size
  would take half a day to a day of CPU at about 75%.

**What already exists, and why it is not evidence.** Outer events are not
used in the group fit. Their residual relative to the frame fitted on
central events was, in units of the one-sweep displacement at lag `N^2/2`:

| N | outer-event residual (one-sweep displacements) | chart units |
|---|---|---|
| 128 | 6.95 [6.11, 7.78] | about `20 / N` |
| 192 | 7.39 [6.68, 8.10] | about `21 / N` |
| 256 | 7.45 [6.96, 7.94] | about `22 / N` |
| 512 | 8.46 [7.98, 8.95] | about `25 / N` |

- Count-like transport would double this ratio from N = 128 to 512.
  Cell-like would leave it flat. It rose by a factor of about 1.2.
- But it mixes transport with the outer events' own looser cells.
- No description was fixed for it.
- And it was read after the records were seen.

So it motivates the design and cannot be its evidence.

**Also available for Route C.** The referee's finite-N test of the local
form ([referee_report.md](referee_report.md), last section). It would turn
the conditional stop into a computed one. It is cheaper than the transport
measurement, but it cannot move the bracket.

Neither has been drafted in full. Drafting, freezing and running each need
the user's approval.

## 4. Outcome of the transport draft (2026-11-26)

The user chose the U2 transport measurement as the next draft. It was
written ([../route_d/transport_design_draft.md](../route_d/transport_design_draft.md))
and reviewed
([../route_d/transport_referee_report.md](../route_d/transport_referee_report.md)).
It is **not to be frozen as written**. Its descriptions would label
cell-scale correlation and the composition of the event set as transport,
and its readout control could not fail. Details are in
[../route_d/RECORD.md](../route_d/RECORD.md) Section 5.

The review also showed that the existing record already speaks to U2 inside
the central set. A transport of the heuristic's size is not seen there by
this instrument, at the lags reached. That is one-sided: it does not move
either route's branch. Its bearing on Route C's joint candidate is the
heuristic's "U2 none" case: count-like resolution, about `alpha > 8/3`.
That case is not supported at the central scale. The cell-like case, which
would make the relevant row 4/3, is not shown either.
