# O4 incident — protocol abort at G3, no scientific verdict

This document records an **event**, not a result. The machine-readable
record is `p14_o4_incident.json`. Nothing here upgrades, grades, or
partially reports the O4 audit: the campaign produced no verdict, and
no gate may be described as having passed.

## 1. What happened

The approved campaign ran from the clean exact checkout of the freeze
merge `1eb9461`, claimed its streams on `refs/o4/reservation`
(`c4da162`), completed G1's 26,200,000 points and G2's function call,
entered G3, and stopped on the frozen fail-closed path:

```
fail-closed: causal_relation returned undecided at a G3 stress point
             -- no result is published
EXIT=1
```

Operator record: launched 2026-08-12 10:40:35 +09:00, exited 22:48:22
+09:00 (43,667 s). Neither frozen cap bound — the 24 h wall cap and the
80 M call cap were both far off. **The termination was the protocol
working, not the budget running out.**

## 2. What is unavailable, and why that matters

G1 and G2 finished *executing*. Their statistics were computed in
memory and **never persisted**, because the v1 runner writes its single
artifact only after all three gates return. So:

- there is no G1 mean, variance, half-width or interval,
- there is no G2 leakage bound,
- and therefore **no gate has a status** — not pass, not fail, not
  inconclusive.

An earlier report in the working session said G2 "appears to have
passed". That is retracted: function completion is not gate passage,
and the distinction is exactly what an audit of this kind exists to
respect.

The failing G3 stress point is also unavailable. The v1 runner records
no coordinates, no per-leg `dt`/`T_min`/`err`, and no returned
tri-state, so the *identity* of the undecided point cannot be recovered
from this run. That is an observability defect of the runner, listed as
such, and it is the first thing the replay diagnostic fixes.

**The lesson the runner must carry:** fail-closed has to mean *no
scientific verdict is published*, never *no record is left*. A run that
consumed 12 hours and retired two campaign scalars must leave a
write-once incident artifact of its own accord. This one had to be
reconstructed by hand from a launch-wrapper log.

## 3. First reading of the cause — a G3 contract defect, derived from
the frozen code

This section is **algebra on the frozen sources**, not analysis of run
data. No execution is involved, so it is stated with certainty about
the *conditions*, and says nothing about how often they occur.

`causal_relation` is deliberately **tri-state** (`s1_schwarzschild_cost.py`):

```
t_min, err = flight_time(p_r, q_r, dpsi, m, tol)
if abs(dt - t_min) <= err:
    return None                     # undecided BY DESIGN
return bool(dt > t_min)
```

`run_g3` probes each cluster at two times and treats **any** `None` as
an abort. Substituting the probe definitions gives the undecided
condition in closed form. With `lo = T1`, `hi = dt - T2`, `L = hi - lo`:

| probe | leg | `|dt - T_min|` | undecided iff |
|---|---|---|---|
| midpoint `(lo+hi)/2` | `p → x` | `L / 2` | `L/2 ≤ err₁` |
| midpoint `(lo+hi)/2` | `x → q` | `L / 2` | `L/2 ≤ err₂` |
| outside `hi + 1e-6` | `p → x` | `L + 1e-6` | `L + 1e-6 ≤ err₁` |
| **outside `hi + 1e-6`** | **`x → q`** | **`1e-6` exactly** | **`err₂ ≥ 1e-6`** |

The last row is the important one, and it must be stated precisely.
What is geometry-independent is the **decision distance**: the outside
probe's second leg sits a hard-coded `1e-6` from the boundary it is
testing, whatever the cluster looks like, so unlike every other row it
does not shrink with the window. The undecided condition reduces to
`err₂ ≥ 1e-6` — a comparison between a fixed offset and the solver's
reported bound, with no `L` in it.

That is not the same as saying the outcome is geometry-independent.
`err₂` is the solver's own reported bound at that leg, inflated by the
angle-mismatch term `|ang - dpsi| · b` and by `tol`, and it therefore
depends on the path and hence on the geometry. The precise statement
is: **the fixed offset does not adapt to the local `err₂`; the
boundary gap is fixed regardless of geometry, but whether that gap is
resolved is decided by an `err₂` that geometry controls.** G3 never
consults `err` at all, so it cannot notice either half of this.

The midpoint rows add a second, geometry-dependent path: G3's clusters
are selected by `L_S1 > 0` alone, which admits arbitrarily thin
windows, and a thin enough window makes `L/2 ≤ err` on its own.

So the frozen requirement — zero undecided among 100,000 clusters — was
never guaranteed by the sample definition it was built on. The
fail-closed device fired correctly; what it caught was a **contract
mismatch between a boolean-only stress test and a deliberately
tri-state predicate**, not evidence about S1's correctness. Whether S1
is defective remains open and is not addressed by this run.

## 4. Diagnosis plan (replay only, no new observation)

To be committed after this record merges. It replays the **already
spent** scalar `40_000_281` via the ledger's explicit `replay_scalar`,
which is a reproduction and never a new draw:

1. no G1 statistic is computed or reported;
2. stop as soon as the leading 100,000 positive-window points are
   collected;
3. at every stress point record: index, `(r, θ)`, `T₁`, `T₂`, `L`,
   probe kind (midpoint / outside), each leg's `dt`, `T_min`, `err`,
   `dt - T_min`, the actual `causal_relation` return value, and a cause
   classification — `midpoint-in-error-band`,
   `outside-in-error-band`, `boolean-mismatch`;
4. report **counts per cause**, not merely the first failure;
5. a tolerance ladder, if used, is diagnostic only and is never
   promoted into a gate.

## 5. Redesign direction (not yet frozen)

Splitting G3 into two contracts is the clean form:

- **G3a wrapper consistency** — the expected value is drawn from
  `{True, False, None}` and compared to what is returned. A legitimate
  `None` is not a failure.
- **G3b decisiveness** — the undecided-cluster fraction is reported
  separately, or judged against a pre-specified margin.

If a boolean-only stress is wanted, eligibility must use the returned
errors at each point — admit only clusters whose
`[T₁ + err₁ + η, Δt - T₂ - err₂ - η]` is non-empty, and place the
outside probe at least `err₂ + η` beyond the edge. That changes the
sampling distribution, so it reopens the preregistration rather than
being a patch.

## 6. Consequence for the audit

After the diagnosis, O4 must be re-run **in full, under a new freeze,
with new G1/G2 scalars**. This run's raw samples and statistics were not
preserved, and reconstructing them by replay and then promoting them
into a scientific result is not permitted: the streams are retired, and
a replay is a reproduction, not an observation.
