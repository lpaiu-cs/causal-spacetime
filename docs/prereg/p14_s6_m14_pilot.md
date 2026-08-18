# S6 M=1.4 ambiguity pilot — freeze (fixed-n tri-state indicator scan)

Freeze document (2026-08-18, quota-delegated S6 ruling). **No
results.** The pilot executes ONCE after this freeze merges and the
exact-checkout preflight passes, from the freeze branch head — and
never twice. No rerun, no seed reuse, no sample extension, no cap
raise.

## 1. What it measures, and why

The M = 1.4 rung's count-stage power sentence needs a certified
ceiling on the undecided-event count U, exactly as the central rung's
did (docs/prereg/p14_o5_amb_pilot.md). The pilot draws a fixed-n
sample from the SAME rung box measure the M = 1.4 count campaign will
sprinkle, asks the full tri-state `causal_relation` **at m = 1.4**
about both legs of each event (early exit), and counts the events
where any asked leg returned None. Deeper curvature is exactly where
solver ambiguity could rise — this pilot is that gate: INCONCLUSIVE
stops this RUNG for a PI ruling and indicts nothing; the other rungs
proceed independently (verdict separation).

## 2. Frozen configuration — sized at 96-bit from the ACTUAL rung endpoints

| key | value |
| --- | --- |
| rung | m = 1.4, μ = 0.18666666666666665, dt = 9.070565190742672 |
| box | r ∈ [11.367185829820095, 18.7053462469963], ψ ≤ 0.623439673163968 (import-gated vs the one exact path) |
| `SCALE_rung` | 18141.08004658374 |
| `n_events` | **6,024,777** (fixed; never extended — the smallest n keeping the k = 1 boundary FEASIBLE with a positive margin at A_provisional·SCALE_rung, mirroring the central pilot's boundary choice; margin 3.98e-9 of tail) |
| `u_max` | 30 |
| `alpha_pilot` | 0.01 (exact one-sided Clopper–Pearson) |
| `tail_budget` | 0.001 |
| `a_provisional` | 830.0 (from the rung oracle's actual endpoints V = [64.336198, 64.982744]; final A at the count freeze) |
| seed | `s6_m14_pilot = 40,000,461` (FRESH; spent by the run) |
| `tol` / `chunk` | 1e-08 / 65,536 |
| `max_calls` | **12,056,071 = 2·n_events + 6,517** (the rung's own deterministic full-wrapper G3a preflight at m = 1.4 AND the rung geometry — measured, not copied from M = 1's 6,537) |
| `max_wall_s` | 86,400 |

Caps are frozen — **no auto-raise, ever**. Reference rows of the
general-k rule at this n (the RULE decides, this table only
illustrates): k=0 → tail 1.4749e-06, FEASIBLE; k=1 → 9.99996e-04,
FEASIBLE; k=2 → 2.4261e-02, INCONCLUSIVE.

## 3. The frozen decision rule — general in k

The central pilot's rule verbatim, with this rung's constants; the
deciding arithmetic runs at 96-bit MPFR precision (the same engine,
digest-pinned):

```
p_upper  = exact one-sided CP upper at confidence 1 − alpha_pilot from (k, n)
lambda_U = a_provisional · SCALE_rung · p_upper
tail     = exact P(Poisson(lambda_U) > u_max)
verdict  = FEASIBLE  iff tail ≤ tail_budget, else INCONCLUSIVE
```

## 4. Execution order and failure discipline

The O5 pilot's discipline verbatim (write-once artifacts, publication
commit decided from the record, publish-time re-verify uncertainty
under its own key, clean pre-push refusals filing no incident, chunk
atomic non-verdict checkpoints, incident naming its failure_point
with per-observation live tallies), with this rung's identities:
reservation `refs/s6m14pilot/reservation` with **retained refs × 4**
(o4, o4b, o5pilot, o5) verified before any claim; G3a wrapper
preflight at the rung mass AND rung geometry (the foundation's
one-rung invariant), 6,517 metered calls re-asserted at run time.
FEASIBLE and INCONCLUSIVE are published identically, as computed.
Nothing here reruns, extends, or relaxes; **Wall-cap risk is NOT
included** in any power statement.
