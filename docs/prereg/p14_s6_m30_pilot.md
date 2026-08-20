# S6 M=3.0 ambiguity pilot — freeze (fixed-n tri-state indicator scan)

Freeze document (2026-08-20, strong-curvature ruling: A = 440,
n = 1,843,669 and the seed allocation approved by the PI). **No
results.** The pilot executes ONCE after this freeze merges and the
exact-checkout preflight passes, from the freeze branch head — and
never twice. No rerun, no seed reuse, no sample extension, no cap
raise.

## 1. What it measures, and why

The M = 3.0 rung's count-stage power sentence needs a certified
ceiling on the undecided-event count U, exactly as every shallower
rung's did (docs/prereg/p14_o5_amb_pilot.md, p14_s6_m14_pilot.md,
p14_s6_m18_pilot.md). The pilot draws a fixed-n sample from the SAME
rung box measure the M = 3.0 count campaign will sprinkle, asks the
full tri-state `causal_relation` **at m = 3.0** about both legs of
each event (early exit), and counts the events where any asked leg
returned None. This is the ladder's deepest curvature — the inner
anchor at r = 4M, between the photon sphere and the ISCO — which is
exactly where solver ambiguity could rise; this pilot is that gate:
INCONCLUSIVE stops this RUNG for a PI ruling and indicts nothing (the
executed rungs are unaffected — verdict separation).

## 2. Frozen configuration — sized at 96-bit from the ACTUAL rung endpoints

| key | value |
| --- | --- |
| rung | m = 3.0, μ = 0.4, dt = 12.442423039673733 |
| box | r ∈ [11.442907037276841, 18.769059682298728], ψ ≤ 0.39909090507371286 (import-gated vs the one exact path) |
| `SCALE_rung` | 10472.024224793164 |
| `n_events` | **1,843,669** (fixed; never extended — the smallest n keeping the k = 1 boundary FEASIBLE with a positive margin at A_provisional·SCALE_rung, the same convention the M = 1.4 and M = 1.8 pilots used; margin 1.65e-8 of tail) |
| `u_max` | 30 |
| `alpha_pilot` | 0.01 (exact one-sided Clopper–Pearson) |
| `tail_budget` | 0.001 |
| `a_provisional` | 440.0 (from the rung oracle's actual endpoints V = [121.225021, 122.443280] — the smallest 10-step A whose count-stage worst endpoint power clears 0.93 (0.931493), keeping the ladder's power profile uniform; final A at the count freeze) |
| seed | `s6_m30_pilot = 40,000,501` (FRESH; spent by the run) |
| `tol` / `chunk` | 1e-08 / 65,536 |
| `max_calls` | **3,693,804 = 2·n_events + 6,466** (the rung's own deterministic full-wrapper G3a preflight at m = 3.0 AND the rung geometry — measured, not copied; m14 was 6,517, m18 was 6,461) |
| `max_wall_s` | 86,400 |

Caps are frozen — **no auto-raise, ever**. Reference rows of the
general-k rule at this n (the RULE decides, this table only
illustrates): k=0 → tail 1.4749e-06, FEASIBLE; k=1 → 9.99983e-04,
FEASIBLE; k=2 → 2.4260e-02, INCONCLUSIVE.

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
reservation `refs/s6m30pilot/reservation` with **retained refs × 8**
(o4, o4b, o5pilot, o5, s6m14pilot, s6m14, s6m18pilot, s6m18) verified
before any claim; G3a wrapper preflight at the rung mass AND rung
geometry (the foundation's one-rung invariant), 6,466 metered calls
re-asserted at run time. FEASIBLE and INCONCLUSIVE are published
identically, as computed. Nothing here reruns, extends, or relaxes;
**Wall-cap risk is NOT included** in any power statement. The count
seed `s6_m30_count = 40,000,511` is NOT registered by this freeze —
it enters the ledger at the count freeze; `40,000,301` remains
unallocated and unspent.
