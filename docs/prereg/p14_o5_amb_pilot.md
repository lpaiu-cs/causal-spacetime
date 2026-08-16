# O5 ambiguity pilot — freeze (fixed-n tri-state indicator scan)

Freeze document (2026-08-16, delegation ruling). **No results.** The
pilot executes ONCE after this PR merges, from the clean exact checkout
of the freeze branch head (the merge commit's second parent), with no
further per-step approval — and never twice. No rerun, no seed reuse,
no sample extension, no cap raise.

## 1. What it measures, and why

The O5 Poisson-count power sentence ("worst-case ≥ 0.90") holds only if
the campaign's undecided-event count U stays below the acceptance-window
shrink it causes. The pilot certifies the frozen ambiguity budget: it
draws a fixed-n sample of events from the SAME 4D box measure the O5
campaign will sprinkle, asks the full tri-state `causal_relation` about
both legs of each event (early exit), and counts the events where any
asked leg returned None.

## 2. Frozen configuration

| key | value |
| --- | --- |
| `n_events` | **10,736,965** (fixed; never extended) |
| `u_max` | 30 |
| `alpha_pilot` | 0.01 (exact one-sided Clopper–Pearson) |
| `tail_budget` | 0.001 |
| `a_provisional` | 940.0 (provisional; re-frozen at the O5 freeze) |
| seed | `o5_amb_pilot = 40,000,441` (FRESH; spent by the run) |
| `tol` | 1e-08 (the campaign predicate's) |
| `chunk` | 65,536 |
| `max_calls` | **21,480,467 = 2·n_events + 6,537** (the hard two-leg bound plus the deterministic full-wrapper G3a preflight, charged to the same budget) |
| `max_wall_s` | 86,400 |

Caps are frozen — **no auto-raise, ever**. The wall cap is a stop rule,
not part of any power statement.

## 3. The frozen decision rule — general in k

For WHATEVER ambiguous count k the fixed-n scan produces (nothing is
hardcoded to small k):

```
p_upper  = exact one-sided CP upper at confidence 1 − alpha_pilot from (k, n)
lambda_U = a_provisional · SCALE · p_upper        (Poisson thinning)
tail     = exact P(Poisson(lambda_U) > u_max)     (finite pmf sum)
verdict  = FEASIBLE  iff tail ≤ tail_budget, else INCONCLUSIVE
```

No normal approximation. The deciding functions are cross-checked by
contract tests against the frozen P14 P2 exact-Garwood engine and an
independent 96-bit gmpy2 engine. Reference evaluations at the frozen n
(the RULE decides, this table only illustrates): k=0 → λ_U 11.5092,
tail 1.475e-06, FEASIBLE; k=1 → λ_U 16.5905, tail 9.99999e-04,
FEASIBLE; k=2 → λ_U 21.0081, tail 2.426e-02, INCONCLUSIVE.

## 4. Budget allocation and the joint bound it certifies

With the sizing audit's exact arithmetic at the actual O3′ endpoints
(`V′ = [56.492959, 57.060667]`, A provisional 940):

```
power(U=30)                      0.911855
− tail_budget                    0.001
− alpha_pilot                    0.010
− campaign call-cap Chernoff     exp(−179,491)   (exact INEQUALITY upper
  bound — not the exact tail; 60M cap minus 6,537 deterministic preflight)
= joint lower bound              0.900855 ≥ 0.90
```

Wall-cap risk is NOT included in this bound (host-dependent, not
certified). Final A, acceptance range, and joint power are re-frozen at
the O5 freeze from the actual O3′ endpoints and this pilot's bound.

## 5. Execution order and failure discipline

static preflight (digests, environment lock, clean tree at the exact
full-40-hex `--freeze-rev`, artifact/incident/checkpoint absence,
fresh-seed assertion, retained-ref verification for `refs/o4` AND
`refs/o4b`, namespace probe) → metered G3a full-wrapper preflight
(6,537 deterministic calls, re-asserted at run time; a failure here has
spent nothing) → reservation claim on `refs/o5pilot/reservation` (the
claim RETURNS its object; nonce-unique attempt; any uncertain push is
SEED POSSIBLY SPENT, fail-closed) → RNG construction → fixed-n chunked
scan with an atomic non-verdict checkpoint per chunk (rng_stream named;
chunk provenance carried) → reservation re-read → write-once
publication. Every failure path — KeyboardInterrupt included — files a
write-once incident naming its `failure_point` and preserving the
partial tallies. FEASIBLE and INCONCLUSIVE are published identically,
as computed.

## 6. Boundaries

- **Fixed n.** A larger-than-hoped k publishes INCONCLUSIVE as-is; the
  sample is never extended, the thresholds never relaxed, the seed
  never reused. An INCONCLUSIVE pilot stops the arc (redesign is a PI
  ruling, never automatic).
- **Lineage independence.** The pilot stream (`40,000,441`) is separate
  from every campaign scalar; no O5 campaign seed, reservation, or
  freeze exists until this pilot's result is merged. A FEASIBLE pilot
  feeds one number forward — the certified ambiguity bound — and the
  pilot data never enters the O5 estimator.
- The verdict is about the FEASIBILITY of the O5 power sentence at
  `u_max`, not about the instrument: INCONCLUSIVE indicts nothing.
- Nothing in Paper A changes with this pilot.
