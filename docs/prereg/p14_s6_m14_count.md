# S6 M=1.4 Poisson-count campaign — freeze (sprinkled tri-state membership count)

Freeze document (2026-08-18, quota-delegated S6 ruling). **No
results.** The campaign executes ONCE after this freeze merges and
the exact-checkout preflight passes, from the freeze branch head —
and never twice. No rerun, no seed reuse, no sample extension, no
cap raise. Any exception (incident, DISCORDANT, INCONCLUSIVE,
inconsistency) stops this rung for a PI ruling.

## 1. What it measures, and why

The prediction-anchored test the oracle arc was built for: a genuine
4D Poisson process of frozen intensity `A` per unit certified 4-volume
is sprinkled into the L4 box, and the count of points falling causally
between the frozen anchors is compared against `A · V_rung` for the
rung's certified interval (M = 1.4, μ = 0.18666666666666665). The **estimand is `V_op − V_true`**, where
`V_op = A⁻¹·E[K_op]` is the measure of the wrapper-accepted region:
the campaign asks whether the operational instrument (sampler +
tri-state predicate) realises the certified volume. `K ~ Poisson(A·V_true)`
is assumed ONLY under the defect-free alternative, which is what every
power number below conditions on — the verdict itself assumes nothing.

This is the program's first *preregistered equivalence gate* between a
sprinkled count and a reproducible directed-rounding volume enclosure.
Poisson-sprinkled count verification on Schwarzschild causal sets has
prior art (Homšak–Veroni, arXiv:2404.11670 §V.1); the novelty here is
the preregistered TOST verdict against a certified MPFR interval whose
endpoints were frozen before any count data existed. O4b's recovered
CONCORDANT was an audit of `sampler + S1 volume response + oracle`
agreement; this campaign confronts a sprinkled COUNT with the
certified prediction `ρ·V`.

## 2. Frozen configuration

| key | value |
| --- | --- |
| `a_intensity` | **830.0** per unit 4-volume (final; re-frozen from the actual rung endpoints) |
| `tau` | **0.025** (equivalence band `B = tau · V_ref = 1.6164867722650458`) |
| `alpha` | 0.05 (exact Garwood, alpha/2 per side, ONE outer interval) |
| `v_lo`, `v_hi` | **64.33619774000562, 64.98274404119803** (the rung artifact verbatim, target-met at ratio 0.0049996) |
| `v_ref` | **64.65947089060182** (the rung artifact's recommendation, adopted HERE as deferred by PR #94) |
| `u_power_ceiling` | 30 (certified by the rung pilot, FEASIBLE at k = 1 — the pre-sized boundary case) |
| seed | `s6_m14_count = 40,000,471` (FRESH; spent by the run) |
| `tol` | 1e-08 (the campaign predicate's) |
| `chunk` | 65,536 |
| `scan_call_cap` | **34,000,000** |
| `max_calls` | **34,006,517 = scan cap + 6,517** (the rung-measured G3a preflight at m = 1.4 and the rung geometry, charged to the same budget) |
| `max_wall_s` | 86,400 |

Caps are frozen — **no auto-raise, ever**. The wall cap is a stop
rule, not part of any power statement. `E[calls] ≤ 2·A·SCALE_rung =
30,114,193` (~7.7 h hard; less at the observed early-exit rate —
far under the wall).

The frozen endpoints are not merely asserted: the module refuses to
import if `p14_s6_m14_volume.json` does not carry exactly these
serialized endpoints with `target-met`, or if `p14_s6_m14_pilot.json`
does not carry FEASIBLE at `u_max = 30`. The manifest pins those
artifacts' bytes; the import gate pins their meaning.

## 3. The frozen verdict rule — general in (K_certain, U_amb)

The scan classifies each sprinkled point by the tri-state membership
(leg1 `P → x`; None → ambiguous, False → out with leg2 never asked,
True → ask leg2 `x → Q`). For WHATEVER tallies it produces:

```
L        = exact Garwood lower at alpha/2 from K_certain
U        = exact Garwood upper at alpha/2 from K_certain + U_amb
C        = [L/A, U/A]                          (volume units)
D        = [C_lo − V_hi,  C_hi − V_lo]         (identified discrepancy)
B        = tau · V_ref
verdict  = CONCORDANT     if D ⊆ [−B, B]
           DISCORDANT     if D ∩ [−B, B] = ∅
           INCONCLUSIVE   otherwise
```

ONE outer interval handles the tri-state: every ambiguous point may or
may not be a member, so `[L(K_certain), U(K_certain + U_amb)]` covers
`A·V_op` at the frozen coverage whatever the ambiguous points are. No
normal approximation. **The deciding arithmetic runs at 96-bit MPFR
precision end to end** — the Garwood bisections, the division by `A`,
the band, every comparison — with no double-precision step on the
verdict path (the pilot's boundary lesson, applied fully). All three
outcomes are published identically, as computed; INCONCLUSIVE and
DISCORDANT are results, not failures.

Ambiguity-measure lineage: the campaign's ambiguity rule is **a
fortiori narrower** than the pilot's (the pilot always asked leg2
after a decided leg1; the campaign early-exits on leg1 False), so the
pilot's certified ceiling `P(U > 30) ≤ 10⁻³` at 99% confidence
transfers to the campaign's U on the same measure.

## 4. The power sentence (design; conditions on the alternative)

At the frozen `A`, `tau`, `alpha` and the actual O3′ endpoints, three
independent engines (NR-style double gamma, the frozen P14 P2 engine,
96-bit gmpy2) agree exactly on the acceptance integers and to 6
decimals on the powers; the contract tests re-derive them from the
runner's own 96-bit rule:

```
acceptance on K_true (U = 0)     [53,045, 54,282]
   (equivalently: decide(k, 0) = CONCORDANT exactly for k in that window)
endpoint powers                   0.937571 @ V_lo   0.932156 @ V_hi
worst power                       0.932156  (≥ the 0.93 proposal margin)
power at U = 30 (window shrink)   0.913550
joint completion-and-concordance  ≥ 0.913550 − 0.001 (pilot tail budget)
                                    − 0.010 (pilot alpha)
                                    − exp(−120,283) (call-cap Chernoff)
                                  = 0.902550 ≥ 0.90
```

The Chernoff term is an exact INEQUALITY upper bound on
`P(2N > 34,000,000)` under `N ~ Poisson(A·SCALE_rung)`, not the exact
tail. **Wall-cap risk is NOT
included** in this bound (host-dependent, not certified).

## 5. Execution order and failure discipline

static preflight (digests, environment lock, clean tree at the exact
full-40-hex `--freeze-rev`, artifact/incident/checkpoint absence,
fresh-seed assertion, retained-ref verification for `refs/o4`, `refs/o4b`, `refs/o5pilot`, `refs/o5` AND
`refs/s6m14pilot`, namespace probe on `refs/s6m14/`) →
metered G3a full-wrapper preflight (6,517 deterministic calls,
re-asserted at run time; a failure here has spent nothing) →
reservation claim on `refs/s6m14/reservation` (the claim RETURNS its
object; nonce-unique attempt; any uncertain push is SEED POSSIBLY
SPENT, fail-closed) → RNG construction → `N ~ Poisson(A·SCALE_rung)` drawn
once from the stream → chunked scan with an atomic non-verdict
checkpoint per chunk (rng_stream named; chunk provenance carried) →
reservation re-read → write-once publication. Every failure path —
KeyboardInterrupt included — files a write-once incident naming its
`failure_point` and preserving the partial tallies. A publish-time
re-verify uncertainty is recorded under its OWN key
(`publish_reverify_uncertainty`) with the stage as the failure point:
a confirmed claim stays a confirmed claim (the pilot's PR #85
lesson, built in from the start here).

**The publication commit is decided from the record, not from control
flow** (the O4b R26–R29 pattern): a receipt is stamped the instant
`os.link` returns, and the artifact is asked whether it names this
run's nonce-unique claim object. A failure after the commit exits
without filing an incident beside the published result; a write-once
refusal against a foreign artifact still files one, because a spent
seed with neither result nor record is the outlawed state.

## 6. Boundaries

- **One run.** Whatever comes out — CONCORDANT, DISCORDANT,
  INCONCLUSIVE, or an incident — is published as computed and the
  seed is spent. A DISCORDANT or INCONCLUSIVE outcome stops the arc
  for a PI ruling; nothing retries, extends, or relaxes.
- **Execution gate.** This freeze does NOT authorize execution. The
  run starts only after the integrated approval names the freeze
  branch head; the runner additionally refuses any checkout that is
  not exactly that SHA, clean, with matching digests, environment,
  fresh seed, retained refs and an unclaimed `refs/s6m14/reservation`.
- **Lineage.** The campaign stream (`40,000,471`) is fresh and
  lineage-independent of the pilot stream (retired, OBSERVED) and of
  every prior scalar; the pilot data never enters the estimator — it
  contributed one number, the certified U ceiling, to the DESIGN.
- **Claim discipline.** A CONCORDANT outcome is a statement that the
  operational instrument realises the rung's certified volume within
  `tau` at the frozen coverage — the program's first preregistered
  equivalence gate against a reproducible directed-rounding enclosure
  (Poisson-sprinkle verification on Schwarzschild has prior art in
  Homšak–Veroni 2404.11670 §V.1; see `literature_priority_oracle_arc.md`).
  It does not retroactively upgrade S4/S5, and nothing in Paper A
  changes with this freeze.
