# O5 Poisson-count campaign — freeze (sprinkled tri-state membership count)

Freeze document (2026-08-17, delegation ruling). **No results.**
**Execution is NOT approved by this freeze**: per the delegation, the
arc STOPS after this PR merges, and the campaign executes only after
ONE integrated PI approval — once, from the clean exact checkout of
the freeze branch head (the merge commit's second parent), and never
twice. No rerun, no seed reuse, no sample extension, no cap raise.

## 1. What it measures, and why

The prediction-anchored test the oracle arc was built for: a genuine
4D Poisson process of frozen intensity `A` per unit certified 4-volume
is sprinkled into the L4 box, and the count of points falling causally
between the frozen anchors is compared against `A · V′` for the O3′
certified interval. The **estimand is `V_op − V_true`**, where
`V_op = A⁻¹·E[K_op]` is the measure of the wrapper-accepted region:
the campaign asks whether the operational instrument (sampler +
tri-state predicate) realises the certified volume. `K ~ Poisson(A·V_true)`
is assumed ONLY under the defect-free alternative, which is what every
power number below conditions on — the verdict itself assumes nothing.

This is the first *prediction-anchored* stage: O4b's recovered
CONCORDANT was an audit of `sampler + S1 volume response + oracle`
agreement; this campaign confronts a sprinkled COUNT with the
certified prediction `ρ·V`.

## 2. Frozen configuration

| key | value |
| --- | --- |
| `a_intensity` | **940.0** per unit 4-volume (final; re-frozen from the actual O3′ endpoints) |
| `tau` | **0.025** (equivalence band `B = tau · V_ref′ = 1.419420321666095`) |
| `alpha` | 0.05 (exact Garwood, alpha/2 per side, ONE outer interval) |
| `v_lo`, `v_hi` | **56.49295916680885, 57.06066656647874** (the O3′ artifact verbatim, target-met at ratio 0.0049995) |
| `v_ref` | **56.776812866643795** (the O3′ artifact's recommendation, adopted HERE as deferred by the O3′ results PR) |
| `u_power_ceiling` | 30 (certified by the ambiguity pilot, FEASIBLE at k = 0) |
| seed | `o5_campaign = 40,000,451` (FRESH; spent by the run) |
| `tol` | 1e-08 (the campaign predicate's) |
| `chunk` | 65,536 |
| `scan_call_cap` | **60,000,000** |
| `max_calls` | **60,006,537 = scan cap + 6,537** (the deterministic full-wrapper G3a preflight, charged to the same budget) |
| `max_wall_s` | 86,400 |

Caps are frozen — **no auto-raise, ever**. The wall cap is a stop
rule, not part of any power statement. `E[calls] ≤ 2·A·SCALE =
53,667,552` (~11.5 h at the design price of 768 µs/pair; the pilot
measured ~19% slower, projecting ~13.7 h — both far under the wall).

The frozen endpoints are not merely asserted: the module refuses to
import if `p14_o3p_volume.json` does not carry exactly these
serialized endpoints with `target-met`, or if `p14_o5_amb_pilot.json`
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
D        = [C_lo − V′_hi,  C_hi − V′_lo]       (identified discrepancy)
B        = tau · V_ref′
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
acceptance on K_true (U = 0)     [52,752, 53,980]
   (equivalently: decide(k, 0) = CONCORDANT exactly for k in that window)
endpoint powers                   0.936665 @ V′_lo   0.930794 @ V′_hi
worst power                       0.930794  (≥ the 0.93 proposal margin)
power at U = 30 (window shrink)   0.911855
U where the floor would break     ≥ 47 (0.900286 at 46, 0.899527 at 47)
joint completion-and-concordance  ≥ 0.911855 − 0.001 (pilot tail budget)
                                    − 0.010 (pilot alpha)
                                    − exp(−179,855) (call-cap Chernoff)
                                  = 0.900855 ≥ 0.90
```

The Chernoff term is an exact INEQUALITY upper bound on
`P(2N > 60,000,000)` under `N ~ Poisson(A·SCALE)`, not the exact
tail. (The pilot-era quote exp(−179,491) assumed the 6,537 preflight
calls charged inside a 60M total; this freeze allocates
`max_calls = 60M + 6,537` so the scan's cap is exactly 60M, which
tightens the exponent. Both are astronomically beyond relevance and
the joint bound is unchanged at 6 decimals.) **Wall-cap risk is NOT
included** in this bound (host-dependent, not certified).

## 5. Execution order and failure discipline

static preflight (digests, environment lock, clean tree at the exact
full-40-hex `--freeze-rev`, artifact/incident/checkpoint absence,
fresh-seed assertion, retained-ref verification for `refs/o4` AND
`refs/o4b` AND `refs/o5pilot`, namespace probe on `refs/o5/`) →
metered G3a full-wrapper preflight (6,537 deterministic calls,
re-asserted at run time; a failure here has spent nothing) →
reservation claim on `refs/o5/reservation` (the claim RETURNS its
object; nonce-unique attempt; any uncertain push is SEED POSSIBLY
SPENT, fail-closed) → RNG construction → `N ~ Poisson(A·SCALE)` drawn
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
  fresh seed, retained refs and an unclaimed `refs/o5/reservation`.
- **Lineage.** The campaign stream (`40,000,451`) is fresh and
  lineage-independent of the pilot stream (retired, OBSERVED) and of
  every prior scalar; the pilot data never enters the estimator — it
  contributed one number, the certified U ceiling, to the DESIGN.
- **Claim discipline.** A CONCORDANT outcome is a statement that the
  operational instrument realises the O3′ certified volume within
  `tau` at the frozen coverage — the first prediction-anchored
  sentence of the program. It does not retroactively upgrade S4/S5,
  and nothing in Paper A changes with this freeze.
