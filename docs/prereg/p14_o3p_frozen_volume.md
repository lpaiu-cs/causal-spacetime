# O3′ freeze — re-certified diamond volume at target ratio 0.005

Freeze document (2026-08-16). **No results.** Execution is NOT approved
by this freeze: after the freeze PR review converges, the campaign run
needs a separate PI approval and an exact-checkout preflight. O3′ is
deterministic certified integration — no seed, no sprinkling, no
stochastic draw anywhere.

## 1. Why O3′, and why r = 0.005

The O5 Poisson-count stage is infeasible against the current O3
interval: its full width `W = 1.135282` sets an equivalence floor
`W/V_ref = 2.0%`, and the exact three-engine-reproduced sizing puts
even `tau = 3.0%` beyond the O5 campaign call caps. The adopted target
(PI ruling on the O5 sizing audit) is the back-computed minimum that
strengthens the O5 allowed sentence within caps:

- `r = 0.005` lowers the oracle floor to 1.0% and makes
  `tau_O5 = 2.5%` feasible at ~5.4e7 expected predicate calls;
- `r = 0.0075` would only preserve `tau = 3.0%` (no sentence gain);
- `r = 0.0025` would allow `tau = 2.0%` but its own projected cost
  exceeds every existing oracle cap under every fit window.

## 2. Frozen configuration

| key | value | grade |
| --- | --- | --- |
| anchors | (12, 18, 8.5)M, M = 1 | unchanged (the frozen configuration) |
| `target_ratio` | **0.005** | adopted target (PI ruling) |
| `n_sub` / `k_micro` / `d_switch` | 16 / 4 / 0.25 | inherited from O3 unchanged — a freeze choice, not optimality |
| `init_rho` × `init_psi` | 12 × 12 | inherited from O3 unchanged |
| `max_depth` | 18 | O3 value unchanged (O3 reached depth 12) |
| `max_calls` | **900,000** | frozen cap — **no auto-raise, ever** |
| `max_wall_s` | **86,400** | frozen cap — **no auto-raise, ever** |

Runner: `experiments/oracle/o3p_frozen_volume.py`. Freeze identity is
content-addressed (`p14_o3p_freeze_manifest.json`: raw SHA-256 of every
protocol file plus the python/gmpy2/MPFR/GMP environment lock),
verified at entry and exit; `--preflight` additionally requires a clean
tree, the absence of the O3′ result, and the presence of the O3 base.

## 3. Cost projection — two curves, fitted separately

Model `calls ~ C0 · (1/ratio)^p`, OLS on log-log points per window;
recomputed by `experiments/oracle/o3p_projection.py` and pinned by a
contract test. **Projections, not certifications.**

| source | window | p | projected calls at r = 0.005 |
| --- | --- | --- | --- |
| neighbor ladder | all 5 points | 1.793 | ~477,688 |
| neighbor ladder | last 4 | 1.846 | ~495,806 |
| neighbor ladder | last 3 | 1.871 | ~504,459 |
| neighbor ladder | last 2 | 2.435 | ~745,501 |
| O3 own curve (1,078 samples) | full curve | 1.796 | ~467,365 |
| O3 own curve | ratio ≤ 0.05 | 1.938 | ~545,681 |
| O3 own curve | ratio ≤ 0.03 | 1.944 | ~549,243 |
| O3 own curve | ratio ≤ 0.02 | 2.276 | ~749,784 |
| central (p = 2 through the O3 point) | — | 2.000 | ~551,501 |

**Range: ~467k–750k calls (~12.6–20.2 h at the O3 call rate).** The
neighbor-ladder fits are anchored through the executed O3 point
(137,958 calls at ratio 0.009997); the O3-curve fits carry their own
intercepts. Two independently measured curves give the same range,
which is why the caps are set at 900,000 calls / 24 h — +20% above the
worst fitted window — and stay frozen whatever happens.

## 4. What the run publishes

Write-once artifact `docs/prereg/p14_o3p_volume.json` (atomic no-clobber
`os.link`, the O3 mechanism), carrying:

1. the **standalone O3′ interval** with outward binary64 endpoints,
   `status`, `termination_reason`, calls/cells/wall, mode decomposition,
   and the full refinement curve;
2. the **intersection with the immutable O3 interval** (the O3 artifact
   is digest-pinned in this freeze's manifest and is never written).
   Both certifications claim to enclose the same true volume, so an
   **empty intersection is a certification inconsistency**: it indicts
   at least one of the two, the artifact records it machine-readably,
   and **no downstream stage (O5 sizing, V_ref′) may consume either
   interval until it is resolved**. On consistency, downstream consumes
   the **standalone O3′ interval** — the intersection is the
   consistency check, not the estimator;
3. a `v_ref_prime_recommendation`, **gated on the consistency check**:
   on consistency, **V_ref′ = the midpoint of the standalone O3′
   interval** — a recommendation adopted, with every O5 number (A,
   acceptance range, powers) re-derived from the actual O3′ endpoints,
   at the O5 freeze. On a certification inconsistency the value is
   `null` and the status reads BLOCKED, so a downstream consumer that
   reads the artifact directly cannot lift a midpoint out of a failed
   record — the blocking invariant lives in the producer;
4. **`target-not-met` is published as-is**: the certified interval
   reached at the moment a cap fired, the exact termination reason, and
   the same provenance (config, environment, git state, curve). A cap
   firing is a fact about the budget, not about the geometry, and the
   partial certification is still a valid enclosure at its achieved
   ratio.

## 5. Boundaries

- The existing O3 artifacts (`p14_o3_volume.json`,
  `p14_o3_freeze_manifest.json`, `p14_o3_executed_freeze_manifest.json`)
  are **immutable**; O3′ has its own write-once path and its own
  manifest.
- This freeze does not touch the seed ledger — there is nothing to
  touch; the integrator draws nothing.
- Nothing in Paper A changes with this freeze; the O5 stage that would
  consume the O3′ interval has its own pending design ruling, and its
  freeze PR is separately gated.
- Execution approval is a separate PI ruling after this PR's review
  converges: approved SHA, exact clean checkout, `--preflight` PASS.
