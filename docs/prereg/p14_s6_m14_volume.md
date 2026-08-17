# S6 rung M = 1.4 — oracle freeze (certified volume at target 0.005)

Freeze document (2026-08-17, S6 ruling). **No results. Execution is
NOT approved by this freeze**: the run starts only after the M = 1.4
**integrated execution approval** (which also covers, in order, this
rung's ambiguity pilot and count campaign, each behind its own frozen
preflight), once, from the clean exact checkout of the freeze branch
head — the merge commit's second parent, never the merge commit.

## 1. The rung

The mass ladder's first new rung. The certified shell `[10, 20]`, the
cap and the anchors `(12, 18)` stay **fixed in absolute coordinates**;
only the mass changes — nothing is rescaled with M, so this is not an
isometric copy of the central rung:

| quantity | value | source |
| --- | --- | --- |
| `m` | **1.4** | S6 ruling |
| `mu = 2M/r_c` (r_c = 15) | **0.18666666666666665** | the pre-frozen indicator (central rung: 0.1333…) |
| `T_min(1.4)` | 7.405857442632261 | `rho(18) − rho(12)` at m = 1.4 |
| `dt` | **9.070565190742672** | the ladder convention `8.5 · T_min(M)/T_min(1)`, one exact path (`s6_rungs`) |
| box `r` | [11.367185829820095, 18.7053462469963] | containment certificate, outward binary64 |
| `psi_max` | 0.623439673163968 | `dt / (2 w_glob)` |
| `SCALE` | 18141.08004658374 | `dt · B_box`, the o4_sizing exact path |

The runner refuses to **import** unless (a) `dt`/`mu` equal the
ladder constants re-derived through the one exact path, and (b) the
rung's full certified lemma table passes — exterior, K < 0,
w-monotone, Q > 0 via the perihelion floor, the L2a patch bound, all
four L4 margins, and the L5 winding cross-check (margins
+1.66/+3.71/+3.03/+5.48 for L4, +19.46 for winding). Nothing is
inherited silently from M = 1 (certification addendum, S6 section).

## 2. Frozen configuration

| key | value |
| --- | --- |
| `target_ratio` | **0.005** (ladder-uniform; the floor arithmetic that opens τ = 2.5%) |
| `n_sub` / `k_micro` / `d_switch` | 16 / 4 / 0.25 (inherited from O3/O3′ unchanged — a freeze choice, not optimality) |
| `init_rho` / `init_psi` / `max_depth` | 12 / 12 / 18 |
| `max_calls` | **900,000** |
| `max_wall_s` | **86,400** |

Caps are frozen — **no auto-raise, ever**. If they fire, the artifact
is the certified interval reached at that moment plus
`target-not-met` plus the exact termination reason, **published as
is**. Cost projection (not a certification): O3′ measured 565,478
calls / 14.6 h at M = 1; the design survey converged faster per call
at higher mass, so the ~467k–750k band carries over as a projection.

## 3. What the artifact publishes

- The **standalone** certified interval, outward-serialized. There is
  NO intersection base at this mass: no prior certified interval
  exists, the design survey's coarse estimate is not a certification
  and never becomes a gate. Cross-rung consistency is a statement the
  S6 integration stage may make later — per-rung verdicts stay
  separate, no joint primary.
- `v_ref_rung` as the standalone midpoint, **recommendation only**:
  adoption and every count-stage number (A, acceptance window, pilot
  n, powers — all provisional in the design survey) are re-derived at
  96-bit from the actual endpoints at this rung's count freeze. In
  particular the survey's A = 820 and pilot n ≈ 5.9M are NOT frozen
  anywhere.
- The rung identity block (m, mu, SCALE, ladder) and the certified
  lemma-table snapshot the import gate verified.

## 4. Discipline

Deterministic certified integration: **no seed, no sprinkling, no
reservation** — the write-once artifact is the only exclusion needed.
`verify_freeze` runs at entry AND exit; the exit lineage compares to
the **approved SHA directly** (the PR #81 R3 rule); `--freeze-rev`
must be the full 40-hex SHA because the manifest cannot certify
itself; a `target-not-met` or any incident stops the rung arc for a
PI ruling — never a retry, never a cap raise.
