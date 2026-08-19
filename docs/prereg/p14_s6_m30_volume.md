# S6 rung M = 3.0 — oracle freeze (certified volume at target 0.005)

Freeze document (2026-08-19, strong-curvature ruling). **No results.**
Execution of this rung **is approved** (PI ruling: push the ladder once
more toward strong curvature, one rung, seeds allocated for the pilot
and count stages only). The run starts once, from the clean exact
checkout of the freeze branch head — the merge commit's second parent,
never the merge commit.

## 1. The rung

The ladder's **strong-curvature rung and its deep end**. The certified
shell `[10, 20]`, the cap and the anchors `(12, 18)` stay **fixed in
absolute coordinates**; only the mass changes — nothing is rescaled
with M, so this is not an isometric copy of any other rung:

| quantity | value | source |
| --- | --- | --- |
| `m` | **3.0** | strong-curvature ruling |
| `mu = 2M/r_c` (r_c = 15) | **0.4** | the pre-frozen indicator (central rung: 0.1333…; ladder spans a factor of three) |
| `T_min(3.0)` | 10.158883083359672 | `rho(18) − rho(12)` at m = 3.0 |
| `dt` | **12.442423039673733** | the ladder convention `8.5 · T_min(M)/T_min(1)`, one exact path (`s6_rungs`) |
| box `r` | [11.442907037276841, 18.769059682298728] | containment certificate, outward binary64 |
| `psi_max` | 0.39909090507371286 | `dt / (2 w_glob)` |
| `SCALE` | 10472.024224793164 | `dt · B_box`, the o4_sizing exact path |

At this compactness the inner anchor sits at `r = 4M`, between the
photon sphere (`3M = 9`) and the ISCO (`6M = 18`), with the shell
floor `R_MIN = 10` above the photon sphere.

**Why this is the deep end.** Every lemma row certifies exactly on
`M ∈ [0.92, 3.33]`: below, the L5 winding cross-check fails; above,
L1/L2c fails because the photon sphere `3M` reaches the shell floor
`R_MIN = 10` at `M = 10/3`. `w_monotone` is the only margin that
shrinks with depth and here still clears by **1.0** — a full tenth of
`R_MIN` — where M = 3.2 / 3.3 would leave 0.4 / 0.1 for at most 10%
more `mu`. Every other margin grows with depth.

The runner refuses to **import** unless (a) `dt`/`mu` equal the ladder
constants re-derived through the one exact path, and (b) the rung's
full certified lemma table passes — exterior, K < 0, w-monotone,
Q > 0 via the perihelion floor, the L2a patch bound, all four L4
margins, and the L5 winding cross-check (margins
+2.28/+6.58/+3.57/+18.73 for L4, +70.20 for winding). Nothing is
inherited silently from M = 1 or from any shallower rung
(certification addendum, S6 section).

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
is**. Cost projection (not a certification): the executed rungs cost
565,478 calls / 14.6 h at M = 1, ~507k / 14.3 h at M = 1.4 and ~503k
/ 12.4 h at M = 1.8. This rung's box is smaller again (SCALE 10472.02
against 13679.09), so a comparable count is the projection — but
stronger curvature may refine harder, and that risk is carried by the
caps, not by a raise.

## 3. What the artifact publishes

- The **standalone** certified interval, outward-serialized. There is
  NO intersection base at this mass: no prior certified interval
  exists, and the design extrapolation is not a certification and
  never becomes a gate. Cross-rung consistency is a statement the S6
  integration stage may make later — per-rung verdicts stay separate,
  no joint primary.
- `v_ref_rung` as the standalone midpoint, **recommendation only**:
  adoption and every count-stage number (A, acceptance window, pilot
  n, powers) are re-derived at 96-bit from the actual endpoints at
  this rung's count freeze. In particular the design extrapolation's
  `V ≈ 114`, `A ≈ 469` and `n ≈ 4.9M` are **NOT frozen anywhere** —
  they are a cost projection only, and the extrapolation is over a
  wider mass gap than any executed rung.
- The rung identity block (m, mu, SCALE, ladder) and the certified
  lemma-table snapshot the import gate verified.

## 4. Discipline

Deterministic certified integration: **no seed, no sprinkling, no
reservation** — the write-once artifact is the only exclusion needed.
The two allocated seeds (`s6_m30_pilot = 40,000,501`,
`s6_m30_count = 40,000,511`) belong to the later stages and are
untouched here; `40,000,301` remains unallocated and unspent.
`verify_freeze` runs at entry AND exit; the exit lineage compares to
the **approved SHA directly** (the PR #81 R3 rule); `--freeze-rev`
must be the full 40-hex SHA because the manifest cannot certify
itself; a `target-not-met` or any incident stops the rung arc for a
PI ruling — never a retry, never a cap raise.
