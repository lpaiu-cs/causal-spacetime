# P14 probe P1 — results (exploratory)

**Status: EXPLORATORY.** One run, seed `20260808`, no gate, no
threshold, no verdict. These numbers inform the P2–P4 designs and the
eventual freeze; they decide nothing by themselves (design §8.2).

Reproduce:

    python experiments/positive_control/p14_probe_p1.py 20260808

Machinery pinned by `tests/test_p14_probe_p1.py`; the numbers below are
this run's, not test-pinned. Convention: `w = 1`, `N = 300` elements
per sprinkling at every operating point — equal compute, not equal
density; density is the free knob (`rho = N / V_box`). 24 sprinklings
per point; MC volumes from 200k shared samples.

## The four answers

### 1. Eligible fraction — thin everywhere the effect lives

| point | `a` | `w·dv` | `w·dx/2, w·dy/2` | f_vol | f_elem | f_pair | elig. pairs/sprinkling |
|---|---|---|---|---|---|---|---|
| slice-a0.3 | 0.3 | 0.2 | 3, 3 | 0.0877 | 0.0821 | 0.0067 | 301 |
| slice-a0.6 | 0.6 | 0.2 | 3, 3 | 0.0508 | 0.0443 | 0.0019 | 86 |
| slice-a1.0 | 1.0 | 0.2 | 3, 3 | 0.0277 | 0.0232 | 0.0005 | 24 |
| aniso-a1.0 | 1.0 | 1.0 | 1, 3 | 0.0258 | 0.0235 | 0.0006 | 26 |
| roomy-a0.2 | 0.2 | 16.0 | 3, 3 | 0.3148 | 0.3058 | 0.0936 | 4200 |
| high-a2.0 | 2.0 | 0.5 | 1, 1 | 1.4e-4 | 1.4e-4 | ~0 | 0.0 |
| edge-a2.4 | 2.4 | 0.2 | 0.6, 0.4 | 2.4e-4 | 3.0e-4 | ~0 | 0.0 |

Measured element fractions track the analytic volume fraction (the
guard's own maximized quantity), which closes that consistency loop.
The pair fraction is its square, as the factorized rule implies.

### 2. Feasibility — the effect lives at large `a`, and external
anchors reach it

The effect for an on-axis diamond is a function of `a = w·Δu` ALONE:
slicing the diamond, `Δv` factors out of the volume ratio, and

    V_A / V_0 = R(a),   R(a) − 1 = a⁴/252 + O(a⁸),

computed by quadrature (never the series — at `a = 2.4` the true
`R − 1 = 0.181` against the series' `0.132`, and it diverges toward
the conjugate point `a = π`). The quadrature reproduces the design
check script's pinned values to the digit (`1.00400047` at `a = 1`,
`1.07300802` at `a = 2`), which retroactively identifies that script's
`T` as the u-separation.

Fattest eligible axis diamond per point; `n90` = sprinklings for 90%
power at two-sided 95% (house convention, §6.1), with the **predicted**
`δ` as signal — verifying the prediction is P2's charter:

| point | `a` | δ = R−1 | λ_A | V_dis/V_A | n90 marginal | n90 paired |
|---|---|---|---|---|---|---|
| slice-a0.3 | 0.3 | 3.2e-5 | 0.14 | ~0 | 7.1e10 | — |
| slice-a0.6 | 0.6 | 5.2e-4 | 0.26 | 0.012 | 1.5e8 | 1.8e6 |
| slice-a1.0 | 1.0 | 4.0e-3 | 0.42 | 0.050 | 1.6e6 | 7.6e4 |
| aniso-a1.0 | 1.0 | 4.0e-3 | 5.80 | 0.037 | 1.1e5 | 4.2e3 |
| roomy-a0.2 | 0.2 | 6.4e-6 | 12.9 | 0.002 | 2.0e10 | 3.5e7 |
| high-a2.0 | 2.0 | 7.3e-2 | 0.07 | 0.104 | 2.7e4 | 3.2e3 |
| edge-a2.4 | 2.4 | 1.8e-1 | 0.38 | 0.329 | 834 | **386** |

`λ_A` is the expected count in the diamond at `N = 300`;
`V_dis` is the volume where the two arms' predicates disagree,
measured by shared-sample MC — it is the exact per-realization
variance of the same-points paired count, `Var(N_A − N_0) = ρ V_dis`.

**Reading.** The product (effect)² × (count) is maximized at the
largest `a` the guard admits, despite the tiny boxes: δ grows as `a⁴`
and beyond, while the count only needs to be O(1) per sprinkling.
`edge-a2.4` needs ~4e2 sprinklings paired (~8e2 marginal); `roomy`
and the `slice` points are dead at any affordable n. Density is a
further free lever this table does not use: P2-style counting inside
an external-anchor diamond is O(N) per anchor, not O(N²), so raising
`N` at the small boxes cuts n linearly at modest cost.

**The split that matters downstream:**

- **External-anchor Class C (P2's volume check, anchors fixed and not
  counted, per the approval condition): FEASIBLE**, at large `a`.
- **Sprinkled-pair Class C (§5 items 2, 3b — chain vectors, MM on
  sprinkled intervals): starved exactly where the effect lives.**
  At `a ≥ 2` the run saw 0.0 eligible pairs per sprinkling; the
  viable middle is `aniso-a1.0`-like (26 pairs/sprinkling, δ = 4e-3).
  Whether that middle carries a usable statistic is §8 P3's question,
  with these numbers as its input.

### 3. Order-invariant candidate (§4.6.2) — does not reproduce the guard

Candidate: interior iff ≥ `k` elements below AND above. Element-level
agreement with the coordinate guard / candidate-only fraction (the
dangerous direction — diamonds that may leave the box):

| point | k=1 | k=2 | k=4 | k=8 | k=16 |
|---|---|---|---|---|---|
| slice-a0.3 | 0.90/0.02 | 0.92/0.00 | 0.92/0.00 | 0.92/0.00 | 0.92/0.00 |
| slice-a1.0 | 0.81/0.16 | 0.90/0.08 | 0.96/0.02 | 0.98/0.00 | 0.98/0.00 |
| aniso-a1.0 | 0.53/0.46 | 0.67/0.31 | 0.81/0.17 | 0.92/0.06 | 0.97/0.00 |
| edge-a2.4 | 0.53/0.47 | 0.68/0.32 | 0.82/0.18 | 0.95/0.05 | 1.00/0.00 |

High agreement appears only at `k` large enough that the candidate
admits almost nothing — the agreement is then mostly on "neither",
not on a shared eligible set. At any `k` that admits a usable set,
the candidate-only fraction (unguarded diamonds) is 2–47%. The
reason is structural: below/above counts measure interiority along
the CAUSAL depth (`u`), while the guard's binding conditions are
transverse; the proxy cannot see the direction that matters.

**Consequence for §7:** the chain vectors do NOT recover pure-poset
standing from this candidate. The strictly single-poset claim stays on
item 3a (global relation fraction) unless a different order-invariant
rule is proposed — one that would have to encode transverse
interiority, which no degree count does.

### 4. Ambiguity and cost — a non-issue at probe densities

~7.5M generic pairs across all points and sprinklings: **zero
escalations, zero ambiguous pairs**. Generic decision ~2.0 µs;
a constructed on-cone pair escalates at ~55 µs (~25×). The
escalation path exists for adversarial/boundary inputs, not for
sprinkled ones; `ambiguous_fraction` at these densities is 0 and
§5.1's interval `[D_lower, D_upper]` collapses to a point.

## Caveats, all real

- `n90` columns take §4.4's predicted δ as the signal — this is a
  reconnaissance of *whether P2 can afford to check the prediction*,
  not a check of it.
- `λ_A < 1` at the large-`a` points: per-sprinkling counts are far
  from normal; the n90 arithmetic is aggregate-Poisson and holds only
  summed over hundreds of sprinklings. Order-of-magnitude, not sizing.
- MC error at `edge-a2.4`: ~250 accepted samples for `V_A` (±6%),
  ~85 for `V_dis` (±11%); n90 there is ±25%-ish.
- One seed, 24 sprinklings per point for the fractions; the
  eligibility numbers have ~binomial error, visible in the f_vol vs
  f_elem gaps.
- The `a > 2.4` region is unexplored here; admitting boxes exist
  (guard tests pin 2.6, 2.8) but their δ and λ were not measured.

## What this probe caught in itself

Two defects found and fixed during the run, both recorded because the
mechanism generalizes:

1. **δ was computed as `(wτ)⁴/252` with `τ² = 2ΔuΔv`** — importing
   `Δv` into an effect `Δv` cancels out of (24,000× overstated at
   `roomy`, 14× understated at `a = 1`). Caught by the identity
   `|V_A − V_0| ≤ V_dis`, which the probe now asserts on every run —
   the bound, not a reader, is what noticed.
2. **The MC betweenness test never read the sample's `v`** — it
   tested the column condition `c1 + c2 ≤ Δv`, measuring the
   diamond's shadow × the full `v` depth (3× the true volume in the
   test case). Caught by the flat diamond's closed answer
   `π s²Δv²/6`, pinned in the test suite.
