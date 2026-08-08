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

Fattest eligible axis diamond per point. **Two sizings, answering
different questions (review R1.1):**

- `n90 detect` — sprinklings to reject a NO-SHIFT null with the raw
  same-points difference `N_A − N_0` (`Var = ρ V_dis`, exact: shared
  points cancel where the predicates agree), taking the **predicted**
  shift `ρ V_0 δ` as the alternative. This is detectability
  reconnaissance, **not** P2 feasibility.
- `sd_Z` — the per-sprinkling standard deviation of P2's preregistered
  residual `Z = N_A − r N_0`, `r` predicted, from
  `Var(Z) = ρ(V_A + r²V_0 − 2rV_∩)` with all three volumes measured.
  Under the prediction `E[Z] = 0`: nothing is detected, and after `n`
  sprinklings the measured ratio carries a standard error of
  `sd_Z / (ρV_0 √n)`. P2 sizes its own tolerance against that; this
  probe only reports the number.

| point | `a` | δ = R−1 | λ_A | V_dis hits | n90 marginal | n90 detect | sd_Z |
|---|---|---|---|---|---|---|---|
| slice-a0.3 | 0.3 | 3.2e-5 | 0.14 | 0¹ | 7.1e10 | [2.3e6, 2.2e9]¹ | ~0¹ |
| slice-a0.6 | 0.6 | 5.2e-4 | 0.26 | 2¹ | 1.5e8 | [1.8e6, 5.6e6]¹ | 0.055 |
| slice-a1.0 | 1.0 | 4.0e-3 | 0.42 | ~250 | 1.6e6 | 1.2e5 | 0.145 |
| aniso-a1.0 | 1.0 | 4.0e-3 | 5.80 | ~600 | 1.1e5 | 4.8e3 | 0.463 |
| roomy-a0.2 | 0.2 | 6.4e-6 | 12.9 | ~340 | 2.0e10 | 5.4e7 | 0.150 |
| high-a2.0 | 2.0 | 7.3e-2 | 0.07 | 5¹ | 2.7e4 | [3.2e3, 6.8e3]¹ | 0.090 |
| edge-a2.4 | 2.4 | 1.8e-1 | 0.38 | ~85 | 834 | **463** | 0.385 |

`λ_A` is the expected count in the diamond at `N = 300`; `V_dis` is
the volume where the two arms' predicates disagree; `V_∩` their
overlap — all from shared-sample MC, and the partition
`V_A + V_0 = V_dis + 2V_∩` holds exactly per sample (test-pinned).

¹ **A low MC count is a bound, never a value (review R2.1, R3.1):**
three rows have too few disagreeing samples for a point estimate. The
`n90 detect` there is a bracket, from the analytic floor
`|V_A − V_0| = δV_0` (best case) to the **exact one-sided 95% Poisson
upper limit** for the observed count times the quantum (worst case) —
`2.996 q` at 0 hits (the rule of three), `6.30 q` at 2, `10.5 q` at 5.
The earlier `max(point, 3q)` under-stated this for every positive
count and called a one-hit row resolved; the estimator now stores the
raw count and computes the real limit (`_poisson_upper_95`,
test-pinned against the standard values). All three unresolved rows
are dead either way. `slice-a0.3`'s `sd_Z` is likewise a floor, not
an estimate. The resolved rows carry hundreds of hits and their
point estimates stand.
Worked example at `edge-a2.4`: `sd_Z = 0.385` per sprinkling and
`ρV_0 ≈ 0.32`, so 100 sprinklings verify the ratio to
`±0.12` (1σ) against a predicted shift of `0.18` — P2's tolerance
choice decides whether that is enough or `n` must grow.

**Reading.** The product (effect)² × (count) is maximized at the
largest `a` the guard admits, despite the tiny boxes: δ grows as `a⁴`
and beyond, while the count only needs to be O(1) per sprinkling.
`edge-a2.4` needs ~4e2 sprinklings to detect the shift (~8e2
marginal); `roomy`
and the `slice` points are dead at any affordable n. Density is a
further free lever this table does not use: P2-style counting inside
an external-anchor diamond is O(N) per anchor, not O(N²), so raising
`N` at the small boxes cuts n linearly at modest cost.

**The split that matters downstream:**

- **External-anchor Class C (P2's volume check, anchors fixed and not
  counted, per the approval condition):** what this probe establishes
  is that the geometry is ADMISSIBLE at large `a` and the predicted
  shift is DETECTABLE there (~4e2 sprinklings at `edge-a2.4`), and
  that P2's validation precision is quantified by `sd_Z`
  (±0.12 on the ratio per 100 sprinklings at `edge-a2.4`).
  **P2 feasibility itself is conditional on the tolerance P2 has not
  yet chosen** — this probe cannot pre-decide it (review R2.2).
- **Sprinkled-pair Class C (§5 items 2, 3b — chain vectors, MM on
  sprinkled intervals): starved exactly where the effect lives.**
  At `a ≥ 2` the run saw 0.0 eligible pairs per sprinkling; the
  viable middle is `aniso-a1.0`-like (26 pairs/sprinkling, δ = 4e-3).
  Whether that middle carries a usable statistic is §8 P3's question,
  with these numbers as its input.

### 3. Order-invariant candidate (§4.6.2) — does not reproduce the guard

Candidate: interior iff ≥ `k` elements below AND above. §8 P1 asks
for the PAIRS each rule admits that the other rejects (review R1.3) —
both rules factorize over endpoints, so pair admissions are derived
combinatorially from the element categories (test-pinned against
materialized pair masks). Pair-level agreement / candidate-only pair
fraction of `C(N,2)` (the dangerous direction — admitted pairs whose
diamonds carry no containment certificate):

Both exclusive directions are reported (review R3.2), as
agree / candidate-only / guard-only, all as fractions of `C(N,2)`.
Candidate-only is the dangerous direction (admitted pairs with no
containment certificate); guard-only is lost coverage — and reporting
it is what tells an "ideal candidate" (`1.00/0/0`) apart from an
"empty candidate that rejects everything the guard admits":

| point | k=2 | k=4 | k=8 | k=16 |
|---|---|---|---|---|
| slice-a1.0 | 0.993/0.006/0.001 | 0.999/0.000/0.001 | 0.999/0.000/0.001 | 0.999/0.000/0.001 |
| aniso-a1.0 | 0.901/0.098/0.001 | 0.970/0.029/0.001 | 0.996/0.004/0.001 | 0.999/0.000/0.001 |
| edge-a2.4 | 0.896/0.104/0.000 | 0.968/0.032/0.000 | 0.997/0.003/0.000 | 1.000/0.000/0.000 |

**Every integer `k` from 1 to 20 was swept (review R3.3)** — the mask
is monotone, so an intermediate threshold cannot hide between sampled
values. For each point, the first `k` that admits **zero**
candidate-only (unsafe) pairs, and what the candidate still admits
there:

| point | first safe `k*` | guard-only at `k*` | guarded pairs the candidate keeps |
|---|---|---|---|
| slice-a1.0 | 8 | 0.001 | **0.0** |
| aniso-a1.0 | 19 | 0.001 | **0.0** |
| roomy-a0.2 | 6 | 0.094 | **0.0** |
| edge-a2.4 | 17 | 0.000 | **0.0** |

**This is the finding, and it is sharper than "does not reproduce".**
At every operating point, the smallest `k` that removes the unsafe
admissions also removes essentially all the useful ones: the candidate
keeps ~0 guard-admitted pairs at `k*`. There is no intermediate
threshold that is both safe and non-empty. The reason is structural:
below/above counts measure interiority along the CAUSAL depth (`u`),
while the guard binds transversally; a `k` high enough to reject the
transversally-shallow pairs rejects nearly all of them. (Element-level,
same run: agreement 0.53–0.98, candidate-only up to 47%.)

**Consequence for §7:** the chain vectors do NOT recover pure-poset
standing from this candidate. The strictly single-poset claim stays on
item 3a (global relation fraction) unless a different order-invariant
rule is proposed — one that would have to encode transverse
interiority, which no degree count does.

### 4. Ambiguity and cost — a non-issue at probe densities

§5.1 defines a pair as ambiguous when EITHER arm is undecided, so
both arms are censused and the union taken (review R1.2): ~7.5M
generic pairs across all points and sprinklings, **zero escalations
in either arm, zero ambiguous pairs in the union**. Generic decision
~2.0 µs; a constructed on-cone pair escalates at ~55 µs (~25–30×). The
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
  ~85 for `V_dis` (±11%); n90 there is ±25%-ish. Rows with fewer
  hits carry their own flags (footnotes 1–2).
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

Review round R1 then corrected three more, folded into the tables
above: the paired sizing presented a raw-difference detectability
number as P2 feasibility (different estimands — P2's `Z` has
`E[Z] = 0` under the prediction and its variance uses `V_∩`); the
ambiguity census read only the curved arm where §5.1 is defined over
both; and the candidate comparison was reported at element level
where §8 P1 asks for pair admissions, which compound.

Round R2 corrected two follow-ons: a zero-hit MC estimate of `V_dis`
flowed into the detection sd as an exact zero, manufacturing a
zero-noise instrument (footnote 1 — the estimator now returns its
quantum and the sizing consumes the bound); and the downstream
conclusion "external-anchor Class C: FEASIBLE" restated the estimand
conflation R1.1 had just removed — narrowed to what the probe
establishes: admissible, detectable, precision quantified, tolerance
P2's own.

Round R3 corrected three more: the `V_dis` upper bound was
`max(point, 3q)`, which is the 95% limit only at zero hits — replaced
by the exact one-sided Poisson limit for the observed count (`6.30 q`
at 2 hits, not `3 q`), which also fixed `disagree_resolved` calling a
one-hit row resolved; the candidate table omitted the guard-only
direction §8 P1 asks for, so `1.00/0.000` could not be told from an
empty candidate — both directions reported now; and the `k` sweep
skipped `9..15`, so an intermediate safe-and-non-empty threshold had
not been ruled out — every integer `k ≤ 20` is swept now, and the
`k*` table shows there is none.
