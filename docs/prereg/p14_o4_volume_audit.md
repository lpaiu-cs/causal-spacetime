# P14 O4 — auditing the S1 predicate's volume response against the O3 certified oracle

Preregistration. Frozen before execution; the campaign runs only after
a separate execution approval, from a clean exact checkout of the
freeze commit.

## 1. What this stage is, and what it is not

O3 produced a certified volume for the frozen Schwarzschild diamond
`(r_in, r_out, dt) = (12, 18, 8.5) M`:

    V_true in [56.212737, 57.348019],   half-width 0.9997%.

O4 confronts that certified value with the **operational instrument
stack** that every S-series result rests on: the sprinkling measure
and the S1 causal predicate (Gauss-Legendre flight time at
`tol = 1e-8`, whose reported `err` is a stopping heuristic, not a
proven bound). The stage is therefore a **three-way consistency
audit** — `sampler + S1 flight-time volume response + oracle` — at
certified precision.

It is deliberately **not** a Poisson-count audit. Quantising the
real-valued `L_S1` into Bernoulli acceptances in order to recover a
Poisson label would discard information already paid for and buy a
property of our own random-number construction rather than a property
of the instrument. A genuine Poisson-count audit (4D sprinkling with
`causal_relation`, comparing an order-interval cardinality to `rho V`)
is a separate, larger stage and is not begun here.

**Nothing in Paper A is upgraded by any outcome of this stage.**
Section 6.7 remains operationally anchored. A concordant result says
the instrument stack agrees with the certified continuum volume within
the frozen margin at one frozen configuration; a discordant result
indicts the stack — sampler, predicate, or certification — and not
causal set theory, and not the S4/S5 verdicts, though it would
trigger an investigation of whether the S1 predicate carries a
boundary bias that those verdicts inherited.

## 2. Estimand

The audited quantity is a **deterministic functional of the S1
implementation**, not a random variable:

    L_S1(X) = [ dt - T1_S1(X) - T2_S1(X) ]_+ ,
    V_S1    = integral over the certified box of L_S1 dmu_3 ,

with `T1_S1 = flight_time(r_in, r_X, psi_X)` and
`T2_S1 = flight_time(r_X, r_out, psi_X)` taken as S1's **deterministic
point output**. The solver's `err` never enters the inference: it is
not a certified enclosure, so an interval built from it would not be a
conservative bound. `err` is used for nothing in this stage except
`causal_relation`'s own undecided rule inside G3, where an undecided
answer aborts the run.

The comparison target is the discrepancy

    Delta = V_S1 - V_true ,   V_true in [V_lo, V_hi] (certified).

`mu_3` is the `r^2 sin(theta)` measure. Because the anchors are
radially aligned, `L` depends only on `(r, theta)`; the `phi` integral
contributes the `2 pi` already inside `B_box`.

## 3. Sampling region and constants

All geometry is derived in ONE path (`experiments/oracle/o4_sizing.py`)
from `certified_flight_time.containment_certificate`, so a boundary
edit cannot desynchronise the sizing from the run.

| Constant | Value | Origin |
| --- | --- | --- |
| sampling box `r` | `[11.3536422001, 18.6949600539]` | L4 certificate, widened OUTWARD to binary64 |
| sampling box `psi` | `<= 0.817912881352` | L4 certificate, rounded up |
| `B_box` | `3358.420056` | `2 pi (1 - cos psi)(r_hi^3 - r_lo^3)/3`, MPFR |
| `B_out` | `3381.100139` | S1 patch shell minus `B_box` |
| `scale = dt * B_box` | `28546.5705` | `V = scale * E[Z]` |
| `L_max` | `1.5599927415085` | `dt - T_cert(P,Q)`, **upward** bound |
| `V_ref` | `56.780378` | `(V_lo + V_hi)/2`, the frozen band scale |

Sampling the certified box **widened** outward is safe: it still
contains the diamond, and the extra shell contributes `L = 0`, costing
a little variance and no bias.

`V_ref` is frozen as the band scale. Using `V_true` as the denominator
instead would move the worst-case floor to `W/V_lo`; the sentence
below is stated against `V_ref` throughout.

## 4. The floor: why `tau` cannot approach 1%

The identified set for `Delta` given the data is

    [ C_lo - V_hi ,  C_hi - V_lo ] ,

whose width is at least `W = V_hi - V_lo` no matter how much is
sampled. At the defect-free alternative `V_S1 = V_true = v`, as the
sample grows the set tends to `[v - V_hi, v - V_lo]`, and `v` may sit
anywhere in the certified interval. Containment in `+/- tau V_ref` for
every admissible `v` therefore requires

    tau >= W / V_ref = 1.999426% .

The oracle's **full** certified width, not its half-width, is the
floor. `tau = 1.5%` is strictly unattainable at any exposure;
`tau = 2.0%` clears the floor by 0.06 pp and is unattainable **within
any frozen budget**. Buying a tighter `tau` requires a tighter oracle,
not more sampling. `tau = 3.0%` is frozen for this round.

## 5. Gates

Familywise budget `0.04 + 0.01 = 0.05`, so the composed G1+G2 sentence
carries **simultaneous coverage >= 95%**.

### G1 — primary: volume equivalence inside the box

Interval: Maurer–Pontil (2009) Theorem 4, empirical Bernstein, applied
to `Z = L_S1/dt in [0,1]` and to `1 - Z` (which share `V_n`), union
bound, `delta = 0.02` per side → two-sided **96%**:

    V_n        = sum_{i<j} (Z_i - Z_j)^2 / (n(n-1))
    half_width = sqrt(2 V_n ln(2/delta)/n) + 7 ln(2/delta)/(3(n-1))
    V_S1 in scale * [ Zbar - half_width , Zbar + half_width ]

Coverage is bought by the theorem, not by simulation. A fixed-seed
coverage simulation exists as a regression diagnostic and is never
called a certification. Preconditions are fail-closed: `Z` outside
`[0,1]` aborts (values are **not** clipped, since clipping biases the
mean and voids the hypothesis), `n >= 2`, every value finite. The
i.i.d. unit is the **sample point**; stream partitioning is an
implementation detail and must not enter the statistic.

Decision, with `band = tau * V_ref = 1.703411`:

| identified set vs `[-band, band]` | G1 |
| --- | --- |
| contained | concordant |
| disjoint | discordant |
| otherwise | inconclusive |

### G2 — co-primary: leakage outside the box

The certificate proves the *true* diamond lies inside the box; it says
nothing about where **S1** thinks the diamond is. Sampling only inside
the box would structurally hide a false-positive branch defect, so the
complement is sampled as a co-primary stratum: the S1 patch shell
minus the box, over `t in [t_p, t_q]`.

Both ends of the decision carry one shared budget, `alpha = 0.005`
each (an earlier draft charged only the upper end):

- upper: Clopper–Pearson on the leaking fraction →
  `V_leak <= dt * B_out * p_upper` (using `L_out <= dt`)
- lower: empirical Bernstein on `E[L_out]` → `V_leak >= B_out * dt * lb`

Budget `0.25% * V_ref = 0.141951`. Concordant if the upper bound is
within budget, discordant if the lower bound exceeds it, otherwise
inconclusive.

### G3 — instrumentation prerequisite (not a scientific verdict)

Does `causal_relation` accept exactly the time window `L_S1` implies?
The independent unit is the **spatial point (cluster)**: the two
stress times at one point share that point's flight times and are not
independent Bernoulli trials. A cluster mismatches if either stress
time disagrees.

Frozen stress distribution: **G1's sampling measure conditioned on
`L_S1 > 0`**, probed at the window midpoint (both relations must hold)
and at the upper edge `+1e-6` (the second relation must fail). The
reported rate is over that distribution, not a generic 4D one.

`100,000` clusters, one-sided 95% Clopper–Pearson. Zero mismatching
clusters bounds the rate at `2.995687e-05`. G3 is **never composed
into the 3.25% sentence**; its failure invalidates the stage.

### Stage verdict

| condition | verdict |
| --- | --- |
| G3 invalid | `INVALID` |
| G1 or G2 discordant | `DISCORDANT` |
| G1 underpowered (frozen `n` not reached) | `INCONCLUSIVE` |
| G1 and G2 both concordant | `CONCORDANT` |
| otherwise | `INCONCLUSIVE` |

A `DISCORDANT` reading survives an underpowered run, because
disjointness is a coverage statement; concordance is not, so it
requires the frozen `n`.

## 6. Sizing — analytic worst-case, no distributional model

`n_G1 = 26,200,000` is certified by the following chain, all
inequalities:

1. `Z_i^2 <= c Z_i` for `Z_i in [0, c]` gives
   `V_n <= c Zbar n/(n-1)` **deterministically**, so the half-width is
   a monotone function of the observed mean and needs no separate
   concentration argument. The bound is nearly lossless:
   `(c - m)/c ~ 0.989`.
2. Monotonicity turns `pass` into `Chat in [c_lo, c_hi]`; using the
   upper bound for `H` shrinks the window, so the resulting
   probability is a **lower** bound on the true pass probability.
3. `Var(Z) <= m(c - m)` for `Z in [0, c]` with mean `m` — the only
   distributional input, and an inequality rather than a model.
4. Bernstein's inequality bounds the two tails.

`c = L_max/dt = 0.183528557825` holds under the defect-free
alternative, which is the alternative the power is computed at.
Coverage never uses it — the interval needs only `Z in [0,1]`.

The worst case over `V_true` is proved on the **continuum** by
interval-arithmetic branch and bound over `[V_lo, V_hi]` (a finite
grid is not a certification). `exp(-x)` is bounded above by the
reciprocal of a truncated Taylor lower bound on `exp(x)`, so the proof
uses only the operations already in the digest-pinned interval module.

| quantity | value |
| --- | --- |
| accept window for `Chat` | `[55.9774017, 57.5787933]` |
| EB half-width at `V_hi` | `0.3367017` (59.3% of the budget) |
| continuum proof | 3 boxes, 1 split |
| failure upper bound | `0.09666241` |
| **power lower bound** | **`0.90334`** |
| worst point | `V_hi` |
| normal-approximation power (diagnostic only) | `0.9849` |

`n_G1` is rounded up from the smallest size reaching 0.90
(`26,012,722`, whose 8e-7 headroom is too tight for a machine-checked
continuum proof) at a cost of +0.72%.

`n_G2 = 1,072,696`; `G3 = 100,000` clusters.

## 7. Budget and caps

Cost convention, pinned: one `flight_time` call = 768 us
(`p14_s1_cost.json`, tol 1e-8); one spatial point = 2 calls = 1.536 ms;
one G3 cluster = 4 `causal_relation` calls.

Expected: `54,945,392` calls, `11.72 h`.
Caps: `80,000,000` calls (**1.46x** the expectation) and `24 h`
(**2.05x**). The two ratios differ and are stated separately. Raising
either during or after the run is forbidden; if a cap binds, the
termination reason is recorded and the verdict degrades per section 5.

## 8. Fail-closed and diagnostics

Aborts with no published result: any `flight_time` exception,
non-finite or negative return; `L_S1 > dt`; any `causal_relation`
undecided at a G3 stress point; `Z` outside `[0,1]`.

Recorded as a **diagnostic outcome, independent of the verdict**:
points with `L_S1 > L_max`, i.e. where S1 claims to beat the certified
optical-triangle minimum. These are never used to clip the estimand —
an estimand that truncates the defect it is looking for would be
worthless.

Smoke runs are capped at `n = 2000`. Larger smoke on the frozen anchors
would pre-observe `V_S1` at a precision comparable to the equivalence
band and collapse the freeze boundary (the O3 lesson).

## 9. Seeds

`g1_audit = 40_000_281` and `g2_leakage = 40_000_291`, allocated fresh
in `probe_seed_ledger.FRESH_PROBE_SCALARS` and moved to
`OBSERVED_PROBE_SCALARS` by the results commit. Separate scalars
rather than children of one, so a stratum rerun can never re-enter
another stratum's stream. The smoke stream `40_000_311` (with its two
reserved successors) is spent from allocation, because the contract
tests observe its output.

**G3 has no seed of its own.** Its unit is the boundary-stress cluster
drawn from `G1 measure | L_S1 > 0`, and the leading `g3_clusters`
accepted points of an i.i.d. G1 stream are i.i.d. draws from exactly
that conditional law — so a separate stream would buy no independence
the Clopper–Pearson bound uses. The v1 draft allocated `g3_wrapper`
and then never read it, which would have written a spent seed into the
artifact as provenance for samples G1 produced; the scalar is
withdrawn **unspent** and returns to the unallocated pool.

**Reservation.** `assert_fresh_scalar` reads the ledger and returns; it
creates no persistent state. The campaign therefore publishes a
write-once `p14_o4_reservation.json` **before its first draw**, and
preflight refuses when it exists. This is what makes "observing any
output spends the seed" enforceable rather than aspirational: a run
that aborts on a fail-closed path leaves no result artifact, and
without the reservation the already-observed streams would read as
fresh on the next attempt. Publication is by `os.link`, so of two
concurrent runs at most one may ever open these streams.

## 10. Freeze identity

Content-addressed manifest `p14_o4_freeze_manifest.json` over the
certified surface, the statistics modules, the sizing certification,
the runner, this document, and the O3 result the audit is measured
against; plus an environment lock (python, gmpy2, MPFR, GMP, numpy).
Verified at entry and exit; `--preflight` additionally requires a
clean tree, the absence of the reservation and of the result artifact,
and fresh campaign seeds, and observes nothing. Publication is atomic
and write-once.

**The manifest cannot certify itself.** A commit made after the
approved freeze that edits a protocol file and re-pins the manifest in
the *same* commit passes every digest check, so digest verification
alone would let a drifted-but-self-consistent tree spend the campaign
seeds; recording the start SHA in the artifact makes that traceable
afterwards but does not prevent it. Both `--preflight` and the
campaign therefore take `--freeze-rev`, the full 40-hex commit named
in the execution approval, and refuse fail-closed unless `HEAD` equals
it exactly. This deliberately replaces the O3 "no-argument official
run" convention: the one argument is the approval itself, and a SHA
the approval does not name is a protocol violation recorded in the
artifact.
