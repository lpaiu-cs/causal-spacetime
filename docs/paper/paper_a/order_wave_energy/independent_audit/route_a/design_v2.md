# Route A, second design: wall participation from an order-preserving chain

2026-11-24. Amendment to [design.md](design.md), **written after its pilot
was seen**. Section 1 says what the first design did and why its main stage
was not run. Section 2 lists every number seen before this freeze. Nothing
here is a result of the second design. The question, the definitions and
the two hypotheses are those of design.md Section 1.

## 1. What happened to the first design

- **Control stage: passed.** Stamp `741ddfc`. In 1+1 dimensions the searches
  found 750 of 750 known realisable flips from the sampled configuration and
  certified none outside the rule (Wilson 95% lower end 0.9949 against the
  floor 0.99). [evidence/control](evidence/control/summary.json).
- **Pilot stage: ran, and showed the instrument cannot carry the main
  stage.** Stamp `228cc84`, N = 96, 8 orders, 64 tested events, 3125
  candidates. [evidence/pilot](evidence/pilot/summary.json).
  - The completeness gate would fail. Of 2187 candidates certified by the
    chain or the second start, the searches from the sampled configuration
    had found 1398: `0.639`, Wilson 95% `[0.619, 0.659]`, against the floor
    `0.90`.
  - The reason is the search budget, not a stall: 1646 joint searches
    failed after a mean of 34.2 of the 40 allowed programmes, with the
    target margin still falling (median final over initial 0.455).
  - Cost was 440 s per order against a planned 60 to 120 s. At that rate
    N = 512 alone is above 30 CPU hours.
- **Main stage: not run.** Running it as frozen would have spent that time
  to return "unresolved: instrument gate failed". The pilot was declared
  for sizing only; it is not read as a verdict on U1 here or anywhere else.
  `evidence/pilot/sizing.json` (12 orders per size) is retained and unused.

The frozen design.md and its instrument are left as they are.

## 2. What was seen before this freeze

Everything in design.md Section 9, and from the pilot at N = 96:

| per tested event | mean |
|---|---|
| candidates `F` | 48.8 |
| walls at the sampled configuration `S` | 15.3 |
| certified from the sampled configuration `R` | 23.1 |
| certified by any search `U` | 35.4 |

- Of the walls at one chain state, the fraction still walls 250 to 2000
  sweeps later was 0.64 to 0.59, flat in the lag. The participation
  estimate of design.md from those overlaps was 25.8.
- Quartiles of `N |m|` over all candidates: 12.5, 28.3, 50.3; over
  certified ones: 9.0, 19.5, 38.4.
- The chain's root-mean-square displacement was 0.32 after 250 sweeps and
  0.38 after 2000, in a diamond of proper height 1.48. The order was checked
  exact at every checkpoint.
- Development runs of the second instrument in 1+1 dimensions gave
  participation 3.77 to 3.99 on six orders (known value 4, Section 4), and
  no hit outside the candidates. In 3+1 only its time per sweep was printed
  (0.020, 0.052, 0.149 s at N = 512, 1024, 2048).

**What the pilot changes in the analysis of design.md Section 2.** The
scaling argument stands; its constants were not written down, and they
matter at these sizes.

- Link scale. With `rho V = N r^2 s^2 / TOTAL^2` for a pair at spatial
  separation r and null margin s, links have `|s| <~ 1.09 / (r sqrt(N))`,
  about 2 to 5 over `sqrt(N)`.
- Counting scale. The pilot has about one candidate per event per unit of
  `N |m|`, so rows with margin below w number about `N^2 w / 2`, and the
  counting threshold `8N` of Section 2.4 is reached at `w* = 16 / N`.

The ratio of the two is about `(0.12 to 0.3) sqrt(N)`: of order one at
N = 96, where "most candidates are realisable" (`U/F = 0.73`) is what both
hypotheses predict. It reaches 4 to 10 only near N = 1000 to 2000. The
first design's grid stopped at 512. A probe of U1 has to reach further,
which an instrument built on linear programmes over all pairs cannot afford
and the one below can.

## 3. Instrument

[hit_probe.py](hit_probe.py). The chain of the first instrument, unchanged
in law: one sweep moves each event once along a uniformly random line; the
proposal is uniform on the segment of that line on which the order is
exactly C; acceptance is the density ratio. The start is the sampled
configuration, an exact draw from the conditional law given C, so every
state is.

Each end of the segment lies on the cone of one other event, and just
beyond it exactly that one relation has changed. So the partner at a segment
end is a wall of the moved event at the current state, found by exact
geometry. Segment ends on a tip's cone are exits from the diamond; they are
counted separately and excluded.

For a tracked event j let `p_j(i)` be the probability, under the
conditional law given C and a uniform direction, that the partner is i. Its
participation is `1 / sum_i p_j(i)^2`. The chain is cut into 8 blocks. With
`f_b` the partner frequencies of block b,

`Q_j(L) = mean over block pairs at lag L of sum_i f_b(i) f_b'(i)`

is unbiased for `sum_i p_j(i)^2` when the two blocks are independent,
because counting noise in different blocks is independent. If the blocks
are still correlated `Q_j` is too large and the participation is
underestimated. The size-level estimate is

`K(N) = 1 / mean over tracked events and orders of Q_j(far lags 4 to 7)`.

**Relation to the wall density of design.md.** `p_j(i)` weights facet i of
the single-move cell by its area times a chord length; `rho_e` weights it
by area over volume. Within one cell these differ by a shape factor of
order one, so the two participations have the same scaling in N. They are
not equal, and no statement here depends on their ratio.

**Internal exactness check.** A wall is a single-relation flip that keeps
the order a partial order, so every partner must be a candidate of
design.md (a link or a critical pair in either orientation). The count of
hits outside the candidates is recorded and must be zero.

## 4. Campaign

CPU only, clean-tree preflight, commit stamp, no overwriting. Seeds
`SeedSequence([20261125, dimension, N, r])`. Per order: 64 tracked events
drawn uniformly, `16 N` sweeps in 8 blocks of `2 N`, exact order comparison
after every block.

| stage | geometry | N | orders per N |
|---|---|---|---|
| control | flat 1+1 diamond | 64, 128, 256 | 8 |
| main_a | background q = 1 | 128, 256, 512, 1024 | 12 |
| main_b | background q = 1 | 2048 | 12 |

The chain length grows with N on purpose. If relative positions were free
at the link scale `3 / sqrt(N)` and single moves have width `c / N`,
diffusion to equilibrium takes of order N sweeps; `16 N` leaves room for
that in the scenario where it is needed. main_a and main_b are one campaign
split only so that the long size can run on its own.

**Control, known answer.** In 1+1 dimensions with a flat density the walls
of an event are its four rank neighbours, and by the exchangeability of the
spacings and the symmetry of the direction each has `p = 1/4`. So the
participation of an event with four real neighbours is exactly 4.

## 5. Gates

**Instrument gates.**

1. *Control.* The 95% two-sided interval of the participation of interior
   events, pooled over the control sizes (bootstrap over orders), must lie
   inside `4 +/- 0.3`, and no hit may fall outside the candidates. The
   tolerance was chosen after the development runs of Section 2, which sat
   slightly below 4; the likely cause is the slow collective modes that 1+1
   dimensions has and 3+1 should not, and the gate is wide enough to pass
   them and narrow enough to fail a wrong partner rule (3 or 5).
2. *Exactness in 3+1.* Zero hits outside the candidates.
3. *The grid can discriminate.* The 95% interval of the log-log slope of
   the mean of `F_j` over the gate sizes lies above `0.25`.

**Decision gate.** Gate sizes N = 256, 512, 1024, 2048: N = 128 is run and
reported but left out of the slope, because Section 2 puts it where the two
scales coincide. With `kappa_K` and `kappa_F` the least-squares slopes of
`log K(N)` and `log mean F_j` on `log N`,

`D = kappa_K - kappa_F / 2`,

95% two-sided percentile interval from a cluster bootstrap over orders
within each size, 10000 draws, seed 20261125.

- interval entirely above 0: **count-limited participation**. This reading
  does not need the mixing gate, since an unmixed chain can only lower K.
- interval entirely below 0 **and** mixing shown at every gate size:
  **cell-limited participation**.
- interval entirely below 0 without mixing shown: **unresolved: flat but
  not shown mixed**.
- otherwise: **undecided**.

*Mixing.* At each gate size the ratio of the mean of `Q_j` at lags 1 to 2
to its mean at lags 4 to 7 must have a 95% interval upper end of at most
`1.10`. A flat ratio cannot exclude a mode slower than the whole chain;
Section 6 says what that leaves open.

**Sizing for 90% power.** As in design.md Section 6: drift allowance
`0.05`, design value `kappa_F = 0.40`, requirement `SE(D) <= 0.0463`, which
on four sizes a factor 2 apart is a standard error of `log K` of at most
`0.0717` per size. Planning value for the between-order relative standard
deviation: `0.174`, the pilot's value for 8-event means of U; means over 64
tracked events should vary less. Twelve orders give `0.174 / sqrt(12) =
0.050`, inside the requirement with a margin of 1.4. After the run the
realised `SE(D)` is compared with `(kappa_F_hat / 2 - 0.05) / 3.2415`, and
a larger value makes the verdict **unresolved: underpowered**.

If main_b is not complete the analysis labels its output provisional and no
verdict is recorded.

## 6. Limits

- The verdicts are about the scaling of wall participation between N = 256
  and 2048 on one background and one diamond. The scale ratio of Section 2
  is 2 to 5 at the low end and 5 to 14 at the top. A slope over that range
  is evidence about the asymptotic alternative, not a proof of it.
- Single-event moves. A collective mode slower than `16 N` sweeps would
  leave `Q` flat in the lag and K too low. The displacement of the chain at
  every block end is recorded so that a still-growing displacement can be
  seen. Exact symmetries of the order (the conformal maps preserving the
  diamond) move events without changing any wall, so they do not bias K
  whether or not they have relaxed.
- `p_j` is a hit distribution, not the wall density (Section 3).
- No wave, no likelihood, no total-variation statement. What a verdict
  triggers is in design.md Section 7, with "participation" in place of
  "realisability" and `K <= U` no longer needed.
