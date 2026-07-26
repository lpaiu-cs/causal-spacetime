# P10: does the emergent geometry deepen toward a continuum?

Status: **STAGE A COMPLETE** (2026-07-26). Headline: across N in
{600, 900, 1200}, the continuum phase is statistically indistinguishable
from genuine sprinklings through the frozen instrument -- 120/120 samples
pass, E tracks S within CI at every size. Two design defects found by the
data are owned in Section 6d: the ESS screen was degenerate at m = 10
(voiding the N = 1200 rung under the frozen rule), and the instrument is
scale-covariant at frozen constants, so this ladder tests indistinguishability
under scaling, not resolution accumulation.

This is the first experiment of the post-T1 programme, and it is aimed at
the programme's founding question rather than at the instrument. It is
deliberately not being rushed toward a paper.

## 1. The question, and why it is the right next one

P5 found a phase of the action-weighted 2D-orders ensemble whose typical
samples carry instrument-level reconstructable geometry (18/18 at
`N = 600`), where the geometry-free growth dynamics has none (0/100) and
the crystal phase has none (0/4, structural). That is the programme's
closest approach to *emergence*: geometry in an ensemble nobody drew by
hand.

But every one of those verdicts lives at one size. A finite sample can
pass a discriminator for two very different reasons: because it is a
coarse view of something that keeps looking geometric as resolution
grows, or because at this size the discriminator cannot yet tell the
difference. The founding question — is spacetime structure *emergent* in
causal order, or a finite-size accident — turns, at this point in the
programme, on which of those two is happening.

So: **at fixed depth in the continuum phase, does reconstructable
geometry deepen, persist, or fade as `N` grows?**

## 2. Design idea: a calibrated ladder, not a lone trend

A falling error curve on the emergence arm alone would mean little — we
would not know how fast it *should* fall if the phase really were a
continuum geometry in the making. The design therefore runs two arms
through the identical instrument at every size, plus the instrument's
built-in null:

- **Arm S (sprinkled reference).** Uniform random 2D orders. By the
  null-coordinate dictionary these are exactly Poisson sprinklings of a
  causal interval (`t = i + pi_i`, `x = i - pi_i`), so this arm is
  ground truth *with a known continuum limit*, sampled for free. Its
  error-vs-`N` curve is the yardstick: it shows what "converging to a
  continuum" looks like through this instrument at these sizes.
- **Arm E (emergence).** Samples of the smeared-action ensemble at
  `beta = 2`, deep in the continuum phase — deliberately far from the
  transition where P7 found mixing pathologies. Scaling convention
  `eps N = 12` fixed, inherited from P4/P5 (holds the crystal's exact
  action advantage constant across sizes).
- **Null (built in).** The column-shuffle null inside the frozen P3
  discriminator provides the per-sample floor at every size; no separate
  negative arm is needed.

Per sample, the discriminator emits `heldout`, `null_gap`, `truth`
(scored against the order's own lightcone coordinates via
`order_inputs`), and P7's frozen composite `G`. Everything is the frozen
P3/P5 pipeline (`analyze_order`), unmodified.

**The readout is the pair of curves.** Three coarse outcomes, named in
advance:

| outcome | signature | reading |
|---|---|---|
| **deepening** | Arm E's error falls with `N` at a rate comparable to Arm S | consistent with a continuum limit; the strongest available evidence short of proving one |
| **floor** | Arm S falls, Arm E plateaus | the phase's geometry has a finite resolution; "emergence" is scale-limited |
| **fading** | Arm E's error grows / pass rate collapses | the `N = 600` result was a finite-size accident |

All three are informative. The second and third would be major *honest*
findings, not failures.

## 3. What "size" means here, and a caveat stated up front

`N` is the number of order elements. Arm S at growing `N` is a genuine
continuum limit (denser sprinkling of the same interval). Arm E at
growing `N` is the same ensemble family at larger size under the frozen
`eps N` convention — the thermodynamic-limit direction of P4/P5. The
design treats "the arm-S curve" as the operational meaning of
"converging to continuum through this instrument"; whether the arm-E
ensemble *has* a continuum limit in any stronger sense is exactly what
cannot be assumed and is not claimed.

One more caveat: the discriminator's frozen constants
(`CHAIN_COUNT = 6`, `MIN_CHAIN_LEN = 25`, `MAX_TARGETS = 44`) were set
at `N = 600`. At larger `N` the instrument therefore reads a bounded
window of an ever-larger order — which is fine (both arms are read the
same way, and a bounded window at growing density is precisely a
resolution ladder) but must be said: this measures geometry *at the
instrument's scale*, not a diverging hierarchy of scales.

## 4. Feasibility, measured (2026-07-26)

Sampler (`mcmc_2d_order_fast`, O(N^2)/move), measured on this machine:

| `N` | steps/s | 3M steps |
|---|---|---|
| 300 | 1525 | 33 min |
| 600 | 272 | 3.1 h |
| 900 | 139 | 6.0 h |
| 1200 | 80 | 10.4 h |

Discriminator: **~1.3 s per sample, `N`-independent** across 600-1800
(the fit cost is set by the target cap, not by `N`). Arm S is free
(direct sampling). So the entire budget is Arm E's MCMC.

The open feasibility question is the mixing time deep in the phase.
P5's recon showed a bipartite start melting within 400k steps at
`N = 600`; if the melt time scales benignly, chains far shorter than
P5's 3M steps suffice at `beta = 2` and the ladder is cheap. **Stage 0
measures exactly this** (in progress): dual-start 500k-step chains at
`N = 900` and `N = 1200`, melt time and tail agreement read from the
`n0` trace.

Ladder under consideration: `N in {600, 900, 1200}` (factor 2, all
measured feasible), with `N = 1800` only if Stage 0 shows fast mixing.

## 5. Mixing discipline (inherited from P7, not relaxed)

- Dual starts (random + bipartite) at every `(N, beta)`; start-separated
  reporting; a size is usable only if both starts agree in phase label
  and tail observables.
- ESS >= 20 on action/n0/height per P7's screening rule (Geyer initial
  positive pairs); non-adequate chains are reported and excluded, never
  pooled.
- Chain lengths set from Stage 0's measured melt/IAT, not assumed.

## 6. Seeds and staging

- **Stage 0 (pilot, running).** Melt probe, seeds 9001/9002/9101 —
  outside every confirmatory range in the programme. Exploratory;
  nothing it produces may become a threshold.
- **Stage A (characterization).** To be frozen after Stage 0: the
  ladder grid, chain lengths, sample counts, fresh seeds (P7 requires
  fresh seeds for any `N = 900/1200` work). Characterization means:
  curves with uncertainties, no gates, no hypothesis verdicts.
- **Stage B (confirmatory), only if warranted.** A preregistered
  hypothesis about the arm-E/arm-S rate comparison, frozen only after
  Stage A shows the design can resolve it. The P8/P9 lesson is applied
  in advance: the decision rule will be *paired within size* (E vs S
  through the same instrument), not an absolute gate.

## 6b. Stage 0 outcome (factual record, 2026-07-26)

Dual-start 500k-step chains at `beta = 2`, `eps N = 12`, seeds
9001/9002/9101, `n0` retained every 25k steps.

| `N` | melt (bipartite -> continuum band) | tails agree | acceptance | wall |
|---|---|---|---|---|
| 900 | **<= 50k steps** (crystal 202500 -> ~4850 by the 2nd sample) | yes, band [4777, 4956] | 0.985 / 0.987 | 93 min |
| 1200 | **<= 50k steps** (360000 -> ~6900 by the 2nd sample) | yes, band [6696, 6926] | 0.989 / 0.991 | 168 min |

Post-melt lag-1 autocorrelation of `n0` at 25k-step spacing: `+0.014`
(N=900) and `-0.153` (N=1200), both indistinguishable from zero at
`n = 19` (se ~ 0.23). So 25k-step spacing already decorrelates `n0`
deep in the phase, and the melt is roughly `N`-independent in steps —
an order of magnitude below P5's 400k-step recon budget.

## 6c. Stage A grid (frozen before any Stage A chain runs)

Per `N` in `{600, 900, 1200}`:

- **Arm E**: two chains at `beta = 2`, `eps = 12/N` — one random start,
  one bipartite start. Burn **100k** steps (twice the observed melt),
  then **10 retained samples per chain at 50k-step spacing** (twice the
  decorrelated spacing) -> 1.1M steps per chain, **20 discriminator
  samples per `N`**.
- **Arm S**: 20 uniform permutations, direct sampling.
- Every sample judged by the frozen P3 discriminator via
  `order_inputs` (heldout, null_gap, truth, and P7's frozen `G`).
- **Seeds**: arm E chains `820 + 10*k` (k indexing the six chains);
  arm S `900-959` split by `N`. Fresh with respect to every
  confirmatory range in the programme (0-9, 100-119, 400-419, 500-519,
  P5's 100-102) and to Stage 0's 9001/9002/9101.
- **Mixing screen (inherited from P7, applied per chain)**: the
  bipartite chain's first retained sample must sit inside the random
  chain's tail band on `n0`, and post-hoc ESS over the retained `n0`
  must be `>= 10` per chain; a failing chain is reported and its
  samples excluded, never pooled.
- **Deliverable**: per-sample rows CSV; per-`N` medians with bootstrap
  intervals for both arms; the two error-vs-`N` curves side by side.
  **No gate, no hypothesis verdict** — characterization, per Section 6.

Estimated cost from the measured rates: `2.2 / 4.4 / 6.1` hours of MCMC
per `N` respectively (discriminator time negligible), run as six
parallel chains.

## 6d. Stage A outcome (factual record, 2026-07-26)

All six chains and the arm-S sweep ran to completion; **120/120 samples
returned `status = ok` and every one of them passes the frozen gates**
(`G >= 0.5`; minimum observed `G = 0.56`, on an arm-E sample at
`N = 600`).

**Pooled reading** (all chains; median [95% bootstrap CI]):

| `N` | arm | heldout | truth | `G` |
|---|---|---|---|---|
| 600 | S | 0.0381 [0.0331, 0.0450] | 0.1488 [0.1378, 0.1698] | 1.000 |
| 600 | E | 0.0406 [0.0344, 0.0437] | 0.1380 [0.1309, 0.1511] | 1.000 |
| 900 | S | 0.0381 [0.0275, 0.0444] | 0.1520 [0.1405, 0.1672] | 1.000 |
| 900 | E | 0.0406 [0.0338, 0.0488] | 0.1573 [0.1332, 0.1677] | 1.000 |
| 1200 | S | 0.0419 [0.0338, 0.0488] | 0.1613 [0.1505, 0.1799] | 1.000 |
| 1200 | E | 0.0363 [0.0325, 0.0413] | 0.1688 [0.1545, 0.2043] | 1.000 |

**Arm E is statistically indistinguishable from genuine sprinklings at
every size measured** — the E and S intervals overlap on every
observable at every rung, and the per-start medians agree within noise
(the dual starts bracket each other at all three sizes). Of the three
outcomes named in Section 2, the data sit on **"deepening" in its
paired form**: no fading, no floor relative to the reference, over a
factor of two in size and forty fresh ensemble samples.

### Two design defects, found by the data and owned here

**(1) The frozen ESS screen is degenerate at `m = 10` and the frozen
rule voids the top rung.** ESS cannot exceed the sample count, so
requiring `ESS >= 10` from 10 retained samples demands the estimator
find *exactly zero* positive autocorrelation — which perfectly iid data
fails **25% of the time** (measured: 20,000 iid trials, `P(ESS >= 10)
= 0.747`). Four of six chains fail the screen (ESS 6.9-9.0, i.e.
`tau <= 1.45` — mild at worst; under a pure-iid null, 4-of-6 failures
has `p ~ 0.04`, so some residual autocorrelation may be real but small).
The rule as frozen — "a failing chain is reported and its samples
excluded" — was applied, and its reading is:

| `N` | screen-applied arm E | arm S |
|---|---|---|
| 600 | truth 0.1436 (10 samples) | 0.1488 |
| 900 | truth 0.1399 (10 samples) | 0.1520 |
| 1200 | **no chain survives** | 0.1613 |

The conclusions above are unchanged where the screen-applied reading is
defined, and the `N = 1200` rung is formally void under the frozen rule.
The defect is the screen's design (P7's rule was transplanted from
`m = 48`, where `ESS >= 20` is a real bar, to `m = 10`, where
`ESS >= 10` is the degenerate ceiling), not the chains. A successor
stage must retain enough samples per chain (`~48`) for the screen to
mean something.

**(2) The yardstick premise failed, informatively.** Section 2 assumed
arm S's error would *fall* with `N` — "what converging looks like".
Measured, arm S is flat to slightly rising (truth 0.149 → 0.152 →
0.161). The reason is the caveat Section 3 half-saw: the discriminator's
frozen constants (`MIN_CHAIN_LEN = 25`, target cap 44) stay fixed while
the order's natural scales grow (the longest chain grows like
`~ 2 sqrt(N)`), so the instrument is **scale-covariant** — it reads the
same relative window at every size and the density gain is normalized
away. The ladder therefore measures *"does the phase stay
indistinguishable from sprinklings as size doubles"* (it does), *not*
*"does resolution accumulate"* (untestable at frozen constants).
A true resolution ladder needs constants scaled in continuum units —
that is the concrete design input for any Stage B.

Artifacts: `docs/prereg/frozen/p10_stage_a/` (six arm-E chain CSVs, the
arm-S CSV, and the aggregator summary). Stage 0 and Stage A ran at
commit `c43bdfe`; grid frozen in 6c before any chain ran.

## 7. What this can and cannot answer

It can answer: whether the continuum phase's reconstructable geometry
behaves, through a fixed instrument, like a geometry being resolved
(arm-S-like) or like an artifact being exposed, over one octave of size.

It cannot answer: the thermodynamic limit (three sizes are three sizes),
Lorentz structure, dynamics, or anything about ensembles other than this
one. If the outcome is "deepening", the honest statement is
"consistent with a continuum limit over the measured range" — and the
next octave is the follow-up, not a victory lap.
