# P10: does the emergent geometry deepen toward a continuum?

Status: **CLOSED** (2026-07-27, at Stage B'0 — see Section 9.5). The
programme's founding continuum-limit question is UNDECIDABLE BY THIS
INSTRUMENT FAMILY: through the frozen discriminator, truth accuracy in
continuum units does not improve — and measurably degrades — as the
order grows, under every combination tried of instrument scaling
(frozen / continuum-scaled) and question restriction (full mix /
margin-cutoff / mix-matched bins — see 9.5-9.6). Three gates, three
failures; every per-size result in the programme stands unchanged.
What a continuum-limit test now requires is a reconstruction whose
continuum accuracy grows with density — a different instrument.
The Stage A record below is preserved as written.

Dates in this record are local dates (UTC+9, Asia/Seoul), matching the
repository's commit timestamps, which carry their +09:00 offset — e.g.
a section dated 2026-07-27 may sit in a commit reading
2026-07-26T18:14Z when normalized to UTC. The chronology is verifiable
from commit order; no date in this document is prospective.

Stage A record: **STAGE A COMPLETE** (2026-07-26; readings revised same day when
review found the ESS re-implementation diverged from the inherited P7
diagnostic — see correction notes in 6d). Headline, primary reading
(frozen mixing screen, computed with the shared P7 diagnostic): one
chain per rung survives, and **no E − S difference is detected at any
size** — median-truth differences −0.005 [−0.030, +0.024] at N = 600,
+0.006 [−0.011, +0.026] at N = 900, and +0.036 [−0.012, +0.052] at
N = 1200 — with 120/120 samples passing the frozen gates. The N = 1200
interval covers zero but leans positive with a point estimate near the
detection sensitivity, so the failure to detect is least informative
exactly at the top rung; that is stated rather than smoothed. "No
detectable difference" is a failure to detect at this sample size, NOT
demonstrated equivalence: Stage A prespecified no equivalence margin,
deliberately, and any equivalence claim is deferred to a Stage B with
one. Two design defects found by the data are owned in Section 6d: the
ESS screen is degenerate at m = 10 (iid data fails it 40% of the time
under the inherited implementation), and the instrument is
scale-covariant at frozen constants, so this ladder tests behaviour
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
  then **10 retained samples per chain at 50k-step spacing** ->
  **600k steps per chain** *(corrected: this line originally said
  "1.1M", a wrong sum of its own stated inputs — `100k + 10 x 50k =
  600k`. The runner implements the stated inputs exactly; the archived
  CSVs record retained samples at steps 100k through 550k. The inputs,
  not the mis-summed total, were always the frozen content)*, **20
  discriminator samples per `N`**.
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

## 6d. Stage A outcome (factual record, 2026-07-26; corrected same day
after review — see the correction note at the end of this section)

All six chains and the arm-S sweep ran to completion; **120/120 samples
returned `status = ok` and every one of them passes the frozen gates**
(`G >= 0.5`). Minimum observed `G`: **0.5125, on an arm-S sample at
`N = 900`** (the arm-E minimum is 0.5625, at `N = 600`).

**Primary reading — the frozen mixing screen applied** (failing chains
excluded per prereg 6c; ESS from the shared P7 diagnostic; the E column
at each rung is the surviving chain):

| `N` | surviving chain | arm E truth | arm S truth | E − S diff [95% CI] |
|---|---|---|---|---|
| 600 | random | 0.1436 (10) | 0.1488 | −0.005 [−0.030, +0.024] |
| 900 | bipartite | 0.1581 (10) | 0.1520 | +0.006 [−0.011, +0.026] |
| 1200 | random | 0.1971 (10) | 0.1613 | +0.036 [−0.012, +0.052] |

The `N = 1200` row deserves its own sentence: the interval covers zero,
but the point estimate (+0.036) sits near the detection sensitivity and
the surviving chain is the higher of the two (the excluded bipartite
chain's median is 0.1596, close to arm S). So at the top rung the
primary reading neither detects a difference nor gives much comfort —
the pooled secondary reading there (+0.008 [−0.019, +0.046]) is
friendlier, and the honest statement is that the top rung is where a
longer-chain successor stage matters most.

**Secondary reading — pooled, ignoring the screen** (kept only because
the screen itself is documented below as degenerate; median [95% CI]):

| `N` | arm | heldout | truth | E − S truth diff |
|---|---|---|---|---|
| 600 | S | 0.0381 [0.0331, 0.0450] | 0.1488 [0.1378, 0.1698] | |
| 600 | E | 0.0406 [0.0344, 0.0437] | 0.1380 [0.1309, 0.1511] | −0.011 [−0.032, +0.009] |
| 900 | S | 0.0381 [0.0275, 0.0444] | 0.1520 [0.1405, 0.1672] | |
| 900 | E | 0.0406 [0.0338, 0.0488] | 0.1573 [0.1332, 0.1677] | +0.005 [−0.022, +0.025] |
| 1200 | S | 0.0419 [0.0338, 0.0488] | 0.1613 [0.1505, 0.1799] | |
| 1200 | E | 0.0363 [0.0325, 0.0413] | 0.1688 [0.1545, 0.2043] | +0.008 [−0.019, +0.046] |

**What may be claimed, stated with its statistical character.** At every
rung where the primary reading is defined, and at every rung of the
secondary reading, the E − S difference CI covers zero: **no difference
was detected** between the emergent phase and genuine sprinklings, at
sample sizes able to detect a median-truth shift of roughly 0.03-0.05.
That is a *failure to detect*, not equivalence — Stage A prespecified no
equivalence margin (deliberately: a characterization stage has nothing
to tune, and a margin chosen now would be chosen after seeing the
spread). An equivalence claim needs a Stage B with a preregistered
margin, exactly as P9 treated its no-improvement leg.

**Against the frozen outcome table of Section 2: unclassifiable.** All
three preregistered signatures share a premise — that arm S's error
falls with `N` — and that premise failed (S runs 0.149 → 0.152 →
0.161). Neither arm falls, so "deepening"'s signature was not met;
arm S does not fall, so "floor"'s was not; arm E neither grows beyond
the reference nor loses passes, so "fading"'s was not. An earlier
version of this record called the result "the paired form of
deepening", which was a post-hoc fourth category wearing the friendliest
frozen label; that classification is withdrawn. What stands is the
descriptive result above — no E − S difference detected at any size, at
the stated sensitivity, on **60 fresh ensemble samples in the pooled
reading (30 surviving the frozen screen in the primary)**, with the dual
starts bracketing each other — and the design finding that the frozen
outcome table was built on a yardstick premise the instrument does not
satisfy at frozen constants. A successor stage with constants scaled in
continuum units is where the Section 2 table becomes decidable.

### Two design defects, found by the data and owned here

**(1) The frozen ESS screen is degenerate at `m = 10`.** ESS cannot
exceed the sample count, so requiring `ESS >= 10` from 10 retained
samples demands the estimator find *exactly zero* positive
autocorrelation — which perfectly iid data fails **40.4% of the time**
under the inherited P7 implementation (measured: 20,000 iid trials; an
earlier figure of 25% in this note was computed with a mis-paired
re-implementation of the estimator, since corrected — note (ix)).
Three of six chains fail the screen (ESS 9.3-9.9, `tau <= 1.08` —
negligible; under a pure-iid null, 3-of-6 failures at a 40.4% rate has
`p ~ 0.46`, i.e. **the observed failures are entirely consistent with
the degenerate screen operating on well-mixed chains**). The rule as
frozen — "a failing chain is reported and its samples excluded" — is
the **primary** reading above, which keeps exactly one chain per rung;
the pooled table is retained as secondary because the screen is
defective in the way just quantified. The defect is the screen's design
(P7's rule was transplanted from `m = 48`, where `ESS >= 20` is a real
bar, to `m = 10`, where `ESS >= 10` is the degenerate ceiling), not the
chains. A successor stage must retain enough samples per chain (`~48`)
for the screen to mean something.

**(2) The yardstick premise failed, informatively.** Section 2 assumed
arm S's error would *fall* with `N` — "what converging looks like".
Measured, arm S is flat to slightly rising (truth 0.149 → 0.152 →
0.161). The reason is the caveat Section 3 half-saw: the discriminator's
frozen constants (`MIN_CHAIN_LEN = 25`, target cap 44) stay fixed while
the order's natural scales grow (the longest chain grows like
`~ 2 sqrt(N)`), so the instrument is **scale-covariant** — it reads the
same relative window at every size and the density gain is normalized
away. The ladder therefore measures *"is a difference from sprinklings
detected as size doubles"* (none was, at the sensitivity stated above —
which is a failure to detect, not equivalence), *not* *"does resolution
accumulate"* (untestable at frozen constants).
A true resolution ladder needs constants scaled in continuum units —
that is the concrete design input for any Stage B.

Artifacts: `docs/prereg/frozen/p10_stage_a/` (six arm-E chain CSVs, the
arm-S CSV, and the aggregator summary). Stage 0 and Stage A ran at
commit `c43bdfe`; grid frozen in 6c before any chain ran.

**Correction note (same day, after automated review of PR #23).** Four
review findings were applied to this record; the measurement rows are
untouched. (i) The first version of this section claimed "statistically
indistinguishable" from overlapping marginal intervals — an equivalence
claim that overlap cannot support. Rewritten as "no difference
detected", with explicit E − S difference CIs and the failure-to-detect
character stated. (ii) The first aggregator pooled all chains despite
the frozen screen; the screen-applied reading is now computed by the
code and is primary, with the `N = 1200` rung null. (iii) The
aggregator's bootstrap seeds came from Python's per-process salted
`hash()`, making the archived intervals unreproducible; seeds now come
from CRC32 of fixed labels, verified identical across processes, and
the frozen summary was regenerated accordingly (the earlier archived
intervals differ in the third decimal at most). (iv) The minimum-`G`
record misattributed the global minimum to arm E at `N = 600`
(`0.5625`); the global minimum is `0.5125`, on an arm-S sample at
`N = 900`. The gate-pass claim (120/120 above 0.5) was and remains
true.

**Second correction note (same day, second review round).** (v) The
6c chain-length line said "1.1M steps per chain" — a wrong sum of its
own stated inputs (`100k + 10 x 50k = 600k`); the runner and the
archived step columns implement the stated inputs, and 6c now carries
the arithmetic correction inline. (vi) The aggregator's `screen_pass`
used the ESS criterion alone, so an unmelted chain with clean ESS
would have entered the primary reading; it now requires melt AND ESS
(no archived verdict changes — all six chains melted). (vii) Undefined
statistics were serialized as bare `NaN` tokens, which are outside
strict JSON; they are now `null`, and the frozen summary was
regenerated and verified against a strict parser. (viii) One sentence
in defect note (2) still concluded "indistinguishable … (it does)";
it now uses the same failure-to-detect wording as the rest of the
record.

**Third correction note (same day, third review round; this one changes
readings, not only wording).** (ix) The aggregator's ESS was a
re-implementation of Geyer pairing with the indexing off by one — it
paired `rho[1]+rho[2], rho[3]+rho[4]` where the inherited P7 diagnostic
(`integrated_autocorrelation`) pairs `rho[0]+rho[1], rho[2]+rho[3]`
with the `−1` adjustment. The screen *claims* to inherit P7's
convention, so P7's code is the convention; the aggregator now calls it
directly. On the archived rows this flips three of six chain verdicts:
the surviving chains are now 600-random, 900-bipartite, 1200-random —
in particular **the `N = 1200` rung, called void in the first two
versions of this record, is defined** with the random chain's 10
samples, and its E − S truth difference is +0.036 [−0.012, +0.052].
Every table and the status header were recomputed; the earlier
correction notes' statements that the 1200 rung was void are superseded
by this note and left in place as the record of what was believed when.
The iid failure rate of the degenerate screen is 40.4% under the
inherited implementation (the earlier 25% was measured with the
mis-paired one). (x) The melt criterion silently widened the frozen
"inside the random tail band" by 10% each way; it is now the strict
band. No archived melt verdict changes — all three bipartite first
samples lie strictly inside their bands.

**Fourth correction note (same day, fourth review round).** (xi) Two
earlier versions of this record classified the outcome as "the paired
form of deepening". Under the frozen Section 2 table that was wrong:
all three preregistered signatures presuppose a falling arm-S curve,
which did not occur, so **the ladder is unclassifiable under its own
frozen outcome table** and the descriptive difference result is what
stands. The friendly label was a post-hoc category and is withdrawn —
the same discipline P8 applied when it declined to invent a passing
rule after the fact. (xii) The archived summary's `code_version`
stamped a commit (`1ae9f61`) whose aggregator still carried the
mis-paired ESS, so checking out the recorded version could not
reproduce the recorded reading. The summary is regenerated at a commit
that contains the corrected aggregator, and reproduction is: check out
the stamped commit, run `--aggregate`, compare.

## 7. What this can and cannot answer

It can answer: whether the continuum phase's reconstructable geometry
behaves, through a fixed instrument, like a geometry being resolved
(arm-S-like) or like an artifact being exposed, over one octave of size.

It cannot answer: the thermodynamic limit (three sizes are three sizes),
Lorentz structure, dynamics, or anything about ensembles other than this
one. If the outcome is "deepening", the honest statement is
"consistent with a continuum limit over the measured range" — and the
next octave is the follow-up, not a victory lap.

## 8. Stage B: the resolution ladder (designed 2026-07-26; **B0 ran
2026-07-27 and its gate FAILED — Stage B stopped at 8.5; no B1 chain
ever ran**)

Stage A's findings dictate this design item by item, and every review
lesson from PR #23's five rounds is pre-applied rather than rediscovered.

### 8.1 The scaled instrument family

One pipeline definition — `analyze_order`, now parameterized with
keyword defaults that ARE the frozen constants (the defaults are pinned
byte-identical against archived Stage A rows) — evaluated at operating
points scaled in continuum units:

| constant | law | 600 | 900 | 1200 |
|---|---|---|---|---|
| `min_chain_len` | `25 sqrt(N/600)` | 25 | 31 | 35 |
| `max_targets` | `44 N/600` | 44 | 66 | 88 |
| `min_targets` | `20 N/600` | 20 | 30 | 40 |
| `train_c` | `3000 N/600` | 3000 | 4500 | 6000 |
| `heldout_c` | `800 N/600` | 800 | 1200 | 1600 |
| `chain_count` | fixed | 6 | 6 | 6 |

Rationale: the longest chain grows like `~ 2 sqrt(N)`, so a fixed
`min_chain_len` is a shrinking continuum window — `sqrt(N)` scaling
holds it fixed. Target and constraint counts scale together because P8
measured what happens when they do not: a fixed budget over more targets
dilutes the fit. `chain_count` stays 6 — observer count is a design
choice, not a resolution property. **The `N = 600` anchor is exactly
the frozen instrument**, so Stage B connects continuously to everything
the programme has measured. P7's `G` margins are not scaled: the score
definition is frozen.

### 8.2 Stage B0 — the yardstick pilot, and why it gates everything

Stage A *assumed* the sprinkled arm would fall with `N` and it did not.
B0 tests exactly that premise before anything else is spent: arm S
(direct sampling, cheap) through the scaled instrument at all three
sizes, 20 samples each, seeds 1100-1159 (fresh).

**Gate for proceeding** (this is a design prerequisite, not a science
hypothesis): the `N = 1200` minus `N = 600` median-truth difference must
have its entire 95% bootstrap CI below zero. If it does not, the scaled
instrument has not produced a resolution ladder either, Stage B stops,
and that is the reportable finding. No E chain runs and no hypothesis
is frozen until B0 passes.

### 8.3 Stage B1 — frozen only after B0, stated now so the shape is on record

- Six chains as in Stage A (dual starts, `beta = 2`, `eps N = 12`),
  seeds **30000-30050** *(corrected: the original table said 1000-1050
  and was not fresh -- P6-B's uniform reference orders are exactly
  `default_rng(1000..1019).permutation(600)`, so the (600, random)
  start would have literally reused one; see the 8.5 correction
  note)*, burn 100k, **48 retained samples per chain
  at 25k spacing** (Stage 0 measured 25k decorrelated) -> 1.3M steps
  per chain. `m = 48` restores the inherited screen to P7's actual bar
  (`ESS >= 20` at `m = 48`), ending the Stage A degeneracy.
- Mixing screen: the corrected Stage A aggregator's rules verbatim —
  shared P7 ESS diagnostic, strict random-tail band, both criteria
  combined per chain, failing chains reported and excluded.
- Uncertainty under dependence *(design amendment 2026-07-27, before
  any B1 run)*: arm-E intervals are computed by a within-chain circular
  block bootstrap, block length `ceil(2 tau)` from each chain's OWN
  truth-series autocorrelation (shared Geyer diagnostic); blocks never
  cross a chain boundary; arm S, being direct sampling, resamples
  ordinarily. Rationale: the frozen `n0` screen admits `tau <= 2.4`,
  and an iid bootstrap under that dependence narrows intervals by
  roughly `sqrt(tau)` — enough to falsely support an equivalence
  hypothesis. The screen itself is unchanged; the dependence is
  handled where the uncertainty is computed. Block lengths are
  reported per chain.
- **H-TRACK (equivalence, TOST).** At each rung, the E − S median-truth
  difference must have its 95% CI inside `[-0.05, +0.05]`. The margin
  is the programme's standing saturation tolerance (P2 and P9's frozen
  0.05), chosen because it is the difference the programme already
  treats as not meaningful in truth-error units — an a-priori anchor,
  not a fit. Supported only if every rung passes.
- **H-DEEPEN.** Arm E's own `N = 1200` minus `N = 600` median-truth
  difference has its 95% CI entirely below zero — the emergent phase's
  geometry *sharpens* under the resolution ladder, not merely tracks.
- Read as a conjunction for "consistent with a continuum limit over the
  measured range"; each reported separately; H-TRACK failing at exactly
  the top rung is the Stage A tension resolving into a detection and
  would be a major honest finding, not a failure.
- Exact hypothesis values above are frozen by this section now; only
  B0's pass/fail decides *whether* B1 runs, never *what* it tests.

### 8.4 Costs, measured basis

B0: minutes (direct sampling; the scaled fit at `N = 1200` costs a few
seconds per sample). B1: 1.3M steps per chain at the measured rates —
1.3/2.6/4.5 hours per chain at 600/900/1200 — about 17 chain-hours,
run as six parallel processes, wall ~5 h.

### 8.5 Stage B0 outcome (factual record, 2026-07-27): the gate FAILED
and Stage B stops here

Seeds 1100-1159, 60/60 samples ok. Under the scaled instrument the
sprinkled arm's truth error does not fall — it **rises, significantly**:

| `N` | median truth [95% CI] | targets | chain floor |
|---|---|---|---|
| 600 | 0.1599 [0.1335, 0.1807] | 44 | 25 |
| 900 | 0.1623 [0.1428, 0.1816] | 66 | 31 |
| 1200 | 0.1907 [0.1807, 0.2011] | 88 | 35 |

Top-minus-bottom: **+0.031 [+0.004, +0.058]** — the entire CI above
zero, the opposite sign of the gate's requirement. Per 8.2, **no B1
chain runs and the frozen 8.3 hypotheses are not evaluated.** Artifacts:
`docs/prereg/frozen/p10_stage_b/`.

**What this is a finding about.** Two instrument families have now been
tried against the same premise and both refuse it: frozen constants
leave the sprinkled curve flat (Stage A — scale-covariance), and
continuum-scaled constants make it rise (B0). **Why it rises is not
established.** A first version of this note offered a mechanism — that
densifying targets shrinks typical pair separations, so the instrument
asks harder questions as fast as it improves — and that mechanism is
**wrong for this observable and is withdrawn**: the truth error samples
comparison pairs *uniformly over all targets*
(`embedding_distance_order_error`), so the typical separation of a
sampled pair is set by the region, not by the target count, and does
not shrink as targets densify (review finding). What survives as a
candidate, still labelled post hoc and untested: the fit-optimization
budget is unscaled in this family (`STEPS = 1500`, `RESTARTS = 5` fixed
while free parameters and constraints double), so larger fits may
simply be less converged. Distinguishing that from an intrinsic effect
is a diagnostic for any successor stage, not something this record can
assert.

**What would make the Section 2 table decidable, and what it would
take.** A *scale-fixed* observable — e.g. discordance restricted to
pairs whose true continuum separation exceeds a threshold held fixed
across `N` — asks the same question at every size while the resolution
improves. That motivation stands on its own construction and does not
depend on the withdrawn mechanism. It is a new observable, so it
requires its own design section with the same discipline (gate first on
arm S, hypotheses frozen before any E chain), not a quiet swap inside
this one; and its arm-S gate is a genuine test, since this stage has
now twice failed to produce a falling yardstick by rescaling alone.
Whether to open that Stage B′ is a programme decision, recorded here as
the open question this stage ends on.

*Correction note (2026-07-27, review round on PR #24).* (i) The
mechanism paragraph above originally asserted "scale-referenced by
construction" via the pair-separation argument; withdrawn as described
— what is measured is flat-then-rising, and the explanation is open.
(ii) The `--stage b1` path originally trusted the documentation to
prevent post-gate runs; it now reads the frozen B0 record and refuses
unless `yardstick_falls` is true, the same enforcement pattern as
P8-B's constants guard. (iii) The B1 seed table 1000-1050 collided
with P6-B's uniform reference seeds 1000-1019 — the (600, random)
start would have reused a P6-B reference order exactly — and is
replaced by 30000-30050, verified clean including the derived streams;
and the aggregator now validates every chain against the frozen table
(membership and a uniform matching `chain_seed`), so a stale, extra,
or wrong-seed CSV exits instead of entering the screen.

## 9. Stage B': the scale-fixed observable (designed 2026-07-27;
**B'0 ran the same day and its gate FAILED — see 9.5; P10 is closed**)

Stage B ended on a diagnosis: both instrument families kept the same
observable — sign discordance over uniformly sampled target pairs —
whose difficulty mix is fixed by the continuum distribution of
comparison margins and therefore never eases as resolution improves.
B' inverts the design: **the instrument is the frozen P3 discriminator,
verbatim and unscaled** (which also removes B0's unscaled-fit-budget
confound: same constants, same convergence regime at every size), and
only the truth *scoring* changes.

### 9.1 The observable, and its two anchors

**Margin-restricted discordance**: the same pair-pair comparisons as
the frozen scorer, restricted to quadruples whose TRUE comparison
margin `||x_a - x_b| - |x_c - x_d||` (continuum units, `x = (i -
pi_i)/N`) is at least `DELTA_MARGIN`. Every size is asked the same
fixed question set while per-question precision improves with density.

Two anchors keep it honest:

- **To the frozen observable**: at `delta = 0`, same seed and count,
  the scorer reproduces the frozen scorer *bit for bit* (same pair
  stream, same reference-tie exclusion) — pinned by regression. It is
  a restriction of the frozen observable, not a second definition.
- **To pure geometry**: `DELTA_MARGIN = 0.1399` is the 25th percentile
  of the comparison-margin distribution of the exact continuum model
  (uniform iid points in the unit `(u, v)` square), computed from
  geometry alone — 4M draws, seed 2718281828, reproducible to 3.5e-6
  across seeds. **No experiment data was consulted in choosing it.**
  Eligible fraction is 75% by construction under the continuum model;
  the realized fraction is reported per sample with a floor of 8,000
  eligible comparisons out of 32,000 sampled.

### 9.2 Stage B'0 — the gate, third attempt at the yardstick

Arm S through the frozen instrument, both scorings per sample (the
unrestricted scoring rides along as the replication of Stage A's flat
curve on the very same fits). 20 samples per size, seeds 40000-40059,
fresh against every range in the programme including P6-B's 1000-1019
and the B1 envelope 30000-30379 (pinned by regression).

**Gate**: the restricted discordance's top-minus-bottom CI entirely
below zero, with full completeness at every rung and the eligible floor
met on every sample — completeness-conditioned verdicts are not
verdicts (the B0 lesson, inherited into `bprime_gate_verdict`).

If this gate fails too, the scale-referenced behaviour is not a
property of the question mix either, three distinct mechanisms have
been eliminated (instrument scale-covariance, instrument rescaling,
question-mix fixity), and P10 closes on that three-sided instrument
finding.

### 9.3 Stage B'1 — frozen now, runnable only on a B'0 pass

Identical in structure to the (never-run) 8.3, with the restricted
observable substituting for the unrestricted one end to end: six
chains at fresh seeds (to be drawn clear of every used range and
listed before any chain runs), 48 retained samples at 25k spacing,
the corrected Stage B machinery verbatim — code-enforced B'0 gate,
frozen-chain validation, index-set completeness, measurement
completeness, dual-start presence, within-chain block bootstrap on the
restricted truth series, ESS >= 20 at m = 48. Hypotheses: **H-TRACK'**
(TOST at margin 0.05 per rung on restricted E - S differences; the
margin is P2/P9's standing tolerance, same provenance as 8.3) and
**H-DEEPEN'** (arm E's own top-minus-bottom restricted CI entirely
below zero); conjunction for "consistent with a continuum limit over
the measured range". The B'1 runner is deliberately NOT implemented
until the gate passes — a 17-chain-hour path with nothing to authorize
it is an invitation.

### 9.4 Cost

B'0: minutes (direct sampling; the dual scoring adds one 32k-comparison
pass per sample). B'1 if authorized: the 8.4 budget unchanged.

### 9.5 Stage B'0 outcome (factual record, 2026-07-27): the gate FAILED
a third time, and P10 closes on its own rule

Seeds 40000-40059, 60/60 samples ok, every eligible-comparison floor
met (minimum 17,057 of 32,000 against a floor of 8,000).

| `N` | unrestricted (replication) | restricted (gates) |
|---|---|---|
| 600 | 0.1462 [0.1310, 0.1608] | 0.0520 [0.0304, 0.0787] |
| 900 | 0.1367 [0.1236, 0.1553] | 0.0457 [0.0306, 0.0672] |
| 1200 | 0.1703 [0.1599, 0.1836] | 0.0892 [0.0843, 0.0981] |

Top-minus-bottom, restricted: **+0.037 [+0.009, +0.062]** — the entire
CI above zero, again the wrong sign. (Unrestricted, same fresh seeds:
+0.024 [+0.005, +0.045], consistent in direction with Stage A's flat
curve and here crossing into significance; the rise is driven by the
`N = 1200` rung, with `N = 900` dipping below `N = 600` on both
scorings.) Per the frozen 9.2 rule, **B'1 is not implemented, its
hypotheses are never evaluated, and P10 closes.**

**The three-sided finding.** The margin restriction did what it was
designed to do at every individual size — the cutoff-restricted
questions are about three times easier than the full mix (0.05-0.09
against 0.14-0.17) — and the improvement with `N` still never came.
**A scope correction from review, applied before closing:** a common
cutoff restricts the question set but does not FIX its distribution —
the surviving margin mix above the cutoff can drift with the targets'
spatial configuration across sizes, and the realized eligible rates
(62.6% at `N = 600` to 66.2% at `N = 1200`, against the continuum 75%)
show that it does. The frozen gate's verdict is unaffected — the gate
was defined on the cutoff-restricted observable and failed as frozen —
but the mechanism-elimination claim needed the mix drift addressed,
which the post-hoc stratified diagnostic of 9.6 does. Three
measurements so far:

| attempt | instrument | question mix | sprinkled curve |
|---|---|---|---|
| Stage A | frozen (scale-covariant) | full | flat |
| Stage B0 | continuum-scaled | full | rising |
| Stage B'0 | frozen | **cutoff-restricted** | **rising** |

The positive characterization, which is what P10 ends on: **through
this discriminator family, truth accuracy in continuum units does not
improve — and measurably degrades — as the order grows**, under every
combination tried of instrument scaling and question fixing. The
discriminator is a fixed-relative-resolution instrument. The Section 2
outcome table is therefore undecidable by this instrument family, full
stop: a continuum-limit test requires a reconstruction whose continuum
accuracy grows with density, which is a different instrument, not a
different protocol around this one.

Why the fixed-question error rises is **not established**. The
unscaled-fit-budget candidate from 8.5 does not directly apply (the
constants, budgets and convergence regime are identical at every size
here); what grows with `N` at frozen constants is the selected chains'
lengths and hence the raw delay-rank magnitudes feeding the
dissimilarity. That observation is recorded as a lead, labelled post
hoc and untested, not as a mechanism.

What this does NOT touch: every per-size result in the programme —
PC-V1 through P9 and P10 Stage A's tracking result — is a
fixed-size statement and stands unchanged. What closed is the road
from this instrument to a continuum limit.

Artifacts: `docs/prereg/frozen/p10_stage_bprime/`.

### 9.6 Post-hoc stratified diagnostic (2026-07-27): the rise survives
mix-matching, bin by bin

Review of 9.5 identified a scope error: a common margin cutoff
restricts the question set but does not fix its distribution, so the
B'0 rise could in principle have been drift in the surviving-margin mix
rather than degrading accuracy. This diagnostic settles that. The same
60 samples (identical seeds, deterministic) were re-scored per a-priori
margin bin — conditional quartiles of the continuum margin given
`margin >= delta`, edges `0.2633 / 0.4160 / 0.6367` from the same
pure-geometry stream as `DELTA_MARGIN` — and the equal-weight bin mean
is the mix-matched error, immune to mix drift by construction.

| `N` | mix-matched [95% CI] | bin 1 | bin 2 | bin 3 | bin 4 |
|---|---|---|---|---|---|
| 600 | 0.0392 [0.0189, 0.0609] | 0.1234 | 0.0340 | 0.0011 | 0.0000 |
| 900 | 0.0346 [0.0226, 0.0561] | 0.1037 | 0.0259 | 0.0010 | 0.0000 |
| 1200 | 0.0750 [0.0696, 0.0908] | 0.1930 | 0.0898 | 0.0175 | 0.0012 |

Top-minus-bottom, mix-matched: **+0.036 [+0.012, +0.060]** — the rise
survives, and at `N = 1200` the error is higher than at `N = 600` in
every one of the four bins. *(Corrected on second review: four broad
bins still admit WITHIN-bin drift, and error varies sharply with
margin, so this subsection alone reduces rather than eliminates the
mix-drift explanation. The elimination is done properly by the
pointwise curve of 9.7.)*

Status of this subsection: post hoc, run after the frozen 9.2 gate had
already failed; it rewords the closure's mechanism claim and cannot
reopen the gate. Artifacts alongside the B'0 record in
`docs/prereg/frozen/p10_stage_bprime/`.

### 9.7 Post-hoc pointwise curve (2026-07-27): domination across the
margin axis, which no reweighting objection can touch

Second review escalation on 9.6: equal-weighting four broad bins leaves
within-bin drift, and the steep dependence of error on margin makes
that a live concern. The decisive form is the pointwise
error-versus-margin curve — nothing is aggregated across margins, so
there is nothing to reweight. Same 60 samples, 16 a-priori fine bins
(continuum 16-quantiles of the full margin range, same geometry stream,
cross-seed delta <= 5e-4), per-bin medians over the 20 samples with
bootstrap intervals.

| margin bin | E(600) | E(900) | E(1200) | 1200 above 600 | CI-sep |
|---|---|---|---|---|---|
| [0.000, 0.033) | 0.4548 | 0.4550 | 0.4691 | Y | |
| [0.033, 0.067) | 0.3558 | 0.3600 | 0.3956 | Y | Y |
| [0.067, 0.103) | 0.2761 | 0.2826 | 0.3383 | Y | |
| [0.103, 0.140) | 0.2165 | 0.2110 | 0.2855 | Y | Y |
| [0.140, 0.179) | 0.1608 | 0.1435 | 0.2351 | Y | Y |
| [0.179, 0.220) | 0.1241 | 0.0995 | 0.1907 | Y | Y |
| [0.220, 0.263) | 0.0879 | 0.0603 | 0.1581 | Y | Y |
| [0.263, 0.310) | 0.0531 | 0.0374 | 0.1196 | Y | |
| [0.310, 0.361) | 0.0318 | 0.0290 | 0.0905 | Y | Y |
| [0.361, 0.416) | 0.0135 | 0.0113 | 0.0625 | Y | |
| [0.416, 0.478) | 0.0023 | 0.0025 | 0.0308 | Y | |
| [0.478, 0.550) | 0.0006 | 0.0006 | 0.0154 | Y | |
| [0.550, 0.637) | 0.0000 | 0.0000 | 0.0060 | Y | |
| [0.637, 0.747) | 0.0000 | 0.0000 | 0.0022 | Y | |
| [0.747, 0.914) | 0.0000 | 0.0000 | 0.0000 | | |
| [0.914, inf) | 0.0000 | 0.0000 | 0.0000 | | |

*(Bin labels here and in 9.9 are display-rounded to 3 decimals; the
exact a-priori edges — 0.4784, 0.5503, 0.6367, 0.7473, 0.9144, … — are
in the runner's constant block and the artifact. Claim-bearing
thresholds in the text use the exact values.)*

**The `N = 1200` curve sits above the `N = 600` curve in 14 of 16
bins — every bin in which either curve is nonzero — with CI separation
in 6.** The two exceptions are the widest-margin bins where both
curves are exactly zero, leaving no room to be above. (The `N = 900`
curve dips at or below `N = 600` through the mid bins, consistent with
the 900 dip seen on every aggregate; the rise is a 1200 phenomenon,
now visible pointwise.) Sparse-bin caveat: the widest bins carry few
comparisons per sample for bracketed targets, and per-bin minimum
counts are recorded in the artifact; the domination conclusion rests
on the well-populated bins.

This settles the mechanism claim in its strongest available form:
**at fixed question difficulty — pointwise in margin — the frozen
instrument answers less accurately at `N = 1200` than at `N = 600`.**
Per-question accuracy in continuum units degrades with `N`; no
mix-drift account survives a pointwise comparison. Status as 9.6:
post hoc, after the frozen gate failed; rewords the mechanism claim,
cannot reopen the gate. Artifact:
`docs/prereg/frozen/p10_stage_bprime/p10_bprime_margincurve_summary.json`.

### 9.8 Third escalation answered by bound and by zero (2026-07-27)

Review pressed once more: any finite bin averages within itself, so
within-bin drift remains logically possible at any granularity. Two
answers, one quantitative and one absolute, both now in the artifact:

**The drift bound.** The most a within-bin margin shift can move a bin
average is `|mean-margin shift| x local slope`, and both factors are
measured. The realized within-bin mean-margin shifts between `N = 600`
and `N = 1200` are 0.00001-0.00125 — inside a fine bin the margin
distribution has almost no room to move — so the maximal drift
contribution is at most 0.0007 in every bin, while the observed rises
are 0.002-0.074: **the rise exceeds the maximal drift contribution in
all 14 nonzero bins, by one to three orders of magnitude** (bin 4:
+0.0743 against a bound of 0.00003). *(Corrected on fourth review:
this is not a bound — equal mean margins can hide mass reallocation
within a bin, and a neighbour-bin slope is no Lipschitz constant for
the within-bin error curve. The "maximally 0.0007" figure and the
"eliminated" conclusion drawn from it are withdrawn; the quantity is
retained in the artifact as a descriptive indicator only. See 9.9.)*

**The zero bins.** In bins 12-13 (margins 0.478-0.637... and the
adjacent 0.637-0.747 in the medians table) the `N = 600` error is
identically zero: an average of nonnegative indicators equal to zero
means zero discordance at every sampled margin in the bin, and no
reweighting of an identically-zero function is positive. `N = 1200` is
positive there. At those margins, fixed-difficulty degradation holds
with no binning assumption of any kind. *(Corrected on fourth review:
those table entries are MEDIANS of 20 per-sample rates, not pooled
averages — a zero median only means at least half the per-sample rates
are zero. The pooled counts show 499 discordant `N = 600` comparisons
in this region; the identically-zero claim is withdrawn as stated. The
per-comparison zero-fact that does hold — 0/3,871 in the last bin — is
established from pooled counts in 9.9.)*

**Provenance correction (same round).** Both diagnostic artifacts had
been produced from then-uncommitted implementations and stamped a clean
HEAD that did not contain them — the same class of error as correction
(xii) of Section 6d, recurring. Both were regenerated at a commit that
contains their implementation, with medians verified IDENTICAL to the
prior freezes; and diagnostic stamps are now dirty-aware
(`-dirty` suffix), so an uncommitted implementation can never again
stamp a commit it is not in.

### 9.9 Fourth round: both 9.8 arguments corrected, and the exact
pooled analysis that replaces them (2026-07-27)

Review corrected both of 9.8's arguments, and both corrections are
right. The "drift bound" was not a bound — a mean-margin shift times a
neighbour-bin slope bounds nothing, since equal means can hide mass
reallocation and the neighbour slope is no Lipschitz constant; it is
demoted to a labelled descriptive indicator. And a zero MEDIAN of
per-sample rates does not mean zero discordance — concretely so: the
pooled counts reveal 499 discordant comparisons at `N = 600` in the
high-margin region whose bin medians read 0.0000. The 9.8 "zero bins"
claim is withdrawn as stated.

What replaces them is exact and assumption-free: **pooled discordant
counts over every sampled comparison**, per rung per bin —

| margin bin | N=600 | N=900 | N=1200 | rate ratio 1200/600 |
|---|---|---|---|---|
| [0.416, 0.478) | 772/36167 (.0213) | 754/36232 (.0208) | 1885/37905 (.0497) | 2.3x |
| [0.478, 0.550) | 347/30786 (.0113) | 416/33377 (.0125) | 1114/34930 (.0319) | 2.8x |
| [0.550, 0.637) | 121/26060 (.0046) | 189/28930 (.0065) | 538/32159 (.0167) | 3.6x |
| [0.637, 0.747) | 26/20414 (.0013) | 46/23151 (.0020) | 200/27596 (.0072) | 5.7x |
| [0.747, 0.914) | 5/13333 (.0004) | 17/15594 (.0011) | 84/22518 (.0037) | 9.9x |
| **[0.914, inf)** | **0/3871** | **0/4799** | **15/11184** | — |

Three statements now carry the mechanism claim, in decreasing order of
absoluteness:

1. **The last bin is a true zero-fact at the per-comparison level.**
   Every one of the 3,871 sampled comparisons at margins `>= 0.9144`
   was answered correctly at `N = 600` (and all 4,799 at `N = 900`);
   `N = 1200` errs there 15 times. No median is involved and no
   within-bin reweighting story survives an empirical error of exactly
   zero across the whole sampled bin. *(Corrected on fifth review:
   "the whole sampled bin" is the pool's sampled support, not the bin —
   the rungs sample different comparison sets, and the artifact records
   different last-bin mean margins (0.9756 vs 1.0045). The zero is an
   empirical pooled sample fact; it does not establish an
   identically-zero response function over `[0.9144, inf)`. Scoped
   properly in 9.10.)*
2. **Every tail bin's pooled rate rises 2.3-9.9x** from `N = 600` to
   `N = 1200`, on counts in the tens of thousands. Pooled comparisons
   share fits, so these are exact sample facts rather than fresh
   inference — the inferential statements remain the per-sample
   bootstrap results of 9.7.
3. The within-bin mean-margin shifts (0.00001-0.00125) remain recorded
   as a **descriptive indicator only**.

Also recorded rather than hidden: the m* construction of the previous
commit was uninformative as designed — `m* = 0.4792`, essentially the
region edge, because `N = 600` errs near the edge too; the per-bin
pooled table above supersedes it. Artifact regenerated at a commit
containing its implementation (stamp `3a9edc3`, medians verified
identical to the prior freeze), under the now-structural dirty-aware
provenance rule.

**The closure claim, in its final calibrated form.** The frozen gates
failed as frozen (unchanged throughout). The mechanism reading is:
per-question accuracy in continuum units degrades with `N` — exactly,
at the highest margins, where `N = 600` is empirically error-free over
thousands of sampled comparisons and `N = 1200` is not; and by 2-10x
pooled-rate rises across the tail, with the pointwise CI separations of
9.7 in the mid bins. What a mix-drift alternative would now require has
no support in any recorded quantity and cannot operate at all in the
last bin. *(Corrected on fifth review: "cannot operate at all"
overstates — a reweighting account confined to within-bin regions
unrepresented in the `N = 600` sample is not excluded by a
sample-scoped zero. The terminal scoped form is 9.10.)*

### 9.10 Fifth round: the zero-fact scoped to its sample (2026-07-27)

Review pressed on 9.9's statement 1, and the point is right: the
`N = 600` and `N = 1200` pools sample different target geometries and
different comparisons, and the artifact itself records different
last-bin mean margins (0.9756 versus 1.0045). Zero observed errors
over the 3,871 comparisons the `N = 600` pool sampled is an empirical
pooled sample result — it does not establish an identically-zero
response function over `[0.9144, inf)`, and a reweighting account in
which the `N = 1200` errors sit in within-bin regions unrepresented in
the `N = 600` sample is not excluded by it. "Cannot operate at all in
the last bin" is withdrawn.

Two recorded quantities bear on that account, both descriptive: the
only measured within-bin location statistic, the mean margin, moves
toward larger margins at `N = 1200` (0.9756 -> 1.0045) — recorded as
description, since a mean does not localize mass; and the degradation
is CI-separated per sample in six mid bins (9.7), where the claim is
inferential rather than support-dependent.

The mechanism claim in its terminal form, claiming nothing beyond
sampled support: **the frozen gates failed as frozen; per-question
accuracy in continuum units degrades with `N` — inferentially in the
mid bins (9.7's bootstrap-separated pointwise rises), and as exact
sample facts in the tail (pooled rates rising 2.3-9.9x on tens of
thousands of counts; `N = 600` error-free over its 3,871 sampled
last-bin comparisons while `N = 1200` errs 15 times in 11,184).**
Each escalation of this review cycle has now been driven to a quantity
that is exact over its stated scope, and to no claim beyond it.

### 9.11 Sixth round: the scorer-seed namespace collision, conceded
and measured (2026-07-27)

Review found a seed-hygiene violation in the B'0 design itself, and it
is real: the derived scorer seeds `base + k + 9` sit inside the
consecutive chain-seed spans, so **51 of the 60 scorer streams
coincide with some sample's chain stream** (values 40009-40059 — e.g.
40020 selects scoring pairs for `N = 600` sample 11 AND drives
`N = 900` sample 0's chain; 40009 scores `N = 600` sample 0 and drives
its sample 9). The freshness test checked derived seeds against every
EARLIER range and never against the B'0 chain seeds themselves — the
same class of hole as the B1 collision of Section 8, one level deeper.
Rows entering the gate bootstrap as independent replicates therefore
shared deterministic RNG streams. The fits are untouched by this: the
chain seeds are mutually distinct, and the only colliding randomness
is the scorer's pair SELECTION.

The seedfix stage (`--stage seedfix`, regression-pinned collision map)
replays the same 60 fits and re-scores **every published B' quantity**
under both namespaces — the original (an exact replay check) and a
disjoint one, 41000-41059, fresh against everything including
ourselves. *(Corrected on seventh review: the first seedfix run copied
the unrestricted truth from the instrument's own output — which is
scored internally at `seed + 9`, the same colliding stream — so
exactly one published quantity had escaped this sweep while this
paragraph said "every". Re-scored in 9.12.)* Frozen artifacts:
`docs/prereg/frozen/p10_stage_bprime/p10_bprime_seedfix.csv` and
`_summary.json` (stamp `b5b0af5`). Replay: the old-namespace medians
reproduce the frozen B'0 summary to 0.0 and the pooled curve counts
match the frozen artifact exactly, so the comparison is
apples-to-apples. The results, old vs clean:

| quantity | old namespace | clean namespace |
|---|---|---|
| gate: restricted top-minus-bottom | +0.0372 [+0.0091, +0.0630] | +0.0379 [+0.0097, +0.0624] |
| gate verdict (falls) | NO (rise) | NO (rise) |
| eligible floors | met (>= 16847) | met (>= 16847) |
| 9.6 mix-matched rise | +0.0358 [+0.0125, +0.0589] | +0.0340 [+0.0129, +0.0577] |
| 9.7 pointwise: bins 1200 above 600 | 14/16, 6 CI-separated | 14/16, 6 CI-separated |
| 9.9 last bin, N = 600 | 0/3871 | 0/3967 |
| 9.9 last bin, N = 900 | 0/4799 | 2/4749 |
| 9.9 last bin, N = 1200 | 15/11184 | 28/11066 |

Row-level maxima of |clean − old| are 0.0036 / 0.0044 / 0.0069
(restricted, by rung) and at most 0.0044 (mix-matched) — the
pair-sampling noise scale, an order below the between-sample spread
that drives every interval. Every direction, verdict, and CI sign is
unchanged; the pooled tail rates rise 2.5-19x under clean seeds
(sparser top bins swing the ratio; the direction is identical).

Owned rather than smoothed: under fresh pair draws `N = 900`'s last
bin reads 2/4749 where the old draw read 0/4799 — a pooled zero is a
sample fact that need not replicate, exactly as 9.10 scoped it; the
`N = 600` zero DID replicate (0/3967), and `N = 1200`'s last-bin rate
is ten times `N = 900`'s under the same clean namespace. Status: post
hoc, labelled; the frozen record stands as run, the gate stays
failed-as-frozen, and the terminal mechanism claim of 9.10 is
unchanged — now measured to be robust to the scorer-seed namespace,
rather than argued to be.

### 9.12 Seventh round: the unrestricted arm, the one quantity the
sweep had missed (2026-07-27)

Review caught 9.11 doing what it had just corrected others for: the
first seedfix run copied `row["truth"]` from the instrument — which
scores truth INTERNALLY at `seed + 9`, the same colliding stream — so
the unrestricted arm (the Stage A replication reading) was the one
published quantity outside the robustness sweep, under a paragraph
saying "every". Conceded; the run now re-scores it under both
namespaces, with an in-run assertion that the old-namespace
recomputation equals the instrument's own output bit for bit (it
does; and the old truth medians reproduce the frozen B'0 summary to
0.0). Artifacts refrozen at stamp `702f4a5`; the disjointness test now
also clears the instrument's other derived windows (`seed + 100` fit
learning).

| unrestricted (8,000-comparison budget) | old | clean |
|---|---|---|
| median truth 600 / 900 / 1200 | 0.1462 / 0.1367 / 0.1703 | 0.1513 / 0.1369 / 0.1676 |
| top-minus-bottom | +0.0241 [+0.0051, +0.0461] | +0.0163 [+0.0069, +0.0404] |

The ordering (1200 > 600 > 900, the Stage A shape) and the
top-minus-bottom sign with an all-positive CI hold in both namespaces.
The point estimate moves by 0.008 and the row-level maxima of
|clean − old| are 0.010-0.016 — three to four times the restricted
quantities' shifts, exactly as the 8,000-comparison budget predicts
(per-row pair-selection SE ~ 0.004 against the restricted scorers'
~ 0.002): the unrestricted arm is the noisiest published quantity, and
its namespace sensitivity is pair-sampling noise at that budget, with
no direction or verdict change anywhere.

