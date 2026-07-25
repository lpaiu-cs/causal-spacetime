# P10: does the emergent geometry deepen toward a continuum?

Status: **DESIGN; STAGE 0 PILOT IN PROGRESS; NOTHING FROZEN** (2026-07-26).

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

## 7. What this can and cannot answer

It can answer: whether the continuum phase's reconstructable geometry
behaves, through a fixed instrument, like a geometry being resolved
(arm-S-like) or like an artifact being exposed, over one octave of size.

It cannot answer: the thermodynamic limit (three sizes are three sizes),
Lorentz structure, dynamics, or anything about ensembles other than this
one. If the outcome is "deepening", the honest statement is
"consistent with a continuum limit over the measured range" — and the
next octave is the follow-up, not a victory lap.
