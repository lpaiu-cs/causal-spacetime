# P14 §8 P3-C — confirmation stage: **confirmed** (exploratory chain)

The confirmation/termination stage of the P3 discriminability
question, run under the machinery frozen and certified in the
preflight (PR #47): three-branch verdict, DeLong AUC, full-variance
`s`, midpoint-threshold BA with first-half training; margins
`ε_s = 0.0806, ε_AUC = 0.0233, ε_BA = 0.0285`; both branches
certified at CP-lower ≥ 0.90 before execution. Fresh unpaired seed
streams (curved `20260841`, flat `20260842`), `n = 4800`
sprinklings per arm at aniso-a1.0, `E[N] = 300`.

This is the SECOND campaign block. The first (seeds
`20260831/20260832`) was downgraded to exploratory after review —
see "Protocol correction" below — and its record is preserved in
`p14_probe_p3c_results_exploratory.json`. The verdict below is the
program's confirmation verdict; it comes from seeds never used
before, under the corrected interval rule frozen BEFORE these
streams were opened.

## Verdict

| metric | value | 95% CI | band |
|---|---|---|---|
| s | 10.988237 | [10.827728, 11.148746] | ±0.0806 |
| AUC | 1.000000 | [0.999232, 1.000000] | 0.5 ± 0.0233 |
| BA | 1.000000 | [0.985855, 1.014145] | 0.5 ± 0.0285 |

**Verdict: `confirmed`** — every CI entirely outside its band in the frozen direction.

The two arms' 3a samples are again completely separated in the raw
data (min over 4800 curved values `0.0539969` exceeds max over 4800
flat values `0.0329959`; means `0.074131` vs `0.023804`). At
complete separation the DeLong influence variance collapses, so the
AUC interval switches to the frozen boundary rule: the one-sided
exact lower bound `0.025^(1/m) = 0.999232` for `m = 4800` — the
smallest population AUC not rejected at 2.5% by `m` concordant
paired placements. The interval `[0.999232, 1]` is a valid
finite-sample statement; its lower endpoint clears the band by a
wide margin, so the verdict does not lean on the boundary rule's
tightness. Both arms recorded `ambiguous = 0, escalated = 0`
relation censuses (stored and asserted in the artifact).

## First campaign (downgraded to exploratory)

The original block on seeds `20260831/20260832` produced the same
qualitative picture — complete separation, `s = 11.103838`
`[10.941745, 11.265932]`, BA `1.0 [0.985855, 1.014145]`, means
`0.074154` vs `0.023772` — but its AUC interval was the degenerate
Wald `[1, 1]`: the estimated influence variance is zero at complete
separation, which erases the finite-sample uncertainty about
unobserved population overlap rather than measuring it. The frozen
verdict rule requires valid 95% CIs for all three metrics, so that
block's verdict is **downgraded from `confirmed` to
`exploratory-separation`**. The record is preserved unedited (plus
a `grade` field) in `p14_probe_p3c_results_exploratory.json`; the
seeds are burned and appear in `BURNED_SEEDS`.

## The licensed sentence (frozen scope)

> "aniso-a1.0의 고정된 유한 상자·밀도에서 global relation fraction이 평탄 ensemble과 분리된다"

Nothing wider: not a general-Weyl discriminator claim, not a
box-independent one. The three rules close only what they measure;
the termination sentence was not needed.

## What this decides

§8.3's positive branch: **P14 does not close here.** The actual
preregistration's power section is built from P3's effect sizes and
P4's paired-variance record (gain `1.532× [1.386, 1.698]`, preflight
artifact) — with S1's Schwarzschild price still to be attached to
the eventual scope statement.

Everything above the licensed sentence is RENDERED from the
committed artifact `p14_probe_p3c_results.json` by
`campaign_table`, and a test asserts the doc embeds it verbatim;
the metrics and verdict recompute from the stored raw samples, and
the seed streams reproduce on a prefix (slow marker).

## Changelog

**Protocol correction (PR #48 review R1).** The first campaign's
`confirmed` verdict relied on the degenerate `[1,1]` AUC interval,
and this document originally described that degeneracy as "a
property of the data, not of the interval construction". That
framing was backwards: complete separation is the data property
that CAUSES the Wald interval to degenerate; the degenerate
interval itself is not a valid finite-sample population-AUC CI.
Corrections, in order:

1. `auc_delong` now returns the exact boundary interval
   `[0.025^(1/m), 1]` when the empirical AUC is 1 (symmetric rule
   at 0). Interior behaviour is unchanged; the rule was frozen and
   the preflight re-certified (identical counts: null equivalent
   `18209/20000`, CP-lower `0.90641`; effect `4000/4000`) before
   any new data was drawn.
2. `arm_samples` now records per-arm `ambiguous`/`escalated`
   totals; the campaign asserts both are zero and stores them in
   the artifact, so undecided pairs can never silently enter the
   relation fraction as "unrelated".
3. The first campaign was downgraded, not erased: artifact renamed
   to `p14_probe_p3c_results_exploratory.json` with an explicit
   `grade`, seeds `20260831/20260832` burned.
4. A fresh confirmation block was run on new seeds
   `20260841/20260842` under the corrected rule — the verdict
   table above.

Initial record: first campaign, seeds `20260831/20260832`
(superseded as the confirmation verdict by this correction).
