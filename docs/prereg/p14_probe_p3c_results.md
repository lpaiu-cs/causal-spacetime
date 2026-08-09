# P14 §8 P3-C — confirmation stage: **confirmed** (exploratory chain)

The confirmation/termination stage of the P3 discriminability
question, run under the machinery frozen and certified in the
preflight (PR #47): three-branch verdict, DeLong AUC with the exact
disjoint-pair boundary rule, full-variance `s`,
midpoint-threshold BA with first-half training; margins
`ε_s = 0.0806, ε_AUC = 0.0233, ε_BA = 0.0285`; both branches
certified at CP-lower ≥ 0.90 before execution. Fresh unpaired seed
streams (curved `20260851`, flat `20260852`), `n = 4800`
sprinklings per arm at aniso-a1.0, `E[N] = 300`. The block was
executed from a clean checkout of the freeze commit `6f74ccb` —
which carries the corrected rule and these seeds but NO results —
and recorded in the commit after it, so the freeze-before-execution
ordering is provable from the commit graph (PR #48 review R2).

## Verdict

| metric | value | 95% CI | band |
|---|---|---|---|
| s | 11.348905 | [11.183450, 11.514360] | ±0.0806 |
| AUC | 1.000000 | [0.999232, 1.000000] | 0.5 ± 0.0233 |
| BA | 1.000000 | [0.985855, 1.014145] | 0.5 ± 0.0285 |

**Verdict: `confirmed`** — every CI entirely outside its band in the frozen direction.

The two arms' 3a samples are again completely separated in the raw
data (min over 4800 curved values `0.0535627` exceeds max over 4800
flat values `0.0327214`; means `0.074191` vs `0.023807`). At
complete separation the DeLong influence variance collapses, so the
AUC interval switches to the frozen boundary rule — the exact FIXED
DISJOINT-PAIR bound: index-pairs `(a_i, b_i)`, `i < k = min(m, n)`,
fixed a priori, give iid `Bernoulli(P(a > b))` indicators with
`P(a > b) ≤ AUC`, all successes at complete separation, hence the
exact CP lower bound `0.025^(1/4800) = 0.999232`. (Not a placement
argument — the placements share both arms and are not independent
trials.) The interval `[0.999232, 1]` is a valid finite-sample
statement; its lower endpoint clears the band by a wide margin, so
the verdict does not lean on the boundary rule's tightness. Both
arms recorded `ambiguous = 0, escalated = 0` relation censuses
(stored and asserted in the artifact).

## Earlier blocks (downgraded to exploratory, PR #48 review)

### First block (seeds `20260831/20260832`) — downgraded, R1

Complete separation (`s = 11.103838 [10.941745, 11.265932]`, BA
`1.0 [0.985855, 1.014145]`, means `0.074154` vs `0.023772`), but
its AUC interval was the degenerate Wald `[1, 1]`: the estimated
influence variance is zero at complete separation, which erases the
finite-sample uncertainty about unobserved population overlap
rather than measuring it. The frozen verdict rule requires valid
95% CIs for all three metrics, so the verdict is downgraded from
`confirmed` to `exploratory-separation`. Record:
`p14_probe_p3c_results_exploratory.json`.

### Second block (seeds `20260841/20260842`) — downgraded, R2

Run under the corrected boundary rule and again completely
separated (`s = 10.988237 [10.827728, 11.148746]`, AUC `1.0
[0.999232, 1.0]`, BA `1.0 [0.985855, 1.014145]`, means `0.074131`
vs `0.023804`, both arms `ambiguous = 0, escalated = 0`) — but the
corrected rule, the re-certified preflight, and this block's
results all entered git history in one commit (`5049ce8`), so the
repository's freeze-before-execution ordering (commit the runner,
then execute from a clean worktree) is not provable from the commit
graph. Test-level reproduction shows computational reproducibility,
not temporal ordering; the verdict is downgraded to
`exploratory-separation`. Record:
`p14_probe_p3c_results_exploratory2.json`.

Both blocks' streams are burned (`BURNED_SEEDS`); neither may feed
a confirmation again. That all three blocks — two exploratory, one
confirmation-grade — separate completely in the same direction is
consistent, but only the third carries the verdict.

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

**Protocol correction 2 (PR #48 review R2).**

1. The boundary lower bound's justification is corrected: the `m`
   placements share both arms and are not independent Bernoulli
   trials. The valid argument is the FIXED DISJOINT-PAIR bound
   (described under the verdict table); the exponent is now
   `min(m, n)` — identical at the campaign's equal arm sizes, so
   the number `0.999232` at 4800 is unchanged.
2. The second campaign block was downgraded, not erased, its seeds
   burned, and the operative campaign contract at the top of the
   runner updated to the new streams `20260851/20260852`.
3. The third block executed from a clean checkout of the freeze
   commit `6f74ccb`, results recorded in the subsequent commit —
   the verdict table above.

**Protocol correction 1 (PR #48 review R1).** The first campaign's
`confirmed` verdict relied on the degenerate `[1,1]` AUC interval,
and this document originally described that degeneracy as "a
property of the data, not of the interval construction". That
framing was backwards: complete separation is the data property
that CAUSES the Wald interval to degenerate; the degenerate
interval itself is not a valid finite-sample population-AUC CI.
The boundary rule replaced it, the preflight was re-certified
(identical counts: null equivalent `18209/20000`, CP-lower
`0.90641`; effect `4000/4000`), `arm_samples` began recording
per-arm `ambiguous`/`escalated` totals (asserted zero — undecided
pairs can never silently enter the relation fraction as
"unrelated"), and the first block was downgraded.

Initial record: first campaign, seeds `20260831/20260832`
(superseded as the confirmation verdict by correction 1).
