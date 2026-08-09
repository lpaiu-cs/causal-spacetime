# P14 §8 P3-C — confirmation stage: **confirmed** (exploratory chain)

The confirmation/termination stage of the P3 discriminability
question, run under the machinery frozen and certified in the
preflight (PR #47): three-branch verdict, DeLong AUC, full-variance
`s`, midpoint-threshold BA with first-half training; margins
`ε_s = 0.0806, ε_AUC = 0.0233, ε_BA = 0.0285`; both branches
certified at CP-lower ≥ 0.90 before execution. Fresh unpaired seed
streams (curved `20260831`, flat `20260832`), `n = 4800`
sprinklings per arm at aniso-a1.0, `E[N] = 300`.

## Verdict

| metric | value | 95% CI | band |
|---|---|---|---|
| s | 11.103838 | [10.941745, 11.265932] | ±0.0806 |
| AUC | 1.000000 | [1.000000, 1.000000] | 0.5 ± 0.0233 |
| BA | 1.000000 | [0.985855, 1.014145] | 0.5 ± 0.0285 |

**Verdict: `confirmed`** — every CI entirely outside its band in the frozen direction.

The two arms' 3a samples are COMPLETELY separated in the raw data
(min over 4800 curved values exceeds max over 4800 flat values;
mean `0.074154` vs `0.023772`), which is what the
degenerate AUC CI reflects — test-pinned as a property of the data,
not of the interval construction.

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

Initial record.
