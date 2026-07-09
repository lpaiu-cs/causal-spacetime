# P3: does geometry emerge from a geometry-free dynamics?

Status: FROZEN v1 (2026-07-09). Calibration provenance: P3-A at commit f4837c8,
seeds 0-9. Frozen gates: `docs/prereg/frozen/p3_test_constants.json`
(heldout <= 0.10, null_gap >= 0.10; order-intrinsic discriminator re-validated
on sprinkled vs column-shuffle, effect size d = 10.4). After this point the
gates are locked and P3-B is the confirmatory application to the dynamics. The
first test of the emergence question the whole program builds toward: run the
validated representability discriminator on causal orders produced by a
DYNAMICS that never refers to geometry, rather
than on geometry that was sprinkled in (PC-V1/P2) or hand-diluted (P1). This
document lives in version control.

## 1. Motivation and honest prior

PC-V1 validated the discriminator; P1 showed it tracks graded geometric
consistency; P2 showed it works in 2+1D and selects the right dimension. All
used orders with geometry put in. The open question is whether a geometry-free
dynamics can PRODUCE order that the discriminator reads as geometric. We use
classical sequential growth / transitive percolation (Rideout-Sorkin), the
canonical geometry-free causal-set dynamics.

Honest prior: generic classical sequential growth is known to be
non-manifoldlike (dominated by Kleitman-Rothschild-type orders), so the honest
expectation is that the discriminator BLOCKS across the dynamics parameter — a
meaningful negative result ("geometry does not emerge from this dynamics"), not
a failure of the experiment. A regime that PASSES would be a striking positive
result demanding heavy scrutiny.

## 2. The chain-selection problem and how P3 handles it

A dynamics-generated causal set has no coordinates, so reference chains cannot
be supplied geometrically. P3 selects them intrinsically: the K longest
disjoint chains of the order, with tick labels = position along each chain
(`order_intrinsic.py`). Targets are elements two-sided bracketed by every
selected chain. Because this is a different measurement protocol from PC-V1
(supplied chains), the discriminator is RE-VALIDATED here with order-selected
chains (Section 4) before being applied to the dynamics.

## 3. Avoiding circularity

The gate must not be calibrated using the dynamics it will judge. P3 calibrates
the gate only on (a) SPRINKLED 1+1D causal sets with order-selected chains (the
pass cluster) and (b) a column-shuffle of those profiles (the geometry-free
fail cluster, marginal-preserving, consistency-destroyed) -- exactly the PC-V1
specificity foil, and independent of the dynamics. Transitive percolation is
the APPLICATION subject only (Section 5); it never enters gate construction.

## 4. P3-A calibration (re-validate the order-intrinsic discriminator)

Seeds 0-9. Sprinkle a 1+1D diamond (N = 1500, T = 2.0), select K = 6 longest
disjoint chains (min length 25), select up to 44 bracketed targets, measure
bracket-width parallax profiles, fit d = 1. Compute, for the structured
profiles and their column-shuffle: held-out violation, null gap
(shuffle heldout - structured heldout), and (structured only, against true x)
truth-order error.

Hard floors: HF1 median structured heldout <= 0.10 and median truth-order error
<= 0.35 (order-selected chains are noisier than supplied chains; recovery must
still be clearly sub-chance). HF2 the structured/shuffle separation must give a
positive null gap (Cohen's d >= 1.0).

Midpoint gate placement (as P2 D1):
- gate_heldout = min(0.10, round(0.5*(max structured heldout + min shuffle
  heldout), 0.05))
- gate_nullgap = max(0.10, round(0.5*(min structured null gap + 0), 0.05))
- gate_truth = round(0.5*(max structured truth + 0.5), 0.05)  (0.5 = chance)

A scene PASSES the geometry gate iff heldout <= gate_heldout AND
null_gap >= gate_nullgap. On calibration the structured scenes must pass and
the column-shuffles must block.

## 5. P3-B application (the emergence test)

Fresh seeds 100-119. For each transitive-percolation parameter p in a sweep
{0.006, 0.008, 0.010, 0.014, 0.020} (spanning wide/low-density to
higher-density orders), grow the order, select order-intrinsic chains and
targets, and evaluate the frozen geometry gate. Record for every (p, seed):
achieved density, whether K chains of the minimum length exist (structural
feasibility), held-out, null gap, and the pass/block verdict. No truth gate
(no coordinates).

Decision:
- EMERGENCE-SUPPORTED for a given p iff >= 16 of 20 seeds PASS the frozen
  geometry gate.
- NO-EMERGENCE for a given p iff >= 16 of 20 seeds block (gate fail or
  structural infeasibility).
- Report the full p-sweep; a mixed p is inconclusive at that p.

## 6. Stop rules

Instrument repair (chain selection, gates) only before the P3 freeze, via a
versioned amendment re-running P3-A on fresh seeds. The frozen PC-V1
dissimilarity/fit pipeline is unchanged. No threshold retuning after freeze; no
dropping seeds after outcomes; the p-sweep is fixed before P3-B.

## 7. Claim boundary

A NO-EMERGENCE result means: this specific geometry-free dynamics
(transitive percolation), judged by this validated discriminator at this scale,
does not produce low-dimensional geometry -- consistent with the known
non-manifoldlikeness of generic classical sequential growth. It is not a claim
about all possible dynamics. An EMERGENCE result would be a controlled,
finite-scale positive signal for one dynamics, not a proof of physical
spacetime emergence, and would require replication and scrutiny before any
stronger reading.

## 8. Freeze and provenance

Mirror PC-V1 Section 12: commit P3-A summary and gates to
`docs/prereg/frozen/p3_test_constants.json` with the calibration commit hash;
P3-B refuses to run without it. Every row carries order digests, seed, stage,
parameter, and code version.

## 9. Deviations log

(empty -- no pre-freeze repair was needed; P3-A calibration passed cleanly.)

## 10. Confirmatory outcome (post-freeze factual record)

Result of executing the frozen P3-B decisions (P3 gates at f4837c8; PC-V1
pipeline at 9162e8e), transitive percolation over the frozen p-sweep, seeds
100-119. Registry and per-(p, seed) CSV under `docs/prereg/frozen/`. Changes no
threshold or rule.

- **NO-EMERGENCE across the entire sweep: 20/20 block at every p** (0 passes of
  100 dynamics runs). Two regimes, both blocking:
  - Low p (0.006, 0.008; density 0.19-0.34): all 20 seeds STRUCTURALLY block --
    the order cannot furnish six disjoint chains of >= 25 ticks (short chains,
    wide antichains), a Kleitman-Rothschild-type non-manifoldlike signature.
  - Higher p (0.010-0.020; density 0.46-0.73): chains exist (16-20 valid) but
    the geometry gate fails -- held-out violation ~0.23-0.25 (>> 0.10) and null
    gap ~0 (+0.019, -0.005, +0.004; << 0.10). The profiles are no better than
    their column-shuffle: no consistent latent geometric structure.
- The verdict is interpretable because the SAME order-intrinsic discriminator
  passes on sprinkled geometry and recovers true position (P3-A: structured
  held-out 0.024-0.065, truth ~0.22 << 0.5 chance, separation d = 10.4). The
  block on the dynamics is therefore a property of the dynamics, not an
  instrument artifact.

Conclusion: transitive percolation (classical sequential growth), judged by the
validated discriminator at this scale, does not produce low-dimensional
geometry -- geometry does not emerge from this geometry-free dynamics. This is
consistent with the known non-manifoldlikeness of generic classical sequential
growth. Per Section 7 it is a controlled negative result about this specific
dynamics, not a statement about all possible dynamics; the natural next step is
dynamics families designed or tuned toward manifoldlikeness (e.g., non-trivial
CSG couplings), evaluated by the same frozen discriminator.
