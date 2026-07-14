# P6: Chain-rich hard negatives and diagnostic comparison

Status: **STAGE A PLANNED, nothing frozen** (2026-07-14). Stage B is locked
until the Stage-A artifacts and `docs/prereg/frozen/p6_test_constants.json`
are committed together.

## 1. Question

Does the frozen P3 order-intrinsic discriminator reject a non-manifoldlike
order only after it has successfully extracted the required observer chains
and targets? The P5 crystal control did not answer this because all crystal
configurations blocked structurally before the numerical gates were reached.

## 2. Frozen inherited instrument

P6 does not recalibrate the discriminator. It inherits the P3/P5 protocol:
six disjoint chains of at least 25 elements, 20-44 bracketed targets,
order-intrinsic bracket profiles, one-dimensional ordinal embedding, and the
fixed gates heldout <= 0.10, null_gap >= 0.10, truth <= 0.40.

## 3. Stage-A constructions

All constructions use N = 600 and seeds 0-9.

### Partially layered family

A complete k-layer order is represented by a permutation formed from k
contiguous descending blocks. Layer sizes have a minimum of six and the
remaining elements are assigned by a seeded multinomial draw. Windowed
transpositions then soften the exact layers.

- layer count k in {25, 40, 60}
- transposition count in {20, 100, 500, 2000}
- position window 60

A cell is eligible for Stage B only if at least 8/10 seeds reach the numerical
gates and at least 8/10 are blocked by those gates. This is a calibration rule,
not a confirmatory result.

### Local-shuffle construction audit

The proposed P5-local-shuffle family requires an identifiability audit before
it can be treated as a hard negative. A permutation after any transpositions
still defines an exact 2D order with new coordinates `(t, x) = (i + pi_i,
i - pi_i)`. Stage A therefore applies {30, 150, 600, 2400} transpositions
within a position window of 32 to uniform N = 600 continuum representatives,
then scores one fitted embedding against both the original and current exact
coordinates.

A current-coordinate pass paired with an original-coordinate failure means
that the operation remapped geometry relative to an externally retained
label. It does **not** establish a geometry-free hard negative. Because this
conclusion follows from the representation itself, expensive regeneration of
the frozen P5 beta=2/8 samples is not scientifically informative unless this
audit finds contrary behavior.

## 4. Freeze and Stage B

After Stage A, eligible layered cells and their expected chain-rich block
verdicts are frozen. Confirmatory seeds are 100-119. A frozen cell confirms if
at least 16/20 seeds both reach the numerical gates and block there. Structural
blocks do not count as successes. No grid, inherited gate, or denominator may
move after the freeze.

The local-shuffle arm proceeds to Stage B only if Stage A demonstrates loss of
intrinsic current-coordinate geometry rather than coordinate remapping.

## 5. P6b comparison

Before computing confirmatory comparison metrics, freeze the input labels,
diagnostics, and metrics: P1 epsilon endpoints/sweep, P3 orders, P5 continuum
and crystal configurations, and confirmed P6 negatives; Myrheim-Meyer
dimension, interval-abundance profile distance, and height; ROC/AUC,
P1 truth-rank correlation, and overlap with the P1 H-LAG region. Optional
heavier diagnostics remain outside the primary comparison.

## 6. Deviations log

(empty)
