# Gram Fourier challenge: population signal and unresolved instrument

2026-09-24. The new candidate uses an identity that is pointwise zero in
the effective background, reconstructed from order probabilities relative
to three anchors. Its oscillatory residual in the exact vacuum wave has a
nonzero complex population coefficient. This is a useful alternative to
the previous coordinate-oracle residual, whose linear conditional average
vanished. It is not yet a demonstrated finite-order detector or a sharp law.

## Frozen population calculation

Source and [design](gram_fourier_challenge.md) were frozen at
`3eb68f452256d7777b266ae868b9545a6502d5df` before the clean-tree CPU run.
The run retained all eight cases, q=1, phase 0.37. Here the displayed complex
quantity is S/(lambda^2 E w), with weights and Fourier phases applied to
the **reconstructed** states. Generator coordinates choose the fixed
population quadrature and anchors; they are not a finite-order observation.

| r | Anchor rotation | Complex coefficient | Relative difference from limit |
|---:|---:|---:|---:|
| 64 | 0 | -0.181459 - 0.038384i | 10.6552 |
| 64 | pi/8 | -0.178447 - 0.053546i | 10.5473 |
| 128 | 0 | -0.022764 - 0.068843i | 4.53430 |
| 128 | pi/8 | -0.021970 - 0.083450i | 4.53190 |
| 256 | 0 | 0.015248 + 0.005876i | 0.163440 |
| 256 | pi/8 | 0.016029 - 0.008762i | 0.162631 |
| 512 | 0 | 0.017393 + 0.008220i | 0.00533311 |
| 512 | pi/8 | 0.018164 - 0.006427i | 0.00564849 |

The predicted limits are respectively `0.01729433359+0.00824624955i` and
`0.01805991954-0.00639796158i`. This supports the predicted phase rotation
and nonzero coefficient in the computed high-frequency cases. Low-frequency
disagreements are large and remain in the record; the calculation is not
eight independent successes, classification accuracy or a confidence interval.

For r=512,theta=0 the predeclared tighter primitive and doubled quadrature
changed the normalized complex signal by `6.0475974e-7` relatively, below
the frozen 1% numerical-consistency criterion. The **weight mass**, however,
changed from `4.33143524e-5` to `4.39769377e-5`, about 1.53%. Thus the
normalized coefficient is much better converged than the separate mass;
do not assign the score's tolerance to the unnormalized integral or its mass.
There is no interval-arithmetic enclosure or independent-author computation
for these geometry results. All active anchor/probe causal slacks were positive.

The background algebra control recovered the Gram products with maximum
absolute error `1.2421e-15`; it also checked the zero residual and the four-angle
quadratic average. The run took 206.57 seconds on CPU, including its refinement.
No GPU was used. Its [manifest](gram_fourier_evidence/manifest.json) and
[raw rows](gram_fourier_evidence/rows.jsonl) are preserved byte for byte.
Raw SHA256:
`5c018c32903926b8055710f3769fa225b905540cd80b2fe8ac3c9d20c1849802`.

## What the new noise argument would change

The subsequent [proof draft](gram_fourier_noise.md) expands the empirical
raw moments before calibration, retains repeated reference indices, and
uses bounded variation along background null generators. It proposes

`noise <= C_K[1/(omega sqrt(n)) + 1/n + omega^K/n^((K+1)/2)]`

on the specified frequency range. A canonical-slice coupling handles wave
noise, including tip contributions. Disjoint pilot selection and a block
majority scan remove generator-chosen anchors, frequency and unknown
asymptotic threshold constants at the level of a mathematical construction.
The resulting candidate would separate each fixed power alpha>2 by choosing
a sufficiently large **fixed** K. It is not a uniform growing-K theorem.

The complete implementation with fixed smooth regularization and certified
increasing precision is still required by the challenge. No finite-data
power campaign or independent proof audit has been completed. The programme's
established bracket therefore remains **below 4/3 insufficient; power 3
sufficient under its stated logarithmic condition**. This candidate does not
prove that exponent 2 is necessary and leaves the original sharp-boundary
goal active.

The next decisive task is the order-only instrument and its error analysis,
followed by attacking the remaining lower/upper gap. More population rows by
themselves would not close that gap. The ordinary use of distances and Gram
matrices is established prior art; the draft records primary references and
does not treat search non-retrieval as novelty certification.
