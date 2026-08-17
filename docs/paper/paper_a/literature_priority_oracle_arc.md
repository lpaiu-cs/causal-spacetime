# Literature priority assessment — oracle arc and volume-count verification

Search date: 2026-08-17.
Scope: claims in the oracle arc (S3 volume certification → O3/O3′ endpoints
→ O4b instrument audit → O5 Poisson-count campaign) that could overlap
existing causal-set or mathematical-relativity literature.

## Method

Searched arXiv (gr-qc, math-ph, hep-th), INSPIRE-HEP, Zenodo, and Google
Scholar for: causal-set volume certification, diamond volume enclosure,
Schwarzschild causal relation, Poisson sprinkle verification, causal-set
Hauptvermutung, and related terms. Directly verified PDFs / TeX sources for
ambiguous cases (Alfyorov 2026, Braun 2025).

## Claim clusters

### 1. Causal-order-and-number/volume reconstruction

**Status: PRIOR ART EXISTS**

The Hauptvermutung — that a Poisson sprinkling into a Lorentzian manifold
determines the manifold up to isometry — is a foundational conjecture of the
causal-set program (Bombelli et al. 1987; Myrheim 1978). Recent progress:

- **Braun 2025** (arXiv:2507.01907): proves that i.i.d. samples together
  with the chronological order determine a smooth isometry (weakened variant
  of the Hauptvermutung). Uses "causal-order-and-number" reconstruction, not
  "order-only".

Our program does not claim to prove or extend the Hauptvermutung. The oracle
arc verifies that a *specific* certified volume interval is consistent with
a *specific* Poisson count — a measurement, not a reconstruction theorem.

### 2. Schwarzschild sprinkle count verification

**Status: PARTIAL OVERLAP — prior art for Poisson-sprinkle counting**

- **Homšak–Veroni 2024** (arXiv:2404.11670, "Boltzmannian state counting
  for black hole entropy in Causal Set Theory", §V.1): Poisson-sprinkled
  count vs volume on Schwarzschild causal sets, with a Poissonity test.

Our novelty: the preregistered TOST equivalence gate against an MPFR
directed-rounding volume enclosure whose endpoints were frozen before any
count data existed. The statistical framework (Garwood intervals, frozen
equivalence band, tri-state ambiguity handling) and the certification
pipeline are new.

### 3. Certified volume oracle / directed-rounding enclosure

**Status: NO DIRECT PRIOR**

- **Berthiere–Gibbons–Solodukhin 2016** (arXiv:1507.03619): analytic volume
  comparison theorems for causal diamonds in curved spacetime. These are
  asymptotic/analytic results, not finite-anchor MPFR enclosures with
  directed rounding.
- **Wang 2019** (arXiv:1904.01034, "Geometry of small causal diamonds",
  Phys. Rev. D 100:064020): causal-diamond volume expansion including
  Weyl curvature contribution. Again analytic, not a numerical
  certification.

We provide a reproducible directed-rounding enclosure (not "the first
certified numerical enclosure"): MPFR flight-time certification, uniform
anchor-diamond containment proof, and certified cell-refinement integrator.

### 4. Schwarzschild causal-relation algorithm

**Status: PARTIAL OVERLAP — He–Rideout is prior art**

- **He–Rideout 2009** (Song He and David Rideout, arXiv:0811.4235,
  CQG 26:125015): Schwarzschild causal-relation algorithm for causal-set
  sprinkling. Subsequent corrections to the inside-horizon branch exist
  in the literature.

Our implementation uses a different reduction (static-spacetime radar
coordinates) but He–Rideout is a mandatory citation for Schwarzschild
causal-relation work.

## Alfyorov–Shnyukov resolution

Flagged as a potential priority conflict due to concurrent causal-set work on
curved spacetimes.

- **Paper 7** ("Weyl curvature from the Hasse diagram", Zenodo v3.1,
  2026-03-31, doi:10.5281/zenodo.19364212): CJ estimator = stratified
  covariance of Hasse path-count observables. Bridge formula:
  ⟨CJ⟩ = C₀·N^{8/9}·E_{ij}E^{ij}·T⁴. Verified on pp-wave diamonds.
  105 Lean 4 theorems. Not on arXiv; indexed on INSPIRE (record 3142005).
  CQG submission inferred from cover letter files.

  **Assessment: orthogonal.** Paper 7 measures *curvature* (Weyl tensor
  E_{ij}E^{ij}) from Hasse diagrams; our program verifies *volume*
  (sprinkled count vs certified 4D measure). Different observables, methods,
  spacetimes, and claim structures. No priority conflict.

- **Paper 8** ("Non-perturbative spectral gravity measure", INSPIRE record
  3143309): functional measure theory for spectral action. Zero overlap
  with our work.

## Summary

| Cluster | Prior art | Our novelty scope |
|---------|-----------|-------------------|
| Order-and-number reconstruction | Braun 2025, Myrheim 1978, BLMS 1987 | Not claimed — we do measurement, not reconstruction |
| Schwarzschild sprinkle count | Homšak–Veroni §V.1 | Preregistered TOST gate against certified MPFR enclosure |
| Directed-rounding enclosure | None direct (Berthiere et al. analytic only) | Reproducible directed-rounding enclosure pipeline |
| Schwarzschild causal relation | He–Rideout 2009 | Different reduction; He–Rideout is mandatory citation |
| Curved-spacetime causal-set (Alfyorov) | Orthogonal (curvature, not volume) | No conflict |
