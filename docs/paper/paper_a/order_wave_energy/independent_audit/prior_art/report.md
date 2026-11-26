# Prior-art search: order-only resolution law for a high-frequency gravitational plane wave

Search date: 2026-11-24. Produced by a search sub-agent that was given only a
description of the claimed bracket. It returned this report as text; the
audit orchestrator saved it unchanged apart from this header. The orchestrator
did not re-fetch the sources. The two preprint scope statements in Section 2
agree with what an earlier session of this programme recorded after reading
the same preprints.

Reading conventions: "read" = full text fetched and summarised; "abstract" =
arXiv abstract page only; "snippet" = search snippets or secondary
descriptions. Theorem statements are paraphrases; verify against the PDFs
before citing. All fetches went through a summariser. It gave 1/28 for
Madsen's d = 4 exponent; the formula gives 1/32.

A search that fails to retrieve a result is not a proof of novelty. The
coverage limits are in Section 5.

## 0. Bottom line

1. No DIRECT item was found. Nothing retrieved states a threshold, minimax
   rate, or total-variation / Hellinger / KL indistinguishability law in N
   versus a wavelength for sprinkled causal orders. Nothing treats an
   Isaacson/Burnett regime in a causal-set setting.
2. Closest by subject is Ibrahim-Bombelli 2611.09056; closest by form are
   Bombelli gr-qc/0002053 and Mueller 2503.01719 Theorem 3; closest by
   technique is the random-geometric-graph versus Erdos-Renyi line.
3. Existing quantitative causal-set results (Madsen 2607.05840, 2610.36695)
   are upper bounds for a fixed smooth geometry, with geometry-dependent
   constants. They say nothing uniform along a family whose curvature grows
   like 1/lambda.
4. The standard lore is that structure below the mean spacing rho^(-1/d) is
   not encoded. In d = 4 a naive reading puts the threshold at alpha = 4. The
   claimed bracket lies well below that, i.e. it asserts sub-spacing
   resolution of a coherent wave. No retrieved paper states or refutes this.

## 1. Table of items

Classes: DIRECT, ADJACENT (ADJACENT* = closest), BACKGROUND.

### 1.1 Causal sets: distinguishing geometries, Hauptvermutung, closeness

| # | Citation | URL | Year | What it establishes (how read) | Class |
|---|---|---|---|---|---|
| A1 | N. Ibrahim, L. Bombelli, "The Weyl tensor for gravitational-wave causal sets: thickness of geodesic analogs and geodesic deviation" | https://arxiv.org/abs/2611.09056 | 2026 | Heuristic estimate that about 1e8 or more elements are needed to see one Weyl component of a linearised plane wave via a longest-chain geodesic-deviation analog. No lower bound. (read) | ADJACENT* |
| A2 | L. Bombelli, "Statistical Lorentzian geometry and the closeness of Lorentzian manifolds", J. Math. Phys. 41 (2000) 6944 | https://arxiv.org/abs/gr-qc/0002053 | 2000 | Pseudo-distance d_n = Bhattacharyya-type angle between n-point sprinkled-poset laws; monotone in n; computed only for n = 2, 3 on a flat 2D rectangle versus cylinder. No scaling of n with feature size. (read) | ADJACENT* |
| A3 | O. Mueller, "On the Hauptvermutung of causal set theory" (v4) | https://arxiv.org/abs/2503.01719 | 2025-26 | Theorem 3: for every K, eps, D and unit-volume Cauchy slab X there is Y at Lorentzian GH distance > D whose K-point order law is within eps in TV. Theorem 5: same with flat cylinders. Theorem 6: positive result for another metric, K > 8s/delta. Theorem 2: countable version true. (read in full) | ADJACENT* |
| A4 | M. Braun, "Spacetime reconstruction by order and number" | https://arxiv.org/abs/2507.01907 | 2025 | Chronological adjacency matrices of i.i.d. samples have the same law iff the spacetimes are smoothly isometric. Qualitative. (abstract) | ADJACENT |
| A5 | N. Madsen, "On the uniqueness of embeddings of causal sets" | https://arxiv.org/abs/2607.05840 | 2026 | Two well-conditioned embeddings are related by an explicit eps-approximate isometry; Poisson sprinklings qualify with probability 1 - (rho V)^(-K'). Upper bound only. (first 100k chars read) | ADJACENT |
| A6 | N. Madsen, "Reconstruction of the Riemann tensor from causal order and element counts" | https://arxiv.org/abs/2610.36695 | 2026 | See Section 2.1. (read) | ADJACENT |
| A7 | S. Surya, "A closeness function on coarse grained Lorentzian geometries" | https://arxiv.org/abs/2510.19403 | 2025 | L1 distance between interval-abundance spectra; numerics at n = 20000; notes finite-n non-separable pairs and that an internal factor of volume below vol/n is invisible. No law in n versus curvature scale. (mostly read) | ADJACENT |
| A8 | L. Bombelli, J. Noldus, CQG 21 (2004) 4429; J. Noldus, gr-qc/0308074 | https://arxiv.org/abs/gr-qc/0402049 | 2003-04 | Lorentzian Gromov-Hausdorff distance. No sampling statement. (abstract) | BACKGROUND |
| A9 | F. Dowker, J. Butterfield, "Recovering general relativity from a Planck scale discrete theory of quantum gravity" | https://arxiv.org/abs/2106.01297 | 2021 | Sec. 4.1: a Lorentzian simplicial complex cannot tell Minkowski from a plane-fronted TT wave burst hidden in a vertex-free slab. Concerns lattices, not sprinklings; no count or wavelength law. (read) | ADJACENT |
| A10 | C. A. Hu, D. A. Meyer, E. J. Q. Meyer, "Reconstructing Minkowski geometry from causal separations" | https://arxiv.org/abs/2601.03280 | 2025 | Timelike separations determine spacelike ones. Qualitative. (abstract) | BACKGROUND |
| A11 | S. Major, D. Rideout, S. Surya, "On recovering continuum topology from a causal set" | https://arxiv.org/abs/gr-qc/0604124 | 2006 | Homology from thickened antichains at high enough density; threshold not quantified. (abstract) | BACKGROUND |

### 1.2 Causal sets: curvature estimators, fluctuations, manifoldlikeness

| # | Citation | URL | Year | What it establishes (how read) | Class |
|---|---|---|---|---|---|
| B1 | D. Benincasa, F. Dowker, "The scalar curvature of a causal set" | https://arxiv.org/abs/1001.2725 | 2010 | Nonlocal d'Alembertian; fluctuations grow with N unless a mesoscale is introduced. No required-density law. (snippet) | ADJACENT |
| B2 | Dowker-Glaser 1305.2588; Aslanbeigi-Saravani-Sorkin 1403.1622; Belenchia-Benincasa-Dowker 1510.04656 | https://arxiv.org/abs/1305.2588 | 2013-15 | Generalised operators; continuum limit of the mean. (abstract/snippet) | BACKGROUND |
| B3 | H. Moradi, Y. Yazdi, M. Zilhao, "Fluctuations and correlations in causal set theory" | https://arxiv.org/abs/2407.03395 | 2024 | Variance of the action decomposed into interval correlations. No density-to-see-curvature conclusion. (abstract) | ADJACENT |
| B4 | L. Machet, J. Wang, continuum limit of the BDG action | https://arxiv.org/abs/2007.13192 | 2020 | Mean only. (abstract) | BACKGROUND |
| B5 | M. Roy, D. Sinha, S. Surya, "The discrete geometry of a small causal diamond" | https://arxiv.org/abs/1212.0631 | 2012 | Curvature from k-chain abundances; Sec. 6 gives estimator error scaling rho^(-1/2), no constants, no numerics. No lower bound. (read) | ADJACENT |
| B6 | L. Glaser, S. Surya, "Towards a definition of locality in a manifoldlike causal set" | https://arxiv.org/abs/1309.3403 | 2013 | Interval-abundance flatness test; no power analysis. (snippet) | ADJACENT |
| B7 | J. Kambor, Nomaan X, "Manifold properties from causal sets using chains" | https://arxiv.org/abs/2007.03835 | 2020 | Chain-based estimators, numerics. (snippet) | ADJACENT |
| B8 | de Brito, Eichhorn, Pfeiffer, higher-order curvature operators | https://arxiv.org/abs/2301.13525 | 2023 | Higher curvature invariants. (abstract line) | BACKGROUND |
| B9 | Eichhorn, Mack, Le, Wagner, "Charting causal set configuration space with graph observables" | https://arxiv.org/abs/2605.27514 | 2026 | Graph observables separate nine classes; notes curvature invariants fluctuate strongly. (abstract) | ADJACENT |
| B10 | Barton, Borza, Roehrig, "Ollivier-Ricci curvature for causal sets" | https://arxiv.org/abs/2606.04910 | 2026 | Transport-based curvature; v2 scaling remark not read. (abstract) | ADJACENT |
| B11 | S. Shuman, "The geometry of tangent spaces on causal sets" | https://arxiv.org/abs/2405.11805 | 2024 | Estimators approach flat values with density; no rate. (abstract) | BACKGROUND |
| B12 | D. Reid, dimension tests in conformally flat spacetimes | https://arxiv.org/abs/gr-qc/0207103 | 2002 | (snippet) | BACKGROUND |
| B13 | Bolognesi-Lamb 1407.1649; Aghili-Bombelli-Pilgrim 1807.08701 | https://arxiv.org/abs/1407.1649 | 2014-18 | Manifoldlikeness indicators; no test theory. (snippet) | BACKGROUND |
| B14 | M. Boguna, D. Krioukov, "Measuring spatial distances in causal sets via causal overlaps" | https://arxiv.org/abs/2401.17376 | 2024 | Spacelike distance from causal overlap; scaling argument for relative overlap error (t_P/tau_c)^((d+2)/2). Flat only; no lower bound. (read) | ADJACENT |
| B15 | M. Boguna, D. Krioukov, "Sampling biases in sprinkled causal sets" | https://arxiv.org/abs/2610.03247 | 2026 | Finite-region bias of link observables. (abstract) | BACKGROUND |
| B16 | Rideout-Wallden 0810.1768; Eichhorn-Surya-Versteegen 1809.06192; Saravani-Aslanbeigi 1403.6429 | https://arxiv.org/abs/0810.1768 | 2008-18 | Spatial distance from order. (abstract/snippet) | BACKGROUND |
| B17 | Brightwell-Henson-Surya 0706.0375; Brightwell-Luczak 1510.05612 | https://arxiv.org/abs/1510.05612 | 2007-15 | A flat 2D sprinkling is exactly a random 2-dimensional order; a conformal factor correlates the null coordinates. (snippet) | BACKGROUND |

### 1.3 Plane waves and the high-frequency limit (continuum)

| # | Citation | URL | Year | Content | Class |
|---|---|---|---|---|---|
| C1 | R. Penrose, Rev. Mod. Phys. 37 (1965) 215 | (journal) | 1965 | Plane waves are not globally hyperbolic. (prior knowledge) | BACKGROUND |
| C2 | Marolf-Ross hep-th/0208197; Hubeny-Rangamani hep-th/0211195; Flores-Sanchez gr-qc/0211086 | https://arxiv.org/abs/hep-th/0211195 | 2002 | Causal structure of plane / pp-waves. No sampling content. | BACKGROUND |
| C3 | Burnett, J. Math. Phys. 30 (1989) 90; Huneau-Luk 1706.09501, survey 2404.07659; Luk-Rodnianski 2009.08968; Touati 2402.17530 | https://arxiv.org/abs/2404.07659 | 1989-2024 | High-frequency limits of vacuum spacetimes are null-dust solutions and conversely. Continuum only. (snippet) | BACKGROUND |

No item combining a causal set with a pp-wave or exact plane wave other than
A1 was retrieved; A1 lists exact pp-waves as future work.

### 1.4 Statistics / probability

| # | Citation | URL | Year | What it establishes (how read) | Class |
|---|---|---|---|---|---|
| D1 | Bubeck, Ding, Eldan, Racz, "Testing for high-dimensional geometry in random graphs" | https://arxiv.org/abs/1411.5713 | 2014 | Hard-threshold spherical RGG versus G(n,p), dense: threshold d ~ n^3; signed triangles; Wishart versus GOE. (snippet) | ADJACENT* |
| D2 | Bubeck-Ganguly 1509.03258; Racz-Richey 1611.05838; Brennan-Bresler-Huang 2103.14011 | https://arxiv.org/abs/1611.05838 | 2015-21 | Other proofs of n^3; limiting TV in the critical window; masked Wishart. (snippet) | ADJACENT |
| D3 | Brennan, Bresler, Nagaraj, "Phase transitions for detecting latent geometry in random graphs" | https://arxiv.org/abs/1910.14167 | 2019 | Sparse regime via couplings and information inequalities. (snippet) | ADJACENT |
| D4 | Liu, Mohanty, Schramm, Yang, "Testing thresholds for high-dimensional sparse random geometric graphs" | https://arxiv.org/abs/2111.11316 | 2021 | Sparse hard-threshold RGG indistinguishable once d >= polylog n; belief propagation / cavity for the posterior. (abstract) | ADJACENT* |
| D5 | Du, Mao, Sun, Wu, Xu, "Sharp detection thresholds for random geometric graphs with polynomial average degrees" | https://arxiv.org/abs/2610.08517 | 2026 | KL -> 0 for d >> (np log(1/p))^3; truncated chi-square for the last vertex, Fourier expansion of the conditional likelihood ratio, cavity recurrences. (abstract) | ADJACENT* |
| D6 | Mao, Wu, Xu, smooth kernels and a spectral conjecture | https://arxiv.org/abs/2602.14998 | 2026 | Critical d = n^(3/4); conjectured criterion n^3 tr^2(kappa^3) = 1; posterior-overlap method. (abstract) | ADJACENT |
| D7 | Huang, Li, "Detecting random geometric graphs with general kernel" | https://arxiv.org/abs/2610.12765 | 2026 | Thresholds for non-differentiable kernels; signed quadrilaterals can beat triangles. (abstract) | ADJACENT |
| D8 | Bangachev-Bresler 2310.14501, 2402.12589, 2305.04802, 2408.00995 | https://arxiv.org/abs/2310.14501 | 2023-24 | Signed 4-cycles optimal among low-degree tests for L_infinity RGG; cluster expansion; Fourier coefficients of hard-threshold models; sandwich coupling. (abstracts) | ADJACENT* |
| D9 | Liu-Racz 2110.15886; Duchemin-de Castro survey 2203.15351 | https://arxiv.org/abs/2110.15886 | 2021-22 | Soft RGG thresholds; survey. | BACKGROUND |
| D10 | Mao-Wein-Zhang 2402.00305; Cai-Liang-Rakhlin 1604.06474 | https://arxiv.org/abs/2402.00305 | 2016-24 | Thresholds for a planted one-dimensional geometry; noisy relation, product null. (snippet) | ADJACENT |
| D11 | Kunisky, Wein, Bandeira, low-degree likelihood ratio notes | https://arxiv.org/abs/1907.11636 | 2019 | The low-degree method. | BACKGROUND |
| D12 | Arias-Castro, Channarond, Pelletier, Verzelen, "On the estimation of latent distances using graph distances" | https://arxiv.org/abs/1804.10611 | 2018 | Indicator link: matching lower bound proved by two latent configurations with the same graph. (read) | ADJACENT* |
| D13 | Dani-Diaz-Hayes-Moore 2107.14323; Huang-Jiradilok-Mossel 2511.05434 | https://arxiv.org/abs/2511.05434 | 2021-25 | Position / distance reconstruction from the graph; volumetric lower bounds. (abstract/snippet) | ADJACENT |
| D14 | Kleindessner-von Luxburg COLT 2014; Arias-Castro 1501.02861; Jain-Jamieson-Nowak 1606.07081 | https://arxiv.org/abs/1501.02861 | 2014-17 | Ordinal embedding consistency and first rates; optimal rates open. (snippet) | ADJACENT |
| D15 | Janson 0902.0306; Hladky-Mathe-Patel-Pikhurko 1211.2473 | https://arxiv.org/abs/0902.0306 | 2009-12 | Kernel representation of exchangeable random posets. No testing rates. (snippet) | BACKGROUND |
| D16 | Winkler, Order 1 (1985); Bollobas-Brightwell, TAMS 324 (1991) 59 | https://doi.org/10.1090/S0002-9947-1991-0986685-9 | 1985-91 | Random k-dimensional orders, box-spaces. No two-model testing result found. (snippet) | BACKGROUND |
| D17 | Gruebel 2206.12153; Kim-Balakrishnan-Wasserman 2003.13208; Berrett-Kontoyiannis-Samworth 2001.05513 | https://arxiv.org/abs/2206.12153 | 2020-22 | Permutons; minimax independence testing. Relevant to the 1+1 analog only. (snippet) | ADJACENT |
| D18 | Ingster-Suslina (2003); Arias-Castro-Pelletier-Saligrama 1607.08156 | https://arxiv.org/abs/1607.08156 | 2003-16 | Minimax detection of density perturbations, coordinates observed. (snippet) | BACKGROUND |
| D19 | Ghoshdastidar et al. 1707.00833, 1705.06168; Maugis et al. 1701.00505 | https://arxiv.org/abs/1707.00833 | 2017 | Two-sample graph tests with subgraph counts; noisy edges. (snippet) | ADJACENT |
| D20 | Aamari, Levrard, curvature estimation rates | https://arxiv.org/abs/1705.00989 | 2017 | Minimax n^(-(k-2)/d) for curvature from points with coordinates. (snippet) | ADJACENT |
| D21 | van der Hoorn et al. 2008.01209; Lubold-Chandrasekhar-McCormick 2012.10559 | https://arxiv.org/abs/2008.01209 | 2020 | Graph curvature converges to Ricci curvature; tests of latent curvature type. Consistency only. (snippet) | ADJACENT |

### 1.5 Works the group already cites (known)

| # | Citation | URL | Threshold or lower bound of the kind in question? |
|---|---|---|---|
| E1 | Henson, gr-qc/0601069 | https://arxiv.org/abs/gr-qc/0601069 | No. 2D, N = 10000; residual null-coordinate discrepancy about 0.4-0.8 fundamental units; qualitative remark that many sprinklings give the same causal set. 20 citers (Semantic Scholar) screened by title: none. (read) |
| E2 | Johnston, 2502.09701 (and 2111.09331) | https://arxiv.org/abs/2502.09701 | No. Simulation metrics for d = 2, 3, 4 at N = 500-1500; no error law. One INSPIRE citation (by elimination, Madsen 2610.36695). (read) |
| E3 | Nenashev-Baranovskii, 2302.12209 (Part II 2312.03487) | https://arxiv.org/abs/2302.12209 | No. Sixteen causal-relation measurements certify non-conformal-flatness in 4D. No random events or sample-size statement in the abstract; treating it as a deterministic chosen configuration is the searcher's inference from the Part II description. One INSPIRE citation (Part II). (abstract) |
| E4 | Gine-Latala-Zinn, math/0003228 | https://arxiv.org/abs/math/0003228 | No. Probability tool; citers not screened. (abstract) |

## 2. Scope statements

### 2.1 Madsen, arXiv:2610.36695 (27 Oct 2026, 55 pp.)

- **Setting.** Fixed smooth globally hyperbolic Lorentzian d-manifold
  (d = spacetime dimension, d >= 2). Data = causal order plus counts in a
  Poisson sample of intensity rho.
- **Structure.** Sec. 2 geometry and counting (2.2 Poisson realization, 2.6
  scales); Sec. 3 normal coordinates from order and interval counts; Sec. 4
  transport, holonomy, curvature (4.7 the reconstruction theorem); Sec. 5
  alternative interval estimators.
- **Deterministic input.** Definition 2.1 (V1, V2): relative count error in
  diamonds of height s at most K sqrt(L_rho/(rho s^d)), plus a global count
  condition. Theorem 2.4: Poisson samples satisfy these with failure
  probability at most rho^(-N). Corollary 2.5: almost surely for large rho.
- **Intermediate.** Theorem 2.16: proper time from counts, relative error
  C (l^2 + delta_t + t/l).
- **Main result.** Theorem 4.9 (deterministic) and Corollary 4.10
  (probability >= 1 - rho^(-N), eventually almost surely): with
  r = (L_rho/rho)^(1/(3d+4)) and H = sqrt(r), uniform Riemann error at most
  C (L_rho/rho)^(1/[2(3d+4)]). The Conclusion uses L_rho = log(e+rho).
  Exponent: 1/20 (d=2), 1/26 (d=3), 1/32 (d=4). Corollary 4.11: same rate for
  scalar curvature. Theorem 5.2: same bound for general sampled charts.
- **Constants.** Depend on the fixed geometry, a volume threshold, the count
  allowance and construction thresholds. No explicit dependence on curvature
  derivatives was found.
- **Lower bound / threshold in N versus wavelength / high-frequency regime:**
  none. The Conclusion only says sharper rates may be sought.
- **Bearing on the bracket.** None directly. The searcher's extrapolation,
  not in the paper: requiring H below the wavelength needs
  N >> lambda^(-2(3d+4)), i.e. alpha > 32 in d = 4.
- **Not seen verbatim.** Definition 2.1, the definition of L_rho,
  Proposition 3.1, Appendix A.

### 2.2 Ibrahim and Bombelli, arXiv:2611.09056 (6 Nov 2026, 24 pp.)

- **Spacetime.** Linearised plus-polarised TT plane wave,
  h_xx = h0 cos[omega(t - z)]. The causal-relation formula is exact for that
  metric, which solves the vacuum equations only to first order in h0.
  h0 = 1e-2.
- **Causal relations (Sec. 4.1).** Sign of a conserved geodesic tangent norm
  (twice Synge's world function), closed forms in Appendix A, valid within a
  convex neighbourhood of size about (1/omega) sqrt(2/h0).
- **What is simulated.** Only Minkowski sprinklings in d = 2, 3, 4 (up to
  about 4e6 elements). Thickness/H = A_d N^(-beta_d), beta about 0.149,
  0.143, 0.111, all below 1/d (Sec. 3).
- **Feasibility (Sec. 6, Appendix B).** Uncertainty of the second derivative
  of tube separation must be below (1/2) omega^2 h0; three-point stencil with
  step lambda/4; ideal distance estimator assumed. Result:
  thickness <= (pi^2/16)(h0/a_s) L and
  N_total ~ 7.3e9 (a_s/a_h)^7 a_w^2 a_L^(-6), about 1e8-1e9. Probe height
  lambda/2, separation a_L lambda, needs L << lambda.
- **Lower bound / threshold / high-frequency regime:** none. The authors
  state that "our result does not say that there is no other way" (Ibrahim
  and Bombelli). Only the long-wavelength regime is discussed.
- **Bearing on the bracket.** A requirement for one estimator at wavelength
  comparable to the probed region. The searcher's extrapolation, not in the
  paper: N_local ~ h0^(-1/beta_4) ~ h0^(-9); with h0 ~ lambda and
  N_local ~ N lambda^4 this is alpha of order 13.
- **Citations.** It cites Madsen 2610.36695; INSPIRE lists none of it yet.

## 3. The abstract statistical problem

### 3.1 Special case of a solved problem? Not as far as found

- **Abstract form.** N i.i.d. latent points, a deterministic 0/1 kernel, two
  (kernel, measure) pairs differing by a small high-frequency perturbation.
  In Janson's language this is a test between two kernel-random posets with
  random-free kernels. No testing or TV theory for poset limits was found.
- **RGG thresholds (D1-D8).** These test a geometric model against a product
  measure with dimension as the parameter. Here both hypotheses are geometric
  and the dimension is fixed.
- **Graphon two-sample results (D19).** These assume noisy, conditionally
  independent edges. A deterministic relation makes the conditional
  likelihood degenerate.
- **One genuine reduction (searcher's observation, supported by B17).** In
  1+1 dimensions the sprinkled order is exactly the rank permutation of a
  sample from the copula of the volume measure in null coordinates. The 2D
  analog is therefore a rank test between copulas (D17). It does not cover
  the plane wave, since Weyl vanishes in 2D. No rank-test result at a
  growing-frequency alternative matching the exponents was located.

### 3.2 Techniques most likely to transfer

1. **Vertex-exposure KL chain rule with posterior/cavity analysis and
   truncated chi-square.**
   - References: 1910.14167; 2111.11316; 2610.08517; 2602.14998 for the
     two-replica variant.
   - For: it is the only method with sharp exponents for deterministic
     hard-threshold relations. It needs posterior localisation of latent
     points given the relation, the same quantity that controls a statistic
     built from order-reconstructed coordinates. The 2026 device of
     Fourier-expanding the conditional likelihood ratio and truncating high
     frequencies is suggestive.
   - Against: all worked cases have a product null. Here one must compare two
     posteriors; no example was found.
2. **Signed motif counts, kernel spectrum, low-degree method.**
   - References: 1411.5713; 2310.14501; 2402.12589; 2610.12765; 1907.11636.
   - For: centred sub-poset counts are the analog and are U-statistics of the
     kind Gine-Latala-Zinn controls. They are indifferent to determinism. The
     lesson that the optimal motif depends on kernel regularity bears on
     choosing the order observable.
   - Against: as a lower bound it only excludes low-degree tests and needs a
     product null for the orthogonal basis. Fixed small motifs average over
     the wave's phase, so this route may give sufficiency bounds that are not
     sharp.
3. **Coupling through the latent configuration (re-embedding, sandwiching,
   collisions).**
   - References: 1411.5713 with 1611.05838; 2408.00995; 1804.10611 Theorem 2.
   - For: native to deterministic relations. TV between order laws is bounded
     by TV between latent laws after an almost order-preserving re-embedding,
     which is how a bound better than the coordinates-revealed exponent can
     arise. Mueller's Theorem 3 and Henson's remark are the same mechanism.
   - Against: data-processing bounds are lossy; in the RGG line this route is
     sharp only in the dense regime. There is no explicit sufficient
     statistic playing the role of the Gram matrix.

## 4. Lore the claim would sharpen

- Bombelli describes the finite-n distance as comparing manifolds only down
  to a volume scale V/n. Surya says the spectrum characterises a spacetime up
  to the discreteness scale. Madsen 2607.05840 labels rho^(-1/d) the smallest
  meaningful scale.
- None of these is a detectability theorem, and none treats a coherent
  oscillation with wavelength below the mean spacing extending over the whole
  region.
- Boguna-Krioukov and Henson indicate that order pins down null coordinates
  more finely than the mean spacing in a favourable frame. Neither states a
  law.

## 5. Queries run and access limits

**Fetched.**
- Abstract pages: 2610.36695, 2611.09056, 2607.05840, 2503.01719,
  2510.19403, 2507.01907, 2302.12209, 2407.03395, 2106.01297, 2605.27514,
  2405.11805, 2606.04910, 2610.03247, 2610.08517, 2610.12765, 2602.14998,
  2408.00995, 2511.05434, 2111.11316, 2305.04802, 2310.14501.
- Full text: 2610.36695 (three windows), 2611.09056 (both windows),
  2607.05840 (first 100k chars), 2503.01719 (complete PDF), 2510.19403 (two
  windows), gr-qc/0002053, gr-qc/0601069, 2502.09701, 1212.0631, 1804.10611,
  2401.17376, 2106.01297.

**arXiv API (abstract-field boolean search).**
- "causal set" AND (plane wave / pp-wave / gravitational wave): 5 hits, only
  2611.09056 on point.
- "causal set(s)" AND (total variation / minimax / indistinguishable /
  hypothesis test / lower bound / sample complexity / convergence rate /
  error bound): 5 hits, only 2610.36695 relevant.
- "causal set(s)" AND (Weyl / Riemann / geodesic deviation / curvature): 34
  hits, all screened.
- Hauptvermutung AND (causal / Lorentzian): 5 hits.
- ("causal set(s)" / sprinkling) AND (wavelength / high-frequency / Fourier /
  resolution / Isaacson): 12 hits, none relevant.
- (random geometric graph / latent geometry / detection threshold) AND
  (Lorentzian / partial order / causal set / Minkowski / de Sitter /
  spacetime): 15 hits, none a threshold result.
- Random partial order / random orders / poset limits / random poset: 639
  hits, polluted by streaming papers; first 43 screened, rest not.

**Web searches (paraphrased topics; no project or personal data).**
- Causal sets and waves: plane-wave / pp-wave sprinkling; high-frequency /
  Isaacson / Burnett; whether a causal set can detect a wave shorter than the
  discreteness scale.
- Causal-set thresholds: density to detect curvature; BDG variance and
  mesoscale; quantitative Hauptvermutung with total variation;
  manifoldlikeness test power; Myrheim-Meyer variance; Roy-Sinha-Surya
  errors; sub-discreteness structure; Fisher information / Cramer-Rao;
  Bombelli closeness; Bombelli-Noldus; Braun.
- Known background: Henson; Johnston; Nenashev-Baranovskii;
  Gine-Latala-Zinn; Boguna-Krioukov overlaps.
- Statistics: RGG versus Erdos-Renyi; Bangachev-Bresler;
  Brennan-Bresler-Nagaraj and Liu-Racz; Wishart-GOE; low-degree method; RGG
  position recovery; ordinal embedding; box orders and random orders; Janson
  poset limits; two-sample graphon tests; testing two random-free graphons;
  permuton / copula rank tests; 2D orders; Ollivier-Ricci convergence; latent
  curvature tests; Aamari-Levrard; planted dense cycles; Ingster-Suslina;
  math.ST / math.PR Lorentzian random orders.
- Continuum: Burnett conjecture; plane-wave causal structure.

**Citation lookups.**
- INSPIRE by record id: 2302.12209 (1), 2610.36695 (1), 2502.09701 (1),
  2611.09056 (0). The citing records are 2312.03487, 2611.09056, 2610.36695.
- Semantic Scholar: 20 citers of gr-qc/0601069; empty for 2302.12209; HTTP
  429 for 2502.09701.

**Could not access or did not do.**
- Google Scholar directly.
- INSPIRE `refersto:arxiv` queries (not honoured).
- The arXiv web search form with booleans (returned zero for known
  positives).
- Body-text search of any kind.
- Local PDF page rendering (long PDFs were read via HTML/ar5iv).
- Full text of 2607.05840 beyond the first part, 2507.01907, 2407.03395,
  2606.04910 v2, 2605.27514, and all statistics papers except 1804.10611.
- Paywalled items (Penrose 1965, Burnett 1989, Winkler 1985,
  Bollobas-Brightwell 1991, Kleindessner-von Luxburg 2014, Ingster-Suslina
  2003).
- Citers of Gine-Latala-Zinn.
- Theses, talks and non-arXiv venues systematically.
- Non-English literature.
