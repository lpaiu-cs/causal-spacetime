# Literature priority assessment — S6 mass ladder (curvature generalization)

Search date: 2026-08-19.
Scope: the claims specific to the S6 arc, which generalizes the oracle-arc
Poisson-count verification from a single Schwarzschild mass to a preregistered
mass ladder. This document is the companion to
`literature_priority_oracle_arc.md` (the single-mass O5 arc) and does not
restate its findings; it assesses only what the mass generalization adds.

The arc under review verifies, at three preregistered Schwarzschild masses
`M ∈ {1.0, 1.4, 1.8}` — a dimensionless curvature ladder `μ = 2M/r_c`,
`r_c = 15`, i.e. `μ ∈ {0.1333, 0.1867, 0.2400}` — that a Poisson sprinkling's
certified-membership count inside a diamond fixed in **absolute** coordinates
(shell `[10, 20]`, anchors `(12, 18)`) is consistent with an independently
certified continuum 4-volume `V(M)` from an MPFR directed-rounding enclosure
whose endpoints were frozen before any count data existed. The gate is
preregistered: a conservative outer count interval is formed from the exact
Garwood lower limit at `K_certain` and the upper limit at `K_certain + U_amb`
(the ambiguous points inflating only the upper end), each tail at level 0.025;
rescaled by the intensity to a volume interval `C = [L/A, U/A]`, its identified
discrepancy `D = [C_lo − V_hi, C_hi − V_lo]` against `V(M)` is required to lie
within a frozen `τ = 2.5%` band (`B = τ·V_ref`) — a tri-state test (contained →
CONCORDANT, disjoint → DISCORDANT, else INCONCLUSIVE). Because the diamond
is held in absolute coordinates while only `M` varies, no rung is an isometric
copy of another. Each rung was decided independently (all three CONCORDANT;
verdicts recorded in `docs/prereg/p14_s6_m14_count.json`,
`p14_s6_m18_count.json`, and — for the central rung — the O5 artifact).

## Method

Searched arXiv (gr-qc, math-ph, hep-th), INSPIRE-HEP, Zenodo, Research Square,
and Crossref for: mass/curvature-parameterized causal-set sprinkling,
number-vs-volume checks in curved and cosmological backgrounds, causal-set
curvature indicators and actions, equivalence / TOST testing for model and
count validation, exact Poisson interval methodology, and the full
Alfyorov–Shnyukov corpus (version history and siblings). Every candidate cited
below was confirmed by fetching its arXiv abstract page, an INSPIRE/Crossref
record, or a DOI resolver; bibliographic metadata is locked in
`citations/references.bib`.

## Claim clusters

### 1. Mass/curvature-parameterized count-vs-certified-volume across a ladder

**Status: NO DIRECT PRIOR (the mass-ladder framing is the arc's own).**

The nearest works each miss the generalization on a different axis:

- **Homšak–Veroni 2024** (arXiv:2404.11670 §V.1; catalogued in the oracle-arc
  document): a Poisson-sprinkled count vs volume in Schwarzschild with a
  Poissonity test — but at essentially one setting, against an analytic
  volume, with no preregistered equivalence band.
- **Roy–Sinha–Surya 2013** (arXiv:1212.0631, Phys. Rev. D 87:044046): the
  expected causal-set element/chain count inside a *small* causal diamond in
  curved spacetime, as a Riemann-normal-coordinate power series in the local
  curvature. This is the closest analytic ancestor of "compare a sprinkled
  count to a curved-diamond volume", but it is a single-setting expansion, its
  reference volume is an analytic series rather than a certified enclosure, and
  it runs no equivalence test on a realized Poisson count.
- **Barton–Counsell–Dowker–Gould–Jubb–Taylor 2019** (arXiv:1909.08620,
  Phys. Rev. D 100:126008): a discrete causal-set count ("horizon molecules")
  reproduces a continuum geometric quantity in a black-hole spacetime — but the
  target is horizon *area* (a 2-surface), not the certified 4-volume of a fixed
  diamond, and there is no mass sweep.

None of these varies curvature across a preregistered ladder of non-isometric
rungs, uses a directed-rounding certified volume as the reference, or applies a
frozen equivalence band with a tri-state verdict. The mass generalization —
novelty (a) — has no direct prior. Its novelty is *not* that a sprinkled count
can track a curved-spacetime volume (that is old; see cluster 3) but that the
tracking is gated, rung by rung, against a certified enclosure across a
curvature axis chosen a priori to make the rungs distinct.

### 2. The dimensionless curvature indicator μ = 2M/r as the ladder axis

**Status: NO OVERLAP — μ is a continuum label, not a discrete estimator.**

A reader who sees "curvature indicator" in a causal-set paper will expect the
causal-set curvature machinery, and it must be cited precisely to *disclaim*
that we propose a new one:

- **Benincasa–Dowker 2010** (arXiv:1001.2725, Phys. Rev. Lett. 104:181301):
  the canonical discrete scalar curvature / Benincasa–Dowker action, read off
  the sprinkled causal set by order-theoretic combinatorics converging to `R`.
- **Dowker–Glaser 2013** (arXiv:1305.2588, Class. Quantum Grav. 30:195016):
  the dimension-generic family of these operators.
- **Belenchia–Benincasa–Dowker 2016** (arXiv:1510.04656,
  Class. Quantum Grav. 33:245018): proof that the 4D causal-set d'Alembertian
  mean converges to `□ − R/2` in curved spacetime — the theorem that
  legitimizes reading curvature off a sprinkling at all.
- **Barton–Borza–Röhrig 2026** (arXiv:2606.04910): the newest such indicator,
  an Ollivier–Ricci (Lorentzian optimal transport) curvature recovered from
  high-density sprinklings on Minkowski/de Sitter/anti-de Sitter diamonds.

All four *estimate* a discrete curvature from the causal set. Our `μ = 2M/r_c`
is the opposite construct: a continuum geometric label naming which
Schwarzschild spacetime a rung inhabits, fixed a priori to guarantee
non-isometry and to order the rungs, and never estimated from the count. There
is no overlap of object or method; these are cited as the machinery we are
*not* using, so that "curvature indicator" is not misread as a discrete
curvature estimate.

### 3. Number-volume correspondence in curved / cosmological backgrounds

**Status: PRIOR ART for the correspondence; the certified gate is ours.**

That a sprinkled count tracks `ρV` in a genuinely curved region is old and
foundational — it is the mechanism of the causal-set cosmological-constant
program:

- **Sorkin 1997** (arXiv:gr-qc/9706002, Int. J. Theor. Phys. 36:2759):
  the number-volume fluctuation argument, `Λ ∼ 1/√V`.
- **Ahmed–Dodelson–Greene–Sorkin 2004** (arXiv:astro-ph/0209274,
  Phys. Rev. D 69:103523, "Everpresent Λ"): turns `N = ρV +` Poisson
  fluctuation in an expanding FRW spacetime into a concrete, testable
  cosmological model.

This line establishes `N ∼ ρV` in curved spacetime as settled physics, and it
is cited as the conceptual backdrop so the S6 correspondence is not overclaimed
as new. The distinction is method, not correspondence: Everpresent-Λ *exploits*
the number-volume fluctuation as physics (a fluctuating dark energy) against an
analytic FRW volume, whereas S6 *gates* a realized count against a
directed-rounding certified enclosure, rung by rung, with a preregistered
equivalence band. Novelty (a) is the certified gate across a curvature ladder,
not the underlying `N ∼ ρV` fact.

### 4. The equivalence-testing statistical framework

**Status: METHODOLOGY TRANSFER — the apparatus is standard; the application
is the novelty.**

The count stage's framework is assembled from established statistics and must
be cited as such, so that the abstract statistical machinery is not presented
as a from-scratch invention:

- **Schuirmann 1987** (J. Pharmacokinet. Biopharm. 15:657): the two one-sided
  tests (TOST) procedure. The S6 gate is a robust confidence-interval relative
  of it, each Garwood tail at level 0.025 (a 95% level — *more conservative*
  than the conventional 0.05 TOST, whose 90% interval it replaces), not a
  single-count symmetric CI: the outer count interval takes its lower limit at
  `K_certain` and its upper limit at `K_certain + U_amb`; rescaled to a volume
  interval `C`, the identified discrepancy `D = [C_lo − V_hi, C_hi − V_lo]` of
  `C` against the certified enclosure is required to lie within `±B` (tri-state:
  contained → CONCORDANT, disjoint → DISCORDANT, else INCONCLUSIVE). This
  describes the frozen rule; it does not change it.
- **Lakens 2017** (Soc. Psychol. Personal. Sci. 8:355) and
  **Lakens–Scheel–Isager 2018** (Adv. Methods Pract. Psychol. Sci. 1:259):
  the modern primer/tutorial, including the "smallest effect size of interest"
  (SESOI) framing that corresponds to fixing the tolerance band `τ` a priori.
- **Garwood 1936** (Biometrika 28:437): the exact ("Garwood") Poisson interval
  used to bound the certified-membership count; **Clopper–Pearson 1934**
  (Biometrika 26:404): the exact binomial interval used in the pilot and power
  certifications; **DeLong–DeLong–Clarke-Pearson 1988** (Biometrics 44:837):
  the nonparametric AUC interval used in the Type-D C2 stage.
- **Robinson–Froese 2004** (Ecol. Modelling 176:349, "Model validation using
  equivalence tests"): the seminal statement of equivalence testing *as a
  model-validation tool* — require the discrepancy to fall inside a pre-set
  region rather than fail to reject an NHST null. The S6 gate is a transfer of
  this idea to a causal-set count.
- **Siebert–Ellenberger 2020** (Transportation 47:3031, arXiv:1802.03341):
  the nearest applied precedent for equivalence-testing a *count* against a
  reference tolerance; **Stucke–Kieser 2013** (Biom. J. 55:203): equivalence /
  noninferiority testing specialized to Poisson counts (two-sample; the S6
  stage is the one-sample analogue against a certified expected value).

The framework is therefore a domain transfer, not a new statistical
apparatus. Its novelty in this program is the *preregistered application*: a
frozen equivalence gate (TOST-family, each Garwood tail at level 0.025) that
requires the identified discrepancy of the conservative outer Garwood interval
(lower at `K_certain`, upper at `K_certain + U_amb`) against a directed-rounding
certified enclosure to lie within `±B`, with a tri-state (concordant /
discordant / inconclusive) verdict, executed across a curvature ladder.

## Alfyorov–Shnyukov corpus — version-history due diligence

The oracle-arc document assessed "Weyl curvature from the Hasse diagram"
(Zenodo, INSPIRE 3142005) as orthogonal. The mass generalization prompted a
fuller sweep of the same group, because a curvature-ladder framing is closer to
their subject matter. Findings:

- **No distinct January-2026 version is verifiable.** An internally flagged
  "2026-01 ResearchGate preprint" of the Weyl-curvature work could not be
  confirmed as a separate item; the earliest fetchable predecessor of the
  Zenodo record is a 2026-03 code/data deposit. The catalogued entry cites a
  real, resolvable **version-of-record** DOI (10.5281/zenodo.19364212,
  2026-04-01); Zenodo additionally mints a **concept DOI** that resolves to the
  latest version. Automated reads of the concept-DOI value disagreed, so the
  bibliography retains the verified version-of-record DOI rather than asserting
  an unconfirmed concept-DOI number. This is a citation-hygiene note, not a
  priority finding.
- **A new same-group sibling exists and is also orthogonal.**
  **Alfyorov 2026** (Research Square rs-9721442, doi:10.21203/rs.3.rs-9721442,
  2026-05-18): a causal-set sprinkling of the Schwarzschild *interior* arguing
  singularity resolution via a bounded discrete curvature proxy and a
  de Sitter (Hayward) core. This is interior-structure / curvature-proxy
  estimation — no directed-rounding certified volume, no `μ = 2M/r` exterior
  ladder, no TOST/Garwood equivalence framework. Orthogonal.
- **Corpus verdict: firmly orthogonal.** Across the group's Papers 1–8 plus the
  de Sitter-core preprint, the causal-set items are curvature/Weyl *estimation*
  (Paper 7) and interior singularity resolution (de Sitter core); the remainder
  are spectral-action / noncommutative-geometry field theory. None performs a
  mass/curvature-parameterized sprinkling-count-vs-certified-volume equivalence
  verification. The oracle-arc orthogonality assessment stands and is
  strengthened.

## Landscape note (orthogonal, surfaced for completeness)

- **Eichhorn–Gamito–Stokes 2026** (arXiv:2605.06813): discrete diagnostics for
  black-hole horizons and geodesic focusing in a 1+1D causal-set toy model.
  Orthogonal to count-vs-volume; note that its "ladders" are geodesic tracer
  chains — a false cognate for the S6 mass ladder, worth flagging so the two
  are not conflated.

## Summary

| Cluster | Prior art | S6 novelty scope |
|---------|-----------|------------------|
| Mass ladder of count-vs-certified-volume | Homšak–Veroni (single mass); Roy–Sinha–Surya (analytic); Barton et al. (area) | No direct prior — certified gate across a non-isometric curvature ladder |
| μ = 2M/r ladder axis | Benincasa–Dowker, Dowker–Glaser, Belenchia et al., Barton–Borza–Röhrig (discrete curvature *estimators*) | No overlap — μ is a continuum label, not an estimator |
| N ∼ ρV in curved spacetime | Sorkin 1997; Ahmed et al. 2004 (Everpresent Λ) | Correspondence is prior art; the certified equivalence gate is ours |
| Equivalence-testing framework | Schuirmann, Lakens ×2, Garwood, Clopper–Pearson, DeLong, Robinson–Froese, Siebert–Ellenberger, Stucke–Kieser | Methodology transfer — novelty is preregistered application vs a certified enclosure |
| Curved-spacetime causal-set (Alfyorov group) | Orthogonal (curvature estimation, interior structure) | No conflict; assessment strengthened |

The allowed generality statement is unchanged by this search and is bounded
accordingly: on each preregistered rung `μ ∈ {0.1333, 0.1867, 0.2400}`
independently, an operational count realized the rung's certified volume within
`τ = 2.5%`. No interpolation, extrapolation, or composite cross-rung statistic
is claimed, and nothing here promotes the ladder into Paper A's claim set — the
manuscript continues to state the prediction-anchored count stage as its own
boundary until the program integrates it.
