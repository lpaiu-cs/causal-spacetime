# T1 literature positioning: what is already known, and what is not

Status: **POSITIONING AUDIT v1.0 (2026-07-25).** Companion to
`t1_parallax_identifiability.md` (v1.4) and `t1_g4c_proof.md` (v1.3).
Bibliography: `citations/t1_references.bib`; verification:
`citations/t1_citation_verification_report.md`.

This document exists because the T1 theory documents carry **zero external
citations**. Every claim in them was developed internally, and until now
nothing established which of those claims are new. This audit answers that,
claim by claim.

Verdict tags used below:

- `[NEW]` — no prior statement found; defensible as a contribution.
- `[KNOWN]` — a published theorem already says this. Cite it; do not claim it.
- `[RESCOPE]` — the claim is partly new, but its current framing overstates
  what is new. The document must narrow it.
- `[RETRACT]` — the claim, as an internal discovery, is superseded by prior
  work and must be reframed as a validation rather than a finding.

Nothing here changes a proof. Every `[PROVED]` tag in the T1 documents stands.
What changes is what may be *claimed* about priority.

---

## 0. Executive summary — three findings that change claims

**Finding 1 `[RETRACT]` — the G2 fluctuation-class results are known theorems.**
The 1+1D Minkowski-sprinkling longest chain is not an analogue of Poissonian
last-passage percolation; in null coordinates `u = t - x`, `v = t + x` it *is*
that object. The causal order becomes coordinatewise order, the linear map has
constant Jacobian so the Poisson intensity is preserved, and an Alexandrov
interval becomes an axis-aligned rectangle. Therefore:

- the order-only chain's `theta = 1/3` and transverse `xi = 2/3` are theorems
  (`bdj1999`, `johansson2000`), and T1's measured `0.322` / `-0.168` are
  numerical confirmations of them;
- the tube-confined chain's Gaussian `theta = 1/2` is also a theorem
  (`deyjosephpeled2024`), **including the confinement mechanism** T1
  describes independently.

This is the single largest correction. The measurements are correct and
valuable — as pipeline validation against exact theory — but they are not
discoveries, and presenting them as such is a direct referee target.

**Finding 2 `[RESCOPE]` — `R >= d+2` is a familiar number in localization,
for a different reason, and T1 must say so explicitly.**
The centered observable `w(x) = M Phi(x)` is exactly a TDOA (range-difference)
measurement with an unknown common range offset. In that setting, with
**known, labeled** receivers, `R = d+1` already gives a generic two-fold
ambiguity and `R >= d+2` gives uniqueness — stated in general `d` by
`calhoun2021`, rigorously analysed for `d = 2` by `compagnoni2014`. T1's
threshold is *not* that one: Section 5b explicitly argues the centering fold
does not bite for targets inside the hull, so T1's `d+2` comes from unlabeled
rigidity rather than from the TDOA fold. **That distinction is currently made
only in passing and only `[MEASURED]` (a dense scan). It must become an
explicit, cited, proved proposition**, or every reader from the localization
community will assume T1 has rediscovered a 1987 fact.

**Finding 3 `[NEW]` — the observer-chain / radar-bracket construction appears
unoccupied in causal set theory.** A search of the spatial-reconstruction
literature found no prior work building observer worldlines of ticks on a
causal set and defining spatial structure from bracket widths. The four
existing families (`rideout2009` sphere distance; `eichhorn2019` mesoscale;
`bogunakrioukov2024` causal overlaps; `major2006`/`major2007` thickened
antichains) each require an externally supplied ingredient — the dimension
`d`, a mesoscale, or a choice of antichain. This is the cleanest novelty in
the program and is currently *understated* in the T1 documents.

---

## 1. Novelty ledger

| # | T1 claim | Prior art | Verdict |
|---|---|---|---|
| 1 | Lemma 1 (continuum radar bracket `2\|dx\|`) | Elementary; standard radar coordinates | `[KNOWN]` — not a claim |
| 2 | Lemma 2 (rank-gap identity `W = N + 1`) | Elementary counting; `campbell1909` for the Model P mean | `[NEW]` but elementary. State as a lemma, claim nothing |
| 3 | Lemma 3 (centering removes the shared-scalar gauge only) | Standard gauge algebra | `[KNOWN]` |
| 4 | Lemma 4a-c (strict Robinson ⇒ order up to reversal, `R >= 2`) | `robinson1951`; `atkins1998`; `preafortin2014`; `carmona2025` ("flat" Robinson space = exactly two compatible orders) | `[RESCOPE]` — the *decoding* is textbook seriation. The contribution is showing the parallax observable **is** strictly Robinson |
| 5 | Lemma 4f (`D` determines order, not spacings) | No counterpart found | `[NEW]` |
| 6 | Theorem 1' (multilateration; `d+1` affinely independent observers) | `coope2000`; `larsson2025`; `dokmanic2015` | `[KNOWN]` — textbook; cite, do not claim |
| 7 | Theorem 2 (Bernstein concentration for order recovery) | Standard technique (`boucheron2013`); the shared-region cancellation is a neat but routine step | `[NEW]` as an application; technique `[KNOWN]` |
| 8 | Same-slice pathwise monotonicity, tie probability `exp(-4 lambda g)` | No counterpart found | `[NEW]` |
| 9 | **`R >= d+2` threshold** | `calhoun2021` / `compagnoni2014` give `d+2` for *labeled, known-anchor* TDOA uniqueness — a different mechanism | `[RESCOPE]` — see Finding 2 |
| 10 | **`n = d^2 + 3d + 1` threshold** | Zero occurrences found in rigidity, distance geometry or localization | `[NEW]` |
| 11 | Cone-point argument (finite singular set forbids continuous ambient symmetry) | Folklore Riemannian/orbifold reasoning; the frameworks-on-surfaces program; `cruickshank2023` as the general framework | `[RESCOPE]` — a new *instance* in a known framework, not new machinery |
| 12 | Unlabeled recovery of targets **and** observers up to congruence | `gtt2019` and successors treat unlabeled distances *among the points themselves*; no result found for distances among **images of a map with unknown parameters** | `[NEW]`, with the reconciliation in §2 required |
| 13 | G2: order-only chain is KPZ (`theta = 0.322`) | `bdj1999`, `johansson2000` — theorems about the same object | `[RETRACT]` |
| 14 | G2: tube-confined chain is Poisson (`theta = 0.487`) | `deyjosephpeled2024` — theorem, same regime, same mechanism | `[RETRACT]` |
| 15 | G2: confinement destroys transverse optimization | `deyjosephpeled2024`'s mechanism | `[RETRACT]` |
| 16 | Observer-chain / radar-bracket construction on a causal set | None found | `[NEW]` — the strongest novelty |
| 17 | Chain-fluctuation measurements at `d >= 3` | `brightwellluczak2015` explicitly flags the variance as open for general `d`; only Talagrand-type upper bounds exist | **open** — a real opportunity, not yet done |

---

## 2. Unlabeled rigidity: reconciling with Gortler-Theran-Thurston

`gtt2019` (*Forum of Mathematics, Sigma* **7** (2019) e21) proves, for `d >= 2`
and `n >= d+2` generic points, that **labels have no effect on generic
uniqueness**: a generic configuration is determined by the unlabeled multiset
of point-pair lengths (given `d` and `n`) iff it is determined by the labeled
edge lengths. Follow-ups: `garamvolgyi2021` (when generic global rigidity
fails), `garamvolgyi2022` (globally rigid ⇒ fully reconstructible),
`connelly2024` (`d = 1`, algorithmic), `gkioulekas2024` (unlabeled path/loop
lengths; greedy trilateration seeded on `d+2` points).

**Why T1 is not a corollary of this.** In every one of those papers the
measurements are Euclidean lengths *among the unknown points themselves*, and
"unlabeled" is a combinatorial quotient (which pair produced which number).
T1's observable is the distance structure of the **images** `w_j = w(x_j)`
under a map whose parameters — the observer positions — are themselves
unknown. T1's permutation-invariance is automatic, not hard-won: `D` depends
on `w` only through Euclidean distances, so the `O(R)` quotient is built in.
GTT's hard content is therefore *not* what T1 proves.

**Three things this obliges the paper to do.**

1. **Cite `gtt2019` and state the distinction in the introduction.** A
   math-aware referee will reach for it immediately. Not citing it is,
   per the survey, the single most likely objection.
2. **Never write "global rigidity."** T1 establishes *infinitesimal*
   rigidity (nullity equals the `d(d+1)/2` gauge) in the exact model. GTT
   proves *global* generic uniqueness. T1's statement is strictly weaker and
   the vocabulary must not blur that. `t1_parallax_identifiability.md`
   already scopes this correctly in Section 5b; the paper must keep that
   discipline in the abstract and title too.
3. **Position within the right framework.** The nearest structural home is
   identifiability under algebraic constraints (`cruickshank2023`) and the
   frameworks-on-surfaces program — which, crucially, always assumes the
   surface is *known and fixed*. T1's surface is parametrized by the hidden
   observers. That is the actual gap being filled, and it should be stated
   that way rather than as "a new rigidity theorem."

**Threshold arithmetic.** `d^2 + 3d + 1 = (d+1)(d+2) - 1` returned zero hits
across rigidity, distance geometry and localization searches. `>= d+2` is of
course ubiquitous, but always as a bound on the number of *configuration
points* (`K_{d+2}` is the smallest generically globally rigid graph in `R^d`),
never on the number of *anchors*.

---

## 3. Localization and TDOA: the `d+2` that is not T1's `d+2`

The centered profile `w = M Phi` is TDOA up to an unknown common offset
(`r = w + c 1`). The relevant known facts:

- **`d+1` affinely independent anchors determine a point** — textbook
  multilateration (`coope2000`, `dokmanic2015`); `larsson2025` catalogues the
  degeneracies when the anchors fail to span (collinear senders in 2D give two
  solutions; in 3D, a whole circle). T1's Theorem 1' restates this. Cite it.
- **Range *differences* cost one more.** With known labeled receivers,
  `R = d+1` gives a generic two-fold TDOA ambiguity and `R >= d+2` gives
  uniqueness (`calhoun2021` states the general-`d` form; `compagnoni2014` is
  the rigorous `d = 2` treatment, with the bifurcation locus mapped
  explicitly).
- **Network localizability ⟺ generic global rigidity** (`aspnes2006`,
  `eren2004`) — the standard bridge between the two literatures.
- **Anchor-free self-calibration** — unknown receivers *and* unknown sources —
  is an active applied area (`thrun2005`, `pollefeysnister2008`,
  `crocco2012`, `kuang2013`, `ferranti2021`, `elbadawy2023`). It publishes
  *minimal solvers for specific small `(d, R, n)`*, not general-`d`
  identifiability theorems, and it always keeps labels.

**The required edit.** T1 Section 5b currently disposes of the TDOA fold in
one bullet: the fold "is *not* the obstruction… for targets in the hull: a
dense scan finds no second zero." That is a `[MEASURED]` claim standing
against a published theorem that says two-valuedness is generic. It needs to
become a proposition with the hull hypothesis doing explicit work, cited
against `compagnoni2014`. Until then, the paper reads as if it is unaware of
the TDOA result — and the coincidence that both thresholds are `d+2` makes
that reading almost inevitable.

**The genuine gap.** No source found combines (a) unknown observer positions,
(b) unlabeled measurements, and (c) a congruence-recovery theorem with a sharp
`R` bound in general `d`. That combination is T1's contribution, and it is
best stated in exactly those three clauses.

---

## 4. Causal set spatial reconstruction

State of the art, all order-intrinsic in principle, each with an external
ingredient:

| proposal | mechanism | supplied ingredient |
|---|---|---|
| `rideout2009` sphere distance | `d`-element equidistant sets | dimension `d`; existence of such sets only conjectured |
| `eichhorn2019` | induced geometry on a Cauchy surface | a mesoscale cutoff |
| `bogunakrioukov2024` | volume overlap of light cones | a dimension-dependent constant |
| `major2006`, `major2007` | thickened antichains, homology | choice of antichain |

`bogunakrioukov2024` (PRD **110** (2024) 024008) is the most recent and
contains an explicit critique of `rideout2009` (no null-surface elements;
ambiguous Cauchy-surface choice without an embedding). **It postdates the T1
work's implicit baseline and must be engaged with.**

**Where T1 sits.** The observer chain supplies the surface *and* the scale
internally, which is a different trade rather than a strictly better one — the
chain itself is supplied protocol structure (axiom A4 in this program). The
honest framing: T1 replaces "supply a dimension / mesoscale / antichain" with
"supply an observer congruence," and then proves what that buys. No prior work
makes that trade.

**Prior identifiability theorems to distinguish.** `humeyermeyer2025` (JMP
**66** (2025) 122502) shows distances between causally related pairs determine
distances between spacelike pairs in any dimension — removing the conformal
ambiguity. This is the nearest prior *identifiability* result and **must be
cited and distinguished**: it assumes exact continuum separations, whereas T1
works from finite order data through a measured instrument. `madsen2026` and
`muller2025` bear on embedding uniqueness / the Hauptvermutung.

**"Order-intrinsic" has an established name.** Surya's Living Review
(`surya2019`) §4 defines an **order invariant** as a function independent of
the labelling of the causal set. That is labelling-invariance, which is weaker
than T1's informal "no embedded coordinates" sense — the sense that matters
for the G2 order-only harvest. **The paper must define its own stronger
notion and anchor it to `surya2019` §4**, rather than using "order-intrinsic"
as if it were standard at that strength.

---

## 5. Longest-chain fluctuations: what must be retracted

The identification (§0, Finding 1) is not new either — `bollobasbrightwell1991`
states it as motivation, and `bachmat2007` builds the LIS ↔ 1+1D Lorentzian
dictionary explicitly.

Consequences, in order of severity:

1. **`theta = 1/3` for the order-only chain is `bdj1999`.** Tracy-Widom GUE,
   scale `N^{1/6}`. T1's `0.322` with Poisson excluded at `12.8` se is a
   clean numerical confirmation — report it as such.
2. **`xi = 2/3` transverse is `johansson2000`.** T1's `-0.168 ≈ -1/6` is
   consistent with `xi = 2/3` **only under a stated normalization**
   (transverse `~ N^{1/3}` against box side `~ N^{1/2}`, i.e. `rho^{-1/6}` in
   density units). Currently the normalization is implicit. A probabilist
   referee will read `-1/6` as contradicting `2/3` unless it is spelled out.
   *This is a presentation bug that looks like an error.*
3. **The tube result is `deyjosephpeled2024`.** Increasing paths in
   `[0,n]^2` confined to a strip of width `n^gamma`, `gamma < 2/3`, have
   **Gaussian** fluctuations with variance `n^{1 - gamma/2 + o(1)}`. At
   `gamma = 0` (fixed-width tube) that is variance `~ n`, i.e. `sd ~ L^{1/2}`
   — exactly T1's `theta = 0.487`. Their mechanism is T1's mechanism.
4. **Do not conflate with thin rectangles.** `bodineaumartin2005` and
   `baiksuidan2005` show LPP to `(n, n^a)` in a thin *rectangle* still
   converges to Tracy-Widom. Thinness alone does not destroy KPZ; constraining
   the *path* while the endpoints stay on the diagonal does. T1's tube must be
   described in the Dey-Joseph-Peled sense, explicitly.

**What survives as a contribution.** Two things, and they are real:

- **The naming gap is documented.** `brightwellluczak2015` — a dedicated
  survey, "The mathematics of causal sets" — cites `bdj1999` yet contains zero
  occurrences of "Tracy-Widom", "KPZ", "Kardar", "last-passage" or
  "Hammersley", and asks openly for tight variance bounds on the height in
  general `d`. Importing the KPZ vocabulary into causal set theory, with the
  dictionary stated, is a legitimate service contribution.
- **`d >= 3` is open.** The same survey flags the general-`d` variance as
  unknown, with only Talagrand-type upper bounds. T1's pipeline can measure
  there. **This is where the scientific claim should move.**

---

## 6. Seriation

`robinson1951` (with `brainerd1951`) is the origin; the modern line is
`chepoifichet1997` (first recognition algorithm), `atkins1998` (spectral /
Fiedler-vector seriation), `preafortin2014` (optimal recognition; the set of
compatible orders is a PQ-tree), `laurentseminaroti2017` (similarity-first
search), and `carmona2024siam` / `carmona2025` (modules and PQ-trees;
**"flat" Robinson space** = exactly two compatible orders, an order and its
reverse — the modern name for T1's "up to reversal").

**On Lemma 4c's decoding step.** "A strictly Robinson dissimilarity determines
the underlying linear order up to reversal" is treated as routine in that
literature, but it could **not** be located as a numbered theorem for the
strict linear case. `armstrong2021` defines "strictly linear Robinson" and
proves only the circular analogue. So:

- keep Lemma 4c as T1's own lemma (the max-distance pair must be the
  extremes; sorting the anchor row is then forced) — this is correct and
  short;
- attribute the *framework* to `atkins1998` (closest explicit statement:
  with no repeated Fiedler values the two monotone permutations are the only
  compatible orders), `preafortin2014` (PQ-tree; strictness collapses it to
  one Q-node), and `carmona2025` (the "flat" terminology);
- **do not** present the seriation structure itself as a discovery. The
  contribution is that the parallax observable *has* that structure.

---

## 7. Required edits to the T1 documents, ranked

1. **Retract the G2 novelty framing.** Rewrite Section 5's count-class
   subsection and gap G2's closure note to present the measurements as
   validation against `bdj1999` / `johansson2000` / `deyjosephpeled2024`.
   Add the null-coordinate identification explicitly, citing
   `bollobasbrightwell1991` and `bachmat2007`.
2. **State the transverse-exponent normalization** so `-1/6` visibly agrees
   with `xi = 2/3`.
3. **Upgrade the hull-vs-TDOA-fold argument** from a dense scan to a
   proposition, cited against `compagnoni2014` / `calhoun2021`, and say
   plainly that T1's `d+2` and TDOA's `d+2` are different facts.
4. **Add the `gtt2019` reconciliation paragraph** and enforce
   infinitesimal-vs-global vocabulary discipline throughout.
5. **Cite and distinguish `humeyermeyer2025`**, the nearest prior
   identifiability theorem.
6. **Engage `bogunakrioukov2024`**, the current competing spatial-distance
   proposal.
7. **Define "order-intrinsic" explicitly**, anchored to `surya2019` §4's
   "order invariant" but acknowledged as stronger.
8. **Reattribute Lemma 4c's decoding** to the seriation literature; keep the
   lemma, drop any implicit priority.
9. **Promote the observer-chain construction** in the framing — it is the
   clearest novelty and is currently buried under the theorems it supports.

Items 1-3 are corrections. Items 4-8 are citations. Item 9 is framing.
None requires new mathematics.

---

## 8. Residual verification gaps

Two searches were started and **not exhausted**; both bear on priority.

1. **The TOA/TDOA self-calibration corpus** (Lund group and successors).
   It solves the unknown-receivers-and-unknown-sources problem, but publishes
   minimal solvers rather than identifiability theorems, and it is large,
   partly conference-only, and poorly indexed. Estimated residual risk that a
   general-`d` count matching `(R = d+2, n = d^2+3d+1)` already exists there:
   **~40%.** *Do this hand-check before any priority claim on items 9-12 of
   the ledger.*
2. **Causal-set chain statistics 2015-2026.** `brightwellluczak2015` is good
   evidence for the naming gap but is a 2015 survey. A full sweep of the
   intervening literature was not completed.

A DOF count performed during the audit (**not found in the literature — treat
as unverified derivation**) is suggestive: for TDOA self-calibration with `m`
receivers and `n` transmitters, measurements `(m-1)n` against unknowns
`d(m+n) - d(d+1)/2` make `m = d+1` infeasible for every `n`, while `m = d+2`
requires `n >= d(d+3)/2` — giving `(4, 5)` in 2D and `(5, 9)` in 3D. The
corresponding TOA count permits `m = d+1`. So the `d+1 → d+2` step on passing
from ranges to range-differences reappears in the self-calibration setting.
This is *consistent with* T1's threshold but arrives at different target
counts than `d^2 + 3d + 1`, which is worth understanding before submission.

---

## 9. Self-calibration corpus: resolution of the residual gap

Status: **CLOSED for the priority question, PARTIALLY CLOSED for the corpus
sweep** (2026-07-25). Residual risk revised from **~40% to ~10%**.

### 9.1 The verdict

**No prior art found for T1's count.** The self-calibration corpus solves a
*different, strictly easier* problem, and the discrepancy §8 flagged as
"worth understanding before submission" is the fingerprint of exactly that
difference. It is now understood, and it is not a defect in either count.

### 9.2 The DOF discrepancy, resolved

§8's derivation is **correct**, and it is the count for the **labeled**
problem. T1's Theorem 2a is the count for the **unlabeled** one. They differ
by a term that can be named exactly.

Recall the T1 decomposition (`t1_g4c_proof.md`, Lemma C):

    nullity = dim ker L + dim(Im L  cap  F)

where `L` is the derivative of scene → centered-profile cloud, and `F` is the
flex space of that cloud. The two halves are two different observables:

- **Labeled** (self-calibration): the receiver index is known, so the
  centered profile *vector* `w_j` is observed. The flex space is `ker L`
  alone.
- **Unlabeled** (T1): only the pairwise distances `D(j,k) = ||w_j - w_k||`
  are observed, so the cloud is known only up to an isometry of the profile
  space `V`. The flex space picks up `Im L cap F` as well.

Writing the two counts out with `m = R - 1`:

    labeled     n(m - d)  >=  dR - d(d+1)/2
    unlabeled   n(m - d)  >=  dR + m(m+1)/2 - d(d+1)/2

**The numerators differ by exactly `dim F = m(m+1)/2 = R(R-1)/2`** — the
dimension of the isometry group of the profile space, which is precisely
what the unlabeled observable cannot see. Nothing else differs.

`(4, 5)` in 2D and `(5, 9)` in 3D are therefore the *labeled* thresholds, and
`(4, 11)` / `(5, 19)` the *unlabeled* ones. The gaps are `6 = 3·4/2` and
`10 = 4·5/2`, as the formula requires.

**Verified numerically, not just algebraically.** Measuring the smallest `n`
at which `dim ker L` falls to the rigid-motion gauge:

| `d` | `R` | labeled count | measured | unlabeled (T1 2a) |
|---|---|---|---|---|
| 2 | 4 | 5 | **5** | 11 |
| 2 | 5 | 4 | **4** | 9 |
| 3 | 5 | 9 | **9** | 19 |
| 3 | 6 | 6 | **6** | 14 |
| 4 | 6 | 14 | **14** | 29 |

Five cells, exact agreement. §8's count is not merely plausible — it is the
labeled threshold, attained.

So the third of §8's candidate explanations is the right one: **the two are
counting different things.** T1's Theorem 2a is not a strengthening of the
naive count and does not contradict it.

### 9.3 What the corpus actually publishes

Searched: the Lund line (`kuang2013`, the EUSIPCO stratified and
minimal-problem papers), `ferranti2021`, `elbadawy2023`, `crocco2012`, the
homotopy-continuation line, and the TOA joint sensor-source localization
line including arXiv:2410.19772.

Three properties hold across all of it, and each is independently fatal to
the prior-art hypothesis:

1. **Everything is labeled.** The measurement is a TDOA or TOA *matrix*
   indexed by `(receiver, transmitter)`. No paper found discards the receiver
   index, which is the whole content of T1's observable.
2. **The published configurations are far smaller, exactly as the labeled
   count predicts.** 2D TDOA minimal problems appear at `(4, 1)`, `(3, 2)`,
   `(5, 2)`; 3D far-field TDOA at four receivers and nine transmitters; 3D
   TOA joint sensor-source at `(4, 6)` and `(6, 4)` in the prior art, reduced
   to `(4, 4)` by arXiv:2410.19772. Nothing in this range resembles `11`,
   `19`, `29`, `41`.
3. **The corpus is solver-oriented, not theorem-oriented.** It enumerates
   minimal problems at fixed small `(d, m, n)` and counts their solutions.
   **No general-`d` requirement of `d+2` receivers, and no transmitter-count
   formula of any kind, was found in it.**

### 9.4 What remains open, stated plainly

- `ferranti2021`'s tables were **not** fully extracted. PDF text extraction
  failed here as it did for the original audit; the 2D entries above came
  through, the 3D ones did not. This sub-item is not closed, and the
  published version should be read by hand before submission.
- The EUSIPCO minimal-configuration tables were located by title but their
  contents were not extracted.
- Patent-adjacent and conference-only material remains unswept by
  construction.

None of these can plausibly overturn 9.2, which is an algebraic identity
rather than a search result: whatever those tables contain, they are counts
for a labeled observable and therefore cannot coincide with T1's.

### 9.5 Consequence for the ledger

Items 9-12 are **unblocked**. Item 10 (`n = d^2+3d+1`) stands as `[NEW]`.
Item 9 (`R >= d+2`) stands as `[RESCOPE]` for the reason Finding 2 already
gives — the coincidence with the TDOA fold threshold is real and must be
disclaimed — and §9.2 now supplies the sharp way to say it: *the labeled
problem needs `d+2` receivers too, and needs far fewer targets; T1's extra
targets are the price of discarding the receiver index.*

D1 (promoting the observer-chain construction) is likewise unblocked.
