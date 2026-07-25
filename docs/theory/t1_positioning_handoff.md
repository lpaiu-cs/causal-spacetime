# Handoff: T1 literature positioning — verification and document repair

Status: **HANDOFF v1.0, issued 2026-07-25.** For a worker session. The
issuing session was a review/audit session and deliberately performed no
edits to `t1_parallax_identifiability.md`, `t1_g4c_proof.md`, or
`docs/paper/paper_b/manuscript.md`.

**Source of this work order:** `docs/theory/t1_literature_positioning.md`
(positioning audit v1.0). Read that document first — it carries the
evidence, the per-claim verdicts, and the reasoning. This file carries only
the instructions.

**Prerequisites for the worker:** read, in order,
`docs/theory/t1_literature_positioning.md`,
`docs/theory/citations/t1_citation_verification_report.md`, then the two
target documents. Do not start editing before finishing item A.

---

## 0. Ground rules

These are hard constraints. A change that violates one of them is a defect
even if it looks like an improvement.

1. **No proof changes.** Every argument in `t1_parallax_identifiability.md`
   and `t1_g4c_proof.md` stands as written. This work order changes what is
   *claimed about priority*, not what is proved.
2. **No proof-status tag may move.** No `[PROVED]` becomes `[CONJECTURED]`
   or vice versa. The one exception is item B3, which *adds* a new
   proposition and may tag it — see that item.
3. **No measured number may change.** Exponents, standard errors, confidence
   intervals, thresholds, seed counts: all frozen. This includes the numbers
   in `docs/theory/*.json`.
4. **No frozen artifact is touched.** Nothing under `docs/prereg/frozen/`,
   no gate, no threshold, no PC-V1 scene definition.
5. **No experiment is re-run and no new measurement is made.** If an item
   seems to need one, stop and report rather than improvising.
6. **Follow the repository's citation discipline.** Verify every entry
   against an authoritative source; never fabricate a field; omit a field
   you cannot confirm and record the omission in the verification report.
   The existing reports are the format to match.
7. **Revision-note convention.** `t1_parallax_identifiability.md` records
   every revision as a numbered note in its "Revision notes" section
   (currently notes 1-15, ending at v1.0). This work is **note 16**, and the
   status header at the top of the file gets a corresponding version bump.
   Do not skip this — the document's audit trail is load-bearing.

---

## A. Blocking research task: the TDOA self-calibration corpus

**Why this is first.** The audit rates the residual risk at **~40%** that a
general-`d` identifiability count matching T1's `(R = d+2, n = d^2+3d+1)`
already exists in the acoustic/RF sensor self-calibration literature. That
literature solves exactly the "unknown receivers *and* unknown sources"
problem. It publishes minimal solvers for small fixed `(d, R, n)` rather
than general-`d` theorems, is partly conference-only, is patent-adjacent,
and is poorly indexed — which is why the audit's search did not settle it.
**No priority claim on ledger items 9-12 may be written until this is
closed.**

### Scope

Search and read, at minimum:

- The Lund group line: Kuang, Åström, Oskarsson, Burgess, Larsson. Start
  from `kuang2013` (ICASSP 2013, microphone position self-calibration) and
  follow citations forward and backward. Include the EUSIPCO papers
  (Kuang & Åström, stratified TDOA self-calibration; Ask, Kuang & Åström,
  minimal problems in collinear and planar TDOA self-calibration) — these
  were named in the audit but their minimal-configuration tables were not
  extracted.
- `ferranti2021` (ICASSP 2021, arXiv:2005.10298) — **the audit could not
  extract its tables; PDF text extraction failed. Retrieve them from the
  published version.**
- `elbadawy2023` (IEEE TSP 71, arXiv:2102.03565).
- `crocco2012`, `thrun2005`, `pollefeysnister2008`.
- arXiv:2410.19772 (Cao et al.) — reports prior minimal TOA joint
  sensor-source configurations in 3D as 4 sensors/6 sources and 6
  sensors/4 sources. Verify and follow its reference list.
- Any survey of "anchor-free localization", "sensor network
  self-calibration", "blind calibration", "unlabeled TDOA".

### The question to answer

For each relevant paper, record:

1. Are receiver positions unknown? Are source positions unknown?
2. Are measurements labeled or unlabeled?
3. Is the result an *identifiability theorem* (for which `(d, R, n)` is the
   configuration determined, and up to what group?) or a *minimal solver*
   (an algorithm for one fixed small case)?
4. Does the paper state, in any form, a general-`d` requirement of
   `d+2` receivers, or a target/source count of `d^2+3d+1`, or
   `(R, n) = (4, 11)` at `d = 2`, `(5, 19)` at `d = 3`, `(6, 29)` at
   `d = 4`, `(7, 41)` at `d = 5`?

### Also settle this

The audit performed a degrees-of-freedom count that it explicitly marks as
an **unverified derivation, not found in the literature**: for TDOA
self-calibration with `m` receivers and `n` transmitters, measurements
`(m-1)n` against unknowns `d(m+n) - d(d+1)/2` make `m = d+1` infeasible for
every `n`, while `m = d+2` requires `n >= d(d+3)/2` — giving `(4, 5)` in 2D
and `(5, 9)` in 3D. **Those target counts differ from T1's
`d^2 + 3d + 1` (11 and 19).** Determine why. Candidate explanations to test:
the naive DOF count is not tight (T1's Theorem 2a bound is a genuine
strengthening); or T1's unlabeled quotient adds constraints the count
misses; or the two are counting different things. This discrepancy is
currently unexplained and a referee will find it.

### Deliverable

A new section appended to `docs/theory/t1_literature_positioning.md`
(§9, "Self-calibration corpus: resolution of the residual gap") stating:

- the verdict: prior art found / not found, with the specific papers checked;
- if found, exactly what it says and how T1's statement differs;
- the resolution of the DOF discrepancy above;
- an updated residual-risk figure replacing the ~40%.

Plus new verified entries in `docs/theory/citations/t1_references.bib` and
their rows in the verification report.

**If prior art IS found:** stop, do not proceed to section B items C1/D1,
and report. The novelty ledger in the audit needs rewriting first, and the
journal decision changes.

---

## B. Document repair

Two target files. Line numbers are as of commit `06cb15f`; quoted anchor
text is authoritative if lines have drifted.

- **T1** = `docs/theory/t1_parallax_identifiability.md` (1697 lines)
- **PB** = `docs/paper/paper_b/manuscript.md` (1709 lines), Section 8
  (lines 837-1350)

Note on PB: its Section 8 currently carries only two pandoc citation keys
across ~500 lines (`[@boucheron2013]`, `[@robinson1951]`). Every key added
to PB must also exist in `docs/paper/paper_b/citations/references.bib` with
a row in `docs/paper/paper_b/citations/citation_verification_report.md`.
The keys in `docs/theory/citations/t1_references.bib` were chosen to be
collision-free with Paper B's, so entries can be copied across verbatim.

### Group B — corrections (highest priority; these are defects)

#### B1. Retract the G2 fluctuation-class novelty framing

**Problem.** The 1+1D Minkowski-sprinkling longest chain *is* Poissonian
last-passage percolation, not an analogue of it: in null coordinates
`u = t - x`, `v = t + x` the causal order becomes coordinatewise order, the
linear map has constant Jacobian so the Poisson intensity is preserved, and
an Alexandrov interval becomes an axis-aligned rectangle. Consequently both
measured fluctuation classes are **published theorems**:

- order-only chain, `theta = 1/3` and transverse `xi = 2/3` — `bdj1999`,
  `johansson2000`;
- tube-confined chain, Gaussian `theta = 1/2` — `deyjosephpeled2024`
  (increasing paths in `[0,n]^2` confined to a strip of width `n^gamma`,
  `gamma < 2/3`, have Gaussian fluctuations with variance
  `n^{1-gamma/2+o(1)}`; at `gamma = 0` that is `sd ~ L^{1/2}`),
  **including the confinement mechanism** the documents derive
  independently.

**Locations.**
- T1 §5, subsection `### The count-fluctuation class (v1.0): ...`
  (line 720 onward), including the mechanism paragraph and the
  tube-width table around lines 742-765.
- T1 §6, gap G2's closure note (lines 936-980), specifically the
  "settled in v1.0, and the answer is that there is no single class"
  passage at lines 965-978.
- PB §8.4, lines 1044-1060 ("fluctuates like a longest chain should
  (Tracy-Widom, exponent 1/3 ...)" through "... swamped by the wandering
  it deliberately allows").

**Required change.** Reframe from discovery to **validation against exact
theory**. Concretely:

1. Add the null-coordinate identification explicitly, as its own short
   paragraph, citing `bollobasbrightwell1991` (states it as motivation) and
   `bachmat2007` (explicit LIS ↔ 1+1D Lorentzian dictionary).
2. State each measured exponent as confirming a named theorem, with the
   citation adjacent to the number.
3. Attribute the confinement mechanism to `deyjosephpeled2024`.
4. Add the non-conflation warning: `bodineaumartin2005` and
   `baiksuidan2005` show LPP to `(n, n^a)` in a thin *rectangle* still
   converges to Tracy-Widom. Thinness alone does not destroy KPZ;
   constraining the *path* while endpoints stay on the diagonal does. T1's
   tube must be described in the Dey-Joseph-Peled sense explicitly.
5. Preserve what survives as contribution, and say it plainly: (a) the
   naming gap — `brightwellluczak2015`, a dedicated survey, cites `bdj1999`
   yet contains zero occurrences of "Tracy-Widom", "KPZ", "Kardar",
   "last-passage" or "Hammersley", and asks openly for tight variance
   bounds in general `d`; (b) `d >= 3` is genuinely open, and that is where
   a scientific claim can live.

**Do not** delete the measurements, weaken the standard-error statements,
or remove the `count_class_status` provenance note.

#### B2. State the transverse-exponent normalization

**Problem.** T1 line 636 reports "transverse RMS exponent `-0.168` (95% CI
`[-0.192, -0.144]`), within `0.002` of the KPZ wandering value `-1/6`".
Johansson's theorem gives `xi = 2/3`. The two agree only under a
normalization that is currently implicit (transverse `~ N^{1/3}` against box
side `~ N^{1/2}`, giving `rho^{-1/6}` in density units). As written a
probabilist referee reads `-1/6` as contradicting `2/3`. **This is a
presentation defect that looks like an arithmetic error.**

**Locations.** T1 lines 633-647 and 659-661; PB lines 1020-1027 and the
figure legend at line 1314.

**Required change.** Add one explicit sentence giving the normalization and
the conversion, adjacent to the first occurrence in each document, citing
`johansson2000`. Every later occurrence of `-1/6` should be readable
against it.

#### B3. Upgrade the hull-vs-TDOA-fold argument to a proposition

**Problem.** The centered observable `w = M Phi` is a TDOA (range-difference)
measurement with an unknown common offset. In that setting, with **known,
labeled** receivers, `R = d+1` gives a generic two-fold ambiguity and
`R >= d+2` gives uniqueness — `calhoun2021` states the general-`d` form,
`compagnoni2014` is the rigorous `d = 2` treatment with the bifurcation
locus mapped. T1's `d+2` is **a different fact**, arising from unlabeled
rigidity, and T1 says so — but only in one bullet, and only `[MEASURED]`:

> T1 lines 869-874: "`Phi~(x') = Phi~(x)` is the TDOA condition ... and
> three-receiver TDOA is famously two-valued — so one expects centering to
> cost an observer. It does not, for targets in the hull: a dense scan finds
> no second zero ..."

A dense numerical scan standing against a published genericity theorem will
not survive review, and the coincidence that both thresholds are `d+2`
makes the misreading ("they rediscovered a 1987 fact") nearly automatic.

**Locations.** T1 lines 863-876 (the G4b bullet list in §5b) and the
corresponding revision note at lines 1720-1727. PB §8.5, around line 1160.

**Required change.**

1. Promote the hull claim to a numbered proposition with the hull
   hypothesis doing explicit work. **If it can be proved, tag `[PROVED]`
   and write the argument. If it cannot, tag it `[MEASURED]` and say
   explicitly that it is a numerical finding standing against a genericity
   theorem, naming the tension.** Do not leave it as an unlabeled bullet
   either way. *This is the one item permitted to add a tag; see ground
   rule 2.*
2. Cite `compagnoni2014` and `calhoun2021` at the point where TDOA is
   invoked.
3. Add one explicit sentence: T1's `R >= d+2` and TDOA's `R >= d+2` are
   different facts with different mechanisms. Say which is which.

### Group C — citations (positioning; none requires new mathematics)

#### C1. Gortler-Theran-Thurston reconciliation and vocabulary discipline

`gtt2019` (*Forum of Mathematics, Sigma* **7** (2019) e21) proves that for
`d >= 2` and `n >= d+2` generic points, **labels have no effect on generic
uniqueness**. Per the audit, not citing it is the single most likely
objection from a math-aware referee.

**Required change**, in T1 §5b and §6 (G4b/G4c) and PB §8.5-8.6:

1. Add a reconciliation paragraph. The distinction to make: in `gtt2019`
   and its successors (`garamvolgyi2021`, `garamvolgyi2022`,
   `connelly2024`, `gkioulekas2024`) the measurements are Euclidean lengths
   *among the unknown points themselves*, and "unlabeled" is a
   combinatorial quotient. T1's observable is the distance structure of the
   **images** `w_j = w(x_j)` under a map whose parameters — the observer
   positions — are themselves unknown. T1's permutation-invariance is
   automatic (`D` depends on `w` only through Euclidean distances), so
   GTT's hard content is not what T1 proves.
2. **Enforce infinitesimal-vs-global vocabulary throughout.** T1
   establishes *infinitesimal* rigidity in the exact model; GTT proves
   *global* generic uniqueness. **The phrase "global rigidity" must not
   appear describing T1's result** — in the documents, and especially in
   any abstract or title. T1 §5b already scopes this correctly; extend the
   same discipline to §6 and to PB.
3. Position within the right framework: identifiability under algebraic
   constraints (`cruickshank2023`) and the frameworks-on-surfaces program,
   which always assumes the surface is *known and fixed*. T1's surface is
   parametrized by hidden observers. That is the gap being filled.
4. Note the arithmetic finding: `d^2+3d+1 = (d+1)(d+2)-1` returned zero
   hits across rigidity, distance geometry and localization; `>= d+2` is
   ubiquitous but always bounds the number of *configuration points*, never
   *anchors*.

#### C2. Cite and distinguish the nearest prior identifiability theorem

`humeyermeyer2025` (*J. Math. Phys.* **66** (2025) 122502) shows distances
between causally related pairs determine distances between spacelike pairs
in any dimension, removing the conformal ambiguity. **This is the closest
prior identifiability result in the causal-set literature and is currently
uncited.** Add it to T1 §1 or §8 and to PB §8.2, with the distinction: it
assumes exact continuum separations; T1 works from finite order data
through a measured instrument. Also add `madsen2026` and `muller2025` where
embedding uniqueness / the Hauptvermutung is discussed.

#### C3. Engage the current competing spatial-distance proposal

`bogunakrioukov2024` (*Phys. Rev. D* **110** (2024) 024008) measures spatial
distances via causal overlaps, is explicitly intrinsic to the causal-set
graph, and **contains an explicit critique of `rideout2009`**. It postdates
this program's implicit baseline. Add a short related-work paragraph to T1
§1 and PB §8, placing T1's trade against the four existing families:

| proposal | supplied ingredient |
|---|---|
| `rideout2009` | dimension `d`; equidistant-set existence only conjectured |
| `eichhorn2019` | mesoscale cutoff |
| `bogunakrioukov2024` | dimension-dependent constant |
| `major2006`, `major2007` | choice of antichain |

T1 replaces those with "supply an observer congruence" — a different trade,
not a strictly better one, since the chain is itself supplied protocol
structure (axiom A4). State it that way.

#### C4. Define "order-intrinsic" explicitly

`surya2019` §4 defines an **order invariant** as a function independent of
the labelling of the causal set. T1 uses "order-intrinsic" in a stronger
sense — no embedded coordinates — which matters for the G2 order-only
harvest and its coordinate-tube contrast. **Define the stronger notion
explicitly in T1 §2, anchor it to `surya2019` §4, and note that the
stronger sense is used informally in the literature rather than
axiomatized.** Supporting phrasings to cite: `rideout2009` ("defined
intrinsically to the causal set"), `bogunakrioukov2024` ("without any
reference to an embedding continuous manifold").

#### C5. Reattribute Lemma 4c's decoding step

T1 Lemma 4 (line 286) and its part (c) (lines 305-312) derive order
recovery from the strict Robinson structure. That decoding is textbook
seriation. **Keep Lemma 4c as T1's own lemma** — it is correct and short,
and the audit could not locate the strict linear case as a numbered theorem
anywhere (`armstrong2021` defines "strictly linear Robinson" but proves
only the circular analogue). But attribute the *framework*:

- `atkins1998` — closest explicit statement (no repeated Fiedler values ⇒
  the two monotone permutations are the only compatible orders);
- `preafortin2014` — set of compatible orders is a PQ-tree; strictness
  collapses it to one Q-node;
- `carmona2025` — "flat" Robinson space is the modern name for
  "unique up to reversal";
- `chepoifichet1997`, `laurentseminaroti2017` for the recognition line;
- `brainerd1951` alongside the existing `robinson1951`.

The contribution to claim is that the parallax observable **has** that
structure, not the seriation itself.

### Group D — framing

#### D1. Promote the observer-chain construction

Per the audit this is the **clearest novelty in the program** — no prior
work builds observer worldlines of ticks on a causal set and defines
spatial structure from bracket widths — and it is currently buried under
the theorems it supports. Nearest neighbours, none of which occupy the
space: `sivaramakrishnan2025` (a "lightbulb clock", but continuum
perturbative quantum gravity, not causal sets); `bevilacqua2026` (the
closest observer-centric causal-set paper, different construction).

**Required change.** Surface it in T1 §1 (Goal) and in PB §8's opening, as
a stated construction-level contribution with the negative-search result
recorded. Keep the honesty discipline: the chain is supplied protocol
structure, so this is a new *route*, not a derivation from order alone.

**Blocked on item A** — do not write a priority claim here until the
self-calibration corpus is settled.

---

## C. Bibliography wiring

1. Copy the entries needed by PB from `docs/theory/citations/t1_references.bib`
   into `docs/paper/paper_b/citations/references.bib`. Keys are already
   collision-free; keep them identical across both files.
2. Add the corresponding rows to
   `docs/paper/paper_b/citations/citation_verification_report.md`, matching
   its existing table format and stating the confirming source.
3. Two known bibliography defects to fix while there, both recorded in
   `docs/theory/citations/t1_citation_verification_report.md` §5:
   - `armstrong2021`'s `year` field is the arXiv preprint year, not the
     journal year. Correct it.
   - `aghili2018`'s title string was confirmed only at DOI-record level.
     Re-read it off the article page and correct if needed.
4. Do not add an entry you have not verified yourself. Do not add
   Bombelli-Noldus-Tafoya (arXiv:1212.0601) — withdrawn by its authors.

---

## D. Acceptance criteria

- [ ] Item A closed with a verdict and an updated residual-risk figure.
- [ ] All B/C/D items applied to **both** T1 and PB where the item names both.
- [ ] `t1_parallax_identifiability.md` status header version-bumped and
      **revision note 16** written, in the style of notes 1-15.
- [ ] Every pandoc key used in PB resolves in
      `docs/paper/paper_b/citations/references.bib`.
- [ ] Every new bib entry has a row in the relevant verification report,
      with its confirming source and any omitted field recorded.
- [ ] `pytest` passes.
- [ ] `ruff` clean (no code should change, but the repo gates on it).
- [ ] Diff review: no measured number changed, no proof-status tag moved
      except B3's new proposition, no frozen artifact touched.
- [ ] The phrase "global rigidity" does not describe T1's result anywhere.

## E. Out of scope for the worker session

- Any new experiment, measurement, or re-run.
- Any change to `docs/theory/*.json` results tables.
- Any change to gates, thresholds, or `docs/prereg/frozen/`.
- The journal/venue decision, and the Paper B ↔ T1 split decision. Both are
  open questions for the issuing session; the worker should not
  pre-empt them in prose. In particular, **do not restructure Section 8 of
  Paper B on the assumption that T1 will be split out.**
- Writing a T1 standalone manuscript.

## F. Provenance

- Audit and this handoff issued from a review session on 2026-07-25, at
  repository commit `06cb15f`.
- The audit's literature searches were performed by four parallel research
  passes covering: unlabeled rigidity/uDGP; TDOA and localization
  identifiability; causal-set spatial reconstruction; KPZ/LPP and seriation.
  Their residual gaps are recorded in
  `docs/theory/citations/t1_citation_verification_report.md` §4 and
  `docs/theory/t1_literature_positioning.md` §8.
- No file under `docs/theory/` other than the three new ones
  (`t1_literature_positioning.md`, `citations/t1_references.bib`,
  `citations/t1_citation_verification_report.md`) and this handoff was
  modified by the issuing session.
