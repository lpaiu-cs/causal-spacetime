# Cover letter (Classical and Quantum Gravity)

*Final text, rev. 5 — four rounds of review of the letter itself. This
round (colleague review): rewritten in plain language for a desk editor
who is likely not a causal-structure specialist — every term of art is
glossed at first use — and the concept figure added. Submission is
performed by the author.*

Dear Editors of *Classical and Quantum Gravity*,

Please consider the enclosed manuscript, "Spacetime quantities from
causal order: an operational ladder with preregistered validations at
finite density," for publication in *Classical and Quantum Gravity* as a
Paper.

A standard result in Lorentzian geometry says that the causal structure
of a spacetime determines its conformal geometry — the pattern of light
cones — while volume information fixes the remaining scale. The
manuscript asks a more practical version of that statement: given only a
finite set of events and the causal order between them, which spacetime
quantities can actually be reconstructed, and what additional
information is needed for each one?

![The manuscript's central object: each row of the reconstruction
hierarchy adds one declared ingredient and names the spacetime quantity
it unlocks — and what that ingredient set still does not
supply.](figures/png300/fig1_ladder.png){width=5in}

We answer with a step-by-step reconstruction hierarchy, summarized in
the figure above.
Starting from causal order alone, we add one declared ingredient at a
time — an overall event density, an observer's clock, an orientation
reference, overlapping observer charts, a local volume profile — and
measure what becomes recoverable at each step, with its error at finite
sampling density. Just as importantly, we show what cannot be recovered
while an ingredient is missing: causal order alone determines neither
the overall scale nor how the volume element varies from place to place;
a single observer cannot tell left from right (a reflection ambiguity);
and a finite signal speed by itself does not produce Lorentzian
geometry.

We then test the framework in curved spacetime, in numerical experiments
that were preregistered: every decision rule was fixed and archived
before the data were generated. The first test uses a vacuum plane wave
built so that the volume element and the random placement of events (a
Poisson "sprinkling") are exactly the same as in flat spacetime, and
only the light cones change. Even at finite sampling density, a
statistic computed from the causal order alone separates the flat and
curved ensembles — and since nothing about the volume differs, the
signal cannot be attributed to volume effects. Matching tests on a
Schwarzschild exterior region confirm the same light-cone effect in a
paired design and show that a single sampled set already discriminates
the two geometries above chance. A final four-mass series then turns to
the volume half of the standard result: event counts from sprinkling
agree, within a preregistered 2.5% tolerance, with continuum
four-volumes certified in advance by interval arithmetic (computer
arithmetic that carries guaranteed error bounds) — so the comparison
runs against a certified prediction rather than a tuned margin.

We believe the manuscript is well suited to *Classical and Quantum
Gravity* because it directly addresses the relation between causal
structure and spacetime geometry, and the reconstruction of continuum
geometry from discrete data. Its contribution is to make that
reconstruction operational: what can be recovered from finite causal
data, what extra information each step requires, and what ambiguities
remain when that information is absent — with the curved-spacetime tests
showing how these statements behave at finite density rather than only
in an ideal continuum limit. The computational claims are fully
auditable: the decision rules were archived before execution, the
numerical evidence is version-controlled in a public repository, and
automated checks rebuild the reported figures from the archived data.

The manuscript is original, has not been published previously, and is
not under consideration elsewhere. There are no conflicts of interest.
Supporting data are cited in the Data Availability Statement. The use of
AI assistance in this work is disclosed in the Acknowledgements, per IOP
policy.

Thank you for your consideration.

Sincerely,

Juneyoung Kim\
Independent researcher\
lpaiu.cs@gmail.com

---

## Submission-form extras (not part of the letter)

**Article type:** Paper.

**Related manuscripts:** none under consideration anywhere. Form-ready
sentence, if asked: "A companion study building a discriminator on this
foundation is in preparation and unsubmitted; the present manuscript is
self-contained and does not depend on it." (The manuscript's Discussion
names it once.)

**Suggested reviewers.** Before entering them on the form: verify each
person's current affiliation and e-mail from their institutional page,
and check for recent co-authorship, shared institutions, or other
conflicts. Final selection is editorial.

1. Sumati Surya — author of the *Living Reviews in Relativity* survey of
   causal set theory and of the small-diamond discrete-geometry results
   the count stage is positioned against.
2. Fay Dowker — co-author of the Lorentz-invariance/discreteness argument
   and of the Benincasa-Dowker curvature estimator the manuscript names
   as the complementary channel.
3. David Rideout — co-author of the Schwarzschild causal-relation
   algorithm (He-Rideout) the manuscript's exact predicate is checked
   against, and of large-scale causal set numerics.
4. Lisa Glaser — causal set numerics and discrete d'Alembertians;
   well placed to judge the estimator-side claims.
5. (optional, if the form allows a fifth) Renate Loll — causal dynamical
   triangulations; a neighbouring discrete Lorentzian programme, well
   placed to judge the manuscript's broader spacetime-reconstruction
   positioning rather than its causal-set internals.

**Opposed reviewers:** none.
