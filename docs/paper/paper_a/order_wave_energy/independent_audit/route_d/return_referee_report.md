# Referee report on the far-start draft

Saved by the orchestrator on 2026-11-25 from the text the referee returned.
The referee is an AI agent in a fresh context. It reviewed
`return_design_draft.md` at commit `1fa4fa0`, was told not to start a chain
from a configuration pushed along a non-family field in 3+1, and wrote its
scripts outside the repository; they are not archived here. Everything
below the line is the referee's text.

---

**Outcome: do not freeze as drafted.**

As written, the Section 3 path never gets far, so the design can only
return "Not far enough". No chain was started from a pushed 3+1 start. The
repository is untouched (`git status --porcelain` empty).

## Verdicts

| # | Item | Verdict |
|---|---|---|
| 1 | Section 1 logic | GAP |
| 2 | Push fields, `P`, `M` | CORRECT for fields, fit and noise; GAP: the built start is not along F |
| 3 | Path | ERROR: 0 of 36 paths reach 2 plateaus |
| 4 | Controls | ERROR (a) cannot be built at N = 128; GAP (b) underpowered; neither exercises `M` |
| 5 | Descriptions | ERROR as written; GAP after repair (power, mislabelled cases) |
| 6 | Cost | CORRECT for the path; GAP for stage 0 |
| 7 | Other | GAP |

## Evidence

**1. Logic.**

- (a) and (c) are kept out of "forgotten" by the distance requirement, but
  they land in "neither" although the chain did not return.
- (b) The one exactly free direction diffuses 0.20 rad, about 0.5 plateau,
  by `N^2` sweeps. Against a 3-plateau push that reads as remembered;
  against the spec path's 0.3 to 0.6 plateau along F it would be noise.
- (c) For the shear, a co-moving push loses projection under soft family
  drift: `M` = 0.96, 0.83, 0.37 for a tilt of 0.3, 0.6, 1.2 rad, with
  distance unchanged. Reference tilts at lag `N^2` average 0.33; one chain
  reached 1.5. The quadratic field is immune (1.00).
- (d) Absorption by the fit is negligible: 0.98 to 1.00 for a dilation or
  boost of 0.08.
- In 1+1, where the law is count-scale broad, the pipeline still returns
  "forgotten". So "forgotten" means narrow relative to the push, not "cell
  scale" as Section 1's first bullet says.
- Per order, `M` tends to minus x0's own offset over the push, not 0.
  "Neighbourhood of x0" should read "the set the reference chain samples",
  and for central events only.

**2. Fields and statistic.**

- Mean square surviving orthogonalisation: shear 0.89 (min 0.68) at N = 128
  and 0.92 at 256; quadratic 0.95 and 0.98.
- For an ideal pinned push, `P(x*)` equals the push exactly and fitted
  parameters are below 1e-9.
- Reference noise of `P` is 0.13 (shear) and 0.12 (quadratic) plateau, flat
  from `N^2/4` to `N^2` (7 orders, N = 128). A clean 3-plateau push gives a
  12-order half-width of about 0.03 on `M`, which is ample.

**3. Path, exactly as specified** (24 paths at N = 128, 12 at N = 256).

| | N = 128 | N = 256 |
|---|---|---|
| distance, plateaus (shear / quadratic) | 0.91 / 0.69 | 0.77 / 0.71 |
| range | 0.28 to 1.59 | 0.36 to 1.24 |
| along F, plateaus | 0.60 / 0.33 | 0.50 / 0.34 |
| seconds per path | about 1 | 2 to 4 |

- The start is closer to x0 than a reference state (1.1 plateau), and 0.28
  to 1.59 plateaus is 1.2 to 6.8 `s1`.
- 35 paths stop on the 1% rule and one on "no theta", after about 32
  programmes.
- The theta rule binds: 47% of programmes have no acceptable theta, and
  half the accepted steps use 1/8.
- Changing the floor (0.1, 0.6), the box (2/N, 16/N), the stall rule or the
  solver gives means of 0.60 to 0.86; none reaches 2.
- Cause: link margins are concave, box-vertex steps move linked events
  against each other, and a floor relative to the current margin lets
  active rows shrink until they sit below solver tolerance (minimum margin
  about 1e-7).
- It is not a corner: 0.1 to 0.2% of rows fall below 1% of their starting
  margin, so the Section 8 statement is wrong. Outer events move 6 to 7/N
  against 8 to 11/N for central ones.
- Repair tested, reach only: floor at 10% of the starting margin, plus a
  per-link bound on the relative spatial step. Theta = 1 in 968 of 972
  rounds. All 24 paths pass 2 plateaus and 14 reach 3 at N = 128; 4 of 6
  reach 3 at 256; 5 to 37 s.
- Even then the F component is only 2.0 (shear) and 1.5 (quadratic)
  plateaus and saturates after 120 rounds. Single outer events are thrown
  74 to 156/N.

**4. Controls.**

- (a) Rotation is an exact symmetry of order, density and proposal, so the
  control chain is a rotated reference chain. It can fail only through the
  fit.
- (a) `30 s1` is 0.68 chart units at N = 128, but a rotation by pi moves
  central events only 0.46. It was constructible in 0 of 19 orders. At 256
  it needs 81 to 117 degrees.
- (a) Substituting 3 plateaus (1.2 to 1.5 rad): memory 0.89 [0.53, 1.26]
  over 5 direct chains, 0.98 [0.85, 1.11] over 7 emulated.
- (a) The fitted "rotation parameter" mixes with tilts. One order read 0.82
  while the composed element kept the rotation to 0.03 rad.
- (b) The start is valid: order identical in 19 of 19, log density lower by
  27 to 42.

| Stiff memory | `N^2/4` | `N^2/2` | `N^2` |
|---|---|---|---|
| N = 64, 12 orders | 0.10 [-0.08, 0.27] | 0.01 [-0.11, 0.13] | -0.15 [-0.45, 0.15] |
| N = 128, 7 orders | 0.03 [-0.19, 0.25] | -0.05 [-0.25, 0.15] | -0.10 [-0.38, 0.18] |

- (b) It passes, with a relaxation time near `0.1 N^2`. Per-order spread is
  0.26, so the chance of passing when the truth is 0 is 0.95 with 12 orders
  but 0.62 with 6.
- Neither control touches `M`, the path or the distance ratio.
- Same chains scored through the main pipeline (six-parameter fit, F the
  released generator): rotation `M` = 0.93 [0.82, 1.04] with distance ratio
  2.7, remembered. Stiff `M` = -0.10 [-0.37, 0.17] with ratio 1.01
  [0.86, 1.17], "neither" at 7 orders.
- 1+1 end to end (N = 64, 13 orders, exact equilibrium by direct sampling,
  push 2.8 plateaus): `M` = 0.82, 0.36, 0.08, -0.02 [-0.11, 0.07] at `N`,
  `N^2/8`, `N^2/4`, `N^2`; ratio 0.92 [0.79, 1.04]; forgotten, as known.
- A smooth-map start with the same `P` gives 0.87, 0.37, 0.19, 0.04, so the
  corner transient is small there. Its ratio [0.72, 1.11] returns "neither"
  by 0.01.

**5. Descriptions.**

- As written, only "Not far enough" can be returned.
- After repair, the ratio of two equilibrated chains has a per-order spread
  of 0.14 to 0.19. The chance its upper end is at most 1.10 when the truth
  is exactly 1 is about 0.4 to 0.65 with 12 orders, so "neither" is roughly
  as likely as "forgotten" under a perfect return. At spread 0.15 it is
  0.88 with 24 orders and 0.99 with 48.
- The "remembered" consequence overclaims for a slow return.

**6. Cost.** Chains ran at 4.0 to 4.5 ms per sweep alone at N = 128.
Stage 0 at N = 256 is two chains of 65,536 sweeps at 11.9 to 13.2 ms, 26
to 29 minutes, not 20. Stage 1 is about 41 to 50 minutes.

**7. Other.**

- Undefined: "distance" in the stop rule, the base of the 1%, whether
  rejected rounds count, control signs, "wave axis".
- `s1` from a single sweep has a spread of 33% per order.
- There is no power section, which the house rule on interval gates
  requires.
- Moving single events within their own cells would be undone by
  construction and is no substitute for a collective push.

## Changes before freezing

1. **MUST.** Replace the path step rule (start-relative floor and per-link
   trust region, or exact cones), penalise non-central motion, and pilot it
   chain-free. Define "far" by `|P(x*)|` against the reference spread of
   `P`, and report the share not along F.
2. **MUST.** Score controls through the main pipeline: the six-parameter
   versions above and the 1+1 run. Drop `30 s1`; define rotation memory by
   group composition.
3. **MUST.** Add a power section. Raise orders to 24 to 48 or average late
   lags for the distance ratio, and give the stiff control 12 orders at
   N = 256.
4. **MUST.** Add a "not returned" description decided on distance alone.
   Soften the "remembered" consequence, limit "forgotten" to central
   events, and reword Section 1.
5. **SHOULD.** Project on F transported by the fitted element, or use
   time-component fields. Report the spread of `P` by lag and the
   outer-event distance.
6. **SHOULD.** Average `s1` over many sweeps, define the terms in item 7,
   correct Section 8 and the stage 0 time, and add more or random fields.

## Not checked

- Any chain from a non-family pushed start in 3+1, so the 3+1 corner
  transient is unknown.
- Controls at N = 256.
- Whether the repaired path gives sound starts; only its reach was tested.
- More than 7 orders at N = 128. Rotation at N = 64 is degenerate (angle
  pi).
- 1+1 beyond one field at N = 64, with a two-parameter stand-in family.
- Twelve-worker timings, and the chain law itself.
