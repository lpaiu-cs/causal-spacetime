# Referee report on the joint-slack draft

Saved by the orchestrator on 2026-11-25 from the text the referee returned.
The referee is an AI agent in a fresh context. It reviewed
`slack_design_draft.md` at commit `3e8fc31`, was told not to run the
measurement itself in 3+1, and wrote its scripts outside the repository;
they are not archived here. Everything below the line is the referee's
text.

---

**Outcome: do not freeze as drafted.** The mathematics of Sections 2 to 3
is right, but the decision logic cannot return its positive description,
its positive control cannot pass, and the positive description would not
mean what Section 6 says. The repository is untouched
(`git status --porcelain` is empty).

## Verdicts

| # | Item | Verdict |
|---|---|---|
| 1 | Linear model | CORRECT (signs; "not a bound"). GAP: rigorous bounds exist and are unused |
| 2 | Pin | CORRECT, not exploitable. GAP: nearly blind to a common shift of the central cluster |
| 3 | Positive control | ERROR: cannot pass as written |
| 4 | Certification rounds | ERROR: under-certifies where slack is known to exist |
| 5 | Maximum vs typical, Section 6 | ERROR in the stated consequence; GAP in the fit form and the box |
| 6 | Cost | ERROR: 1 to 3 s at N = 1024 is not realistic |
| 7 | Definitions | GAP: several undefined terms |

## Evidence

**1. Linear model.**

- The margin is concave, so tangent link rows relax and tangent critical
  rows restrict. The draft has this right.
- Upper bound: an LP with link rows only is a global relaxation, so its
  optimum is a rigorous upper bound on the pinned slack inside the box.
  Whether it is informative is not known (links are 1.2 to 4.4 per variable
  from N = 128 to 1024).
- Lower bound: exact link cones plus tangent critical rows is a convex
  region inside the true cell, so every solution is a certificate with no
  theta search.
- In 1+1, where the cone is polyhedral, that inner model certified 0.99 to
  1.00 of the exact slack in one programme, with the order identical in
  every case.
- In 3+1 it needs a second-order cone solver, and none is installed in the
  venv.

**2. Pin.**

- The pin is the exact stationarity condition of `fit_group(x0 -> x)` at
  the identity, for any displacement size. Cumulative use with generators
  at `x0` is therefore sound, provided the fit runs from `x0`.
- For pinned fields up to box scale the fit returned parameters below 2e-8,
  and multi-start found no lower minimum in 96 trials per size (N = 128,
  512). The final fit is a no-op for joint programmes.
- A pure common shift of the central events survives a central-set pin at
  2 to 3% mean amplitude (7 to 11% rms over events). It survives a pin over
  non-extremal events at 57 to 77%.
- Section 7's "the rows ... are what stops such a motion" is therefore
  untested by this programme: the pin forbids the motion before the rows
  can.

**3. Positive control.**

- `g_c . G_a` is not zero (up to 0.74). It equals `alpha_c m_c` with
  `|alpha| <= 1.00` for the dilation and boosts, 0.088 for the tilts and 0
  for the axis rotation, tip rows included. So the model is feasible only
  for `|c| < 1`; at N = 128 that binds together with the box.
- A pure dilation inside the box gives the wave objective a median 0.65 of
  `4/sqrt(N)` (the corner is 1.41) and y1 a median 0.025. "Box reached in
  every released programme" is not guaranteed.
- The seven-parameter fit removes 97.7 to 99.8% of a straight dilation
  step, and an exact dilation to 1e-16. So
  `w_cert(released) / w_cert(joint)` is about 1 by construction.
- In 1+1 that ratio had medians of 0.68 to 1.86, never reliably above 2.

**4. Certification.**

- 3+1 single-event programmes work: theta is typically 1/2 five times,
  reaching 97 to 100% of a long-run reference.
- In 1+1 the exact pinned slack is a linear programme on rank-adjacent null
  differences and is at box scale. Against it, the draft loop does poorly:

| N | certified / exact, wave | certified / exact, z | ended "no theta", wave | ended "no theta", z |
|---|---|---|---|---|
| 64 | 0.35 | 0.34 | 14% | 19% |
| 128 | 0.63 | 0.65 | 17% | 17% |
| 256 | 0.19 | 0.21 | 44% | 58% |
| 512 | 0.40 | 0.36 | 49% | 49% |

- `J_cert < 1`, impossible for the true quantity, occurred for 11 to 33% of
  events.
- For the z direction, B for `J_cert` was 7.0 +/- 9.0 while the exact B was
  18.3 +/- 8.7. The draft would return "model slack not certified" where
  box-scale slack exists by construction.
- In 3+1, a straight step along an exactly free generator keeps the order
  only to an rms central displacement of 3.4 to 5.7/N (dilation), 2.5 to
  4.9/N (boosts) or 1.0 to 1.6/N (rotations), at every N. The box is 45/N
  to 128/N.
- Six straight rounds therefore cannot certify count-scale motion with
  relative transverse displacement. This is an inference; I did not observe
  a 3+1 joint loop.

**5. Maximum vs typical.**

- A ball-like cell in 4N dimensions with the chain's typical 12.4/N has a
  one-sided extent of about 8.8/sqrt(N), twice the box, and
  `J^2` of about `(4N+2)/3`.
- So "count-scale slack certified" is the expected outcome for an isotropic
  cell that is typically rigid. "Its plateau cannot be read as rigidity"
  overclaims, and contradicts the draft's own "or they carry little weight".
- A box-saturated model gives `J_model` proportional to sqrt(N)
  automatically, so B > 0 for `J_model` carries no information once the box
  is reached.
- Widths add linearly, so `J = A + B sqrt(N)` is the natural form.
- At N = 128 the wave box corner lies outside the diamond for 55 to 62% of
  central events (18 to 20% at N = 256).

**6. Cost.**

| N | rows | central events | default HiGHS | interior point |
|---|---|---|---|---|
| 128 | 3.6k | 21 | 0.6 to 0.75 s | not timed |
| 256 | 9.6k | 35 | 4.5 to 13 s | not timed |
| 512 | 28k | 71 | 22 to 184 s | 1.3 to 2.2 s |
| 1024 | 76k | 136 | unfinished at 300 s | 9 to 12 s |

- The single-event row set is fine either way: all relations and the
  event's own links and critical pairs gave the same width within 0.5%.

**7. Definitions.** Undefined or inconsistent:

- "reached the box";
- which 16 events are tested, when one N = 128 sample had exactly 16
  central events;
- the per-order J (mean of ratios and ratio of means differed up to 3 times
  in 1+1);
- the base of the 5% stop;
- the fit direction and the residual sign;
- what happens if the control fails;
- Section 6 says "three sizes" where Section 5 lists four.

## Changes before freezing

1. **MUST.** Replace the certification, or prove it on a control. Options
   are the inner model with exact link cones, or a path-following loop
   whose round budget is set by a control. As drafted, "model slack not
   certified" is the expected outcome whether or not slack exists.
2. **MUST.** Replace the positive control and gate the descriptions on it:
   - (a) 1+1 end-to-end against the exact null-coordinate programme,
     requiring certified/exact above a stated fraction;
   - (b) 3+1 released generator, scored on raw displacement or a
     six-parameter fit, with the objective along the released generator;
   - (c) a planted free block: remove the rows between event i's neighbours
     and the rest, so rigid translation is the exact known answer.
3. **MUST.** Rewrite the "certified" consequence in Section 6 as existence
   of extremal configurations only. Add the cheap discriminator: start the
   existing chain at the certificate and see whether the residual relaxes
   to about 12/N.
4. **MUST.** Fix the solver to interior point and re-estimate. The default
   trips the stop rule at N = 512.
5. **MUST.** Define the terms in item 7. For "reached the box", use "the
   optimum changes when the box is doubled".
6. **SHOULD.** Add the links-only programme per tested event as a rigorous
   upper bound; it would turn a negative into a proof inside the box.
7. **SHOULD.** Add a programme pinned over non-extremal events, and reword
   Section 7.
8. **SHOULD.** Report `w_single` as raw displacement. The fit shrinks it 11
   to 13% at N = 128 but 3% at N = 512, which biases J with N.
9. **SHOULD.** Fit `J = A + B sqrt(N)` with a bootstrap over orders, and
   treat N = 128 box statistics separately.
10. **SHOULD.** Put a clearance on certificate margins, floor `w_cert` at
    the single-event certificate, and add a secondary objective so the
    other events' displacement is well defined.

## Not checked

- No 3+1 joint programme was solved with a displacement objective. Joint
  timings used a random dense objective or one outer event, with only
  status and seconds printed.
- Not run in 3+1: the inner model, the links-only bound, the planted-block
  control.
- Interior-point timing is from 3 to 4 programmes per size in one process,
  not twelve workers. I did not check whether it returns a vertex.
- The 1+1 sample is 6 orders x 6 events per size up to N = 512.
- I did run one thing beyond your literal list: an exact-order check along
  straight pure-generator lines in 3+1, with no optimisation.
- Two runs exceeded the 6-minute guideline (9.6 and 7.2 minutes).
