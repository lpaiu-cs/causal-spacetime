# Nested landmarks: precision gain with an unresolved coverage warning

The [frozen eight-pair pilot](nested_pilot_protocol.md), generated at
**c59cdb3edd359d4c20d2ffddb35d26210e7f810f**, completed both arms on the
same 4096 selection and 262144 measurement events per case. The candidate
retained every baseline landmark and used **122--131** interior landmarks,
versus 64 in the baseline. Both use the identical corrected fifth-degree
model, with ten free Weyl and 24 free gradient components.

Actual full-tensor error decreased in **all eight pairs**, including all six
curved cases. The median paired error ratio is **0.5628**: a 43.7% reduction.
The median largest-variance ratio is **0.1837** and the median magnetic-
variance ratio is **0.2036**. This is a substantial gain in these new data,
unlike the earlier wider-center or curvature-square-only comparisons.
These are descriptive pilot findings, not population precision guarantees.

| A | w | Baseline relative error | Nested relative error | Nested/baseline largest variance | Nested continuum bias |
|---:|---:|---:|---:|---:|---:|
| .25 | 0 | 96.60% | 44.00% | 0.2635 | 0.0041% |
| .75 | 0 | 29.70% | 7.23% | 0.1632 | 0.0425% |
| .25 | .5 | 119.30% | 74.25% | 0.1885 | 1.2362% |
| .75 | .5 | 38.70% | 18.20% | 0.2157 | 1.2244% |
| .25 | -.5 | 94.47% | 88.35% | 0.1440 | 1.2801% |
| .75 | -.5 | 26.75% | 22.88% | 0.2096 | 1.2379% |

The flat cases' absolute errors decrease from 0.67287 to 0.33864 and from
0.35830 to 0.24102. Relative error at zero remains undefined. All errors
include all ten tidal components in the shared internal frame modulo the
remaining reflection, not just eigenvalues or a fitted wave amplitude.

## Coverage warning and limit of the result

Both arms satisfy the unchanged 2% continuum-bias condition in every case.
The baseline nominal balls contain truth in eight cases, but the nested
balls contain truth in **seven of eight**. Case **6** (A=.25,w=-.5) has
absolute error **0.44176**, exceeding radius **0.41321**. Its point error
improves by only 6.5% even though its nominal maximum variance falls by 85.6%.
Keep this noncoverage; do not enlarge its radius or remove it after the fact.
Eight cases do not establish or refute 95% population coverage.

The three strong cases have point errors below 25%, but their confidence
radii are 0.459--0.538, exceeding the 0.375 relative-accuracy target radius
at truth norm 1.5. Thus **no 25% confidence-contained accuracy claim is made**.
Weak-case errors remain large. The objective is not completed by this pilot.

The data support retaining more relational information as a precision lever.
They do not yet support the new uncertainty construction's reliability.
The next step is independent coverage validation of this unchanged candidate,
including zero-clock, bias and fit failures. Any later accuracy study still
needs explicit regularity/bias scope, independent high-count calibration and
90% power sizing before its data. No unrestricted-frequency or global result
follows, and the old fifth-degree failure remains unchanged.

## Resource cost and audit

The total sample is **2,129,920 events plus 16 supplied tips**; events were
shared, not doubled. The nested fits use **3727--4185** interval clocks versus
**948--1089**. Querying both directional membership bits requires roughly
65.0--69.7 million landmark/event predicates per nested case versus 34.6
million for the baseline; these counts are not CPU action-evaluation counts.
The actual paired generator queries only the union once. Compressed order
archives are approximately 4.24--4.49 MB versus 2.58--2.60 MB.

Computational cost also increases. One sequential pure-estimator replay of
the same case 5 took **4.37 seconds** for the baseline and **60.03 seconds**
for the nested fit on this host. No concurrent research job was running in
that replay, but it is one descriptive timing, not a controlled benchmark.
The sum of pilot per-pair times, including generation, both finite fits and
both continuum evaluations under variable load, is **1116.27 seconds**.

All 32 fits have Weyl rank ten, full nuisance rank for their actual landmark
count, successful optimizers and no boundary flags. The audit verified all
11 historical source hashes, all 48 artifact hashes, all eight shared-order
and shared-measurement bit comparisons, preserved frame prefixes, and all
saved error/variance ratios. Separate pure-estimator replays of **both**
case-5 arms match every saved fit field exactly without generator/profile/
truth imports. The additional coordinates are inferred by each estimator.

[Recorded results and audit](nested_pilot_summary.json) retain every case,
including case 6. Raw files are in `outputs/intrinsic_weyl_nested_pilot_v1`;
the raw summary SHA256 is
`981b88fee30c04ec5e3ca223ffbe5f60d182f60968c4f7ec196a4c36e480ecb7`.
The pure-replay timing record is `pure_replay_audit.json` in that folder.
