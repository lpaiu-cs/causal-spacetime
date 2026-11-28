# Wider central landmarks: the tested lever did not work

2026-09-22. The [paired protocol](design_pilot_protocol.md) was frozen at
`1aeb142651050e36db8134c537c0dae136bec573` before fresh data. All eight
pairs completed with 64 landmarks per arm and no selection, empty-clock,
rank, boundary or optimizer failures. Both arms saw exactly the same
4096-event selection order and 262144-event measurement stream per case.

**Do not adopt the wider selector as an improvement.** Magnetic variance
increased in five of eight pairs, with median ratio **1.0941**. The largest
full-tidal covariance eigenvalue had median ratio **1.1448**. Actual tensor
error increased in seven pairs, with median ratio **1.1189**. These are
descriptive paired findings; the study was not powered to establish a
population improvement or an accuracy guarantee.

## Complete paired outcomes

All ratios are `wide/connected`; below one favors the wider arrangement.
Errors compare all ten signed Weyl components in the SAME internal frame,
modulo the common remaining reflection.

| Case | A | w | Magnetic variance ratio | Largest variance ratio | Actual tensor-error ratio | Continuum relative bias, connected -> wide |
|---|---:|---:|---:|---:|---:|---|
| 0 | 0 | 0 | .7757 | .8859 | 1.0938 | numerical zero -> zero |
| 1 | 0 | .5 | 1.0974 | 1.1991 | 1.1441 | numerical zero -> zero |
| 2 | .25 | 0 | 1.1645 | 1.3378 | .7604 | .2673% -> .3907% |
| 3 | .75 | 0 | 1.2277 | 1.3277 | 1.0042 | 1.2208% -> 1.5844% |
| 4 | .25 | .5 | .6993 | .6897 | 1.0303 | 4.9067% -> 3.6180% |
| 5 | .75 | .5 | .8490 | .8994 | 1.3460 | 5.2353% -> 6.0080% |
| 6 | .25 | -.5 | 1.1939 | 1.1334 | 1.5177 | 3.0489% -> 3.5557% |
| 7 | .75 | -.5 | 1.0908 | 1.1561 | 1.4018 | 4.3372% -> 5.3353% |

Electric variance improved in seven pairs, but only modestly: median ratio
.9562. This does not address the dominant magnetic uncertainty. All four
varying-profile cases in BOTH arms violated the unchanged 2% continuum-bias
allowance. Their broad nominal balls nevertheless contained the targets;
that is not transferred coverage or useful-precision recovery.

Expanding a coordinate-proxy radius does increase the raw magnetic quartic
response. The experiment shows that this alone is insufficient after the
changed clock set, correlations and unknown positions are accounted for.
It does not disprove every possible geometric/information-based design.
The tested `.30` selector is retained only for reproduction, and is not
the default. No new radius is tuned against these same outcomes.

## Verification and provenance

Raw files are in `outputs/intrinsic_weyl_design_pilot_v1`. The committed
`design_pilot_summary.json` is the complete summary with LF line endings.
Raw summary SHA-256:
`0e3aff4e3c3c6e0adb9533910949f6400034296cf57b33d7dbcbce98b8717f6e`.
LF summary SHA-256:
`25c8533fe51a638507bedc913a5adc27edd7adf5e7af931793ce9bfd843fe356`.

The post-run check verified all **48 artifact hashes** and **10 source
hashes** against the generating commit and the unchanged source files.
All paired selection orders and first two frame IDs matched. Measurement
query bytes matched on all **362 shared landmarks across the eight pairs**.
The oracle canonical tensors agreed to at most `1.23e-15` between arms.
Variance traces and reported ratios were independently recalculated.
The pure estimator replay of case 7, `wide`, reproduced coefficients,
coordinates and covariance exactly (maximum difference **0.0**) without
opening a truth archive or importing a generator.

Python 3.14.0, NumPy 2.5.3 and SciPy 1.18.1 were used. Sum of case runtimes
was 304.19 seconds. There were 2,129,920 distinct sampled events plus tips;
the common measurement stream is not counted twice.

## What changes next

The radius-only noise lever is rejected for this implementation. A separate
[gradient diagnosis](gradient_model_diagnosis.md) identifies a specific
omitted term behind much of the varying-profile bias and checks its local
information cost without supplying a wave profile. The next substantive
step is a blind joint fit of curvature and its allowed first derivatives,
followed by fresh validation and powered precision. Neither that fit nor
the ultimate accuracy claim has been completed by this comparison.
