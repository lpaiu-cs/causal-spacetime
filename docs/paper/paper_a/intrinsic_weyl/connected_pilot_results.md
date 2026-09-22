# Connected selection is retained; Weyl accuracy remains open

2026-09-22. All six cases of the [fixed connected pilot](connected_pilot_protocol.md)
ran from clean pre-data commit `efa8624`, with no seed replacement or rule
changes. Each used 4096 selection events and 262144 measurement events.
The complete records are [the result](connected_pilot_result.json) and
[the retention and spectra audit](connected_retention_audit.json).

Every case retained **21 past, 22 central, and 21 future landmarks**. The
minimum central-to-past and central-to-future qualified degrees were,
respectively, `(8,9), (8,7), (9,7), (7,6), (7,7), (8,9)`. Thus the specific
central-retention failure in the preceding pilot is repaired in these six
new samples. This is an observed instrument-construction result, not a
population reliability or accuracy statement. Selection events increased
fourfold; the measurement count was unchanged.

## All observed spectra

Evaluation truth for BOTH spectra is `(-A,0,A)`. The estimator receives
neither that form nor the amplitude.

| Case | A | Electric eigenvalues | Magnetic eigenvalues |
|---|---:|---|---|
| 0 | 0 | (-0.1139, -0.0291, 0.1429) | (-0.2112, -0.1150, 0.3262) |
| 1 | 0 | (-0.2126, 0.0192, 0.1934) | (-0.2932, -0.1596, 0.4528) |
| 2 | 0.25 | (-0.4206, 0.0566, 0.3641) | (-0.4073, -0.2308, 0.6381) |
| 3 | 0.25 | (-0.2001, -0.0518, 0.2519) | (-0.2654, -0.1626, 0.4280) |
| 4 | 0.75 | (-0.6354, 0.0366, 0.5988) | (-0.6879, -0.0626, 0.7505) |
| 5 | 0.75 | (-0.8604, 0.0131, 0.8472) | (-0.6583, 0.1402, 0.5181) |

All fits reported optimizer success, full coordinate rank, projected Weyl
rank ten and no inferred-time boundary excess. The generation-and-fit times
were 6.0--8.0 seconds per case; all ambiguity escalation counts were zero.
Some measurement intervals contained only **8--25 events** at each case's
minimum, so the local Gaussian covariance approximation needs explicit
coverage checking. The flat spectra and the weak-wave magnetic spectra
still show substantial noise. No accuracy or improvement verdict is assigned.

All six source/archival hashes were checked against the reachable pre-data
commit. Retention and minimum central degrees were independently recomputed
from each saved selection order. A fresh estimator process replayed case 4
to `1e-11` in all ten coefficients with no generator imports. Its archive
contains only the previously specified six allowed fields.

## Full-tensor audit, not just separate spectra

Correct separate E and B eigenvalues would not suffice: their relative
directions also matter. The post-pilot evaluator now supports comparison of
all ten components against the true tensor in an aligned internal frame.
It also reuses the continuum-volume audit on these SAME selected landmarks.
These are labelled truth-using diagnostics, never estimator inputs or a
second independent recovery result. The earlier audit record is preserved.

For the constant Brinkmann metric, the midpoint is `(u,v,x,y)=(1/2,0,0,0)`.
Write a tangent as `(U,V,X,Y)`, with geodesic parameter s in [0,1]. The
transverse geodesics obey

```text
x(s) = X s sinhc(sqrt(A) U s),
y(s) = Y s sinc(sqrt(A) U s),
v(1) = V - U A integral_0^1 (x(s)^2 - y(s)^2) ds.
```

Here sinc(z)=sin(z)/z and sinhc(z)=sinh(z)/z, with their continuous values
at zero. Thus endpoint data determine U, X, Y and V. Dimensionless normal
coordinates in the tip observer's tetrad are
`((U-V)/2, X/sqrt(2), Y/sqrt(2), (U+V)/2)`. The implementation evaluates
the smooth integral by 32-point Gauss quadrature. An independent geodesic
ODE round trip and an action-based squared-distance check cover both signs
of U, U=0, flat space and two nonzero amplitudes.

The evaluator aligns true and inferred spatial clouds by orthogonal
Procrustes. It permits O(3), including reflections: causal order does not
provide a preferred spatial orientation. Transforming the full covariant
tensor preserves the corresponding magnetic parity automatically. With
`C_uiuj=-diag(A,-A)` and squared tip time 2, the error is
`sqrt(||E_hat-E||_F^2 + ||B_hat-B||_F^2)`; its relative version divides by
`2A` for nonzero A. This evaluates signed components and relative direction,
not merely two eigenvalue lists. It is a point-error diagnostic, not a
confidence bound.

The audit is invoked explicitly with input and NEW output paths; it refuses
to overwrite earlier audit files. Independent accuracy/coverage validation,
truncation control beyond these geometries, and varying-profile finite-order
tests remain outstanding. The active goal is not complete.
