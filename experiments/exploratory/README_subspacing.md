# Exploratory subspacing detector

Run from the repository root with the existing environment:

```powershell
E:/lab/causal-spacetime/.venv/Scripts/python.exe experiments/exploratory/run_subspacing_pilot.py
E:/lab/causal-spacetime/.venv/Scripts/pytest.exe tests/test_subspacing_detector.py
E:/lab/causal-spacetime/.venv/Scripts/ruff.exe check experiments/exploratory tests/test_subspacing_detector.py
E:/lab/causal-spacetime/.venv/Scripts/ruff.exe format --check experiments/exploratory tests/test_subspacing_detector.py
```

The JSON is `outputs/subspacing_pilot_summary.json`. A rerun preserves earlier
campaigns inside `previous_campaigns`. Each new campaign records source hashes,
seeds, actual split sizes, timings, spectra, calibration residuals and failures.
The runner budgets 285 seconds and records cells it does not attempt. The default
full grid includes both q modes and three background repetitions per Stage 2 cell.
Stage 2 uses N=2^17, 2^19, 2^21. Stage 3 uses N=2^20 and 2^22
at q=.25, 1, 4 for background r=128 and waves r=64, 96, 128, 192,
followed by five N=2^23 cells (q=1 background/64; q=4 background/64/96).
Each cell runs known and unknown q. Four worker processes overlap CPU fitting
with GPU counting; fixed job-index seeds and output ordering are preserved.
Unattempted cells remain explicit. Per-cell GPU pool bytes describe that worker
process and include reusable allocations.

## Boundary between truth and observation

`subspacing_wave_sim.py` owns physical coordinates, cumulative quadrature,
physical-volume sprinkling, continuum interval volumes and CUDA kernels.
`subspacing_observer.py` imports none of that code. Its `observe` entry point
receives event indices, a boolean strict-order callback and optionally an
accelerated interval-count callback. Coordinates in its background model and
calibration result are inferred nuisance parameters, not supplied truth.
Only Stage 1 and the unit tests use continuum truth to check equations.

The pilot, reference and probe pool are disjoint. Pilot counts alone select
anchors; independent reference counts calibrate them. All reference events
contribute to counts. The pilot uses min(16384,N/4) events; the rest of the
first half supplies probes, and the second half supplies references. Four compact
anchors per side select comparable probes across the entire pool in tiles of
65536. Up to 32768 comparable probes are counted. Up to ten ordered anchors per
side calibrate the background; auxiliary anchors never restrict probe support.
The core maximizes volume-profile spread times common pilot support cubed.
Auxiliary pairs also require pilot interval mass >.0005 to avoid empty counts.
Count and window filters further reduce the sample; all reductions are reported.

## Numerical and scientific limits

- The local Gauss polynomial primitive and CDF inversion approximate the exact
  integrals at floating-point accuracy. Agreement is tested against the existing
  DOP853 `WaveFixture`, independent volume formulas and analytic Jacobians.
  No uniform error bound is claimed as q approaches pi^2.
- The ellipse area is **2 pi / sqrt(det K)**. The slice volume after integrating
  the v width is **pi / sqrt(det K)**, which supplies the u marginal.
- GPU time for 10 million primitive evaluations excludes host-to-device input
  transfer; download time is reported separately. CUDA compilation uses CuPy's
  in-memory cache to avoid inaccessible user-cache writes in the sandbox.
- Calibration scans both O(3) chiralities on 64 Fibonacci directions, with
  separately ranked q=.05,.5,1.5,4,7 buckets. Gauge-equivalent duplicate seeds
  are removed. Population fits use up to 16 seeds per bucket (Stage 1: 64);
  finite-N fits use one per bucket to bound cost. None certifies a global optimum.
  Fits use inverse binomial-SD weights, but reported RMS is unweighted. Small q
  is ill-conditioned; at q=0 the preferred null axis is unidentifiable.
- The probe initializer solves the lifted linear system for (v,y1,y2,|y|^2)
  on 33 u values, followed by damped Gauss-Newton. Stage 1 checks every sampled
  active point with initial=None after aligning the transverse O(2) gauge.
- The background interval coefficient uses its exact trigonometric identity
  and a cancellation-free small-argument series, checked against hyp0f1.
  The standalone runner defaults BLAS to one thread for small calibration SVDs.
- Core bands are p in [.0005,.04], f in [.25,.80] and their late dual.
  Auxiliary bands extend to [.0005,.12] and [.12,.80]. Earlier campaigns remain
  in JSON history. These adaptations are exploratory, not a frozen design.
- The shear coefficient is the exact interval differential to endpoint STF
  shear. Factoring out a Fourier carrier uses the leading high-frequency
  endpoint symbol, rather than an exact finite-lambda likelihood template.
- Stage 1 uses exact population interval probabilities on a finite collection
  of interior test points. Its coordinate-noise diagnostic holds population
  anchors fixed and does not include anchor-estimation error.
- Finite-probe shot noise need not follow 1/(omega sqrt(n_R)). The runner measures
  slopes and records background peaks instead of asserting that noise law.
- Spectra divide the summed matched-filter numerator by the global weighted
  symbol norm. No per-point symbol is inverted. Both signed frequencies are
  scanned: transverse parity conjugates the circular symbol and swaps the two
  amplitudes. The reported envelope is invariant; the signed peak is gauge-dependent.
- Off-peak RMS excludes frequencies within 20 of the peak; its high-frequency
  variant additionally requires omega>=44. Target SNR uses the same exclusion
  around the injected frequency. Noise L2 slopes also use omega>=48, and the
  omega slope is compared with -1 (inverse frequency) and 0 (flat shot noise).
- Off-peak SNR, Z and two-block differences are exploratory diagnostics, not
  p-values. Blocks share the reference sample and anchor calibration. Peaks,
  q point errors and finite alpha>2 cases establish neither a detection error
  rate nor asymptotic reconstruction accuracy. No accuracy gate is passed by a
  sampled point estimate alone.

## Observed Step 2 campaign (2026-11-28)

RTX 5090, four worker processes: **228.621 s**; all 36 Stage 2
and 70 Stage 3 mode-specific results completed, with no failures or budget skips.
This meets the 285 s cap, but exceeds the suggested 140-220 s range.
Full spectra, seeds, source hashes, timing, memory-pool observations and previous
campaigns are retained in `outputs/subspacing_pilot_summary.json`. Earlier runs
were runtime diagnostics and may contain explicit budget skips; they are not
additional independent repetitions. The rejected cross-product GPU experiment
is not part of the delivered order kernel.

Validation: 35 tests passed, plus Ruff check and format check. All 657 active
diamond test points across q=0,.25,1,4 converged from initial=None; maximum
coordinate error was 1.471e-15. Stage 1 uses independently generated population
probabilities and fitted anchors: maximum coordinate error 2.973e-13
over all 204 active points at q=.05,.25,1,4, after transverse gauge alignment.
Maximum forward-map error: 6.523e-16.
The 96 overlapping results from the last runtime diagnostic have a maximum
amplitude difference of 0.000e+00 under four-process execution.

### Stage 2 measured slopes

| q | q mode | L2 vs n_ref | High L2 vs n_ref | Mean high amplitude vs omega |
|---:|:---|---:|---:|---:|
| 0 | known | -0.719 | -0.888 | -0.393 |
| 0 | unknown | -0.718 | -0.865 | -0.275 |
| 1 | known | -0.668 | -0.801 | -0.248 |
| 1 | unknown | -0.786 | -0.965 | -0.229 |

Each row uses nine observations (three repetitions at each N). The high-frequency
omega slopes are closer to flat shot noise (0) than inverse frequency (-1),
but neither a flat law nor an asymptotic n_ref law is established. Probe counts
also grow with N before reaching the cap; these slopes do not isolate reference noise.

### Anchor support and conditioning

| q | Comparable pool fraction | Known-q condition range | Unknown-q condition range |
|---:|---:|---:|---:|
| 0 | 3.28-7.79% | 1.1e+16-1.53e+16 | 1.48e+03-1.42e+16 |
| 0.25 | 2.84-6.00% | 2.01e+03-4.26e+03 | 2.55e+03-3.4e+04 |
| 1 | 2.89-6.93% | 440-2.09e+03 | 428-2.18e+03 |
| 4 | 2.06-4.93% | 107-953 | 134-792 |

The 150-600 condition range is not universal. Small q remains poorly conditioned;
q=0 has no identifiable preferred null axis. Conditions are from the weighted
calibration Jacobian, whereas calibration RMS remains unweighted.

### Stage 3: all observed target scores

Every pair below is **known / unknown q**. Full and high SNR use the injected
target frequency, not necessarily the observed global peak. Background rows
use r=128 only as a diagnostic target. Peak columns are unsigned envelope peaks.

| N | q | Scene | r | Target SNR | Target SNR high | Global peak | Unknown q_hat |
|---:|---:|:---|---:|---:|---:|---:|---:|
| 2^20 | 0.25 | background | 128 | 0.50 / 0.73 | 0.73 / 1.25 | 24 / 24 | 0.2036 |
| 2^20 | 0.25 | circular_wave | 64 | 0.69 / 0.78 | 1.72 / 1.83 | 24 / 24 | 0.2104 |
| 2^20 | 0.25 | circular_wave | 96 | 0.58 / 0.15 | 0.76 / 0.22 | 24 / 24 | 0.1903 |
| 2^20 | 0.25 | circular_wave | 128 | 0.53 / 0.29 | 0.84 / 0.59 | 24 / 32 | 0.2774 |
| 2^20 | 0.25 | circular_wave | 192 | 0.35 / 0.37 | 0.63 / 0.60 | 24 / 24 | 0.2991 |
| 2^20 | 1 | background | 128 | 1.16 / 1.25 | 1.53 / 1.63 | 24 / 32 | 0.9118 |
| 2^20 | 1 | circular_wave | 64 | 0.75 / 1.22 | 1.40 / 1.91 | 24 / 24 | 0.8989 |
| 2^20 | 1 | circular_wave | 96 | 1.27 / 1.50 | 2.03 / 1.89 | 24 / 24 | 0.9521 |
| 2^20 | 1 | circular_wave | 128 | 1.01 / 0.38 | 1.02 / 0.46 | 24 / 44 | 1.1473 |
| 2^20 | 1 | circular_wave | 192 | 0.29 / 0.29 | 0.49 / 0.53 | 28 / 24 | 0.9489 |
| 2^20 | 4 | background | 128 | 1.12 / 1.09 | 1.57 / 1.51 | 24 / 28 | 3.9952 |
| 2^20 | 4 | circular_wave | 64 | 1.05 / 0.99 | 2.30 / 2.10 | 24 / 24 | 3.9954 |
| 2^20 | 4 | circular_wave | 96 | 0.85 / 1.16 | 1.05 / 1.40 | 24 / 24 | 3.9772 |
| 2^20 | 4 | circular_wave | 128 | 0.67 / 0.59 | 0.86 / 0.90 | 24 / 24 | 4.0275 |
| 2^20 | 4 | circular_wave | 192 | 0.60 / 0.56 | 1.03 / 0.93 | 24 / 24 | 3.9820 |
| 2^22 | 0.25 | background | 128 | 0.69 / 1.00 | 0.99 / 1.36 | 32 / 24 | 0.1991 |
| 2^22 | 0.25 | circular_wave | 64 | 1.34 / 0.77 | 1.93 / 1.36 | 24 / 24 | 0.2271 |
| 2^22 | 0.25 | circular_wave | 96 | 0.40 / 0.53 | 0.75 / 0.91 | 28 / 28 | 0.2637 |
| 2^22 | 0.25 | circular_wave | 128 | 0.45 / 0.32 | 0.57 / 0.45 | 28 / 28 | 0.2731 |
| 2^22 | 0.25 | circular_wave | 192 | 0.20 / 0.22 | 0.42 / 0.49 | 24 / 24 | 0.2436 |
| 2^22 | 1 | background | 128 | 0.64 / 0.71 | 0.78 / 1.00 | 24 / 24 | 0.9581 |
| 2^22 | 1 | circular_wave | 64 | 1.95 / 1.26 | 3.80 / 3.12 | 24 / 24 | 1.0179 |
| 2^22 | 1 | circular_wave | 96 | 0.98 / 1.00 | 1.21 / 1.28 | 24 / 24 | 1.0056 |
| 2^22 | 1 | circular_wave | 128 | 0.39 / 0.54 | 0.65 / 0.75 | 24 / 24 | 1.0016 |
| 2^22 | 1 | circular_wave | 192 | 0.27 / 0.25 | 0.50 / 0.47 | 24 / 24 | 0.9955 |
| 2^22 | 4 | background | 128 | 0.23 / 0.24 | 0.48 / 0.49 | 24 / 24 | 3.9995 |
| 2^22 | 4 | circular_wave | 64 | 4.01 / 3.35 | 8.16 / 6.99 | 60 / 24 | 3.9907 |
| 2^22 | 4 | circular_wave | 96 | 0.42 / 0.52 | 0.79 / 0.94 | 28 / 32 | 4.0904 |
| 2^22 | 4 | circular_wave | 128 | 1.58 / 1.62 | 2.18 / 2.16 | 24 / 24 | 4.0020 |
| 2^22 | 4 | circular_wave | 192 | 0.72 / 0.59 | 1.33 / 1.38 | 24 / 24 | 4.0209 |
| 2^23 | 1 | background | 128 | 0.34 / 0.32 | 0.63 / 0.64 | 36 / 32 | 0.9840 |
| 2^23 | 1 | circular_wave | 64 | 1.99 / 1.99 | 5.15 / 5.14 | 24 / 24 | 0.9991 |
| 2^23 | 4 | background | 128 | 0.45 / 0.42 | 0.80 / 0.79 | 24 / 28 | 3.9937 |
| 2^23 | 4 | circular_wave | 64 | 6.29 / 8.45 | 10.98 / 14.55 | 60 / 64 | 3.9935 |
| 2^23 | 4 | circular_wave | 96 | 1.70 / 1.50 | 4.27 / 3.74 | 24 / 24 | 3.9909 |

At N=2^23, q=4, r=64 is the strongest listed response (high target SNR
10.98/14.55; peaks 60/64). At q=1,r=64 and q=4,r=96 the high target scores
are elevated, but their full-band global peaks remain at 24. Those cells do
not establish blind recovery of the injected frequency.

The r=128/192 and q=.25 cases remain inconsistent or weak. For example, at
N=2^22,q=1 the r=128/192 high target SNRs are .65/.75 and .50/.47;
at q=.25 they are .57/.45 and .42/.49. q=4,r=128 reaches 2.18/2.16,
but this is a single exploratory realization, not a certified detection.
No false-positive rate, power, confidence-interval accuracy gate or
asymptotic detectability claim follows from these scores.
