# Preregistration S2: Two-Component Huber-Whittle Confirmatory Campaign

**Status:** Frozen post-diagnostic continuation of S1 (`P2 -> F2 -> R2` provenance chain)  
**Target Branch:** `codex/subspacing-order-detector`  
**Natural Units:** `c = 1` throughout; strictly blind order-only observer (`q = None`).

## 1. Motivation from Audited Campaign S1 (Commit `36e551a`)

In the audited confirmatory run **S1** (`f88d637` / `36e551a`), the single-parameter
`1 + c_1 / omega^2` Huberized Whittle background model passed `16/18` frozen gates
(including `35/36` detection and `36/36` frequency lock on `W1`, `33/36` detection
and `36/36` frequency lock on Poisson/random-phase `W1p`, `27/32` detection and
`32/32` lock on `W2`, `20/20` lock on `W3`, and all `G7a/b/c` crossover checks),
but **failed two gates (`G1`, `G6a`)**:

1. **Low-frequency background tail accumulation (`G1`):** Across all `60` null runs
   in `S1`, `36/60` (`60.0%`) peaked in the lowest bins `omega <= 32` (`5.1%` of bins),
   producing `4/24` false positives at `omega in {24, 28}` in `N23` (`4/48` pooled,
   Wilson upper95 `= 0.1955 > 0.18`).
2. **Inflation of `W3` baseline (`G6a`):** The upward `N23` null drift (`mean = 2.686`
   vs `2.402` in `N22`) shrank the `W3 - N23` separation CI lower bound to `+0.319 < +0.40`.

Diagnostic inspection on disjoint exploratory seeds (`759000+`) established that smooth
window-tapered triangulation drift falls faster than `omega^-2` (roughly `omega^-3` to
`omega^-4`) across `omega in [24, 64]`. Upgrading `spectrum()` in
`experiments/exploratory/subspacing_observer.py` to the two-component scale-invariant
Huberized Whittle model

    sigma_bg(omega)^2 = var_shot + c_1 / omega^2 + c_2 / omega^4   (c_1, c_2 >= 0)

with multi-start analytic-gradient L-BFGS-B + boundary minimization eliminated
low-frequency peak pile-up on `44/44` fresh exploratory nulls (`0/44` peaks at
`omega <= 32`; `mean Z_max = 2.008`, `sd = 0.202`, `q95 = 2.277`, `max = 2.584`).

Rather than altering S1 retroactively, S1 remains preserved verbatim at `36e551a` and
**S2** executes a completely new preregistered campaign in fresh seed decade
`52_000_000..52_999_999`.

## 2. Frozen Instrument and Protocol Constants

| Parameter | Frozen S2 Value |
|---|---|
| Background model in `spectrum()` | `var_shot + c_1 / omega^2 + c_2 / omega^4` (non-negative, scale-invariant) |
| Frequency scan `FREQUENCIES` | `np.arange(24.0, 257.0, 4.0)` (59 uniform bins) |
| Primary detection statistic | `whitened_peak_z` (`Z_max`, Eq. (11)) |
| Primary threshold `Z_TH` | `2.50` (positive iff `Z_max >= 2.50`; `1/44 = 2.27%` in exploratory Step 2b) |
| Wavenumber lock window | `abs(whitened_peak_frequency - r) <= 8.0` (5/59 = 8.47% chance) |
| Calibration escalation | `calibrate_progressive` (`n_ref = len(references)`, multi-start `4` + cores `5, 4`) |
| Split & cap | `N_pilot = min(16384, N // 4)`, `N_ref = N // 2`, `MAX_PROBES = 65536`, `q = None` |
| Weyl-jet crossover baseline | `intrinsic_weyl_joint.fit_joint(..., starts=1, model_order=5, clock_mask=two)` at `N = 2^21`, `q = 0.25` |

## 3. Disjoint S2 Confirmatory Seed Allocation (`52_000_000..52_999_999`)

Certified disjoint via `seed_windows.assert_windows_disjoint_and_fresh` against
`0..999_999` (exploratory), `1_000_000..49_999_999` (P11--P14), and
`51_000_000..51_999_999` (S1):

| Arm | Model / Parameters | Count `n` | Seed Block (Inclusive) | Role |
|---|---|---:|---|---|
| `N22` | `background`, `N=2^22`, `q=4.0` | `28` | `52_000_000..52_000_027` | Null control (`2^22`) |
| `N23` | `background`, `N=2^23`, `q=4.0` | `28` | `52_100_000..52_100_027` | Null control (`2^23`) |
| `W1` | `circular_wave`, `N=2^22`, `q=4.0`, `r=64`, `phi=0.37` | `36` | `52_200_000..52_200_035` | Primary sub-discreteness (`alpha=3.667`) |
| `W1p` | `circular_wave`, `N~Poisson(2^22)`, `q=4.0`, `r=64`, `phi~U(0,2pi)` | `36` | `52_300_000..52_300_035` (+ phase `52_310_000..52_310_035`) | Poisson & random phase |
| `W2` | `circular_wave`, `N=2^23`, `q=4.0`, `r=96`, `phi=0.37` | `32` | `52_400_000..52_400_031` | Deeper sub-discreteness (`alpha=3.493`) |
| `W3` | `circular_wave`, `N=2^23`, `q=4.0`, `r=128`, `phi=0.37` | `28` | `52_500_000..52_500_027` | Deep frontier (`alpha=3.286`, full TPR + lock) |
| `X_null_q1` | `background`, `N=2^22`, `q=1.0` | `16` | `52_600_000..52_600_015` | Crossover null (`q=1.0`) |
| `X_sub_q4_r32` | `circular_wave`, `N=2^22`, `q=4.0`, `r=32`, `phi=0.37` | `24` | `52_610_000..52_610_023` | Crossover short wave (`q=4, alpha=4.40`) |
| `X_sub_q1_r48` | `circular_wave`, `N=2^23`, `q=1.0`, `r=48`, `phi=0.37` | `24` | `52_620_000..52_620_023` | Crossover intermediate (`q=1, alpha=4.12`) |
| `X_weyl_r4` | Weyl jet `N=2^21`, `q=0.25`, `r=4` | `8` | `52_630_000..52_630_007` (+ lm `52_635_000..52_635_007`) | Long-wave Weyl jet accuracy |
| `X_weyl_r32` | Weyl jet `N=2^21`, `q=0.25`, `r=32` | `8` | `52_640_000..52_640_007` (+ lm `52_645_000..52_645_007`) | Short-wave Weyl jet collapse |

## 4. Two-Sided 95% CI Gates and >=90%-Power Sizing (`AGENTS.md` Section 6)

| Gate | Estimand | Frozen Two-Sided 95% CI Pass Criterion | Sizing Basis (Step 2b Exploratory) | Power |
|---|---|---|---|---:|
| **G1** | Pooled Null FPR `P(Z >= 2.50)` (`N22+N23`, `n=56`) | Wilson upper `<= 0.18` (`<= 4/56`) | `1/44=0.023`; conservative `p=0.03` | **97.4%** |
| **G1b** | Pooled Null low-freq pile-up `P(pf <= 32)` (`n=56`) | Wilson upper `<= 0.20` (`<= 5/56`) | Pilot `0/44`; chance `3/59=0.051` | **94.0%** |
| **G2a** | Null calibration completion & residual (`n=56`) | `56/56` valid, `c_ref > 0`, `cal_rms` upper95 `<= 1.2e-4` | Pilot mean `4.4e-5` | **>99.9%** |
| **G2b** | Null blind `mean(|q_hat - 4|)` (`N22`, `N23`) | `N22` upper95 `<= 0.035`; `N23` upper95 `<= 0.025` | Pilot `N22`: `0.016`, `N23`: `0.007` | **>99.9%** |
| **G3a** | `W1` (`r=64`, `2^22`, `n=36`) TPR | Wilson lower `>= 0.50` (`>= 24/36`) | Conservative `0.80` (`6/6` pilot) | **98.2%** |
| **G3b** | `W1` `Z_max` sep vs `N22` | Welch 95% CI lower `>= +0.65` | Pilot `+1.29`, `SE=0.10` (`> 6 sigma`) | **>99.9%** |
| **G3c** | `W1` freq lock & `q` accuracy | Lock lower `>= 0.65` (`>= 30/36`), `|q_err|` upper `<= 0.030` | Pilot `6/6`, `q_err=0.011` | **93.7%+** |
| **G4a** | `W1p` (`Poisson`, random `phi`, `n=36`) TPR | Wilson lower `>= 0.50` (`>= 24/36`) | Conservative `0.80` (`6/6` pilot) | **98.2%** |
| **G4b** | `W1p` `Z_max` sep vs `N22` | Welch 95% CI lower `>= +0.50` | Pilot `+0.98`, `SE=0.07` | **>99.9%** |
| **G4c** | `W1p` freq lock & `q` accuracy | Lock lower `>= 0.65` (`>= 30/36`), `|q_err|` upper `<= 0.030` | Pilot `6/6`, `q_err=0.010` | **93.7%+** |
| **G5a** | `W2` (`r=96`, `2^23`, `n=32`) TPR | Wilson lower `>= 0.50` (`>= 22/32`) | Conservative `0.85` (`6/6` pilot) | **99.5%** |
| **G5b** | `W2` `Z_max` sep vs `N23` | Welch 95% CI lower `>= +1.20` | Pilot `+2.68`, `SE=0.19` | **>99.9%** |
| **G5c** | `W2` freq lock & `q` accuracy | Lock lower `>= 0.65` (`>= 27/32`), `|q_err|` upper `<= 0.025` | Pilot `6/6`, `q_err=0.010` | **90.6%+** |
| **G6a** | `W3` (`r=128`, `2^23`, `n=28`) frontier TPR | Wilson lower `>= 0.50` (`>= 20/28`) | Conservative `0.80` (`6/6` pilot) | **91.0%** |
| **G6b** | `W3` `Z_max` sep vs `N23` | Welch 95% CI lower `>= +0.90` | Pilot `+1.67`, `SE=0.135` | **>99.9%** |
| **G6c** | `W3` freq lock & `q` accuracy | Lock lower `>= 0.60` (`>= 22/28`), `|q_err|` upper `<= 0.025` | Pilot `6/6` (`98.2%` at `0.90`) | **98.2%** |
| **G7a** | Crossover Weyl `r=4` recovery (`n=8`) | Student-`t` lower95 `>= 0.65` | S1 held-out `0.738` | **>99.9%** |
| **G7b** | Crossover Weyl `r=32` collapse (`n=8`) | Student-`t` upper95 `<= 0.10` | S1 held-out `0.010` | **>99.9%** |
| **G7c** | Crossover Subspacing `X_sub_q4_r32` (`n=24`) vs `N22` | TPR lower `>= 0.50`, Lock lower `>= 0.55`, Welch diff lower `>= +0.55` | Pilot `mean=3.14` (`3/3 >= 2.91`) | **91.1%+** |
| **G7d** | Crossover Subspacing `X_sub_q1_r48` (`n=24`) vs `X_null_q1` | TPR lower `>= 0.50`, Lock lower `>= 0.55`, Welch diff lower `>= +0.40` | Pilot `mean=2.83` (`3/3 >= 2.69`) | **91.1%+** |

**Overall Campaign Verdict:** `CONFIRMED` iff all `20` gates `G1..G7d` pass simultaneously on their frozen two-sided 95% CIs (`MIXED_OR_FAIL` otherwise; zero retroactive alterations).
