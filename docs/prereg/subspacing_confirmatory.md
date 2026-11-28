# Preregistration S1: Held-Out Confirmatory Sub-Discreteness Detection and Wavelength Crossover

**Status:** Frozen prior to confirmatory execution (`P -> F -> R` provenance chain)  
**Target Branch:** `codex/subspacing-order-detector`  
**Natural Units:** `c = 1` throughout; causal order relations only (`q = None`) for the spectral detector.

## 1. Scientific Claims and Scope

This preregistration tests whether the blind, unknown-`q` anchor-triangulated spectral detector (`experiments/exploratory/subspacing_observer.py`) operationally detects coherent circularly polarized vacuum plane gravitational waves below the discrete sprinkling spacing `ell = N^{-1/4}` (`alpha = log(N) / log(r) < 4`) without false-positive inflation, without prior knowledge of focusing `q`, carrier phase `phi`, or probe coordinates, and with a quantitative crossover against the local fifth-degree Weyl jet (`experiments/positive_control/intrinsic_weyl_joint.py`).

Specifically, S1 evaluates:
1. **Null calibration & false-positive bound (H1--H2):** Background controls (`N = 2^22` and `N = 2^23`, `q = 4.0`) control false positives of `Z_max >= Z_th` and recover blind background focusing `q_hat`.
2. **Primary sub-discreteness detection (`N = 2^22`, `alpha = 3.6667 < 4`) (H3):** Coherent waves at `r = 64` achieve high single-run TPR, clean mean `Z_max` separation over matched `N = 2^22` background controls, accurate peak frequency lock (`|f_peak - r| <= 8`), and blind `q` accuracy.
3. **Poisson count & random-phase invariance (`Poisson(2^22)`, `phi ~ U[0, 2pi)`) (H4):** Allowing total event counts `N` to fluctuate as `Poisson(2^22)` and drawing carrier phases independently across `[0, 2pi)` preserves detection, separation, frequency lock, and `q` accuracy.
4. **Deeper sub-discreteness detection (`N = 2^23`, `alpha = 3.4928 < 4`) (H5):** At `r = 96`, blind order-only detection remains clearly positive with `q` and frequency accuracy.
5. **Deep frontier ensemble signal (`N = 2^23`, `alpha = 3.2857 < 4`) (H6):** At `r = 128` (where single-run `Z_th` exceedance sits near 50% in the Step 2 pilot), the ensemble mean `Z_max` remains elevated above matched `N = 2^23` background controls and the global spectral peak locks onto `r = 128` far above chance (`5/59 = 8.47%`). We explicitly pre-specify that `W3` tests ensemble elevation and wavenumber locking rather than high single-run exceedance.
6. **Held-out finite-`N` Weyl-jet crossover (`r = 4` vs `r = 32`) (H7):** Local fifth-degree Weyl jets (`N = 2^21`, `q = 0.25`) recover `> 0.65` of crest electric Weyl curvature at `r = 4` and collapse `< 0.10` at `r = 32`, whereas the blind spectral detector (`N = 2^22`, `q = 1.0`) switches on strongly at `r = 32` relative to matched `q = 1.0` background controls.

**What S1 does NOT claim:** S1 does not claim to reconstruct the full 10-component Weyl tensor from the spectral detector, does not address non-stationary bursts or non-plane-wave geometries without separate calibration, and does not break the asymptotic `alpha -> 2^+` phase-decoherence barrier documented in Section 3.

## 2. Frozen Instrument and Protocol Constants

All constants are frozen before confirmatory runs begin:

| Parameter | Frozen Value |
|---|---|
| Frequency scan `FREQUENCIES` | `np.arange(24.0, 257.0, 4.0)` (59 uniform bins) |
| Primary detection statistic | `whitened_peak_z` (`Z_max`, Eq. (11)) |
| Primary threshold `Z_TH` | `3.30` (positive iff `Z_max >= 3.30`; 0/44 Step 2 null FP) |
| Wavenumber lock window | `abs(whitened_peak_frequency - r) <= 8.0` (5 of 59 bins = 8.47%) |
| Pilot split `N_pilot` | `min(16384, N // 4)` |
| Reference split `N_ref` | `N // 2` |
| Probe pool cap `MAX_PROBES` | `65536` |
| Focusing parameter `q` passed to detector | `None` (strictly blind) |
| Calibration escalation | `calibrate_progressive` (`n_ref = len(references)`, multi-start `4` + cores `5, 4`) |
| Weyl-jet crossover baseline | `intrinsic_weyl_joint.fit_joint(..., starts=1, model_order=5, clock_mask=two)` at `N = 2^21`, `q = 0.25` |

## 3. Strictly Disjoint Confirmatory Seed Allocation (`51_000_000..51_999_999`)

Enforced via `seed_windows.assert_windows_disjoint_and_fresh` against all house spent ledgers (`1_000_000..39_999_999`), P14 (`40_000_000..49_999_999`), and Step 1--3 exploratory seeds (`0..999_999`, including Step 2 `750_000+` and Step 3 `760_000+`):

| Arm | Model / Parameters | Count `n` | Seed Block (Inclusive) | Role |
|---|---|---:|---|---|
| `N22` | `background`, fixed `N=2^22`, `q=4.0` | `24` | `51_000_000..51_000_023` | Null control (`2^22`) |
| `N23` | `background`, fixed `N=2^23`, `q=4.0` | `24` | `51_100_000..51_100_023` | Null control (`2^23`) |
| `W1` | `circular_wave`, fixed `N=2^22`, `q=4.0`, `r=64`, `phi=0.37` | `36` | `51_200_000..51_200_035` | Primary sub-discreteness (`alpha=3.667`) |
| `W1p` | `circular_wave`, `N~Poisson(2^22)`, `q=4.0`, `r=64`, `phi~U(0,2pi)` | `36` | `51_300_000..51_300_035` (N/sprinkle), `51_310_000..51_310_035` (phase) | Poisson-$N$ & phase invariance |
| `W2` | `circular_wave`, fixed `N=2^23`, `q=4.0`, `r=96`, `phi=0.37` | `32` | `51_400_000..51_400_031` | Deeper sub-discreteness (`alpha=3.493`) |
| `W3` | `circular_wave`, fixed `N=2^23`, `q=4.0`, `r=128`, `phi=0.37` | `20` | `51_500_000..51_500_019` | Frontier ensemble elevation (`alpha=3.286`) |
| `X_null_q1` | `background`, fixed `N=2^22`, `q=1.0` | `12` | `51_600_000..51_600_011` | Crossover null (`q=1.0`) |
| `X_sub_q1_r32` | `circular_wave`, fixed `N=2^22`, `q=1.0`, `r=32`, `phi=0.37` | `16` | `51_610_000..51_610_015` | Crossover spectral switch-on (`alpha=4.40`) |
| `X_weyl_r4` | Weyl jet `N=2^21`, `q=0.25`, `r=4` | `8` | `51_620_000..51_620_007` (+ lm `51_625_000..51_625_007`) | Long-wave Weyl jet accuracy |
| `X_weyl_r32` | Weyl jet `N=2^21`, `q=0.25`, `r=32` | `8` | `51_630_000..51_630_007` (+ lm `51_635_000..51_635_007`) | Short-wave Weyl jet collapse |

## 4. Two-Sided 95% CI Gates and Pre-Data 90%-Power Sizing (`AGENTS.md` Section 6)

Per house policy (`AGENTS.md`), **every** accuracy, rate, separation, and calibration-error gate is decided strictly by its **two-sided 95% confidence interval** (`Wilson score CI` with `z = Phi^{-1}(0.975)` for binomial rates; `Welch t-CI` with Satterthwaite df for two-sample mean differences; `Student t-CI` for one-sample mean accuracy/error).

| Gate | Estimand | Frozen Pass Criterion (Two-Sided 95% CI) | Sizing Assumption (Step 2 Pilot) | Pre-Data Power |
|---|---|---|---|---:|
| **G1** | Pooled Null FPR `P(Z_max >= 3.30)` (`N22+N23`, `n=48`) | Wilson upper `<= 0.18` (needs `<= 3/48`) | Pilot `0/44`; conservative true `p_0 = 0.025` (at `0.03`: `94.5%`) | **96.8%** |
| **G2a** | Null calibration completion (`N22+N23`, `n=48`) | `48/48` valid fits, `mean(cal_rms)` upper95 `<= 1.2e-4` | Pilot mean `~3.8e-5` | **>99.9%** |
| **G2b** | Null blind `mean(|q_hat - q|)` (`N22` and `N23`) | `N22` upper95 `<= 0.035`; `N23` upper95 `<= 0.025` | Pilot `N22`: `0.013`, `N23`: `0.010` | **>99.9%** |
| **G3a** | `W1` (`r=64`, `2^22`) single-run TPR | Wilson lower `>= 0.50` (needs `>= 24/36`) | Conservative `p = 0.75` (`90.8%`); pilot `0.8125` | **98.9%** |
| **G3b** | `W1` mean `Z_max` separation over `N22` | Welch 95% CI lower `>= +1.00` | Pilot diff `+2.15`, `SE = 0.175` (`> 6 sigma`) | **>99.9%** |
| **G3c** | `W1` frequency lock `P(|f_peak - 64| <= 8)` | Wilson lower `>= 0.65` (needs `>= 29/36`) | Conservative `0.90` (`98.0%`); pilot `16/16` | **>98.0%** |
| **G3d** | `W1` blind `mean(|q_hat - 4|)` | Student-`t` upper95 `<= 0.030` | Pilot `0.010` | **>99.9%** |
| **G4a** | `W1p` (`Poisson`, random `phi`) TPR | Wilson lower `>= 0.50` (needs `>= 24/36`) | Conservative `0.75` (`90.8%`); pilot `0.8125` | **98.9%** |
| **G4b** | `W1p` mean `Z_max` separation over `N22` | Welch 95% CI lower `>= +1.00` | Pilot diff `+2.10`, `SE = 0.175` | **>99.9%** |
| **G4c** | `W1p` freq lock & `q` error upper95 `<= 0.030` | Lock Wilson lower `>= 0.65`; `q_err` upper95 `<= 0.030` | Pilot lock `16/16`, `q_err ~ 0.014` | **>98.0%** |
| **G5a** | `W2` (`r=96`, `2^23`) single-run TPR | Wilson lower `>= 0.50` (needs `>= 22/32`) | Conservative `p = 0.80` (`95.9%`); pilot `4/4` | **95.9%** |
| **G5b** | `W2` mean `Z_max` separation over `N23` | Welch 95% CI lower `>= +0.90` | Pilot diff `+2.05`, `SE = 0.172` | **>99.9%** |
| **G5c** | `W2` freq lock & `q` error | Lock Wilson lower `>= 0.60` (needs `>= 25/32`); `q_err` upper95 `<= 0.025` | Pilot lock `4/4` (`p=0.90`: `98.8%`); `q_err=0.006` | **98.8%** |
| **G6a** | `W3` (`r=128`, `2^23`) ensemble `Z_max` sep vs `N23` | Welch 95% CI lower `>= +0.40` | Pilot diff `+0.98`, `SE = 0.13` (`4.25 sigma`) | **98.8%** |
| **G6b** | `W3` frontier freq lock `P(|f_peak - 128| <= 8)` | Wilson lower `>= 0.45` (needs `>= 14/20`; vs chance `0.085`) | Pilot `4/4`; at `p = 0.85` | **97.8%** |
| **G7a** | Crossover Weyl `r=4` recovery ratio | Student-`t` 95% CI lower `>= 0.65` | Step 3 finite `0.806` | **>99.9%** |
| **G7b** | Crossover Weyl `r=32` collapse ratio | Student-`t` 95% CI upper `<= 0.10` | Step 3 finite `0.011` | **>99.9%** |
| **G7c** | Crossover Subspacing `X_sub_q1_r32` vs `X_null_q1` | TPR Wilson lower `>= 0.50`, freq lock lower `>= 0.50`, Welch diff lower `>= +2.00` | Pilot `mean ~ 6.08` vs `2.04` (`SE=0.31`) | **>92.1%** |

**Overall Campaign Verdict:** `CONFIRMED` iff **ALL** gates `G1..G7c` pass simultaneously on their frozen two-sided 95% CIs (`MIXED` / `FAIL` reported verbatim if any gate fails; no post hoc threshold relaxation, seed addition, or re-scoring permitted).

## 5. Provenance and Deviation Policy

1. **P -> F -> R Commit Chain (`AGENTS.md` Merge/Cleanliness Contract):**
   - Commit `P`: Freeze this document, `experiments/subspacing_confirmatory.py`, `experiments/check_subspacing_confirmatory.py`, and `tests/test_subspacing_confirmatory.py` on a clean worktree.
   - Run `--mode preflight` on clean checkout `P` -> records `preflight_code_version = P` and source SHA-256 digests into `docs/prereg/frozen/subspacing_confirmatory/freeze_manifest.json`.
   - Commit `F`: Commit ONLY the freeze manifest (`P -> F` adds no code modifications).
   - Run `--mode campaign` on clean checkout `F` (`HEAD == F`, empty `git status --porcelain`, all frozen source SHA-256 digests verified against `P`).
   - Run clean-room independent verifier `experiments/check_subspacing_confirmatory.py` on saved `results.json` and commit outputs as `R` (`P < F < R` verified by `git merge-base --is-ancestor`).
2. **Interrupt / Crash Recovery:** An atomic temporary checkpoint in gitignored `outputs/subspacing_confirmatory_ckpt.json` resumes only unexecuted deterministic `(arm, seed)` entries without changing seed sequence or completed data; no extra seeds may ever be added on failure.
