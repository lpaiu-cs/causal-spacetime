"""Power and cost for the finite-N test of the five-window correction (design aid).

Deterministic except for seeded simulations; no new event is drawn. Inputs:

- the multi-window study's window fits (continuum, the benchmark run's clock
  covariances, cases 0 and 1; outputs/intrinsic_weyl_window_study_v1) for the
  continuum readings and for the variance of the corrected reading relative to
  one window;
- the 63 completed fifth-degree covariances of the flat calibration at
  N = 2^24 (outputs/intrinsic_weyl_calibration_v1) for the spread of the
  one-window reading over selections and orientations.

The reading of an estimate c_hat against the truth c is the tidal projection
<c, c_hat> / <c, c> (full tensor, unweighted). For R replicates the gate is
decided on the 95% two-sided t-interval of the mean reading; per-replicate
readings are simulated as the continuum reading plus Gaussian noise whose
standard deviation is drawn from the selection-by-orientation distribution,
scaled by sqrt(kappa) for a possible variance understatement.
"""

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SPECS = {
    "single": ("single_window", 0.0, "none"),
    "three_point_plane_wave": ("three_point_first", 0.5, "plane_wave_estimated_k"),
    "five_point_plane_wave": ("five_point_second", 0.5, "plane_wave_estimated_k"),
    "five_point_minimum_norm": ("five_point_second", 0.5, "minimum_norm_time_time"),
}
TOLERANCES = (0.03, 0.04, 0.05, 0.06)
REPLICATES = (
    8,
    12,
    16,
    20,
    24,
    28,
    32,
    40,
    48,
    56,
    64,
    80,
    96,
    112,
    128,
    160,
    192,
    256,
)
LEVELS = (2**24, 2**25, 2**26)
TARGET_POWER = 0.9
# Cost model per window, scaled from one flat calibration replicate at 2^24
# (131 landmarks, about 3.8e3 clocks, two BLAS threads): event generation and
# relation bits 55.8 s and clock moments 15.7 s on the GPU, fifth-degree fit
# 64.7 s on the CPU, peak 14.8 GB. Side windows carry about 1.25 times the
# landmarks and 1.5-1.6 times the clocks (frame-link study), assumed to cost
# 1.25x generation, 1.6x moments and 3x fit time. To be checked by the pilot.
GENERATION_SECONDS = 55.8
MOMENT_SECONDS = 15.7
FIT_SECONDS = 64.7
PEAK_GB_AT_2_24 = 14.8


def replicate_cost(level):
    """GPU seconds, CPU fit seconds and peak memory of one five-window replicate."""
    scale = level / 2**24
    gpu = (1 + 4 * 1.25) * GENERATION_SECONDS * scale + (
        1 + 4 * 1.6
    ) * MOMENT_SECONDS * scale
    cpu = (1 + 4 * 3.0) * FIT_SECONDS
    return {
        "gpu_seconds": gpu,
        "cpu_fit_seconds": cpu,
        "peak_gb": PEAK_GB_AT_2_24 * 1.25 * scale,
    }


def reading(estimate, covariance, c):
    """Tidal projection of the estimate on c and its variance."""
    from intrinsic_weyl_joint import TIDAL_GRAM

    norm2 = c @ TIDAL_GRAM @ c
    value = (c @ TIDAL_GRAM @ estimate) / norm2
    variance = (c @ TIDAL_GRAM @ covariance @ TIDAL_GRAM @ c) / norm2**2
    return float(value), float(variance)


def continuum_inputs(window_study, amplitude, orientations):
    """Continuum readings and corrected-to-single variance ratios, cases 0, 1."""
    import intrinsic_weyl_window_study as ws

    windows = json.loads((window_study / "summary_windows.json").read_text())
    truth = amplitude * ws.unit_truth()
    result = {}
    for case in (0, 1):
        chosen = [w for w in windows if w["case"] == case]
        _, context = ws.analyse_case(case, chosen, window_study / "cache", amplitude)
        window_set = context["window_sets"][2.0]
        specs = ws.correction_specs(
            sorted(window_set),
            context["alias_n"],
            np.array([1.0, 0, 0, 1]),
            context["estimated_k"][2.0],
            context["minimum_norm"],
        )
        combined = {}
        for label, key in SPECS.items():
            _, _, _, shifts, weights = next(s for s in specs if s[:3] == key)
            combined[label] = ws.combine([window_set[s] for s in shifts], weights)
        row = {}
        for label, (estimate, covariance) in combined.items():
            value, variance = reading(estimate, covariance, truth)
            ratios = []
            for rotation in orientations:
                c = ws.weyl_rotation(rotation) @ truth
                ratios.append(
                    reading(estimate, covariance, c)[1]
                    / reading(*combined["single"], c)[1]
                )
            row[label] = {
                "continuum_reading": value,
                "reading_sd_at_2_24_realized": float(np.sqrt(variance)),
                "variance_ratio_to_single_median": float(np.median(ratios)),
                "variance_ratio_to_single_5_95": np.quantile(
                    ratios, [0.05, 0.95]
                ).tolist(),
            }
        result[case] = row
    return result


def single_window_spread(calibration, amplitude, orientations):
    """One-window reading variances over flat calibration selections x orientations."""
    import intrinsic_weyl_window_study as ws

    truth = amplitude * ws.unit_truth()
    variances = []
    for path in sorted(glob.glob(str(calibration / "case_*_n16777216.json"))):
        record = json.loads(Path(path).read_text())
        if record.get("evaluation", {}).get("kind") != "flat":
            continue
        pure = json.loads(Path(path.replace(".json", "_pure.json")).read_text())
        covariance = np.asarray(pure["fits"]["order_5"]["nominal_component_covariance"])
        for rotation in orientations:
            variances.append(
                reading(np.zeros(10), covariance, ws.weyl_rotation(rotation) @ truth)[1]
            )
    return np.array(variances)


def gate_power(mean, sds, replicates, tolerance, rng, draws=20000):
    """P(95% t-interval of the mean reading lies inside [1 - t, 1 + t])."""
    from scipy.stats import t as student

    quantile = student.ppf(0.975, replicates - 1)
    sigma = rng.choice(sds, size=(draws, replicates))
    values = mean + sigma * rng.normal(size=(draws, replicates))
    centre = values.mean(axis=1)
    half = quantile * values.std(axis=1, ddof=1) / np.sqrt(replicates)
    return float(
        np.mean((centre - half > 1 - tolerance) & (centre + half < 1 + tolerance))
    )


def below_power(mean, sds, replicates, bound, rng, draws=20000):
    """P(95% t-interval of the mean reading lies below ``bound``)."""
    from scipy.stats import t as student

    quantile = student.ppf(0.975, replicates - 1)
    sigma = rng.choice(sds, size=(draws, replicates))
    values = mean + sigma * rng.normal(size=(draws, replicates))
    centre = values.mean(axis=1)
    half = quantile * values.std(axis=1, ddof=1) / np.sqrt(replicates)
    return float(np.mean(centre + half < bound))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--window-study",
        type=Path,
        default=ROOT / "outputs/intrinsic_weyl_window_study_v1",
    )
    parser.add_argument(
        "--calibration",
        type=Path,
        default=ROOT / "outputs/intrinsic_weyl_calibration_v1",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from intrinsic_weyl_bombelli_benchmark import benchmark_amplitude
    from intrinsic_weyl_landmark_design import code_state
    from scipy.spatial.transform import Rotation

    amplitude = float(benchmark_amplitude())
    commit, clean = code_state()
    orientations = list(Rotation.random(400, random_state=2026112107).as_matrix())
    continuum = continuum_inputs(args.window_study, amplitude, orientations)
    single = single_window_spread(args.calibration, amplitude, orientations[:100])
    result = {
        "design_aid": "power and cost for the finite-N five-window test; no new events",
        "code_version": commit,
        "worktree_clean": clean,
        "benchmark_amplitude": amplitude,
        "continuum": {str(k): v for k, v in continuum.items()},
        "single_window_reading_sd_at_2_24": {
            "selections_x_orientations": int(len(single)),
            "rms": float(np.sqrt(single.mean())),
            "median": float(np.sqrt(np.median(single))),
            "5_95": np.sqrt(np.quantile(single, [0.05, 0.95])).tolist(),
        },
        "designs": {},
    }
    rng = np.random.default_rng(2026112108)
    for label in (
        "five_point_minimum_norm",
        "five_point_plane_wave",
        "three_point_plane_wave",
    ):
        ratio = max(
            continuum[c][label]["variance_ratio_to_single_median"] for c in (0, 1)
        )
        mean = float(
            np.mean([continuum[c][label]["continuum_reading"] for c in (0, 1)])
        )
        sds = np.sqrt(single * ratio)
        design = {
            "variance_ratio_used": ratio,
            "continuum_reading_mean": mean,
            "rows": [],
        }
        levels = LEVELS if label != "three_point_plane_wave" else LEVELS[:1]
        for level in levels:
            # Linear regime: the reading variance scales as 1/N.
            scaled = sds * np.sqrt(2**24 / level)
            for kappa in (1.0, 1.2):
                for shift in (0.0, -0.01, -0.02):
                    for tolerance in TOLERANCES:
                        powers = {
                            r: gate_power(
                                mean + shift, scaled * np.sqrt(kappa), r, tolerance, rng
                            )
                            for r in REPLICATES
                        }
                        needed = next(
                            (r for r in REPLICATES if powers[r] >= TARGET_POWER), None
                        )
                        row = {
                            "events_per_window": level,
                            "kappa": kappa,
                            "mean_reading": mean + shift,
                            "tolerance": tolerance,
                            "power_by_replicates": powers,
                            "replicates_for_90pct": needed,
                        }
                        if needed is not None:
                            cost = replicate_cost(level)
                            row["cost_at_needed"] = {
                                "measurement_events": 5 * level * needed,
                                "gpu_hours": cost["gpu_seconds"] * needed / 3600,
                                "cpu_fit_hours": cost["cpu_fit_seconds"]
                                * needed
                                / 3600,
                                "peak_gb_per_window": cost["peak_gb"],
                            }
                        design["rows"].append(row)
                        print(
                            label,
                            level,
                            kappa,
                            round(mean + shift, 4),
                            tolerance,
                            needed,
                            flush=True,
                        )
        result["designs"][label] = design
    single_mean = float(
        np.mean([continuum[c]["single"]["continuum_reading"] for c in (0, 1)])
    )
    result["control"] = {
        "single_reading_mean": single_mean,
        "below_1_minus_t_power": {
            str(tolerance): {
                r: below_power(
                    single_mean, np.sqrt(single * 1.2), r, 1 - tolerance, rng
                )
                for r in (16, 32, 64)
            }
            for tolerance in TOLERANCES
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
