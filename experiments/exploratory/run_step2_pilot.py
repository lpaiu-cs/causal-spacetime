"""Exploratory GPU pilot for null tail, q/r frontier, phase, and Poisson checks.

Strictly exploratory seeds (750000+); disjoint from any preregistered confirmatory run.
Outputs incremental JSON checkpoints to outputs/exploratory_step2_pilot.json.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exploratory"))

from subspacing_observer import observe  # noqa: E402
from subspacing_wave_sim import OrderOracle, WaveGeometry  # noqa: E402

FREQUENCIES = np.arange(24.0, 257.0, 4.0)


def json_safe(obj):
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [json_safe(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer, np.bool_)):
        return obj.item()
    return obj


def build_jobs():
    jobs = []
    # 1. Null distribution at N=2^22 and N=2^23 (q=4.0) + q=0.5, 1.0 at N=2^22
    for i in range(16):
        jobs.append(
            dict(
                group="null_22_q4",
                n=2**22,
                q=4.0,
                mode="background",
                r=64,
                phi=0.37,
                poisson=False,
                seed=750000 + i,
            )
        )
    for i in range(16):
        jobs.append(
            dict(
                group="null_23_q4",
                n=2**23,
                q=4.0,
                mode="background",
                r=64,
                phi=0.37,
                poisson=False,
                seed=750100 + i,
            )
        )
    for q in (0.5, 1.0):
        base = 750200 if q == 0.5 else 750210
        for i in range(4):
            jobs.append(
                dict(
                    group=f"null_22_q{q:g}",
                    n=2**22,
                    q=q,
                    mode="background",
                    r=64,
                    phi=0.37,
                    poisson=False,
                    seed=base + i,
                )
            )

    # 2. Frontier & power sweep at N=2^22 across q in (0.5, 1.0, 4.0)
    # and r in (32, 64, 96, 128) + r=192 for q=4.0
    seed_cursor = 751000
    for q in (0.5, 1.0, 4.0):
        r_list = (32, 64, 96, 128, 192) if q == 4.0 else (32, 64, 96, 128)
        for r in r_list:
            for _ in range(4):
                jobs.append(
                    dict(
                        group=f"wave_22_q{q:g}_r{r}",
                        n=2**22,
                        q=q,
                        mode="circular_wave",
                        r=r,
                        phi=0.37,
                        poisson=False,
                        seed=seed_cursor,
                    )
                )
                seed_cursor += 1

    # 3. High-N sub-discreteness frontier at N=2^23
    seed_cursor = 752000
    for q, r in ((4.0, 96), (4.0, 128), (4.0, 192), (1.0, 128)):
        for _ in range(4):
            jobs.append(
                dict(
                    group=f"wave_23_q{q:g}_r{r}",
                    n=2**23,
                    q=q,
                    mode="circular_wave",
                    r=r,
                    phi=0.37,
                    poisson=False,
                    seed=seed_cursor,
                )
            )
            seed_cursor += 1

    # 4. Phase invariance (q=4.0, r=64, N=2^22)
    seed_cursor = 753000
    for phi in (0.0, 1.15, 2.11, 3.85):
        for _ in range(2):
            jobs.append(
                dict(
                    group=f"phase_q4_r64_phi{phi:.2f}",
                    n=2**22,
                    q=4.0,
                    mode="circular_wave",
                    r=64,
                    phi=phi,
                    poisson=False,
                    seed=seed_cursor,
                )
            )
            seed_cursor += 1

    # 5. Poisson-N realism check (N_mean=2^22, q=4.0)
    for i in range(4):
        jobs.append(
            dict(
                group="poisson_22_q4_bg",
                n=2**22,
                q=4.0,
                mode="background",
                r=64,
                phi=0.37,
                poisson=True,
                seed=754000 + i,
            )
        )
    for i in range(4):
        jobs.append(
            dict(
                group="poisson_22_q4_r64",
                n=2**22,
                q=4.0,
                mode="circular_wave",
                r=64,
                phi=0.37,
                poisson=True,
                seed=754100 + i,
            )
        )

    return jobs


def run_one(job, max_probes=65536):
    import cupy as cp

    started = time.perf_counter()
    rng = np.random.default_rng(job["seed"])
    n_draw = int(rng.poisson(job["n"])) if job["poisson"] else int(job["n"])
    geometry = WaveGeometry(job["q"], job["mode"], 1.0 / job["r"], phi=job["phi"])
    points = geometry.sprinkle(n_draw, rng, "gpu")
    oracle = OrderOracle(geometry, points, "gpu")
    gen_s = time.perf_counter() - started

    pilot_n = min(16384, n_draw // 4)
    probe_n = n_draw // 2 - pilot_n
    pilot = np.arange(pilot_n)
    probes = np.arange(pilot_n, pilot_n + probe_n)
    references = np.arange(n_draw // 2, n_draw)
    cache = {}

    def counts(refs, at, early=(), late=(), _oracle=oracle):
        key = tuple(
            hashlib.blake2b(
                np.asarray(g, dtype=np.int64).tobytes(), digest_size=16
            ).digest()
            for g in (refs, at, early, late)
        )
        if key not in cache:
            cache[key] = _oracle.counts(refs, at, early, late)
        return cache[key]

    row = {
        **job,
        "n_actual": n_draw,
        "alpha": float(np.log(n_draw) / np.log(job["r"])),
        "generation_seconds": gen_s,
    }
    try:
        res = observe(
            oracle.query,
            pilot,
            references,
            probes,
            FREQUENCIES,
            q=None,
            count_query=counts,
            max_probes=max_probes,
        )
        idx = int(np.argmin(np.abs(FREQUENCIES - job["r"])))
        whitened = res["whitened_amplitude"]
        off = np.abs(FREQUENCIES - job["r"]) > 20.0
        w_off = (
            float(np.sqrt(np.mean(whitened[off] ** 2))) if off.any() else float("nan")
        )
        row.update(
            status="completed",
            q_hat=float(res["q_hat"]),
            q_err=float(res["q_hat"] - job["q"]),
            calibration_rms=float(res["calibration_rms"]),
            calibration_rms_scaled=float(
                res["calibration_rms"] * np.sqrt(res["n_ref"])
            ),
            calibration_condition=float(res["calibration_condition"]),
            dtheta_norm=float(res["dtheta_norm"]),
            whitened_peak_z=float(res["whitened_peak_z"]),
            whitened_snr=float(res["whitened_snr"]),
            whitened_peak_frequency=float(res["whitened_peak_frequency"]),
            whitened_peak_signed_frequency=float(res["whitened_peak_signed_frequency"]),
            target_z=float(whitened[idx]),
            target_snr=float(whitened[idx] / max(w_off, 1e-30)),
            shot_sd=float(res["shot_sd"]),
            c_ref_estimate=float(res["c_ref_estimate"]),
            block_relative_difference=float(res["block_relative_difference"]),
            tangent_leakage=float(res["tangent_leakage"]),
            accepted=int(res["accepted"]),
            n_ref=int(res["n_ref"]),
            n_probe=int(res["n_probe"]),
            comparable=int(res["comparable"]),
            total_seconds=time.perf_counter() - started,
        )
    except Exception as exc:
        row.update(
            status="failed",
            error=repr(exc),
            total_seconds=time.perf_counter() - started,
        )
    finally:
        del oracle, points
        cp.get_default_memory_pool().free_all_blocks()
    return row


def summarize(rows):
    by_group = {}
    for r in rows:
        by_group.setdefault(r["group"], []).append(r)
    summary = {}
    for grp, items in by_group.items():
        ok = [x for x in items if x.get("status") == "completed"]
        if not ok:
            summary[grp] = {"n": len(items), "completed": 0}
            continue
        pz = np.array([x["whitened_peak_z"] for x in ok])
        snr = np.array([x["whitened_snr"] for x in ok])
        rms = np.array([x["calibration_rms"] for x in ok])
        qerr = np.array([x["q_err"] for x in ok])
        pf = np.array([x["whitened_peak_frequency"] for x in ok])
        hit8 = (
            float(np.mean(np.abs(pf - ok[0]["r"]) <= 8.0))
            if ok[0]["mode"] == "circular_wave"
            else None
        )
        summary[grp] = {
            "n": len(items),
            "completed": len(ok),
            "alpha": ok[0]["alpha"],
            "pz_min": float(pz.min()),
            "pz_med": float(np.median(pz)),
            "pz_max": float(pz.max()),
            "snr_min": float(snr.min()),
            "snr_med": float(np.median(snr)),
            "snr_max": float(snr.max()),
            "rms_max": float(rms.max()),
            "q_abs_err_max": float(np.max(np.abs(qerr))),
            "freq_hit8_rate": hit8,
            "c_ref_zero_count": int(sum(x["c_ref_estimate"] <= 1e-30 for x in ok)),
        }
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path, default=ROOT / "outputs" / "exploratory_step2_pilot.json"
    )
    parser.add_argument("--max-probes", type=int, default=65536)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    existing = {}
    if args.output.exists():
        data = json.loads(args.output.read_text(encoding="utf-8"))
        for r in data.get("runs", []):
            if r.get("status") == "completed":
                existing[(r["group"], r["seed"])] = r

    jobs = build_jobs()
    runs = []
    for idx, job in enumerate(jobs, 1):
        key = (job["group"], job["seed"])
        if key in existing:
            row = existing[key]
        else:
            row = run_one(job, max_probes=args.max_probes)
        runs.append(row)
        print(
            f"[{idx}/{len(jobs)}] {job['group']} seed={job['seed']} "
            f"st={row['status']} pz={row.get('whitened_peak_z', float('nan')):.2f} "
            f"snr={row.get('whitened_snr', float('nan')):.2f} "
            f"pf={row.get('whitened_peak_frequency', float('nan')):.0f} "
            f"rms={row.get('calibration_rms', float('nan')):.2e} "
            f"qh={row.get('q_hat', float('nan')):.3f} "
            f"dt={row.get('total_seconds', 0.0):.1f}s",
            flush=True,
        )
        payload = {
            "max_probes": args.max_probes,
            "runs": runs,
            "summary": summarize(runs),
        }
        tmp = args.output.with_suffix(".tmp")
        tmp.write_text(json.dumps(json_safe(payload), indent=2), encoding="utf-8")
        tmp.replace(args.output)


if __name__ == "__main__":
    main()
