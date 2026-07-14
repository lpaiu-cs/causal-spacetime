"""P7 N=600 dense-grid characterization with start-separated diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from p3_dynamics import analyze_order
from p5_two_orders_emergence import order_inputs
from pc_common import git_describe, write_rows_csv

from causal_spacetime_lab.positive_control.geometry_score import (
    geometry_order_parameter,
)
from causal_spacetime_lab.positive_control.mcmc_diagnostics import (
    IAT_ESTIMATOR,
    classify_phase,
    integrated_autocorrelation,
)
from causal_spacetime_lab.positive_control.two_orders import bipartite_perm

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "positive_control"
FROZEN = ROOT / "docs" / "prereg" / "frozen"
CONSTANTS = FROZEN / "p7_n600_stage_a_constants.json"


def _config(constants: dict, chain: int) -> dict:
    matches = [row for row in constants["chains_per_beta"] if row["chain"] == chain]
    if len(matches) != 1:
        raise SystemExit("chain must be 0, 1, or 2")
    return matches[0]


def _seed(base: int, beta: float, chain: int) -> int:
    return int(base + 1000 * int(beta) + chain)


def _score(result: dict) -> float:
    return geometry_order_parameter(
        status=result["status"],
        heldout=result.get("heldout"),
        null_gap=result.get("null_gap"),
        truth_error=result.get("truth"),
    )


def run_chain(beta: float, chain: int) -> None:
    from causal_spacetime_lab.positive_control.accelerated_two_orders import (
        mcmc_2d_order_replay_accelerated,
    )

    constants = json.loads(CONSTANTS.read_text(encoding="utf-8"))
    if beta not in constants["new_betas"]:
        raise SystemExit(f"new beta must be one of {constants['new_betas']}")
    config = _config(constants, chain)
    if config["start"] == "bipartite":
        initial = bipartite_perm(constants["n_elements"])
    else:
        initial_seed = _seed(config["initial_seed_base"], beta, chain)
        initial = np.random.default_rng(initial_seed).permutation(
            constants["n_elements"]
        )
    proposal_seed = _seed(constants["proposal_seed_base"], beta, chain)
    samples, acceptance, permutations = mcmc_2d_order_replay_accelerated(
        initial,
        beta=beta,
        eps=constants["eps"],
        steps=constants["steps"],
        seed=proposal_seed,
        sample_every=constants["sample_every"],
        burn_frac=constants["burn_frac"],
        collect_perms=True,
        resync_every=constants["resync_every"],
    )
    if not constants["minimum_sample_count"] <= len(samples) <= constants[
        "nominal_sample_count"
    ]:
        raise RuntimeError(f"unexpected retained sample count: {len(samples)}")

    phase = classify_phase(samples, constants["phase_reference"])
    version = git_describe()
    tag = f"b{beta:g}_c{chain}"
    chain_rows = [
        {
            "beta": beta,
            "chain": chain,
            "start": config["start"],
            "proposal_seed": proposal_seed,
            "sample": index,
            "acceptance": acceptance,
            "phase": phase,
            "code_version": version,
            **sample,
        }
        for index, sample in enumerate(samples)
    ]
    write_rows_csv(OUT / f"p7_n600_chain_{tag}.csv", chain_rows)

    instrument_rows = []
    for sample_index in (0, len(samples) // 2, len(samples) - 1):
        causal, times, coords = order_inputs(permutations[sample_index])
        result = analyze_order(
            causal,
            times,
            coords,
            seed=proposal_seed + sample_index,
            want_truth=True,
        )
        result["G"] = _score(result)
        instrument_rows.append(
            {
                "beta": beta,
                "chain": chain,
                "start": config["start"],
                "proposal_seed": proposal_seed,
                "sample": sample_index,
                "acceptance": acceptance,
                "phase": phase,
                "code_version": version,
                **result,
            }
        )
    write_rows_csv(OUT / f"p7_n600_instrument_{tag}.csv", instrument_rows)
    mean_g = float(np.mean([row["G"] for row in instrument_rows]))
    print(
        f"beta={beta:g} chain={chain} start={config['start']} phase={phase} "
        f"acceptance={acceptance:.4f} G={mean_g:.3f}",
        flush=True,
    )


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _paths(beta: float, chain: int, reused: bool) -> tuple[Path, Path]:
    prefix = "p6_near" if reused else "p7_n600"
    base = FROZEN if reused else OUT
    tag = f"b{beta:g}_c{chain}"
    return base / f"{prefix}_chain_{tag}.csv", base / f"{prefix}_instrument_{tag}.csv"


def _binder(values: list[float]) -> float:
    array = np.asarray(values, dtype=float)
    centered = array - array.mean()
    second = float(np.mean(centered**2))
    if second == 0.0:
        return float("nan")
    return float(1.0 - np.mean(centered**4) / (3.0 * second**2))


def _instrument_score(row: dict[str, str]) -> float:
    if row.get("G") not in (None, ""):
        return float(row["G"])
    return geometry_order_parameter(
        status=row["status"],
        heldout=float(row["heldout"]) if row.get("heldout") else None,
        null_gap=float(row["null_gap"]) if row.get("null_gap") else None,
        truth_error=float(row["truth"]) if row.get("truth") else None,
    )


def _validate_shard_contents(
    *,
    beta: float,
    chain: int,
    start: str,
    chain_rows: list[dict[str, str]],
    instrument_rows: list[dict[str, str]],
    minimum_samples: int,
    nominal_samples: int,
) -> None:
    """Require one complete P7 chain and its three instrument snapshots."""

    identity = (float(beta), int(chain), start)
    for label, rows in (("chain", chain_rows), ("instrument", instrument_rows)):
        try:
            identities = {
                (float(row["beta"]), int(row["chain"]), row["start"])
                for row in rows
            }
        except (KeyError, TypeError, ValueError) as error:
            raise SystemExit(
                f"malformed P7 {label} row for beta={beta:g} chain={chain}: {error}"
            ) from error
        if identities != {identity}:
            raise SystemExit(
                f"P7 {label} identity mismatch for beta={beta:g} chain={chain}: "
                f"observed={sorted(identities)}"
            )

    if not minimum_samples <= len(chain_rows) <= nominal_samples:
        raise SystemExit(
            f"P7 chain beta={beta:g} chain={chain} has {len(chain_rows)} samples; "
            f"expected {minimum_samples}-{nominal_samples}"
        )
    try:
        sample_indices = sorted(int(row["sample"]) for row in chain_rows)
        instrument_indices = sorted(int(row["sample"]) for row in instrument_rows)
    except (KeyError, TypeError, ValueError) as error:
        raise SystemExit(
            f"malformed P7 sample index for beta={beta:g} chain={chain}: {error}"
        ) from error
    if sample_indices != list(range(len(chain_rows))):
        raise SystemExit(
            f"P7 chain beta={beta:g} chain={chain} has missing or duplicate "
            "sample indices"
        )
    expected_instrument = [0, len(chain_rows) // 2, len(chain_rows) - 1]
    if instrument_indices != expected_instrument:
        raise SystemExit(
            f"P7 instrument beta={beta:g} chain={chain} has sample indices "
            f"{instrument_indices}; expected {expected_instrument}"
        )


def _group_chains_by_start(chains: list[dict]) -> dict[str, list[dict]]:
    """Return the required two-random/one-bipartite start groups."""

    grouped = {
        start: sorted(
            (row for row in chains if row["start"] == start),
            key=lambda row: int(row["chain"]),
        )
        for start in ("random", "bipartite")
    }
    unexpected = sorted(
        {row["start"] for row in chains} - {"random", "bipartite"}
    )
    if len(grouped["random"]) != 2 or len(grouped["bipartite"]) != 1 or unexpected:
        raise SystemExit(
            "P7 aggregation requires exactly two random starts and one "
            f"bipartite start; got random={len(grouped['random'])}, "
            f"bipartite={len(grouped['bipartite'])}, unexpected={unexpected}"
        )
    return grouped


def _write_figure(beta_rows: list[dict]) -> None:
    betas = [row["beta"] for row in beta_rows]
    figure, axes = plt.subplots(1, 3, figsize=(11.2, 3.5))
    for start, label, marker in (
        ("random", "random starts", "o"),
        ("bipartite", "bipartite start", "s"),
    ):
        axes[0].plot(
            betas,
            [row[f"mean_S_{start}"] for row in beta_rows],
            marker=marker,
            label=label,
        )
        axes[1].plot(
            betas,
            [row[f"mean_n0_{start}"] for row in beta_rows],
            marker=marker,
        )
        axes[2].plot(
            betas,
            [row[f"mean_G_{start}"] for row in beta_rows],
            marker=marker,
        )
    axes[0].set_ylabel("Mean action")
    axes[1].set_ylabel("Mean n0")
    axes[2].set_ylabel("Mean G")
    axes[2].axhline(0.5, color="black", linestyle="--", linewidth=1)
    for axis in axes:
        axis.set_xlabel("beta")
        axis.grid(alpha=0.2)
    axes[0].legend(frameon=False)
    figure.tight_layout()
    figure.savefig(OUT / "p7_n600_stage_a_observables.png", dpi=180)
    figure.savefig(OUT / "p7_n600_stage_a_observables.pdf")
    plt.close(figure)


def aggregate() -> None:
    constants = json.loads(CONSTANTS.read_text(encoding="utf-8"))
    chain_summaries = []
    beta_summaries = []
    all_chain_rows: dict[tuple[float, int], list[dict[str, str]]] = {}
    all_instrument_rows: dict[tuple[float, int], list[dict[str, str]]] = {}

    for beta in constants["all_betas"]:
        reused = beta in constants["reused_p6_betas"]
        for config in constants["chains_per_beta"]:
            chain = int(config["chain"])
            chain_path, instrument_path = _paths(beta, chain, reused)
            if not chain_path.exists() or not instrument_path.exists():
                raise SystemExit(
                    f"missing P7 chain inputs for beta={beta:g} chain={chain}"
                )
            rows = _read(chain_path)
            instrument = _read(instrument_path)
            _validate_shard_contents(
                beta=beta,
                chain=chain,
                start=config["start"],
                chain_rows=rows,
                instrument_rows=instrument,
                minimum_samples=constants["minimum_sample_count"],
                nominal_samples=constants["nominal_sample_count"],
            )
            all_chain_rows[(beta, chain)] = rows
            all_instrument_rows[(beta, chain)] = instrument
            summary = {
                "beta": beta,
                "chain": chain,
                "start": rows[0]["start"],
                "phase": rows[0]["phase"],
                "acceptance": float(rows[0]["acceptance"]),
                "n_samples": len(rows),
                "mean_G": float(
                    np.mean([_instrument_score(row) for row in instrument])
                ),
                "n_G": len(instrument),
                "source": "P6-reuse" if reused else "P7-new",
            }
            for field in ("S", "n0", "height"):
                values = [float(row[field]) for row in rows]
                tau, ess = integrated_autocorrelation(values)
                summary[f"mean_{field}"] = float(np.mean(values))
                summary[f"iat_{field}"] = tau
                summary[f"ess_{field}"] = ess
            chain_summaries.append(summary)

    for beta in constants["all_betas"]:
        chains = [row for row in chain_summaries if row["beta"] == beta]
        grouped = _group_chains_by_start(chains)
        phases = [row["phase"] for row in chains]
        min_ess = min(
            row[f"ess_{field}"] for row in chains for field in ("S", "n0", "height")
        )
        row = {
            "beta": beta,
            "phase_random_0": grouped["random"][0]["phase"],
            "phase_random_1": grouped["random"][1]["phase"],
            "phase_bipartite": grouped["bipartite"][0]["phase"],
            "phase_agreement": len(set(phases)) == 1,
            "minimum_ess": min_ess,
            "mixing_screen_adequate": len(set(phases)) == 1
            and min_ess >= constants["ess_screening_min"],
        }
        for start, start_chains in grouped.items():
            chain_ids = [int(item["chain"]) for item in start_chains]
            samples = [
                sample
                for chain in chain_ids
                for sample in all_chain_rows[(beta, chain)]
            ]
            instruments = [
                instrument
                for chain in chain_ids
                for instrument in all_instrument_rows[(beta, chain)]
            ]
            n0 = [float(sample["n0"]) for sample in samples]
            row[f"mean_S_{start}"] = float(
                np.mean([float(sample["S"]) for sample in samples])
            )
            row[f"mean_n0_{start}"] = float(np.mean(n0))
            row[f"susceptibility_n0_{start}"] = float(
                np.var(n0, ddof=1) / constants["n_elements"]
            )
            row[f"binder_n0_{start}"] = _binder(n0)
            row[f"mean_height_{start}"] = float(
                np.mean([float(sample["height"]) for sample in samples])
            )
            row[f"mean_G_{start}"] = float(
                np.mean([_instrument_score(item) for item in instruments])
            )
        row["dual_start_delta_S"] = (
            row["mean_S_bipartite"] - row["mean_S_random"]
        )
        row["dual_start_delta_n0"] = (
            row["mean_n0_bipartite"] - row["mean_n0_random"]
        )
        row["dual_start_delta_G"] = (
            row["mean_G_bipartite"] - row["mean_G_random"]
        )
        beta_summaries.append(row)

    write_rows_csv(OUT / "p7_n600_stage_a_chain_summary.csv", chain_summaries)
    write_rows_csv(OUT / "p7_n600_stage_a_beta_summary.csv", beta_summaries)
    result = {
        "stage": "P7-N600-A",
        "code_version": git_describe(),
        "characterization_only": True,
        "iat_estimator": IAT_ESTIMATOR,
        "chain_summaries": chain_summaries,
        "beta_summaries": beta_summaries,
        "mixing_adequate_betas": [
            row["beta"] for row in beta_summaries if row["mixing_screen_adequate"]
        ],
        "inference_warning": (
            "Beta points failing the frozen mixing screen are start-separated "
            "finite-run characterizations, not equilibrium estimates."
        ),
    }
    (OUT / "p7_n600_stage_a_summary.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    _write_figure(beta_summaries)
    print(json.dumps(result, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("chain", "aggregate"), required=True)
    parser.add_argument("--beta", type=float)
    parser.add_argument("--chain", type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage == "chain":
        if args.beta is None or args.chain is None:
            raise SystemExit("--stage chain requires --beta and --chain")
        run_chain(args.beta, args.chain)
    else:
        aggregate()


if __name__ == "__main__":
    main()
