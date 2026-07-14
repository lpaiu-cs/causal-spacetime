"""P6a near-critical N=600 characterization and P7 foundation."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from p3_dynamics import analyze_order
from p5_two_orders_emergence import order_inputs
from pc_common import git_describe, write_rows_csv

from causal_spacetime_lab.positive_control.accelerated_two_orders import (
    mcmc_2d_order_replay_accelerated,
)
from causal_spacetime_lab.positive_control.mcmc_diagnostics import (
    classify_phase,
    integrated_autocorrelation,
)
from causal_spacetime_lab.positive_control.two_orders import bipartite_perm

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "positive_control"
FROZEN = ROOT / "docs" / "prereg" / "frozen" / "p6_near_critical_constants.json"


def _gate_pass(row: dict, gates: dict) -> bool:
    return bool(
        row["status"] == "ok"
        and row["heldout"] <= gates["heldout_max"]
        and row["null_gap"] >= gates["null_gap_min"]
        and row["truth"] <= gates["truth_max"]
    )


def _chain_config(constants: dict, beta: float, chain: int) -> dict:
    if beta not in constants["betas"]:
        raise SystemExit(f"beta must be one of {constants['betas']}")
    configs = [
        entry for entry in constants["chains_per_beta"] if entry["chain"] == chain
    ]
    if len(configs) != 1:
        raise SystemExit("chain must be 0, 1, or 2")
    return configs[0]


def run_chain(beta: float, chain: int) -> None:
    constants = json.loads(FROZEN.read_text(encoding="utf-8"))
    config = _chain_config(constants, beta, chain)
    beta_index = constants["betas"].index(beta)
    if config["start"] == "bipartite":
        pi0 = bipartite_perm(constants["n_elements"])
    else:
        initial_seed = config["initial_seed_offset"] + 10 * beta_index + chain
        pi0 = np.random.default_rng(initial_seed).permutation(constants["n_elements"])
    proposal_seed = constants["proposal_seed_base"] + 100 * beta_index + chain
    samples, acceptance, permutations = mcmc_2d_order_replay_accelerated(
        pi0,
        beta=beta,
        eps=constants["eps"],
        steps=constants["steps"],
        seed=proposal_seed,
        sample_every=constants["sample_every"],
        burn_frac=constants["burn_frac"],
        collect_perms=True,
        resync_every=constants["resync_every"],
    )
    if (
        not constants["minimum_sample_count"]
        <= len(samples)
        <= constants["nominal_sample_count"]
    ):
        raise RuntimeError(
            f"expected {constants['minimum_sample_count']}-"
            f"{constants['nominal_sample_count']} samples, got {len(samples)}"
        )

    version = git_describe()
    phase = classify_phase(samples, constants["phase_reference"])
    rows = []
    for sample_index, row in enumerate(samples):
        rows.append(
            {
                "beta": beta,
                "chain": chain,
                "start": config["start"],
                "proposal_seed": proposal_seed,
                "sample": sample_index,
                "acceptance": acceptance,
                "phase": phase,
                "code_version": version,
                **row,
            }
        )
    tag = f"b{beta:g}_c{chain}"
    write_rows_csv(OUT / f"p6_near_chain_{tag}.csv", rows)

    instrument_rows = []
    instrument_indices = [0, len(samples) // 2, len(samples) - 1]
    for sample_index in instrument_indices:
        causal, times, coords = order_inputs(permutations[sample_index])
        result = analyze_order(
            causal,
            times,
            coords,
            seed=proposal_seed + sample_index,
            want_truth=True,
        )
        result["gate_pass"] = float(_gate_pass(result, constants["instrument_gates"]))
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
    write_rows_csv(OUT / f"p6_near_instrument_{tag}.csv", instrument_rows)
    print(
        f"beta={beta:g} chain={chain} start={config['start']} phase={phase} "
        f"acceptance={acceptance:.4f}",
        flush=True,
    )


def _read_rows(pattern: str) -> list[dict[str, str]]:
    rows = []
    for path in sorted(OUT.glob(pattern)):
        with path.open(newline="", encoding="utf-8") as handle:
            rows.extend(csv.DictReader(handle))
    return rows


def aggregate() -> None:
    constants = json.loads(FROZEN.read_text(encoding="utf-8"))
    expected = {
        (float(beta), int(config["chain"]))
        for beta in constants["betas"]
        for config in constants["chains_per_beta"]
    }
    chain_rows = _read_rows("p6_near_chain_b*_c*.csv")
    instrument_rows = _read_rows("p6_near_instrument_b*_c*.csv")
    observed = {(float(row["beta"]), int(row["chain"])) for row in chain_rows}
    if observed != expected:
        missing = sorted(expected - observed)
        raise SystemExit(f"near-critical aggregation missing chains: {missing}")

    chains = []
    for beta, chain in sorted(expected):
        subset = [
            row
            for row in chain_rows
            if float(row["beta"]) == beta and int(row["chain"]) == chain
        ]
        instrument = [
            row
            for row in instrument_rows
            if float(row["beta"]) == beta and int(row["chain"]) == chain
        ]
        summary = {
            "beta": beta,
            "chain": chain,
            "start": subset[0]["start"],
            "phase": subset[0]["phase"],
            "acceptance": float(subset[0]["acceptance"]),
            "n_samples": len(subset),
            "instrument_valid": sum(row["status"] == "ok" for row in instrument),
            "instrument_gate_pass": sum(float(row["gate_pass"]) for row in instrument),
        }
        for field in ("S", "n0", "height"):
            values = [float(row[field]) for row in subset]
            tau, ess = integrated_autocorrelation(values)
            summary[f"mean_{field}"] = float(np.mean(values))
            summary[f"iat_{field}"] = tau
            summary[f"ess_{field}"] = ess
        chains.append(summary)
    write_rows_csv(OUT / "p6_near_critical_chain_summary.csv", chains)
    result = {
        "code_version": git_describe(),
        "characterization_only": True,
        "chains": chains,
        "dual_start_phase_by_beta": {
            f"{beta:g}": [row["phase"] for row in chains if row["beta"] == beta]
            for beta in constants["betas"]
        },
    }
    (OUT / "p6_near_critical_summary.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
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
