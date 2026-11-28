"""Verify the frozen input hashes and summarize the first deterministic challenge."""

import argparse
import hashlib
import json
import subprocess
from collections import Counter
from pathlib import Path
from statistics import NormalDist

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def summarize(folder):
    manifest = json.loads((folder / "manifest.json").read_text())
    raw = (folder / "rows.jsonl").read_bytes()
    assert manifest["status"] == "diagnostic_completed"
    assert hashlib.sha256(raw).hexdigest() == manifest["rows_sha256"]
    stamp = manifest["code_version"]
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", stamp, "HEAD"], cwd=ROOT, check=True
    )
    for key, path in (
        ("design_sha256", "docs/paper/paper_a/order_wave_energy/design.md"),
        ("source_sha256", "experiments/positive_control/order_wave_energy.py"),
    ):
        content = subprocess.check_output(["git", "show", f"{stamp}:{path}"], cwd=ROOT)
        # Git and the working tree both use LF in these new source files.
        assert hashlib.sha256(content).hexdigest() == manifest[key]
    rows = [json.loads(line) for line in raw.splitlines()]
    counts = Counter(row["stage"] for row in rows)
    assert counts == {
        "limit": 45,
        "transfer": 36,
        "refinement": 3,
        "analytic_transfer": 160,
    }
    limits = [row for row in rows if row["stage"] == "limit"]
    transfers = [row for row in rows if row["stage"] == "transfer"]
    refined = [row for row in rows if row["stage"] == "refinement"]
    theory = [row for row in rows if row["stage"] == "analytic_transfer"]
    ratios = []
    for row in transfers:
        if row["q"] == 0.00125:
            high = next(
                x
                for x in transfers
                if x["q"] == 0.0025
                and x["r"] == row["r"]
                and x["profile"] == row["profile"]
            )
            ratios.append(row["residual"] / high["residual"])
    refinement_differences = []
    for row in refined:
        if row["q"] == 4:
            original = next(
                x
                for x in limits
                if x["q"] == row["q"]
                and x["profile"] == row["profile"]
                and x["lambda"] == 1 / row["r"]
            )
            difference = row["p2_values"][0] - original["p2"]
            unit = "relation_probability"
        else:
            original = next(
                x
                for x in transfers
                if x["q"] == row["q"]
                and x["profile"] == row["profile"]
                and x["r"] == row["r"]
            )
            difference = row["mean_readout_ratio"] - original["mean"]
            unit = "energy_readout_ratio"
        refinement_differences.append(
            {
                "profile": row["profile"],
                "q": row["q"],
                "r": row["r"],
                "difference": difference,
                "unit": unit,
            }
        )
    normal = NormalDist()
    power_factor = normal.inv_cdf(0.975) + normal.inv_cdf(0.95)
    summary = {
        "source_commit": stamp,
        "input_rows_sha256": manifest["rows_sha256"],
        "rows": counts,
        "elapsed_seconds": manifest["elapsed_seconds"],
        "max_limit_error_at_lambda_1_over_128": max(
            abs(x["residual"]) for x in limits if x["lambda"] == 1 / 128
        ),
        "max_limit_error_relative_to_curvature_contrast": max(
            abs(x["residual"] / (0.1 - x["background_p2"]))
            for x in limits
            if x["lambda"] == 1 / 128
        ),
        "max_relative_transfer_error": {
            str(q): max(
                abs(x["residual"] / x["prediction"]) for x in transfers if x["q"] == q
            )
            for q in (0.0025, 0.00125)
        },
        "residual_ratio_on_halving_q": [min(ratios), max(ratios)],
        "prediction_grid_maximum": max(theory, key=lambda x: x["sine"]),
        "refinement_differences": refinement_differences,
        "minimum_rosen_a": min(x["minimum_a"] for x in limits),
        "flat_limit_variance": "Var(P2hat)=(0.04*N+0.10)/(N*(N-1))",
        "weak_field_normal_planning_only": {
            "coverage": 0.95,
            "sidedness": "two-sided",
            "target_power": 0.90,
            "relative_accuracy": 0.25,
            "q_L_squared": 0.00125,
            "N_approx": (power_factor * 52.5 / (0.25 * 0.00125)) ** 2,
            "warning": "Raw P2 only; not a campaign or a bound for all estimators",
        },
    }
    (folder / "summary.json").write_text(json.dumps(summary, indent=2))

    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout="constrained")
    for profile, color, marker in (("sine", "#2463ad", "o"), ("mixed", "#b85320", "s")):
        axes[0].semilogx(
            [x["r"] for x in theory],
            [x[profile] for x in theory],
            color=color,
            label=f"{profile}: derived response",
        )
        data = [x for x in transfers if x["profile"] == profile and x["q"] == 0.00125]
        axes[0].scatter(
            [x["r"] for x in data],
            [x["mean"] for x in data],
            color=color,
            marker=marker,
            s=30,
            facecolors="none",
            label=f"{profile}: exact vacuum, phase mean",
        )
    axes[0].axhline(1, color="#777777", linestyle="--", linewidth=1)
    axes[0].set(
        xlabel=r"window size / inverse frequency, $r=\omega L$",
        ylabel=r"order readout / effective Ricci, $\langle E\rangle/(2q)$",
        title=r"Resolution response ($qL^2=0.00125$)",
        ylim=(-0.02, 1.47),
        xlim=(0.1, 128),
    )
    axes[0].legend(fontsize=8, loc="upper left")
    for profile, color in (
        ("sine", "#2463ad"),
        ("mixed", "#b85320"),
        ("circular", "#437c50"),
    ):
        data = [x for x in limits if x["profile"] == profile and x["q"] == 1]
        axes[1].loglog(
            [1 / x["lambda"] for x in data],
            [abs(x["residual"] / (0.1 - x["background_p2"])) for x in data],
            marker="o",
            color=color,
            label=profile,
        )
    axes[1].set(
        xlabel=r"frequency parameter, $L/\lambda$",
        ylabel="error / background curvature contrast",
        title=r"Approach to the effective geometry ($qL^2=1$)",
    )
    axes[1].legend(fontsize=9)
    for ax in axes:
        ax.grid(alpha=0.15, which="both")
    fig.savefig(folder / "response.png", dpi=180)
    fig.savefig(folder / "response.pdf")
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    summarize(parser.parse_args().folder.resolve())
