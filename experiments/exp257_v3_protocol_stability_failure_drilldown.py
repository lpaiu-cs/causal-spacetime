"""Drill down v3 protocol restart and latent-order instability."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import save_metric_bar
from v3_protocol_blocked_v4_experiment_helpers import (
    data_path,
    read_csv_rows,
    to_float_rows,
)

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_stability_drilldown import (
    summarize_v3_protocol_stability_failures,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 protocol stability drilldown.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    return summarize_v3_protocol_stability_failures(
        to_float_rows(
            read_csv_rows(
                data_path(
                    config.output_dir,
                    "v3_protocol_manifest_family_restart_stability.csv",
                )
            )
        ),
        to_float_rows(
            read_csv_rows(
                data_path(
                    config.output_dir,
                    "v3_protocol_manifest_family_latent_order_stability.csv",
                )
            )
        ),
        default_cross_family_robustness_criteria(),
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_csv(
        rows,
        data_path(config.output_dir, "v3_protocol_stability_failure_drilldown.csv"),
        ["family_name", "restart_std", "latent_order_disagreement"],
    )
    figure = save_metric_bar(
        rows,
        label_key="family_name",
        value_key="latent_order_disagreement",
        path=config.output_dir / "figures" / "v3_protocol_stability_failures.png",
        ylabel="latent-order disagreement",
    )
    print(f"Wrote v3 protocol stability drilldown: {path}")
    print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()

