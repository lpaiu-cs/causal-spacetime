"""Audit metadata on production protocol-invariant v3 manifests."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_v3_protocol_diagnostics import (
    manifest_metadata_audit_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="Audit v3 protocol metadata.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def main() -> None:
    config = parse_args()
    rows = manifest_metadata_audit_rows(config.output_dir / "manifests_v3")
    path = write_csv(
        rows,
        data_path(config.output_dir, "v3_protocol_manifest_metadata_audit.csv"),
        ["manifest_id", "family_name", "family_kind"],
    )
    print(f"Wrote v3 protocol manifest metadata audit: {path}")


if __name__ == "__main__":
    main()
