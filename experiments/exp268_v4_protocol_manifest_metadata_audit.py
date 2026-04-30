"""Audit metadata on production v4 protocol manifests."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_v4_diagnostics import (
    manifest_v4_metadata_audit_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="Audit v4 manifest metadata.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def main() -> None:
    config = parse_args()
    rows = manifest_v4_metadata_audit_rows(config.output_dir / "manifests_v4")
    path = write_csv(
        rows,
        data_path(config.output_dir, "v4_protocol_manifest_metadata_audit.csv"),
        ["manifest_id", "family_name", "family_kind"],
    )
    print(f"Wrote v4 protocol manifest metadata audit: {path}")


if __name__ == "__main__":
    main()
