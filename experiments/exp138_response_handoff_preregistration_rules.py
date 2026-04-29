"""Write response handoff preregistration rules."""

from __future__ import annotations

import csv
from pathlib import Path

from causal_spacetime_lab.state_change_response_preregistration import (
    preregistration_rule_table,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def run_experiment() -> list[dict[str, str]]:
    """Return preregistration rule rows."""

    return preregistration_rule_table()


def write_outputs(
    rows: list[dict[str, str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write preregistration rule table."""

    path = output_dir / "data" / "response_handoff_preregistration_rules.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote response handoff preregistration rules: {output_path}")


if __name__ == "__main__":
    main()
