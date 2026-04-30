from __future__ import annotations

import csv
import json
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_remediation_plan import (
    build_preregistered_remediation_plan,
    default_new_manifest_family_specs_v2,
    remediation_plan_to_jsonable,
)
from causal_spacetime_lab.state_change_manifest_v2_spec import (
    load_v2_family_specs_from_design_csv,
    load_v2_family_specs_from_remediation_plan,
    mark_v2_family_specs_executed,
)


def _write_design_csv(path: Path) -> None:
    rows = default_new_manifest_family_specs_v2()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_load_v2_family_specs_from_design_csv(tmp_path: Path) -> None:
    path = tmp_path / "new_manifest_family_design_v2.csv"
    _write_design_csv(path)

    specs = load_v2_family_specs_from_design_csv(path)

    assert {spec.family_name for spec in specs} >= {
        "rank_gap_more_protocol_columns_v2",
        "failed_controls_v2",
    }


def test_load_v2_family_specs_from_remediation_plan(tmp_path: Path) -> None:
    plan = build_preregistered_remediation_plan([], [])
    path = tmp_path / "remediation_plan_m36.json"
    path.write_text(
        json.dumps(remediation_plan_to_jsonable(plan), sort_keys=True),
        encoding="utf-8",
    )

    specs = load_v2_family_specs_from_remediation_plan(path)

    assert specs
    assert all(spec.execution_status == "planned_only" for spec in specs)


def test_mark_v2_family_specs_executed(tmp_path: Path) -> None:
    path = tmp_path / "new_manifest_family_design_v2.csv"
    _write_design_csv(path)
    specs = load_v2_family_specs_from_design_csv(path)

    executed = mark_v2_family_specs_executed(specs)

    assert all(spec.execution_status == "executed" for spec in executed)
    assert all(spec.execution_status == "planned_only" for spec in specs)
