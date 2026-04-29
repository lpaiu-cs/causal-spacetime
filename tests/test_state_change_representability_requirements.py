from __future__ import annotations

from causal_spacetime_lab.state_change_representability_requirements import (
    default_response_representability_ladder,
    representability_ladder_table,
)


def test_representability_ladder_table_includes_expected_levels() -> None:
    rows = representability_ladder_table()
    levels = {row["level"] for row in rows}

    assert "scalar_response_rank" in levels
    assert "multi_reference_response_profile" in levels
    assert "calibrated_metric_representation" in levels


def test_scalar_level_denies_target_target_distance_order() -> None:
    scalar_level = default_response_representability_ladder()[0]

    assert scalar_level.level == "scalar_response_rank"
    assert "target-target distance order" in scalar_level.what_it_does_not_imply
