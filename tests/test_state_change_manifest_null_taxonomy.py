from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_null_taxonomy import (
    classify_null_type,
    default_null_taxonomy,
    null_taxonomy_table,
)


def test_classify_null_type_maps_shuffled_sides_to_destructive_null() -> None:
    assert classify_null_type("shuffled_sides") == "destructive_null"


def test_classify_null_type_maps_permuted_targets_to_symmetry_control() -> None:
    assert classify_null_type("permuted_targets") == "symmetry_control"


def test_default_null_taxonomy_contains_expected_entries() -> None:
    entries = default_null_taxonomy()
    table = null_taxonomy_table()

    assert any(entry.null_type == "random_same_marginals" for entry in entries)
    assert any(row["null_type"] == "shuffle_delays" for row in table)
