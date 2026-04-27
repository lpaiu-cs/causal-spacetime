from __future__ import annotations

import numpy as np

from causal_spacetime_lab.observer_slice_pipeline import (
    reconstruct_stationary_oriented_slices_1p1,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond


def test_reconstruct_stationary_oriented_slices_returns_required_keys() -> None:
    events = sprinkle_1p1_causal_diamond(40, T=2.0, seed=1)

    result = reconstruct_stationary_oriented_slices_1p1(
        events,
        T=2.0,
        tick_count=16,
        beacon_separation=0.15,
        bin_width=2,
    )

    required = {
        "accessible",
        "reconstructed_T",
        "reconstructed_X",
        "predecessor_tick_positions",
        "successor_tick_positions",
        "radar_time_rank",
        "radar_distance_rank",
        "slice_labels",
    }
    assert required <= set(result)
    assert all(np.asarray(result[key]).shape == (40,) for key in required)
    assert np.all(np.asarray(result["slice_labels"])[~result["accessible"]] == -1)
