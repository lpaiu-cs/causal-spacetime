"""Definability table for order-first cross-slice predicates."""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class PredicateRequirement:
    """Protocol requirements for a spatial or dynamical predicate."""

    name: str
    requires_causal_order: bool
    requires_observer_protocol: bool
    requires_slice_selection: bool
    requires_persistence: bool
    requires_transport: bool
    requires_metric_calibration: bool
    defined_without_transport: bool
    description: str


def default_predicate_requirements() -> list[PredicateRequirement]:
    """Return the default conservative predicate-definability table."""

    return [
        PredicateRequirement(
            name="same_slice_distance_order",
            requires_causal_order=True,
            requires_observer_protocol=True,
            requires_slice_selection=True,
            requires_persistence=False,
            requires_transport=False,
            requires_metric_calibration=False,
            defined_without_transport=True,
            description=(
                "Within one observer-selected slice, pair-distance order can be "
                "compared without cross-slice transport."
            ),
        ),
        PredicateRequirement(
            name="cross_slice_same_position",
            requires_causal_order=True,
            requires_observer_protocol=True,
            requires_slice_selection=True,
            requires_persistence=False,
            requires_transport=True,
            requires_metric_calibration=False,
            defined_without_transport=False,
            description=(
                "Same-position across slices requires a transport, anchor, or "
                "equivalent identification rule."
            ),
        ),
        PredicateRequirement(
            name="same_direction",
            requires_causal_order=True,
            requires_observer_protocol=True,
            requires_slice_selection=True,
            requires_persistence=False,
            requires_transport=True,
            requires_metric_calibration=False,
            defined_without_transport=False,
            description=(
                "Cross-slice direction comparison requires a chosen orientation "
                "and transport rule."
            ),
        ),
        PredicateRequirement(
            name="object_persistence",
            requires_causal_order=False,
            requires_observer_protocol=False,
            requires_slice_selection=False,
            requires_persistence=True,
            requires_transport=False,
            requires_metric_calibration=False,
            defined_without_transport=True,
            description=(
                "Object identity across slices is additional persistence "
                "structure, not derived here from causal order."
            ),
        ),
        PredicateRequirement(
            name="pair_distance_order_history",
            requires_causal_order=True,
            requires_observer_protocol=True,
            requires_slice_selection=True,
            requires_persistence=True,
            requires_transport=False,
            requires_metric_calibration=False,
            defined_without_transport=True,
            description=(
                "Persistent object labels plus slice-local distance order define "
                "pair-distance order histories without same-position transport."
            ),
        ),
        PredicateRequirement(
            name="relational_shape_change",
            requires_causal_order=True,
            requires_observer_protocol=True,
            requires_slice_selection=True,
            requires_persistence=True,
            requires_transport=False,
            requires_metric_calibration=False,
            defined_without_transport=True,
            description=(
                "Changes in ordinal shape signatures are persistence-dependent "
                "but transport-independent relational statements."
            ),
        ),
        PredicateRequirement(
            name="coordinate_velocity",
            requires_causal_order=True,
            requires_observer_protocol=True,
            requires_slice_selection=True,
            requires_persistence=True,
            requires_transport=True,
            requires_metric_calibration=True,
            defined_without_transport=False,
            description=(
                "Coordinate velocity needs object persistence, time calibration, "
                "and cross-slice spatial transport."
            ),
        ),
        PredicateRequirement(
            name="constant_velocity",
            requires_causal_order=True,
            requires_observer_protocol=True,
            requires_slice_selection=True,
            requires_persistence=True,
            requires_transport=True,
            requires_metric_calibration=True,
            defined_without_transport=False,
            description=(
                "Constant-velocity judgments compare transported positions over "
                "calibrated slice times."
            ),
        ),
        PredicateRequirement(
            name="metric_spatial_evolution",
            requires_causal_order=True,
            requires_observer_protocol=True,
            requires_slice_selection=True,
            requires_persistence=True,
            requires_transport=True,
            requires_metric_calibration=True,
            defined_without_transport=False,
            description=(
                "Metric spatial evolution requires a calibrated metric "
                "representation and transport across slices."
            ),
        ),
        PredicateRequirement(
            name="quantitative_dynamics",
            requires_causal_order=True,
            requires_observer_protocol=True,
            requires_slice_selection=True,
            requires_persistence=True,
            requires_transport=True,
            requires_metric_calibration=True,
            defined_without_transport=False,
            description=(
                "Quantitative dynamics is an effective representation-layer "
                "statement, not a primitive order predicate."
            ),
        ),
    ]


def predicate_requirement_table() -> list[dict[str, str | bool]]:
    """Return predicate requirements as serializable dictionaries."""

    return [asdict(requirement) for requirement in default_predicate_requirements()]
