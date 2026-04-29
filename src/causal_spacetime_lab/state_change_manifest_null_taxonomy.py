"""Null taxonomy for frozen-manifest family comparisons."""

from __future__ import annotations

from dataclasses import asdict, dataclass

NULL_TAXONOMY_CLASSES = {
    "destructive_null",
    "symmetry_control",
    "marginal_preserving_null",
}


@dataclass(frozen=True)
class NullTaxonomyEntry:
    """Classification of a representation null or control."""

    null_type: str
    taxonomy_class: str
    description: str
    expected_behavior: str

    def __post_init__(self) -> None:
        if self.taxonomy_class not in NULL_TAXONOMY_CLASSES:
            allowed = ", ".join(sorted(NULL_TAXONOMY_CLASSES))
            raise ValueError(f"taxonomy_class must be one of: {allowed}")


def default_null_taxonomy() -> list[NullTaxonomyEntry]:
    """Return the default null taxonomy for family comparisons."""

    return [
        NullTaxonomyEntry(
            "shuffled_sides",
            "destructive_null",
            "Randomly swaps constraint sides.",
            "Usually worsens representation fit.",
        ),
        NullTaxonomyEntry(
            "random_constraints",
            "destructive_null",
            "Samples unrelated response-comparison constraints.",
            "Usually worsens representation fit.",
        ),
        NullTaxonomyEntry(
            "permuted_targets",
            "symmetry_control",
            "Applies a target-label isomorphism to all constraints.",
            "May behave similarly to structured constraints.",
        ),
        NullTaxonomyEntry(
            "shuffle_delays",
            "marginal_preserving_null",
            "Shuffles delay ranks within protocol columns.",
            "Can be a harder marginal-preserving baseline.",
        ),
        NullTaxonomyEntry(
            "shuffle_reachability",
            "marginal_preserving_null",
            "Shuffles reachability masks within protocol columns.",
            "Tests sensitivity to reachability structure.",
        ),
        NullTaxonomyEntry(
            "random_same_marginals",
            "marginal_preserving_null",
            "Randomizes profiles while preserving per-column marginals.",
            "Can preserve coarse distributions while disrupting profiles.",
        ),
    ]


def classify_null_type(null_type: str) -> str:
    """Return taxonomy class for a null type, or ``unknown``."""

    for entry in default_null_taxonomy():
        if entry.null_type == null_type:
            return entry.taxonomy_class
    return "unknown"


def null_taxonomy_table() -> list[dict[str, str]]:
    """Return the default null taxonomy as dictionaries."""

    return [asdict(entry) for entry in default_null_taxonomy()]
