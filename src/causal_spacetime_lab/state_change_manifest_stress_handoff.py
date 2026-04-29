"""Future stress-test handoff planning from carry-forward registries."""

from __future__ import annotations

from dataclasses import asdict, dataclass

from causal_spacetime_lab.state_change_manifest_carry_forward import (
    CarryForwardRegistry,
)


@dataclass(frozen=True)
class StressTestPlanEntry:
    """Planned future stress-test allowance for one family."""

    family_name: str
    decision: str
    allowed: bool
    allowed_tests: list[str]
    required_caveats: list[str]
    stop_if_failed: bool


def default_allowed_future_stress_tests(decision: str) -> list[str]:
    """Return allowed future stress tests for a carry-forward decision."""

    if decision == "carry_forward":
        return [
            "constraint_sparsity",
            "harder_nulls",
            "optimizer_stability",
            "protocol_transfer",
        ]
    if decision == "provisional":
        return ["constraint_sparsity", "optimizer_stability"]
    return []


def build_stress_test_plan(
    registry: CarryForwardRegistry,
) -> list[StressTestPlanEntry]:
    """Build future stress-test plan entries from a registry."""

    entries: list[StressTestPlanEntry] = []
    for record in registry.records:
        allowed_tests = default_allowed_future_stress_tests(record.decision)
        caveats: list[str] = []
        if record.decision == "provisional":
            caveats.append("provisional_family_requires_explicit_caveats")
        if record.decision in {"blocked", "report_only", "failed_control"}:
            caveats.append("report_only_no_future_stress_tests")
        entries.append(
            StressTestPlanEntry(
                family_name=record.family_name,
                decision=record.decision,
                allowed=bool(allowed_tests),
                allowed_tests=allowed_tests,
                required_caveats=caveats,
                stop_if_failed=record.decision == "carry_forward",
            )
        )
    return entries


def stress_test_plan_table(
    entries: list[StressTestPlanEntry],
) -> list[dict[str, float | str]]:
    """Convert stress-test plan entries to CSV rows."""

    rows: list[dict[str, float | str]] = []
    for entry in entries:
        row = asdict(entry)
        row["allowed"] = float(entry.allowed)
        row["allowed_tests"] = ";".join(entry.allowed_tests)
        row["required_caveats"] = ";".join(entry.required_caveats)
        row["stop_if_failed"] = float(entry.stop_if_failed)
        rows.append(row)
    return rows
