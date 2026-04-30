"""Check theory docs for risky overclaim language."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

BANNED_PHRASES = (
    "time is information-transfer speed",
    "space is information-transfer speed",
    "spacetime is information-transfer speed",
    "causal order implies global discrete time",
    "distance equals c times time slice",
    "distance = c times time slice",
    "causal order alone gives metric geometry",
    "this derives quantum mechanics",
    "this replaces relativity",
    "space is made of absolute cells",
    "we proved spacetime emergence",
    "generation order is physical time",
    "event id is time",
    "local finiteness implies global time",
    "state-change network derives relativity",
    "state-change network derives quantum mechanics",
    "observer-like chain",
    "observer quality",
    "observer quality score",
    "observer selection ambiguity",
    "true observer",
    "best observer",
    "observer derived uniquely",
    "observer is derived uniquely",
    "best chain is the true observer",
    "chain score proves observerhood",
    "chain interval profile is physical time",
    "chain interval count is seconds",
    "observer-chain selection derives metric",
    "observer-like ordered protocol",
    "bracket width is distance",
    "rank width is distance",
    "radar rank is metric time",
    "bracket rank gives meters",
    "chain brackets derive spatial geometry",
    "order-level bracket proves observerhood",
    "echo delay is distance",
    "echo rank gives distance",
    "echo order gives spatial geometry",
    "finite-speed spatial geometry is implemented",
    "same-emission echo proves distance",
    "echo protocol derives metric",
    "planted delay is distance",
    "echo motif gives distance",
    "response motif derives geometry",
    "echo-response motif is finite-speed geometry",
    "planted delay is physical time",
    "shortcut is metric noise",
    "background interference is metric noise",
    "return spectrum gives distance",
    "shortcut depth is distance",
    "shortcut depth is speed",
    "echo interference derives geometry",
    "shortcut classification proves spatial distance",
    "coarse-graining derives geometry",
    "reference subsampling calibrates time",
    "edge thinning is metric noise",
    "shortcut depth is metric error",
    "response order is distance",
    "response signature is spatial geometry",
    "stable response order proves geometry",
    "echo response order derives metric",
    "stable response core is distance",
    "response layer is physical space",
    "R_rank",
    "echo distance",
    "response distance",
    "D_echo is distance",
    "D_echo is duration",
    "D_echo is time",
    "gated echo fixes shortcut",
    "stable response order is distance order",
    "scalar representability proves geometry",
    "scalar rank is metric",
    "single response order gives distance",
    "response order defines pairwise distance",
    "scalar response order is distance order",
    "multi-reference profile is metric",
    "response profile proves geometry",
    "representability ladder proves metric",
    "scalar rank representation is spatial geometry",
    "pairwise response dissimilarity is distance",
    "pairwise response comparison is distance order",
    "response-profile dissimilarity is spatial distance",
    "response profile recovers pairwise distance",
    "null baselines prove geometry",
    "pairwise response protocol derives space",
    "validated constraints are distance",
    "validated constraints are distance-order constraints",
    "constraint pool recovers geometry",
    "validation gates prove geometry",
    "held-out agreement proves metric",
    "null baseline separation proves space",
    "pre-embedding validation reconstructs distance",
    "handoff proves geometry",
    "handoff manifest proves representability",
    "eligible pool is distance order",
    "eligible constraints are spatial distances",
    "manifest contains metric reconstruction",
    "handoff exports embedding",
    "failed handoff is ignored",
    "fitted coordinates are spatial coordinates",
    "manifest embedding recovers geometry",
    "ordinal representation proves metric",
    "heldout violation proves space",
    "latent representation is physical distance",
    "dimension curve gives spacetime dimension",
    "null fit proves geometry",
    "manifest fit reconstructs distance",
    "family comparison proves geometry",
    "manifest family recovers space",
    "dimension selected is physical dimension",
    "permuted target null destroys geometry",
    "target permutation null proves structure",
    "failed manifests ignored",
    "stricter criteria prove metric",
    "carry-forward proves geometry",
    "family robustness proves geometry",
    "family-level fit recovers space",
    "carry-forward family recovers distance",
    "stress-test eligibility proves metric",
    "selected family is spatial geometry",
    "blocked families ignored",
    "threshold sensitivity chooses best threshold",
    "failure analysis retunes thresholds",
    "missing metric counts as pass",
    "blocked family can be stress-tested anyway",
    "stop condition ignored",
    "remediation recovers geometry",
    "failure decomposition proves metric",
    "thresholds changed after failure",
    "failed families dropped",
    "remediation plan fixes geometry",
    "remediation recovers metric",
    "planned v2 family is current result",
    "future run can execute now",
    "thresholds lowered after failure",
    "blocked family stress-tested anyway",
    "remediation plan exports new fit",
    "v2 manifests prove geometry",
    "diagnostic complete means successful",
    "v2 generation evaluates carry-forward",
    "v2 generation runs stress tests",
    "v2 manifests recover metric",
    "diagnostic-complete bundle proves space",
    "v2 family passes by construction",
    "v2 carry-forward proves geometry",
    "v2 family recovers metric",
    "v2 carry-forward recovers space",
    "v2 decision runs stress tests",
    "v2 decision reruns fits",
    "diagnostic completeness means carry-forward",
    "carry-forward means physical coordinates",
    "threshold sensitivity chooses threshold",
    "v3 design fixes geometry",
    "v3 family will pass",
    "structural counterfactual changes decision",
    "ignore manifest count failure",
    "thresholds adjusted for v3",
    "v3 preregistration executes manifests",
    "blocked v2 family can be stress-tested",
    "measured failure is only structural",
    "mixed protocol profile",
    "profile distance",
    "response-profile distance",
    "measurement-rule-mixed dissimilarity",
    "D_echo distance",
    "protocol-mixed profile is admissible",
    "mixed protocols can be concatenated",
    "measurement protocol variation is reference variation",
    "v3 protocol patch proves geometry",
    "protocol-invariant profile is distance",
    "patched v3 family will pass",
    "v3 protocol manifests prove geometry",
    "protocol metadata proves geometry",
    "parameter-complete protocol is distance",
    "top-down handoff proves geometry",
    "top-down manifest is evidence",
    "top-down design recovers metric",
    "handoff provenance proves theory",
    "v3 protocol generation means success",
    "v3 diagnostic complete means carry-forward",
    "v3 protocol generation evaluates carry-forward",
    "v3 protocol generation runs stress tests",
    "v3 manifests recover metric",
    "protocol-invariant family passes by construction",
    "parameter metadata recovers metric",
    "v3 protocol carry-forward proves geometry",
    "v3 protocol family recovers metric",
    "v3 carry-forward recovers space",
    "top-down carry-forward proves geometry",
    "hybrid handoff proves theory",
    "provenance integrity proves metric",
    "protocol metadata means physical coordinates",
    "v3 decision runs stress tests",
    "v3 decision reruns fits",
    "v4 design fixes geometry",
    "v4 family will pass",
    "thresholds adjusted for v4",
    "blocked v3 ignored",
    "top-down provenance proves structure",
    "hybrid provenance proves theory",
    "counterfactual changes decision",
    "report-only counterfactual justifies threshold",
    "v3 failure proves theory false",
    "v4 remediation proves metric",
    "v4 generation means success",
    "v4 diagnostic complete means carry-forward",
    "v4 manifests prove geometry",
    "v4 manifests recover metric",
    "v4 protocol manifests prove geometry",
    "v4 execution evaluates carry-forward",
    "v4 generation runs stress tests",
    "v4 family passes by construction",
    "v4 fixes theory",
    "diagnostic-complete v4 proves space",
)


@dataclass(frozen=True)
class LanguageViolation:
    """A risky phrase occurrence in a theory-facing document."""

    path: Path
    line_number: int
    phrase: str
    line: str


def default_markdown_files(root: Path) -> list[Path]:
    """Return the theory-facing markdown files checked by default."""

    files = [root / "README.md"]
    files.extend(sorted((root / "docs" / "theory").glob("*.md")))
    return [path for path in files if path.exists()]


def default_python_files(root: Path) -> list[Path]:
    """Return source and experiment Python files for optional language checks."""

    files = sorted((root / "src" / "causal_spacetime_lab").glob("*.py"))
    files.extend(sorted((root / "experiments").glob("*.py")))
    files.extend(sorted((root / "scripts").glob("*.py")))
    files.extend(sorted((root / "docs").glob("*.md")))
    files.append(root / "README.md")
    return [path for path in files if path.exists()]


def _heading_level(line: str) -> int | None:
    stripped = line.lstrip()
    if not stripped.startswith("#"):
        return None
    return len(stripped) - len(stripped.lstrip("#"))


def _allowed_rejected_contexts(lines: list[str]) -> list[bool]:
    """Mark lines inside explicit rejected-language sections."""

    allowed = [False] * len(lines)
    active_level: int | None = None
    for index, line in enumerate(lines):
        level = _heading_level(line)
        if level is not None:
            if active_level is not None and level <= active_level:
                active_level = None
            lower = line.lower()
            if "rejected language" in lower or "deprecated terminology" in lower:
                active_level = level
        allowed[index] = active_level is not None
    return allowed


def _allowed_python_forbidden_contexts(path: Path, lines: list[str]) -> list[bool]:
    """Mark lines inside explicit Python forbidden-interpretation helpers."""

    allowed = [False] * len(lines)
    if path.suffix != ".py":
        return allowed
    if path.name == "check_theory_language.py":
        active_tuple = False
        for index, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith("BANNED_PHRASES = ("):
                active_tuple = True
            allowed[index] = active_tuple
            if active_tuple and stripped == ")":
                active_tuple = False
        return allowed
    active_indent: int | None = None
    for index, line in enumerate(lines):
        stripped = line.lstrip()
        indent = len(line) - len(stripped)
        if active_indent is not None and stripped and indent <= active_indent:
            active_indent = None
        if stripped.startswith("def forbidden_"):
            active_indent = indent
        allowed[index] = active_indent is not None
    return allowed


def find_language_violations(
    files: list[Path],
    banned_phrases: tuple[str, ...] = BANNED_PHRASES,
) -> list[LanguageViolation]:
    """Find risky phrases outside allowed rejected-language contexts."""

    violations: list[LanguageViolation] = []
    for path in files:
        lines = path.read_text(encoding="utf-8").splitlines()
        allowed_context = _allowed_rejected_contexts(lines)
        python_forbidden_context = _allowed_python_forbidden_contexts(path, lines)
        path_is_rejected_claims = path.name == "rejected_claims.md"
        for line_number, line in enumerate(lines, start=1):
            lower = line.lower()
            for phrase in banned_phrases:
                pattern = (
                    r"(?<![a-z0-9_])"
                    + re.escape(phrase.lower())
                    + r"(?![a-z0-9_])"
                )
                if re.search(pattern, lower):
                    if (
                        path_is_rejected_claims
                        or allowed_context[line_number - 1]
                        or python_forbidden_context[line_number - 1]
                    ):
                        continue
                    violations.append(
                        LanguageViolation(
                            path=path,
                            line_number=line_number,
                            phrase=phrase,
                            line=line.strip(),
                        )
                    )
    return violations


def main() -> int:
    """Run the language guard on README and docs/theory."""

    parser = argparse.ArgumentParser(description="Check theory-facing language.")
    parser.add_argument(
        "--include-python",
        action="store_true",
        help="also scan src/causal_spacetime_lab/*.py and experiments/*.py",
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    files = default_markdown_files(root)
    if args.include_python:
        files.extend(default_python_files(root))
    violations = find_language_violations(files)
    if violations:
        print("Theory language guard found risky phrases:")
        for violation in violations:
            rel_path = violation.path.relative_to(root)
            print(
                f"{rel_path}:{violation.line_number}: "
                f"{violation.phrase!r}: {violation.line}"
            )
        return 1
    print("Theory language guard passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
