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
    files.extend(sorted((root / "docs").glob("*.md")))
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


def find_language_violations(
    files: list[Path],
    banned_phrases: tuple[str, ...] = BANNED_PHRASES,
) -> list[LanguageViolation]:
    """Find risky phrases outside allowed rejected-language contexts."""

    violations: list[LanguageViolation] = []
    for path in files:
        lines = path.read_text(encoding="utf-8").splitlines()
        allowed_context = _allowed_rejected_contexts(lines)
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
                    if path_is_rejected_claims or allowed_context[line_number - 1]:
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
