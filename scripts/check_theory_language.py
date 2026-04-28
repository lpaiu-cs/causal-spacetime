"""Check theory docs for risky overclaim language."""

from __future__ import annotations

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
            if "rejected language" in line.lower():
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
                if phrase in lower:
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

    root = Path(__file__).resolve().parents[1]
    violations = find_language_violations(default_markdown_files(root))
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
