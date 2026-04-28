from __future__ import annotations

from pathlib import Path

from scripts.check_theory_language import (
    default_markdown_files,
    find_language_violations,
)


def test_banned_phrase_is_detected_in_temporary_doc(tmp_path: Path) -> None:
    doc = tmp_path / "bad.md"
    doc.write_text("Best chain is the true observer.\n", encoding="utf-8")

    violations = find_language_violations([doc])

    assert len(violations) == 1
    assert violations[0].phrase == "best chain is the true observer"


def test_rejected_language_section_is_allowed(tmp_path: Path) -> None:
    doc = tmp_path / "allowed.md"
    doc.write_text(
        "# Test\n\n"
        "## Rejected Language\n\n"
        '- "causal order alone gives metric geometry"\n',
        encoding="utf-8",
    )

    assert find_language_violations([doc]) == []


def test_rejected_claims_file_is_allowed(tmp_path: Path) -> None:
    doc = tmp_path / "rejected_claims.md"
    doc.write_text('- "this replaces relativity"\n', encoding="utf-8")

    assert find_language_violations([doc]) == []


def test_real_docs_pass_language_guard() -> None:
    root = Path(__file__).resolve().parents[1]

    assert find_language_violations(default_markdown_files(root)) == []
