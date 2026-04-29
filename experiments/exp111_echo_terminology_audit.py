"""Deterministic terminology audit for echo-order documentation and code."""

from __future__ import annotations

import argparse
import csv
import sys
from importlib import import_module
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

DEFAULT_OUTPUT = Path("outputs/data/echo_terminology_audit.csv")

DEPRECATED_TERMS = (
    "observer-like " + "chain",
    "observer-like " + "ordered protocol",
    "observer " + "quality",
    "observer " + "quality score",
    "observer " + "selection ambiguity",
    "best " + "observer",
    "true " + "observer",
    "R" + "_rank",
)


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Echo terminology audit.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / "echo_terminology_audit.csv"


def run_experiment(root: Path | None = None) -> list[dict[str, float | str]]:
    """Run deterministic terminology checks."""

    guard = import_module("scripts.check_theory_language")
    project_root = root or Path(__file__).resolve().parents[1]
    markdown_files = guard.default_markdown_files(project_root)
    python_files = guard.default_python_files(project_root)
    markdown_violations = guard.find_language_violations(markdown_files)
    python_violations = guard.find_language_violations(markdown_files + python_files)
    rows: list[dict[str, float | str]] = [
        {
            "check": "markdown_language_guard",
            "violation_count": float(len(markdown_violations)),
            "passed": float(len(markdown_violations) == 0),
            "detail": "README.md and docs/theory/*.md",
        },
        {
            "check": "python_language_guard",
            "violation_count": float(len(python_violations)),
            "passed": float(len(python_violations) == 0),
            "detail": "markdown plus src/ and experiments/",
        },
    ]
    for phrase in DEPRECATED_TERMS:
        phrase_violations = guard.find_language_violations(
            markdown_files + python_files,
            banned_phrases=(phrase,),
        )
        rows.append(
            {
                "check": f"deprecated_phrase:{phrase}",
                "violation_count": float(len(phrase_violations)),
                "passed": float(len(phrase_violations) == 0),
                "detail": phrase,
            }
        )
    banned_lower = {phrase.lower() for phrase in guard.BANNED_PHRASES}
    missing_from_guard = [
        phrase for phrase in DEPRECATED_TERMS if phrase.lower() not in banned_lower
    ]
    rows.append(
        {
            "check": "deprecated_terms_registered",
            "violation_count": float(len(missing_from_guard)),
            "passed": float(len(missing_from_guard) == 0),
            "detail": ";".join(missing_from_guard),
        }
    )
    return rows


def write_outputs(
    rows: list[dict[str, float | str]],
    output_path: Path = DEFAULT_OUTPUT,
) -> Path:
    """Write audit rows."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def main() -> None:
    output_path = write_outputs(run_experiment(), parse_args())
    print(f"Wrote echo terminology audit: {output_path}")


if __name__ == "__main__":
    main()
