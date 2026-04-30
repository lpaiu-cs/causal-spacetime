"""Run language guard checks for protocol-invariance terminology."""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="Protocol-invariance language audit.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def _load_banned_phrases() -> tuple[str, ...]:
    path = Path("scripts/check_theory_language.py")
    spec = importlib.util.spec_from_file_location("check_theory_language", path)
    if spec is None or spec.loader is None:
        return ()
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return tuple(getattr(module, "BANNED_PHRASES", ()))


def run_experiment() -> list[dict[str, float | str]]:
    banned_phrases = _load_banned_phrases()
    doc_guard = subprocess.run(
        [sys.executable, "scripts/check_theory_language.py"],
        check=False,
        capture_output=True,
        text=True,
    )
    python_guard = subprocess.run(
        [sys.executable, "scripts/check_theory_language.py", "--include-python"],
        check=False,
        capture_output=True,
        text=True,
    )
    phrase_a = "mixed " + "protocol " + "profile"
    phrase_b = "profile " + "distance"
    phrase_c = "patched " + "v3 family " + "will pass"
    return [
        {
            "check": "doc_language_guard_passed",
            "passed": float(doc_guard.returncode == 0),
            "detail": doc_guard.stdout.strip() or doc_guard.stderr.strip(),
        },
        {
            "check": "python_language_guard_passed",
            "passed": float(python_guard.returncode == 0),
            "detail": python_guard.stdout.strip() or python_guard.stderr.strip(),
        },
        {
            "check": "profile_mixing_phrases_configured",
            "passed": float(
                phrase_a in banned_phrases
                and phrase_b in banned_phrases
                and phrase_c in banned_phrases
            ),
            "detail": "risk phrase configuration checked",
        },
    ]


def main() -> None:
    config = parse_args()
    rows = run_experiment()
    path = write_csv(
        rows,
        data_path(config.output_dir, "protocol_invariance_language_audit.csv"),
        ["check", "passed", "detail"],
    )
    print(f"Wrote protocol-invariance language audit: {path}")


if __name__ == "__main__":
    main()
