"""Run every audit check script and record the outcome with a commit stamp.

Usage, from the repository root, on a clean worktree:

    python docs/paper/paper_a/order_wave_energy/independent_audit/run_all_checks.py

The runner refuses to start unless `git status --porcelain` is empty
(untracked files count as dirty), so that the `code_version` it records is a
commit that exists and contains exactly the scripts that ran. It writes
`run_logs/summary.json` and one stdout log per script next to this file.
Those outputs are committed afterwards, in a later commit than the stamp.

All scripts are CPU only. One script (proof2 s6) imports the author's module
and needs SciPy, so it runs under a second interpreter given by
`--scipy-python`.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import json
import pathlib
import re
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
HELPERS = {
    "audit_common.py",
    "audit_geometry.py",
    "toy_common.py",
    "wavegeom.py",
    "popsignal.py",
}
NEEDS_SCIPY = {("proof2_shared_reference", "s6_author_routine_crosscheck.py")}
FOLDERS = (
    "proof1_lower_bound",
    "proof2_shared_reference",
    "proof3a_null_noise",
    "proof3b_wave_signal",
)
LINE = re.compile(r"^\s*(PASS|FAIL)\b")


def git(*args: str) -> str:
    out = subprocess.run(
        ["git", *args], cwd=HERE, capture_output=True, text=True, check=True
    )
    return out.stdout.strip()


def discover() -> list[tuple[str, str]]:
    jobs = [("", "orchestrator_checks.py")]
    for folder in FOLDERS:
        for path in sorted((HERE / folder).glob("*.py")):
            if path.name not in HELPERS:
                jobs.append((folder, path.name))
    return jobs


def run_one(job: tuple[str, str], python: str, scipy_python: str) -> dict:
    folder, name = job
    cwd = HERE / folder if folder else HERE
    exe = scipy_python if job in NEEDS_SCIPY else python
    start = time.perf_counter()
    try:
        proc = subprocess.run(
            [exe, name],
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=900,
            encoding="utf-8",
            errors="replace",
        )
        code, out, err = proc.returncode, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired as exc:
        code, out, err = -9, exc.stdout or "", "TIMEOUT"
    seconds = time.perf_counter() - start
    marks = [m.group(1) for m in map(LINE.match, out.splitlines()) if m]
    label = f"{folder}__{name}" if folder else name
    log = HERE / "run_logs" / (label.removesuffix(".py") + ".log")
    log.write_text(out + ("\n[stderr]\n" + err if err.strip() else ""), newline="\n")
    return {
        "folder": folder or ".",
        "script": name,
        "interpreter": "scipy" if job in NEEDS_SCIPY else "system",
        "exit_code": code,
        "pass": marks.count("PASS"),
        "fail": marks.count("FAIL"),
        "seconds": round(seconds, 1),
        "log": log.relative_to(HERE).as_posix(),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument(
        "--scipy-python",
        default=sys.executable,
        help="interpreter with NumPy and SciPy for the author-module cross-check",
    )
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    dirty = git("status", "--porcelain")
    if dirty:
        print("refusing to run: worktree is not clean (untracked counts as dirty)")
        print(dirty)
        return 2
    stamp = git("rev-parse", "HEAD")

    (HERE / "run_logs").mkdir(exist_ok=True)
    jobs = discover()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        rows = list(
            pool.map(lambda job: run_one(job, args.python, args.scipy_python), jobs)
        )

    def version(exe: str) -> str:
        code = (
            "import sys, importlib\n"
            "mods = ['numpy', 'sympy', 'mpmath', 'scipy']\n"
            "out = [sys.version.split()[0]]\n"
            "for m in mods:\n"
            "    try:\n"
            "        out.append(m + ' ' + importlib.import_module(m).__version__)\n"
            "    except Exception:\n"
            "        out.append(m + ' absent')\n"
            "print('; '.join(out))\n"
        )
        return subprocess.run(
            [exe, "-c", code], capture_output=True, text=True
        ).stdout.strip()

    summary = {
        "code_version": stamp,
        "code_version_short": stamp[:7],
        "worktree_clean_at_start": True,
        "run_date": dt.date.today().isoformat(),
        "device": "cpu",
        "interpreters": {
            "system": version(args.python),
            "scipy": version(args.scipy_python),
        },
        "totals": {
            "scripts": len(rows),
            "pass": sum(r["pass"] for r in rows),
            "fail": sum(r["fail"] for r in rows),
            "nonzero_exit": sum(1 for r in rows if r["exit_code"] != 0),
        },
        "rows": rows,
    }
    text = json.dumps(summary, indent=2) + "\n"
    (HERE / "run_logs" / "summary.json").write_text(text, newline="\n")
    for row in rows:
        print(
            f"{row['folder']:<26} {row['script']:<36} exit={row['exit_code']} "
            f"PASS={row['pass']} FAIL={row['fail']} {row['seconds']}s"
        )
    print(json.dumps(summary["totals"]))
    bad = summary["totals"]["fail"] or summary["totals"]["nonzero_exit"]
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
