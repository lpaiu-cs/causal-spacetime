"""Verify the frozen correlation challenge and summarize its primary inverse."""

import hashlib
import json
import subprocess
import sys
from pathlib import Path, PureWindowsPath

import numpy as np

root = Path(__file__).resolve().parents[2]
folder = Path(sys.argv[1])
manifest = json.loads((folder / "manifest.json").read_text())
raw = (folder / "rows.jsonl").read_bytes()
assert manifest["status"] == "completed"
assert hashlib.sha256(raw).hexdigest() == manifest["rows_sha256"]
stamp = manifest["code_version"]
subprocess.run(
    ["git", "merge-base", "--is-ancestor", stamp, "HEAD"], cwd=root, check=True
)
for name, digest in manifest["files"].items():
    path = PureWindowsPath(name).as_posix()
    data = subprocess.check_output(["git", "show", f"{stamp}:{path}"], cwd=root)
    assert hashlib.sha256(data).hexdigest() == digest
rows = [json.loads(line) for line in raw.splitlines()]
assert len(rows) == 13
assert sum(row["stage"] == "exact_inverse" for row in rows) == 6
table = [
    dict(
        case=row["case"],
        q=row["q"],
        rmse=score["normalized_rmse"],
        max_error=score["normalized_max_error"],
    )
    for row in rows
    if row["stage"] == "exact_inverse"
    for score in row["scores"]
    if score["degree"] == 18
]
main = next(
    row
    for row in rows
    if row["stage"] == "exact_inverse"
    and row["case"] == "rotating"
    and row["q"] == 0.00125
)
refined = next(row for row in rows if row["stage"] == "refinement")
delta_p = np.max(np.abs(np.array(main["difference"]) - np.array(refined["difference"])))
delta_s = (
    np.max(
        np.abs(
            np.array(main["scores"][1]["recovered"])
            - np.array(refined["scores"][1]["recovered"])
        )
    )
    / 0.00125
)
conditioning = next(
    row for row in rows if row["stage"] == "conditioning" and row["degree"] == 18
)
summary = {
    "source_commit": stamp,
    "raw_sha256": manifest["rows_sha256"],
    "elapsed_seconds": manifest["elapsed_seconds"],
    "primary_degree": 18,
    "table": table,
    "refinement_max_p2_difference": float(delta_p),
    "refinement_max_structure_difference_over_q": float(delta_s),
    "refinement_scores": refined["scores"],
    "noise_row_norm_range": [
        min(conditioning["row_norms"]),
        max(conditioning["row_norms"]),
    ],
    "pointwise_planning_events_per_scale_range": [
        min(conditioning["normal_planning_N_per_scale_for_each_lag"]),
        max(conditioning["normal_planning_N_per_scale_for_each_lag"]),
    ],
}
(folder / "summary.json").write_text(json.dumps(summary, indent=2))
print(
    json.dumps({k: v for k, v in summary.items() if k != "refinement_scores"}, indent=2)
)
