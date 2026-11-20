"""Post hoc: are the check's single-threaded refits the stored fits up to axis signs?

The frozen check (check_intrinsic_weyl_calibration.py, f4854be) compared the
refits with the stored pure fits field by field and failed. The gauge fixes the
first landmark in the (t, x) plane and the second in the (t, x, y) hyperplane
but not the signs of the three spatial axes, which `framed` fixes only when the
statistic is formed. This script, written after that outcome, finds for each
refit the diagonal sign pattern that best maps the stored coordinates onto the
new ones and compares every tensor field under that map. It decides nothing.
"""

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from check_intrinsic_weyl_bombelli_benchmark import framed, statistic
from check_intrinsic_weyl_weak_accuracy import read
from intrinsic_weyl_joint import axis_signs
from replay_intrinsic_weyl_threads import write_json


def compare(new, old):
    a, b = (np.asarray(f["coordinates_inferred_from_counts"]) for f in (new, old))
    gap, signs = min(
        (float(abs(a - b * np.r_[1, s]).max()), s)
        for s in itertools.product((1, -1), repeat=3)
    )
    m = np.asarray(axis_signs(list(signs)), dtype=float)
    c_new, c_old = (np.asarray(f["coefficients"]) for f in (new, old))
    s_new, s_old = (np.asarray(f["nominal_component_covariance"]) for f in (new, old))
    framed_new, framed_old = framed(new), framed(old)
    mirror = np.asarray(axis_signs([1, 1, -1]), dtype=float)
    return {
        "axis_signs": list(signs),
        "coordinates_max_difference": gap,
        "coefficients_max_relative_difference": float(
            abs(c_new - m * c_old).max() / abs(c_old).max()
        ),
        "covariance_max_difference": float(abs(s_new - np.outer(m, m) * s_old).max()),
        "framed_coefficients_max_difference_up_to_z": float(
            min(
                abs(framed_new[0] - framed_old[0]).max(),
                abs(framed_new[0] - mirror * framed_old[0]).max(),
            )
        ),
        "q_new": statistic(new),
        "q_stored": statistic(old),
        "curved_cost_relative_difference": abs(new["curved_cost"] - old["curved_cost"])
        / old["curved_cost"],
    }


def main(folder, receipt, output):
    audit = read(receipt)
    rows = []
    for refit in audit["refits"]:
        r, n = refit["replicate"], refit["measurement_n"]
        target = folder.with_name(folder.name + "_audit")
        new = read(sorted(target.glob(f"case_{r}_n{n}_*_pure.json"))[-1])["fits"]
        old = read(folder / f"case_{r}_n{n}_pure.json")["fits"]
        rows.append(
            {
                "replicate": r,
                "measurement_n": n,
                **{order: compare(new[order], old[order]) for order in old},
            }
        )
    write_json(
        output,
        {
            "analysis": (
                "post hoc gauge comparison of the check refits; decides nothing"
            ),
            "receipt_passed": audit["passed"],
            "rows": rows,
        },
    )
    for row in rows:
        print(json.dumps(row))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", type=Path)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    main(args.folder.resolve(), args.receipt, args.output)
