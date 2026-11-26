# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of quantized_order_bound_r2.md.
"""Section-by-section diff of quantized_order_bound.md (orig) vs _r2.md."""

import difflib
import re
import sys
from pathlib import Path

here = Path(__file__).parent
orig = (here / "orig.md").read_text(encoding="utf-8").splitlines()
r2 = (here / "r2.md").read_text(encoding="utf-8").splitlines()


def sections(lines):
    out, name, buf = {}, "preamble", []
    for ln in lines:
        m = re.match(r"^## (\d+)\.", ln)
        if m:
            out[name] = buf
            name, buf = m.group(1), []
        buf.append(ln)
    out[name] = buf
    return out


so, sr = sections(orig), sections(r2)
which = sys.argv[1:] or sorted(set(so) | set(sr))
for key in which:
    a, b = so.get(key, []), sr.get(key, [])
    sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
    print(f"===== section {key}: orig {len(a)} lines, r2 {len(b)} lines, ratio {sm.ratio():.3f}")
    if a == b:
        print("IDENTICAL")
        continue
    for line in difflib.unified_diff(a, b, "orig", "r2", n=0, lineterm=""):
        print(line)
