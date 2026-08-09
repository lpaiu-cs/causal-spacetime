"""Contract test for the Paper A evidence manifest.

The manifest's hash tables are the paper's provenance claim; this test
recomputes every SHA-256 (and the row counts of the 19-table package)
against the committed files, so the manifest and the evidence cannot
drift apart silently -- the pattern the P14 capstone bundle
established, applied to the legacy package.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / "docs" / "paper" / "paper_a" / "artifact_manifest.md"
FIG_DATA = REPO / "docs" / "paper" / "paper_a" / "figures" / "data"


def _sha(path: Path) -> str:
    """SHA-256 of the LF-normalized content -- the committed blob
    bytes. Working-tree line endings vary by platform (git's
    autocrlf materializes CRLF on Windows and normalizes back to LF
    on commit), so hashing raw working-tree bytes would make the
    contract checkout-dependent."""

    return hashlib.sha256(
        path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def test_the_19_table_package_matches_its_manifest_hashes():
    """Every legacy-table row (name | rows | sha256): the file exists
    under figures/data/, its digest matches, and its data row count
    matches."""

    text = MANIFEST.read_text(encoding="utf-8")
    rows = re.findall(
        r"^\| `([^`/]+\.csv)` \| (\d+) \| `([0-9a-f]{64})` \|",
        text, re.M)
    assert len(rows) == 19, [r[0] for r in rows]
    assert len({r[0] for r in rows}) == 19
    for name, n_rows, digest in rows:
        path = FIG_DATA / name
        assert path.exists(), name
        assert _sha(path) == digest, name
        with path.open(encoding="utf-8") as f:
            data_rows = sum(1 for _ in f) - 1
        assert data_rows == int(n_rows), name


def test_the_p14_bundle_matches_its_manifest_hashes():
    """Every Section 6 capstone bundle row (repo-relative path |
    sha256): the file exists and its digest matches."""

    text = MANIFEST.read_text(encoding="utf-8")
    rows = re.findall(
        r"^\| `((?:docs|experiments|tests)/[^`]+)` \| `([0-9a-f]{64})` \|",
        text, re.M)
    assert len(rows) == 13, [r[0] for r in rows]
    for rel, digest in rows:
        path = REPO / rel
        assert path.exists(), rel
        assert _sha(path) == digest, rel
