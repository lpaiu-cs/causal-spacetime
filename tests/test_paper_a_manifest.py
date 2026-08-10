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
    """SHA-256 of the raw file bytes. Every hashed artifact is pinned
    to LF at the storage boundary (.gitattributes, eol=lf), so the
    working-tree bytes equal the committed blob bytes on every platform
    and this digest is what plain sha256sum reports on any checkout.
    Hashing raw bytes keeps the contract sensitive to every byte,
    including a line-ending change."""

    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    """Every Section 6 capstone bundle row and every Section 6.7
    extension bundle row (repo-relative path | sha256): the file
    exists and its digest matches. 13 capstone + 11 extension rows."""

    text = MANIFEST.read_text(encoding="utf-8")
    rows = re.findall(
        r"^\| `((?:docs|experiments|tests)/[^`]+)` \| `([0-9a-f]{64})` \|",
        text, re.M)
    assert len(rows) == 24, [r[0] for r in rows]
    assert len({r[0] for r in rows}) == 24
    for rel, digest in rows:
        path = REPO / rel
        assert path.exists(), rel
        assert _sha(path) == digest, rel
