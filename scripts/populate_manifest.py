"""Populate genomes_manifest.json from manifest_values.tsv.

Reads release-assets/v1.0.0/manifest_values.tsv (produced when staging
release assets) and fills the url/sha256/size_mb fields of
blast_align_tree/data/genomes_manifest.json for each genome entry.

Run from the repo root:
    python scripts/populate_manifest.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TSV = REPO / "release-assets" / "v1.0.0" / "manifest_values.tsv"
MANIFEST = REPO / "blast_align_tree" / "data" / "genomes_manifest.json"


def main() -> int:
    with TSV.open(encoding="utf-8", newline="") as f:
        rows = {r["name"]: r for r in csv.DictReader(f, delimiter="\t")}

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    missing = [n for n in manifest["genomes"] if n not in rows]
    if missing:
        print(f"ERROR: no TSV row for: {', '.join(missing)}")
        return 1

    for name, entry in manifest["genomes"].items():
        row = rows[name]
        entry["url"] = row["url"]
        entry["sha256"] = row["sha256_uncompressed"]
        entry["size_mb"] = float(row["size_mb_gz"])

    MANIFEST.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"Updated {len(manifest['genomes'])} entries in {MANIFEST.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
