"""Download bundled genome FASTA files from external hosting into ./genomes/.

Genome files are too large for PyPI, so they are hosted separately (GitHub
Releases / Zenodo) and fetched on demand. The manifest lives in package data.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import sys
import urllib.request
from importlib.resources import files as _pkg_files
from pathlib import Path


_MANIFEST_PATH = _pkg_files("blast_align_tree") / "data" / "genomes_manifest.json"


def load_manifest() -> dict:
    return json.loads(_MANIFEST_PATH.read_text(encoding="utf-8"))


def sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as resp, tmp.open("wb") as out:
        shutil.copyfileobj(resp, out)
    tmp.replace(dest)


def decompress_gz(src: Path, dest: Path) -> None:
    with gzip.open(src, "rb") as gz_in, dest.open("wb") as out:
        shutil.copyfileobj(gz_in, out)


def fetch_one(name: str, entry: dict, dest_dir: Path) -> bool:
    tag = f"{entry.get('emoji', '  ')} {name}"
    url = entry.get("url")
    if not url:
        print(f"  [{tag}] skipped — no url in manifest", file=sys.stderr)
        return False

    filename = entry.get("filename") or f"{name}.fa"
    final_path = dest_dir / filename

    if final_path.exists():
        expected = entry.get("sha256")
        if expected and sha256_of(final_path) == expected:
            print(f"  [{tag}] already present, checksum OK")
            return True
        if not expected:
            print(f"  [{tag}] already present (no checksum to verify)")
            return True
        print(f"  [{tag}] present but checksum mismatch — re-downloading")

    print(f"  [{tag}] downloading {url}")
    if url.endswith(".gz"):
        tmp_gz = dest_dir / (filename + ".gz")
        download(url, tmp_gz)
        decompress_gz(tmp_gz, final_path)
        tmp_gz.unlink()
    else:
        download(url, final_path)

    expected = entry.get("sha256")
    if expected:
        got = sha256_of(final_path)
        if got != expected:
            final_path.unlink(missing_ok=True)
            print(f"  [{tag}] ERROR: sha256 mismatch (got {got}, expected {expected})",
                  file=sys.stderr)
            return False
    print(f"  [{tag}] done -> {final_path}")
    return True


def cmd_list(manifest: dict) -> None:
    genomes = manifest["genomes"]
    print(f"Available genomes ({len(genomes)}):")
    for name, entry in genomes.items():
        emoji = entry.get("emoji") or "  "
        default = " (default)" if entry.get("default") else ""
        size = entry.get("size_mb")
        size_str = f" ~{size} MB" if size else ""
        url_state = "" if entry.get("url") else " [not yet hosted]"
        print(f"  {emoji} {name}{default}{size_str}{url_state}")
        if entry.get("description"):
            print(f"      {entry['description']}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="blast-align-tree-fetch",
        description="Download bundled genome databases into ./genomes/",
    )
    ap.add_argument("names", nargs="*",
                    help="Specific genome names to fetch (default: the 'default' set)")
    ap.add_argument("--all", action="store_true", help="Fetch every genome listed in the manifest")
    ap.add_argument("--list", action="store_true", help="List available genomes and exit")
    ap.add_argument("--dest", default="genomes",
                    help="Destination directory (default: ./genomes)")
    args = ap.parse_args(argv)

    manifest = load_manifest()

    if args.list:
        cmd_list(manifest)
        return 0

    genomes = manifest["genomes"]
    if args.names:
        unknown = [n for n in args.names if n not in genomes]
        if unknown:
            print(f"Unknown genomes: {', '.join(unknown)}", file=sys.stderr)
            return 2
        selected = {n: genomes[n] for n in args.names}
    elif args.all:
        selected = genomes
    else:
        selected = {n: e for n, e in genomes.items() if e.get("default")}
        if not selected:
            print("No default genomes configured in manifest; "
                  "pass names or --all.", file=sys.stderr)
            return 2

    dest_dir = Path(args.dest).resolve()
    print(f"Fetching {len(selected)} genome(s) into {dest_dir}")
    failures = 0
    for name, entry in selected.items():
        if not fetch_one(name, entry, dest_dir):
            failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
