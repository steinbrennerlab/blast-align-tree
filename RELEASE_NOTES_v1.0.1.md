# v1.0.1 — Auto-index fetched genomes; add Klepikova + cowpea defaults

Quality-of-life release that makes `blast-align-tree-fetch` fully ready-to-run:
downloaded genomes are indexed automatically, and the default set now includes
every file the bundled tutorial examples need.

## Highlights

- **Automatic `makeblastdb` on fetch.** `blast-align-tree-fetch` now runs
  `makeblastdb -in <fa> -dbtype <nucl|prot> -parse_seqids` after each successful
  download. Nucleotide vs. protein is auto-detected. Pass `--no-index` to skip.
  Existing `.nhr`/`.phr` indexes are detected and reused. If `makeblastdb` is
  not on `PATH`, a warning is printed and the download still succeeds.
- **Cowpea (`Vung469cds`) promoted to the default set** so the bundled
  tutorial/example commands work out of the box.
- **Klepikova Arabidopsis expression atlas subset** is now a release asset and
  a default-fetched dataset. It lands in `./datasets/` (not `./genomes/`) so
  existing overlay paths in tutorials resolve directly.
- **Per-entry `dest_dir` in the fetch manifest.** Non-genome assets (like the
  Klepikova TSV) can be routed to a different directory. Non-FASTA extensions
  skip BLAST indexing automatically.

## Upgrade

```bash
pip install --upgrade blast-align-tree
blast-align-tree-fetch          # default set now includes cowpea + Klepikova
```

## Full changelog

- `fetch.py`: automatic `makeblastdb` with `--no-index` opt-out
- `fetch.py`: per-entry `dest_dir` override; skip indexing for non-FASTA files
- Manifest: `Vung469cds` → `default: true`
- Manifest: new `klepikova_atlas` entry (default, routed to `./datasets/`)
- Release asset: `klepikova_atlas_subset.tsv`
