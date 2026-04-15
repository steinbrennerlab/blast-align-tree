# v1.0.0 — Python package release

First pip-installable release. Major rewrite of the pipeline into a proper Python package, with bundled genome databases hosted as release assets.

## Highlights

- **pip-installable package** (`blast-align-tree`) with three CLIs:
  - `blast-align-tree` — the main pipeline
  - `bat-genome-selector` — Tk GUI for picking genome sets
  - `blast-align-tree-fetch` — downloads bundled genome FASTAs from this release into `./genomes/`
- **Bundled genome databases** hosted as release assets (see Downloads below). Default set covers Arabidopsis (TAIR10 CDS + proteins), common bean (Pvul218), and *N. benthamiana* (Niben261 proteins). Opt-in animal/model genomes: human, mouse, rat, chimp, zebrafish, fly, worm, yeast.
- **Tutorial rewrite**: SOBIR1 example for *N. benthamiana*, Klepikova expression overlay, MAFFT + RAxML-NG walkthrough.

## What's new since v0.1.1

### Packaging & distribution
- Package as installable `blast-align-tree` for PyPI (setuptools, entry points, package-data manifest)
- README rewrite for the pip-packaged pipeline; documented conda install
- Genome fetch manifest + sha256 verification for reproducible downloads

### Pipeline
- Integrated MAFFT7 and RAxML-NG alternatives alongside the original aligner/tree builder
- HMM motif scanner with visualization overlays
- RAxML bootstrap support
- Generalized tree dataset overlays and heatmaps
- Klepikova expression data overlay with eFP-browser-compatible headers
- Tree files moved into per-run assets folder for cleaner output

### Genome data
- Added Niben261 as a default genome (promoted from opt-in)
- Added animal/model genomes (human, mouse, rat, chimp, zebrafish, fly, worm, yeast)
- SOBIR1 used as the tobacco example query

### UI / UX
- Genome-selector GUI with column-header alignment fix
- Suffix parser for custom run names
- Cleaner terminal output and PDF legends
- `--check_env` for environment diagnostics

### Bug fixes
- UnicodeDecodeError on Windows when reading FASTA files
- Quote bundled Rscript path; suffix-aware redraw filenames
- MSA plot for wrong sequence counts
- Pipe character in identifiers
- Single-sequence hit handling
- gene_symbols version conflict

## Downloads

Fetch the default genome set after installing:

```bash
pip install blast-align-tree
blast-align-tree-fetch              # downloads the default set
blast-align-tree-fetch --all        # downloads every genome below
blast-align-tree-fetch --list       # see what's available
```

Assets are downloaded, sha256-verified, and decompressed automatically into `./genomes/` (animals sort into `./genomes/animals/`).
