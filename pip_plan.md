# PyPI Packaging Plan for blast-align-tree

## Goal

Publish `blast-align-tree` to PyPI so users can `pip install blast-align-tree` and
get two CLI commands: `blast-align-tree` (the pipeline) and `bat-genome-selector`
(the tkinter GUI). The R script and helper scripts ship as bundled package data.

---

## Current State

```
bat-april/                          (repo root)
    blast_align_tree.py             (1417 lines, main pipeline, has main())
    genome_selector.py              (1056 lines, tkinter GUI, has main())
    visualize_tree.r                (716 lines, called via Rscript subprocess)
    gene_symbols.txt                (data file, ~168 KB)
    scripts/                        (small helper scripts: add_seq, extract_seq, etc.)
    genomes/                        (user genome FASTA files — NOT shipped)
    datasets/                       (user data — NOT shipped)
    hmm_files/                      (user HMM files — NOT shipped)
    environments/                   (conda YAML files)
    images/                         (README images)
```

Key behaviors:
- `blast_align_tree.py` runs from the project root and expects `visualize_tree.r`
  in the workdir (i.e. the same directory).
- `genome_selector.py` uses `Path(__file__).parent / "genomes"` to find genomes.
- Both scripts have `def main()` and `if __name__ == "__main__": main()` already.
- External tools (BLAST+, MAFFT, RAxML-NG, trimAl, Rscript, etc.) are called via
  `subprocess` and must be on `PATH` — pip cannot install these.

---

## Target Structure

```
bat-april/                          (repo root)
    pyproject.toml                  (NEW — package metadata + entry points)
    LICENSE                         (NEW — required for PyPI)
    README.md                       (existing — renders as PyPI project page)
    blast_align_tree/               (NEW — Python package directory)
        __init__.py                 (NEW — version string)
        cli.py                      (MOVED from blast_align_tree.py)
        genome_selector.py          (MOVED from genome_selector.py)
        data/                       (NEW — bundled non-Python files)
            visualize_tree.r        (MOVED from visualize_tree.r)
            gene_symbols.txt        (MOVED from gene_symbols.txt)
    scripts/                        (stays as-is, not part of pip package)
    environments/                   (stays as-is, not part of pip package)
    genomes/                        (stays as-is, user data, not shipped)
    datasets/                       (stays as-is, user data, not shipped)
    hmm_files/                      (stays as-is, user data, not shipped)
```

---

## Step-by-step Plan

### Step 1: Add LICENSE file

PyPI requires a license. Pick one (MIT is common for academic bioinformatics tools)
and add a LICENSE file to the repo root.

### Step 2: Create the package directory

```bash
mkdir -p blast_align_tree/data
```

### Step 3: Move source files into the package

| From (current)           | To (new)                              |
|--------------------------|---------------------------------------|
| `blast_align_tree.py`    | `blast_align_tree/cli.py`             |
| `genome_selector.py`     | `blast_align_tree/genome_selector.py` |
| `visualize_tree.r`       | `blast_align_tree/data/visualize_tree.r` |
| `gene_symbols.txt`       | `blast_align_tree/data/gene_symbols.txt` |

### Step 4: Create `blast_align_tree/__init__.py`

```python
"""blast-align-tree: BLAST, align, and build phylogenetic trees."""
__version__ = "1.0.0"
```

### Step 5: Create `pyproject.toml`

```toml
[build-system]
requires = ["setuptools>=64"]
build-backend = "setuptools.backends._legacy:_Backend"

[project]
name = "blast-align-tree"
version = "1.0.0"
description = "BLAST, align, and build phylogenetic trees from genomic databases"
readme = "README.md"
license = "MIT"
requires-python = ">=3.9"
dependencies = ["biopython"]

[project.scripts]
blast-align-tree = "blast_align_tree.cli:main"
bat-genome-selector = "blast_align_tree.genome_selector:main"

[project.urls]
Homepage = "https://github.com/steinbrennerlab/blast-align-tree"

[tool.setuptools.package-data]
blast_align_tree = ["data/*"]
```

### Step 6: Update file path references in code

These are the code changes needed so bundled files are found after install:

#### 6a. `cli.py` — locate `visualize_tree.r` via package data

Current (`blast_align_tree.py:743`):
```python
rscript = workdir / "visualize_tree.r"
```

Change to:
```python
from importlib.resources import files
_PACKAGE_DATA = files("blast_align_tree") / "data"

# In visualize_tree():
rscript = _PACKAGE_DATA / "visualize_tree.r"
```

This finds the R script inside the installed package instead of assuming it's in
the working directory. The `Rscript` call and `cwd=workdir` stay the same.

#### 6b. `genome_selector.py` — keep `genomes/` relative to cwd

Current (`genome_selector.py:21-22`):
```python
PROJ_DIR = Path(__file__).parent
GENOMES_DIR = PROJ_DIR / "genomes"
```

Change to:
```python
PROJ_DIR = Path.cwd()
GENOMES_DIR = PROJ_DIR / "genomes"
```

After install, `__file__` points into site-packages, not the project directory.
The GUI should scan `genomes/` in the current working directory where the user
runs it. (Or accept a `--genomes-dir` CLI argument for flexibility.)

#### 6c. `cli.py` — same for `genomes/` references

The pipeline already takes genome paths as CLI arguments (`-dbs`, `-qdbs`), so
most references are fine. The `workdir / "genomes"` references (lines 505, 586)
resolve relative to the run directory, which is correct.

### Step 7: Add a genome fetch command

Genome files are too large to ship via pip (PyPI has a ~60 MB per-file soft limit,
and the `genomes/` folder is hundreds of megabytes to gigabytes). Instead, host
them externally and provide a CLI command that downloads them on demand.

#### 7a. Host the genome files

Two good options:

- **GitHub Releases** (simplest): Attach genome tarballs as release assets on
  `steinbrennerlab/blast-align-tree`. Each asset can be up to 2 GB. Free.
- **Zenodo** (best for academia): Upload once, get a permanent DOI, citable in
  papers. Free, unlimited academic storage.

Recommendation: start with GitHub Releases for development; move to Zenodo when
you publish.

#### 7b. Create a manifest file

Add `blast_align_tree/data/genomes_manifest.json` listing each available genome:

```json
{
  "genomes": {
    "TAIR10cds": {
      "url": "https://github.com/steinbrennerlab/blast-align-tree/releases/download/v1.0.0/TAIR10cds.fa.gz",
      "sha256": "abc123...",
      "size_mb": 15,
      "default": true,
      "description": "Arabidopsis thaliana CDS (TAIR10)"
    },
    "Vung469cds": {
      "url": "https://github.com/steinbrennerlab/blast-align-tree/releases/download/v1.0.0/Vung469cds.fa.gz",
      "sha256": "def456...",
      "size_mb": 22,
      "default": true,
      "description": "Vigna unguiculata CDS"
    },
    "Niben261_proteins": {
      "url": "...",
      "sha256": "...",
      "size_mb": 180,
      "default": false,
      "description": "Nicotiana benthamiana proteins"
    }
  }
}
```

The `default: true` flag marks the small essential set; everything else is opt-in.

#### 7c. Add a `fetch-genomes` subcommand

Add a new module `blast_align_tree/fetch.py` with `def main()` exposed as an
entry point. Behavior:

```bash
blast-align-tree-fetch                     # downloads default set to ./genomes/
blast-align-tree-fetch --all               # downloads everything
blast-align-tree-fetch TAIR10cds Vung469cds  # picks specific genomes
blast-align-tree-fetch --list              # shows available genomes + sizes
blast-align-tree-fetch --dest ~/my-genomes # custom destination
```

Implementation outline:
- Read manifest from `importlib.resources.files("blast_align_tree") / "data" / "genomes_manifest.json"`
- For each requested genome: stream-download, verify SHA256, decompress if `.gz`,
  write to `./genomes/<name>.fa`
- Skip files that already exist and pass checksum
- Use only stdlib (`urllib.request`) — no new dependency

#### 7d. Register the entry point

Add to `pyproject.toml`:

```toml
[project.scripts]
blast-align-tree = "blast_align_tree.cli:main"
bat-genome-selector = "blast_align_tree.genome_selector:main"
blast-align-tree-fetch = "blast_align_tree.fetch:main"
```

#### 7e. Document in README

```
After installing, run:
  blast-align-tree-fetch          # downloads small default genome set
  blast-align-tree-fetch --all    # or everything

This populates ./genomes/ in the current directory — run blast-align-tree
from that same directory.
```

### Step 8: Verify it works locally

```bash
pip install -e .
blast-align-tree --help
bat-genome-selector --help
```

The `-e` (editable) install lets you test without rebuilding after each change.

### Step 9: Create a PyPI account and upload

```bash
# One-time: create account at pypi.org, generate API token

# Build
python -m build

# Upload
twine upload dist/*
```

After upload, the package is live at `pypi.org/project/blast-align-tree/`.

---

## What does NOT ship in the pip package

- `genomes/` — too large for PyPI; fetched on demand via `blast-align-tree-fetch`
  (see Step 7). Hosted as GitHub Release assets or on Zenodo.
- `datasets/` — user data
- `hmm_files/` — user HMM profiles
- `environments/` — conda YAML files (still useful for installing external tools)
- `scripts/` — standalone helper scripts (could be added later if desired)
- BLAST indices (`.nhr`, `.nin`, etc.)

---

## What users still need to install separately

All external CLI tools. The README should document:

```
blast-align-tree requires these tools on your PATH:
- BLAST+ (blastdbcmd, tblastn/blastp, makeblastdb)
- MAFFT (>=7)
- trimAl
- RAxML-NG or FastTree
- Rscript + R packages (ggtree, ape, phytools, ggplot2, optparse, treeio, tidytree)
- Optionally: Clustal Omega, PRANK, HMMER

Install via conda:
  conda env create -f environments/bat-environment-linux.yml
```

---

## Future: Bioconda recipe

Once the PyPI package exists, a Bioconda recipe can use `url: pypi` as its source
and add all the non-Python dependencies. This gives users the single-command
`conda install blast-align-tree` experience. The `pyproject.toml` and package
structure from this plan are prerequisites for that.

---

## Risks and Notes

- **`importlib.resources` requires Python 3.9+** for the `files()` API.
  Python 3.8 is EOL, so 3.9+ is reasonable.
- **tkinter**: Ships with most Python installs but not all (e.g. some minimal
  Docker images). Not a pip dependency — document as a requirement.
- **gene_symbols.txt (168 KB)**: Small enough to bundle. If it grows significantly,
  consider a separate download.
- **Existing users**: Anyone cloning the repo and running `python blast_align_tree.py`
  directly can still do so after restructuring, as long as they run from the repo
  root. But `pip install -e .` is the recommended dev workflow after this change.
