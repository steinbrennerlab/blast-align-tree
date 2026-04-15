#!/usr/bin/env python3
"""
Python rewrite of tblastn-align-tree.sh with all helper scripts inlined as functions.

- Orchestrates BLAST+ CLI tools via subprocess
- Replaces awk/xargs/sed with native Python
- **All previous helper scripts (`extract_translate.py`, `translate_db.py`, `pull_id_fasta*.py`, `remove_stop.py`, `add_translations.py`) are now built-in**
- Parallelizes (query, database) jobs
- Requires selected BLAST/alignment/tree tools, Rscript/trimal, Biopython, and R packages; use --check-env to validate

Usage example:
python blast_align_tree.py   -q ENTRY Q2   -qdbs ENTRYDB.fa Q2DB.fa   -n 50 50   -hdr '>' 'gene='   -dbs Vung469.cds.fa TAIR10cds.fa   -add OUTGROUP1   -add_db OUTGROUP_DB.fa   -aa 10 200
"""
from __future__ import annotations
import argparse
import concurrent.futures as cf
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Iterable, List, Tuple, Optional
import re
import tempfile
from datetime import datetime
from importlib.resources import files as _pkg_files

_PACKAGE_DATA = Path(str(_pkg_files("blast_align_tree") / "data"))

# Biopython
try:
    import Bio
    from Bio import SeqIO
    from Bio.Seq import Seq, translate
    BIO_IMPORT_ERROR = None
except ImportError as exc:
    Bio = None
    SeqIO = None
    Seq = None
    translate = None
    BIO_IMPORT_ERROR = exc

# -----------------------
# Utilities
# -----------------------

def check_tool(name: str):
    if shutil.which(name) is None:
        raise SystemExit(f"Required tool not found in PATH: {name}")

def run(cmd: List[str], cwd: Optional[Path] = None, capture: bool = False) -> subprocess.CompletedProcess:
    kw = dict(cwd=str(cwd) if cwd else None, text=True)
    if capture:
        kw["capture_output"] = True
    else:
        kw["stdout"] = subprocess.PIPE
        kw["stderr"] = subprocess.PIPE
    result = subprocess.run(cmd, **kw)
    if result.returncode != 0:
        if not capture and result.stderr:
            sys.stderr.write(result.stderr)
        raise subprocess.CalledProcessError(result.returncode, cmd)
    return result

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def read_lines(fp: Path) -> List[str]:
    if not fp.exists():
        return []
    with fp.open("r", encoding="utf-8", errors="ignore") as f:
        return f.read().splitlines()

def write_text(fp: Path, s: str):
    ensure_dir(fp.parent)
    with fp.open("w", encoding="utf-8") as f:
        f.write(s)

def append_text(fp: Path, s: str):
    ensure_dir(fp.parent)
    with fp.open("a", encoding="utf-8") as f:
        f.write(s)

def merge_files(in_paths: Iterable[Path], out_path: Path):
    ensure_dir(out_path.parent)
    with out_path.open("w", encoding="utf-8") as out:
        for p in in_paths:
            if p.exists():
                with p.open("r", encoding="utf-8", errors="ignore") as f:
                    shutil.copyfileobj(f, out)

def dedup_fasta_by_id(in_path: Path, out_path: Path):
    """Deduplicate FASTA entries by the header token up to first whitespace."""
    ensure_dir(out_path.parent)
    seen = set()
    with in_path.open("r", encoding="utf-8", errors="ignore") as fin, \
         out_path.open("w", encoding="utf-8") as fout:
        write = False
        for line in fin:
            if line.startswith(">"):
                tok = line.strip().split()[0]  # header token
                if tok not in seen:
                    seen.add(tok)
                    write = True
                    fout.write(line)
                else:
                    write = False
            else:
                if write:
                    fout.write(line)

def count_fasta_records(fp: Path) -> int:
    if not fp.exists():
        return 0
    return sum(1 for _ in SeqIO.parse(str(fp), "fasta"))

def prepend_header_line(fp: Path, header: str):
    """Prepend a single header line to a file (like sed -i '1s/^/.../')."""
    if not fp.exists():
        return
    content = fp.read_text(encoding="utf-8", errors="ignore")
    fp.write_text(header + content, encoding="utf-8")

def archive_run(entry_root: Path, timestamp: str) -> Path:
    """
    Move all current contents of ./ENTRY into ./ENTRY/runs/<timestamp>/,
    excluding the 'runs' folder itself. No symlinks are created.
    Returns the final run directory path.
    """
    runs_root = entry_root / "runs"
    run_dir = runs_root / timestamp
    ensure_dir(run_dir)

    for item in list(entry_root.iterdir()):
        if item.name == "runs":
            continue
        target = run_dir / item.name
        ensure_dir(target.parent)
        shutil.move(str(item), str(target))

    print(f"[info] archived to {run_dir}")
    return run_dir

def bt_suffix(blast_type: str) -> str:
    # Used in filenames to distinguish modes
    return "tblastn" if blast_type == "tblastn" else "blastp"

def db_label(db: str) -> str:
    """Flatten a db path (e.g. 'subdir/file.fa') into a safe filename component."""
    return db.replace("/", ".").replace("\\", ".")

def outbase(workdir: Path, entry: str, q: str, db: str, blast_type: str) -> Path:
    return workdir / entry / f"{q}.{db_label(db)}.seq.{bt_suffix(blast_type)}"

def move_old_files(entry_root: Path, timestamp: Optional[str] = None):
    """
    If entry_root contains items other than 'runs' or 'old_files',
    move them into 'old_files/<timestamp>/'.
    - Uses current time if timestamp is not provided.
    - Ensures the target folder is unique if it already exists.
    """
    items = [
        item for item in entry_root.iterdir()
        if item.name not in ("runs", "old_files")
    ]
    if not items:
        return  # nothing to move beyond runs/old_files

    # Build a timestamped folder under old_files/
    if timestamp is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    base_old_dir = entry_root / "old_files" / timestamp
    old_dir = base_old_dir

    # Ensure uniqueness if same-second invocations occur
    suffix = 1
    while old_dir.exists():
        old_dir = entry_root / "old_files" / f"{timestamp}_{suffix}"
        suffix += 1

    ensure_dir(old_dir)

    # Move everything except 'runs' and 'old_files' into the timestamped folder
    for item in items:
        target = old_dir / item.name
        ensure_dir(target.parent)
        shutil.move(str(item), str(target))

    print(f"[info] moved existing files to {old_dir}")


def _move_file_if_exists(src: Path, dst: Path) -> bool:
    if not src.is_file():
        return False
    ensure_dir(dst.parent)
    if dst.exists():
        dst.unlink()
    shutil.move(str(src), str(dst))
    return True


def _unlink_file_if_exists(fp: Path) -> bool:
    if not fp.is_file():
        return False
    fp.unlink()
    return True


def _write_deduped_fasta(in_paths: Iterable[Path], out_path: Path) -> bool:
    parts = [p for p in in_paths if p.is_file()]
    if not parts:
        return False
    ensure_dir(out_path.parent)
    tmp = out_path.with_name(out_path.name + ".tmp")
    merge_files(parts, tmp)
    dedup_fasta_by_id(tmp, out_path)
    _unlink_file_if_exists(tmp)
    return True


def cleanup_run_root(entry_dir: Path, entry: str, queries: List[str], databases: List[str], blast_type: str):
    """
    Keep reports at the run root, move merged/hit sequence material into
    hits/, then remove derivable BLAST/merge staging files before archiving.
    """
    bt = bt_suffix(blast_type)
    hits_dir = entry_dir / "hits"
    ensure_dir(hits_dir)

    moved = 0
    generated = 0
    removed = 0

    for q in queries:
        if _move_file_if_exists(entry_dir / f"{q}.seq.fa", hits_dir / f"{q}.seq.fa"):
            moved += 1

    for db in databases:
        dbl = db_label(db)
        db_dir = hits_dir / dbl
        ensure_dir(db_dir)

        for q in queries:
            q_label = db_label(q)
            base = entry_dir / f"{q}.{dbl}.seq.{bt}"
            if _move_file_if_exists(base, db_dir / f"{q_label}.hit_ids.txt"):
                moved += 1
            full = base.with_suffix(base.suffix + ".full")
            if _move_file_if_exists(full, db_dir / f"{q_label}.blast_report.txt"):
                moved += 1

        aa_out = db_dir / "all_hits.aa.fa"
        db_rmdup = entry_dir / f"{dbl}.parse.merged.rmdup.fa"
        if _move_file_if_exists(db_rmdup, aa_out):
            moved += 1
        else:
            if blast_type == "tblastn":
                aa_parts = sorted(entry_dir.glob(f"*.{dbl}.seq.{bt}.blastdb.translate.fa.parse.fa"))
            else:
                aa_parts = sorted(entry_dir.glob(f"*.{dbl}.seq.{bt}.blastdb.fa.parse.fa"))
            if _write_deduped_fasta(aa_parts, aa_out):
                generated += 1

        if blast_type == "tblastn":
            nt_parts = sorted(entry_dir.glob(f"*.{dbl}.seq.{bt}.blastdb.fa.parse.fa"))
            if _write_deduped_fasta(nt_parts, db_dir / "all_hits.nt.fa"):
                generated += 1

        removed += int(_unlink_file_if_exists(entry_dir / f"{dbl}.parse.merged.fa"))

    if blast_type == "tblastn":
        staging_patterns = [
            f"*.seq.{bt}.blastdb.stop.fa",
            f"*.seq.{bt}.blastdb.fa",
            f"*.seq.{bt}.blastdb.fa.parse.fa",
            f"*.seq.{bt}.blastdb.translate.fa",
            f"*.seq.{bt}.blastdb.translate.fa.parse.fa",
            f"*.seq.{bt}.blastdb.translate.fa.coding.txt",
        ]
        staging_files = [
            entry_dir / f"{entry}.seq.{bt}.blastdb.merged.fa",
            entry_dir / f"{entry}.seq.{bt}.blastdb.nt.parse.merged.fa",
            entry_dir / f"{entry}.seq.{bt}.blastdb.translate.merged.fa",
            entry_dir / f"{entry}.seq.{bt}.blastdb.translate.fa.parse.merged.fa",
            entry_dir / f"{entry}.merged.fa",
            entry_dir / f"{entry}.nt.merged.fa",
        ]
    else:
        staging_patterns = [
            f"*.seq.{bt}.blastdb.fa",
            f"*.seq.{bt}.blastdb.fa.parse.fa",
            f"*.seq.{bt}.blastdb.fa.coding.txt",
        ]
        staging_files = [
            entry_dir / f"{entry}.seq.{bt}.blastdb.merged.fa",
            entry_dir / f"{entry}.seq.{bt}.blastdb.fa.parse.merged.fa",
            entry_dir / f"{entry}.merged.fa",
        ]

    for pattern in staging_patterns:
        for fp in entry_dir.glob(pattern):
            removed += int(_unlink_file_if_exists(fp))

    staging_files.extend([
        Path(str(entry_dir / f"{entry}.parse.merged.aligned.fa") + ".nwk"),
        entry_dir / "mafft.log",
        entry_dir / "merged_coding.txt",
    ])
    for fp in staging_files:
        removed += int(_unlink_file_if_exists(fp))

    for fp in entry_dir.glob(f"{entry}.raxmlng.raxml.*"):
        removed += int(_unlink_file_if_exists(fp))

    final_hit_files = [
        entry_dir / f"{entry}.parse.merged.fa",
        entry_dir / f"{entry}.parse.merged.aligned.fa",
    ]
    if blast_type == "tblastn":
        final_hit_files.append(entry_dir / f"{entry}.nt.parse.merged.fa")
    for fp in final_hit_files:
        if _move_file_if_exists(fp, hits_dir / fp.name):
            moved += 1

    print(f"[cleanup] hits files moved: {moved}; generated: {generated}; staging files removed: {removed}")



# -----------------------
# Inlined helpers (formerly separate scripts)
# -----------------------

def _extract_translate_tblastn(genomes_dir: Path, entry_dir: Path, q: str, qdb: str, aa_slice: Optional[List[int]]):
    """
    Former: extract_translate.py
    - Find nucleotide FASTA record with id==q in genomes/<qdb>
    - Translate to AA (optionally slice aa_start:aa_end)
    - Write to <entry>/<q>.seq.fa
    """
    src = genomes_dir / qdb
    dest = entry_dir / f"{q}.seq.fa"
    ensure_dir(dest.parent)

    aa_start = None
    aa_end = None
    if aa_slice and len(aa_slice) >= 2:
        aa_start, aa_end = aa_slice[0], aa_slice[1]

    found = False
    with open(dest, "w", encoding="utf-8") as out:
        for rec in SeqIO.parse(str(src), "fasta"):
            if rec.id == q:
                aa = str(translate(rec.seq))
                if aa_start is not None:
                    aa = aa[aa_start:aa_end]
                out.write(f">{rec.id}\n{aa}\n")
                found = True
                break

    if not found:
        raise SystemExit(f"Query id '{q}' not found in nucleotide FASTA '{src}' for translation")

def _extract_copy_blastp(genomes_dir: Path, entry_dir: Path, q: str, qdb: str):
    """
    blastp mode: copy protein sequence with header==q from genomes/<qdb> into <entry>/<q>.seq.fa
    """
    src = genomes_dir / qdb
    dest = entry_dir / f"{q}.seq.fa"
    ensure_dir(dest.parent)

    found = False
    with open(src, "r", encoding="utf-8", errors="ignore") as fin, \
         open(dest, "w", encoding="utf-8") as fout:
        write = False
        for line in fin:
            if line.startswith(">"):
                header_tok = line.strip()[1:].split()[0]
                write = (header_tok == q)
                if write:
                    fout.write(f">{q}\n")
                    found = True
            else:
                if write:
                    fout.write(line)
    if not found:
        raise SystemExit(f"Query id '{q}' not found in protein FASTA '{src}'")

def _remove_stop_codons(in_fa: Path, out_fa: Path):
    """
    Former: remove_stop.py
    Remove TAG/TGA/TAA anywhere in-frame across the sequence.
    """
    stops = {"TAG", "TGA", "TAA"}
    my_records = []
    for record in SeqIO.parse(str(in_fa), "fasta"):
        seq_list = list(str(record.seq))
        i = 0
        while i + 2 < len(seq_list):
            codon = "".join(seq_list[i:i+3]).upper()
            if codon in stops:
                del seq_list[i:i+3]
                # do not advance i; next codon now at same index
            else:
                i += 3
        record.seq = Seq("".join(seq_list))
        my_records.append(record)
    ensure_dir(out_fa.parent)
    SeqIO.write(my_records, str(out_fa), "fasta")

def _translate_fasta(in_fa: Path, out_fa: Path):
    """
    Former: translate_db.py
    For each nucleotide record in in_fa, write translated AA with the original description as header.
    """
    ensure_dir(out_fa.parent)
    with open(out_fa, "w", encoding="utf-8") as out:
        for rec in SeqIO.parse(str(in_fa), "fasta"):
            seq = rec.seq

            #enforce full codons AFTER stop removal
            seq = seq[:len(seq) - (len(seq) % 3)]

            aa = translate(seq)  # keep semantics unchanged

            out.write(">" + rec.description + "\n")
            out.write(str(aa) + "\n")

def _parse_header_token(description: str, headerword: str, fallback_id: str, suffix: str = "") -> str:
    """
    Former: pull_id_fasta*.py logic
    If headerword == 'id' → return fallback_id (record.id).
    Else: find substring after the first occurrence of headerword and take until next space.
    If headerword not found, return fallback_id.
    If suffix is provided, strip it from the end of the token.
    """
    if headerword == "id":
        token = fallback_id
    elif headerword not in description:
        token = fallback_id
    else:
        # Split on the headerword and take the part immediately following, up to first space
        try:
            part = description.split(headerword, 1)[1]
            token = part.split(" ", 1)[0]
        except Exception:
            token = fallback_id
    if suffix and token.endswith(suffix):
        token = token[:-len(suffix)]
    return token

def _parse_fasta_headers(in_fa: Path, out_fa: Path, headerword: str, suffix: str = ""):
    """
    Write a new FASTA where each header is the parsed token.
    """
    ensure_dir(out_fa.parent)
    with open(out_fa, "w", encoding="utf-8") as out:
        for rec in SeqIO.parse(str(in_fa), "fasta"):
            token = _parse_header_token(rec.description, headerword, rec.id, suffix)
            out.write(f">{token}\n{str(rec.seq)}\n")

def _coding_table(in_fa: Path, out_txt: Path, headerword: str, db_name: str, suffix: str = ""):
    """
    Create <parsed_token>\t<db_name> lines for each record.
    """
    ensure_dir(out_txt.parent)
    with open(out_txt, "a", encoding="utf-8") as out:
        for rec in SeqIO.parse(str(in_fa), "fasta"):
            token = _parse_header_token(rec.description, headerword, rec.id, suffix)
            out.write(f"{token}\t{db_name}\n")

def _add_translation_from_db(genomes_dir: Path, entry_dir: Path, db: str, seq_id: str):
    """
    Former: add_translations.py
    Append translation of seq_id from genomes/<db> to:
      - <entry>/<entry>.parse.merged.fa
      - <entry>/merged_coding.txt
    """
    src = genomes_dir / db
    out_fa = entry_dir / f"{entry_dir.name}.parse.merged.fa"
    out_txt = entry_dir / "merged_coding.txt"
    ensure_dir(out_fa.parent)
    ensure_dir(out_txt.parent)

    for rec in SeqIO.parse(str(src), "fasta"):
        if rec.id == seq_id:
            with open(out_fa, "a", encoding="utf-8") as fa:
                fa.write(f">{rec.id}\n{str(translate(rec.seq))}\n")
            with open(out_txt, "a", encoding="utf-8") as txt:
                txt.write(f"\n{rec.id}\t{db}")
            print(f"[add] added {seq_id} from {db}")
            return
    raise SystemExit(f"Sequence id '{seq_id}' not found in {src}")

# -----------------------
# Core pipeline steps
# -----------------------

def extract_and_translate(
    entry: str,
    q: str,
    qdb: str,
    slice_args: Optional[List[str]],
    workdir: Path,
    blast_type: str
):
    entry_dir = workdir / entry
    genomes_dir = workdir / "genomes"
    aa_slice = None
    if slice_args and len(slice_args) >= 2:
        # Expecting AA start end
        aa_slice = [int(slice_args[0]), int(slice_args[1])]

    if blast_type == "tblastn":
        _extract_translate_tblastn(genomes_dir, entry_dir, q, qdb, aa_slice)
    else:
        _extract_copy_blastp(genomes_dir, entry_dir, q, qdb)

def blast_and_post(entry: str, q: str, db: str, max_targets: str, workdir: Path, blast_type: str) -> Tuple[str, str, str]:
    """
    Run tblastn or blastp for a (query, db) pair and return:
      - path to sseqid list (tab list)
      - path to full report
      - path to fetched sequences FASTA (nt for tblastn; aa for blastp)
    """
    q_fa = workdir / entry / f"{q}.seq.fa"
    db_path = workdir / "genomes" / db
    out_base = outbase(workdir, entry, q, db, blast_type)

    # IDs only
    run([
        "tblastn" if blast_type == "tblastn" else "blastp",
        "-query", str(q_fa), "-db", str(db_path),
        "-max_target_seqs", str(max_targets), "-max_hsps", "1",
        "-outfmt", "6 sseqid", "-out", str(out_base)
    ])
    # full report (pairwise)
    full = out_base.with_suffix(out_base.suffix + ".full")
    run([
        "tblastn" if blast_type == "tblastn" else "blastp",
        "-query", str(q_fa), "-db", str(db_path),
        "-max_hsps", "1", "-outfmt", "1", "-out", str(full)
    ])
    prepend_header_line(full, "hit query_id	subject_id	pct_identity	aln_length	n_of_mismatches	gap_openings	q_start q_end	s_start   s_end	e_value bit_score\n")

    # Fetch sequences from the BLAST DB for the hit list
    if blast_type == "tblastn":
        stop_fa = Path(str(out_base) + ".blastdb.stop.fa")
        nt_fa   = Path(str(out_base) + ".blastdb.fa")
        run(["blastdbcmd", "-db", str(db_path), "-entry_batch", str(out_base), "-out", str(stop_fa)])
        _remove_stop_codons(stop_fa, nt_fa)
        return str(out_base), str(full), str(nt_fa)
    else:
        prot_fa = Path(str(out_base) + ".blastdb.fa")
        run(["blastdbcmd", "-db", str(db_path), "-entry_batch", str(out_base), "-out", str(prot_fa)])
        return str(out_base), str(full), str(prot_fa)

def translate_and_parse_headers(
    entry: str,
    q: str,
    db: str,
    header_rule: str,
    workdir: Path,
    blast_type: str,
    header_suffix: str = ""):
    bt = bt_suffix(blast_type)
    entry_dir = workdir / entry
    dbl = db_label(db)

    if blast_type == "tblastn":
        # translate_db produces *.seq.tblastn.blastdb.translate.fa
        in_nt = entry_dir / f"{q}.{dbl}.seq.{bt}.blastdb.fa"
        out_aa = entry_dir / f"{q}.{dbl}.seq.{bt}.blastdb.translate.fa"
        _translate_fasta(in_nt, out_aa)

        # pull_id_fasta for nt and translated
        _parse_fasta_headers(in_nt, entry_dir / f"{q}.{dbl}.seq.{bt}.blastdb.fa.parse.fa", header_rule, header_suffix)
        _parse_fasta_headers(out_aa, entry_dir / f"{q}.{dbl}.seq.{bt}.blastdb.translate.fa.parse.fa", header_rule, header_suffix)

        # coding table on translated AA
        _coding_table(out_aa, entry_dir / f"{q}.{dbl}.seq.{bt}.blastdb.translate.fa.coding.txt", header_rule, db, header_suffix)
    else:
        # blastp: only AA present → only pull/parse once from *.blastdb.fa
        in_aa = entry_dir / f"{q}.{dbl}.seq.{bt}.blastdb.fa"
        _parse_fasta_headers(in_aa, entry_dir / f"{q}.{dbl}.seq.{bt}.blastdb.fa.parse.fa", header_rule, header_suffix)
        _coding_table(in_aa, entry_dir / f"{q}.{dbl}.seq.{bt}.blastdb.fa.coding.txt", header_rule, db, header_suffix)

def optional_add_translations(entry: str, add_dbs: List[str], add_seqs: List[str], workdir: Path):
    genomes_dir = workdir / "genomes"
    entry_dir = workdir / entry
    for db, seq in zip(add_dbs, add_seqs):
        _add_translation_from_db(genomes_dir, entry_dir, db, seq)

def align_and_build_tree(entry: str, workdir: Path, aligner: str, tree_builder: str, threads: int, mafft_mode: str):
    entry_dir = workdir / entry
    in_fa = entry_dir / f"{entry}.parse.merged.fa"
    seq_count = count_fasta_records(in_fa)

    if seq_count < 2:
        raise SystemExit(
            f"Need at least 2 sequences to align after merge/dedup, but found {seq_count} in {in_fa}. "
            "This usually means BLAST returned only one unique hit for the query."
        )

    # Canonical alignment filename (independent of aligner)
    aln_fa = entry_dir / f"{entry}.parse.merged.aligned.fa"

    print(f"→ Aligning sequences with {aligner} ...")
    align_start = datetime.now()


    if aligner == "clustalo":
        # Clustal Omega writes directly to output file
        run([
            "clustalo",
            "-i", str(in_fa),
            "-o", str(aln_fa),
            "--force"
        ])

    elif aligner == "mafft":
        cmd = ["mafft", "--thread", str(threads)]

        if mafft_mode == "auto":
            cmd.append("--auto")
        elif mafft_mode == "linsi":
            cmd += ["--localpair", "--maxiterate", "1000"]
        elif mafft_mode == "einsi":
            cmd += ["--genafpair", "--maxiterate", "1000"]
        elif mafft_mode == "fftns2":
            cmd += ["--retree", "2"]

        cmd.append(str(in_fa))

        result = subprocess.run(
            cmd,
            text=True,
            capture_output=True
        )

        mafft_log = entry_dir / "mafft.log"
        mafft_log.write_text(result.stderr, encoding="utf-8")

        if result.returncode != 0:
            sys.stderr.write(result.stderr)
            raise subprocess.CalledProcessError(result.returncode, "mafft")

        strategy = "unknown"
        lines = result.stderr.splitlines()

        for i, line in enumerate(lines):
            line = line.strip()

            # Case 1: single-line strategy
            if line.startswith("Strategy") and len(line.split()) > 1:
                strategy = line
                break

            if line.startswith("Using"):
                strategy = line
                break

            # Case 2: "Strategy:" followed by next line
            if line == "Strategy:" and i + 1 < len(lines):
                strategy = f"Strategy {lines[i + 1].strip()}"
                break

        print(f"  MAFFT strategy: {strategy}")


        # Keep only FASTA output
        fasta_lines = [
            line for line in result.stdout.splitlines()
            if line.startswith(">") or re.match(r"^[A-Za-z\-\*]+$", line)
        ]

        aln_fa.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")


    align_end = datetime.now()
    print(f"  Alignment time: {align_end - align_start}")
    print(f"→ Building tree with {tree_builder} ...")
    tree_start = datetime.now()
    tree_out = Path(str(aln_fa) + ".nwk")

    if tree_builder == "FastTree":
        run(["FastTree", "-out", str(tree_out), str(aln_fa)])

    elif tree_builder == "RAxML":
        prefix = f"{entry}.raxmlng"
        prefix_path = entry_dir / prefix

        # ---- Step 1: ML tree search ----
        print("  [RAxML] Running ML tree search...")
        run([
            "raxml-ng",
            "--msa", str(aln_fa),
            "--model", "LG+G",
            "--threads", str(threads),
            "--blopt", "nr_safe",
            "--prefix", prefix
        ], cwd=entry_dir)

        best_tree = entry_dir / f"{prefix}.raxml.bestTree"
        if not best_tree.exists():
            raise SystemExit(f"RAxML-NG did not produce a bestTree: {best_tree}")

        # ---- Step 2: Bootstrap replicates ----
        print("  [RAxML] Running bootstrap replicates...")
        run([
            "raxml-ng",
            "--bootstrap",
            "--msa", str(aln_fa),
            "--model", "LG+G",
            "--threads", str(threads),
            "--bs-trees", "100",   # <-- adjust as desired
            "--prefix", prefix
        ], cwd=entry_dir)

        bs_trees = entry_dir / f"{prefix}.raxml.bootstraps"
        if not bs_trees.exists():
            raise SystemExit(f"RAxML-NG did not produce bootstrap trees: {bs_trees}")

        # ---- Step 3: Map bootstrap support onto ML tree ----
        print("  [RAxML] Mapping bootstrap support...")
        run([
            "raxml-ng",
            "--support",
            "--tree", str(best_tree),
            "--bs-trees", str(bs_trees),
            "--prefix", prefix
        ], cwd=entry_dir)

        support_tree = entry_dir / f"{prefix}.raxml.support"
        if not support_tree.exists():
            raise SystemExit(f"RAxML-NG did not produce support tree: {support_tree}")

        # Final output: bootstrap-supported tree (like FastTree output)
        shutil.copyfile(support_tree, tree_out)

def visualize_tree(entry: str, queries: List[str], workdir: Path, datasets: Optional[str] = None):
    """
    Run visualize-tree.r on the combinedtree.nwk
    By default, --write argument is set to the first query name.
    """
    rscript = _PACKAGE_DATA / "visualize_tree.r"

    if not rscript.exists():
        raise SystemExit(f"R script not found: {rscript}")

    write_arg = queries[0]  # default = first -q
    cmd = [
        "Rscript", str(rscript),
        "--entry", entry,
        "--write", write_arg
    ]
    if datasets:
        cmd.extend(["--datasets", datasets])
    run(cmd, cwd=workdir)

# -----------------------
# Motif & HMMER feature detection
# -----------------------

def _prosite_to_regex(pat: str) -> str:
    """
    Minimal PROSITE→regex converter:
      - x or X → .
      - - is removed
      - {P}  → [^P]
      - [ST] stays [ST]
      - x(2) or  → .{2} / [ST]{2}
      - x(2,4) → .{2,4}
      - '<' anchor start, '>' anchor end (optional)
    Examples:
      C-x(2)-C       -> C.{2}C
      H-x(2)-H       -> H.{2}H
      G-[ST]-x(2)-E  -> G[ST].{2}E
              -> [DE]{2}
      <M-x(3)-K>     -> ^M.{3}K$
    """
    s = pat.strip()
    # anchors
    anchor_start = s.startswith("<")
    anchor_end   = s.endswith(">")
    s = s.lstrip("<").rstrip(">")

    # basic substitutions
    s = s.replace("-", "")
    s = re.sub(r"[xX]", ".", s)
    s = re.sub(r"\{([^}]+)\}", r"[^\1]", s)  # {P} -> [^P]

    # (n) and (m,n) repeat counts on the previous token
    # previous token can be: literal AA, '.', bracket class [...]
    def _rep_sub(m):
        token = m.group(1)
        count = m.group(2)
        return f"{token}{{{count}}}"

    def _rep_rng_sub(m):
        token = m.group(1)
        a, b = m.group(2), m.group(3)
        return f"{token}{{{a},{b}}}"

    s = re.sub(r"(\.|\[[^\]]+\]|[A-Za-z])\((\d+)\)", _rep_sub, s)
    s = re.sub(r"(\.|\[[^\]]+\]|[A-Za-z])\((\d+),(\d+)\)", _rep_rng_sub, s)

    if anchor_start:
        s = "^" + s
    if anchor_end:
        s = s + "$"
    return s

def _parse_named_motifs(raw: list[str], syntax: str) -> list[tuple[str, re.Pattern]]:
    """
    raw items can be:
      NAME=PATTERN   (feature name + pattern)
      PATTERN        (feature name becomes the original string)
    syntax: "regex" or "prosite"
    """
    motifs = []
    for item in raw:
        if "=" in item:
            name, pat = item.split("=", 1)
        else:
            name, pat = item, item
        if syntax == "prosite":
            pat = _prosite_to_regex(pat)
        compiled = re.compile(pat)
        motifs.append((name, compiled))
    return motifs

def _build_alignment_maps(aln_fa: Path) -> dict[str, list[int]]:
    """
    For each aligned sequence, build a map:
      map[ungapped_aa_index] = aligned_column_index   (1-based, inclusive)
    Useful to translate unaligned hits -> alignment coordinates.
    """
    maps = {}
    for rec in SeqIO.parse(str(aln_fa), "fasta"):
        seq = str(rec.seq)
        ungapped = 0
        m = {}
        for i, ch in enumerate(seq, start=1):
            if ch != "-":
                ungapped += 1
                m[ungapped] = i
        maps[rec.id] = m
    return maps

def _scan_motifs(seq_fa: Path, motifs: list[tuple[str, re.Pattern]], allow_overlap: bool) -> list[tuple[str, int, int, str]]:
    """
    Returns list of (label, start_aa_unaligned, end_aa_unaligned, feature)
    """
    hits = []
    for rec in SeqIO.parse(str(seq_fa), "fasta"):
        s = str(rec.seq)
        for feat_name, pat in motifs:
            if allow_overlap:
                # Use lookahead to capture overlaps
                la = re.compile(f"(?=({pat.pattern}))")
                for m in la.finditer(s):
                    span = m.span(1)
                    if span[0] == -1:
                        continue
                    start = span[0] + 1
                    end   = span[1]
                    hits.append((rec.id, start, end, feat_name))
            else:
                for m in pat.finditer(s):
                    start = m.start() + 1
                    end   = m.end()
                    hits.append((rec.id, start, end, feat_name))
    return hits

def _run_hmmscan(hmm_file: Path, seq_fa: Path, tmpdir: Path) -> list[tuple[str, int, int, str]]:
    """
    Run hmmscan --domtblout on seq_fa against hmm_file.
    Returns list of (label, start_aa_unaligned, end_aa_unaligned, feature=HMM name).
    """
    domtbl = tmpdir / (hmm_file.stem + ".domtblout")
    run(["hmmscan", "--noali", "--cut_ga", "--domtblout", str(domtbl), str(hmm_file), str(seq_fa)])
    out = []
    if not domtbl.exists():
        return out
    with domtbl.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 23:
                continue
            hmm_name   = parts[0]    # target (HMM) name
            query_name = parts[3]    # sequence id
            ali_from   = int(parts[17])
            ali_to     = int(parts[18])
            start, end = (ali_from, ali_to) if ali_from <= ali_to else (ali_to, ali_from)
            out.append((query_name, start, end, f"HMM:{hmm_name}"))
    return out

def _to_alignment_coords(hits_unaligned: list[tuple[str,int,int,str]],
                         aln_maps: dict[str, list[int]]) -> list[tuple[str,int,int,str]]:
    aligned = []
    for label, s, e, feat in hits_unaligned:
        m = aln_maps.get(label)
        if not m:
            # label not present in alignment; skip
            continue
        if s in m and e in m:
            aligned.append((label, m[s], m[e], feat))
        else:
            # If edges fall off (rare), try best-effort clamp inward
            # Find nearest mapped positions
            s2 = next((m[i] for i in range(s, 0, -1) if i in m), None)
            e2 = next((m[i] for i in range(e, 0, -1) if i in m), None)
            if s2 is not None and e2 is not None and s2 <= e2:
                aligned.append((label, s2, e2, feat))
    return aligned

def _write_features_tsv(rows: list[tuple[str,int,int,str]], out_fp: Path):
    ensure_dir(out_fp.parent)
    rows_sorted = sorted(rows, key=lambda r: (r[0], r[1], r[2], r[3]))
    with out_fp.open("w", encoding="utf-8") as w:
        w.write("label\taa_start\taa_end\tfeature\n")
        for label, a, b, feat in rows_sorted:
            w.write(f"{label}\t{a}\t{b}\t{feat}\n")

def annotate_features(entry: str,
                      workdir: Path,
                      motifs_raw: list[str],
                      motif_syntax: str,
                      motif_overlap: bool,
                      hmm_files: list[str]) -> Path:
    """
    Generate features.txt with **raw (unaligned)** AA coordinates
    taken directly from the unaligned per-sequence FASTA.
    """
    entry_dir = workdir / entry
    seq_fa = entry_dir / f"{entry}.parse.merged.fa"   # unaligned AA with tree tip labels
    out_fp = entry_dir / "features.txt"

    all_hits_unaligned: list[tuple[str,int,int,str]] = []

    # Motifs
    if motifs_raw:
        motifs = _parse_named_motifs(motifs_raw, motif_syntax)
        all_hits_unaligned.extend(
            _scan_motifs(seq_fa, motifs, allow_overlap=motif_overlap)
        )

    # HMMER
    if hmm_files:
        check_tool("hmmscan")
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            for hmm in hmm_files:
                all_hits_unaligned.extend(_run_hmmscan(Path(hmm), seq_fa, tdp))

    # Write **unaligned** coordinates directly
    _write_features_tsv(all_hits_unaligned, out_fp)
    return out_fp


# -----------------------
# Environment checks
# -----------------------

REQUIRED_R_PACKAGES = [
    "ggplot2",
    "treeio",
    "phytools",
    "ggtree",
    "optparse",
    "tidytree",
    "ape",
    "broom",
    "Biostrings",
]


def selected_external_tools(args: argparse.Namespace) -> list[str]:
    tools = [
        "blastdbcmd",
        "Rscript",
        "python",
        "trimal",
        "tblastn" if args.blast_type == "tblastn" else "blastp",
    ]

    if args.aligner == "clustalo":
        tools.append("clustalo")
    elif args.aligner == "mafft":
        tools.append("mafft")

    if args.tree_builder == "FastTree":
        tools.append("FastTree")
    elif args.tree_builder == "RAxML":
        tools.append("raxml-ng")

    if args.hmms:
        tools.append("hmmscan")

    return list(dict.fromkeys(tools))


def _one_line(text: str) -> str:
    return " ".join(text.strip().split())


def active_env_roots() -> list[Path]:
    roots = [Path(sys.prefix).resolve()]
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        roots.append(Path(conda_prefix).resolve())
    return list(dict.fromkeys(roots))


def _is_under_any_root(path: str, roots: list[Path]) -> bool:
    exe = Path(path).resolve()
    for root in roots:
        try:
            exe.relative_to(root)
            return True
        except ValueError:
            continue
    return False


def _format_roots(roots: list[Path]) -> str:
    return "; ".join(str(root) for root in roots)


def _check_path_python_biopython() -> tuple[bool, str]:
    if shutil.which("python") is None:
        return False, "python not found on PATH"

    code = (
        "import sys, Bio; "
        "print(getattr(Bio, '__version__', 'unknown')); "
        "print(sys.executable)"
    )
    result = subprocess.run(["python", "-c", code], text=True, capture_output=True)
    if result.returncode == 0:
        lines = result.stdout.splitlines()
        version = lines[0] if lines else "unknown"
        executable = lines[1] if len(lines) > 1 else "python"
        return True, f"{version} ({executable})"

    detail = _one_line(result.stderr or result.stdout)
    return False, detail or "python could not import Bio"


def _check_r_packages() -> tuple[bool, str]:
    if shutil.which("Rscript") is None:
        return False, "Rscript not found on PATH"

    package_vector = ", ".join(f'"{pkg}"' for pkg in REQUIRED_R_PACKAGES)
    expr = (
        f"pkgs <- c({package_vector}); "
        "missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; "
        "if (length(missing)) { cat(paste(missing, collapse = '\\n')); quit(status = 1) }; "
        "cat(paste(pkgs, collapse = ', '))"
    )
    result = subprocess.run(["Rscript", "-e", expr], text=True, capture_output=True)
    if result.returncode == 0:
        return True, _one_line(result.stdout)

    missing = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if missing:
        return False, "missing: " + ", ".join(missing)

    detail = _one_line(result.stderr or result.stdout)
    return False, detail or "R package check failed"


def collect_environment_checks(args: argparse.Namespace) -> tuple[list[tuple[str, str, str]], list[str], list[str]]:
    checks: list[tuple[str, str, str]] = []
    missing: list[str] = []
    warnings: list[str] = []
    roots = active_env_roots()

    def record(status: str, label: str, detail: str):
        checks.append((status, label, detail))
        if status == "MISSING":
            missing.append(f"{label}: {detail}")
        elif status == "WARN":
            warnings.append(f"{label}: {detail}")

    record("OK", "active env root", _format_roots(roots))
    record("OK", "current Python", sys.executable)
    if BIO_IMPORT_ERROR is None:
        record("OK", "Biopython (current Python)", getattr(Bio, "__version__", "unknown"))
    else:
        record("MISSING", "Biopython (current Python)", str(BIO_IMPORT_ERROR))

    for tool in selected_external_tools(args):
        path = shutil.which(tool)
        if path:
            if _is_under_any_root(path, roots):
                record("OK", f"tool {tool}", path)
            else:
                record(
                    "WARN",
                    f"tool {tool}",
                    f"{path} (outside active env: {_format_roots(roots)})",
                )
        else:
            record("MISSING", f"tool {tool}", "not found on PATH")

    ok, detail = _check_path_python_biopython()
    record("OK" if ok else "MISSING", "Biopython (PATH python)", detail)

    ok, detail = _check_r_packages()
    record("OK" if ok else "MISSING", "R packages", detail)

    return checks, missing, warnings


def print_environment_report(args: argparse.Namespace) -> int:
    checks, missing, warnings = collect_environment_checks(args)

    print("Environment check")
    print("=================")
    for status, label, detail in checks:
        print(f"{status:<7} {label}: {detail}")

    strict_failures = warnings if args.strict_env else []
    if missing or strict_failures:
        print("\nEnvironment check failed:")
        for item in missing:
            print(f"  - {item}")
        for item in strict_failures:
            print(f"  - strict-env violation: {item}")
        return 1

    if warnings:
        print("\nEnvironment check passed with warnings:")
        for item in warnings:
            print(f"  - {item}")
        print("Use --strict-env to treat these warnings as failures.")
        return 0

    print("\nEnvironment check passed.")
    return 0


def ensure_environment(args: argparse.Namespace):
    _, missing, warnings = collect_environment_checks(args)
    strict_failures = warnings if args.strict_env else []
    if missing or strict_failures:
        items = missing + [f"strict-env violation: {item}" for item in strict_failures]
        lines = "\n".join(f"  - {item}" for item in items)
        raise SystemExit(
            "Environment check failed before running the pipeline:\n"
            f"{lines}\n"
            "Run with --check-env for the full report."
        )


def validate_pipeline_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> list[str]:
    required = [
        ("queries", "-q/--queries"),
        ("query_databases", "-qdbs/--query_databases"),
        ("seqs", "-n/--seqs"),
        ("header", "-hdr/--header"),
        ("database", "-dbs/--database"),
    ]
    missing = [flag for attr, flag in required if getattr(args, attr) is None]
    if missing:
        parser.error(
            "the following arguments are required for a pipeline run: "
            + ", ".join(missing)
        )

    if len(args.query_databases) != len(args.queries):
        parser.error("queries and query_databases lengths must match")
    if len(args.seqs) != len(args.database):
        parser.error("seqs (-n) must have same length as database (-dbs)")
    if len(args.header) != len(args.database):
        parser.error("header (-hdr) must have same length as database (-dbs)")
    if args.add_seqs and (len(args.add_seqs) != len(args.add_dbs)):
        parser.error("add_seqs and add_dbs must have the same length")

    if args.header_suffix is not None:
        if len(args.header_suffix) != len(args.database):
            parser.error("header_suffix (-hdr_sfx) must have same length as database (-dbs)")
        return [("" if s.lower() == "none" else s) for s in args.header_suffix]

    return [""] * len(args.database)


# -----------------------
# Main
# -----------------------

def main():
    ap = argparse.ArgumentParser(description="Python rewrite of tblastn-align-tree.sh (single-file, helpers inlined)")

    ap.add_argument("--aligner", choices=["clustalo", "mafft"], default="clustalo", help="Multiple sequence aligner to use (default=clustalo)")
    ap.add_argument("--mafft_mode", choices=["auto", "linsi", "einsi", "fftns2"], default="auto", help="MAFFT alignment strategy (default: auto)")
    ap.add_argument("--tree_builder", choices=["FastTree", "RAxML"], default="FastTree", help="Tree builder to use (default=FastTree)")
    ap.add_argument("--check-env", action="store_true", help="check selected Python, R, and external CLI dependencies, then exit")
    ap.add_argument("--strict-env", action="store_true", help="treat tools found outside the active environment as failures")
    ap.add_argument("-q", "--queries", nargs="+", help="fasta ID for query sequences")
    ap.add_argument("-qdbs", "--query_databases", nargs="+", help="fasta files containing the queries")
    ap.add_argument("-n", "--seqs", nargs="+", help="-max_target_seqs per search database (align with -dbs)")
    ap.add_argument("-hdr", "--header", nargs="+", help="header parsing rules per db (align with -dbs)")
    ap.add_argument("-hdr_sfx", "--header_suffix", nargs="+", default=None,
                    help="suffix to trim from parsed header tokens, per db (use 'none' for no trimming)")
    ap.add_argument("-dbs", "--database", nargs="+", help="blast databases to search (filenames or subfolder paths under ./genomes)")
    ap.add_argument("-add", "--add_seqs", nargs="*", default=[], help="additional sequences (optional)")
    ap.add_argument("-add_db", "--add_dbs", nargs="*", default=[], help="databases for additional sequences (optional)")
    ap.add_argument("-aa", "--slice", nargs="*", default=[], help="AA slice start end, optional")
    ap.add_argument("--threads", type=int, default=max(1, os.cpu_count() // 2), help="parallel jobs for BLAST steps")
    ap.add_argument("--workdir", default=".", help="working directory (default=.)")
    ap.add_argument("--datasets", default=None,
                    help="Comma- or semicolon-separated dataset files/directories to pass to visualize_tree.r")
    ap.add_argument(
        "--blast_type",
        choices=["tblastn", "blastp"],
        default="tblastn",
        help="Choose between tblastn (DNA DBs) or blastp (protein DBs). Default=tblastn"
    )
    ap.add_argument("--motif", dest="motifs", nargs="*", default=[],
                    help="AA motif patterns. Accepts NAME=PATTERN to name a motif. Default syntax is regex.")
    ap.add_argument("--motif_syntax", choices=["regex", "prosite"], default="regex",
                    help="Syntax for --motif patterns. Use 'prosite' for simple PROSITE-like forms (e.g., C-x(2)-C).")
    ap.add_argument("--motif_overlap", action="store_true",
                    help="Allow overlapping motif matches (uses regex lookahead).")
    ap.add_argument("--hmm", dest="hmms", nargs="*", default=[],
                    help="One or more HMMER profile files (.hmm). Scans unaligned AA sequences with hmmscan.")
    args = ap.parse_args()
    blast_type = args.blast_type

    if args.check_env:
        raise SystemExit(print_environment_report(args))

    header_suffixes = validate_pipeline_args(args, ap)
    ensure_environment(args)

    workdir = Path(args.workdir).resolve()
    entry = args.queries[0]
    entrydb = args.query_databases[0]  # kept for parity with bash; extract uses both

    # Create directories
    entry_dir = workdir / entry

    ensure_dir(entry_dir)

    # Check if directory has leftover files besides 'runs'
    move_old_files(entry_dir)

    print(f"\n{'='*50}")
    print(f"  blast-align-tree")
    print(f"{'='*50}")
    print(f"  Entry:      {entry}")
    print(f"  Queries:    {', '.join(args.queries)}")
    print(f"  Databases:  {', '.join(args.database)}")
    print(f"  Blast type: {blast_type}")
    print(f"  Max seqs:   {', '.join(args.seqs)}")
    print(f"{'='*50}")

    # Step 1: extract and translate each query (optionally AA slice)
    print(f"\n→ Extracting query sequences")
    for q, qdb in zip(args.queries, args.query_databases):
        extract_and_translate(entry, q, qdb, args.slice if args.slice else None, workdir, blast_type)

    # Step 2: parallel BLAST per (query, db)
    print(f"  Parallel BLAST Jobs {args.threads}")
    print(f"\n→ Running BLAST searches")
    jobs: List[Tuple[str, str, str, str, str]] = []  # (q, db, n, header_rule, header_suffix)
    for q in args.queries:
        for db, n, hdr, sfx in zip(args.database, args.seqs, args.header, header_suffixes):
            jobs.append((q, db, n, hdr, sfx))

    def _one(job: Tuple[str, str, str, str, str]):
        q, db, n, hdr, sfx = job
        out_base, full, fetched_fa = blast_and_post(entry, q, db, n, workdir, blast_type)
        translate_and_parse_headers(entry, q, db, hdr, workdir, blast_type, sfx)
        return job

    blast_start = datetime.now()

    with cf.ThreadPoolExecutor(max_workers=args.threads) as ex:
        list(ex.map(_one, jobs))

    blast_end = datetime.now()
    blast_dt = blast_end - blast_start
    print(f"  BLAST time: {blast_dt}")

    print(f"\n→ Merging results")
    # Step 3: combine coding txt → merged_genome_mapping.txt and prepend header
    coding_txts = sorted(entry_dir.glob("*.coding.txt"))
    merged_genome_mapping = entry_dir / "merged_genome_mapping.txt"
    merge_files(coding_txts, merged_genome_mapping)
    prepend_header_line(merged_genome_mapping, "taxa\tgenome\n")

    # Step 4: merge and dedup FASTAs (conditional on blast_type)
    bt = bt_suffix(blast_type)

    if blast_type == "tblastn":
        # NT FASTAs (only for tblastn)
        all_nt = sorted(entry_dir.glob(f"*.seq.{bt}.blastdb.fa"))
        nt_merged = entry_dir / f"{entry}.seq.{bt}.blastdb.merged.fa"
        merge_files(all_nt, nt_merged)
        nt_dedup = entry_dir / f"{entry}.nt.merged.fa"
        dedup_fasta_by_id(nt_merged, nt_dedup)

        all_nt_parse_sources = sorted(entry_dir.glob(f"*.seq.{bt}.blastdb.fa.parse.fa"))
        nt_parse_merged_src = entry_dir / f"{entry}.seq.{bt}.blastdb.nt.parse.merged.fa"
        merge_files(all_nt_parse_sources, nt_parse_merged_src)
        nt_parse_dedup = entry_dir / f"{entry}.nt.parse.merged.fa"
        dedup_fasta_by_id(nt_parse_merged_src, nt_parse_dedup)

        # Translated AA from tblastn
        all_trans = sorted(entry_dir.glob(f"*.seq.{bt}.blastdb.translate.fa"))
        trans_merged = entry_dir / f"{entry}.seq.{bt}.blastdb.translate.merged.fa"
        merge_files(all_trans, trans_merged)
        trans_dedup = entry_dir / f"{entry}.merged.fa"
        dedup_fasta_by_id(trans_merged, trans_dedup)

        all_trans_parse = sorted(entry_dir.glob(f"*.seq.{bt}.blastdb.translate.fa.parse.fa"))
        trans_parse_merged = entry_dir / f"{entry}.seq.{bt}.blastdb.translate.fa.parse.merged.fa"
        merge_files(all_trans_parse, trans_parse_merged)
        parse_merged = entry_dir / f"{entry}.parse.merged.fa"
        dedup_fasta_by_id(trans_parse_merged, parse_merged)

    else:
        # blastp: AA are directly in *.seq.blastp.blastdb.fa → parse creates *.parse.fa
        all_aa = sorted(entry_dir.glob(f"*.seq.{bt}.blastdb.fa"))
        aa_merged = entry_dir / f"{entry}.seq.{bt}.blastdb.merged.fa"
        merge_files(all_aa, aa_merged)
        aa_dedup = entry_dir / f"{entry}.merged.fa"
        dedup_fasta_by_id(aa_merged, aa_dedup)

        all_aa_parse = sorted(entry_dir.glob(f"*.seq.{bt}.blastdb.fa.parse.fa"))
        aa_parse_merged = entry_dir / f"{entry}.seq.{bt}.blastdb.fa.parse.merged.fa"
        merge_files(all_aa_parse, aa_parse_merged)
        parse_merged = entry_dir / f"{entry}.parse.merged.fa"
        dedup_fasta_by_id(aa_parse_merged, parse_merged)

    # Step 5: per-database merges → Orthofinder-ready copies
    orthof = entry_dir / "hits" / "orthofinder-input"
    ensure_dir(orthof)
    if blast_type == "tblastn":
        pattern = f"*.{{dbl}}.seq.{bt}.blastdb.translate.fa.parse.fa"
    else:
        pattern = f"*.{{dbl}}.seq.{bt}.blastdb.fa.parse.fa"

    for db in args.database:
        dbl = db_label(db)
        db_parts = sorted(entry_dir.glob(pattern.format(dbl=dbl)))
        db_merge = entry_dir / f"{dbl}.parse.merged.fa"
        merge_files(db_parts, db_merge)
        db_rmdup = entry_dir / f"{dbl}.parse.merged.rmdup.fa"
        dedup_fasta_by_id(db_merge, db_rmdup)
        shutil.copyfile(db_rmdup, orthof / dbl)

    # Step 6: optional add translations
    if args.add_seqs:
        optional_add_translations(entry, args.add_dbs, args.add_seqs, workdir)

    # Step 7: alignment and tree building
    print(f"Alignment & Tree Threads: {args.threads}")
    align_and_build_tree(entry, workdir, args.aligner, args.tree_builder, args.threads, args.mafft_mode)

    aln_fa = entry_dir / f"{entry}.parse.merged.aligned.fa"

    tree_src = Path(str(aln_fa) + ".nwk")
    tree_dst = entry_dir / "combinedtree.nwk"

    # Make sure downstream code knows the aligned FASTA path
    aligned_fasta_path = aln_fa

    if not aln_fa.exists():
        raise SystemExit(f"Aligned FASTA not found: {aln_fa}")
    if tree_src.exists():
        shutil.copyfile(tree_src, tree_dst)
    else:
        raise SystemExit(f"Expected tree file not found: {tree_src}")


    # Step 7.5: annotate motifs/HMMs (optional)
    hmm_dir = workdir / "hmm_files"

    hmm_paths = [str(hmm_dir / hmm_name) for hmm_name in args.hmms]

    if args.motifs or args.hmms:
        print(f"\n→ Annotating features (motifs/HMMs)")
        feats_path = annotate_features(entry, workdir,
                                       motifs_raw=args.motifs,
                                       motif_syntax=args.motif_syntax,
                                       motif_overlap=args.motif_overlap,
                                       hmm_files=hmm_paths)
        print(f"  Wrote {feats_path}")

    print(f"\n→ Generating tree visualizations (R/ggtree)")


    # Step 8: run visualize-tree.r
    visualize_tree(entry, args.queries, workdir, args.datasets)

    # Step 9: compact run-root outputs before archiving
    print(f"\n→ Cleaning run-root intermediates")
    cleanup_run_root(entry_dir, entry, args.queries, args.database, blast_type)

    # Step 10: archive into runs/<timestamp>
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    run_dir = archive_run(entry_dir, timestamp)

    print(f"\n{'='*50}")
    print(f"  Done.")
    print(f"{'='*50}")
    print(f"  Alignment: {run_dir / 'hits' / f'{entry}.parse.merged.aligned.fa'}")
    print(f"  Tree:      {run_dir / 'combinedtree.nwk'}")
    print(f"  Mapping:   {run_dir / 'merged_genome_mapping.txt'}")
    print(f"  Reports:   {run_dir}")
    subdir = f"runs/{timestamp}"
    write_arg = args.queries[0]
    print(f"\n  To re-draw trees (e.g. with a subnode):")
    rscript_path = _PACKAGE_DATA / "visualize_tree.r"
    redraw_name = f"{write_arg}_redraw"
    redraw_cmd = (
        f'  Rscript "{rscript_path}" -e {entry} -b {redraw_name} '
        f'--subdir "{subdir}" -n <NODE>'
    )
    if args.datasets:
        redraw_cmd += f' --datasets "{args.datasets}"'
    print(redraw_cmd)

if __name__ == "__main__":
    main()
