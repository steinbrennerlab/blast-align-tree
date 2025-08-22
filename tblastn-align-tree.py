#!/usr/bin/env python3
"""
Python rewrite of tblastn-align-tree.sh

- Orchestrates BLAST+ CLI tools via subprocess
- Replaces awk/xargs/sed with native Python
- Parallelizes (query, database) jobs
- Requires: tblastn, blastdbcmd, clustalo, FastTree, Python helper scripts in ./scripts

Usage example:
python tblastn_align_tree.py \
  -q ENTRY Q2 \
  -qdbs ENTRYDB.fa Q2DB.fa \
  -n 50 50 \
  -hdr '>' 'gene=' \
  -dbs Vung469.cds.fa TAIR10cds.fa \
  -add OUTGROUP1 \
  -add_db OUTGROUP_DB.fa \
  -aa 10 200
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

# -----------------------
# Utilities
# -----------------------



def check_tool(name: str):
    if shutil.which(name) is None:
        raise SystemExit(f"Required tool not found in PATH: {name}")

def run(cmd: List[str], cwd: Optional[Path] = None, capture: bool = False) -> subprocess.CompletedProcess:
    print("RUN:", " ".join(cmd))
    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        check=True,
        text=True,
        capture_output=capture
    )

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
        cur_header = None
        for line in fin:
            if line.startswith(">"):
                # Header token: everything until whitespace
                tok = line.strip().split()[0]
                if tok not in seen:
                    seen.add(tok)
                    write = True
                    cur_header = tok
                    fout.write(line)
                else:
                    write = False
                    cur_header = None
            else:
                if write:
                    fout.write(line)

def prepend_header_line(fp: Path, header: str):
    """Prepend a single header line to a file (like sed -i '1s/^/.../')."""
    if not fp.exists():
        return
    content = fp.read_text(encoding="utf-8", errors="ignore")
    fp.write_text(header + content, encoding="utf-8")

from datetime import datetime
from pathlib import Path
import shutil

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



# -----------------------
# Core pipeline steps
# -----------------------

def extract_and_translate(entry: str, q: str, qdb: str, slice_args: Optional[List[str]], scripts_dir: Path, workdir: Path):
    cmd = ["python", str(scripts_dir / "extract_translate.py"), entry, q, qdb]
    if slice_args and len(slice_args) >= 2:
        cmd += slice_args[:2]
    run(cmd, cwd=workdir)

def tblastn_and_post(entry: str, q: str, db: str, max_targets: str, workdir: Path) -> Tuple[str, str, str]:
    """
    Run tblastn for a (query, db) pair and return paths to:
      - sseqid list
      - full report
      - blastdb fasta (with stops removed)
    """
    q_fa = workdir / entry / f"{q}.seq.fa"
    db_path = workdir / "genomes" / db
    out_base = workdir / entry / f"{q}.{db}.seq.tblastn"
    # outfmt 6 list of sseqids
    run([
        "tblastn", "-query", str(q_fa), "-db", str(db_path),
        "-max_target_seqs", str(max_targets), "-max_hsps", "1",
        "-outfmt", "6 sseqid", "-out", str(out_base)
    ])
    # full report
    full = out_base.with_suffix(out_base.suffix + ".full")
    run([
        "tblastn", "-query", str(q_fa), "-db", str(db_path),
        "-max_hsps", "1", "-outfmt", "1", "-out", str(full)
    ])
    # Add header line like the sed command in bash version
    prepend_header_line(full, "hit query_id\tsubject_id\tpct_identity\taln_length\tn_of_mismatches\tgap_openings\tq_start q_end\ts_start   s_end\te_value bit_score\n")

    # blastdbcmd to pull sequences; then remove stop codons
    stop_fa = Path(str(out_base) + ".blastdb.stop.fa")
    nt_fa   = Path(str(out_base) + ".blastdb.fa")
    run(["blastdbcmd", "-db", str(db_path), "-entry_batch", str(out_base), "-out", str(stop_fa)])
    run(["python", str(workdir / "scripts" / "remove_stop.py"), str(stop_fa), str(nt_fa)])

    return str(out_base), str(full), str(nt_fa)

def translate_and_parse_headers(entry: str, q: str, db: str, header_rule: str, scripts_dir: Path, workdir: Path):
    # translate_db.py
    run(["python", str(scripts_dir / "translate_db.py"), entry, q, db], cwd=workdir)
    # pull_id_fasta for nt and translated
    run(["python", str(scripts_dir / "pull_id_fasta.py"), entry, q, db, header_rule, ".seq.tblastn.blastdb.fa"], cwd=workdir)
    run(["python", str(scripts_dir / "pull_id_fasta.py"), entry, q, db, header_rule, ".seq.tblastn.blastdb.translate.fa"], cwd=workdir)
    # coding table
    run(["python", str(scripts_dir / "pull_id_fasta_coding.py"), entry, q, db, header_rule, ".seq.tblastn.blastdb.translate.fa"], cwd=workdir)

def optional_add_translations(entry: str, add_dbs: List[str], add_seqs: List[str], scripts_dir: Path, workdir: Path):
    for db, seq in zip(add_dbs, add_seqs):
        run(["python", str(scripts_dir / "add_translations.py"), entry, db, seq], cwd=workdir)

def clustalo_and_fasttree(entry: str, workdir: Path):
    in_fa  = workdir / entry / f"{entry}.parse.merged.fa"
    aln_fa = workdir / entry / f"{entry}.parse.merged.clustal.fa"
    run(["clustalo", "-i", str(in_fa), "-o", str(aln_fa)])
    tree_out = Path(str(aln_fa) + ".nwk")
    print("Running FastTree")
    run(["FastTree", "-out", str(tree_out), str(aln_fa)])
    shutil.copyfile(tree_out, workdir / entry / "combinedtree.nwk")
    
def visualize_tree(entry: str, queries: List[str], workdir: Path):
    """
    Run visualize-tree.r on the combinedtree.nwk 
    By default, --write argument is set to the first query name.
    """
    rscript = workdir / "visualize-tree.r"

    if not rscript.exists():
        raise SystemExit(f"R script not found: {rscript}")

    write_arg = queries[0]  # default = first -q
    cmd = [
        "Rscript", str(rscript),
        "--entry", entry,
        "--write", write_arg
    ]
    run(cmd, cwd=workdir)


# -----------------------
# Main
# -----------------------

def main():
    ap = argparse.ArgumentParser(description="Python rewrite of tblastn-align-tree.sh")
    ap.add_argument("-q", "--queries", nargs="+", required=True, help="fasta ID for query sequences")
    ap.add_argument("-qdbs", "--query_databases", nargs="+", required=True, help="fasta files containing the queries")
    ap.add_argument("-n", "--seqs", nargs="+", required=True, help="-max_target_seqs per search database (align with -dbs)")
    ap.add_argument("-hdr", "--header", nargs="+", required=True, help="header parsing rules per db (align with -dbs)")
    ap.add_argument("-dbs", "--database", nargs="+", required=True, help="blast databases to search (filenames under ./genomes)")
    ap.add_argument("-add", "--add_seqs", nargs="*", default=[], help="additional sequences (optional)")
    ap.add_argument("-add_db", "--add_dbs", nargs="*", default=[], help="databases for additional sequences (optional)")
    ap.add_argument("-aa", "--slice", nargs="*", default=[], help="AA slice start end, optional")
    ap.add_argument("--threads", type=int, default=max(1, os.cpu_count() // 2), help="parallel jobs for BLAST steps")
    ap.add_argument("--workdir", default=".", help="working directory (default=.)")
    args = ap.parse_args()

    # Sanity checks
    if len(args.query_databases) != len(args.queries):
        raise SystemExit("queries and query_databases lengths must match")
    if len(args.seqs) != len(args.database):
        raise SystemExit("seqs (-n) must have same length as database (-dbs)")
    if len(args.header) != len(args.database):
        raise SystemExit("header (-hdr) must have same length as database (-dbs)")
    if args.add_seqs and (len(args.add_seqs) != len(args.add_dbs)):
        raise SystemExit("add_seqs and add_dbs must have the same length")

    # External tool checks
    for tool in ["tblastn", "blastdbcmd", "clustalo", "FastTree"]:
        check_tool(tool)

    workdir = Path(args.workdir).resolve()
    scripts_dir = workdir / "scripts"

    entry = args.queries[0]
    entrydb = args.query_databases[0]  # kept for parity with bash; extract_translate uses both

    # Create directories
    entry_dir = workdir / entry
    entry_root = workdir / entry
    ensure_dir(entry_dir / "output")
    print(f"Making directory based on first query: {entry}")
    print(f"First database to search, entrydb: {entrydb}")
    print("All queries:", *args.queries, sep="\n  ")
    print("All source databases for queries:", *args.query_databases, sep="\n  ")
    print("-dbs, Queries will be BLASTed against all Databases:", *args.database, sep="\n  ")
    print("Number of subject seqs to pull from tblastn search:", *args.seqs, sep="\n  ")
    print(f"dbs length is {len(args.database)}")

    # Step 1: extract and translate each query (optionally AA slice)
    for q, qdb in zip(args.queries, args.query_databases):
        extract_and_translate(entry, q, qdb, args.slice if args.slice else None, scripts_dir, workdir)

    # Step 2: parallel BLAST per (query, db)
    jobs: List[Tuple[str, str, str, str]] = []  # (q, db, n, header_rule)
    for q in args.queries:
        for db, n, hdr in zip(args.database, args.seqs, args.header):
            jobs.append((q, db, n, hdr))

    def _one(job: Tuple[str, str, str, str]):
        q, db, n, hdr = job
        out_base, full, nt_fa = tblastn_and_post(entry, q, db, n, workdir)
        translate_and_parse_headers(entry, q, db, hdr, scripts_dir, workdir)
        return job

    with cf.ThreadPoolExecutor(max_workers=args.threads) as ex:
        list(ex.map(_one, jobs))

    # Step 3: combine coding txt → merged_genome_mapping.txt and prepend header
    coding_txts = sorted(entry_dir.glob("*.txt"))  # created by pull_id_fasta_coding.py
    merged_genome_mapping = entry_dir / "merged_genome_mapping.txt"
    merge_files(coding_txts, merged_genome_mapping)
    prepend_header_line(merged_genome_mapping, "taxa\tgenome\n")

    # Step 4: merge and dedup FASTAs (nt and translated + parsed variants)
    # Original script produces multiple merged intermediates; we mirror file names.
    # All original header translated fasta files → $ENTRY.seq.tblastn.blastdb.merged.fa → $ENTRY.nt.merged.fa (dedup)
    all_nt = sorted(entry_dir.glob("*.seq.tblastn.blastdb.fa"))
    nt_merged = entry_dir / f"{entry}.seq.tblastn.blastdb.merged.fa"
    merge_files(all_nt, nt_merged)
    nt_dedup = entry_dir / f"{entry}.nt.merged.fa"
    dedup_fasta_by_id(nt_merged, nt_dedup)

    # All parsed nt fastas → merged → dedup
    all_nt_parse = sorted(entry_dir.glob("*.seq.tblastn.blastdb.nt.parse.merged.fa"))  # may not exist yet; build below
    # But in bash the sources are "*.seq.tblastn.blastdb.fa.parse.fa"
    all_nt_parse_sources = sorted(entry_dir.glob("*.seq.tblastn.blastdb.fa.parse.fa"))
    nt_parse_merged_src = entry_dir / f"{entry}.seq.tblastn.blastdb.nt.parse.merged.fa"
    merge_files(all_nt_parse_sources, nt_parse_merged_src)
    nt_parse_dedup = entry_dir / f"{entry}.nt.parse.merged.fa"
    dedup_fasta_by_id(nt_parse_merged_src, nt_parse_dedup)

    # All translated header fastas → merged → dedup → simplify
    all_trans = sorted(entry_dir.glob("*.seq.tblastn.blastdb.translate.fa"))
    trans_merged = entry_dir / f"{entry}.seq.tblastn.blastdb.translate.merged.fa"
    merge_files(all_trans, trans_merged)
    trans_dedup = entry_dir / f"{entry}.merged.fa"
    dedup_fasta_by_id(trans_merged, trans_dedup)

    # All translated parsed → merged → dedup → final parse.merged.fa for alignment
    all_trans_parse = sorted(entry_dir.glob("*.seq.tblastn.blastdb.translate.fa.parse.fa"))
    trans_parse_merged = entry_dir / f"{entry}.seq.tblastn.blastdb.translate.fa.parse.merged.fa"
    merge_files(all_trans_parse, trans_parse_merged)
    parse_merged = entry_dir / f"{entry}.parse.merged.fa"
    dedup_fasta_by_id(trans_parse_merged, parse_merged)

    # Step 5: per-database merges → orthofinder-input copies
    orthof = entry_dir / "output" / "orthofinder-input"
    ensure_dir(orthof)
    for db in args.database:
        # gather all *.{db}.seq.tblastn.blastdb.translate.fa.parse.fa
        db_parts = sorted(entry_dir.glob(f"*.{db}.seq.tblastn.blastdb.translate.fa.parse.fa"))
        db_merge = entry_dir / f"{db}.parse.merged.fa"
        merge_files(db_parts, db_merge)
        db_rmdup = entry_dir / f"{db}.parse.merged.rmdup.fa"
        dedup_fasta_by_id(db_merge, db_rmdup)
        shutil.copyfile(db_rmdup, orthof / db)

    # Step 6: optional add translations
    if args.add_seqs:
        optional_add_translations(entry, args.add_dbs, args.add_seqs, scripts_dir, workdir)

    # Step 7: clustalo and FastTree
    clustalo_and_fasttree(entry, workdir)

    print("\nDone.")
    print(f"Key outputs under: {entry_dir}")
    print(f"  Alignment: {entry_dir / (entry + '.parse.merged.clustal.fa')}")
    print(f"  Tree:      {entry_dir / 'combinedtree.nwk'}")
    print(f"  Mapping:    {merged_genome_mapping}")
    
    # Step 8: run visualize-tree.r
    visualize_tree(entry, args.queries, workdir)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    archive_run(entry_root, timestamp)



if __name__ == "__main__":
    main()
