#!/usr/bin/env python3
"""
Extract sequences from FASTA file based on a list of gene IDs.

Usage:
  extract_seq.py <entry_dir> <gene> <list.csv> [aa_start] [aa_end]  - extract NT from .nt.parse.merged.fa
  extract_seq.py <entry_dir> <gene> <list.csv> --aa [aa_start] [aa_end]  - extract AA from .parse.merged.aligned.fa
"""
import argparse
from pathlib import Path
from Bio import SeqIO

RUN_ASSETS_DIRNAME = "genes_alignments_trees"

def existing_path(*candidates: Path) -> Path:
    for path in candidates:
        if path.exists():
            return path
    return candidates[0]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("entry", help="Entry directory path")
    parser.add_argument("gene", help="Gene name (used for input filename)")
    parser.add_argument("list", help="File with one gene ID per line")
    parser.add_argument("--aa", action="store_true", help="Extract aligned amino acids instead of nucleotides")
    parser.add_argument("aa_start", type=int, nargs='?', default=0)
    parser.add_argument("aa_end", type=int, nargs='?')
    args = parser.parse_args()

    entry_dir = Path(args.entry)
    hits_dir = entry_dir / "hits"
    assets_hits_dir = entry_dir / RUN_ASSETS_DIRNAME / "hits"

    if args.aa:
        fasta_path = existing_path(
            entry_dir / f"{args.gene}.parse.merged.aligned.fa",
            assets_hits_dir / f"{args.gene}.parse.merged.aligned.fa",
            hits_dir / f"{args.gene}.parse.merged.aligned.fa",
        )
        out_ext = ".aa.fa"
    else:
        fasta_path = existing_path(
            entry_dir / f"{args.gene}.nt.parse.merged.fa",
            assets_hits_dir / f"{args.gene}.nt.parse.merged.fa",
            hits_dir / f"{args.gene}.nt.parse.merged.fa",
        )
        out_ext = ".cds.fa"

    list_path = entry_dir / args.list
    out_path = entry_dir / f"{args.list}{out_ext}"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Read gene IDs
    with open(list_path) as f:
        genes = [line.strip() for line in f if line.strip()]

    # Index FASTA for efficient lookup
    seq_index = SeqIO.index(str(fasta_path), "fasta")

    # Write matching records in order
    with open(out_path, "w") as out_f:
        for gid in genes:
            if gid in seq_index:
                record = seq_index[gid]
                if args.aa_start > 0:
                    record.seq = record.seq[args.aa_start:args.aa_end]
                SeqIO.write(record, out_f, "fasta")

if __name__ == "__main__":
    main()
