#!/usr/bin/env python3
"""
Extract sequences from FASTA file based on a list of gene IDs.

Usage:
  extract_seq.py <entry_dir> <gene> <list.csv> [aa_start] [aa_end]  - extract NT from .nt.parse.merged.fa
  extract_seq.py <entry_dir> <gene> <list.csv> --aa [aa_start] [aa_end]  - extract AA from .parse.merged.clustal.fa
"""
import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("entry", help="Entry directory path")
    parser.add_argument("gene", help="Gene name (used for input filename)")
    parser.add_argument("list", help="File with one gene ID per line")
    parser.add_argument("--aa", action="store_true", help="Extract aligned amino acids instead of nucleotides")
    parser.add_argument("aa_start", type=int, nargs='?', default=0)
    parser.add_argument("aa_end", type=int, nargs='?')
    args = parser.parse_args()

    if args.aa:
        fasta_path = f"{args.entry}/{args.gene}.parse.merged.clustal.fa"
        out_ext = ".aa.fa"
    else:
        fasta_path = f"{args.entry}/{args.gene}.nt.parse.merged.fa"
        out_ext = ".cds.fa"

    list_path = f"{args.entry}/output/{args.list}"
    out_path = f"{args.entry}/output/{args.list}{out_ext}"

    # Read gene IDs
    with open(list_path) as f:
        genes = [line.strip() for line in f if line.strip()]

    # Index FASTA for efficient lookup
    seq_index = SeqIO.index(fasta_path, "fasta")

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
