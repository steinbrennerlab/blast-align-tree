#!/usr/bin/env python3
"""
Append a raw nucleotide sequence (no translation) to the merged FASTA and coding table.

Usage: python add_seq.py <entry> <database> <seq_id>
"""
import argparse
import os
from Bio import SeqIO

def first_existing(*paths):
    for path in paths:
        if os.path.exists(path):
            return path
    return paths[0]

def main():
    parser = argparse.ArgumentParser(description="Append a sequence to merged outputs")
    parser.add_argument("entry", help="Entry folder name")
    parser.add_argument("database", help="Database filename in genomes/")
    parser.add_argument("addseq", help="Sequence ID to add")
    args = parser.parse_args()

    cwd = os.getcwd()
    entry_dir = os.path.join(cwd, args.entry)
    hits_dir = os.path.join(entry_dir, "hits")
    db_path = os.path.join(cwd, "genomes", args.database)
    fa_out = first_existing(
        os.path.join(entry_dir, f"{args.entry}.parse.merged.fa"),
        os.path.join(hits_dir, f"{args.entry}.parse.merged.fa"),
    )
    coding_out = first_existing(
        os.path.join(entry_dir, "merged_genome_mapping.txt"),
        os.path.join(entry_dir, "merged_coding.txt"),
    )

    for record in SeqIO.parse(db_path, "fasta"):
        if record.id == args.addseq:
            os.makedirs(os.path.dirname(fa_out), exist_ok=True)
            with open(fa_out, "a") as fa:
                fa.write(f">{record.id}\n{record.seq}\n")
            with open(coding_out, "a") as txt:
                txt.write(f"\n{record.id}\t{args.database}")
            print(f"Added {args.addseq} from {args.database}")
            return

    print(f"Sequence {args.addseq} not found in {db_path}")
    exit(1)

if __name__ == "__main__":
    main()
