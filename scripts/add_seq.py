#!/usr/bin/env python3
"""
Append a raw nucleotide sequence (no translation) to the merged FASTA and coding table.

Usage: python add_seq.py <entry> <database> <seq_id>
"""
import argparse
import os
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Append a sequence to merged outputs")
    parser.add_argument("entry", help="Entry folder name")
    parser.add_argument("database", help="Database filename in genomes/")
    parser.add_argument("addseq", help="Sequence ID to add")
    args = parser.parse_args()

    cwd = os.getcwd()
    db_path = os.path.join(cwd, "genomes", args.database)
    fa_out = os.path.join(cwd, args.entry, f"{args.entry}.parse.merged.fa")
    coding_out = os.path.join(cwd, args.entry, "merged_coding.txt")

    for record in SeqIO.parse(db_path, "fasta"):
        if record.id == args.addseq:
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
