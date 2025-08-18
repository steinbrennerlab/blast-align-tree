#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("entry", help="name of gene folder")
    parser.add_argument("list", help="file with one gene ID per line")
    parser.add_argument("aa_start", type=int, nargs='?', default=0)
    parser.add_argument("aa_end", type=int, nargs='?')
    args = parser.parse_args()

    fasta_path = f"{args.entry}/{args.entry}.nt.parse.merged.fa"
    list_path  = f"{args.entry}/output/{args.list}"
    out_path   = f"{args.entry}/output/{args.list}.cds.fa"

    # 1) Read gene IDs
    with open(list_path) as f:
        genes = [line.strip() for line in f if line.strip()]

    # 2) Index FASTA
    seq_index = SeqIO.index(fasta_path, "fasta")

    # 3) Open output once and write matching records
    with open(out_path, "w") as out_f:
        for gid in genes:
            if gid in seq_index:
                record = seq_index[gid]
                seq = record.seq
                if args.aa_start > 0:
                    subseq = seq[args.aa_start:args.aa_end]
                    record = record[:0]  # drop existing seq
                    record.seq = subseq
                SeqIO.write(record, out_f, "fasta")

if __name__ == "__main__":
    main()
