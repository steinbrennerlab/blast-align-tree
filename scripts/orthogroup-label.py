#!/usr/bin/env python3
"""
Map sequence IDs to OrthoFinder orthogroup assignments.

Usage: python orthogroup-label.py <entry> <list.txt>
"""
import argparse
import os
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Map sequences to OrthoFinder orthogroups")
    parser.add_argument("entry", help="Entry folder name")
    parser.add_argument("list", help="Text file with one gene ID per line")
    args = parser.parse_args()

    cwd = os.getcwd()
    og_dir = os.path.join(cwd, args.entry, "output", "Orthogroup_Sequences")
    seq_list = os.path.join(cwd, args.entry, "output", args.list)
    output = os.path.join(cwd, args.entry, "output", f"{args.list}.orthogroups.txt")

    # Read gene IDs
    with open(seq_list) as f:
        genes = [line.strip() for line in f if line.strip()]

    # Build mapping: gene_id -> orthogroup filename
    og_map = {}
    for og_file in os.listdir(og_dir):
        og_path = os.path.join(og_dir, og_file)
        for record in SeqIO.parse(og_path, "fasta"):
            if record.id in genes:
                og_map[record.id] = og_file

    # Write output
    with open(output, "w") as out:
        out.write("taxa\tOG\tOG_type\n")
        for gene in genes:
            if gene in og_map:
                out.write(f"{gene}\t{og_map[gene]}\tnonunique-OG\n")

if __name__ == "__main__":
    main()
