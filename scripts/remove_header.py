#!/usr/bin/env python3
"""Remove FASTA descriptions, keeping only sequence names."""
import sys
from Bio import SeqIO

if len(sys.argv) != 3:
    print("Usage: remove_header.py <input.fa> <output.fa>")
    sys.exit(1)

input_file, output_file = sys.argv[1], sys.argv[2]

with open(output_file, "w") as out:
    for record in SeqIO.parse(input_file, "fasta"):
        record.description = ""
        SeqIO.write(record, out, "fasta")
