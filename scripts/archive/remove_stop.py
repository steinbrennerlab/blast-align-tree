from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import io
import os
import re

#Remove stop codon from all sequences in the entry file by looking for stop codons in the last three nucleotides
parser = argparse.ArgumentParser()
parser.add_argument("entry", help="name of gene given to original blast scripts")
parser.add_argument("output")
args = parser.parse_args()

codon_stop_array=["TAG","TGA","TAA"]

my_record=[]
print("Removing stop from "+args.entry)

for record in SeqIO.parse(args.entry, "fasta"):
	#print(record.seq)
	tempRecordSeq = list(record.seq)
	for index in range(0, len(record.seq), 3):
		codon = record.seq[index:index+3]
		if codon in codon_stop_array:
			del tempRecordSeq[index:index+3]
	record.seq = Seq("".join(tempRecordSeq))
	my_record.append(record)
	

SeqIO.write(my_record,args.output,"fasta")