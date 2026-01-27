#Trying to pull fasta nucleotide entries based on a text list, with optional aa start and end points

#Takes "entry" the entry gene and subfolder and "list" a filename containing the cds of interest.
#Example: python extract_seq.py gene_ID list.txt

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import translate
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import argparse
import io
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument("entry", help="name of gene given to original blast scripts")
parser.add_argument("list", help="name of the genes list, a text file with one gene on each line.  As of 9/8/18, tree.R can generate this")
parser.add_argument("aa_start", type=int, nargs='?', default=0)
parser.add_argument("aa_end", type=int, nargs='?')
args = parser.parse_args()
cwd = os.getcwd()

print('\n')
print("Extract_seqs is using {} as entry".format(args.entry))
filename = cwd + "/" + str(args.entry) + "/" + str(args.entry) + ".nt.parse.merged.fa"
seqlist = cwd + "/" + str(args.entry) + "/output/" + str(args.list)
output = cwd + "/" + str(args.entry) + "/output/" + str(args.list) + ".cds.fa"
print("Finding and translating this gene in the query db: {}".format(filename))

print(seqlist)

with open(seqlist) as f:
	for line2 in f:
		line = line2.rstrip('\n') #Need to strip the line break
		#print(line)
		for seq_record in SeqIO.parse(filename, "fasta"):
			#print(seq_record.id)
			save_file = open(output, 'a+')

			if line == seq_record.id:
				if args.aa_start > 0:
					save_file = open(output, 'a')
					save_file.write('>')
					save_file.write(seq_record.id)
					save_file.write('\n') #line break
					temp = str((seq_record.seq))
					save_file.write(temp[args.aa_start:args.aa_end])
					save_file.write('\n')
					print("Found entry seq! " + args.entry)
					print(seq_record.description)
					print(seq_record.seq())
					save_file.close()
				else:
					save_file = open(output, 'a')
					save_file.write('>')
					save_file.write(seq_record.id)
					save_file.write('\n')
					save_file.write(str((seq_record.seq)))
					save_file.write('\n')

					save_file.close()

