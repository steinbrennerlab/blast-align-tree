#April 2019: editing extract_translate.py to pull a list from a file
#Usage: python extract_translate.py $ENTRY.fa $ENTRYDB ${SLICE[0]} ${SLICE[1]}

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
parser.add_argument("database")
parser.add_argument("aa_start", type=int, nargs='?', default=0)
parser.add_argument("aa_end", type=int, nargs='?')
args = parser.parse_args()
cwd = os.getcwd()

print('\n')
print("Extract_Translate is using {} as entry".format(args.entry))
print("Extract_Translate is using {} as database".format(args.database))
filename = cwd + "/" + str(args.entry)
database = cwd + "/genomes/" + str(args.database)
output = cwd + "/" + str(args.entry) + ".seq.fa"
print("Finding and translating this gene in the query db: {}".format(filename))

#takes in the "entry" file as a tabular sequence file; two columns, first is name second is sequence (leave it blank if you want)
for seq_record in SeqIO.parse(filename, "tab"):

	#open file in append mode
	save_file = open(output, 'a+')
	print(seq_record.id)
	n=0
	for seq_record2 in SeqIO.parse(database, "fasta"):
		if seq_record.id == seq_record2.id:
			if args.aa_start > 0:
				save_file = open(output, 'a')
				save_file.write('>')
				save_file.write(seq_record.id)
				save_file.write('\n') #line break
				temp = str((translate(seq_record.seq)))
				save_file.write(temp[args.aa_start:args.aa_end])
				save_file.write('\n')
				print("Found entry seq! " + args.entry)
				print(seq_record.description)
				print(seq_record.translate())
			else:
				print("Found:")
				print(seq_record.id)
				save_file = open(output, 'a')
				save_file.write(seq_record2.format("fasta")) 

				
				print(seq_record.description)
				print(seq_record.translate())
				n=1
	if n==0:
		print("Could not find:")
		print(seq_record.id)
		save_file.write('>')
		save_file.write(seq_record.id)
		save_file.write('\n') #line break

