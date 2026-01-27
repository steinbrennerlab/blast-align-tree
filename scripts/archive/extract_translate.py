# Adam Steinbrenner
# astein10@uw.edu
# UW Seattle, Dept of Biology

# Find a sequence in a fasta file in subfolder "genomes" with matching fasta header to argument "entry". Optional arguments will pull a specified substring of amino acids

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

#Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("entry", help="name of gene given to original blast scripts")
parser.add_argument("query", help="query")
parser.add_argument("database")
parser.add_argument("aa_start", type=int, nargs='?', default=0)
parser.add_argument("aa_end", type=int, nargs='?')
args = parser.parse_args()

#Define directory
cwd = os.getcwd()

#Command line output of steps
print('\n')
print("Extract_Translate is using {} as entry".format(args.entry))
print("Extract_Translate is using {} as database".format(args.database))

#Use arguments to define files to search
filename = cwd + "/genomes/" + str(args.database)
output = cwd + "/" + str(args.entry) + "/" + str(args.query) + ".seq.fa"
print("Finding and translating this gene in the query db: {}".format(filename))


#Loop through all sequences in the specified fasta file database, and write the correct sequence to "query.seq.fa"
for seq_record in SeqIO.parse(filename, "fasta"):
	
	save_file = open(output, 'w+') #w+ for write, and plus for create if doesn't exist

	if args.query == seq_record.id:
		if args.aa_start > 0:
			save_file = open(output, 'a')
			save_file.write('>')
			save_file.write(seq_record.id)
			save_file.write('\n') #line break
			temp = str((translate(seq_record.seq)))
			save_file.write(temp[args.aa_start:args.aa_end])
			save_file.write('\n')
			print("Found query seq! " + args.query)
			print(seq_record.description)
			print(seq_record.translate())
			exit()
		else:
			save_file = open(output, 'a')
			save_file.write('>')
			save_file.write(seq_record.id)
			save_file.write('\n') #line break
			save_file.write(str((translate(seq_record.seq))))
			save_file.write('\n')
			print("Found query seq! " + args.query)
			print(seq_record.description)
			print(seq_record.translate())
			exit()

