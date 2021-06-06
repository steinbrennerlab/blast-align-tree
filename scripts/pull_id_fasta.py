#Trying to parse a fasta file to give whatever input header thing that I want

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
parser.add_argument("headerword", help="string right before what I'm trying to pull from a fasta header")
parser.add_argument("trailing", help="trailing filename to merge")
args = parser.parse_args()
cwd = os.getcwd()

print("Pull ID script is using {} as entry".format(args.entry))
print("Pull ID script is using {} as database".format(args.database))
print("Pull ID script is looking for word number {} in fasta headers".format(args.headerword))
filename = cwd + "/" + str(args.entry) + "/" + str(args.entry) + "." + str(args.database) + ".seq.tblastn.blastdb.fa"
filename = cwd + "/" + str(args.entry) + "/" + str(args.entry) + "." + str(args.database) + ".seq.tblastn.blastdb.translate.fa"
filename = cwd + "/" + str(args.entry) + "/" + str(args.entry) + "." + str(args.database) + str(args.trailing)
print("Filename is {}".format(filename))
#filename = "sbic.fa" #tester entry
output = cwd + "/" + str(args.entry) + "/" + str(args.entry) + "." + str(args.database) + ".seq.tblastn.blastdb.translate.parse.fa"
output = cwd + "/" + str(args.entry) + "/" + str(args.entry) + "." + str(args.database) + str(args.trailing) + ".parse.fa"
#output = "sbic_pull.fa" #tester output
print("Output is {}".format(output))
print("Parsing {}".format(filename))
print(args.headerword)
header=str(args.headerword) #makes a string out of the argument headerword
print(header)

for seq_record in SeqIO.parse(filename, "fasta"):

	save_file = open(output, 'a')
	save_file.write('>')

	if header == "id":
		save_file.write(seq_record.id)

	elif header not in str(seq_record.description):
		save_file.write(seq_record.id)
	else:
		word = ((re.split(header,seq_record.description))[1]) #Takes above, this becomes a string which can be split
		save_file.write(word.split(' ',1)[0]) #Splits the description to write only what's between header (args.headerword) and the next space
	save_file.write('\n') #line break
	save_file.write(str(seq_record.seq))
	save_file.write('\n')