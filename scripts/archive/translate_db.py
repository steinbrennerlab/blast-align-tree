from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import translate
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import argparse
import io
import os

#This script takes an argument -- gene name -- and formats it to receive the output from his blast pipeline (genename.seq.tblastn.blastdb.fa) and outputs a translated fasta


parser = argparse.ArgumentParser()
parser.add_argument("entry", help="name of gene given to original blast scripts")
parser.add_argument("query", help="query")
parser.add_argument("database")
args = parser.parse_args()
cwd = os.getcwd()

print("Translate script is using {} as entry".format(args.entry))
print("Translate script is using {} as database".format(args.database))
filename = cwd + "/" + str(args.entry) + "/" + str(args.query) + "." + str(args.database) + ".seq.tblastn.blastdb.fa"
output = cwd + "/" + str(args.entry) + "/" + str(args.query) + "." + str(args.database) + ".seq.tblastn.blastdb.translate.fa"
print("Translating {}".format(filename))

for seq_record in SeqIO.parse(filename, "fasta"):
	#print(seq_record.translate())
	#print(len(seq_record))
	#save_file = open("file{}.xml".format(len(seq_record), 'w'))
	#s = seq_record.id
	save_file = open(output, 'a')
	save_file.write('>')
	#save_file.write(seq_record.id)
	save_file.write(seq_record.description)

	#save_file.write('\n') #line break
	#save_file.write(str(seq_record.translate()))

	save_file.write('\n') #line break
	save_file.write(str(translate(seq_record.seq)))
	#word = str(seq_record.seq)
	#save_file.write(word[-100:])
	save_file.write('\n')
	#save_file.write('/n')
	save_file.close()
	#print(save_file)
print(seq_record.id)